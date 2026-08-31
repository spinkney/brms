library(brms)
library(cmdstanr)

root <- Sys.getenv("FM_BENCHMARK_DIR")
if (!nzchar(root)) {
  stop("Set FM_BENCHMARK_DIR to a writable benchmark working directory.")
}
repo <- Sys.getenv("NY_R_TALK_DIR", file.path(root, "ny_r_talk"))
draw_dir <- file.path(root, "draws")

load(file.path(repo, "simulated_data.rdata"))

prepare_dataset <- function(sim) {
  d <- sim[[1]]
  X <- d$X
  if (d$N < d$J) {
    X <- X[, 2:1, drop = FALSE]
    tmp <- d$N
    d$N <- d$J
    d$J <- tmp
  }
  X <- unname(X)
  list(
    dat = data.frame(
      y = d$y,
      group_1 = factor(X[, 1], levels = seq_len(d$N)),
      group_2 = factor(X[, 2], levels = seq_len(d$J))
    ),
    source = d
  )
}

datasets <- list(
  simulation_1 = prepare_dataset(fm_simulation_1),
  simulation_2 = prepare_dataset(fm_simulation_2)
)
seeds <- c(simulation_1 = 123, simulation_2 = 216)

# These lognormal parameters approximate the full-table RMS induced by the
# source model's centered iid N(0, 1) rank-five factors at 100 x 20. The
# spectrum prior is necessarily different: the exact source law is a
# product-Wishart spectrum, whereas fm() uses ordered uniform energy gaps.
rms_log_mean <- 0.7678554
rms_log_sd <- 0.08103256
d0 <- datasets[[1]]$source
intercept_scale <- sqrt(
  d0$beta_sigma^2 * (1 / d0$N + 1 / d0$J) +
    d0$K / (d0$N * d0$J)
)

matched_prior <- c(
  prior_string(
    sprintf("normal(0, %.17g)", intercept_scale),
    class = "Intercept"
  ),
  # Canonical source marginal SDs are 3.026136 and 2.932149. A common scale
  # of three is a close, symmetric approximation, not an exact joint match.
  prior(constant(3), class = "sdfm_main", group = "group_1"),
  prior(constant(3), class = "sdfm_main", group = "group_2"),
  prior_string(
    sprintf("lognormal(%.17g, %.17g)", rms_log_mean, rms_log_sd),
    class = "sdfm", coef = "group_1:group_2"
  ),
  prior(constant(1), class = "sigma")
)

for (nm in names(datasets)) {
  dir.create(
    file.path(draw_dir, paste0("semi_", nm)),
    showWarnings = FALSE, recursive = TRUE
  )
}

message("Compiling and sampling semi-orthogonal FM: simulation_1")
semi_1 <- brm(
  y ~ fm(group_1, group_2, k = 5),
  data = datasets$simulation_1$dat,
  family = gaussian(),
  prior = matched_prior,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  seed = unname(seeds["simulation_1"]),
  control = list(adapt_delta = 0.8, max_treedepth = 20),
  refresh = 100,
  save_pars = save_pars(all = TRUE),
  output_dir = file.path(draw_dir, "semi_simulation_1"),
  silent = 0
)
saveRDS(semi_1, file.path(root, "fit_semi_simulation_1.rds"))

message("Reusing the semi-orthogonal executable for simulation_2")
semi_2 <- update(
  semi_1,
  newdata = datasets$simulation_2$dat,
  recompile = FALSE,
  seed = unname(seeds["simulation_2"]),
  output_dir = file.path(draw_dir, "semi_simulation_2"),
  file = NULL,
  refresh = 100,
  silent = 0
)
saveRDS(semi_2, file.path(root, "fit_semi_simulation_2.rds"))

saveRDS(
  list(
    seeds = seeds,
    rms_log_mean = rms_log_mean,
    rms_log_sd = rms_log_sd,
    intercept_scale = intercept_scale,
    cmdstan_version = as.character(cmdstan_version())
  ),
  file.path(root, "semi_benchmark_config.rds")
)

message("Semi-orthogonal fits complete")
