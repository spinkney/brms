library(cmdstanr)

root <- Sys.getenv("FM_BENCHMARK_DIR")
if (!nzchar(root)) {
  stop("Set FM_BENCHMARK_DIR to a writable benchmark working directory.")
}
repo <- Sys.getenv("NY_R_TALK_DIR", file.path(root, "ny_r_talk"))
source_file <- file.path(repo, "simulated_data.rdata")
if (!file.exists(source_file)) {
  stop("Cannot find simulated_data.rdata in NY_R_TALK_DIR.")
}
load(source_file)

prepare_dataset <- function(sim) {
  d <- sim[[1]]
  X <- d$X
  if (d$N < d$J) {
    X <- X[, 2:1, drop = FALSE]
    tmp <- d$N
    d$N <- d$J
    d$J <- tmp
  }
  list(
    N = d$N,
    J = d$J,
    K = d$K,
    X = unname(X),
    y = d$y,
    beta_sigma = d$beta_sigma,
    y_sigma = d$y_sigma
  )
}

datasets <- list(
  simulation_1 = prepare_dataset(fm_simulation_1),
  simulation_2 = prepare_dataset(fm_simulation_2)
)
seeds <- c(simulation_1 = 123, simulation_2 = 216)
draw_dir <- file.path(root, "draws")
dir.create(draw_dir, showWarnings = FALSE, recursive = TRUE)
exe_dir <- file.path(root, "executables")
dir.create(exe_dir, showWarnings = FALSE, recursive = TRUE)

model_file <- file.path(
  "analysis", "fm-s2z-diagnostics", "benchmark", "stan_fm_1_modern.stan"
)
model <- cmdstan_model(
  model_file,
  exe_file = file.path(exe_dir, "stan_fm_1_modern")
)
for (nm in names(datasets)) {
  output_dir <- file.path(draw_dir, paste0("repo_", nm))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  fit <- model$sample(
    data = datasets[[nm]],
    seed = unname(seeds[nm]),
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.8,
    max_treedepth = 20,
    refresh = 100,
    save_warmup = FALSE,
    output_dir = output_dir
  )
  fit$save_object(file.path(root, paste0("fit_repo_", nm, ".rds")))
}
