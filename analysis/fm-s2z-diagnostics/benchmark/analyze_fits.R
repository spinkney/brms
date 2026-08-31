library(cmdstanr)
library(loo)
library(matrixStats)
library(posterior)

options(mc.cores = 4L)

root <- Sys.getenv("FM_BENCHMARK_DIR")
if (!nzchar(root)) {
  stop("Set FM_BENCHMARK_DIR to the benchmark working directory.")
}
repo <- Sys.getenv("NY_R_TALK_DIR", file.path(root, "ny_r_talk"))

load(file.path(repo, "simulated_data.rdata"))

prepare_dataset <- function(sim) {
  d <- sim[[1]]
  X <- d$X
  swapped <- d$N < d$J
  if (swapped) {
    X <- X[, 2:1, drop = FALSE]
    tmp <- d$N
    d$N <- d$J
    d$J <- tmp
  }
  X <- unname(X)
  theta_true_obs <- as.numeric(
    sim[[2]]$linear_predictor + sim[[2]]$factor_terms
  )
  cell_index <- X[, 1] + (X[, 2] - 1L) * d$N
  theta_true <- numeric(d$N * d$J)
  theta_true[cell_index] <- theta_true_obs
  list(
    N = d$N,
    J = d$J,
    K = d$K,
    X = X,
    y = as.numeric(d$y),
    cell_index = cell_index,
    theta_true = theta_true,
    theta_true_obs = theta_true_obs,
    swapped = swapped
  )
}

datasets <- list(
  simulation_1 = prepare_dataset(fm_simulation_1),
  simulation_2 = prepare_dataset(fm_simulation_2)
)

read_csv_fit <- function(kind, dataset) {
  csv <- Sys.glob(file.path(root, "draws", paste0(kind, "_", dataset), "*.csv"))
  stopifnot(length(csv) == 4L)
  as_cmdstan_fit(sort(csv))
}

draw_matrix <- function(fit, variable) {
  as.matrix(fit$draws(variables = variable, format = "draws_matrix"))
}

chain_layout <- function(fit) {
  x <- as.data.frame(fit$draws(variables = "lp__", format = "draws_df"))
  list(chain = x$.chain, iteration = x$.iteration)
}

as_derived_draws <- function(x, layout, prefix) {
  chains <- sort(unique(layout$chain))
  niter <- max(layout$iteration)
  out <- array(
    NA_real_,
    dim = c(niter, length(chains), ncol(x)),
    dimnames = list(
      iteration = as.character(seq_len(niter)),
      chain = as.character(chains),
      variable = paste0(prefix, "[", seq_len(ncol(x)), "]")
    )
  )
  for (j in seq_along(chains)) {
    take <- which(layout$chain == chains[j])
    take <- take[order(layout$iteration[take])]
    out[, j, ] <- x[take, , drop = FALSE]
  }
  as_draws_array(out)
}

diag_summary <- function(x) {
  tib <- summarise_draws(x, rhat, ess_bulk, ess_tail)
  finite_rhat <- tib$rhat[is.finite(tib$rhat)]
  finite_bulk <- tib$ess_bulk[is.finite(tib$ess_bulk)]
  finite_tail <- tib$ess_tail[is.finite(tib$ess_tail)]
  c(
    n = nrow(tib),
    rhat_max = max(finite_rhat),
    rhat_q99 = unname(quantile(finite_rhat, 0.99)),
    rhat_gt_1_01 = sum(finite_rhat > 1.01),
    ess_bulk_min = min(finite_bulk),
    ess_bulk_q10 = unname(quantile(finite_bulk, 0.10)),
    ess_bulk_median = median(finite_bulk),
    ess_tail_min = min(finite_tail),
    ess_tail_q10 = unname(quantile(finite_tail, 0.10)),
    ess_tail_median = median(finite_tail)
  )
}

canonical_components <- function(theta, N, J) {
  ndraw <- nrow(theta)
  grand <- numeric(ndraw)
  main1 <- matrix(NA_real_, ndraw, N)
  main2 <- matrix(NA_real_, ndraw, J)
  interaction <- matrix(NA_real_, ndraw, N * J)
  for (s in seq_len(ndraw)) {
    m <- matrix(theta[s, ], nrow = N, ncol = J)
    grand[s] <- mean(m)
    r <- rowMeans(m) - grand[s]
    c <- colMeans(m) - grand[s]
    h <- sweep(sweep(m, 1L, r, "-"), 2L, c, "-") - grand[s]
    main1[s, ] <- r
    main2[s, ] <- c
    interaction[s, ] <- as.vector(h)
  }
  list(
    intercept = matrix(grand, ncol = 1L),
    main1 = main1,
    main2 = main2,
    interaction = interaction
  )
}

build_repo_theta <- function(fit, dat) {
  a <- draw_matrix(fit, "group_1_betas")
  b <- draw_matrix(fit, "group_2_betas")
  g <- draw_matrix(fit, "gammas")
  d <- draw_matrix(fit, "deltas")
  theta <- matrix(NA_real_, nrow(a), dat$N * dat$J)
  for (s in seq_len(nrow(a))) {
    gm <- matrix(g[s, ], nrow = dat$N, ncol = dat$K)
    dm <- matrix(d[s, ], nrow = dat$J, ncol = dat$K)
    m <- tcrossprod(gm, dm)
    m <- sweep(m, 1L, a[s, ], "+")
    m <- sweep(m, 2L, b[s, ], "+")
    theta[s, ] <- as.vector(m)
  }
  theta
}

build_semi_theta <- function(fit, dat) {
  alpha <- draw_matrix(fit, "Intercept")[, 1L]
  u <- draw_matrix(fit, "Qfm_1_1")
  v <- draw_matrix(fit, "Qfm_1_2")
  singular <- draw_matrix(fit, "fm_singular_1")
  sdfm <- draw_matrix(fit, "sdfm_1")[, 1L]
  z1 <- draw_matrix(fit, "zfm_main_1_1")
  z2 <- draw_matrix(fit, "zfm_main_1_2")
  c1 <- sqrt(dat$N / (dat$N - 1))
  c2 <- sqrt(dat$J / (dat$J - 1))
  cint <- sqrt(dat$N * dat$J)
  theta <- matrix(NA_real_, length(alpha), dat$N * dat$J)
  for (s in seq_along(alpha)) {
    um <- matrix(u[s, ], nrow = dat$N, ncol = dat$K)
    vm <- matrix(v[s, ], nrow = dat$J, ncol = dat$K)
    m <- tcrossprod(sweep(um, 2L, singular[s, ], "*"), vm)
    m <- sdfm[s] * cint * m
    m <- sweep(m, 1L, 3 * c1 * z1[s, ], "+")
    m <- sweep(m, 2L, 3 * c2 * z2[s, ], "+")
    theta[s, ] <- as.vector(m + alpha[s])
  }
  theta
}

raw_diagnostics <- function(fit, kind) {
  variables <- if (kind == "repo") {
    c("gammas", "deltas")
  } else {
    c("zfm_frame_1_1", "zfm_frame_1_2")
  }
  x <- fit$summary(variables = variables)
  rhat <- x$rhat[is.finite(x$rhat)]
  bulk <- x$ess_bulk[is.finite(x$ess_bulk)]
  tail <- x$ess_tail[is.finite(x$ess_tail)]
  c(
    n = nrow(x),
    rhat_max = max(rhat),
    rhat_q99 = unname(quantile(rhat, 0.99)),
    rhat_gt_1_01 = sum(rhat > 1.01),
    ess_bulk_min = min(bulk),
    ess_bulk_median = median(bulk),
    ess_tail_min = min(tail),
    ess_tail_median = median(tail)
  )
}

sampler_summary <- function(fit) {
  basic <- fit$diagnostic_summary()
  diagnostics <- as.matrix(fit$sampler_diagnostics(format = "draws_matrix"))
  time <- fit$time()
  chain_time <- time$chains
  c(
    divergences = sum(basic$num_divergent),
    treedepth_hits = sum(basic$num_max_treedepth),
    ebfmi_min = min(basic$ebfmi),
    ebfmi_median = median(basic$ebfmi),
    accept_mean = mean(diagnostics[, "accept_stat__"]),
    leapfrog_median = median(diagnostics[, "n_leapfrog__"]),
    leapfrog_q90 = unname(quantile(diagnostics[, "n_leapfrog__"], 0.90)),
    leapfrog_q99 = unname(quantile(diagnostics[, "n_leapfrog__"], 0.99)),
    leapfrog_mean = mean(diagnostics[, "n_leapfrog__"]),
    warmup_max_chain_seconds = max(chain_time$warmup),
    sampling_max_chain_seconds = max(chain_time$sampling),
    total_max_chain_seconds = max(chain_time$total)
  )
}

truth_components <- function(dat) {
  canonical_components(matrix(dat$theta_true, nrow = 1L), dat$N, dat$J)
}

quality_summary <- function(theta, components, dat, layout, seed) {
  theta_obs <- theta[, dat$cell_index, drop = FALSE]
  theta_mean <- colMeans(theta)
  interaction_mean <- colMeans(components$interaction)
  truth <- truth_components(dat)
  theta_q <- colQuantiles(theta, probs = c(0.05, 0.95), drop = FALSE)
  interaction_q <- colQuantiles(
    components$interaction,
    probs = c(0.05, 0.95),
    drop = FALSE
  )

  set.seed(seed)
  yrep <- theta_obs + matrix(rnorm(length(theta_obs)), nrow = nrow(theta_obs))
  pred_q <- colQuantiles(yrep, probs = c(0.05, 0.95), drop = FALSE)
  ppc_mse <- rowMeans((yrep - rep(dat$y, each = nrow(yrep)))^2)
  bayes_r2 <- apply(theta_obs, 1L, var)
  bayes_r2 <- bayes_r2 / (bayes_r2 + 1)

  log_lik <- -0.5 * (
    theta_obs - rep(dat$y, each = nrow(theta_obs))
  )^2 - 0.5 * log(2 * pi)
  r_eff <- relative_eff(exp(log_lik), chain_id = layout$chain)
  loo_fit <- loo(log_lik, r_eff = r_eff, cores = 4)

  metrics <- c(
    rmse_observed = sqrt(mean((colMeans(theta_obs) - dat$y)^2)),
    mae_observed = mean(abs(colMeans(theta_obs) - dat$y)),
    rmse_truth_theta = sqrt(mean((theta_mean - dat$theta_true)^2)),
    cor_truth_theta = cor(theta_mean, dat$theta_true),
    coverage90_truth_theta = mean(
      dat$theta_true >= theta_q[, 1L] & dat$theta_true <= theta_q[, 2L]
    ),
    mean_sd_theta = mean(colSds(theta)),
    mean_ci90_width_theta = mean(theta_q[, 2L] - theta_q[, 1L]),
    rmse_truth_interaction = sqrt(mean(
      (interaction_mean - as.numeric(truth$interaction))^2
    )),
    cor_truth_interaction = cor(
      interaction_mean,
      as.numeric(truth$interaction)
    ),
    coverage90_truth_interaction = mean(
      as.numeric(truth$interaction) >= interaction_q[, 1L] &
        as.numeric(truth$interaction) <= interaction_q[, 2L]
    ),
    mean_sd_interaction = mean(colSds(components$interaction)),
    predictive_coverage90_y = mean(
      dat$y >= pred_q[, 1L] & dat$y <= pred_q[, 2L]
    ),
    predictive_width90_y = mean(pred_q[, 2L] - pred_q[, 1L]),
    ppc_mse_median = median(ppc_mse),
    ppc_mse_q05 = unname(quantile(ppc_mse, 0.05)),
    ppc_mse_q95 = unname(quantile(ppc_mse, 0.95)),
    bayes_r2_median = median(bayes_r2),
    elpd_loo = loo_fit$estimates["elpd_loo", "Estimate"],
    elpd_loo_se = loo_fit$estimates["elpd_loo", "SE"],
    p_loo = loo_fit$estimates["p_loo", "Estimate"],
    looic = loo_fit$estimates["looic", "Estimate"],
    pareto_k_gt_0_7 = sum(pareto_k_values(loo_fit) > 0.7),
    pareto_k_max = max(pareto_k_values(loo_fit))
  )

  list(
    metrics = metrics,
    loo = loo_fit,
    theta_mean = theta_mean,
    theta_sd = colSds(theta),
    interaction_mean = interaction_mean,
    interaction_sd = colSds(components$interaction),
    intercept_mean = mean(components$intercept),
    intercept_sd = sd(components$intercept),
    main1_mean = colMeans(components$main1),
    main2_mean = colMeans(components$main2)
  )
}

geometry_summary <- function(theta, components, fit, kind, dat) {
  max_interaction_margin <- 0
  max_sixth_sv <- 0
  take <- seq_len(min(100L, nrow(theta)))
  for (s in take) {
    h <- matrix(components$interaction[s, ], dat$N, dat$J)
    max_interaction_margin <- max(
      max_interaction_margin,
      abs(rowSums(h)),
      abs(colSums(h))
    )
    sv <- svd(h, nu = 0, nv = 0)$d
    if (length(sv) > dat$K) {
      max_sixth_sv <- max(max_sixth_sv, sv[dat$K + 1L])
    }
  }

  out <- c(
    interaction_draws_checked = length(take),
    interaction_margin_max = max_interaction_margin,
    sixth_singular_value_max = max_sixth_sv
  )
  if (kind == "semi") {
    u <- draw_matrix(fit, "Qfm_1_1")
    v <- draw_matrix(fit, "Qfm_1_2")
    singular <- draw_matrix(fit, "fm_singular_1")
    sdfm <- draw_matrix(fit, "sdfm_1")[, 1L]
    z1 <- draw_matrix(fit, "zfm_main_1_1")
    z2 <- draw_matrix(fit, "zfm_main_1_2")
    matrix_margin <- 0
    orthogonality_error <- 0
    spectrum_order_error <- 0
    spectrum_norm_error <- 0
    rms_error <- 0
    for (s in seq_len(nrow(u))) {
      um <- matrix(u[s, ], dat$N, dat$K)
      vm <- matrix(v[s, ], dat$J, dat$K)
      matrix_margin <- max(
        matrix_margin,
        abs(colSums(um)), abs(colSums(vm))
      )
      orthogonality_error <- max(
        orthogonality_error,
        abs(crossprod(um) - diag(dat$K)),
        abs(crossprod(vm) - diag(dat$K))
      )
      spectrum_order_error <- max(
        spectrum_order_error,
        pmax(diff(singular[s, ]), 0)
      )
      spectrum_norm_error <- max(
        spectrum_norm_error,
        abs(sum(singular[s, ]^2) - 1)
      )
      h <- sqrt(dat$N * dat$J) * sdfm[s] *
        tcrossprod(sweep(um, 2L, singular[s, ], "*"), vm)
      rms_error <- max(rms_error, abs(sqrt(mean(h^2)) - sdfm[s]))
    }
    out <- c(
      out,
      frame_column_sum_max = matrix_margin,
      frame_orthogonality_error_max = orthogonality_error,
      spectrum_order_error_max = spectrum_order_error,
      spectrum_norm_error_max = spectrum_norm_error,
      interaction_rms_error_max = rms_error,
      native_vector_sum_max = max(abs(rowSums(z1)), abs(rowSums(z2)))
    )
  }
  out
}

analyze_fit <- function(dataset, kind) {
  message("Analyzing ", dataset, " / ", kind)
  dat <- datasets[[dataset]]
  fit <- read_csv_fit(kind, dataset)
  layout <- chain_layout(fit)
  theta <- if (kind == "repo") {
    build_repo_theta(fit, dat)
  } else {
    build_semi_theta(fit, dat)
  }
  components <- canonical_components(theta, dat$N, dat$J)

  invariant <- list(
    theta = diag_summary(as_derived_draws(theta, layout, "theta")),
    interaction = diag_summary(as_derived_draws(
      components$interaction,
      layout,
      "interaction"
    )),
    additive = diag_summary(as_derived_draws(
      cbind(components$intercept, components$main1, components$main2),
      layout,
      "additive"
    ))
  )
  if (kind == "semi") {
    invariant$spectrum <- diag_summary(fit$draws(
      variables = c("sdfm_1", "fm_singular_1"),
      format = "draws_array"
    ))
  }

  quality <- quality_summary(
    theta,
    components,
    dat,
    layout,
    seed = if (dataset == "simulation_1") 1001 else 1002
  )

  list(
    sampler = sampler_summary(fit),
    raw_factors = raw_diagnostics(fit, kind),
    invariants = invariant,
    geometry = geometry_summary(theta, components, fit, kind, dat),
    quality = quality
  )
}

results <- list()
for (dataset in names(datasets)) {
  results[[dataset]] <- list()
  for (kind in c("repo", "semi")) {
    results[[dataset]][[kind]] <- analyze_fit(dataset, kind)
    gc()
  }
}

comparisons <- list()
for (dataset in names(datasets)) {
  r <- results[[dataset]]$repo$quality
  s <- results[[dataset]]$semi$quality
  comparisons[[dataset]] <- c(
    theta_mean_rmse_between = sqrt(mean((s$theta_mean - r$theta_mean)^2)),
    theta_mean_cor_between = cor(s$theta_mean, r$theta_mean),
    interaction_mean_rmse_between = sqrt(mean(
      (s$interaction_mean - r$interaction_mean)^2
    )),
    interaction_mean_cor_between = cor(
      s$interaction_mean,
      r$interaction_mean
    ),
    mean_abs_theta_difference_over_quadrature_sd = mean(
      abs(s$theta_mean - r$theta_mean) /
        sqrt(s$theta_sd^2 + r$theta_sd^2)
    ),
    elpd_loo_semi_minus_repo =
      s$loo$estimates["elpd_loo", "Estimate"] -
      r$loo$estimates["elpd_loo", "Estimate"]
  )
}

saveRDS(
  list(results = results, comparisons = comparisons, datasets = datasets),
  file.path(root, "semi_diagnostic_results.rds")
)

flatten_section <- function(section) {
  rows <- list()
  for (dataset in names(results)) {
    for (kind in names(results[[dataset]])) {
      values <- results[[dataset]][[kind]][[section]]
      if (is.list(values) && section == "invariants") {
        for (component in names(values)) {
          rows[[length(rows) + 1L]] <- data.frame(
            dataset = dataset,
            model = kind,
            component = component,
            metric = names(values[[component]]),
            value = as.numeric(values[[component]])
          )
        }
      } else if (is.list(values) && section == "quality") {
        values <- values$metrics
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = dataset,
          model = kind,
          metric = names(values),
          value = as.numeric(values)
        )
      } else {
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = dataset,
          model = kind,
          metric = names(values),
          value = as.numeric(values)
        )
      }
    }
  }
  do.call(rbind, rows)
}

for (section in c("sampler", "raw_factors", "invariants", "geometry", "quality")) {
  write.csv(
    flatten_section(section),
    file.path(root, paste0("semi_", section, ".csv")),
    row.names = FALSE
  )
}

comparison_rows <- do.call(rbind, lapply(names(comparisons), function(dataset) {
  data.frame(
    dataset = dataset,
    metric = names(comparisons[[dataset]]),
    value = as.numeric(comparisons[[dataset]])
  )
}))
write.csv(
  comparison_rows,
  file.path(root, "semi_comparisons.csv"),
  row.names = FALSE
)

message("Analysis complete")
