library(cmdstanr)
library(matrixStats)
library(posterior)

script_dir <- file.path(
  "analysis", "fm-s2z-diagnostics", "benchmark", "expanded"
)
source(file.path(script_dir, "scenarios.R"))

root <- Sys.getenv("FM_EXPANDED_DIR")
if (!nzchar(root)) {
  stop("Set FM_EXPANDED_DIR to the completed benchmark directory.")
}
scenario_file <- file.path(root, "scenarios.rds")
if (!file.exists(scenario_file)) {
  stop("Cannot find scenarios.rds; run expanded/run.R first.")
}
scenarios <- readRDS(scenario_file)
selected <- Sys.getenv("FM_EXPANDED_SCENARIOS")
if (nzchar(selected)) {
  selected <- trimws(strsplit(selected, ",", fixed = TRUE)[[1L]])
  unknown <- setdiff(selected, names(scenarios))
  if (length(unknown)) {
    stop("Unknown scenarios: ", paste(unknown, collapse = ", "))
  }
  scenarios <- scenarios[selected]
}
results_dir <- file.path(root, "results")
dir.create(results_dir, showWarnings = FALSE)

read_fit <- function(scenario, model) {
  run_dir <- file.path(root, "runs", scenario, model)
  manifest <- readRDS(file.path(run_dir, "manifest.rds"))
  csv_files <- file.path(run_dir, manifest$csv_files)
  if (!all(file.exists(csv_files))) {
    stop("Missing CmdStan CSV files for ", scenario, " / ", model)
  }
  list(fit = as_cmdstan_fit(csv_files), manifest = manifest)
}

draw_matrix <- function(fit, variable) {
  out <- as.matrix(fit$draws(variables = variable, format = "draws_matrix"))
  class(out) <- "matrix"
  out
}

chain_layout <- function(fit) {
  x <- as.data.frame(fit$draws(variables = "lp__", format = "draws_df"))
  list(chain = x$.chain, iteration = x$.iteration)
}

as_derived_draws <- function(x, layout, prefix) {
  x <- as.matrix(x)
  chains <- sort(unique(layout$chain))
  counts <- table(layout$chain)
  stopifnot(length(unique(as.integer(counts))) == 1L)
  niter <- as.integer(counts[1L])
  out <- array(
    NA_real_, dim = c(niter, length(chains), ncol(x)),
    dimnames = list(
      iteration = NULL, chain = as.character(chains),
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
  summary <- summarise_draws(x, rhat, ess_bulk, ess_tail)
  rhat_values <- summary$rhat[is.finite(summary$rhat)]
  bulk_values <- summary$ess_bulk[is.finite(summary$ess_bulk)]
  tail_values <- summary$ess_tail[is.finite(summary$ess_tail)]
  summarize_or_na <- function(values, fun) {
    if (length(values)) fun(values) else NA_real_
  }
  c(
    n = nrow(summary),
    rhat_max = summarize_or_na(rhat_values, max),
    rhat_gt_1_01 = sum(rhat_values > 1.01),
    ess_bulk_min = summarize_or_na(bulk_values, min),
    ess_bulk_median = summarize_or_na(bulk_values, median),
    ess_tail_min = summarize_or_na(tail_values, min),
    ess_tail_median = summarize_or_na(tail_values, median)
  )
}

anchor_frames <- function(left, right) {
  for (k in seq_len(ncol(left))) {
    pivot <- which.max(abs(left[, k]))
    if (left[pivot, k] < 0) {
      left[, k] <- -left[, k]
      right[, k] <- -right[, k]
    }
  }
  list(left = left, right = right)
}

build_surfaces <- function(fit, scenario, model) {
  n1 <- scenario$n1
  n2 <- scenario$n2
  k <- scenario$k
  intercept <- draw_matrix(fit, "Intercept")[, 1L]
  zmain1 <- draw_matrix(fit, "zmain1")
  zmain2 <- draw_matrix(fit, "zmain2")
  ndraw <- length(intercept)
  ncells <- n1 * n2
  interaction <- matrix(NA_real_, ndraw, ncells)
  theta <- matrix(NA_real_, ndraw, ncells)
  rms <- numeric(ndraw)
  spectrum <- matrix(NA_real_, ndraw, k)
  frames <- matrix(NA_real_, ndraw, (n1 + n2) * k)
  anchored_frames <- matrix(NA_real_, ndraw, (n1 + n2) * k)
  component1 <- matrix(NA_real_, ndraw, ncells)
  tie_block <- if (scenario$id == "spectral_tie") {
    matrix(NA_real_, ndraw, ncells)
  } else {
    NULL
  }
  left_projector1 <- right_projector1 <- NULL
  left_projector2 <- right_projector2 <- NULL
  if (scenario$id == "spectral_tie") {
    left_projector1 <- matrix(NA_real_, ndraw, n1 * n1)
    right_projector1 <- matrix(NA_real_, ndraw, n2 * n2)
    left_projector2 <- matrix(NA_real_, ndraw, n1 * n1)
    right_projector2 <- matrix(NA_real_, ndraw, n2 * n2)
  }
  determinant <- rep(NA_real_, ndraw)

  if (model == "centered_raw") {
    left_draws <- draw_matrix(fit, "left")
    right_draws <- draw_matrix(fit, "right")
  } else {
    q1_draws <- draw_matrix(fit, "Q1")
    q2_draws <- draw_matrix(fit, "Q2")
    singular_draws <- draw_matrix(fit, "singular")
    sdfm <- draw_matrix(fit, "sdfm")[, 1L]
  }

  c1 <- scenario$semi_data$C1
  c2 <- scenario$semi_data$C2
  main_sd <- scenario$semi_data$main_sd
  basis1 <- centered_basis(n1)
  basis2 <- centered_basis(n2)
  for (s in seq_len(ndraw)) {
    if (model == "centered_raw") {
      # Stan flattens array[K] vector[N] with the array index varying fastest:
      # left[1,1], ..., left[K,1], left[1,2], ... .  Rebuild K x N first.
      left <- t(matrix(left_draws[s, ], nrow = k, ncol = n1))
      right <- t(matrix(right_draws[s, ], nrow = k, ncol = n2))
      h <- scenario$factor_sd^2 * c1 * c2 * tcrossprod(left, right)
      decomposition <- svd(h, nu = k, nv = k)
      q1 <- decomposition$u[, seq_len(k), drop = FALSE]
      q2 <- decomposition$v[, seq_len(k), drop = FALSE]
      singular_values <- decomposition$d[seq_len(k)]
    } else {
      q1 <- matrix(q1_draws[s, ], nrow = n1, ncol = k)
      q2 <- matrix(q2_draws[s, ], nrow = n2, ncol = k)
      normalized <- singular_draws[s, ]
      singular_values <- sdfm[s] * sqrt(n1 * n2) * normalized
      h <- tcrossprod(sweep(q1, 2L, singular_values, "*"), q2)
    }
    normalized <- singular_values / sqrt(sum(singular_values^2))
    aligned <- anchor_frames(q1, q2)
    main1 <- main_sd * c1 * zmain1[s, ]
    main2 <- main_sd * c2 * zmain2[s, ]
    mean_surface <- h + outer(main1, rep(1, n2)) +
      outer(rep(1, n1), main2) + intercept[s]

    interaction[s, ] <- as.vector(h)
    theta[s, ] <- as.vector(mean_surface)
    rms[s] <- sqrt(mean(h^2))
    spectrum[s, ] <- normalized
    frames[s, ] <- c(as.vector(q1), as.vector(q2))
    anchored_frames[s, ] <- c(
      as.vector(aligned$left), as.vector(aligned$right)
    )
    component1[s, ] <- as.vector(
      singular_values[1L] * tcrossprod(q1[, 1L], q2[, 1L])
    )
    if (scenario$id == "spectral_tie") {
      tie_block[s, ] <- as.vector(tcrossprod(
        sweep(q1[, 1:2, drop = FALSE], 2L, singular_values[1:2], "*"),
        q2[, 1:2, drop = FALSE]
      ))
      left_projector1[s, ] <- as.vector(tcrossprod(q1[, 1L]))
      right_projector1[s, ] <- as.vector(tcrossprod(q2[, 1L]))
      left_projector2[s, ] <- as.vector(tcrossprod(q1[, 1:2, drop = FALSE]))
      right_projector2[s, ] <- as.vector(tcrossprod(q2[, 1:2, drop = FALSE]))
    }
    if (n1 - 1L == k && n2 - 1L == k) {
      small <- crossprod(basis1, h) %*% basis2
      determinant[s] <- det(small)
    }
  }

  list(
    theta = theta, interaction = interaction, rms = matrix(rms, ncol = 1L),
    spectrum = spectrum, frames = frames, anchored_frames = anchored_frames,
    component1 = component1, tie_block = tie_block,
    left_projector1 = left_projector1,
    right_projector1 = right_projector1,
    left_projector2 = left_projector2,
    right_projector2 = right_projector2,
    determinant = determinant
  )
}

semi_radii <- function(fit, scenario) {
  radii_one <- function(draws, nsmall, k, special) {
    maxcols <- k - as.integer(special == 1L && nsmall == k)
    if (maxcols == 0L) {
      return(matrix(numeric(0), nrow(draws), 0L))
    }
    lengths <- nsmall - seq_len(maxcols) + 1L
    starts <- cumsum(c(1L, head(lengths, -1L)))
    out <- matrix(NA_real_, nrow(draws), maxcols)
    for (j in seq_len(maxcols)) {
      take <- starts[j]:(starts[j] + lengths[j] - 1L)
      out[, j] <- sqrt(rowSums(draws[, take, drop = FALSE]^2))
    }
    out
  }
  z1 <- draw_matrix(fit, "zframe1")
  z2 <- draw_matrix(fit, "zframe2")
  cbind(
    radii_one(
      z1, scenario$n1 - 1L, scenario$k,
      scenario$semi_data$special1
    ),
    radii_one(
      z2, scenario$n2 - 1L, scenario$k,
      scenario$semi_data$special2
    )
  )
}

sampler_summary <- function(fit, manifest) {
  basic <- fit$diagnostic_summary()
  diagnostics <- as.matrix(fit$sampler_diagnostics(format = "draws_matrix"))
  chain_time <- manifest$time$chains
  c(
    divergences = sum(basic$num_divergent),
    treedepth_hits = sum(basic$num_max_treedepth),
    ebfmi_min = min(basic$ebfmi),
    step_size_mean = mean(diagnostics[, "stepsize__"]),
    leapfrog_mean = mean(diagnostics[, "n_leapfrog__"]),
    leapfrog_q90 = unname(quantile(diagnostics[, "n_leapfrog__"], 0.90)),
    leapfrog_q99 = unname(quantile(diagnostics[, "n_leapfrog__"], 0.99)),
    leapfrog_total = sum(diagnostics[, "n_leapfrog__"]),
    warmup_max_chain_seconds = max(chain_time$warmup),
    sampling_max_chain_seconds = max(chain_time$sampling),
    total_max_chain_seconds = max(chain_time$total)
  )
}

log_mean_exp <- function(x) {
  maximum <- max(x)
  maximum + log(mean(exp(x - maximum)))
}

quality_summary <- function(surface, scenario) {
  theta <- surface$theta
  interaction <- surface$interaction
  theta_mean <- colMeans(theta)
  interaction_mean <- colMeans(interaction)
  theta_quantiles <- colQuantiles(theta, probs = c(0.05, 0.95))
  spectrum_mean <- colMeans(surface$spectrum)
  replicate_lpd <- vapply(seq_along(scenario$truth$theta), function(cell) {
    log_mean_exp(dnorm(
      scenario$y_replicate[cell], theta[, cell],
      sd = scenario$semi_data$y_sd, log = TRUE
    ))
  }, numeric(1))
  train_cells <- scenario$train_cells
  train_predictions <- theta_mean[train_cells]
  observed <- scenario$centered_raw_data$y
  out <- c(
    rmse_observed = sqrt(mean((train_predictions - observed)^2)),
    rmse_truth_theta = sqrt(mean((theta_mean - scenario$truth$theta)^2)),
    cor_truth_theta = cor(theta_mean, scenario$truth$theta),
    coverage90_truth_theta = mean(
      scenario$truth$theta >= theta_quantiles[, 1L] &
        scenario$truth$theta <= theta_quantiles[, 2L]
    ),
    rmse_truth_interaction = sqrt(mean(
      (interaction_mean - as.vector(scenario$truth$interaction))^2
    )),
    cor_truth_interaction = cor(
      interaction_mean, as.vector(scenario$truth$interaction)
    ),
    rms_posterior_mean = mean(surface$rms),
    rms_posterior_sd = sd(surface$rms),
    spectrum_rmse = sqrt(mean(
      (spectrum_mean - scenario$truth$singular)^2
    )),
    replicate_elpd_per_cell = mean(replicate_lpd)
  )
  if (length(scenario$test_cells)) {
    heldout <- scenario$test_cells
    heldout_lpd <- replicate_lpd[heldout]
    out <- c(
      out,
      heldout_rmse_truth = sqrt(mean(
        (theta_mean[heldout] - scenario$truth$theta[heldout])^2
      )),
      heldout_elpd_per_cell = mean(heldout_lpd)
    )
    if (scenario$design == "sparse_zipf") {
      train_unique <- unique(scenario$train_cells)
      train_i1 <- ((train_unique - 1L) %% scenario$n1) + 1L
      train_i2 <- ((train_unique - 1L) %/% scenario$n1) + 1L
      degree1 <- tabulate(train_i1, nbins = scenario$n1)
      degree2 <- tabulate(train_i2, nbins = scenario$n2)
      heldout_i1 <- ((heldout - 1L) %% scenario$n1) + 1L
      heldout_i2 <- ((heldout - 1L) %/% scenario$n1) + 1L
      degree <- pmin(degree1[heldout_i1], degree2[heldout_i2])
      low <- degree <= median(degree)
      out <- c(
        out,
        heldout_rmse_low_degree = sqrt(mean(
          (theta_mean[heldout[low]] - scenario$truth$theta[heldout[low]])^2
        )),
        heldout_rmse_high_degree = sqrt(mean(
          (theta_mean[heldout[!low]] -
             scenario$truth$theta[heldout[!low]])^2
        )),
        heldout_elpd_low_degree = mean(heldout_lpd[low]),
        heldout_elpd_high_degree = mean(heldout_lpd[!low])
      )
    }
  }
  list(
    metrics = out, theta_mean = theta_mean, theta_sd = colSds(theta),
    interaction_mean = interaction_mean,
    interaction_sd = colSds(interaction)
  )
}

flatten_metrics <- function(values, scenario, model, component = NULL) {
  out <- data.frame(
    scenario = scenario, model = model,
    metric = names(values), value = as.numeric(values)
  )
  if (!is.null(component)) {
    out$component <- component
    out <- out[, c("scenario", "model", "component", "metric", "value")]
  }
  out
}

sampler_rows <- diagnostic_rows <- quality_rows <- raw_rows <- list()
component_rows <- comparison_rows <- count_rows <- list()

for (scenario_id in names(scenarios)) {
  scenario <- scenarios[[scenario_id]]
  message("Analyzing ", scenario_id)
  scenario_results <- list()
  for (model in c("centered_raw", "semiorthogonal")) {
    loaded <- read_fit(scenario_id, model)
    fit <- loaded$fit
    layout <- chain_layout(fit)
    surface <- build_surfaces(fit, scenario, model)
    sampler <- sampler_summary(fit, loaded$manifest)
    diagnostics <- list(
      theta = diag_summary(as_derived_draws(
        surface$theta, layout, "theta"
      )),
      interaction = diag_summary(as_derived_draws(
        surface$interaction, layout, "interaction"
      )),
      rms_spectrum = diag_summary(as_derived_draws(
        cbind(surface$rms, surface$spectrum), layout, "rms_spectrum"
      )),
      component1 = diag_summary(as_derived_draws(
        surface$component1, layout, "component1"
      )),
      anchored_frames = diag_summary(as_derived_draws(
        surface$anchored_frames, layout, "anchored_frames"
      ))
    )
    if (model == "semiorthogonal") {
      diagnostics$unanchored_frames <- diag_summary(as_derived_draws(
        surface$frames, layout, "unanchored_frames"
      ))
      raw <- cbind(draw_matrix(fit, "zframe1"), draw_matrix(fit, "zframe2"))
      raw_rows[[length(raw_rows) + 1L]] <- flatten_metrics(
        diag_summary(as_derived_draws(raw, layout, "raw_reflector")),
        scenario_id, model, "raw_reflector"
      )
      radii <- semi_radii(fit, scenario)
      raw_rows[[length(raw_rows) + 1L]] <- flatten_metrics(
        diag_summary(as_derived_draws(radii, layout, "reflector_radius")),
        scenario_id, model, "reflector_radius"
      )
    } else {
      raw <- cbind(draw_matrix(fit, "left"), draw_matrix(fit, "right"))
      raw_rows[[length(raw_rows) + 1L]] <- flatten_metrics(
        diag_summary(as_derived_draws(raw, layout, "raw_factor")),
        scenario_id, model, "raw_factor"
      )
    }
    if (!is.null(surface$tie_block)) {
      diagnostics$tied_block <- diag_summary(as_derived_draws(
        surface$tie_block, layout, "tied_block"
      ))
      diagnostics$left_projector1 <- diag_summary(as_derived_draws(
        surface$left_projector1, layout, "left_projector1"
      ))
      diagnostics$right_projector1 <- diag_summary(as_derived_draws(
        surface$right_projector1, layout, "right_projector1"
      ))
      diagnostics$left_projector2 <- diag_summary(as_derived_draws(
        surface$left_projector2, layout, "left_projector2"
      ))
      diagnostics$right_projector2 <- diag_summary(as_derived_draws(
        surface$right_projector2, layout, "right_projector2"
      ))
    }
    for (component in names(diagnostics)) {
      values <- diagnostics[[component]]
      values <- c(
        values,
        ess_bulk_median_per_million_leapfrogs =
          unname(values["ess_bulk_median"] /
                   sampler["leapfrog_total"] * 1e6),
        ess_bulk_median_per_sampling_second =
          unname(values["ess_bulk_median"] /
                   sampler["sampling_max_chain_seconds"])
      )
      diagnostic_rows[[length(diagnostic_rows) + 1L]] <- flatten_metrics(
        values, scenario_id, model, component
      )
    }
    quality <- quality_summary(surface, scenario)
    sampler_rows[[length(sampler_rows) + 1L]] <- flatten_metrics(
      sampler, scenario_id, model
    )
    quality_rows[[length(quality_rows) + 1L]] <- flatten_metrics(
      quality$metrics, scenario_id, model
    )
    if (scenario_id == "full_rank") {
      sign_positive <- surface$determinant > 0
      for (chain in sort(unique(layout$chain))) {
        take <- layout$chain == chain
        lp <- draw_matrix(fit, "lp__")[take, 1L]
        component_rows[[length(component_rows) + 1L]] <- data.frame(
          scenario = scenario_id, model = model, chain = chain,
          positive_fraction = mean(sign_positive[take]),
          transitions = sum(diff(sign_positive[take]) != 0),
          mean_lp = mean(lp), mean_lp_positive = mean(lp[sign_positive[take]]),
          mean_lp_negative = mean(lp[!sign_positive[take]])
        )
      }
    }
    scenario_results[[model]] <- list(quality = quality)
    rm(fit, surface)
    gc()
  }

  raw_quality <- scenario_results$centered_raw$quality
  semi_quality <- scenario_results$semiorthogonal$quality
  comparisons <- c(
    theta_mean_rmse_between = sqrt(mean(
      (semi_quality$theta_mean - raw_quality$theta_mean)^2
    )),
    theta_mean_cor_between = cor(
      semi_quality$theta_mean, raw_quality$theta_mean
    ),
    interaction_mean_rmse_between = sqrt(mean(
      (semi_quality$interaction_mean - raw_quality$interaction_mean)^2
    )),
    interaction_mean_cor_between = cor(
      semi_quality$interaction_mean, raw_quality$interaction_mean
    ),
    mean_abs_theta_difference_over_quadrature_sd = mean(
      abs(semi_quality$theta_mean - raw_quality$theta_mean) /
        sqrt(semi_quality$theta_sd^2 + raw_quality$theta_sd^2)
    )
  )
  comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
    scenario = scenario_id, metric = names(comparisons),
    value = as.numeric(comparisons)
  )

  centered_interaction <- (scenario$n1 + scenario$n2 - 2L) * scenario$k
  semi_interaction <- scenario$semi_data$M1 + scenario$semi_data$M2 +
    scenario$k
  identifiable_interaction <- scenario$k * (
    scenario$n1 + scenario$n2 - 2L - scenario$k
  )
  additive <- scenario$n1 + scenario$n2 - 1L
  count_rows[[length(count_rows) + 1L]] <- rbind(
    data.frame(
      scenario = scenario_id, model = "centered_raw",
      total = additive + centered_interaction,
      interaction = centered_interaction,
      interaction_excess = centered_interaction - identifiable_interaction
    ),
    data.frame(
      scenario = scenario_id, model = "semiorthogonal",
      total = additive + semi_interaction,
      interaction = semi_interaction,
      interaction_excess = semi_interaction - identifiable_interaction
    )
  )
}

write.csv(do.call(rbind, sampler_rows), file.path(results_dir, "sampler.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, diagnostic_rows),
          file.path(results_dir, "diagnostics.csv"), row.names = FALSE)
write.csv(do.call(rbind, quality_rows), file.path(results_dir, "quality.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, raw_rows), file.path(results_dir, "raw.csv"),
          row.names = FALSE)
write.csv(do.call(rbind, comparison_rows),
          file.path(results_dir, "comparisons.csv"), row.names = FALSE)
write.csv(do.call(rbind, count_rows),
          file.path(results_dir, "parameter_counts.csv"), row.names = FALSE)
if (length(component_rows)) {
  write.csv(do.call(rbind, component_rows),
            file.path(results_dir, "determinant_components.csv"),
            row.names = FALSE)
}

message("Expanded FM analysis complete.")
