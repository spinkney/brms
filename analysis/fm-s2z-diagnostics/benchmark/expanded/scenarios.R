# Scenario generation and prior calibration for the expanded FM benchmark.

centered_basis <- function(n) {
  basis <- stats::contr.helmert(n)
  sweep(basis, 2L, sqrt(colSums(basis^2)), "/")
}

centered_frame <- function(n, k) {
  basis <- centered_basis(n)
  z <- matrix(stats::rnorm((n - 1L) * k), nrow = n - 1L, ncol = k)
  basis %*% qr.Q(qr(z), complete = FALSE)
}

normalize_spectrum <- function(x) {
  x / sqrt(sum(x^2))
}

separated_spectrum <- function(k) {
  normalize_spectrum(exp(-0.12 * (seq_len(k) - 1L)))
}

is_connected_design <- function(i1, i2, n1, n2) {
  neighbors <- vector("list", n1 + n2)
  for (n in seq_along(i1)) {
    a <- i1[n]
    b <- n1 + i2[n]
    neighbors[[a]] <- c(neighbors[[a]], b)
    neighbors[[b]] <- c(neighbors[[b]], a)
  }
  seen <- rep(FALSE, n1 + n2)
  queue <- 1L
  seen[1L] <- TRUE
  while (length(queue)) {
    node <- queue[1L]
    queue <- queue[-1L]
    next_nodes <- unique(neighbors[[node]])
    next_nodes <- next_nodes[!seen[next_nodes]]
    if (length(next_nodes)) {
      seen[next_nodes] <- TRUE
      queue <- c(queue, next_nodes)
    }
  }
  all(seen)
}

sparse_cells <- function(n1, n2, nobs, zipf = FALSE) {
  stopifnot(nobs >= n1 + n2 - 1L, nobs < n1 * n2)
  cover <- unique(c(
    seq_len(n1) + (((seq_len(n1) - 1L) %% n2) * n1),
    1L + ((seq_len(n2) - 1L) * n1)
  ))
  remaining <- setdiff(seq_len(n1 * n2), cover)
  if (zipf) {
    row_weight <- sample(seq_len(n1)^-1.2)
    col_weight <- sample(seq_len(n2)^-0.9)
    i1 <- ((remaining - 1L) %% n1) + 1L
    i2 <- ((remaining - 1L) %/% n1) + 1L
    probability <- row_weight[i1] * col_weight[i2]
  } else {
    probability <- NULL
  }
  c(
    cover,
    sample(
      remaining, nobs - length(cover), replace = FALSE,
      prob = probability
    )
  )
}

calibrate_rms_prior <- function(n1, n2, k, factor_sd, seed,
                                draws = 40000L, batch = 2000L) {
  set.seed(seed)
  rms <- numeric(draws)
  c1 <- sqrt(n1 / (n1 - 1))
  c2 <- sqrt(n2 / (n2 - 1))
  position <- 1L
  while (position <= draws) {
    size <- min(batch, draws - position + 1L)
    left <- stats::rWishart(size, df = n1 - 1L, Sigma = diag(k))
    right <- stats::rWishart(size, df = n2 - 1L, Sigma = diag(k))
    values <- vapply(seq_len(size), function(i) {
      sqrt(sum(left[, , i] * right[, , i]) / (n1 * n2))
    }, numeric(1))
    take <- position:(position + size - 1L)
    rms[take] <- factor_sd^2 * c1 * c2 * values
    position <- position + size
  }
  variance_log <- log1p((stats::sd(rms) / mean(rms))^2)
  list(
    meanlog = log(mean(rms)) - variance_log / 2,
    sdlog = sqrt(variance_log),
    mean = mean(rms),
    sd = stats::sd(rms),
    quantiles = unname(stats::quantile(rms, c(0.05, 0.5, 0.95)))
  )
}

expanded_fm_specs <- function() {
  rank_specs <- lapply(c(1L, 2L, 5L, 8L), function(k) {
    list(
      id = paste0("rank_k", k), purpose = "rank_sweep",
      n1 = 30L, n2 = 20L, k = k, interaction_rms = 2.5,
      spectrum = separated_spectrum(k), design = "dense", reps = 1L
    )
  })
  c(
    rank_specs,
    list(
      list(
        id = "sparse_uniform", purpose = "matrix_completion",
        n1 = 60L, n2 = 40L, k = 5L, interaction_rms = 2,
        spectrum = separated_spectrum(5L), design = "sparse_uniform",
        nobs = 600L, reps = 1L
      ),
      list(
        id = "sparse_zipf", purpose = "degree_imbalance",
        n1 = 60L, n2 = 40L, k = 5L, interaction_rms = 2,
        spectrum = separated_spectrum(5L), design = "sparse_zipf",
        nobs = 600L, reps = 1L
      ),
      list(
        id = "spectral_tie", purpose = "tied_singular_subspace",
        n1 = 40L, n2 = 25L, k = 5L, interaction_rms = 2,
        spectrum = normalize_spectrum(c(1, 1, 0.65, 0.4, 0.25)),
        design = "dense", reps = 1L
      ),
      list(
        id = "weak_interaction", purpose = "weak_signal",
        n1 = 40L, n2 = 25L, k = 5L, interaction_rms = 0.25,
        spectrum = separated_spectrum(5L), design = "dense", reps = 1L
      ),
      list(
        id = "full_rank_fixed", purpose = "fixed_determinant_component",
        n1 = 6L, n2 = 6L, k = 5L, interaction_rms = 1.5,
        spectrum = normalize_spectrum(c(1, 0.82, 0.65, 0.48, 0.32)),
        design = "replicated", reps = 15L
      )
    )
  )
}

make_expanded_fm_scenario <- function(spec, index) {
  data_seed <- 7100L + 37L * index
  prior_seed <- 9100L + 41L * index
  sampler_seed <- 11100L + 43L * index
  set.seed(data_seed)

  n1 <- spec$n1
  n2 <- spec$n2
  k <- spec$k
  q1 <- centered_frame(n1, k)
  q2 <- centered_frame(n2, k)
  if (k == n1 - 1L && k == n2 - 1L) {
    basis1 <- centered_basis(n1)
    basis2 <- centered_basis(n2)
    orientation <- det(crossprod(basis1, q1)) *
      det(crossprod(basis2, q2))
    if (orientation < 0) {
      q2[, 1L] <- -q2[, 1L]
    }
  }
  singular <- normalize_spectrum(spec$spectrum)
  interaction <- spec$interaction_rms * sqrt(n1 * n2) *
    tcrossprod(sweep(q1, 2L, singular, "*"), q2)
  main1 <- stats::rnorm(n1)
  main1 <- main1 - mean(main1)
  main1 <- 0.8 * main1 / sqrt(mean(main1^2))
  main2 <- stats::rnorm(n2)
  main2 <- main2 - mean(main2)
  main2 <- 0.8 * main2 / sqrt(mean(main2^2))
  intercept <- 0.3
  theta <- as.vector(intercept + outer(main1, rep(1, n2)) +
                       outer(rep(1, n1), main2) + interaction)

  all_cells <- seq_len(n1 * n2)
  if (spec$design == "dense") {
    train_cells <- all_cells
  } else if (spec$design == "replicated") {
    train_cells <- rep(all_cells, times = spec$reps)
  } else {
    train_cells <- sparse_cells(
      n1, n2, spec$nobs, zipf = spec$design == "sparse_zipf"
    )
  }
  i1 <- ((train_cells - 1L) %% n1) + 1L
  i2 <- ((train_cells - 1L) %/% n1) + 1L
  test_cells <- setdiff(all_cells, unique(train_cells))
  y_sd <- 1
  y <- theta[train_cells] + stats::rnorm(length(train_cells), sd = y_sd)
  y_replicate <- theta + stats::rnorm(length(theta), sd = y_sd)

  factor_sd <- sqrt(spec$interaction_rms / sqrt(k))
  rms_prior <- calibrate_rms_prior(
    n1, n2, k, factor_sd, seed = prior_seed
  )
  c1 <- sqrt(n1 / (n1 - 1))
  c2 <- sqrt(n2 / (n2 - 1))
  # A square orthogonal factor has two disconnected determinant components.
  # Select SO(K) for every square factor.  In the doubly-square benchmark this
  # deliberately restricts H to its positive-determinant component; the truth
  # is generated in that component above.  Leaving one determinant free is a
  # separate failure mode, not something Euclidean HMC can traverse.
  special1 <- as.integer(k == n1 - 1L)
  special2 <- as.integer(k == n2 - 1L)
  m1 <- as.integer((n1 - 1L) * k - k * (k - 1L) / 2L - special1)
  m2 <- as.integer((n2 - 1L) * k - k * (k - 1L) / 2L - special2)
  common_data <- list(
    N1 = n1, N2 = n2, K = k, Nobs = length(y),
    I1 = i1, I2 = i2, y = y, y_sd = y_sd,
    intercept_sd = 1, main_sd = 1.5, C1 = c1, C2 = c2
  )

  stopifnot(
    max(abs(rowSums(interaction))) < 1e-10,
    max(abs(colSums(interaction))) < 1e-10,
    abs(sqrt(mean(interaction^2)) - spec$interaction_rms) < 1e-10,
    is_connected_design(i1, i2, n1, n2),
    length(intersect(unique(train_cells), test_cells)) == 0L
  )

  list(
    id = spec$id, purpose = spec$purpose, design = spec$design,
    n1 = n1, n2 = n2, k = k, train_cells = train_cells,
    test_cells = test_cells, y_replicate = y_replicate,
    truth = list(
      intercept = intercept, main1 = main1, main2 = main2,
      interaction = interaction, theta = theta,
      rms = spec$interaction_rms, singular = singular
    ),
    factor_sd = factor_sd, rms_prior = rms_prior,
    centered_raw_data = c(common_data, list(factor_sd = factor_sd)),
    semi_data = c(common_data, list(
      Cinteraction = sqrt(n1 * n2), special1 = special1,
      special2 = special2, M1 = m1, M2 = m2,
      rms_meanlog = rms_prior$meanlog, rms_sdlog = rms_prior$sdlog
    )),
    seeds = list(
      data = data_seed, prior = prior_seed, sampler = sampler_seed
    )
  )
}

make_expanded_fm_scenarios <- function() {
  specs <- expanded_fm_specs()
  out <- lapply(seq_along(specs), function(i) {
    make_expanded_fm_scenario(specs[[i]], i)
  })
  names(out) <- vapply(out, `[[`, character(1), "id")
  out
}
