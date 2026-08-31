context("Centering transformations for physical S2Z effects")

# R mirror of inst/chunks/fun_sum_to_zero.stan. Keeping this small reference
# implementation in the tests makes the geometry checks independent of Stan.
.s2z_constrain_r <- function(y) {
  n <- length(y)
  z <- numeric(n + 1L)
  sum_w <- 0
  for (ii in seq_len(n)) {
    i <- n - ii + 1L
    w <- y[i] / sqrt(i * (i + 1))
    sum_w <- sum_w + w
    z[i] <- z[i] + sum_w
    z[i + 1L] <- z[i + 1L] - i * w
  }
  z
}

.s2z_basis <- function(n) {
  stopifnot(n >= 2L)
  vapply(
    seq_len(n - 1L),
    function(j) .s2z_constrain_r(diag(n - 1L)[, j]),
    numeric(n)
  )
}

# Compute a dense log absolute determinant after row and column equilibration.
# This keeps the numerical check useful when the map contains scales many
# orders of magnitude apart. QR is the primary result and singular values are
# retained as an independent stability check.
.s2z_balanced_logabsdet <- function(x) {
  stopifnot(is.matrix(x), nrow(x) == ncol(x), all(is.finite(x)))
  row_scale <- apply(abs(x), 1L, max)
  stopifnot(all(is.finite(row_scale)), all(row_scale > 0))
  balanced <- sweep(x, 1L, row_scale, "/")
  col_scale <- apply(abs(balanced), 2L, max)
  stopifnot(all(is.finite(col_scale)), all(col_scale > 0))
  balanced <- sweep(balanced, 2L, col_scale, "/")
  adjustment <- sum(log(row_scale)) + sum(log(col_scale))

  R <- qr.R(qr(balanced, LAPACK = TRUE))
  qr_value <- sum(log(abs(diag(R)))) + adjustment
  singular_values <- svd(balanced, nu = 0L, nv = 0L)$d
  svd_value <- sum(log(singular_values)) + adjustment
  list(
    qr = qr_value,
    svd = svd_value,
    balanced_condition = max(singular_values) / min(singular_values)
  )
}

# Independent reference for the direct non-centered S2Z coordinates.  The
# centered coordinates are physical contrasts X = Z L', while the direct
# non-centered form samples Z and applies L after the orthonormal group
# rotation.  This helper deliberately constructs the full Jacobian rather
# than assuming its determinant.
.s2z_noncentered_change_case <- function(scales) {
  n <- 7L
  k <- length(scales)
  basis <- .s2z_basis(n)
  z <- matrix(
    sin(seq_len((n - 1L) * k) * 0.41) +
      cos(seq_len((n - 1L) * k) * 0.73),
    nrow = n - 1L, ncol = k
  ) / 2.3
  correlation <- outer(seq_len(k), seq_len(k), function(i, j) {
    0.23^abs(i - j)
  })
  L_cor <- t(chol(correlation))
  L_base <- diag(scales, nrow = k) %*% L_cor

  centered_coordinates <- z %*% t(L_base)
  centered_effects <- basis %*% centered_coordinates
  noncentered_effects <- (basis %*% z) %*% t(L_base)
  whitened_effects <- t(forwardsolve(L_base, t(centered_effects)))

  log_centered <- -0.5 * sum(whitened_effects^2) -
    (n - 1L) * sum(log(diag(L_base))) -
    0.5 * (n - 1L) * k * log(2 * pi)
  log_noncentered <- sum(stats::dnorm(z, log = TRUE))
  coordinate_map <- kronecker(L_base, diag(n - 1L))
  log_jacobian <- as.numeric(
    determinant(coordinate_map, logarithm = TRUE)$modulus
  )

  # Use the same columns as population and group-level predictors.  Shifting
  # a conventional common group mean into the population coefficients must
  # leave the likelihood-scale predictor unchanged in either S2Z form.
  group <- rep(seq_len(n), each = 3L)
  design <- matrix(
    sin(seq_len(length(group) * k) * 0.29) -
      cos(seq_len(length(group) * k) * 0.17),
    nrow = length(group), ncol = k
  )
  group_mean <- seq(-0.45, 0.65, length.out = k)
  population <- seq(0.8, -0.55, length.out = k)
  conventional_effects <- sweep(
    centered_effects, 2L, group_mean, "+"
  )
  finite_population <- population + group_mean
  eta_conventional <- drop(design %*% population) + rowSums(
    design * conventional_effects[group, , drop = FALSE]
  )
  eta_centered <- drop(design %*% finite_population) + rowSums(
    design * centered_effects[group, , drop = FALSE]
  )
  eta_noncentered <- drop(design %*% finite_population) + rowSums(
    design * noncentered_effects[group, , drop = FALSE]
  )

  list(
    n = n, k = k, L_base = L_base,
    centered_effects = centered_effects,
    noncentered_effects = noncentered_effects,
    log_centered = log_centered,
    log_noncentered = log_noncentered,
    log_jacobian = log_jacobian,
    eta_conventional = eta_conventional,
    eta_centered = eta_centered,
    eta_noncentered = eta_noncentered
  )
}

# Independent reference for the group-specific partially centered change of
# variables.  Each row uses D_j = diag(rho_j) L + diag(1 - rho_j), followed
# by A_j = L D_j^-1, and is projected back onto the physical zero-sum
# subspace.  The dense restricted Jacobian is constructed explicitly rather
# than using the determinant formula exercised by the generated Stan code.
.s2z_partial_change_case <- function(L, rho) {
  n <- nrow(rho)
  k <- ncol(rho)
  stopifnot(
    n >= 2L, nrow(L) == k, ncol(L) == k,
    all(rho >= 0), all(rho <= 1)
  )
  basis <- .s2z_basis(n)
  z <- matrix(
    sin(seq_len((n - 1L) * k) * 0.53) -
      cos(seq_len((n - 1L) * k) * 0.31),
    nrow = n - 1L, ncol = k
  ) / 2.1

  transform <- function(z) {
    raw <- basis %*% z
    transformed <- matrix(NA_real_, nrow = n, ncol = k)
    for (j in seq_len(n)) {
      D_j <- diag(rho[j, ], nrow = k) %*% L +
        diag(1 - rho[j, ], nrow = k)
      transformed[j, ] <- drop(
        L %*% forwardsolve(D_j, raw[j, ], upper.tri = FALSE)
      )
    }
    sweep(transformed, 2L, colMeans(transformed), "-")
  }

  effects <- transform(z)
  dimensions <- (n - 1L) * k
  restricted_map <- vapply(
    seq_len(dimensions),
    function(i) {
      unit <- numeric(dimensions)
      unit[i] <- 1
      as.vector(crossprod(
        basis,
        transform(matrix(unit, nrow = n - 1L, ncol = k))
      ))
    },
    numeric(dimensions)
  )
  dense_log_jacobian <- .s2z_balanced_logabsdet(restricted_map)
  numeric_log_jacobian <- dense_log_jacobian$qr

  D <- lapply(seq_len(n), function(j) {
    diag(rho[j, ], nrow = k) %*% L +
      diag(1 - rho[j, ], nrow = k)
  })
  D_bar <- Reduce("+", D) / n
  formula_log_jacobian <- (n - 1L) * sum(log(diag(L))) -
    sum(vapply(D, function(x) sum(log(diag(x))), numeric(1))) +
    sum(log(diag(D_bar)))
  determinant_correction <-
    -sum(vapply(D, function(x) sum(log(diag(x))), numeric(1))) +
    sum(log(diag(D_bar)))

  whitened <- t(forwardsolve(L, t(effects)))
  physical_log_kernel <- -0.5 * sum(whitened^2) -
    (n - 1L) * sum(log(diag(L)))

  list(
    n = n, k = k, basis = basis, z = z, rho = rho, L = L,
    effects = effects, restricted_map = restricted_map,
    numeric_log_jacobian = numeric_log_jacobian,
    numeric_log_jacobian_svd = dense_log_jacobian$svd,
    balanced_condition = dense_log_jacobian$balanced_condition,
    formula_log_jacobian = formula_log_jacobian,
    determinant_correction = determinant_correction,
    physical_log_kernel = physical_log_kernel,
    transformed_log_kernel = physical_log_kernel + numeric_log_jacobian,
    expected_transformed_log_kernel =
      -0.5 * sum(whitened^2) + determinant_correction
  )
}

# Full explicit-mean reference for an S2Z chart with a non-Gaussian
# population prior. Besides the restricted contrast map, this includes the
# standardized omitted mean, the unit-determinant population-coordinate
# shear, and optional conditional Student-t group scales.
.s2z_explicit_logistic_change_case <- function(L, rho, mixers = NULL) {
  n <- nrow(rho)
  k <- ncol(rho)
  stopifnot(
    n >= 2L, nrow(L) == k, ncol(L) == k,
    all(rho >= 0), all(rho <= 1)
  )
  if (is.null(mixers)) {
    mixers <- rep(1, n)
  }
  stopifnot(length(mixers) == n, all(is.finite(mixers)), all(mixers > 0))
  basis <- .s2z_basis(n)
  contrast_dim <- (n - 1L) * k

  transform <- function(x) {
    z <- matrix(x[seq_len(contrast_dim)], nrow = n - 1L, ncol = k)
    v <- x[contrast_dim + seq_len(k)]
    theta <- x[contrast_dim + k + seq_len(k)]
    raw <- basis %*% z
    delta <- matrix(NA_real_, nrow = n, ncol = k)
    for (j in seq_len(n)) {
      D_j <- diag(rho[j, ], nrow = k) %*% L +
        diag(1 - rho[j, ], nrow = k)
      delta[j, ] <- drop(
        L %*% forwardsolve(D_j, raw[j, ], upper.tri = FALSE)
      )
    }
    delta <- sweep(delta, 2L, colMeans(delta), "-")
    mean_effect <- drop(L %*% v / sqrt(n))
    group_effect <- sweep(delta, 2L, mean_effect, "+")
    population <- theta - mean_effect
    c(as.vector(group_effect), population)
  }

  input_dim <- contrast_dim + 2L * k
  x <- c(
    sin(seq_len(contrast_dim) * 0.37) / 1.9,
    seq(-0.45, 0.35, length.out = k),
    seq(0.8, -0.25, length.out = k)
  )
  jacobian <- vapply(seq_len(input_dim), function(i) {
    unit <- numeric(input_dim)
    unit[i] <- 1
    transform(unit)
  }, numeric(input_dim))
  numeric_log_jacobian <- .s2z_balanced_logabsdet(jacobian)$qr

  D <- lapply(seq_len(n), function(j) {
    diag(rho[j, ], nrow = k) %*% L +
      diag(1 - rho[j, ], nrow = k)
  })
  D_bar <- Reduce("+", D) / n
  determinant_correction <-
    -sum(vapply(D, function(x) sum(log(diag(x))), numeric(1))) +
    sum(log(diag(D_bar)))
  formula_log_jacobian <- n * sum(log(diag(L))) +
    determinant_correction

  transformed <- transform(x)
  group_effect <- matrix(
    transformed[seq_len(n * k)], nrow = n, ncol = k
  )
  population <- tail(transformed, k)
  whitened <- t(forwardsolve(L, t(group_effect))) / mixers
  physical_log_density <- sum(dnorm(whitened, log = TRUE)) -
    n * sum(log(diag(L))) - k * sum(log(mixers)) +
    sum(dlogis(population, location = -0.4, scale = 1.3, log = TRUE))
  chart_log_density <- sum(dnorm(whitened, log = TRUE)) -
    k * sum(log(mixers)) + determinant_correction +
    sum(dlogis(population, location = -0.4, scale = 1.3, log = TRUE))

  list(
    numeric_log_jacobian = numeric_log_jacobian,
    formula_log_jacobian = formula_log_jacobian,
    determinant_correction = determinant_correction,
    physical_plus_jacobian = physical_log_density + numeric_log_jacobian,
    chart_log_density = chart_log_density
  )
}

# Mirror the diagonal-only Fisher contraction emitted by the Stan generator.
# If C C' = I + L' J L and W = C^-1 L', then the desired posterior covariance
# is W' W, so its diagonal is available without forming either inverse.
.s2z_fisher_reliability_diag <- function(gram, L, obs_prec) {
  k <- nrow(L)
  stopifnot(
    ncol(L) == k, identical(dim(gram), c(k, k)),
    is.finite(obs_prec), obs_prec > 0
  )
  K <- obs_prec * crossprod(L, gram %*% L)
  K <- 0.5 * (K + t(K))
  C <- t(chol(diag(k) + K))
  W <- forwardsolve(C, t(L))
  pmin(1, pmax(0, 1 - colSums(W^2) / rowSums(L^2)))
}

.s2z_fisher_rho_dense <- function(Sigma, information) {
  posterior <- solve(solve(Sigma) + information)
  pmin(1, pmax(0, 1 - diag(posterior) / diag(Sigma)))
}

test_that("direct non-centered S2Z is an exact scaled change of variables", {
  scale_cases <- c(
    lapply(c(1e-5, 1e-3, 0.1, 1, 5), c),
    list(c(1e-5, 1e-3, 0.1, 1, 5))
  )
  for (scales in scale_cases) {
    ans <- .s2z_noncentered_change_case(scales)
    expected_log_jacobian <-
      (ans$n - 1L) * sum(log(diag(ans$L_base)))

    expect_equal(
      ans$noncentered_effects, ans$centered_effects,
      tolerance = 2e-13, scale = 1
    )
    expect_equal(
      ans$log_jacobian, expected_log_jacobian,
      tolerance = 2e-12, scale = 1
    )
    expect_equal(
      ans$log_centered + ans$log_jacobian,
      ans$log_noncentered, tolerance = 3e-11, scale = 1
    )
    expect_equal(
      -(ans$n - 1L) * sum(log(diag(ans$L_base))) +
        ans$log_jacobian,
      0, tolerance = 2e-12, scale = 1
    )
    expect_equal(
      ans$eta_centered, ans$eta_conventional,
      tolerance = 3e-13, scale = 1
    )
    expect_equal(
      ans$eta_noncentered, ans$eta_conventional,
      tolerance = 3e-13, scale = 1
    )
  }
})

test_that("scalar partial S2Z uses the exact restricted Jacobian", {
  n <- 7L
  rho_cases <- list(
    rep(1, n),
    rep(0, n),
    rep(0.37, n),
    c(0, 0.08, 0.21, 0.46, 0.7, 0.91, 1)
  )
  for (scale in c(0.03, 0.4, 1, 3.7)) {
    for (rho in rho_cases) {
      ans <- .s2z_partial_change_case(
        matrix(scale, nrow = 1L), matrix(rho, ncol = 1L)
      )
      d <- rho * scale + 1 - rho
      expected_correction <- -sum(log(d)) + log(mean(d))

      expect_equal(colSums(ans$effects), 0, tolerance = 3e-14)
      expect_equal(
        ans$numeric_log_jacobian, ans$formula_log_jacobian,
        tolerance = 2e-12, scale = 1
      )
      expect_equal(
        ans$determinant_correction, expected_correction,
        tolerance = 2e-14, scale = 1
      )
      expect_equal(
        ans$transformed_log_kernel,
        ans$expected_transformed_log_kernel,
        tolerance = 3e-12, scale = 1
      )
    }
  }

  centered <- .s2z_partial_change_case(
    matrix(2.3, nrow = 1L), matrix(1, nrow = n, ncol = 1L)
  )
  noncentered <- .s2z_partial_change_case(
    matrix(2.3, nrow = 1L), matrix(0, nrow = n, ncol = 1L)
  )
  expect_equal(
    centered$effects, centered$basis %*% centered$z,
    tolerance = 3e-14
  )
  expect_equal(
    noncentered$effects, 2.3 * noncentered$basis %*% noncentered$z,
    tolerance = 3e-14
  )
})

test_that("exact logistic S2Z means compose with every centering chart", {
  n <- 5L
  L_cor <- matrix(c(1.7, -0.45, 0, 0.65), 2L, 2L)
  cases <- list(
    correlated_noncentered = list(
      L = L_cor, rho = matrix(0, n, 2L), mixers = rep(1, n)
    ),
    correlated_partial = list(
      L = L_cor,
      rho = matrix(c(
        0, 0.15, 0.4, 0.75, 1,
        0.9, 0.65, 0.35, 0.1, 0
      ), n, 2L),
      mixers = rep(1, n)
    ),
    correlated_centered = list(
      L = L_cor, rho = matrix(1, n, 2L), mixers = rep(1, n)
    ),
    independent_partial = list(
      L = diag(c(0.35, 2.4)),
      rho = matrix(seq(0.05, 0.95, length.out = 2L * n), n, 2L),
      mixers = rep(1, n)
    ),
    scalar_student_partial = list(
      L = matrix(1.25, 1L, 1L),
      rho = matrix(c(0, 0.2, 0.55, 0.85, 1), n, 1L),
      mixers = c(0.45, 0.8, 1.1, 1.7, 2.3)
    )
  )

  for (name in names(cases)) {
    case <- cases[[name]]
    ans <- .s2z_explicit_logistic_change_case(
      case$L, case$rho, mixers = case$mixers
    )
    expect_equal(
      ans$numeric_log_jacobian, ans$formula_log_jacobian,
      tolerance = 2e-10, scale = 1, info = name
    )
    expect_equal(
      ans$physical_plus_jacobian, ans$chart_log_density,
      tolerance = 2e-10, scale = 1, info = name
    )
  }

  expect_equal(
    .s2z_explicit_logistic_change_case(
      L_cor, matrix(0, n, 2L)
    )$determinant_correction,
    0, tolerance = 2e-14
  )
  expect_equal(
    .s2z_explicit_logistic_change_case(
      L_cor, matrix(1, n, 2L)
    )$determinant_correction,
    -(n - 1L) * sum(log(diag(L_cor))), tolerance = 2e-14
  )
})

test_that("diagonal-only Fisher reliabilities equal the dense contraction", {
  set.seed(72841)
  for (k in c(1L, 2L, 4L)) {
    for (iteration in seq_len(20L)) {
      design <- matrix(rnorm((k + 3L) * k), ncol = k)
      gram <- crossprod(design)
      L <- matrix(0, k, k)
      L[lower.tri(L, diag = TRUE)] <- rnorm(k * (k + 1L) / 2L)
      diag(L) <- exp(rnorm(k, sd = 0.7))
      obs_prec <- exp(rnorm(1L, sd = 1.2))

      optimized <- .s2z_fisher_reliability_diag(gram, L, obs_prec)
      white_cov <- solve(diag(k) + obs_prec * crossprod(L, gram %*% L))
      posterior_cov <- L %*% white_cov %*% t(L)
      prior_cov <- L %*% t(L)
      dense <- 1 - diag(posterior_cov) / diag(prior_cov)

      expect_equal(optimized, dense, tolerance = 2e-12, scale = 1)
    }
  }
})

test_that("Fisher fractions respect rescaling, ordering, and absent slopes", {
  Sigma <- matrix(c(
    1.7, 0.35, -0.2,
    0.35, 0.9, 0.16,
    -0.2, 0.16, 1.25
  ), 3L, 3L)
  design <- matrix(c(
    1, -0.4, 0.2,
    1, 0.1, -0.3,
    1, 0.7, 0.8,
    1, 1.2, -0.1
  ), ncol = 3L, byrow = TRUE)
  information <- crossprod(design)
  reference <- .s2z_fisher_rho_dense(Sigma, information)

  rescale <- diag(c(3.5, 0.4, 1.8))
  scaled_sigma <- solve(rescale) %*% Sigma %*% solve(rescale)
  scaled_information <- rescale %*% information %*% rescale
  expect_equal(
    .s2z_fisher_rho_dense(scaled_sigma, scaled_information),
    reference,
    tolerance = 2e-14
  )

  permutation <- c(3L, 1L, 2L)
  expect_equal(
    .s2z_fisher_rho_dense(
      Sigma[permutation, permutation],
      information[permutation, permutation]
    ),
    reference[permutation],
    tolerance = 2e-14
  )

  independent_sigma <- diag(c(1.2, 0.7, 2.1))
  absent_slope_information <- diag(c(4.5, 0, 1.3))
  absent <- .s2z_fisher_rho_dense(
    independent_sigma, absent_slope_information
  )
  expect_equal(absent[2], 0, tolerance = 1e-15)
  expect_equal(
    .s2z_fisher_rho_dense(independent_sigma, matrix(0, 3L, 3L)),
    numeric(3L),
    tolerance = 1e-15
  )
})

test_that("independent Fisher reliabilities cancel diagonal prior scales", {
  set.seed(4147)
  for (k in c(2L, 4L)) {
    for (iteration in seq_len(20L)) {
      design <- matrix(rnorm((k + 3L) * k), ncol = k)
      gram <- crossprod(design)
      scale <- exp(rnorm(k, sd = 0.8))
      row_var <- exp(rnorm(1L, sd = 0.6))
      obs_prec <- exp(rnorm(1L, sd = 1.1))
      K <- row_var * obs_prec * (scale * gram * rep(scale, each = k))
      C <- t(chol(diag(k) + 0.5 * (K + t(K))))
      relative_post_var <- numeric(k)
      for (j in seq_len(k)) {
        rhs <- numeric(k)
        rhs[j] <- 1
        relative_post_var[j] <- sum(forwardsolve(C, rhs)^2)
      }
      optimized <- 1 - relative_post_var
      dense <- .s2z_fisher_reliability_diag(
        gram, diag(scale), row_var * obs_prec
      )

      expect_equal(optimized, dense, tolerance = 2e-12, scale = 1)
    }
  }
})

test_that("correlated partial S2Z uses the exact restricted Jacobian", {
  n <- 6L
  k <- 4L
  L <- matrix(
    c(
      0.55, 0, 0, 0,
      0.17, 1.2, 0, 0,
      -0.11, 0.28, 0.8, 0,
      0.21, -0.16, 0.19, 1.6
    ),
    nrow = k, byrow = TRUE
  )
  rho <- matrix(
    c(
      0, 0.15, 0.4, 1,
      0.08, 0.3, 0.65, 0.92,
      0.2, 0.45, 0.78, 0.83,
      0.37, 0.62, 0.91, 0.71,
      0.73, 0.86, 0.22, 0.58,
      1, 0.04, 0.53, 0
    ),
    nrow = n, byrow = TRUE
  )
  partial <- .s2z_partial_change_case(L, rho)

  expect_equal(colSums(partial$effects), numeric(k), tolerance = 4e-14)
  expect_equal(
    partial$numeric_log_jacobian, partial$formula_log_jacobian,
    tolerance = 2e-11, scale = 1
  )
  expect_equal(
    partial$transformed_log_kernel,
    partial$expected_transformed_log_kernel,
    tolerance = 2e-11, scale = 1
  )

  centered <- .s2z_partial_change_case(
    L, matrix(1, nrow = n, ncol = k)
  )
  noncentered <- .s2z_partial_change_case(
    L, matrix(0, nrow = n, ncol = k)
  )
  expect_equal(
    centered$effects, centered$basis %*% centered$z,
    tolerance = 3e-14
  )
  expect_equal(
    noncentered$effects,
    (noncentered$basis %*% noncentered$z) %*% t(L),
    tolerance = 3e-14
  )
  expect_equal(
    centered$determinant_correction,
    -(n - 1L) * sum(log(diag(L))), tolerance = 2e-14
  )
  expect_equal(noncentered$determinant_correction, 0, tolerance = 2e-14)
})

test_that("correlated partial Jacobian is stable under extreme Cholesky scales", {
  n <- 6L
  k <- 4L
  # Any lower-triangular matrix with positive diagonal is a valid covariance
  # Cholesky factor. Here the diagonal spans more than 1e16, while large
  # off-diagonal entries imply correlations as high as 0.99995 in magnitude.
  L <- matrix(
    c(
      1e-8, 0, 0, 0,
      2e-1, 2e-3, 0, 0,
      -4e3, 2e3, 5e2, 0,
      2e10, -1.2e10, 4e9, 2e8
    ),
    nrow = k, byrow = TRUE
  )
  rho <- matrix(
    c(
      0, 1, 0.2, 0.999999,
      1, 0, 0.8, 0.000001,
      0.00001, 0.99999, 0.5, 0.35,
      0.65, 0.15, 1, 0,
      0.93, 0.42, 0, 0.77,
      0.31, 0.58, 0.999999, 0.000001
    ),
    nrow = n, byrow = TRUE
  )
  ans <- .s2z_partial_change_case(L, rho)
  implied_cor <- cov2cor(tcrossprod(L))
  relative_zero_sum <- abs(colSums(ans$effects)) /
    pmax(1, apply(abs(ans$effects), 2L, max))

  expect_equal(range(diag(L)), c(1e-8, 2e8))
  expect_gt(max(abs(implied_cor[lower.tri(implied_cor)])), 0.9999)
  expect_true(all(c(0, 1) %in% rho))
  expect_true(any(rho > 0 & rho < 1))
  expect_true(all(is.finite(c(
    ans$effects, ans$restricted_map,
    ans$numeric_log_jacobian, ans$numeric_log_jacobian_svd,
    ans$formula_log_jacobian, ans$determinant_correction
  ))))
  expect_lt(max(relative_zero_sum), 5e-15)
  expect_gt(ans$balanced_condition, 1e7)
  expect_equal(
    ans$numeric_log_jacobian, ans$formula_log_jacobian,
    tolerance = 5e-8, scale = 1
  )
  expect_equal(
    ans$numeric_log_jacobian_svd, ans$formula_log_jacobian,
    tolerance = 5e-8, scale = 1
  )
  expect_equal(
    ans$numeric_log_jacobian, ans$numeric_log_jacobian_svd,
    tolerance = 2e-8, scale = 1
  )
})
