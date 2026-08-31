context("Tests for physical sum-to-zero group-level effects")

expect_match2 <- brms:::expect_match2

s2z_count_fixed <- function(x, pattern) {
  unname(lengths(regmatches(
    x, gregexpr(pattern, x, fixed = TRUE)
  )))
}

s2z_stan_between <- function(x, start, end) {
  start_at <- regexpr(start, x, fixed = TRUE)[1L]
  expect_gt(start_at, 0L)
  end_at <- regexpr(
    end, substring(x, start_at + nchar(start)), fixed = TRUE
  )[1L]
  expect_gt(end_at, 0L)
  substring(
    x, start_at + nchar(start),
    start_at + nchar(start) + end_at - 2L
  )
}

s2z_dat <- local({
  set.seed(1916)
  n <- 72
  data.frame(
    y = rnorm(n),
    x = seq(1.25, 4.75, length.out = n),
    z = rep(c(-0.5, 2), length.out = n),
    f = factor(rep(c("a", "b", "c"), length.out = n)),
    g = factor(rep(seq_len(6), each = 12)),
    h = factor(rep(seq_len(8), each = 9)),
    w = rep(seq(0.8, 1.2, length.out = 6), each = 12)
  )
})

s2z_ten_dat <- local({
  data.frame(
    y = seq(-1, 1, length.out = 80),
    ten = factor(rep(letters[1:10], 8)),
    g = factor(rep(seq_len(8), each = 10))
  )
})

test_that("S2Z centering API defaults compatibly and reaches the reframe", {
  default_form <- y ~ x + (1 + x | gr(g, s2z = TRUE))
  centered_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = TRUE))
  noncentered_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = FALSE))
  partial_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = 0.35))
  auto_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = "auto"))
  fisher_form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = "fisher"))

  default_terms <- brmsterms(default_form)
  centered_terms <- brmsterms(centered_form)
  noncentered_terms <- brmsterms(noncentered_form)
  partial_terms <- brmsterms(partial_form)
  auto_terms <- brmsterms(auto_form)
  fisher_terms <- brmsterms(fisher_form)
  expect_equal(default_terms$dpars$mu$re$gcall[[1]]$s2z_center, 1)
  expect_equal(centered_terms$dpars$mu$re$gcall[[1]]$s2z_center, 1)
  expect_equal(noncentered_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0)
  expect_equal(partial_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0.35)
  expect_equal(auto_terms$dpars$mu$re$gcall[[1]]$s2z_center, 0.5)
  expect_false(default_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_true(auto_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_true(fisher_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto)
  expect_equal(brms:::frame_re(default_terms, s2z_dat)$s2z_center, c(1, 1))
  expect_equal(brms:::frame_re(centered_terms, s2z_dat)$s2z_center, c(1, 1))
  expect_equal(
    brms:::frame_re(noncentered_terms, s2z_dat)$s2z_center, c(0, 0)
  )
  expect_equal(
    brms:::frame_re(partial_terms, s2z_dat)$s2z_center, c(0.35, 0.35)
  )
  expect_true(all(
    brms:::frame_re(auto_terms, s2z_dat)$s2z_center_auto
  ))
  auto_reframe <- brms:::frame_re(auto_terms, s2z_dat)
  expect_identical(brms:::re_s2z_center_mode(auto_reframe), "auto")
  expect_identical(
    stancode(fisher_form, data = s2z_dat),
    stancode(auto_form, data = s2z_dat)
  )

  default_code <- stancode(default_form, data = s2z_dat)
  centered_code <- stancode(centered_form, data = s2z_dat)
  noncentered_code <- stancode(noncentered_form, data = s2z_dat)
  expect_identical(default_code, centered_code)
  expect_false(identical(default_code, noncentered_code))
  expect_identical(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE, center = 1)),
      data = s2z_dat
    ),
    centered_code
  )
  expect_identical(
    stancode(
      y ~ x + (1 + x | gr(g, s2z = TRUE, center = 0)),
      data = s2z_dat
    ),
    noncentered_code
  )
  expect_null(standata(centered_form, data = s2z_dat)$rho_s2z_1)
  expect_null(standata(noncentered_form, data = s2z_dat)$rho_s2z_1)
  expect_null(standata(auto_form, data = s2z_dat)$rho_s2z_1)

  # Reframes made from the original S2Z formula representation did not carry
  # s2z_center. They must retain the original centered behavior.
  legacy_terms <- default_terms
  legacy_terms$dpars$mu$re$gcall[[1]]$s2z_center <- NULL
  legacy_terms$dpars$mu$re$gcall[[1]]$s2z_center_auto <- NULL
  expect_equal(
    brms:::frame_re(legacy_terms, s2z_dat)$s2z_center, c(1, 1)
  )
  expect_false(any(
    brms:::frame_re(legacy_terms, s2z_dat)$s2z_center_auto
  ))

  conventional <- stancode(y ~ x + (1 + x | gr(g)), data = s2z_dat)
  conventional_false <- stancode(
    y ~ x + (1 + x | gr(g, center = FALSE)), data = s2z_dat
  )
  conventional_zero <- stancode(
    y ~ x + (1 + x | gr(g, center = 0)), data = s2z_dat
  )
  expect_identical(conventional_false, conventional)
  expect_identical(conventional_zero, conventional)
  conventional_terms <- brmsterms(y ~ x + (1 + x | gr(g)))
  expect_equal(
    brms:::frame_re(conventional_terms, s2z_dat)$s2z_center, c(0, 0)
  )
  conventional_centered <- stancode(
    y ~ x + (1 + x | gr(g, center = TRUE)), data = s2z_dat
  )
  expect_identical(
    conventional_centered,
    stancode(y ~ x + (1 + x | gr(g, center = 1)), data = s2z_dat)
  )
  expect_false(identical(conventional_centered, conventional))
  conventional_partial <- stancode(
    y ~ x + (1 + x | gr(g, center = 0.35)), data = s2z_dat
  )
  expect_false(identical(conventional_partial, conventional_centered))
  expect_identical(
    stancode(
      y ~ x + (1 + x | gr(g, center = "fisher")), data = s2z_dat
    ),
    stancode(
      y ~ x + (1 + x | gr(g, center = "auto")), data = s2z_dat
    )
  )
  conventional_partial_terms <- brmsterms(
    y ~ x + (1 + x | gr(g, center = 0.35))
  )
  expect_equal(
    brms:::frame_re(
      conventional_partial_terms, s2z_dat
    )$s2z_center,
    c(0.35, 0.35)
  )
  # Old ordinary reframes did not carry centering metadata and must retain
  # their historical non-centered default.
  legacy_conventional <- conventional_terms
  legacy_conventional$dpars$mu$re$gcall[[1]]$s2z_center <- NULL
  legacy_conventional$dpars$mu$re$gcall[[1]]$s2z_center_auto <- NULL
  expect_equal(
    brms:::frame_re(legacy_conventional, s2z_dat)$s2z_center, c(0, 0)
  )
  expect_error(gr(g, s2z = TRUE, center = NA), "center")
  expect_error(
    gr(g, s2z = TRUE, center = c(TRUE, FALSE)), "center"
  )
  for (value in list(-0.01, 1.01, NaN, Inf, "AUTO", "partial")) {
    expect_error(gr(g, s2z = TRUE, center = value), "center")
  }

  mixed_data <- standata(
    y ~ x +
      (1 | gr(g, id = "mixed", s2z = TRUE, center = 0.2)) +
      (0 + x | gr(g, id = "mixed", s2z = TRUE, center = 0.8)),
    data = s2z_dat
  )
  expect_equal(unname(mixed_data$rho_s2z_1[, 1]), rep(0.2, 6))
  expect_equal(unname(mixed_data$rho_s2z_1[, 2]), rep(0.8, 6))

  expect_error(
    standata(
      y ~ x +
        (1 | gr(g, id = "mixed-auto", s2z = TRUE,
                center = "auto")) +
        (0 + x | gr(g, id = "mixed-auto", s2z = TRUE,
                    center = 0.8)),
      data = s2z_dat
    ),
    "must use Fisher centering if any coefficient does",
    fixed = TRUE
  )
})


test_that("stored partial charts follow stable blocks under re_formula", {
  form <- y ~ x +
    (1 | gr(g, s2z = TRUE, center = 0.2)) +
    (0 + x | gr(h, s2z = TRUE, center = 0.8))
  fit <- brm(form, data = s2z_dat, empty = TRUE, cores = 1)

  full <- standata(fit, internal = TRUE)
  expect_equal(unname(full$rho_s2z_1), matrix(0.2, 6L, 1L))
  expect_equal(unname(full$rho_s2z_2), matrix(0.8, 8L, 1L))

  restored <- unserialize(serialize(fit, NULL))
  restored_data <- standata(restored, internal = TRUE)
  expect_equal(unname(restored_data$rho_s2z_1), matrix(0.2, 6L, 1L))
  expect_equal(unname(restored_data$rho_s2z_2), matrix(0.8, 8L, 1L))

  incomplete <- fit
  colnames(
    incomplete$basis$dpars$mu$re_s2z_center$rho_s2z_1
  ) <- NULL
  expect_error(
    standata(incomplete, internal = TRUE),
    "invalid or incomplete grouping-level or coefficient labels"
  )

  retained <- standata(
    fit,
    re_formula = ~(0 + x | gr(h, s2z = TRUE, center = 0.8)),
    internal = TRUE
  )
  expect_equal(retained$M_1, 1L)
  expect_equal(unname(retained$rho_s2z_1), matrix(0.8, 8L, 1L))
  expect_null(retained$rho_s2z_2)
})


test_that("scalar Fisher S2Z hoists exposure and uses closed-form rho", {
  form <- y ~ 1 +
    (1 | gr(g, s2z = TRUE, center = "fisher"))
  scode <- stancode(form, data = s2z_dat)
  tdata <- s2z_stan_between(
    scode, "transformed data {", "\nparameters {"
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )

  expect_match2(
    scode, "matrix<lower=0,upper=1>[N_1, M_1] rho_s2z_1;"
  )
  expect_match2(
    scode, "vector<lower=0,upper=1>[M_1] mean_rho_s2z_1;"
  )
  expect_match2(
    tdata, "vector<lower=0>[N_1] exposure_fisher_s2z_1;"
  )
  expect_match2(tdata, "exposure_fisher_s2z_1 = zeros_vector(N_1);")
  expect_match2(
    tdata,
    "exposure_fisher_s2z_1[J_1[n]] += square(Z_1_1[n]);"
  )
  expect_match2(tpar, "inv_square(sigma)")
  expect_match2(tpar, "rho_s2z_1[j, 1]")
  expect_match2(tpar, "exposure_fisher_s2z_1[j]")
  expect_match2(
    tpar,
    paste0(
      "real scaled_info_fisher_s2z = 1.0 * square(sd_1[1]) * ",
      "obs_prec_fisher_s2z * exposure_fisher_s2z_1[j];"
    )
  )
  expect_match2(
    tpar,
    "rho_s2z_1[j, 1] = 1.0 - inv(1.0 + scaled_info_fisher_s2z);"
  )
  expect_false(grepl("J_1[n]", tpar, fixed = TRUE))
  expect_false(grepl("mdivide_left_spd", scode, fixed = TRUE))
  expect_false(grepl("quad_form(", scode, fixed = TRUE))
  expect_match2(scode, "scale_partial_s2z = 1.0 - rho_s2z_1[, 1]")
  expect_match2(scode, "+ log_det_partial_s2z_1")
  expect_false(grepl("eigenvectors_sym", scode, fixed = TRUE))
  expect_false(grepl("inverse_spd", scode, fixed = TRUE))
  expect_false(grepl("if (mean_rho_s2z_1", scode, fixed = TRUE))

  student_code <- stancode(form, data = s2z_dat, family = student())
  student_tpar <- s2z_stan_between(
    student_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    student_tpar,
    "(nu + 1.0) / (nu + 3.0) * inv_square(sigma)"
  )
  expect_match2(student_code, "+ log_det_partial_s2z_1")
})


test_that("multivariate Fisher S2Z hoists Gram matrices and solves diagonals", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = "fisher"))
  scode <- stancode(form, data = s2z_dat)
  tdata <- s2z_stan_between(
    scode, "transformed data {", "\nparameters {"
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )

  expect_match2(tdata, "array[N_1] matrix[M_1, M_1]")
  expect_match2(tdata, "J_1[n]")
  expect_match2(tdata, "Z_1_1[n]")
  expect_match2(tdata, "Z_1_2[n]")
  expect_false(grepl("J_1[n]", tpar, fixed = TRUE))
  expect_false(grepl("design_fisher_s2z", tpar, fixed = TRUE))
  expect_match2(tpar, "inv_square(sigma)")
  expect_match2(tpar, "cholesky_decompose(")
  expect_match2(tpar, "mdivide_left_tri_low(")
  expect_match2(tpar, "L_post_precision_fisher_s2z")
  expect_match2(tpar, "white_factor_fisher_s2z")
  expect_match2(tpar, "post_var_fisher_s2z = columns_dot_self(")
  expect_false(grepl("mdivide_left_spd", scode, fixed = TRUE))
  expect_equal(s2z_count_fixed(tpar, "quad_form("), 1L)
  expect_false(grepl("identity_matrix(M_1)", scode, fixed = TRUE))
  expect_false(grepl("white_cov_fisher_s2z", scode, fixed = TRUE))
  expect_false(grepl("post_cov_fisher_s2z", scode, fixed = TRUE))
  expect_match2(
    scode, "diag_pre_multiply(rho_s2z_1[j]', L_Sigma_s2z_1);"
  )
  expect_match2(scode, "+ log_det_partial_s2z_1")
  expect_null(standata(form, data = s2z_dat)$rho_s2z_1)

  independent_form <- y ~ x +
    (1 + x || gr(g, s2z = TRUE, center = "fisher"))
  independent_code <- stancode(independent_form, data = s2z_dat)
  independent_tdata <- s2z_stan_between(
    independent_code, "transformed data {", "\nparameters {"
  )
  independent_tpar <- s2z_stan_between(
    independent_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    independent_tdata,
    "array[N_1] matrix[M_1, M_1] gram_fisher_s2z_1;"
  )
  expect_false(grepl("J_1[n]", independent_tpar, fixed = TRUE))
  expect_match2(
    independent_tpar,
    "quad_form_diag(gram_fisher_s2z_1[j], sd_1)"
  )
  expect_match2(independent_tpar, "unit_rhs_fisher_s2z")
  expect_match2(independent_tpar, "dot_self(unit_column_fisher_s2z)")
  expect_equal(
    s2z_count_fixed(independent_tpar, "quad_form("), 0L
  )
  expect_false(grepl("diag_matrix(sd_1)", independent_tpar, fixed = TRUE))
  expect_false(grepl("identity_matrix(M_1)", independent_tpar, fixed = TRUE))
  expect_false(grepl(
    "mdivide_left_spd", independent_code, fixed = TRUE
  ))
  expect_false(grepl(
    "post_cov_fisher_s2z", independent_code, fixed = TRUE
  ))

})


test_that("response-free Fisher rules cover representative likelihoods", {
  dat <- data.frame(
    y = c(0L, 1L, 2L, 4L, 3L, 5L, 1L, 2L, 0L, 4L, 2L, 3L),
    trials = rep(6L, 12),
    g = factor(rep(seq_len(3), each = 4)),
    category = factor(rep(c("a", "b", "c"), 4)),
    p1 = seq(0.15, 0.25, length.out = 12),
    p2 = seq(0.25, 0.35, length.out = 12)
  )
  dat$p3 <- 1 - dat$p1 - dat$p2
  dat$simplex <- I(as.matrix(dat[c("p1", "p2", "p3")]))

  check_code <- function(form, family, markers, prior = empty_prior()) {
    scode <- stancode(form, data = dat, family = family, prior = prior)
    tpar <- s2z_stan_between(
      scode, "transformed parameters {", "\nmodel {"
    )
    sdata <- standata(form, data = dat, family = family, prior = prior)
    expect_false(any(startsWith(names(sdata), "rho_s2z_")))
    expect_false(grepl("Y[n]", tpar, fixed = TRUE))
    for (marker in markers) {
      expect_match2(tpar, marker)
    }
    invisible(tpar)
  }

  # Exact conditional expected information.
  check_code(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    poisson(), "real obs_prec_fisher_s2z = value_fisher_s2z_mu;"
  )

  # Three-bin coarsening: {0}, {trials}, and the interior counts.
  check_code(
    bf(
      y | trials(trials) ~ 1 +
        (1 | gr(g, s2z = TRUE, center = "fisher")),
      phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    beta_binomial(),
    c("pmid_fisher_s2z_bb", "dpmid_fisher_s2z_bb")
  )

  # Structural-atom decomposition for both the count and atom predictors.
  zip_prior <-
    prior(normal(0, 2), class = "Intercept") +
    prior(normal(0, 2), class = "Intercept", dpar = "zi")
  check_code(
    bf(
      y ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
      zi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    zero_inflated_poisson(),
    c("q0_fisher_s2z_zi", "atom_derivative_fisher_s2z_zi"),
    prior = zip_prior
  )

  # Marginal category information from the population-only softmax.
  check_code(
    category ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    categorical(), "prob_fisher_s2z_cat"
  )

  # Simplex-coordinate information from the Dirichlet trigamma identity.
  check_code(
    simplex ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    dirichlet(),
    c(
      "prob_fisher_s2z_dir",
      "alpha_fisher_s2z_dir < 1e-6 ? 1.0 :",
      "square(alpha_fisher_s2z_dir) * trigamma(alpha_fisher_s2z_dir)"
    )
  )

  # COM-Poisson moment information; at shape one the location variance
  # marker reduces to the Poisson mean information.
  check_code(
    bf(
      y ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
      shape ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    brmsfamily("com_poisson"),
    c(
      paste0(
        "log(fmax(value_fisher_s2z_mu, 1e-12)) / ",
        "fmax(value_fisher_s2z_shape, 1e-12);"
      ),
      paste0(
        "mode_fisher_s2z_cmp / ",
        "fmax(value_fisher_s2z_shape, 1e-12));"
      ),
      "real log_factorial_variance_fisher_s2z_cmp =",
      "real obs_prec_fisher_s2z = variance_fisher_s2z_cmp;",
      paste0(
        "square(derivative_fisher_s2z_shape) * ",
        "log_factorial_variance_fisher_s2z_cmp"
      )
    )
  )
})


test_that("response-free Fisher rules guard zero-trial observations", {
  dat <- data.frame(
    y = c(0L, 0L, 1L, 2L, 0L, 4L),
    trials = c(0L, 1L, 2L, 4L, 0L, 5L),
    g = factor(rep(seq_len(2), each = 3))
  )
  beta_binomial_form <- bf(
    y | trials(trials) ~ 1 +
      (1 | gr(g, s2z = TRUE, center = "fisher")),
    phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
  )
  beta_binomial_code <- stancode(
    beta_binomial_form, data = dat, family = beta_binomial()
  )
  beta_binomial_tpar <- s2z_stan_between(
    beta_binomial_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    beta_binomial_tpar,
    "(trials[n] == 0 ? 0.0 : (trials[n] == 1 ?"
  )
  expect_match2(
    beta_binomial_tpar,
    "(trials[n] <= 1 ? 0.0 :"
  )
  expect_match2(
    beta_binomial_tpar,
    paste0(
      "((value_fisher_s2z_mu) * ",
      "(1.0 - (value_fisher_s2z_mu)))"
    )
  )
  expect_match2(
    beta_binomial_tpar,
    paste0(
      "prob_fisher_s2z_bb = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, value_fisher_s2z_mu));"
    )
  )
  beta_binomial_data <- standata(
    beta_binomial_form, data = dat, family = beta_binomial()
  )
  expect_equal(as.vector(beta_binomial_data$trials), dat$trials)
  expect_false(any(startsWith(
    names(beta_binomial_data), "rho_s2z_"
  )))
  expect_false(grepl("Y[n]", beta_binomial_tpar, fixed = TRUE))

  dm_dat <- data.frame(
    trials = c(0L, 3L, 4L, 0L, 5L, 2L),
    g = factor(rep(seq_len(2), each = 3))
  )
  dm_dat$counts <- I(rbind(
    c(0L, 0L, 0L), c(1L, 1L, 1L), c(1L, 2L, 1L),
    c(0L, 0L, 0L), c(2L, 1L, 2L), c(1L, 0L, 1L)
  ))
  dm_form <- bf(
    counts | trials(trials) ~ 1 +
      (1 | gr(g, s2z = TRUE, center = "fisher")),
    phi ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
  )
  dm_code <- stancode(
    dm_form, data = dm_dat, family = dirichlet_multinomial()
  )
  dm_tpar <- s2z_stan_between(
    dm_code, "transformed parameters {", "\nmodel {"
  )
  expect_match2(
    dm_tpar,
    paste0(
      "((trials[n]) * (1.0 + (value_fisher_s2z_phi)) / ",
      "((trials[n]) + (value_fisher_s2z_phi)))"
    )
  )
  expect_match2(dm_tpar, "(trials[n] == 0 ? 0.0 : 0.5 *")
  dm_data <- standata(
    dm_form, data = dm_dat, family = dirichlet_multinomial()
  )
  expect_equal(as.vector(dm_data$trials), dm_dat$trials)
  expect_false(any(startsWith(names(dm_data), "rho_s2z_")))
  expect_false(grepl("Y[n]", dm_tpar, fixed = TRUE))
})


test_that("Fisher link algebra is stable at parameter boundaries", {
  dat <- data.frame(
    y_unit = seq(0.1, 0.9, length.out = 8),
    y_pos = seq(0.3, 2, length.out = 8),
    y_bin = rep(0:1, 4),
    g = factor(rep(seq_len(2), each = 4))
  )

  xbeta_code <- stancode(
    y_unit ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    data = dat, family = xbeta()
  )
  for (marker in c(
    paste0(
      "prob_fisher_s2z_xbeta = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, value_fisher_s2z_mu));"
    ),
    paste0(
      "(1.0 + (value_fisher_s2z_phi)) * (value_fisher_s2z_mu) * ",
      "(1.0 - (value_fisher_s2z_mu))"
    )
  )) {
    expect_match2(xbeta_code, marker)
  }

  frechet_code <- stancode(
    bf(
      y_pos ~ 1,
      nu ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
    ),
    data = dat, family = frechet()
  )
  for (marker in c(
    "real boundary_fraction_fisher_s2z_frechet =",
    "boundary_fraction_fisher_s2z_frechet < 1e-6 ?",
    "real shape_info_fisher_s2z_frechet =",
    "real obs_prec_fisher_s2z = shape_info_fisher_s2z_frechet;"
  )) {
    expect_match2(frechet_code, marker)
  }

  softit_code <- stancode(
    y_bin ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    data = dat, family = bernoulli("softit")
  )
  expect_match2(softit_code, "return log(expm1(-p ./ (p - 1)));")
  expect_match2(softit_code, "return log1p_exp(y) ./ (1 + log1p_exp(y));")
})


test_that("Wiener Fisher rules use exact decision coarsening", {
  dat <- data.frame(
    q = seq(0.5, 1.6, length.out = 12),
    decision = rep(0:1, 6),
    g = factor(rep(seq_len(3), each = 4))
  )
  form <- bf(
    q | dec(decision) ~ 1,
    bs ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    bias ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher")),
    ndt = 0.1
  )
  s2z_prior <-
    prior(normal(0, 2), class = "Intercept", dpar = "bs") +
    prior(normal(0, 2), class = "Intercept", dpar = "bias")
  scode <- stancode(
    form, data = dat, family = wiener(), prior = s2z_prior
  )
  tpar <- s2z_stan_between(
    scode, "transformed parameters {", "\nmodel {"
  )
  for (marker in c(
    "real choice_scale_fisher_s2z_wiener =",
    "if (abs(choice_scale_fisher_s2z_wiener) < 1e-5)",
    "real dp_dscale_fisher_s2z_wiener;",
    "real dp_dbias_fisher_s2z_wiener;",
    paste0(
      "prob_safe_fisher_s2z_wiener = fmin(1.0 - 1e-12, ",
      "fmax(1e-12, p_upper_fisher_s2z_wiener));"
    ),
    paste0(
      "square((derivative_fisher_s2z_bs) * 2.0 * ",
      "(value_fisher_s2z_mu) * dp_dscale_fisher_s2z_wiener)"
    ),
    paste0(
      "square((derivative_fisher_s2z_bias) * ",
      "dp_dbias_fisher_s2z_wiener)"
    )
  )) {
    expect_match2(tpar, marker)
  }
  expect_false(grepl("Y[n]", tpar, fixed = TRUE))
  expect_false(grepl("dec[n]", tpar, fixed = TRUE))
  sdata <- standata(
    form, data = dat, family = wiener(), prior = s2z_prior
  )
  expect_false(any(startsWith(names(sdata), "rho_s2z_")))

  expect_error(
    stancode(
      bf(
        q | dec(decision) ~ 1,
        ndt ~ 1 + (1 | gr(g, s2z = TRUE, center = "fisher"))
      ),
      data = dat, family = wiener()
    ),
    paste0(
      "has no response-free expected-information rule for family ",
      "'wiener' and distributional parameter 'ndt'"
    ),
    fixed = TRUE
  )
})


test_that("partial S2Z log determinants are stable at centering endpoints", {
  rho <- c(0, 1 - 1e-8, 1)
  scale <- rep(1e-20, length(rho))
  stable_log_term <- vapply(seq_along(rho), function(i) {
    if (rho[i] == 1) {
      log(scale[i])
    } else {
      log1p(-rho[i] * (1 - scale[i]))
    }
  }, numeric(1))
  reference <- log((1 - rho) + rho * scale)

  expect_equal(stable_log_term, reference, tolerance = 1e-10)
  expect_true(all(is.finite(stable_log_term)))
  expect_identical(log1p(-(1 - 1e-20)), -Inf)
})


test_that("direct non-centered scalar S2Z scales contrasts exactly", {
  centered <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE)), data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  direct <- stancode(
    y ~ 1 + (1 | gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )

  expect_match2(
    centered,
    "r_s2z_1_1 = sum_to_zero_constrain_brms(z_s2z_1);"
  )
  expect_match2(
    direct,
    paste0(
      "r_s2z_1_1 = sum_to_zero_constrain_brms(",
      "sd_1[1] * z_s2z_1);"
    )
  )
  expect_match2(centered, "- (N_1 - 1) * log(sd_1[1])")
  expect_false(grepl(
    "- (N_1 - 1) * log(sd_1[1])", direct, fixed = TRUE
  ))
  for (term in c(
    "-0.5 * group_quad_s2z_1",
    "- 0.5 * log(D_s2z_1)",
    "+ 0.5 * log(1.0 * N_1)",
    "mean_r_s2z_1 = mhat_s2z_1 + sd_1[1] * std_normal_rng()",
    "q_recovered_s2z_1 = theta_s2z - H_s2z_1 * mean_r_s2z_1",
    "r_1_1 = r_s2z_1_1 + mean_r_s2z_1"
  )) {
    expect_true(grepl(term, direct, fixed = TRUE), info = term)
  }

  student <- stancode(
    y ~ 1 + (1 | gr(
      g, s2z = TRUE, center = FALSE, dist = "student"
    )),
    data = s2z_dat,
    prior = prior(normal(0, 2), class = Intercept)
  )
  expect_match2(
    student,
    paste0(
      "r_s2z_1_1 = sum_to_zero_constrain_brms(",
      "sd_1[1] * z_s2z_1);"
    )
  )
  expect_false(grepl(
    "- (N_1 - 1) * log(sd_1[1])", student, fixed = TRUE
  ))
  for (term in c(
    "dot_product(r_s2z_1_1, group_prec_s2z_1)",
    "- sum(log(group_scale_s2z_1))",
    "- 0.5 * log(D_s2z_1)",
    "+ 0.5 * log(1.0 * N_1)"
  )) {
    expect_true(grepl(term, student, fixed = TRUE), info = term)
  }
})


test_that("direct independent S2Z scales K4 and K10 component-wise", {
  four_centered <- stancode(
    y ~ x * z + (1 + x * z || gr(g, s2z = TRUE)),
    data = s2z_dat
  )
  four_direct <- stancode(
    y ~ x * z +
      (1 + x * z || gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_dat
  )
  ten_direct <- stancode(
    y ~ 0 + ten +
      (0 + ten || gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_ten_dat,
    prior = prior(normal(0, 2), class = b)
  )

  expect_match2(four_centered, "- (N_1 - 1) * sum(log(sd_1))")
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(sd_1))", four_direct, fixed = TRUE
  ))
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(sd_1))", ten_direct, fixed = TRUE
  ))
  for (code_and_k in list(c(four_direct, 4L), c(ten_direct, 10L))) {
    code <- code_and_k[[1]]
    k <- as.integer(code_and_k[[2]])
    for (j in seq_len(k)) {
      expect_match2(
        code,
        sprintf(
          paste0(
            "r_s2z_1_%1$s = sum_to_zero_constrain_brms(",
            "sd_1[%1$s] * segment(z_s2z_1, ",
            "(%1$s - 1) * (N_1 - 1) + 1, N_1 - 1));"
          ),
          j
        )
      )
    }
    for (term in c(
      "- 0.5 * sum(log(D_diag_s2z_1))",
      "- 0.5 * log1p(rank1_info_s2z_1)",
      "q_recovered_s2z_1", "mean_r_s2z_1"
    )) {
      expect_true(grepl(term, code, fixed = TRUE), info = term)
    }
    expect_false(grepl("matrix[M_1, M_1]", code, fixed = TRUE))
    expect_false(grepl("cholesky_decompose(", code, fixed = TRUE))
  }

  student_direct <- stancode(
    y ~ x * z + (1 + x * z || gr(
      g, s2z = TRUE, center = FALSE, dist = "student"
    )),
    data = s2z_dat
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(sd_1))", student_direct, fixed = TRUE
  ))
  for (term in c(
    paste0(
      "sd_1[4] * segment(z_s2z_1, ",
      "(4 - 1) * (N_1 - 1) + 1, N_1 - 1)"
    ),
    "- M_1 * sum(log(group_scale_s2z_1))",
    "- 0.5 * sum(log(D_diag_s2z_1))",
    "- 0.5 * log1p(rank1_info_s2z_1)"
  )) {
    expect_true(grepl(term, student_direct, fixed = TRUE), info = term)
  }
})


test_that("direct correlated S2Z applies the reference Cholesky", {
  centered <- stancode(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)),
    data = s2z_dat
  )
  direct <- stancode(
    y ~ x * z +
      (1 + x * z | gr(g, s2z = TRUE, center = FALSE)),
    data = s2z_dat
  )

  expect_match2(
    centered,
    paste0(
      "r_s2z_1[, k] = sum_to_zero_constrain_brms(segment(z_s2z_1, ",
      "(k - 1) * (N_1 - 1) + 1, N_1 - 1));"
    )
  )
  expect_false(grepl(
    "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';",
    centered, fixed = TRUE
  ))
  expect_match2(
    direct, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';"
  )
  expect_match2(
    centered,
    "- N_1 * sum(log(diagonal(L_Sigma_s2z_1)))"
  )
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    direct, fixed = TRUE
  ))
  for (term in c(
    "-0.5 * group_quad_s2z_1",
    "- sum(log(diagonal(L_P_s2z_1)))",
    "+ 0.5 * M_1 * log(1.0 * N_1)",
    "q_recovered_s2z_1 = theta_s2z - H_s2z_1 * mean_r_s2z_1",
    "r_1 = r_s2z_1;",
    "for (j in 1:N_1) r_1[j] += mean_r_s2z_1';"
  )) {
    expect_true(grepl(term, direct, fixed = TRUE), info = term)
  }

  student <- stancode(
    y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = FALSE, dist = "student"
    )),
    data = s2z_dat
  )
  expect_match2(student, "r_s2z_1 = r_s2z_1 * L_Sigma_s2z_1';")
  expect_false(grepl(
    "- (N_1 - 1) * sum(log(diagonal(L_Sigma_s2z_1)))",
    student, fixed = TRUE
  ))
  for (term in c(
    "contrast_score_s2z_1 = white_s2z * group_prec_s2z_1",
    "- M_1 * sum(log(group_scale_s2z_1))",
    "- sum(log(diagonal(L_P_s2z_1)))",
    "+ 0.5 * M_1 * log(1.0 * N_1)"
  )) {
    expect_true(grepl(term, student, fixed = TRUE), info = term)
  }
})


test_that("partial S2Z Jacobian state follows save_pars", {
  form <- y ~ x +
    (1 + x | gr(g, s2z = TRUE, center = 0.4, dist = "student"))
  default_fit <- brm(form, data = s2z_dat, empty = TRUE, cores = 1)
  saved_fit <- brm(
    form, data = s2z_dat, empty = TRUE,
    save_pars = save_pars(all = TRUE), cores = 1
  )
  default_excluded <- unlist(
    brms:::exclude_pars(default_fit), use.names = FALSE
  )
  saved_excluded <- unlist(
    brms:::exclude_pars(saved_fit), use.names = FALSE
  )

  expect_true("log_det_partial_s2z_1" %in% default_excluded)
  expect_false("log_det_partial_s2z_1" %in% saved_excluded)
  expect_true("contrast_score_s2z_1" %in% default_excluded)
  expect_false("contrast_score_s2z_1" %in% saved_excluded)
})


test_that("Matheron supports overlapping blocks and all centering modes", {
  form <- y ~ x * z +
    (1 + x | gr(g, s2z = TRUE, center = TRUE)) +
    (1 + x * z || gr(h, s2z = TRUE, center = FALSE)) +
    (1 | gr(f, s2z = TRUE, center = 0.4))
  bprior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b)
  scode <- stancode(form, data = s2z_dat, prior = bprior)

  for (term in c(
    "// fast Gaussian Matheron system for S2Z blocks 1, 2, 3",
    "matrix[4, 4] W_matheron_s2z_1;",
    "matrix[4, 4] L_W_matheron_s2z_1;",
    "W_matheron_s2z_1 = add_diag(",
    "square(prior_scale_s2z_1[{1, 2, 3, 4}])",
    "H_active_s2z = H_s2z_1[{1, 2, 3, 4}, ];",
    paste0(
      "theta_difference_s2z = theta_s2z[{1, 2, 3, 4}] - ",
      "prior_mean_s2z_1[{1, 2, 3, 4}];"
    ),
    "- (N_2 - 1) * sum(log(diagonal(L_Sigma_s2z_2)))",
    "+ log_det_partial_s2z_1",
    "- 0.5 * (N_1 - 1) * M_1 * log(2 * pi())",
    "- 0.5 * (N_2 - 1) * M_2 * log(2 * pi())",
    "- 0.5 * (N_3 - 1) * M_3 * log(2 * pi())",
    "- 0.5 * 4 * log(2 * pi())",
    "mean_r_s2z_1 += L_Sigma_s2z_1 *",
    "mean_r_s2z_2 += L_Sigma_s2z_2 *",
    "mean_r_s2z_3 += L_Sigma_s2z_3 *",
    "group_quad_s2z_3 = dot_self(z_s2z_3);",
    "q_recovered_s2z_1 = theta_s2z;"
  )) {
    expect_true(grepl(term, scode, fixed = TRUE), info = term)
  }
  expect_false(grepl(
    "- (N_3 - 1) * sum(log(diagonal(L_Sigma_s2z_3)))",
    scode, fixed = TRUE
  ))
  expect_false(grepl("log_det_partial_s2z_3", scode, fixed = TRUE))
  expect_false(grepl(
    "r_s2z_3 ./ rep_matrix", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "L_Sigma_s2z_3, r_s2z_3", scode, fixed = TRUE
  ))
  expect_false(grepl("matrix[7, 7] P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("L_P_s2z_1", scode, fixed = TRUE))
  expect_false(grepl("H_joint_s2z_1", scode, fixed = TRUE))
  expect_false(grepl(
    "W_matheron_s2z_1 = rep_matrix", scode, fixed = TRUE
  ))
  expect_false(grepl(
    "diag_matrix(square(prior_scale_s2z_1", scode, fixed = TRUE
  ))
  expect_equal(
    s2z_count_fixed(
      scode,
      "W_matheron_s2z_1 += tcrossprod(H_active_s2z * "
    ),
    2L
  )
  expect_equal(s2z_count_fixed(scode, "cholesky_decompose("), 1L)

  selective_prior <- prior(normal(0, 2), class = Intercept) +
    prior(normal(0, 1), class = b, coef = "x:z")
  selective_code <- stancode(form, data = s2z_dat, prior = selective_prior)
  for (term in c(
    "W_matheron_s2z_1 = add_diag(",
    "square(prior_scale_s2z_1[{1, 4}])",
    "H_active_s2z = H_s2z_1[{1, 4}, ];",
    paste0(
      "theta_difference_s2z = theta_s2z[{1, 4}] - ",
      "prior_mean_s2z_1[{1, 4}];"
    )
  )) {
    expect_true(grepl(term, selective_code, fixed = TRUE), info = term)
  }
})


test_that("explicit centered S2Z preserves legacy code for Plan 03 kernels", {
  cases <- list(
    scalar = list(
      data = s2z_dat,
      default = y ~ 1 + (1 | gr(g, s2z = TRUE)),
      explicit = y ~ 1 + (1 | gr(g, s2z = TRUE, center = TRUE))
    ),
    independent = list(
      data = s2z_ten_dat,
      default = y ~ 0 + ten + (0 + ten || gr(g, s2z = TRUE)),
      explicit = y ~ 0 + ten +
        (0 + ten || gr(g, s2z = TRUE, center = TRUE))
    ),
    correlated = list(
      data = s2z_dat,
      default = y ~ x * z + (1 + x * z | gr(g, s2z = TRUE)),
      explicit = y ~ x * z +
        (1 + x * z | gr(g, s2z = TRUE, center = TRUE))
    ),
    student = list(
      data = s2z_dat,
      default = y ~ x +
        (1 + x | gr(g, s2z = TRUE, dist = "student")),
      explicit = y ~ x + (1 + x | gr(
        g, s2z = TRUE, center = TRUE, dist = "student"
      ))
    )
  )
  for (case in cases) {
    expect_identical(
      stancode(case$explicit, data = case$data),
      stancode(case$default, data = case$data)
    )
  }
})

test_that("fixed partial S2Z covers scalar independent correlated and Student", {
  scalar_form <- y ~ 1 + (1 | gr(g, s2z = TRUE, center = 0.35))
  scalar <- stancode(scalar_form, data = s2z_dat)
  scalar_data <- standata(scalar_form, data = s2z_dat)
  expect_equal(
    unname(scalar_data$rho_s2z_1),
    matrix(0.35, nrow = 6L, ncol = 1L)
  )
  expect_match2(scalar, "real log_det_partial_s2z_1;")
  expect_match2(scalar, "vector[N_1] scale_partial_s2z")
  expect_match2(scalar, "+ log_det_partial_s2z_1")

  independent <- stancode(
    y ~ 0 + ten + (0 + ten || gr(g, s2z = TRUE, center = 0.62)),
    data = s2z_ten_dat, prior = prior(normal(0, 2), class = b)
  )
  expect_match2(independent, "rho_s2z_1[, 10]")
  expect_match2(independent, "r_s2z_1_10 = centered_partial_s2z")
  expect_false(grepl("matrix[M_1, M_1] L_partial_s2z", independent,
                     fixed = TRUE))

  correlated <- stancode(
    y ~ x * z + (1 + x * z | gr(g, s2z = TRUE, center = 0.44)),
    data = s2z_dat
  )
  expect_match2(
    correlated,
    "diag_pre_multiply(rho_s2z_1[j]', L_Sigma_s2z_1)"
  )
  expect_match2(correlated, "mean_partial_s2z /= N_1;")
  expect_match2(correlated, "+ log_det_partial_s2z_1")

  student <- stancode(
    y ~ x + (1 + x | gr(
      g, s2z = TRUE, center = 0.44, dist = "student"
    )),
    data = s2z_dat
  )
  expect_match2(student, "group_prec_s2z_1 = inv_square(group_scale_s2z_1)")
  expect_match2(student, "- M_1 * sum(log(group_scale_s2z_1))")
  expect_match2(student, "+ log_det_partial_s2z_1")
})

test_that("Fisher centering stays local and response-free", {
  expect_error(
    stancode(
      y | weights(w) ~ 1 +
        (1 | gr(g, s2z = TRUE, center = "fisher")),
      data = s2z_dat
    ),
    "does not yet support response addition term 'weights'",
    fixed = TRUE
  )

  ordinal_data <- transform(
    s2z_dat,
    yo = ordered(rep(c("low", "middle", "high"), length.out = nrow(s2z_dat)))
  )
  expect_error(
    stancode(
      yo ~ x + (1 | gr(g, s2z = TRUE, center = "fisher")),
      data = ordinal_data, family = cumulative()
    ),
    "not yet supported for ordinal sum-to-zero location predictors",
    fixed = TRUE
  )

  mv_data <- transform(s2z_dat, y2 = y + 0.2 * x)
  mv_form <- bf(
    y ~ 1 + (1 | gr(g, id = "y_block", s2z = TRUE,
                    center = "fisher"))
  ) +
    bf(
      y2 ~ 1 + (1 | gr(g, id = "y2_block", s2z = TRUE,
                       center = "fisher"))
    ) +
    set_rescor(FALSE)
  mv_code <- stancode(mv_form, data = mv_data)
  expect_match2(mv_code, "rho_s2z_1[j, 1]")
  expect_match2(mv_code, "rho_s2z_2[j, 1]")
  tpar <- s2z_stan_between(
    mv_code, "transformed parameters {", "\nmodel {"
  )
  expect_false(grepl("Y[n]", tpar, fixed = TRUE))
})
