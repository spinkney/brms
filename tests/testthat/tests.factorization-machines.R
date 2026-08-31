context("Factorization-machine terms")

fm_test_data <- function() {
  out <- expand.grid(
    user = factor(letters[1:4]),
    item = factor(LETTERS[1:3])
  )
  out <- out[-c(2L, 7L), , drop = FALSE]
  out$y <- seq_len(nrow(out)) / 10
  out
}

fm_matrix_names <- function(name, nrow, ncol) {
  as.vector(outer(
    seq_len(nrow), seq_len(ncol),
    function(i, j) sprintf("%s[%d,%d]", name, i, j)
  ))
}

test_that("fm validates and parses two categorical fields", {
  term <- fm(user, item, k = 2)
  expect_s3_class(term, "fm_term")
  expect_equal(term$term, c("user", "item"))
  expect_equal(term$k, 2L)
  expect_true(term$main)

  expect_error(fm(user, user), "must be different")
  expect_error(fm(log(user), item), "single untransformed variables")
  expect_error(fm(user, item, k = 0), "positive integer")
  expect_error(fm(user, item, k = 1.5), "positive integer")
  expect_error(fm(user, item, k = Inf), "positive integer")
  expect_error(fm(user, item, k = .Machine$integer.max + 1), "positive integer")
  expect_error(fm(user, item, main = NA), "single logical")
  expect_error(
    brmsterms(y ~ (fm(user, item, k = 2) | group)),
    "factorization machines"
  )

  bterms <- brmsterms(y ~ x + fm(user, item, k = 2))
  expect_true(is.formula(bterms$dpars$mu$fm))
  expect_equal(all.vars(bterms$dpars$mu$fm), c("user", "item"))
  expect_false(any(grepl("fm\\(", all_terms(bterms$dpars$mu$fe))))
  expect_equal(
    sort(all.vars(bterms$allvars)),
    sort(c("y", "x", "user", "item"))
  )
})

test_that("fm standata stores reflector sizes and marginal normalizations", {
  dat <- fm_test_data()
  sdata <- standata(y ~ fm(user, item, k = 1), dat)

  expect_equal(sdata$Nfm_1_1, 4L)
  expect_equal(sdata$Nfm_1_2, 3L)
  expect_equal(sdata$Kfm_1, 1L)
  expect_equal(sdata$Sfm_1_1, 0L)
  expect_equal(sdata$Sfm_1_2, 0L)
  expect_equal(sdata$Mfm_1_1, 3L)
  expect_equal(sdata$Mfm_1_2, 2L)
  expect_equal(sdata$Jfm_1_1, as.array(as.integer(dat$user)))
  expect_equal(sdata$Jfm_1_2, as.array(as.integer(dat$item)))
  expect_equal(sdata$Cfm_main_1_1, sqrt(4 / 3))
  expect_equal(sdata$Cfm_main_1_2, sqrt(3 / 2))
  expect_equal(sdata$Cfm_1, sqrt(4 * 3))

  no_main <- standata(y ~ fm(user, item, k = 1, main = FALSE), dat)
  expect_null(no_main$Cfm_main_1_1)
  expect_null(no_main$Cfm_main_1_2)

  tiny <- expand.grid(
    user = factor(letters[1:2]), item = factor(LETTERS[1:2])
  )
  tiny$y <- 0
  tiny_data <- standata(y ~ fm(user, item, k = 1), tiny)
  expect_equal(tiny_data$Sfm_1_1, 1L)
  expect_equal(tiny_data$Sfm_1_2, 0L)
  expect_equal(tiny_data$Mfm_1_1, 0L)
  expect_equal(tiny_data$Mfm_1_2, 1L)

  expect_error(
    standata(y ~ fm(user, item, k = 3), dat),
    "must not exceed 2"
  )
  expect_error(
    standata(
      y ~ fm(user, item, k = 1, main = FALSE) +
        fm(item, user, k = 1, main = FALSE),
      dat
    ),
    "Duplicated factorization-machine field pairs"
  )
  expect_error(
    standata(y ~ user + fm(user, item, k = 1), dat),
    "Set 'main = FALSE'"
  )
  expect_error(
    standata(y ~ factor(user) + fm(user, item, k = 1), dat),
    "Set 'main = FALSE'"
  )
  expect_error(
    standata(y ~ as.factor(item) + fm(user, item, k = 1), dat),
    "Set 'main = FALSE'"
  )
  expect_error(
    standata(
      y ~ user:item + fm(user, item, k = 1, main = FALSE),
      dat
    ),
    "duplicates the factorization-machine interaction"
  )
  expect_error(
    standata(
      y ~ factor(user):as.factor(item) +
        fm(user, item, k = 1, main = FALSE),
      dat
    ),
    "duplicates the factorization-machine interaction"
  )
  dat$x <- seq_len(nrow(dat))
  expect_silent(standata(
    y ~ user:x + fm(user, item, k = 1, main = FALSE), dat
  ))
  expect_silent(standata(
    y ~ I(as.numeric(user) + as.numeric(item)) +
      fm(user, item, k = 1, main = FALSE),
    dat
  ))
  bad <- transform(dat, user = as.numeric(user) + 0.5)
  expect_error(
    standata(y ~ fm(user, item, k = 2), bad),
    "integer-valued identifier"
  )
  bad_infinite <- transform(dat, user = as.numeric(user))
  bad_infinite$user[1] <- Inf
  expect_error(
    suppressWarnings(standata(y ~ fm(user, item, k = 2), bad_infinite)),
    "integer-valued identifier"
  )
})

test_that("conditional-effect grids treat numeric FM identifiers as levels", {
  dat <- expand.grid(
    user = c(10L, 30L, 80L),
    item = c(101L, 305L, 900L)
  )
  dat$y <- seq_len(nrow(dat)) / 10
  fit <- brm(
    y ~ fm(user, item, k = 2), dat,
    empty = TRUE, backend = "mock", cores = 1
  )
  effects <- list("user")
  conditions <- prepare_conditions(fit, effects = effects)
  expect_true(conditions$item %in% unique(dat$item))

  bterms <- brmsterms(fit$formula)
  fm_vars <- get_fm_vars(bterms)
  cond_data <- prepare_cond_data(
    model.frame(fit)[, "user", drop = FALSE],
    conditions = conditions,
    factor_vars = fm_vars
  )
  expect_equal(sort(unique(as.numeric(as.character(cond_data$user)))),
               sort(unique(dat$user)))
  expect_true(all(cond_data$item %in% unique(dat$item)))
  expect_identical(unname(attr(cond_data, "types")), "factor")

  made <- make_conditions(fit, "item")
  expect_equal(sort(as.numeric(as.character(made$item))),
               sort(unique(dat$item)))

  new_frame <- brmsframe(bterms, cond_data, basis = fit$basis)
  expect_silent(data_fm(new_frame$dpars$mu, cond_data))
})

test_that("fm priors and semi-orthogonal Stan code are generated correctly", {
  dat <- fm_test_data()
  priors <- get_prior(y ~ fm(user, item, k = 2), dat)
  expect_true(any(priors$class == "sdfm" & priors$coef == "user:item"))
  expect_true(any(priors$class == "sdfm_main" & priors$group == "user"))
  expect_true(any(priors$class == "sdfm_main" & priors$group == "item"))
  sdata <- standata(y ~ fm(user, item, k = 2), dat)
  expect_equal(sdata$Sfm_1_1, 0L)
  expect_equal(sdata$Sfm_1_2, 1L)
  expect_equal(sdata$Mfm_1_1, 5L)
  expect_equal(sdata$Mfm_1_2, 2L)

  prior <- prior(normal(0, 1), class = "sdfm", coef = "user:item") +
    prior(exponential(2), class = "sdfm_main", group = "user")
  scode <- suppressWarnings(stancode(
    y ~ fm(user, item, k = 2), dat, prior = prior,
    backend = "mock", parse = FALSE
  ))
  expect_match2(
    scode,
    "matrix fm_apply_reflector_brms(real tau, vector v, matrix B)"
  )
  expect_match2(
    scode,
    "vector[Mfm_1_1] zfm_frame_1_1;"
  )
  expect_match2(scode, "simplex[Kfm_1] zfm_spectrum_1;")
  expect_match2(
    scode,
    paste0(
      "matrix[Nfm_1_1, Kfm_1] Qfm_1_1 = ",
      "fm_centered_semiorthogonal_constrain_brms("
    )
  )
  expect_match2(scode, "sum_to_zero_vector[Nfm_1_1] zfm_main_1_1;")
  expect_match2(scode, "dot_product(Qfm_1_1[Jfm_1_1[n]]")
  expect_match2(
    scode,
    "lprior += dirichlet_lpdf(zfm_spectrum_1 | rep_vector(1, Kfm_1))"
  )
  expect_match2(scode, "lprior += std_normal_lpdf(zfm_frame_1_1)")
  expect_match2(scode, "lprior += std_normal_lpdf(zfm_main_1_1)")
  expect_match2(scode, "normal_lpdf(sdfm_1 | 0, 1)")
  expect_match2(scode, "exponential_lpdf(sdfm_main_1_1 | 2)")

  no_main <- stancode(
    y ~ fm(user, item, k = 2, main = FALSE), dat,
    backend = "mock", normalize = FALSE, parse = FALSE
  )
  expect_false(grepl("zfm_main", no_main, fixed = TRUE))
  expect_match2(
    no_main,
    "lprior += std_normal_lupdf(zfm_frame_1_1)"
  )
  expect_match2(no_main, "dirichlet_lupdf(zfm_spectrum_1")

  constant_code <- suppressWarnings(stancode(
    y ~ fm(user, item, k = 2), dat,
    prior = prior(constant(1.2), class = "sdfm", coef = "user:item"),
    backend = "mock", parse = FALSE
  ))
  assignments <- regmatches(
    constant_code,
    gregexpr("sdfm_1 = 1.2;", constant_code, fixed = TRUE)
  )
  expect_length(assignments[[1]], 1L)

  prior_draw_code <- stancode(
    y ~ fm(user, item, k = 1), dat, sample_prior = "yes",
    backend = "mock", parse = FALSE
  )
  expect_false(grepl("prior_zfm_frame", prior_draw_code, fixed = TRUE))
  expect_false(grepl("prior_zfm_main", prior_draw_code, fixed = TRUE))
  expect_false(grepl("prior_zfm_spectrum", prior_draw_code, fixed = TRUE))
  expect_match2(
    prior_draw_code,
    "real prior_sdfm_1 = student_t_rng(3,0,2.5);"
  )
  expect_match2(
    prior_draw_code,
    "real prior_sdfm_main_1_1 = student_t_rng(3,0,2.5);"
  )
  expect_match2(
    prior_draw_code,
    "real prior_sdfm_main_1_2 = student_t_rng(3,0,2.5);"
  )
  expect_match2(
    prior_draw_code,
    "lprior += std_normal_lpdf(zfm_frame_1_1)"
  )
  expect_match2(
    prior_draw_code,
    "lprior += dirichlet_lpdf(zfm_spectrum_1"
  )
  expect_match2(
    prior_draw_code,
    "lprior += student_t_lpdf(sdfm_1 | 3, 0, 2.5)"
  )
  expect_match2(
    prior_draw_code,
    "lprior += student_t_lpdf(sdfm_main_1_1 | 3, 0, 2.5)"
  )
  expect_match2(
    prior_draw_code,
    "lprior += student_t_lpdf(sdfm_main_1_2 | 3, 0, 2.5)"
  )
})

test_that("fm terms support threaded likelihood code", {
  dat <- fm_test_data()
  scode <- stancode(
    y ~ fm(user, item, k = 2), dat, threads = threading(2),
    backend = "mock", parse = FALSE
  )
  expect_match2(scode, "data array[] int Jfm_1_1")
  expect_match2(scode, "matrix Qfm_1_1")
  expect_match2(scode, "vector fm_singular_1")
  expect_match2(scode, "Jfm_1_1[nn]")
  expect_match2(scode, "reduce_sum(partial_log_lik_lpmf")
})

test_that("fm code keeps distributional parameter prefixes distinct", {
  dat <- fm_test_data()
  scode <- stancode(
    bf(
      y ~ fm(user, item, k = 2),
      sigma ~ fm(user, item, k = 1, main = FALSE)
    ),
    dat, backend = "mock", parse = FALSE
  )
  expect_match2(
    scode,
    "vector[Mfm_1_1] zfm_frame_1_1;"
  )
  expect_match2(
    scode,
    paste0(
      "vector[Mfm_sigma_1_1] ",
      "zfm_frame_sigma_1_1;"
    )
  )
  helper_definitions <- regmatches(
    scode,
    gregexpr(
      "matrix fm_apply_reflector_brms(real tau, vector v, matrix B)",
      scode, fixed = TRUE
    )
  )
  expect_length(helper_definitions[[1]], 1L)
  expect_match2(scode, "sigma[n] += sdfm_sigma_1 * Cfm_sigma_1")
})

test_that("fm predictions reconstruct interaction and main effects", {
  dat <- expand.grid(
    user = factor(letters[1:4]),
    item = factor(LETTERS[1:3])
  )
  dat$y <- 0
  bterms <- brmsterms(y ~ fm(user, item, k = 2))
  bframe <- brmsframe(bterms, dat)$dpars$mu
  sdata <- standata(y ~ fm(user, item, k = 2), dat)

  frame1 <- cbind(
    c(1, -1, 0, 0) / sqrt(2),
    c(1, 1, -1, -1) / 2
  )
  frame2 <- cbind(
    c(1, -1, 0) / sqrt(2),
    c(1, 1, -2) / sqrt(6)
  )
  singular <- c(0.8, 0.6)
  main1 <- c(0.5, -0.5, 0, 0)
  main2 <- c(1, -1, 0)
  values <- c(
    2, as.vector(frame1), as.vector(frame2), singular,
    0.7, main1, 0.8, main2
  )
  variables <- c(
    "sdfm_user:item",
    fm_matrix_names("Qfm_1_1", 4, 2),
    fm_matrix_names("Qfm_1_2", 3, 2),
    sprintf("sifm_user:item[%d]", 1:2),
    "sdfm_main_user", sprintf("zfm_main_1_1[%d]", 1:4),
    "sdfm_main_item", sprintf("zfm_main_1_2[%d]", 1:3)
  )
  draws <- posterior::as_draws_matrix(matrix(
    values, nrow = 1, dimnames = list(NULL, variables)
  ))
  fm_prep <- prepare_predictions_fm(bframe, draws, sdata)
  prep <- list(fm = fm_prep, ndraws = 1L, nobs = nrow(dat))
  got <- as.numeric(predictor_fm(prep))

  expected <- 2 * sdata$Cfm_1 * rowSums(
    sweep(
      frame1[sdata$Jfm_1_1, , drop = FALSE] *
        frame2[sdata$Jfm_1_2, , drop = FALSE],
      2L, singular, "*"
    )
  )
  expected <- expected +
    0.7 * sdata$Cfm_main_1_1 * main1[sdata$Jfm_1_1] +
    0.8 * sdata$Cfm_main_1_2 * main2[sdata$Jfm_1_2]
  expect_equal(got, expected)
  expect_equal(
    as.numeric(predictor_fm(prep, i = c(2L, 8L))),
    expected[c(2L, 8L)]
  )

  no_main_terms <- brmsterms(
    y ~ fm(user, item, k = 2, main = FALSE)
  )
  no_main_frame <- brmsframe(no_main_terms, dat)$dpars$mu
  no_main_data <- standata(
    y ~ fm(user, item, k = 2, main = FALSE), dat
  )
  no_main_values <- c(
    2, as.vector(frame1), as.vector(frame2), singular
  )
  no_main_variables <- c(
    "sdfm_user:item",
    fm_matrix_names("Qfm_1_1", 4, 2),
    fm_matrix_names("Qfm_1_2", 3, 2),
    sprintf("sifm_user:item[%d]", 1:2)
  )
  no_main_draws <- posterior::as_draws_matrix(matrix(
    no_main_values, nrow = 1,
    dimnames = list(NULL, no_main_variables)
  ))
  no_main_prep <- list(
    fm = prepare_predictions_fm(
      no_main_frame, no_main_draws, no_main_data
    ),
    ndraws = 1L, nobs = nrow(dat)
  )
  expect_equal(
    as.numeric(predictor_fm(no_main_prep)),
    2 * no_main_data$Cfm_1 * rowSums(
      sweep(
        frame1[no_main_data$Jfm_1_1, , drop = FALSE] *
          frame2[no_main_data$Jfm_1_2, , drop = FALSE],
        2L, singular, "*"
      )
    )
  )

  surface <- frame1 %*% diag(singular) %*% t(frame2)
  expect_equal(rowSums(surface), rep(0, 4))
  expect_equal(colSums(surface), rep(0, 3))
  expect_equal(crossprod(frame1), diag(2), tolerance = 1e-14)
  expect_equal(crossprod(frame2), diag(2), tolerance = 1e-14)
  expect_equal(sum(singular^2), 1)
  expect_true(all(diff(singular) < 0))
  expect_lte(qr(surface)$rank, 2L)
  expect_equal(
    sqrt(mean((2 * sdata$Cfm_1 * surface)^2)), 2,
    tolerance = 1e-14
  )
})

test_that("distributional fm predictions retain prefixes on newdata", {
  train <- expand.grid(
    user = factor(letters[1:4]),
    item = factor(LETTERS[1:3])
  )
  train$y <- 0
  bterms <- brmsterms(bf(
    y ~ 1,
    sigma ~ fm(user, item, k = 2, main = FALSE)
  ))
  fit_frame <- brmsframe(bterms, train)
  basis <- frame_basis(fit_frame, train)
  newdata <- data.frame(
    y = 0,
    user = factor(c("d", "a"), levels = letters[1:4]),
    item = factor(c("C", "B"), levels = LETTERS[1:3])
  )
  new_frame <- brmsframe(bterms, newdata, basis = basis)$dpars$sigma
  sdata <- data_fm(new_frame, newdata)

  frame1 <- cbind(
    c(1, -1, 0, 0) / sqrt(2),
    c(1, 1, -1, -1) / 2
  )
  frame2 <- cbind(
    c(1, -1, 0) / sqrt(2),
    c(1, 1, -2) / sqrt(6)
  )
  singular <- c(0.8, 0.6)
  values <- c(1.5, as.vector(frame1), as.vector(frame2), singular)
  variables <- c(
    "sdfm_sigma_user:item",
    fm_matrix_names("Qfm_sigma_1_1", 4, 2),
    fm_matrix_names("Qfm_sigma_1_2", 3, 2),
    sprintf("sifm_sigma_user:item[%d]", 1:2)
  )
  draws <- posterior::as_draws_matrix(matrix(
    values, nrow = 1, dimnames = list(NULL, variables)
  ))
  prep <- list(
    fm = prepare_predictions_fm(new_frame, draws, sdata),
    ndraws = 1L, nobs = nrow(newdata)
  )
  expected <- 1.5 * sdata$Cfm_sigma_1 * rowSums(sweep(
    frame1[sdata$Jfm_sigma_1_1, , drop = FALSE] *
      frame2[sdata$Jfm_sigma_1_2, , drop = FALSE],
    2L, singular, "*"
  ))
  expect_equal(as.numeric(predictor_fm(prep)), expected)
})

test_that("fm exposes its invariant ordered spectrum under the term label", {
  dat <- fm_test_data()
  bframe <- brmsframe(
    brmsterms(y ~ fm(user, item, k = 2)), dat
  )$dpars$mu
  pars <- c(
    "sdfm_1", "fm_singular_1[1]", "fm_singular_1[2]",
    "Qfm_1_1[1,1]"
  )
  renamed <- rename_fm(bframe, pars)
  fnames <- unlist(lapply(renamed, `[[`, "fnames"), use.names = FALSE)
  expect_true("sdfm_user:item" %in% fnames)
  expect_true(all(sprintf("sifm_user:item[%d]", 1:2) %in% fnames))
})

test_that("fm prediction metadata permits new pairs but rejects new levels", {
  train <- expand.grid(
    user = factor(letters[1:4]),
    item = factor(LETTERS[1:3])
  )
  train <- train[-1L, , drop = FALSE]
  train$y <- 0
  bterms <- brmsterms(y ~ fm(user, item, k = 2))
  fit_frame <- brmsframe(bterms, train)
  basis <- frame_basis(fit_frame, train)

  known_pair <- data.frame(
    y = 0,
    user = factor("a", levels = letters[1:4]),
    item = factor("A", levels = LETTERS[1:3])
  )
  new_frame <- brmsframe(bterms, known_pair, basis = basis)
  new_data <- data_fm(new_frame$dpars$mu, known_pair)
  expect_equal(new_data$Jfm_1_1, as.array(1L))
  expect_equal(new_data$Jfm_1_2, as.array(1L))

  unseen <- data.frame(y = 0, user = "e", item = "A")
  unseen_frame <- brmsframe(bterms, unseen, basis = basis)
  expect_error(
    data_fm(unseen_frame$dpars$mu, unseen),
    "Unknown level.*e"
  )
})

test_that("fm compares CmdStan versions numerically", {
  expect_false(brms:::fm_cmdstan_version_supported("2.9.0"))
  expect_false(brms:::fm_cmdstan_version_supported("2.35.9"))
  expect_true(brms:::fm_cmdstan_version_supported("2.36.0"))
  expect_true(brms:::fm_cmdstan_version_supported("2.39.0"))
  expect_false(brms:::fm_cmdstan_version_supported(character(0)))
  expect_false(brms:::fm_cmdstan_version_supported(NA_character_))
})

test_that("fm enforces its Stan backend requirement", {
  dat <- fm_test_data()
  expect_error(
    stancode(
      y ~ fm(user, item, k = 2), dat,
      backend = "rstan", parse = FALSE
    ),
    "require 'backend = \"cmdstanr\"'"
  )

  skip_if_not_installed("cmdstanr")
  version <- try(cmdstanr::cmdstan_version(), silent = TRUE)
  skip_if(inherits(version, "try-error") || !length(version))
  skip_if(!brms:::fm_cmdstan_version_supported(version))
  expect_silent(stancode(
    bf(
      y ~ fm(user, item, k = 2),
      sigma ~ fm(user, item, k = 1, main = FALSE)
    ),
    dat, threads = threading(2), backend = "cmdstanr", parse = TRUE
  ))

  tiny <- expand.grid(
    user = factor(letters[1:2]),
    item = factor(LETTERS[1:2])
  )
  tiny$y <- c(0, 1, 1, 0)
  expect_silent(stancode(
    y ~ fm(user, item, k = 1, main = FALSE), tiny,
    backend = "cmdstanr", parse = TRUE
  ))
})

test_that("fm terms are not ignored by response-dependent S2Z Fisher rules", {
  dat <- fm_test_data()
  dat$y <- rep(0:1, length.out = nrow(dat))
  dat$group <- factor(rep(1:2, length.out = nrow(dat)))
  expect_error(
    stancode(
      y ~ fm(user, item, k = 1) +
        (1 | gr(group, s2z = TRUE, center = "fisher")),
      dat, family = bernoulli(), backend = "mock", parse = FALSE
    ),
    "factorization-machine terms"
  )
})
