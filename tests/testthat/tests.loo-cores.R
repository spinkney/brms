context("LOO post-processing core defaults")

test_that("LOO and PSIS defaults do not fork when mc.cores is greater than one", {
  # Emulate the Positron Console's restriction without requiring the IDE.
  withr::local_options(mc.cores = 4L, loo.cores = NULL)
  testthat::local_mocked_bindings(
    mcfork = function(...) stop("Forking is disabled in this test"),
    .package = "parallel"
  )
  fit <- rename_pars(brmsfit_example3)
  reference <- suppressWarnings(loo(fit, cores = 1))
  automatic <- suppressWarnings(loo(fit))
  pointwise <- suppressWarnings(loo(fit, pointwise = TRUE))
  expect_equal(automatic$estimates, reference$estimates)
  expect_equal(automatic$pointwise, reference$pointwise)
  expect_equal(pointwise$estimates, reference$estimates)
  expect_equal(pointwise$pointwise, reference$pointwise)

  weights <- suppressWarnings(psis(fit))
  weights_reference <- suppressWarnings(psis(fit, cores = 1))
  expect_equal(weights, weights_reference)

  stored <- suppressWarnings(add_criterion(fit, "loo"))
  expect_equal(stored$criteria$loo$estimates, reference$estimates)
  expect_identical(getOption("mc.cores"), 4L)
  expect_null(getOption("loo.cores"))
})

test_that("explicit cores reach all three post-processing stages", {
  withr::local_options(mc.cores = 8L)
  seen <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    log_lik = function(object, cores, ...) {
      seen$log_lik <- cores
      matrix(0, 10, 2)
    },
    r_eff_log_lik = function(x, cores, ...) {
      seen$r_eff <- cores
      c(1, 1)
    }
  )
  explicit <- prepare_loo_args(list(), NULL, NULL, FALSE, cores = 2L)
  expect_identical(seen$log_lik, 2L)
  expect_identical(seen$r_eff, 2L)
  expect_identical(explicit$cores, 2L)

  automatic <- prepare_loo_args(list(), NULL, NULL, FALSE)
  expect_identical(seen$log_lik, 1L)
  expect_identical(seen$r_eff, 1L)
  expect_identical(automatic$cores, 1L)
  expect_identical(getOption("mc.cores"), 8L)
})
