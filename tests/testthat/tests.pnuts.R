context("Native PNUTS backend")

test_that("external error output is preserved literally", {
  expect_error(.pnuts_stop('Python traceback: {"compile": compile_model}'),
               'Python traceback: {"compile": compile_model}', fixed = TRUE)
})

pnuts_fixture <- function(chains = 1L, save_warmup = FALSE, thin = 1L) {
  directory <- tempfile("pnuts-fixture-")
  dir.create(directory)
  files <- file.path(directory, paste0("chain-", seq_len(chains), ".csv"))
  for (i in seq_along(files)) {
    x <- data.frame(iteration__ = 1:7, warmup__ = c(1, 1, rep(0, 5)),
                    lp__ = -1:-7, stepsize__ = c(.1, .2, rep(.3, 5)),
                    adapting__ = c(1, 1, rep(0, 5)), metric_updated__ = 0,
                    probe_grad_evals__ = c(8, rep(0, 6)),
                    energy_fidelity__ = .95, numerical_failures__ = 0,
                    b.1 = i * 100 + 1:7, b.2 = i * 100 + 11:17,
                    sigma = 1, omit.1 = 0, b_Intercept = i)
    if (!save_warmup) x <- x[x$warmup__ == 0, ]
    write.table(x, files[i], sep = ",", row.names = FALSE, quote = FALSE)
    cat("# seconds_total=1\n", file = files[i], append = TRUE)
  }
  list(files = files, model = list(library = "unused.so"), exclude = "omit",
       iter = 7L, warmup = 2L, thin = thin, seed = 321L,
       control = .pnuts_control(NULL),
       run = list(statuses = data.frame(seconds = rep(1, chains))),
       save_warmup = save_warmup)
}

test_that("engine selects PNUTS without changing unrelated backends", {
  dat <- data.frame(y = 1:5)
  expect_identical(brm(y ~ 1, dat, engine = "pnuts", empty = TRUE, cores = 1)$backend, "pnuts")
  expect_identical(brm(y ~ 1, dat, backend = "pnuts", empty = TRUE, cores = 1)$backend, "pnuts")
  expect_identical(brm(y ~ 1, dat, backend = "mock", engine = "walnuts", empty = TRUE, cores = 1)$stan_args$engine, "walnuts")
  expect_error(brm(y ~ 1, dat, backend = "pnuts", engine = "other", empty = TRUE), "requires engine")
  cached <- tempfile(fileext = ".rds")
  saveRDS(brm(y ~ 1, dat, empty = TRUE, cores = 1), cached)
  expect_error(brm(y ~ 1, dat, engine = "pnuts", file = cached), "cached fit was not sampled with PNUTS")
  expect_error(brm(y ~ 1, dat, backend = "pnuts", file = cached), "cached fit was not sampled with PNUTS")
})

test_that("PNUTS controls reject misleading or malformed options", {
  expect_identical(.pnuts_control(NULL)$geometry_frame, "pilot")
  expect_error(.pnuts_control(list(adapt_delta = .9)), "Unknown PNUTS control")
  expect_error(.pnuts_control(list(max_treedepth = 10)), "Unknown PNUTS control")
  expect_error(.pnuts_control(list(adapt = c("diag", "none"))), "single string")
  expect_error(.pnuts_control(list(adapt = NULL)), "single string")
  expect_error(.pnuts_control(list(gamma = NA)), "finite number")
  expect_error(.pnuts_control(list(chain_timeout = 0)), "positive")
  expect_error(.pnuts_control(list(gamma = 1, gamma = 2)), "uniquely named")
  expect_error(.pnuts_integer(1.5, "thin", 1), "integer")
  expect_error(.pnuts_integer(NA, "seed"), "integer")
})

test_that("PNUTS initializes chains independently and retains partial lists", {
  expect_equal(.pnuts_initializations("0", 2), rep(list(list(kind = "zero", radius = 0)), 2))
  expect_length(.pnuts_initializations(NULL, 3), 3)
  expect_equal(.pnuts_initializations(function(chain_id) list(sigma = chain_id), 2)[[2]]$values$sigma, 2)
  expect_identical(.pnuts_initializations(list(list(b = 0)), 1)[[1]]$kind, "constrained")
  expect_error(.pnuts_initializations(list(list(0)), 1), "unique names")
  expect_error(.pnuts_initializations(c(1, 2), 2), "init must be")
})

test_that("native conversion preserves constrained draws, diagnostics and thinning", {
  for (chains in c(1L, 2L, 5L)) {
    for (saved in c(FALSE, TRUE)) {
      input <- pnuts_fixture(chains, saved, thin = 2L)
      fit <- do.call(.pnuts_read_fit, input)
      draws <- as.array(fit)
      expect_equal(dim(draws), c(3L, chains, 5L))
      expect_equal(as.numeric(draws[, 1, "b[1]"]), c(103, 105, 107))
      expect_true("b_Intercept" %in% fit@sim$fnames_oi)
      expect_false("omit[1]" %in% fit@sim$fnames_oi)
      expect_equal(fit@par_dims$b, 2L)
      expect_length(fit@stan_args, chains)
      expect_equal(fit@sim$warmup2, rep(as.integer(saved), chains))
      expect_identical(fit@stan_args[[chains]]$sampler_t, "PNUTS")
      info <- attr(fit, "pnuts")
      expect_equal(dim(info$diagnostics), c(3L, chains, 9L))
      expect_false(any(c("accept_stat__", "divergent__") %in% posterior::variables(info$diagnostics)))
      object <- brm(y ~ 1, data.frame(y = 1:5), engine = "pnuts", empty = TRUE, cores = 1)
      object$fit <- fit
      expect_equal(dim(pnuts_diagnostics(object, inc_warmup = TRUE))[1], 3L + as.integer(saved))
      expect_error(nuts_params(object), "pnuts_diagnostics")
      expect_equal(control_params(object), input$control)
      expect_equal(elapsed_time(object)$total, rep(1, chains))
      expect_error(combine_models(object, object), "Combining PNUTS")
    }
  }
})

test_that("incomplete or still-adapting ensembles cannot become brms fits", {
  input <- pnuts_fixture()
  text <- readLines(input$files)
  writeLines(text[-length(text)], input$files)
  expect_error(do.call(.pnuts_read_fit, input), "completion record")
  input$files <- paste0(input$files, ".partial")
  expect_error(do.call(.pnuts_read_fit, input), "Refusing partial")

  input <- pnuts_fixture()
  x <- read.csv(input$files, comment.char = "#", check.names = FALSE)
  x$adapting__[1] <- 1
  write.table(x, input$files, sep = ",", row.names = FALSE, quote = FALSE)
  cat("# seconds_total=1\n", file = input$files, append = TRUE)
  expect_error(do.call(.pnuts_read_fit, input), "frozen tuning")

  input <- pnuts_fixture(chains = 2)
  writeLines("# FAILED: aborted\n", input$files[2])
  expect_error(do.call(.pnuts_read_fit, input), "completion record")
})
