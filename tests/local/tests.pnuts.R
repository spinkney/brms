# Opt-in integration test; requires PNUTS, Python BridgeStan 2.9.x and CmdStan.
# Configure brms.pnuts.* options or PNUTS_EXECUTABLE / PNUTS_PYTHON / BRIDGESTAN.
# PNUTS_TEST_OUTPUT may name a fresh directory to preserve results.
library(brms)
library(testthat)
output <- Sys.getenv("PNUTS_TEST_OUTPUT", tempfile("brms-pnuts-validation-"))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(mc.cores = 4)
fit_pnuts <- function(label, ...) {
  fit <- brm(..., engine = "pnuts", chains = 4, cores = 4,
             iter = 2000, warmup = 1000, seed = 1921,
             output_dir = file.path(output, label), save_warmup = TRUE)
  saveRDS(fit, file.path(output, paste0(label, ".rds")))
  fit
}
check_fit <- function(fit, nobs) {
  expect_identical(fit$backend, "pnuts")
  expect_equal(dim(pnuts_diagnostics(fit))[1:2], c(1000, 4))
  expect_equal(dim(pnuts_diagnostics(fit, TRUE))[1:2], c(2000, 4))
  expect_equal(dim(posterior_predict(fit, ndraws = 20)), c(20, nobs))
  expect_equal(dim(posterior_epred(fit, ndraws = 20)), c(20, nobs))
  expect_equal(dim(log_lik(fit, ndraws = 20)), c(20, nobs))
  stats <- posterior::summarise_draws(as_draws_array(fit))
  expect_lt(max(stats$rhat), 1.01)
  expect_gt(min(stats$ess_bulk), 400)
  print(summary(fit))
  stats
}

set.seed(481)
dat <- data.frame(x = rnorm(120))
dat$y <- 1.2 - .7 * dat$x + rnorm(120, sd = .8)
gaussian <- fit_pnuts("gaussian", y ~ x, data = dat)
gaussian_stats <- check_fit(gaussian, nrow(dat))
expect_s3_class(loo(gaussian, cores = 1), "loo")
restored <- readRDS(file.path(output, "gaussian.rds"))
expect_equal(as_draws_array(restored), as_draws_array(gaussian))

# Same target through the PR's existing CmdStan converter.
stan <- brm(y ~ x, data = dat, backend = "cmdstanr", chains = 4, cores = 4,
            iter = 2000, warmup = 1000, seed = 7121, refresh = 0,
            control = list(metric = "dense_e"))
saveRDS(stan, file.path(output, "gaussian-stan.rds"))
compare_stats <- function(fit) {
  posterior::summarise_draws(as_draws_array(fit), "mean", "mcse_mean", "rhat", "ess_bulk")
}
a <- compare_stats(gaussian)
b <- compare_stats(stan)
comparison <- merge(as.data.frame(a), as.data.frame(b), by = "variable", suffixes = c("_pnuts", "_stan"))
comparison <- subset(comparison, variable != "lp__")
comparison$combined_mcse_z <- abs(comparison$mean_pnuts - comparison$mean_stan) /
  sqrt(comparison$mcse_mean_pnuts^2 + comparison$mcse_mean_stan^2)
expect_lt(max(comparison$combined_mcse_z), 4)
write.csv(comparison, file.path(output, "gaussian-comparison.csv"), row.names = FALSE)
print(comparison)

set.seed(482)
multi <- data.frame(x = rnorm(240), group = factor(rep(1:12, each = 20)))
multi$y <- .4 + .8 * multi$x + rep(rnorm(12, sd = .9), each = 20) + rnorm(240, sd = .6)
multilevel <- fit_pnuts("multilevel", y ~ x + (1 | group), data = multi,
                        prior = prior(normal(0, 1), class = sd))
multi_stats <- check_fit(multilevel, nrow(multi))
expect_equal(dim(ranef(multilevel)$group), c(12, 4, 1))

set.seed(483)
binary <- data.frame(x = rnorm(200))
binary$y <- rbinom(200, 1, plogis(-.4 + .8 * binary$x))
bernoulli_fit <- fit_pnuts("bernoulli", y ~ x, data = binary, family = bernoulli(),
                           prior = prior(normal(0, 2), class = b),
                           init = function(chain_id) list(b = 0, Intercept = 0))
binary_stats <- check_fit(bernoulli_fit, nrow(binary))
expect_true(all(posterior_predict(bernoulli_fit, ndraws = 10) %in% 0:1))

# Reuse a compiled model, with fresh output paths and thinning across saved warmup.
updated_data <- dat
updated_data$y <- updated_data$y + .3
updated <- update(gaussian, newdata = updated_data, chains = 2, cores = 2,
                  iter = 1200, warmup = 1000, thin = 2, seed = 901,
                  output_dir = file.path(output, "update"))
expect_equal(dim(as_draws_array(updated))[1:2], c(100, 2))
expect_gt(fixef(updated)["Intercept", "Estimate"] - fixef(gaussian)["Intercept", "Estimate"], .2)
replayed <- update(gaussian, newdata = updated_data, chains = 2, cores = 2,
                   iter = 1200, warmup = 1000, thin = 2, seed = 901)
expect_equal(as_draws_array(replayed), as_draws_array(updated), tolerance = 0)
expect_equal(pnuts_diagnostics(replayed), pnuts_diagnostics(updated), tolerance = 0)
expect_false(identical(attr(replayed$fit, "pnuts")$files, attr(updated$fit, "pnuts")$files))

# Validation refuses unsupported algorithms and never overwrites an earlier run.
expect_error(update(gaussian, algorithm = "meanfield", iter = 10), "algorithm = 'sampling' only")
old_data <- readBin(file.path(output, "gaussian", "data.json"), "raw", n = 1e7)
expect_error(update(gaussian, newdata = updated_data, output_dir = file.path(output, "gaussian")), "already contains a run")
expect_identical(readBin(file.path(output, "gaussian", "data.json"), "raw", n = 1e7), old_data)

# Switch an existing CmdStan fit through the public engine alias.
switched <- update(stan, engine = "pnuts", chains = 2, cores = 2,
                    iter = 1200, warmup = 1000, seed = 901,
                    output_dir = file.path(output, "switch"))
expect_identical(switched$backend, "pnuts")
expect_equal(dim(as_draws_array(switched))[1:2], c(200, 2))

# A real subprocess timeout invalidates the whole fit and preserves its logs.
expect_error(update(gaussian, chains = 2, cores = 2,
                     control = list(chain_timeout = .001),
                     output_dir = file.path(output, "timeout")), "fit failed; no partial ensemble")
failed <- jsonlite::read_json(file.path(output, "timeout", "run.json"), simplifyVector = TRUE)
expect_false(failed$completed)
expect_true(any(failed$statuses$timed_out))

saveRDS(list(gaussian = gaussian_stats, multilevel = multi_stats, bernoulli = binary_stats,
             comparison = comparison), file.path(output, "validation-summary.rds"))
cat("PNUTS integration tests passed. Results:", output, "\n")
