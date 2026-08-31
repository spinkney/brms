# Reproduce the product-Wishart Monte Carlo calibration used to choose a
# lognormal approximation to the source model's centered full-table RMS.

set.seed(1)
S <- 200000L
N <- 100L
J <- 20L
K <- 5L

# Draw factors directly to reproduce the original RNG stream. Their centered
# cross-products are independently Wishart with N - 1 and J - 1 degrees of
# freedom. If H = A B', then ||H||_F^2 = tr((A'A)(B'B)).
rms <- numeric(S)
for (s in seq_len(S)) {
  left <- scale(
    matrix(rnorm(N * K), nrow = N, ncol = K),
    center = TRUE,
    scale = FALSE
  )
  right <- scale(
    matrix(rnorm(J * K), nrow = J, ncol = K),
    center = TRUE,
    scale = FALSE
  )
  rms[s] <- sqrt(sum(crossprod(left) * crossprod(right)) / (N * J))
}

# Moment-match a lognormal distribution on the RMS scale.
rms_mean <- mean(rms)
rms_var <- var(rms)
sigma_log <- sqrt(log1p(rms_var / rms_mean^2))
mean_log <- log(rms_mean) - sigma_log^2 / 2

expected <- c(meanlog = 0.767855432321696507, sdlog = 0.081032556096903921)
actual <- c(meanlog = mean_log, sdlog = sigma_log)
stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-12)))
print(c(mean_rms = rms_mean, sd_rms = sqrt(rms_var), actual))
