library(cmdstanr)
library(posterior)

chunk <- file.path("inst", "chunks", "fun_semiorthogonal_fm.stan")
if (!file.exists(chunk)) {
  stop("Run this script from the brms repository root.")
}

stan_code <- c(
  "functions {",
  "  #include fun_semiorthogonal_fm.stan",
  "}",
  "data {",
  "  int<lower=2> L;",
  "  int<lower=1, upper=L - 1> K;",
  "  int<lower=0, upper=1> special;",
  "}",
  "transformed data {",
  "  int P = fm_semiorthogonal_num_params_brms(L - 1, K, special);",
  "}",
  "generated quantities {",
  "  vector[P] y;",
  "  matrix[L - 1, K] R;",
  "  matrix[L, K] Q;",
  "  if (P > 0) {",
  "    for (p in 1:P) y[p] = std_normal_rng();",
  "  }",
  "  R = fm_semiorthogonal_constrain_brms(y, L - 1, K, special);",
  "  Q = fm_centered_semiorthogonal_constrain_brms(y, L, K, special);",
  "}"
)

stan_file <- tempfile(fileext = ".stan")
writeLines(stan_code, stan_file)
model <- cmdstan_model(
  stan_file,
  include_paths = normalizePath(dirname(chunk))
)

draw_frames <- function(L, K, special, draws = 1000L, seed = 1L) {
  fit <- model$sample(
    data = list(L = L, K = K, special = special),
    seed = seed,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 0,
    iter_sampling = draws,
    fixed_param = TRUE,
    sig_figs = 18,
    refresh = 0
  )
  list(
    R = as.matrix(fit$draws("R", format = "draws_matrix")),
    Q = as.matrix(fit$draws("Q", format = "draws_matrix"))
  )
}

check_geometry <- function(x, L, K, special) {
  orth_error <- center_error <- 0
  det_error <- numeric(nrow(x$R))
  for (s in seq_len(nrow(x$Q))) {
    R <- matrix(x$R[s, ], nrow = L - 1L, ncol = K)
    Q <- matrix(x$Q[s, ], nrow = L, ncol = K)
    orth_error <- max(
      orth_error,
      abs(crossprod(R) - diag(K)),
      abs(crossprod(Q) - diag(K))
    )
    center_error <- max(center_error, abs(colSums(Q)))
    if (special == 1L && K == L - 1L) {
      det_error[s] <- abs(det(R) - 1)
    }
  }
  stopifnot(orth_error < 1e-10, center_error < 1e-10)
  if (special == 1L && K == L - 1L) {
    stopifnot(max(det_error) < 1e-10)
  }
  invisible(c(orth_error = orth_error, center_error = center_error))
}

rectangular <- draw_frames(6L, 3L, 0L, seed = 601L)
check_geometry(rectangular, 6L, 3L, 0L)

# A Haar column in the centered subspace has E[Q[i,k]^2] = 1 / L.
q11 <- rectangular$Q[, "Q[1,1]"]
stopifnot(abs(mean(q11^2) - 1 / 6) < 0.025)

check_geometry(draw_frames(4L, 3L, 1L, seed = 403L), 4L, 3L, 1L)
check_geometry(draw_frames(5L, 4L, 1L, seed = 504L), 5L, 4L, 1L)
check_geometry(draw_frames(2L, 1L, 1L, draws = 10L, seed = 201L), 2L, 1L, 1L)

# The unrestricted 1 x 1 frame retains both determinant components. Its final
# scalar seed maps to a sign step; this checks support, not cross-mode mixing.
ordinary_edge <- draw_frames(2L, 1L, 0L, draws = 2000L, seed = 202L)
positive <- ordinary_edge$R[, "R[1,1]"] > 0
stopifnot(mean(positive) > 0.45, mean(positive) < 0.55)

message("Semi-orthogonal transform checks passed.")
