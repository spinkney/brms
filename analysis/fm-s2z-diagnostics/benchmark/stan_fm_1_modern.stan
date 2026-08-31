// Syntax-only modernization of adamlauretig/ny_r_talk/stan_fm_1.stan.
// The posterior model is unchanged. Posterior predictive draws are omitted
// so timing is comparable to the brms model, which computes them on demand.
data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> K;
  array[N * J, 2] int X;
  vector[N * J] y;
  real<lower=0> beta_sigma;
  real<lower=0> y_sigma;
}

parameters {
  vector[N] group_1_betas;
  vector[J] group_2_betas;
  matrix[N, K] gammas;
  matrix[J, K] deltas;
}

model {
  vector[N * J] linear_predictor;

  group_1_betas ~ normal(0, beta_sigma);
  group_2_betas ~ normal(0, beta_sigma);
  to_vector(gammas) ~ std_normal();
  to_vector(deltas) ~ std_normal();

  for (i in 1:(N * J)) {
    linear_predictor[i] =
      group_1_betas[X[i, 1]] + group_2_betas[X[i, 2]] +
      dot_product(gammas[X[i, 1]], deltas[X[i, 2]]);
  }
  y ~ normal(linear_predictor, y_sigma);
}
