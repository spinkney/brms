data {
  int<lower=2> N1;
  int<lower=2> N2;
  int<lower=1> K;
  int<lower=1> Nobs;
  array[Nobs] int<lower=1, upper=N1> I1;
  array[Nobs] int<lower=1, upper=N2> I2;
  vector[Nobs] y;
  real<lower=0> y_sd;
  real<lower=0> intercept_sd;
  real<lower=0> main_sd;
  real<lower=0> factor_sd;
  real<lower=0> C1;
  real<lower=0> C2;
}

parameters {
  real Intercept;
  sum_to_zero_vector[N1] zmain1;
  sum_to_zero_vector[N2] zmain2;
  array[K] sum_to_zero_vector[N1] left;
  array[K] sum_to_zero_vector[N2] right;
}

model {
  Intercept ~ normal(0, intercept_sd);
  target += std_normal_lupdf(zmain1);
  target += std_normal_lupdf(zmain2);
  for (k in 1:K) {
    target += std_normal_lupdf(left[k]);
    target += std_normal_lupdf(right[k]);
  }

  for (n in 1:Nobs) {
    real interaction = 0;
    for (k in 1:K) {
      interaction += left[k][I1[n]] * right[k][I2[n]];
    }
    y[n] ~ normal(
      Intercept
      + main_sd * C1 * zmain1[I1[n]]
      + main_sd * C2 * zmain2[I2[n]]
      + square(factor_sd) * C1 * C2 * interaction,
      y_sd
    );
  }
}
