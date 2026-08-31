functions {
  #include fun_semiorthogonal_fm.stan
}

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
  real<lower=0> C1;
  real<lower=0> C2;
  real<lower=0> Cinteraction;
  int<lower=0, upper=1> special1;
  int<lower=0, upper=1> special2;
  int<lower=0> M1;
  int<lower=0> M2;
  real rms_meanlog;
  real<lower=0> rms_sdlog;
}

parameters {
  real Intercept;
  sum_to_zero_vector[N1] zmain1;
  sum_to_zero_vector[N2] zmain2;
  vector[M1] zframe1;
  vector[M2] zframe2;
  simplex[K] spectrum_gaps;
  real<lower=0> sdfm;
}

transformed parameters {
  matrix[N1, K] Q1 = fm_centered_semiorthogonal_constrain_brms(
    zframe1, N1, K, special1
  );
  matrix[N2, K] Q2 = fm_centered_semiorthogonal_constrain_brms(
    zframe2, N2, K, special2
  );
  vector[K] singular;
  for (r in 1:K) {
    real energy = 0;
    for (j in r:K) {
      energy += spectrum_gaps[j] / j;
    }
    singular[r] = sqrt(energy);
  }
}

model {
  Intercept ~ normal(0, intercept_sd);
  target += std_normal_lupdf(zmain1);
  target += std_normal_lupdf(zmain2);
  target += std_normal_lupdf(zframe1);
  target += std_normal_lupdf(zframe2);
  target += dirichlet_lupdf(spectrum_gaps | rep_vector(1, K));
  sdfm ~ lognormal(rms_meanlog, rms_sdlog);

  for (n in 1:Nobs) {
    y[n] ~ normal(
      Intercept
      + main_sd * C1 * zmain1[I1[n]]
      + main_sd * C2 * zmain2[I2[n]]
      + sdfm * Cinteraction * dot_product(
          Q1[I1[n]] .* to_row_vector(singular), Q2[I2[n]]
        ),
      y_sd
    );
  }
}
