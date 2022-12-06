data {
  int<lower=0> N;
  vector[N] y;
  real scale_global;
}
transformed data {
  real<lower=1> nu_global = 1;
  real<lower=1> nu_local = 1;
  real<lower=0> slab_scale = 1;
  real<lower=0> slab_df = 1;
}
parameters {
  real<lower=0> sigma;
  vector[N] mu_raw;
  real<lower=0> global_shrinkage;
  vector<lower=0>[N] local_shrinkage; // called lambda in paper
  real<lower=0> caux;
}
transformed parameters {
  vector<lower=0>[N] truncated_local_shrinkage; // called lambda_tilde in paper
  real<lower=0> c_slab;
  vector[N] mu;

  c_slab = slab_scale * sqrt(caux);
  truncated_local_shrinkage = sqrt(c_slab^2 * square(local_shrinkage) ./ (c_slab^2 + global_shrinkage^2 * square(local_shrinkage)));
  mu = mu_raw .* truncated_local_shrinkage * global_shrinkage;
}
model {
  sigma ~ normal(0, 1);
  y ~ normal(mu, sigma);
}

