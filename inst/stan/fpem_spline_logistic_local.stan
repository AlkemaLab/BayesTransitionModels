functions {
#include ./fill_AR.stan
#include ./scale_blocks.stan
#include ./deboor.stan

  // From the Stan User's Guide 2.28: https://mc-stan.org/docs/2_28/stan-users-guide/truncated-random-number-generation.html.
  // Released under CC-BY 4.0 license (https://creativecommons.org/licenses/by/4.0/legalcode)
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }

  // From the Stan User's Guide 2.28: https://mc-stan.org/docs/2_28/stan-users-guide/truncated-random-number-generation.html.
  // Released under CC-BY 4.0 license (https://creativecommons.org/licenses/by/4.0/legalcode)
  real normal_lub_cdf(real x, real mu, real sigma, real lb, real ub) {
    real p_x  = normal_cdf(x  | mu, sigma);
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    return (p_x - p_lb) / p_ub;
  }

  real rate_spline(real P, real P_tilde, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
    return deboor(P / P_tilde, ext_knots, a, spline_degree);
  }
  
  real rate_logistic_logit(real P, real P_tilde, real omega) {
    if(P == 1 && P_tilde == 1) return 0;

    return (P < P_tilde) ? (
      (P - P_tilde) * omega / (P_tilde * (P - 1))
    ) : 0;
  }
}

data {
  int N; // Number of observations
  int T; // Number of time points
  int S; // Number of sources

  int<lower=1, upper=T> t_star;
  int<lower=1, upper=T> t_last;

  array[N] real y;                         // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=0, upper=S> source;   // Source of each observation
  array[N] int<lower=0, upper=1> held_out; // Whether to hold out each observation

  array[N] real<lower=0> s;                // Standard deviation

  int num_knots;
  vector[num_knots] knots;

  int spline_degree;

  int num_grid;
  vector[num_grid] grid;
  matrix[num_knots + spline_degree - 1, num_grid] B;
  vector[num_knots + spline_degree - 1] spline_maxima;

  real P_tilde_mu;
  real<lower=0> P_tilde_sigma;

  real Omega_mu;
  real<lower=0> Omega_sigma;

  real omega_mu;
  real<lower=0> omega_sigma;

  real<lower=0, upper=1> rho;
  real<lower=0> tau;

  // Data model
  array[S] real<lower=0> nonse;

  real a_lower_bound;
  real a_upper_bound;

  int<lower=0, upper=1> smoothing; // 0: no smoothing component. 1: AR(1) smoothing component.
}
transformed data {
  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  int num_constrained_zero = spline_degree + 1;

  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;
}

parameters {
  // P_tilde
  real P_tilde_raw;

  // Omega
  real Omega_raw;

  // omega
  real omega_raw;

  // Smoothing component
  vector[T * smoothing] epsilon_innovation;
}

transformed parameters {
  vector[T] eta;
  vector[num_basis] a;
  vector[T * smoothing] epsilon;

  real P_tilde_star;
  real<lower=0, upper=1> P_tilde;

  real Omega_star;
  real Omega;

  real omega_star;
  real omega;

  // P_tilde
  P_tilde_star = P_tilde_mu + P_tilde_raw * P_tilde_sigma;
  P_tilde = 0.5 + 0.45 * inv_logit(P_tilde_star);

  // Omega
  Omega = Omega_mu + Omega_raw * Omega_sigma;

  omega = 0.5 * inv_logit(omega_mu + omega_raw * omega_sigma);

  vector[T] logit_eta;
  real transition_function;

  for(i in 1:(num_basis - num_constrained_zero)) {
    a[i] = rate_logistic_logit(spline_maxima[i] * P_tilde, P_tilde, omega);
  }

  a[(num_basis - num_constrained_zero + 1):num_basis] = rep_vector(0, num_constrained_zero);

  if(smoothing == 1) {
    epsilon = to_vector(fill_AR(to_row_vector(epsilon_innovation), rho, tau, t_star));
  }

  logit_eta[t_star] = Omega;

  // Additive smoothing
  if(smoothing == 1) {
    for(t in (t_star + 1):T) {
      transition_function = rate_spline(inv_logit(logit_eta[t - 1]), P_tilde, to_row_vector(a), ext_knots, num_basis, spline_degree);
      logit_eta[t] = logit_eta[t - 1] + transition_function + epsilon[t];
    }

    for(q in 1:(t_star - 1)) {
      int t = t_star - q;
      transition_function = rate_spline(inv_logit(logit_eta[t + 1]), P_tilde, to_row_vector(a), ext_knots, num_basis, spline_degree);
      logit_eta[t] = logit_eta[t + 1] - transition_function - epsilon[t + 1];
    }
  }
  // No smoothing
  else if(smoothing == 0) {
    for(t in (t_star + 1):T) {
      transition_function = rate_spline(inv_logit(logit_eta[t - 1]), P_tilde, to_row_vector(a), ext_knots, num_basis, spline_degree);
      logit_eta[t] = logit_eta[t - 1] + transition_function;
    }

    for(q in 1:(t_star - 1)) {
      int t = t_star - q;
      transition_function = rate_spline(inv_logit(logit_eta[t + 1]), P_tilde, to_row_vector(a), ext_knots, num_basis, spline_degree);
      logit_eta[t] = logit_eta[t + 1] - transition_function;
    }
  }

  eta = inv_logit(logit_eta);
}

model {
  // P_tilde
  P_tilde_raw   ~ std_normal();

  // Omega
  Omega_raw   ~ std_normal();

  // omega
  omega_raw ~ std_normal();

  if(smoothing == 1) {
    to_vector(epsilon_innovation) ~ std_normal();
  }

  nonse ~ normal(0, 0.1);

  for(i in 1:N) {
    if(held_out[i] == 0) {
      y[i] ~ normal(eta[time[i]], s[i] + nonse[source[i]]) T[0, 1];
    }
  }
}
generated quantities {
  vector[N] y_pred;
  vector[N] resid;
  vector[N] pit;

  for(i in 1:N) {
    y_pred[i]  = normal_lub_rng(eta[time[i]], s[i] + nonse[source[i]], 0, 1);
    resid[i]   = y_pred[i] - y[i];
    pit[i]     = normal_lub_cdf(y[i] | eta[time[i]], s[i] + nonse[source[i]], 0, 1);
  }
}
