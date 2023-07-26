functions {
#include ./fill_AR.stan
#include ./scale_blocks.stan
#include ./deboor.stan

  // https://mc-stan.org/docs/2_28/stan-users-guide/truncated-random-number-generation.html
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }

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
  int C; // Number of countries
  int S; // Number of sources

  array[C] int<lower=1, upper=T> t_star;
  int<lower=1, upper=T> t_last;

  array[N] real y;                         // Observations
  array[N] int<lower=1, upper=T> time;     // Time of each observation
  array[N] int<lower=1, upper=C> country;  // Country of each observation
  array[N] int<lower=0, upper=S> source;   // Source of each observation
  array[N] int<lower=0, upper=1> held_out;

  array[N] real<lower=0> s;                // Standard deviation

  int num_knots;
  vector[num_knots] knots;

  int spline_degree;

  int num_grid;
  vector[num_grid] grid;
  matrix[num_knots + spline_degree - 1, num_grid] B;
  vector[num_knots + spline_degree - 1] spline_maxima;

  // hierarchical P_tilde
  int<lower=0> P_tilde_n_terms;
  int<lower=0> P_tilde_n_re;
  array[P_tilde_n_re] int<lower=1, upper=P_tilde_n_terms> P_tilde_re_start;
  array[P_tilde_n_re] int<lower=1, upper=P_tilde_n_terms> P_tilde_re_end;
  matrix[C, P_tilde_n_terms] P_tilde_model_matrix;

  // hierarchical omega;
  int<lower=0> omega_n_terms;
  int<lower=0> omega_n_re;
  array[omega_n_re] int<lower=1, upper=omega_n_terms> omega_re_start;
  array[omega_n_re] int<lower=1, upper=omega_n_terms> omega_re_end;
  matrix[C, omega_n_terms] omega_model_matrix;

  // hierarchical Omega;
  int<lower=0> Omega_n_terms;
  int<lower=0> Omega_n_re;
  array[Omega_n_re] int<lower=1, upper=Omega_n_terms> Omega_re_start;
  array[Omega_n_re] int<lower=1, upper=Omega_n_terms> Omega_re_end;
  matrix[C, Omega_n_terms] Omega_model_matrix;


  // hierarchical a
  int<lower=0> a_n_terms;
  int<lower=0> a_n_re;
  array[a_n_re] int<lower=1, upper=a_n_terms> a_re_start;
  array[a_n_re] int<lower=1, upper=a_n_terms> a_re_end;
  matrix[C, a_n_terms] a_model_matrix;

  int<lower=0, upper=1> smoothing; // 0: no smoothing component. 1: AR(1) smoothing component.
}
transformed data {
  vector[rows(csr_extract_w(P_tilde_model_matrix))] P_tilde_model_matrix_w    = csr_extract_w(P_tilde_model_matrix);
  array[size(csr_extract_v(P_tilde_model_matrix))] int P_tilde_model_matrix_v = csr_extract_v(P_tilde_model_matrix);
  array[size(csr_extract_u(P_tilde_model_matrix))] int P_tilde_model_matrix_u = csr_extract_u(P_tilde_model_matrix);

  vector[rows(csr_extract_w(omega_model_matrix))] omega_model_matrix_w     = csr_extract_w(omega_model_matrix);
  array[size(csr_extract_v(omega_model_matrix))] int omega_model_matrix_v  = csr_extract_v(omega_model_matrix);
  array[size(csr_extract_u(omega_model_matrix))] int omega_model_matrix_u  = csr_extract_u(omega_model_matrix);

  vector[rows(csr_extract_w(Omega_model_matrix))] Omega_model_matrix_w     = csr_extract_w(Omega_model_matrix);
  array[size(csr_extract_v(Omega_model_matrix))] int Omega_model_matrix_v  = csr_extract_v(Omega_model_matrix);
  array[size(csr_extract_u(Omega_model_matrix))] int Omega_model_matrix_u  = csr_extract_u(Omega_model_matrix);

  vector[rows(csr_extract_w(a_model_matrix))] a_model_matrix_w             = csr_extract_w(a_model_matrix);
  array[size(csr_extract_v(a_model_matrix))] int a_model_matrix_v          = csr_extract_v(a_model_matrix);
  array[size(csr_extract_u(a_model_matrix))] int a_model_matrix_u          = csr_extract_u(a_model_matrix);


  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  int num_constrained_zero = spline_degree + 1;

  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;
}

parameters {
  // P_tilde
  vector[P_tilde_n_terms] P_tilde_raw;
  vector<lower=0>[P_tilde_n_re - 1] P_tilde_sigma_raw;
  //real P_tilde_mu;

  // omega
  vector[omega_n_terms] omega_raw;
  vector<lower=0>[omega_n_re - 1] omega_sigma_raw;

  // Omega
  vector[Omega_n_terms] Omega_raw;
  vector<lower=0>[Omega_n_re - 1] Omega_sigma_raw;
  //real Omega_mu;

  // Smoothing component
  //matrix<lower=-5, upper=5>[C * smoothing, T * smoothing] epsilon_innovation;
  matrix[C * smoothing, T * smoothing] epsilon_innovation;
  array[smoothing] real<lower=0, upper=1> est_rho;
  array[smoothing] real<lower=0> est_tau;

  // Data model
  array[S] real<lower=0> nonse;
}

transformed parameters {
  matrix[C, T] eta;
  matrix[C, num_basis] a;
  array[smoothing] real rho;
  array[smoothing] real tau;
  matrix[C * smoothing, T * smoothing] epsilon;
	vector[N] scale;

  vector<lower=0, upper=1>[C] P_tilde;
  vector[P_tilde_n_terms] P_tilde_star;

  vector[C] omega;
  vector[omega_n_terms] omega_star;

  vector[C] Omega;
  vector[Omega_n_terms] Omega_star;

  vector<lower=0>[P_tilde_n_re] P_tilde_sigma;
  vector<lower=0>[Omega_n_re] Omega_sigma;
  vector<lower=0>[omega_n_re] omega_sigma;

  omega_sigma[1] = 5;
  if(omega_n_re > 1) omega_sigma[2:omega_n_re] = omega_sigma_raw;

  Omega_sigma[1] = 3;
  if(Omega_n_re > 1) Omega_sigma[2:Omega_n_re] = Omega_sigma_raw;

  P_tilde_sigma[1] = 3;
  if(P_tilde_n_re > 1) P_tilde_sigma[2:P_tilde_n_re] = P_tilde_sigma_raw;

  if(smoothing == 1) {
    rho[1] = est_rho[1];
    tau[1] = est_tau[1];
  }

  // P_tilde
  P_tilde_star = scale_blocks(P_tilde_raw, P_tilde_sigma, P_tilde_re_start, P_tilde_re_end);
  P_tilde = 0.5 + 0.45 * inv_logit(csr_matrix_times_vector(C, P_tilde_n_terms, P_tilde_model_matrix_w, P_tilde_model_matrix_v, P_tilde_model_matrix_u, P_tilde_star));

  // omega
  omega_star = scale_blocks(omega_raw, omega_sigma, omega_re_start, omega_re_end);
  omega = 0.5 * inv_logit(csr_matrix_times_vector(C, omega_n_terms, omega_model_matrix_w, omega_model_matrix_v, omega_model_matrix_u, omega_star));

  // Omega
  Omega_star = scale_blocks(Omega_raw, Omega_sigma, Omega_re_start, Omega_re_end);
  Omega = csr_matrix_times_vector(C, Omega_n_terms, Omega_model_matrix_w, Omega_model_matrix_v, Omega_model_matrix_u, Omega_star);

  // Initialize the non-zero spline coefficients

  for(c in 1:C) {
    row_vector[T] logit_eta;
    real transition_function;

    for(i in 1:(num_basis - num_constrained_zero)) {
      a[c, i] = rate_logistic_logit(spline_maxima[i] * P_tilde[c], P_tilde[c], omega[c]);
    }

    a[c, (num_basis - num_constrained_zero + 1):num_basis] = rep_row_vector(0, num_constrained_zero);

    if(smoothing == 1) {
      epsilon[c,] = fill_AR(epsilon_innovation[c, ], rho[1], tau[1], t_star[c]);
    }

    logit_eta[t_star[c]] = Omega[c];

    // Additive smoothing
    if(smoothing == 1) {
      for(t in (t_star[c] + 1):T) {
        transition_function = rate_spline(inv_logit(logit_eta[t - 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t - 1] + transition_function + epsilon[c, t];
      }

      for(q in 1:(t_star[c] - 1)) {
        int t = t_star[c] - q;
        transition_function = rate_spline(inv_logit(logit_eta[t + 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t + 1] - transition_function - epsilon[c, t + 1];
      }
    }
    // No smoothing
    else if(smoothing == 0) {
      for(t in (t_star[c] + 1):T) {
        transition_function = rate_spline(inv_logit(logit_eta[t - 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t - 1] + transition_function;
      }

      for(q in 1:(t_star[c] - 1)) {
        int t = t_star[c] - q;
        transition_function = rate_spline(inv_logit(logit_eta[t + 1]), P_tilde[c], a[c,], ext_knots, num_basis, spline_degree);
        logit_eta[t] = logit_eta[t + 1] - transition_function;
      }
    }

    eta[c, ] = inv_logit(logit_eta);
  }
	for(i in 1:N) {
    scale[i] = sqrt(square(s[i]) + square(nonse[source[i]]));
  }
}

model {
  // P_tilde
  //P_tilde_mu    ~ std_normal();
	//for(i in 1:(P_tilde_n_re - 1)) {
  //	P_tilde_sigma_raw[i] ~ normal(0, 3) T[0, positive_infinity()];
	//}
	P_tilde_sigma_raw ~ normal(0, 3);
  P_tilde_raw   ~ std_normal();

  // omega
	//for(i in 1:(omega_n_re - 1)) {
  //	omega_sigma_raw[i] ~ normal(0, 5) T[0, positive_infinity()];
	//}
	omega_sigma_raw ~ normal(0, 5);
  omega_raw   ~ std_normal();

  // Omega
  //Omega_mu    ~ std_normal();
	//for(i in 1:(Omega_n_re - 1)) {
  //	Omega_sigma_raw[i] ~ normal(0, 3) T[0, positive_infinity()];
	//}
	Omega_sigma_raw ~ normal(0, 3);
  Omega_raw   ~ std_normal();

  if(smoothing == 1) {
    est_rho[1] ~ {{RHO_PRIOR}} T[0, 1];
    est_tau[1] ~ {{TAU_PRIOR}} T[0, positive_infinity()];
    to_vector(epsilon_innovation) ~ std_normal();
  }

	//for(i in 1:S) {
  //	nonse[i] ~ normal(0, 0.1) T[0, positive_infinity()];
	//}
	nonse ~ normal(0, 0.1);

  for(i in 1:N) {
    if(held_out[i] == 0) {
      y[i] ~ normal(eta[country[i], time[i]], scale[i]) T[0, 1];
    }
  }
}
generated quantities {
  vector[N] y_pred;
  vector[N] resid;
  vector[N] pit;

  for(i in 1:N) {
    y_pred[i]  = normal_lub_rng(eta[country[i], time[i]], scale[i], 0, 1);
    resid[i]   = y_pred[i] - y[i];
    pit[i]     = normal_lub_cdf(y[i] | eta[country[i], time[i]], scale[i], 0, 1);
  }
}
