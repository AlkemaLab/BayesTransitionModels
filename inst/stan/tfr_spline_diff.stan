functions {
  #include ./scale_blocks.stan
  #include ./deboor.stan

  real rate_spline(real P, real P_tilde, real P_tilde2, row_vector a, vector ext_knots, int num_basis, int spline_degree) {
    return deboor((P - P_tilde) / P_tilde2, ext_knots, a, spline_degree);
  }
}

data {
  int N; // Number of observations
  int T; // Number of time points
  int C; // Number of countries

  int<lower=1, upper=T> start_phase2[C];
  int<lower=1, upper=T> end_phase2[C];

  real y[N];                         // Observations
  int<lower=1, upper=T> time[N];     // Time of each observation
  int<lower=1, upper=C> country[N];  // Country of each observation
  int<lower=0, upper=1> held_out[N];
  int t_last;

  int num_knots;
  vector[num_knots] knots;

  int spline_degree;

  int num_grid;
  vector[num_grid] grid;
  matrix[num_knots + spline_degree - 1, num_grid] B;

  // hierarchical a
  int<lower=0> a_n_terms;
  int<lower=0> a_n_re;
  int<lower=1, upper=a_n_terms> a_re_start[a_n_re];
  int<lower=1, upper=a_n_terms> a_re_end[a_n_re];
  matrix[C, a_n_terms] a_model_matrix;

  real a_lower_bound;
  real a_upper_bound;
}
transformed data {
  vector[rows(csr_extract_w(a_model_matrix))] a_model_matrix_w             = csr_extract_w(a_model_matrix);
  int a_model_matrix_v[size(csr_extract_v(a_model_matrix))]                = csr_extract_v(a_model_matrix);
  int a_model_matrix_u[size(csr_extract_u(a_model_matrix))]                = csr_extract_u(a_model_matrix);

  int num_basis = num_knots + spline_degree - 1;
  vector[2 * spline_degree + num_knots] ext_knots;
  int num_constrained_zero = spline_degree + 1;

  int num_epsilon = 0;
  matrix[C, t_last] ymat = rep_matrix(0, C, t_last);


  vector[C] Omega;
  for(c in 1:C) {
    for(i in 1:N) {
      if(country[i] == c && time[i] == start_phase2[c]) {
        Omega[c] = y[i];
      }
    }
  }

  for(i in 1:N) {
    ymat[country[i], time[i]] = y[i];
  }

  ext_knots[1:spline_degree] = rep_vector(knots[1], spline_degree);
  ext_knots[(num_knots + spline_degree + 1):(num_knots + 2 * spline_degree)] = rep_vector(knots[num_knots], spline_degree);
  ext_knots[(spline_degree + 1):(num_knots + spline_degree)] = knots;
}

parameters {
  // Spline rate vs. level function
  matrix[a_n_terms, num_basis - num_constrained_zero] a_raw;
  matrix<lower=0>[a_n_re, num_basis - num_constrained_zero] a_sigma;

  real<lower=0> epsilon_scale;
  //real<lower=0> sigma;
}

transformed parameters {
  matrix[C, T] transition_function = rep_matrix(0, C, T);
  matrix[C, num_basis] a;
  vector[C] P_tilde;
  vector[C] P_tilde2;
  real<lower=0> sigma = 0.25;

  matrix[a_n_terms, num_basis - num_constrained_zero] a_star;

  P_tilde = rep_vector(1, C);

  // Initialize the non-zero spline coefficients
  for(i in 1:(num_basis - num_constrained_zero)) {
    a_star[, i] = scale_blocks(a_raw[, i], a_sigma[, i], a_re_start, a_re_end);
    a[, num_constrained_zero + i] = a_lower_bound + (a_upper_bound - a_lower_bound) * inv_logit(csr_matrix_times_vector(C, a_n_terms, a_model_matrix_w, a_model_matrix_v, a_model_matrix_u, a_star[, i]));
  }

  for(c in 1:C) {
    a[c, 1:num_constrained_zero] = rep_row_vector(0, num_constrained_zero);
    P_tilde2[c] = Omega[c] - P_tilde[c];

    for(t in (start_phase2[c] + 1):end_phase2[c]) {
      transition_function[c, t] = rate_spline(ymat[c, t - 1], P_tilde[c], P_tilde2[c], a[c,], ext_knots, num_basis, spline_degree);
    }
  }
}

model {
  to_vector(a_raw) ~ std_normal();
  to_vector(a_sigma) ~ normal(0, 0.5);

  epsilon_scale ~ std_normal();

  for(i in 1:N) {
    if(held_out[i] == 0 && time[i] > start_phase2[country[i]]) {
      ymat[country[i], time[i]] - ymat[country[i], time[i] - 1] ~ normal(transition_function[country[i], time[i]], epsilon_scale);
    }
  }
}
generated quantities {
  matrix[C, T] eta = rep_matrix(0, C, T);

  for(c in 1:C) {
    eta[c, 1:t_last] = ymat[c, 1:t_last];
    for(t in (end_phase2[c] + 1):T) {
      real error = normal_rng(0, epsilon_scale);
      real transition = rate_spline(eta[c, t - 1], P_tilde[c], P_tilde2[c], a[c,], ext_knots, num_basis, spline_degree);
      eta[c, t] = eta[c, t - 1] + transition + error;
    }
  }
}
