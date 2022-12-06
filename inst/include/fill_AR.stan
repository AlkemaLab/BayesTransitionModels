/* Fill out a partially filled AR process vector of length T
 *
 * @param epsilon_star A matrix which can be partially filled
 *        with an AR process.  The un-filled part of the
 *        matrix must be filled with N(0,1) draws.  Those
 *        draws are used as a source of randomness.  This function
 *        does not generate any random numbers itself.
 * @param rho the autocorrelation of the AR process.  No constraints
 *        are applied to this parameter although an exact value of 1
 *        will generate a 'divide by zero' error.
 * @param tau the standard deviation of the AR process.  No constraints
 *        are applied to this parameter although negative values will
 *        cause the AR process to flip at each step.
 * @param t_star an integer determining where the AR simulation will start.
 */
row_vector fill_AR(row_vector epsilon_star, real rho, real tau, int t_star) {
  int T = cols(epsilon_star);
  row_vector[T] epsilon = epsilon_star;

  epsilon[t_star] = epsilon[t_star] * tau / sqrt(1 - rho^2);

  for (q in 1:(t_star - 1)) {
    int t = t_star - q;
    epsilon[t] = rho * epsilon[t + 1] + epsilon[t] * tau;
  }
  for (t in (t_star + 1):T) {
    epsilon[t] = rho * epsilon[t - 1] + epsilon[t] * tau;
  }
  return epsilon;
}
