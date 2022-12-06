real deboor(real x, vector ext_knots, row_vector a, int degree) {
  int k = degree + 1;
  row_vector[degree + 1] d;
  int n_ext_knots = rows(ext_knots);

  if(x <= ext_knots[1]) return a[1];
  if(x >= ext_knots[n_ext_knots]) return a[cols(a)];

  while(!(ext_knots[k + 1] > x) && k < n_ext_knots - degree - 1) {
    k = k + 1;
  }

  d = a[(k - degree):(k)];

  for(r in 2:(degree + 1)) {
    for(j in (k + r - degree - 1):k) {
      int j2 = (k + r - degree - 1) - j + k;
      real alpha = (x - ext_knots[j2]) / (ext_knots[j2 + 1 + degree - (r - 1)] - ext_knots[j2]);

      d[j2 + 1 - k + degree] = (1 - alpha) * d[j2 - k + degree] + alpha * d[j2 + 1 - k + degree];
    }
  }

  return d[degree + 1];
}
