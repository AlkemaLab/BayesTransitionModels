functions {
  #include ./deboor.stan
}
data {
  int num_knots;
  vector[num_knots] knots;
  int spline_degree;
  int num_grid;
  vector[num_grid] grid;
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
  row_vector[num_knots + spline_degree - 1] spline_coefs;
}
transformed parameters {
  vector[num_grid] y;
  for(n in 1:num_grid) {
    y[n] = deboor(grid[n], ext_knots, spline_coefs, spline_degree);
    //if(spline_degree == 2) {
    //  y[n] = deboor2(grid[n], spline_coefs, ext_knots);
    //}
    //else if(spline_degree == 3) {
    //  y[n] = deboor3(grid[n], spline_coefs, ext_knots);
    //}
  }
}
model {
  to_vector(spline_coefs) ~ normal(0, 1);
}
