data {
  int N;
}

generated quantities {
  //fabs(x) return the absolute value of a float point number
  real<lower=0> lambda = fabs(normal_rng(0, 5.82337));
  int y[N];
  for (n in 1:N) y[n] = poisson_rng(lambda);
}
