data {
  int N;
}

generated quantities {
  // Simulate the mixture weight from a Beta(1, 1)
  real<lower=0, upper=1> theta = beta_rng(1, 1);
  real<lower=0> lambda = fabs(normal_rng(0, 5.82337));
  int y[N] = rep_array(0, N);
  for (n in 1:N)
    if (!bernoulli_rng(theta))
      y[n] = poisson_rng(lambda);
}
