data {
  int<lower=0> N;
  int y[N];
}

parameters {
  real<lower=0> lambda;
}

model {
  lambda ~ normal(0, 5.82337);
  y ~ poisson(lambda);
}

generated quantities {
  int y_ppc[N];
  for (n in 1:N)
    y_ppc[n] = poisson_rng(lambda);
}
