data {
  int<lower=0> N;
  int y[N];
}

parameters {
  real<lower=0, upper=1> theta;
  real<lower=0> lambda;
}

model {
  theta ~ beta(1, 1);
  lambda ~ normal(0, 5.82337);
  for (n in 1:N) {
    real lpdf = poisson_lpmf(y[n] | lambda);
    if (y[n] == 0)
      target += log_mix(theta, 0, lpdf);
    else
      target += log(1 - theta) + lpdf;
  }
}
