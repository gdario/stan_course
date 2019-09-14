// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int y[N];
}

parameters {
  real<lower=0> lambda;
}

model {
  // This is actually a half normal distribution due to the <lower=0> constraint
  lambda ~ normal(0, 5.82337);
  // y is a vector of N components, one per detector.
  y ~ poisson(lambda);
}
