data {
  int<lower=0> d; //dirichlet dimension
  vector[d] alpha;
  # real beta_a1;
  # real beta_a2;
  real lambda;
}

parameters {
  real<lower=0,upper=1> p[d];
  # real<lower=0,upper=1> tau; //inverse of lambda
}

model {

  target += sum( to_vector(log(p)) .* (alpha -1.0));
  target += - lambda * (sum(p) -1)* (sum(p)-1);

  # target += - 1.0./tau * (sum(p) -1)* (sum(p)-1);
  # target +=  beta_lpdf( tau | beta_a1, beta_a2);

  # target +=  -log(tau) + beta_a2 * log(1-tau);
}
