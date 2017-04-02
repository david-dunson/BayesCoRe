data {
  int<lower=0> d; //dirichlet dimension
  real<lower=0> C[d];
  real<lower=0> alpha;
  real beta_a1;
  real beta_a2;
}

parameters {
  real<lower=0,upper=1> p[d];
  real<lower=0,upper=1> tau; //inverse of lambda
}

model {


  for (i in 1:d) {
      target += log(p[i]) * (C[i] + alpha/d) ;
  }
  target += - 1.0./tau * (sum(p) -1)* (sum(p)-1);
  # target +=  beta_lpdf( tau | beta_a1, beta_a2);

  target +=  -log(tau) + beta_a2 * log(1-tau);
}
