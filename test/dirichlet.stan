data {
  int<lower=0> d; //dirichlet dimension
  real<lower=0> C[d];
  real<lower=0> alpha;
  real lambda;
}

parameters {
  real<lower=0,upper=1> p[d];
}

model {
  for (i in 1:d) {
      target += log(p[i]) * (C[i] + alpha/d) ;
  }
  target += - lambda* (sum(p) -1)* (sum(p)-1);
   // target += - lambda* log( (sum(p) -1)* (sum(p)-1));
  
  // target +=  -  log(1+ (sum(p)-1)* (sum(p)-1) * lambda *lambda);
}
