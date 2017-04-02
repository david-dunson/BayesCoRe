data {
  int<lower=0> N;
  int<lower=0> d; //dirichlet dimension
  real y[N];
  real<lower=0> alpha;
  real lambda;
}

parameters {
  # simplex[d] p;
  real<lower=0, upper=1> p[d];
  real mu[d];
  real<lower=0> sigma;
}

model {
  
  for(i in 1:N){
    real logprob[d];
    for (j in 1:d) {
      logprob[j] = log(p[j]) + normal_lpdf(y[i]|mu[j], sigma ); 
    }
    target += log_sum_exp( logprob);
  }

  for(j in 1:d){
     target += log(p[j]) * alpha / d;
  }
  //prior for sigma and mu
  target += normal_lpdf(mu | 0, 1000);

  //simplex constraint
  target += - lambda* (sum(p) -1)* (sum(p)-1);

  for(i in 1:(d-1)){
    if(p[i]<p[i+1]){
      target += - lambda *  (p[i+1] - p[i]);
    }
  }

}
