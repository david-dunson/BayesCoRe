data {
  int<lower=0> N;
  int<lower=0> d; //dirichlet dimension
  real y1[N];
  real y2[N];
  real<lower=0> alpha;
  real lambda1;
  real lambda2;
}

parameters {
  # simplex[d] p;
  real<lower=0, upper=1> p[d];
  real<lower=0, upper=10> mu1[d];
  real<lower=0, upper=10> mu2[d];
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

model {
  
  for(i in 1:N){
    real logprob[d];
    for (j in 1:d) {
      logprob[j] = log(p[j]) + normal_lpdf( y1[i] |  mu1[j], sigma1 ) + normal_lpdf( y2[i] |  mu2[j], sigma2 ); 
    }
    target += log_sum_exp( logprob);
  }

  for(j in 1:d){
     target += log(p[j]) * (alpha -1 );
  }
  //prior for sigma and mu
  target += normal_lpdf(mu1 | 0, 3.16);
  target += normal_lpdf(mu2 | 0, 3.16);

  target += gamma_lpdf( 1 ./sigma1 ./sigma1 | 2, 1);
  target += gamma_lpdf( 1 ./sigma2 ./sigma2 | 2, 1);

  // for(j in 1:d){
    // target += gamma_lpdf( 1 ./sigma[j] ./sigma[j] | 2, 1);
  // }
  
  // order constraint in p
  for(i in 1:(d-1)){
    if(p[i]<p[i+1]){
      target += - lambda1 *  (p[i+1] - p[i]) ;
    }
  }
  
    // order constraint in mu
 //   for(i in 1:(d-1)){
  //   if(mu[i]<mu[i+1]){
  //     target += - lambda1 *  (mu[i+1] - mu[i]) ;
  //   }
 //  }

  //simplex constraint
  target += - lambda2* (sum(p) -1)* (sum(p)-1);


}
