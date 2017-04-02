data {
  int<lower=0> d; //dirichlet dimension
  real<lower=0> C[d];
  real<lower=0> alpha;
  real beta_a1;
  real beta_a2;
  int dist;
}

parameters {
  real<lower=0,upper=1> p[d];
  # real<lower=0,upper=1> tau; //inverse of lambda
}

model {

  real kappa;


  for (i in 1:d) {
      target += log(p[i]) * (C[i] + alpha/d) ;
  }
  # target += - 1.0./beta_a1 * (sum(p) -1)* (sum(p)-1);
 # target += - 1.0./beta_a1 * fabs(sum(p) -1);

 if(dist ==1) // normal
  target += normal_lpdf( (sum(p) -1) | 0 , beta_a1);
  if(dist==2) //laplace
  target += double_exponential_lpdf( (sum(p) -1) | 0 , beta_a1);
  if(dist==3){// horseshoe
    target += normal_lpdf( (sum(p) -1) | 0 , fabs(kappa));
    target += cauchy_lpdf(kappa | 0, beta_a1);
  }
  if(dist==4){// beta
    target += (beta_a1-1) * log(fabs(sum(p) -1)) + (beta_a2-1) * log(1- fabs(sum(p) -1));
  }  

}
