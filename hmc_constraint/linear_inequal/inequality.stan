data {
  int<lower=1> N;
  int<lower=1> p;
  matrix[N,p] X;
  vector[N] y;
  real lambda1;
}
parameters {
  real<lower=0,upper=1> theta[p];
  real<lower=0> tau;
}
model {

    target +=  normal_lpdf( y | X * to_vector(theta), sqrt(tau));
    target +=  gamma_lpdf( 1.0/tau | 2,1);
  
    target += normal_lpdf(to_vector(theta) |0 , 10);

    //constrain
    if(sum(theta)-1>0){
      target += - (sum(theta)-1)/ lambda1;
    }
  
}
