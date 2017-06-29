data {
  int<lower=1> N;
  vector[2] F;
  real y[N,2];
  real lambda1;
}
parameters {
  vector[2] theta;
  real<lower=0> tau;
}
model {

    for(i in 1:N){
      target +=  normal_lpdf(to_vector(y[i,:]) | theta, sqrt(tau));
    }

    target +=  theta' * F;
    target +=  gamma_lpdf( 1.0/tau | 2,1);

    //constrain orthonormality in U
    target += - fabs(theta' * theta -1) / lambda1;
  
}
