data {
  int<lower=1> N;
  int<lower=1> d;
  int y[N,N];
  real lambda1;
}
transformed data {
  vector[N] mu;
  for (i in 1:N) 
    mu[i] = 0;
}
parameters {
  matrix[N,d] U;
  real<lower=0> tau[d];
  real<lower=0> phi;
  matrix[d,d] eta;
  # vector[d] eta;
}
model {

    matrix[N,N] UDU;
    matrix[d,d]  U2;
    matrix[d,d] eta1;


    for(i in 1:d){
      target +=  normal_lpdf(to_vector(eta[i,:]) | 0, sqrt(tau[i]));
      eta1[i,i] =eta[i,i];
      for(j in 1:(i-1)){
        eta1[i,j] = eta[i,j];
        eta1[j,i] = eta[i,j];
      }
    }

    UDU = U * eta1 * U';

    for(i in 1:N){
      for(j in 1:(i-1)){
          target += bernoulli_logit_lpmf(y[i,j] | UDU[i,j]);
      }
    }

    for(i in 1:d){
        target +=  gamma_lpdf( 1.0/tau[i] | 2, 1);
    }

    target +=  normal_lpdf( to_vector(U) | 0, sqrt(phi));
    target +=  gamma_lpdf( 1.0/phi | 2,1);

    //constrain orthonormality in U
    U2 = U' * U - diag_matrix(rep_vector(1.0, d));
    target += - lambda1* max(fabs(U2));
  
}
