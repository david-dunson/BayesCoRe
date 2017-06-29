data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=1> d1;
  int<lower=1> d3;

  int y[N,N,p];
  real lambda1;
  real lambda2;
  real lambda3;
  real a1;
  real a2;
}
parameters {
  matrix[N,d1] U;
  matrix[p,d3] W;


  matrix[N,N] Z;

  matrix[d1,d3] D;
  
  real<lower=0> nu1[d1];
  real<lower=0> nu2[d3];

  real<lower=0> phi1;
  real<lower=0> phi2;
  real<lower=0> phi3;
}
model {

    matrix[N,N] UDU;
    matrix[d1,d1]  U2;
    matrix[N,N] Z1;
    matrix[d1,p] eta;
    
    real tau[d1,d3];
    real nu1_cumprod[d1];
    real nu2_cumprod[d3];


    eta = D * W';

    for(i in 1:N){
        for(j in 1:i){
          target +=  normal_lpdf(Z[i,j]| 0, sqrt(phi3));
          Z1[i,j] = Z[i,j];
          Z1[j,i] = Z[i,j];
        }
    }
    

    //tensor product
    for(l in 1:p){
      UDU = U * diag_matrix(to_vector(eta[:,l])) * U' + Z1; 
      for(i in 1:N){
        for(j in 1:(i-1)){
          target += bernoulli_logit_lpmf(y[i,j,l] | UDU[i,j]);
        }
      }
    }
    
    target +=  gamma_lpdf( 1.0/nu1[1] | a1, 1);
    nu1_cumprod[1] = nu1[1];
    target +=  gamma_lpdf( 1.0/nu2[1] | a1, 1);
    nu2_cumprod[1] = nu2[1];
    
    for(i in 2:d1){
        target +=  gamma_lpdf( 1.0/nu1[i] | a2, 1);
        nu1_cumprod[i] = nu1_cumprod[i-1] * nu1[i];
    }
    for(i in 2:d3){
        target +=  gamma_lpdf( 1.0/nu2[i] | a2, 1);
        nu2_cumprod[i] = nu2_cumprod[i-1] * nu2[i];
    }

    //variance for the core
    for(i in 1:d1){
        for(j in 1:d3){
        tau[i,j] = nu1_cumprod[i] * nu2_cumprod[j];
        target +=  normal_lpdf(D[i,j]| 0, sqrt(tau[i,j]));
        }
    }

    //variance for the factors
    target +=  normal_lpdf( to_vector(U) | 0, sqrt(phi1));
    target +=  normal_lpdf( to_vector(W) | 0, sqrt(phi2));

    target +=  gamma_lpdf( 1.0/phi1 | 2,1);
    target +=  gamma_lpdf( 1.0/phi2 | 2,1);
    target +=  gamma_lpdf( 1.0/phi3 | 2,1);

    //constrain orthonormality in U
    U2 = U' * U - diag_matrix(rep_vector(1.0, d1));
    target += - lambda2* trace( U2 *U2);

  
   

}
