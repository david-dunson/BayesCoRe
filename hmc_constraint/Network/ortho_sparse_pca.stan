data {
  int<lower=1> V;
  int<lower=1> n;
  int<lower=1> d;
  int y[V,V,n];
  real lambda1;

}
parameters {
  matrix[V,d] U0;
  vector[d] D;
  vector[n] z;

  matrix[V,d] etaU;// C(0,1)
  real phiU; //C(0,1)
  real<lower=0> phiD;
  real<lower=0> phiZ;
  real<lower=0> phi3;
}
transformed parameters {
  matrix[V,d] U= qr_Q(U0)[, 1:d]* sqrt(V - 1);
  // matrix[d,d] W= qr_R(U0)[1:d, 1:d] / sqrt(V - 1);
}


model {

    matrix[V,V] UDU;
    matrix[V,V] U2;
    
  
    UDU = U * diag_matrix(to_vector(D)) * U'; 

    // logit bernoulli with random effect
    for(l in 1:n){
      for(i in 1:V){
        for(j in 1:(i-1)){
          target += bernoulli_logit_lpmf(y[i,j,l] | UDU[i,j]+z[l]);
        }
      }
    }

    // horseshoe prior
    target +=  normal_lpdf( to_vector(U) | 0, fabs(to_vector(etaU) * phiU));
    target +=  cauchy_lpdf( to_vector(etaU) | 0, 1);
    target +=  cauchy_lpdf( phiU | 0, 1);

    //variance for the loading
    target +=  normal_lpdf( D | 0, sqrt(phiD));
    target +=  gamma_lpdf( 1.0/phiD | 2,1);

    //variance for random effect
    target +=  normal_lpdf( z | 0, sqrt(phiZ));
    target +=  gamma_lpdf( 1.0/phiZ | 2,1);

    //constrain orthonormality in U
    U2 = U0' * U0 - diag_matrix(rep_vector(1.0, d));
    target += - lambda1* trace( U2 *U2);

    
  
   

}
