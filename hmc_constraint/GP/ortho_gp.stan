data {
  int<lower=1> N;
  int<lower=1> d;
  int<lower=1> p;
  vector[N] x;
  matrix[N,p] y;
  real<lower=0> jitter;
  real lambda1;
  real lambda2;
  real lambda3;
  # vector[d] rho;
}
transformed data {
  vector[N] mu;
  for (i in 1:N) 
    mu[i] = 0;
}
parameters {
  matrix[N,d] g;
  matrix[d,p] eta;

  real<lower=0> tau[d];

  real<lower=0, upper = 10> phi[d];
  real<lower=0, upper = 5> rho[d];
  real<lower=0> sigma_sq;
}
model {


  for(m in 1:d){
  
    matrix[N,N] Sigma;

    //extrinsic prior for modeling constraint
    matrix[d,d]  G2 = g' * g - diag_matrix(rep_vector(1.0, d));
    target += - lambda2* trace( G2 *G2);

    //for(i in 1:d){
    //  target += -lambda2 * pow((g[:,i])' * g[:,i] - 1, 2);
    //}


    // off-diagonal elements
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Sigma[i,j] = phi[m] * exp(- pow(x[i] - x[j],2) / 2/rho[m]/rho[m]);
        Sigma[j,i] = Sigma[i,j];
      }
    }

    // diagonal elements
    for (k in 1:N)
      Sigma[k,k] = phi[m] * (1.0 + jitter);

    g[:,m] ~ multi_normal(mu,Sigma);

    target +=  normal_lpdf( to_vector(eta[m,:]) | 0, sqrt(tau[m]));
    target +=  gamma_lpdf( 1.0/tau[m] | 2, 1);
    target +=  normal_lpdf( rho[m] | 0,1);
    target +=  gamma_lpdf( 1.0/phi[m] | 2,1);

    //constrain first element in each factor to be positive
    if( g[1,m] < 0){
      target +=   -lambda3 * (-g[1,m]);
    }
  }

  for(m in 1:(d-1)){
    if( tau[m] < tau[m+1]){
      target +=   -lambda1 *( tau[m+1] - tau[m]);
    }
  }
  
  target +=  normal_lpdf( to_vector(y) | to_vector(g*eta), sqrt(sigma_sq));

  sigma_sq ~ cauchy(0,5);


}
