data {
  int<lower=1> N;
  int<lower=1> d;
  int<lower=1> p;
  vector[N] x;
  matrix[N,p] y;
  real<lower=0> jitter;
  real lambda;
}
transformed data {
  vector[N] mu;
  for (i in 1:N) 
    mu[i] = 0;
}
parameters {
  matrix[N,d] g;
  matrix[d,p] eta;


  real<lower=0> eta_sq[d];
  real<lower=0> rho_sq[d];
  real<lower=0> sigma_sq;
}
model {


  for(m in 1:d){
  
    matrix[N,N] Sigma;

    //extrinsic prior for modeling constraint
    matrix[d,d]  G2 = g' * g - diag_matrix(rep_vector(1.0, d));
    target += - lambda* trace( G2 *G2);

    // off-diagonal elements
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        Sigma[i,j] = eta_sq[m] * exp(-rho_sq[m] * pow(x[i] - x[j],2));
        Sigma[j,i] = Sigma[i,j];
      }
    }

    // diagonal elements
    for (k in 1:N)
      Sigma[k,k] = eta_sq[m] + jitter;

    g[:,m] ~ multi_normal(mu,Sigma);

    eta_sq[m] ~ cauchy(0,5);
    rho_sq[m] ~ cauchy(0,5);

  }
  
  target +=  normal_lpdf( to_vector(y) | to_vector(g*eta), sqrt(sigma_sq));

  sigma_sq ~ cauchy(0,5);


}