data {
  int<lower=0> N;
  int<lower=0> p;
  int<lower=0> d; // latent dimension
  int<lower=0> Bd; // dimension of B spline basis
  vector[N] x;
  matrix[N,p] Y;
  matrix[N,Bd] B; //B spline basis
  real lambda;
}

parameters {
  real<lower=0, upper = 10> tau;
  matrix[Bd,d] Bcoef;
  matrix[d,p] L;
}

transformed parameters {
  matrix[N,d] F; 
  F = B * Bcoef;
}

model {

  //extrinsic prior for modeling constraint
  matrix[d,d]  F2 = F' * F - diag_matrix(rep_vector(1.0, d));
  target += - lambda* trace( F2 *F2);


  //likelihood
  target += normal_lpdf( to_vector(Y)| to_vector(F*L), sqrt(tau));

  //prior
  target += normal_lpdf( to_vector(Bcoef)| 0, 1000);
  target += normal_lpdf( to_vector(L)| 0, 1000);


}
