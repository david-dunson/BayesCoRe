data {
  # C: n by r
  # theta: r by p 
  int<lower=0> n;
  int<lower=0> r; 
  int<lower=0> p;
  vector[r] w;
  vector[r] g;

  matrix[n,p] Y;
  real lambda1;
  real lambda2;

}

parameters {
  matrix[r,p] alpha;
  matrix[r,p] beta;
  real<lower=0, upper=1> x[n];
  real<lower=0> sigma;
  real<lower=0> tau;
}

model {

  matrix[n, r] xw;
  matrix[n, p] deriX;
  matrix[n, p] deri2X;


  xw = to_vector(x) * w' ;

  for(j in 1:p){
    target+= normal_lpdf( alpha[:,j] | 0, sqrt(tau * g));
    target+= normal_lpdf( beta[:,j]| 0, sqrt(tau * g));
  }
  
  target += normal_lpdf(to_vector(Y) | to_vector( cos(xw) * alpha + sin(xw) * beta), sqrt(sigma));
  
  target += gamma_lpdf(1.0 ./sigma | 2, 1);
  target += gamma_lpdf(1.0 ./tau | 2, 1);

  deriX = -sin(xw) * diag_matrix(w) * alpha + cos(xw) * diag_matrix(w) * beta  ;
  deri2X = -cos(xw) * diag_matrix(w .* w) * alpha - sin(xw) * diag_matrix(w .* w) * beta;

  for(i in 1:n){
    for(j in 1:p){
      if(deriX[i,j]>0){
        target += - lambda1 * deriX[i,j];
      }
      if(deri2X[i,j]>0){
        target += - lambda2 * deri2X[i,j];
      }
    }
  }
  
  
  
}
