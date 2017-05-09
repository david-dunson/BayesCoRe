data {
  int<lower=0> N;
  int<lower=0> d;
  real lambda;
  matrix[N,d] F; 
}

parameters {
  matrix[N,d] X; 
}

model {
  matrix[d,d]  x2= X' * X - diag_matrix(rep_vector(1.0, d));
  target += trace(F' * X);
  target += - lambda* trace( x2 *x2);
}
