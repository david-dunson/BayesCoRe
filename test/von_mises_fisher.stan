data {
  int<lower=0> N;
  real lambda;
  vector[N] F; 
}

parameters {
  vector[N] X; 
}

model {
  real x2 = X' * X;

target += X'*F;
target += - lambda* (x2 -1)* (x2-1);
}
