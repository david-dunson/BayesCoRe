data {
  int<lower=1> N;
  vector[3] F;
  real lambda1;
  real sigma2;
  matrix[3,3] Sigma;
  real m; //degrees of freedom in t
}
parameters {
  vector[3] theta_vmf;
  vector[3] theta_fb;
  vector[3] theta_t;
  // real<lower=0> tau;
}
model {

    target +=  theta_vmf' * F /sigma2;
    target +=  multi_normal_lpdf(theta_fb | F, Sigma*sigma2);
    
    target +=  -(m+3)/2* log(
      1+ (1+F'*F)/m/sigma2
      - 2* theta_t' * F/m/sigma2);

    //constrain orthonormality in U
    target += - fabs(theta_vmf' * theta_vmf -1) / lambda1;
    target += - fabs(theta_fb' * theta_fb -1) / lambda1;
    target += - fabs(theta_t' * theta_t -1) / lambda1;
}
