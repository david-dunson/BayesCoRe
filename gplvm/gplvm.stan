data {
  # C: n by r
  # theta: r by p 
  int<lower=0> n;
  int<lower=0> r; 
  int<lower=0> p;
  vector[r] w;
  vector[r] g;

  matrix[n,p] Y;
}

parameters {
  matrix[r,p] alpha;
  matrix[r,p] beta;
  real<lower=0, upper=1> x[n-1];
  real<lower=0> sigma[p];
}

model {

  matrix[n, r] xw;
  matrix[n, 2*r] C;
  matrix[2*r, p] theta;
  matrix[n, p] A;
  vector[n] x1;

  x1[1]=0.1;
  for(i in 1:(n-1)){
    x1[i+1] = x[i];
  }

  target+= uniform_lpdf(x| 0,1);

  for(i in 1:n){
    xw[i,] = x1[i]* w';
  }

  C = append_col(cos(xw), sin(xw));

  for(j in 1:p){
    target+= normal_lpdf( alpha[,j] | 0, g);
    target+= normal_lpdf( beta[,j]| 0, g);
  }

  theta = append_row(alpha,beta);

  A = Y - C * theta;

  target += - 0.5* trace_quad_form( diag_matrix(1.0 ./ to_vector(sigma)), A');
  target += - 0.5 * sum(log(sigma)) * n;
  
  target += gamma_lpdf(1.0 ./to_vector(sigma) | 0.1, 0.1);

  # sum(diag(A * diag_matrix(1/sigma) * A')

}
