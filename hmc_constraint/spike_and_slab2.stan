data {
  int<lower=0> N;
  int<lower=0> p;
  matrix[N,p] X;
  vector[N] Z;
  real d;
}

parameters {
  vector[p] beta;
  vector[p] y;
  real<lower=0, upper=1> prob;
}

transformed parameters{

  matrix[p,p] W = X'*X;
  matrix[p,p] D = diag_matrix(rep_vector(d, p));
  matrix[p,p] L = cholesky_decompose(-W +D);

  vector[p] ZXbeta= (Z' * X* diag_matrix(beta))';
  vector[p] dbeta2 = d * beta .* beta;
  vector[p] yLBeta = (y' * L* diag_matrix(beta))';
 
 
   vector[p] includ_prob;
  includ_prob =     inv_logit( ZXbeta + yLBeta - dbeta2/2 + log( prob/(1-prob)) )  ;
}


model {

  
  target += -0.5 * y'*y;
  
  target += log_sum_exp( p*log(1-prob), sum( log(prob) +  ZXbeta + yLBeta - dbeta2/2));
  // target+= sum(log(  1 - prob + exp(ZXbeta + yLBeta - dbeta2/2 + log( prob))  ));

  target +=  normal_lpdf( beta | 0,1000);
  target += beta_lpdf(prob| 0.5,0.5);

}


