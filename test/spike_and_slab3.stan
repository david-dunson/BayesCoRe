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

  vector[p] DWy = (D-W)*y;

  vector[p] ZXbeta= (Z' * X* diag_matrix(beta))';
  vector[p] dbeta2 = D * beta .* beta;
  vector[p] yDWBeta = (DWy' * diag_matrix(beta))';
 
 
   vector[p] includ_prob;
  includ_prob =     inv_logit( ZXbeta + yDWBeta - dbeta2/2 + log( prob/(1-prob)) )  ;
}


model {

  // target += 0.5 * log_determinant(D-W)
  target += -0.5 * y'* DWy;
  
  target += log_sum_exp( p*log(1-prob), sum( log(prob) +  ZXbeta + yDWBeta - dbeta2/2));
  // target+= sum(log(  1 - prob + exp(ZXbeta + yDWBeta - dbeta2/2 + log( prob))  ));

  target +=  normal_lpdf( beta | 0,1000);
  target += beta_lpdf(prob| 0.5,0.5);

}


