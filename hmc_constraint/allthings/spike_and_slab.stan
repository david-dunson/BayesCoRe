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
   matrix[N,p] Xbeta= X* diag_matrix(beta);
   vector[p] ZXbeta= Xbeta' * Z;
 
   vector[p] includ_prob;
  includ_prob =    (exp(ZXbeta + y - d/2)) * prob ./ ( 1 - prob + exp(ZXbeta + y - d/2) * prob);  
}


model {
  
  // vector[N] XbetaY= Xbeta * y;
  // real part1 =  y'*y / d ;
  // matrix[N,N] temp =  (Xbeta * Xbeta' /d - diag_matrix(rep_vector(1.0, N)));
  // real part2 = - XbetaY' * ( temp \XbetaY) / d/ d;
  // real det_part =   log_determinant(-temp) + p*log(d);
  // target+=  - 0.5* (part1+part2) - 0.5* det_part;
  
  matrix[p,p] W = Xbeta' * Xbeta;
  matrix[p,p] D = diag_matrix(rep_vector(d, p)); 
  
  target += -0.5 * y' * inverse(- W +D ) * y;
  target += -0.5* log_determinant(-W+D);
  
  target += log_sum_exp(p * log(1-prob), sum( log(prob) + (ZXbeta + y - d/2)) );
  // target+= sum(log(  1 - prob + exp(ZXbeta + y - d/2) * prob));

  
  target +=  normal_lpdf( beta | 0,1000);
  target += beta_lpdf(prob| 1,1);

}

// 
// generated quantities{
// 
// }
