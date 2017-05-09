data {
  int<lower=0> N;
  int<lower=0> d;
  int<lower=0> m; //number of subjects
  int<lower=0> MC_n;

  matrix[N,d] Xsum;
  real lambda;
}

parameters {
  matrix[N,d] F;
  matrix[N*MC_n, d] Y;
}

model {

  real logprob[MC_n];

  for(i in 1:MC_n){
      int startIdx = (i-1)*N + 1;
      int endIdx = i*N;
      matrix[N, d] Y_i = Y[startIdx:endIdx,];
      matrix[d,d]  y2= Y_i' * Y_i - diag_matrix(rep_vector(1.0, d));
      // Monte Carlo integrand & extrinsic prior
      logprob[i] = trace(F' * Y_i) - lambda* trace( y2 *y2); 
  }

  //Used in the denominator
  target += - m * log_sum_exp( logprob);

  //likelihood
  # for(i in 1:m){
  #     int startIdx = (i-1)*N+ 1;
  #     int endIdx = i*N;
  #     matrix[N, d] X_i = X[startIdx:endIdx,];
  #     target += trace(F' * X_i);
  # }
  target += trace(F' * Xsum);
}
