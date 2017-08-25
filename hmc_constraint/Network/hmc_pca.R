setwd("~/git/constrainedBayes/hmc_constraint/Network")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#######################################################

ss_model = stan_model(file= "ortho_sparse_pca.stan")

load("../../tensorA.RDa")
y=A[,,1:21]
N = dim(A)[1]
p = dim(A)[3]

d1 = 10
d2 = 10
  
  #initialize using HOSVD
  # require("tensr")
  # y0 = y
  # y0[y0 == 1] = 0.999
  # y0[y0 == 0] = 0.001
  # y0=log(y0)/log(1-y0)
  # hosvdy0 = hosvd(y0, r= c(d1,d1,d2))
  # 
  # U0 = hosvdy0$U[[1]]
  # V0 = hosvdy0$U[[3]]
  # eta0 = hosvdy0$S 


  lambda1=0 # ordering in tau
  lambda2=1E3  # orthonormality
  lambda3=0  # positive

  d3=d2
  
  input_dat <- list(N=N, p=p,d1=d1,d3=d3,a1=10,a2=10, y=y,lambda1 = lambda1, lambda2 = lambda2,lambda3= lambda3)
  
  n_steps = 1E4
  
  L = 7
  ss_fit <- sampling(ss_model, data = input_dat, iter = n_steps, chains = 1
                     ,algorithm="NUTS",
                     control = list(adapt_engaged = T,  max_treedepth=L))

  save(ss_fit,file="network_real_data_constrained.RDa")

  
  quit()
  
  
  
  
  
load(file="network_real_data_unconstrainted.RDa")

sampling_idx<- c((n_steps/2+1):n_steps)

# post_tau = extract(ss_fit,"tau")
# post_eta = extract(ss_fit,"eta")

post_U = extractPosteriorMat("U",N,d1,"ss_fit")

acf(post_U[sampling_idx,1],lag.max = 40)
acf(post_U[sampling_idx,N+1],lag.max = 40)

ts.plot(post_U[sampling_idx,1])
ts.plot(post_U[sampling_idx,N+1])
ts.plot(post_U[sampling_idx,2])
ts.plot(post_U[sampling_idx,3])
ts.plot(post_U[sampling_idx,4])

d=d1

U1 = matrix(0,N,d1)
for(j in 1:d1){
  U1[,j]=post_U[n_steps/2,((j-1)*N+1):(j*N)]
}
t(U1)%*%U1


post_D = extract(ss_fit,"D")

post_D[[1]][500,,]


post_nu1 = extract(ss_fit,"nu1")
post_nu1[[1]][500,]

post_nu2 = extract(ss_fit,"nu2")
post_nu2[[1]][500,]




1/rgamma(1000,shape = 10,rate = 1 )
