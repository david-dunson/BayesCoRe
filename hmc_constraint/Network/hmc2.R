setwd("~/git/constrainedBayes/hmc_constraint/Network")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#######################################################

extractPosterior<-function(varname, dimen, stan_fit){
  if(dimen==1){
    eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"`",sep="")))
  }else{
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
  }
}

extractPosteriorMat<-function(varname, d1,d2, stan_fit){
  	L = lapply(c(1:d2), function(j){
  		sapply(c(1:d1), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i,",",j, "]`",sep=""))))
	})
	do.call("cbind",L)
}


extractPosterior3D<-function(varname, d1,d2,d3, stan_fit){
  n = length( eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",1,",",1,",",1, "]`",sep=""))))
  L = matrix(0, n, d1*d2*d3)
  for(i in 1:d1){
    for(j in 1:d2){
      for(l in 1:d3){
          idx = (i-1)*d2*d3 + (j-1)*d3+ l;
          L[,idx] = eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i,",",j,",",l, "]`",sep="")))
      }
    }
  }
  L
}


#model
ss_model = stan_model(file= "ortho_tensor_diag.stan")
# ss_model = stan_model(file= "ortho_tensor.stan")

##### Data ###

if(FALSE){
  N = 20
  d1 = 5
  d2 = 5
  p = 20
  
  
  U = matrix(rnorm(N*d1),N,d1)
  # U = qr.Q(qr(U))
  y = array(0, c(N,N,p))
  
  
  for(l in 1:p){
    eta = rnorm(d1)* c(d1:1)
    UDU = U%*%diag(eta)%*%t(U)
    p_mat = 1/(1+exp(-UDU))
    y_l = (matrix(runif(N*N),N)<p_mat)*1
    Lower=lower.tri(y_l)
    y_l[Lower] = t(y_l)[Lower]
    y[,,l] = y_l
  }

}
####

load("tensorA.RDa")
y=A[,,1:40]
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
  lambda2=0  # orthonormality
  lambda3=0  # positive

  d3=d2
  
  input_dat <- list(N=N, p=p,d1=d1,d3=d3,a1=10,a2=10, y=y,lambda1 = lambda1, lambda2 = lambda2,lambda3= lambda3)
  
  n_steps = 1E4
  
  L = 7
  ss_fit <- sampling(ss_model, data = input_dat, iter = n_steps, chains = 1
                     ,algorithm="NUTS",
                     control = list(adapt_engaged = T,  max_treedepth=L))

  save(ss_fit,file="network_real_data_unconstrained.RDa")

  
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
