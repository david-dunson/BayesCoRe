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
ss_model = stan_model(file= "ortho_tensor.stan")

##### Data ###
N = 20
d1 = 2
d2 = 3
p = 5


U = matrix(rnorm(N*d1),N,d1)
U = qr.Q(qr(U))
y = array(0, c(N,N,p))


for(l in 1:p){
  eta = rnorm(d1)
  UDU = U%*%diag(eta)%*%t(U)
  p_mat = 1/(1+exp(-UDU))
  y_l = (matrix(runif(N*N),N)<p_mat)*1
  Lower=lower.tri(y_l)
  y_l[Lower] = t(y_l)[Lower]
  y[,,l] = y_l
}

lambda1=0   # ordering in tau
lambda2=1E3  # orthonormality
lambda3=1E3  # positive

input_dat <- list(N=N, p=p,d1=d1,d2=d2 , y=y,lambda1 = lambda1, lambda2 = lambda2,lambda3= lambda3)

n_steps = 1E4
ss_fit <- sampling(ss_model, data = input_dat, iter = n_steps, chains = 1, algorithm = "NUTS")

# data = list("y"=y,"g"=g,"x"=x,"ss_fit"=ss_fit)
# save(ss_fit,file="network1.RDa")

load(file="network1.RDa")

sampling_idx<- c((n_steps/2+1):n_steps)

post_U = extractPosteriorMat("U",N,d1,"ss_fit")

acf(post_U[sampling_idx,1])
acf(post_U[sampling_idx,N+1])

ts.plot(post_U[sampling_idx,1])
ts.plot(post_U[sampling_idx,N+2])


post_V = extractPosteriorMat("V",p,d2,"ss_fit")

acf(post_V[sampling_idx,1])
acf(post_V[sampling_idx,N+1])

ts.plot(post_V[sampling_idx,1])
ts.plot(post_V[sampling_idx,N+2])

post_eta = extractPosterior3D("eta",d1,d1,d2,"ss_fit")

ts.plot(post_eta[sampling_idx,1])
acf(post_eta[sampling_idx,2])


d=d1

U1 = matrix(0,N,d1)
for(j in 1:d1){
  U1[,j]=post_U[n_steps/2,((j-1)*N+1):(j*N)]
}
t(U1)%*%U1



