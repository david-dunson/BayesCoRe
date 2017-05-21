setwd("~/git/empiricalTensor/hmc_constraint/Network")
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

#model
ss_model = stan_model(file= "ortho_factor.stan")

##### Data ###
N = 10
d = 2
p = 1

U = matrix(rnorm(N*d),N,d)
U = qr.Q(qr(U))

eta = rnorm(d)
UDU = U%*%diag(eta)%*%t(U)
p_mat = 1/(1+exp(-UDU))
y = (matrix(runif(N*N),N)<p_mat)*1
Lower=lower.tri(y)
y[Lower] = t(y)[Lower]

lambda1=0   # ordering in tau
lambda2=1E3  # orthonormality
lambda3=1E3  # positive

input_dat <- list(N=N, d=d,p=p, y=y,lambda1 = lambda1, lambda2 = lambda2,lambda3= lambda3)

init<- list(list(U= U, tau=rep(1,d), phi =1, eta=diag(eta)))


n_steps = 1E4
ss_fit <- sampling(ss_model, data = input_dat, iter = n_steps, chains = 1, algorithm = "NUTS")

# data = list("y"=y,"g"=g,"x"=x,"ss_fit"=ss_fit)
# save(data,file="gp1.RDa")

sampling_idx<- c((n_steps/2+1):n_steps)

post_U = extractPosteriorMat("U",N,d,"ss_fit")

acf(post_U[sampling_idx,1])
acf(post_U[sampling_idx,N+1])

ts.plot(post_U[sampling_idx,1])
ts.plot(post_U[sampling_idx,N+1])

post_eta = extractPosteriorMat("eta",d,d,"ss_fit")

ts.plot(post_eta[sampling_idx,1])
acf(post_eta[sampling_idx,2])


U1 = matrix(0,N,d)
for(j in 1:d){
  U1[,j]=post_U[n_steps/2,((j-1)*N+1):(j*N)]
}
t(U1)%*%U1



eta1 = matrix(0,d,d)
for(j in 1:d){
  eta1[,j]=post_eta[n_steps/2,((j-1)*d+1):(j*d)]
}

acf(post_eta[,1])
acf(post_eta[,2])
ts.plot(post_eta[,3])
acf(post_eta[,4])

