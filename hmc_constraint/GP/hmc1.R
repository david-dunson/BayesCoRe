setwd("~/git/empiricalTensor/hmc_constraint/GP")
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
ss_model = stan_model(file= "ortho_gp.stan")

##### Data ###
N = 20
d = 5
p = 10

x = runif(N)


g= matrix(0,N,d)

g[,1] = sin(x*16)/x
g[,2] = 30*cos(x*5*pi)*(1-x)
g[,3] = 15*cos(x*10*pi)

order_x =order(x)

eta = matrix(rnorm(d*p),d,p)
mu = g%*%eta
y = rnorm(N*p) * 2 + mu


plot(rep(x,3),c(g[,1:3] ), type="n")
for(i in 1:3){
  lines(x[order_x], g[order_x,i])
}


plot(rep(x,p),c(y[,1:p] ), type="n")
for(i in 1:p){
  lines(x[order_x], y[order_x,i])
}


rho = rep(1, d);

lambda1=1E6   # ordering in tau
lambda2=1E3  # orthonormality
lambda3=1E3  # positive

input_dat <- list(N=N, d=d,p=p, x=x, y=y, jitter=1E-2,lambda1 = lambda1, lambda2 = lambda2,lambda3= lambda3, rho = rho)

ss_fit <- sampling(ss_model, data = input_dat, iter = 10000, chains = 1, algorithm = "NUTS")

data = list("y"=y,"g"=g,"x"=x,"ss_fit"=ss_fit)
save(data,file="gp1.RDa")

