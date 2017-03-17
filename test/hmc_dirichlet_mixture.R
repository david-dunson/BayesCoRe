# setwd("~/work/git/empiricalTensor/test/")
setwd("D:/work/git/empiricalTensor/test/")



require("rstan")
rstan_options(auto_write = TRUE)


#######################################################
ss_model = stan_model(file= "dp_mixture.stan")



extractPosterior<-function(varname, dimen, stan_fit){
  if(dimen==1){
    eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"`",sep="")))
  }else{
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
  }
}

##### Data ###
N = 100
d=4;
p = c(0.1,0.2,0.7)
mu<- c(1,2,5)
sigma=0.5
y<- numeric(N)
alpha =0.5

for(i in 1:N){
  y[i] = rnorm(1, sum(rmultinom(1,1,p)*mu),sigma)
}


lambda=1E5

toy_dat <- list(N=N, d=d,y=y, alpha=alpha, lambda=lambda)
init<- list(list(p= rep(1/d,d), mu=rnorm(d),sigma=1))

ss_fit <- sampling(ss_model, data = toy_dat,init=init, iter = 20000, chains = 1)


post_p<- extractPosterior("p", d,"ss_fit")
ts.plot(post_p[,c(3,4)])


save(ss_fit,file="hmc_dp.RDa")

rowSums(post_p)  


post_mu<- extractPosterior("mu", d,"ss_fit")
ts.plot(post_mu[10000:20000,c(3,4)])

post_sigma<- extractPosterior("sigma", 1,"ss_fit")
ts.plot(post_sigma)
