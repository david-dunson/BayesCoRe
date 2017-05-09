setwd("~/git/empiricalTensor/hmc_constraint/")

options(max.print=100)

require("rstan")
rstan_options(auto_write = TRUE)


langevin_model = stan_model(file= "langevin.stan")

N = 2
d = 1

F = matrix(1,N,d)
Xini = matrix(0,N,d)
for(i in 1:d){
  Xini[i,i]<-1
}

input_dat <- list(N=N,d=d, F=F, lambda= 1E3)

iter = 20000
langevin_fit <- sampling(langevin_model , data = input_dat, init = list(list(X= Xini)),iter = iter, chains = 1, algorithm= "HMC")

posterior= numeric()
for(j in 1:d){
  each_col = sapply(c(1:N), function(i)  eval(parse(text=paste("langevin_fit@sim$samples[[1]]$`X[",i,",",j, "]`",sep=""))))
  posterior = cbind(posterior,each_col)
}

X0 =matrix(posterior[iter,],N,d)

t(X0)%*%X0

sample_idx = c((iter/2+1):iter)


acf(posterior[sample_idx = c((iter/2+1):iter),1], lag.max=40)

plot(posterior[sample_idx,c(1,2)][1:1000,],type="l")

plot(posterior[sample_idx,c(1,2)][1:1000,],type="p")



