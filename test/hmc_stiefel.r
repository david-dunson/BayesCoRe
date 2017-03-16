setwd("~/work/git/empiricalTensor/test/")

require("rstan")
rstan_options(auto_write = TRUE)


#######################################################
# one dimensional Langevin
toy_model = stan_model(file= "toy_code.stan")
##### Data ###
N = 2
F = rep(0,N)
toy_dat <- list(N=N, F=F, lambda= 1E3)
toy_fit3 <- sampling(toy_model, data = toy_dat, init = list(list(X=c(1,rep(0,N-1)))),
         iter = 1000, chains = 1)

posterior = sapply(c(1:N), function(i)  eval(parse(text=paste("toy_fit3@sim$samples[[1]]$`X[",i, "]`",sep=""))))

plot(posterior)
###
#
#
langevin_model = stan_model(file= "langevin.stan")

N = 100
d = 5

F = matrix(1,N,d)
Xini = matrix(0,N,d)
for(i in 1:d){
  Xini[i,i]<-1
}
input_dat <- list(N=N,d=d, F=F, lambda= 1E4)

langevin_fit <- sampling(langevin_model , data = input_dat, init = list(list(X= Xini)),
                     iter = 1000, chains = 1)


posterior= numeric()
for(j in 1:d){
  each_col = sapply(c(1:N), function(i)  eval(parse(text=paste("langevin_fit@sim$samples[[1]]$`X[",i,",",j, "]`",sep=""))))
  posterior = cbind(posterior,each_col)
}

X0 =matrix(posterior[100,],N,d)

t(X0)%*%X0

acf(posterior[,4])
