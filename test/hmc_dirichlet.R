# setwd("~/work/git/empiricalTensor/test/")
setwd("D:/work/git/empiricalTensor/test/")



require("rstan")
rstan_options(auto_write = TRUE)


#######################################################
ss_model = stan_model(file= "dirichlet.stan")



extractPosterior<-function(varname, dimen, stan_fit){
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
}

##### Data ###
d=10;
alpha= 0.5
C= c(1000,150,300,rep(0,7))
lambda=1E5

toy_dat <- list(d=d, alpha=alpha, C =C, lambda=lambda)
init<- list(list(p= rep(1/d,d)))

ss_fit <- sampling(ss_model, data = toy_dat, iter = 20000, chains = 1)


post_p<- extractPosterior("p", d,"ss_fit")
ts.plot(post_p[10000:20000,1:10],ylim=c(0,1))

colMeans(post_p[10000:20000,])  


ts.plot(post_p[,1:3])

acf(post_p[,1])
acf(post_p[,2])
acf(post_p[,3],lag.max = 100)
