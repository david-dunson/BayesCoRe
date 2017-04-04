setwd("~/git/empiricalTensor/hmc_constraint/")
# setwd("C:/Users/leo/git/empiricalTensor/hmc_constraint/")


require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


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
d = 5
p = c(0.1,0.3,0.6)
mu<- c(2,5,7)
sigma=0.5
y<- numeric(N)
alpha =0.5

for(i in 1:N){
  y[i] = rnorm(1, sum(rmultinom(1,1,p)*mu),sigma)
}


lambda=1E3

hist(y,breaks=50)


toy_dat <- list(N=N, d=d,y=y, alpha=alpha, lambda=lambda)
init<- list(list(p= rep(1/d,d), mu=rnorm(d),sigma=1))

ss_fit <- sampling(ss_model, data = toy_dat,init=init, iter = 20000, chains = 1, algorithm = "NUTS")

sampling_idx<- c(10001:20000)


post_p<- extractPosterior("p", d,"ss_fit")
ts.plot(post_p[sampling_idx,])

colMeans(post_p)

save(ss_fit,file="hmc_dp.RDa")

rowSums(post_p)  

post_mu<- extractPosterior("mu", d,"ss_fit")
ts.plot(post_mu[sampling_idx,c(1:3)])

post_sigma<- extractPosterior("sigma", 1,"ss_fit")
ts.plot(post_sigma)



post_p = post_p[sampling_idx,]


for(i in 1:(d-1)){
	print(sum(post_p[,i] > post_p[,i+1]))
}

rowSums(post_p)


acf(post_p[9000:10000,1:3],lag.max=100)