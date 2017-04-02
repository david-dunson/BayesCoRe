# setwd("~/work/git/empiricalTensor/test/")
setwd("~/git/empiricalTensor/hmc_constraint")

options(max.print=100)

require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#####################################################
ss_model = stan_model(file= "dirichlet.stan")




extractPosterior<-function(varname, dimen, stan_fit){
	if(dimen >1){
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
	}else{
		eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$",varname,sep="")))
	}
}

##### Data ###
d=10;
alpha= 0.5
C= c(1000,150,300,rep(0,7))


#flexible lambda

lambda=1E8
beta_a1 = 1/lambda
beta_a2 = lambda


toy_dat <- list(d=d, alpha=alpha, C =C,beta_a1=beta_a1, beta_a2= beta_a2)  #, lambda=lambda)
init<- list(list(p= rep(1/d,d)))

iter =2000

ss_fit <- sampling(ss_model, data = toy_dat, iter = iter, chains = 1, algorithm="HMC")
sample_idx = c((iter/2+1):iter)


tau = extractPosterior("tau", 1,"ss_fit")
hist(tau[sample_idx])
post_p<- extractPosterior("p", d,"ss_fit")

ts.plot(post_p[sample_idx,1:10],ylim=c(0,1))

hist(rowSums(post_p[sample_idx,]))

acf(post_p[sample_idx,1])
acf(post_p[sample_idx,2])
acf(post_p[sample_idx,3])


hist(rowSums(post_p[sample_idx,]))

hist(tau[sample_idx])

sum(tau[sample_idx]<1*1E-8)

#fixed lambda

ss_model_fixed_lambda = stan_model(file= "dirichlet_fixed_lambda.stan")

beta_a1 = 1E-8
# beta_a2 = 1E4
toy_dat <- list(d=d, alpha=alpha, C =C,beta_a1=beta_a1, beta_a2= beta_a2, dist=1)  #, lambda=lambda)

init<- list(list(p= rep(1/d,d)))

ss_fit_fixed_lambda <- sampling(ss_model_fixed_lambda, data = toy_dat, iter = iter, chains = 1, algorithm="HMC")

post_p2<- extractPosterior("p", d,"ss_fit_fixed_lambda")

ts.plot(post_p2[sample_idx,1:10],ylim=c(0,1))
acf(post_p2[sample_idx,1])

hist(rowSums(post_p2[sample_idx,]))

acf(post_p2[sample_idx,1])
acf(post_p2[sample_idx,2])
acf(post_p2[sample_idx,3])

sample_p2 = post_p2[sample_idx,]

err = abs(rowSums(sample_p2)-1)
sum(err<1E-5)