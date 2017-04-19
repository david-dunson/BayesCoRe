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
d=3;
alpha= rep(0.01,d)

alpha = c(10,1,1)

#flexible lambda

lambda=1E3

toy_dat <- list(d=d, alpha=alpha, lambda=lambda)
init<- list(list(p= rep(1/d,d)))

iter =20000

ss_fit <- sampling(ss_model, data = toy_dat, iter = iter, chains = 1, algorithm="HMC")
sample_idx = c((iter/2+1):iter)


# tau = extractPosterior("tau", 1,"ss_fit")
# hist(tau[sample_idx])
post_p<- extractPosterior("p", d,"ss_fit")


require("rgl")

plot3d(post_p)

acf(post_p[sample_idx,1:3])
# acf(post_p[sample_idx,2])
# acf(post_p[sample_idx,3])

