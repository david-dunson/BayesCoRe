setwd("~/git/empiricalTensor/hmc_constraint/")
# setwd("C:/Users/leo/git/empiricalTensor/hmc_constraint/")


require("rstan")
require("splines")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

##################################

N = 100
p = 20
d = 4

x= rnorm(N)

F = cbind(sin(x), cos(x), sin(2*x), cos(2*x))


#sort the columns so that they appear decreasing
sigma = sort( rgamma(4,1,1), dec= T)

L = matrix( rnorm(d*p),d)* sqrt(sigma)

Y = F%*% L + rnorm(N*p) * 0.1

plot(x, Y[,1])
plot(x, Y[,2])
plot(x, Y[,4])




Bd = 10
B = as.matrix( bs(x, df = Bd))

#######################################################
ss_model = stan_model(file= "fpca.stan")

extractPosterior<-function(varname, dimen, stan_fit){
  if(dimen==1){
    eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"`",sep="")))
  }else{
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
  }
}

##### Data ###
lambda=1E3
input_dat <- list(N=N, p=p,d=d, Bd= Bd, x=x, Y=Y, B=B ,lambda=lambda)
# init<- list(list(p= rep(1/d,d), mu=rnorm(d),sigma=1))

ss_fit <- sampling(ss_model, data = input_dat, iter = 20000, chains = 1, algorithm = "NUTS")

save(ss_fit,file="fpca.RDa")

load(file="fpca.RDa")

sampling_idx<- c(10001:20000)

post_sample = do.call("cbind",ss_fit@sim$samples[[1]])

post_sample = post_sample[sampling_idx,] 

post_F = post_sample[, substr( colnames(post_sample), 1, 1) =="F"]
post_F_tensor = array(post_F, dim=c(10000,N,d))

F_sample = post_F_tensor[2,,]

t(F_sample) %*% F_sample

plot(x,F_sample[,1])
plot(x,F_sample[,2])
plot(x,F_sample[,3])
plot(x,F_sample[,4])

acf(post_F_tensor[,1,3],lag.max=40)


post_L = post_sample[, substr( colnames(post_sample), 1, 1) =="L"]
post_L_tensor = array(post_L, dim=c(10000,d,p))

image(post_L_tensor[4,,])




post_p = post_p[sampling_idx,]


for(i in 1:(d-1)){
	print(sum(post_p[,i] > post_p[,i+1]))
}

rowSums(post_p)


acf(post_p[9000:10000,1:3],lag.max=100)