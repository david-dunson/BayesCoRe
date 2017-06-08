options(max.print=100)

setwd("~/git/constrainedBayes/gplvm/")

require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#################################################

r = 30
w = seq(0,pi,length.out = r)
n = 10
p = 1

x<- runif(n)
# x[1] = 0.001

spec<- function(rho, phi){
  phi*exp(-1/4 * rho^2 * w^2)
}

g = spec(0.05,1)
min(g)

alpha<- matrix(rnorm(r*p),r) * sqrt(g)
beta<- matrix(rnorm(r*p),r) * sqrt(g)
theta<- rbind(alpha, beta)

Sigma = rep(0.1,p)

xw =outer(x,w,"*")
X = cbind(cos(xw),sin(xw))

y = X%*%theta + matrix(rnorm(n*p),n)%*%diag(Sigma^0.5,1)

# y[,1] = 5*exp(-x) +  rnorm(n) * 0.1



plot(x,y[,1])
# plot(x,y[,2])
# plot(x,y[,3])

#######################################################
ss_model = stan_model(file= "gplvm.stan")

# g = spec(1,1)

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


##### Data ###
lambda1= 1E3
lambda2= 1E3

input_dat <- list(n=n, r=r,p=p, w=w, g=g, Y=as.matrix(y),lambda1= lambda1, lambda2 =lambda2)
ss_fit <- sampling(ss_model, data = input_dat, iter = 2000, chains = 1)


# save(ss_fit,file="gplvm.fit")
# quit()
# load("gplvm.fit")


post_x<- extractPosterior("x", n,"ss_fit")


# post_x[,1]<- 0.001


acf(post_x[1001:2000,1])

ts.plot(ss_fit@sim$samples[[1]]$sigma)

post_alpha<- extractPosteriorMat("alpha", r,1,"ss_fit")
ts.plot(post_alpha)
post_beta<- extractPosteriorMat("beta", r,1,"ss_fit")
post_beta


getFitted<-function(idx){
  x1 = post_x[idx,]
  cos(x1 %*% t(w)) %*% post_alpha[idx,] + sin(x1 %*% t(w)) %*% post_beta[idx,]
}

fittedCurve<- t(sapply(c(1:2000),getFitted))

# dim(fittedCurve)

ts.plot(fittedCurve[1501:2000,1])
ts.plot(fittedCurve[1501:2000,2])
ts.plot(fittedCurve[1501:2000,3])

plot(post_x[1001,],fittedCurve[1001,])
plot(post_x[1002,],fittedCurve[1002,])

plot(post_x[1500,],fittedCurve[1500,])
plot(post_x[2000,],fittedCurve[2000,])

ts.plot(post_x[,4])
