options(max.print=100)

setwd("~/git/empiricalTensor/gplvm/")

require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#################################################

r = 100
w = seq(0,20,length.out = r+1)
w = w[-1]
n = 100
p = 3

x<- runif(n)
x[1] = 0.1

spec<- function(rho, phi){
  phi*exp(-1/4 * rho^2 * w^2)
}

g = spec(0.1,10)
min(g)

alpha<- matrix(rnorm(r*p),r) * sqrt(g)
beta<- matrix(rnorm(r*p),r) * sqrt(g)
theta<- rbind(alpha, beta)

Sigma = rep(0.1,p)

xw =outer(x,w,"*")
X = cbind(cos(xw),sin(xw))

y = X%*%theta + matrix(rnorm(n*p),n)%*%diag(Sigma^0.5)


plot(x,y[,1])
plot(x,y[,2])
plot(x,y[,3])

#######################################################
ss_model = stan_model(file= "gplvm.stan")



extractPosterior<-function(varname, dimen, stan_fit){
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
}

##### Data ###

toy_dat <- list(n=n, r=r,p=p, w=w, g=g, Y=y)

init<- list(list(alpha =alpha, beta =beta, x=x[-1], sigma = Sigma))


ss_fit <- sampling(ss_model, data = toy_dat, init = init,iter = 20000, chains = 1)


save(ss_fit,file="gplvm.fit")

quit()

load("gplvm.fit")


post_x<- extractPosterior("x", 100,"ss_fit")

ts.plot(post_x[10000:20000,1:20])



post_sigma<- extractPosterior("sigma", 3,"ss_fit")

ts.plot(post_sigma[10000:20000,1:3])

post_sigma[,2]

