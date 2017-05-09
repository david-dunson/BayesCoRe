setwd("~/work/git/empiricalTensor/test/")
setwd("D:/work/git/empiricalTensor/test/")



require("rstan")
rstan_options(auto_write = TRUE)


#######################################################
# spike-and-slab
ss_model = stan_model(file= "spike_and_slab3.stan")



extractPosterior<-function(varname, dimen, stan_fit){
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
}

##### Data ###
N = 40
p = 10
p_eff = 3

X = matrix(rnorm(N*p),N)
beta = rep(0,p)
beta[1:p_eff] = c(1:p_eff)*3
Z = rnorm(N, X%*%beta,1)

# beta

d= eigen(t(X)%*%X)$values[1] + 1

# d = 128

beta0 = c(solve( t(X)%*%X + diag(0.1,p), t(X)%*%Z))

toy_dat <- list(N=N, p=p, X=X,Z=Z, d=d)
init<- list(list(beta= beta0, y=rnorm(p),prob= rep(0.5,p)))

ss_fit <- sampling(ss_model, data = toy_dat, iter = 20000, chains = 1)


post_beta<- extractPosterior("beta", p,"ss_fit")
ts.plot(post_beta[,1:3])
colMeans(post_beta)
beta



includ_prob<- extractPosterior("includ_prob", p,"ss_fit")
plot(colMeans(includ_prob))
ts.plot(includ_prob)


prob<- extractPosterior("prob", p,"ss_fit")
plot(colMeans(prob))
ts.plot(includ_prob)


acf(post_beta[,1:5])




#
require("monomvn")
hs_fit= blasso(X,Z, case="hs")
plot(colMeans(hs_fit$beta))




###
Xbeta= X %*%diag(beta);
d= eigen(t(Xbeta)%*%Xbeta,only.values = T)$values

600 - d

eigen(diag(500,p) - t(Xbeta)%*%Xbeta)$values

print(d)
# 
# solve(diag(d,p)-t(Xbeta)%*%Xbeta) 
# 
# XbetaY= Xbeta %*% y;
# ZXbeta= t(Xbeta)%*% Z;
# 
# part1 = t(y)%*%y/d
# temp =  (Xbeta %*% t(Xbeta) /d - diag(1, N));
# part2 = - t(XbetaY) %*%solve ( temp , XbetaY) / d/ d;
# 
# part1 + part2
# t(y)%*% solve( - t(Xbeta)%*%Xbeta + diag(d,p), y)
# 
# det(temp)
# det( (- t(Xbeta)%*%Xbeta + diag(d,p)) /d )
