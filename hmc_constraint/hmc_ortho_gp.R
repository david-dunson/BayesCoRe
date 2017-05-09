setwd("~/git/empiricalTensor/hmc_constraint/")
# setwd("C:/Users/leo/git/empiricalTensor/hmc_constraint/")

install.packages("rstan")

require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#######################################################

extractPosterior<-function(varname, dimen, stan_fit){
  if(dimen==1){
    eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"`",sep="")))
  }else{
  sapply(c(1:dimen), function(i)  eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i, "]`",sep=""))))
  }
}

#model
ss_model = stan_model(file= "ortho_gp.stan")

##### Data ###
N = 20
d = 3
p = 5

x = runif(N)



g= matrix(0,N,d)

g[,1] = exp(-x^3/2)
g[,2] = sin(x*3)
g[,3] = cos(x*8)

eta = matrix(rnorm(d*p),d,p)

mu = g%*%eta

order_x =order(x)
plot(rep(x,p), mu, type="n")
for(i in 1:p){
	lines(x[order_x],mu[order_x,i])
}

y = rnorm(N*p) * 0.1 + mu


lambda=1E6

input_dat <- list(N=N, d=d,p=p, x=x, y=y, jitter=1E-3,lambda = lambda)
init<- list(list(p= rep(1/d,d), mu1=abs(rnorm(d)),mu2=abs(rnorm(d)),sigma=0.1))

ss_fit <- sampling(ss_model, data = input_dat, iter = 10000, chains = 1, algorithm = "NUTS")

save(ss_fit,file="gp.RDa")
quit()


sampling_idx<- c(5001:10000)


post_p<- extractPosterior("p", d,"ss_fit")
post_mu1<- extractPosterior("mu1", d,"ss_fit")
post_mu2<- extractPosterior("mu2", d,"ss_fit")

# post_mu[post_p<0.05]<- NA

# save(ss_fit,file="hmc_dp.RDa")
# save(ss_fit,file="hmc_dp2.RDa")


pdf("../draft/fmm_mu1_hmc.pdf",6,3)
df = data.frame( step =c(1:5000), mu1=c(post_mu1[sampling_idx,]), component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=mu1, col=component),alpha =0.8,size=.8) + theme_bw() +ylim(0,6)
dev.off()




pdf("../draft/fmm_mu2_hmc.pdf",6,3)
df = data.frame( step =c(1:5000), mu2=c(post_mu2[sampling_idx,]),component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=mu2, col=component),alpha =0.8,size=.8) + theme_bw()+ylim(0,6)
dev.off()

pdf("../draft/fmm_w_hmc.pdf",6,3)

post_p = post_p/rowSums(post_p)

sum(post_p[,2]<post_p[,3])

df = data.frame( step =c(1:5000), w=c(post_p[sampling_idx,]), component =as.factor(rep(c(1:3), each=5000)))

ggplot(data=df) + geom_path( aes(x=step, y=w, col=component),alpha =0.8,size=.8) + theme_bw() +ylim(0,0.8)
dev.off()




###### no ordering

lambda1=0
lambda2=1E3

# hist(y)


input_dat <- list(N=N, d=d,y1=y1,y2=y2, alpha=alpha, lambda1=lambda1, lambda2=lambda2)
init<- list(list(p= rep(1/d,d), mu1=abs(rnorm(d)),mu2=abs(rnorm(d)),sigma=0.1))


ss_fit2 <- sampling(ss_model, data = input_dat,init=init, iter = 10000, chains = 1, algorithm = "NUTS")


post_p<- extractPosterior("p", d,"ss_fit2")
post_mu1<- extractPosterior("mu1", d,"ss_fit2")
post_mu2<- extractPosterior("mu2", d,"ss_fit2")

# post_mu[post_p<0.05]<- NA



pdf("../draft/fmm_mu1_hmc_unordered.pdf",6,3)
df = data.frame( step =c(1:5000), mu1=c(post_mu1[sampling_idx,]), component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=mu1, col=component),alpha =0.8,size=.8) + theme_bw() +ylim(0,6)
dev.off()




pdf("../draft/fmm_mu2_hmc_unordered.pdf",6,3)
df = data.frame( step =c(1:5000), mu2=c(post_mu2[sampling_idx,]),component =as.factor(rep(c(1:3), each=5000)))
ggplot(data=df) + geom_path( aes(x=step, y=mu2, col=component),alpha =0.8,size=.8) + theme_bw()+ylim(0,6)
dev.off()

pdf("../draft/fmm_w_hmc_unordered.pdf",6,3)

post_p = post_p/rowSums(post_p)

sum(post_p[,2]<post_p[,3])

df = data.frame( step =c(1:5000), w=c(post_p[sampling_idx,]), component =as.factor(rep(c(1:3), each=5000)))

ggplot(data=df) + geom_path( aes(x=step, y=w, col=component),alpha =0.8,size=.8) + theme_bw() +ylim(0,0.8)
dev.off()

