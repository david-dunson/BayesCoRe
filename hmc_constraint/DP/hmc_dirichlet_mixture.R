setwd("~/git/empiricalTensor/hmc_constraint/DP")
# setwd("C:/Users/leo/git/empiricalTensor/hmc_constraint/")


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
ss_model = stan_model(file= "dp_mixture.stan")

##### Data ###
N = 100
d = 3
p = c(0.6,0.3,0.1)
mu1<- c(1,3,3)
mu2<- c(5,3,5)


sigma=1
y1<- numeric(N)
y2<- numeric(N)

for(i in 1:N){
	c_i = rmultinom(1,1,p)
  y1[i] = rnorm(1, sum(c_i*mu1),sigma)
  y2[i] = rnorm(1, sum(c_i*mu2),sigma)
}

alpha =2


lambda1=1E6
lambda2=1E3

# hist(y1,breaks = 100)


input_dat <- list(N=N, d=d,y1=y1,y2=y2, alpha=alpha, lambda1=lambda1, lambda2=lambda2)
init<- list(list(p= rep(1/d,d), mu1=abs(rnorm(d)),mu2=abs(rnorm(d)),sigma=0.1))

ss_fit <- sampling(ss_model, data = input_dat,init=init, iter = 10000, chains = 1, algorithm = "NUTS")

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

#plot contour for mu's

mu1<- c(1,3,3)
mu2<- c(5,3,5)
sigmas = sigma^2/ (N*p)



mu1,sqrt(sigmas)

sigma=1
y1<- numeric(N)
y2<- numeric(N)

for(i in 1:N){
	c_i = rmultinom(1,1,p)
  y1[i] = rnorm(1, sum(c_i*mu1),sigma)
  y2[i] = rnorm(1, sum(c_i*mu2),sigma)
}


df <- expand.grid(mu1=seq(0,4,by=0.1), mu2=seq(2,6,by=0.1))


sqrt(sigmas)

logDen = function(x1,x2){
	den=numeric(d)
	for(i in 1:d){
		den[i] =  dnorm(x1, mu1[i], sqrt(sigmas[i]), log=T) + dnorm(x2, mu2[i], sqrt(sigmas[i]), log=T)
	}
	(max(den))
}



df$z <- apply(cbind(df$mu1,df$mu2),1,function(x)logDen(x[1],x[2]))

pdf("../../draft/fmm_mu_contour.pdf",4,4)
ggplot(aes(x=mu1, y=mu2, z=z), data = df) + geom_contour(binwidth=4)+
theme_bw()
dev.off()

