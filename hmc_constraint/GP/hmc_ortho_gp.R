setwd("~/git/empiricalTensor/hmc_constraint/GP")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



#model
ss_model = stan_model(file= "ortho_gp.stan")

##### Data ###
N = 20
d = 3
p = 5

x = runif(N)


g= matrix(0,N,d)

g[,1] = sin(x*16)/x
g[,2] = 30*cos(x*5*pi)*(1-x)
g[,3] = 15*cos(x*10*pi)

order_x =order(x)

eta = matrix(rnorm(d*p),d,p)
mu = g%*%eta
y = rnorm(N*p) * 0.1 + mu


plot(rep(x,3),c(g[,1:3] ), type="n")
for(i in 1:3){
  lines(x[order_x], g[order_x,i])
}


lambda1=1E6
lambda2=1E3

input_dat <- list(N=N, d=d,p=p, x=x, y=y, jitter=1E-5,lambda1 = lambda1, lambda2 = lambda2, a= 2, q= 2)

ss_fit <- sampling(ss_model, data = input_dat, iter = 10000, chains = 1, algorithm = "NUTS")


lambda1=0
lambda2=0

input_dat <- list(N=N, d=d,p=p, x=x, y=y, jitter=1E-5,lambda1 = lambda1, lambda2 = lambda2, a= 2, q= 2)

ss_fit2 <- sampling(ss_model, data = input_dat, iter = 10000, chains = 1, algorithm = "NUTS")

data = list("y"=y,"g"=g,"x"=x,"ss_fit"=ss_fit, "ss_fit2"=ss_fit2)

save(data,file="gp.RDa")

quit()

load(file="gp2.RDa")

sampling_idx<- c(8001:10000)


##### Data ###
N = 20
d = 3
p = 5

ss_fit = data$ss_fit
x = data$x
g = data$g
y = data$y
ss_fit2= data$ss_fit2

# trace plot with orthonormality
post_g = extractPosteriorMat("g",N,d,"ss_fit")
post_eta = extractPosteriorMat("eta",d,p,"ss_fit")

trace_mu <- sapply(sampling_idx, function(i){

  g_est =post_g[i,]
  g_est = matrix(g_est, N,d)

  eta_est = post_eta[i,]
  eta_est = matrix(eta_est, d,p)

  mu_est = g_est %*% eta_est
  c(mu_est)
})

trace_mu= t(trace_mu)


# trace plot without orthonormality
post_g2 = extractPosteriorMat("g",N,d,"ss_fit2")
post_eta2 = extractPosteriorMat("eta",d,p,"ss_fit2")

trace_mu2 <- sapply(sampling_idx, function(i){

  g_est =post_g2[i,]
  g_est = matrix(g_est, N,d)

  eta_est = post_eta2[i,]
  eta_est = matrix(eta_est, d,p)

  mu_est = g_est %*% eta_est
  c(mu_est)
})

trace_mu2= t(trace_mu2)




acf(post_g[sampling_idx,N*4+1],lag.max=100)
ts.plot(post_g[sampling_idx,N*4+1])

acf(post_g2[sampling_idx,2],lag.max=100)
ts.plot(post_g2[sampling_idx,5])

acf(post_eta[sampling_idx,5],lag.max=100)
ts.plot(post_eta[sampling_idx,N+2])

acf(post_eta2[sampling_idx,5],lag.max=100)
ts.plot(post_eta2[sampling_idx,N+2])


ts.plot(post_g[sampling_idx,1])
ts.plot(post_g2[sampling_idx,1])

ts.plot(post_eta[sampling_idx,1])
ts.plot(post_eta2[sampling_idx,1])

eta_est = post_eta2[5001,]
eta_est = matrix(eta_est, d,p)
eta_est


1/rgamma(100,2*1,1)
1/rgamma(100,2*10^3,10^2)
1/rgamma(100,2*10^6,10^4)
1/rgamma(100,1*1^9,1^6)


order_x=order(x)

curve1 = matrix(0,N,50)
curve2 = matrix(0,N,50)

for(i in 1:50){
  curve1[,i]=post_g[(5001+i*50),1:N]
  curve2[,i]=post_g2[(5001+i*50),1:N]
}

plot(range(x), range(curve1),type="n",ylim=c(-0.5,0.5))
for(i in 1:50){
lines(x[order_x], curve1[order_x,i])
}

plot(range(x), range(post_g2[sampling_idx,]),type="n")
for(i in 1:50){
lines(x[order_x], curve2[order_x,i])
}
# plot(lowess(x, curve1[,1])$x, lowess(x, curve1[,1])$y)

# plot(loess.smooth(x, curve1[,1])$x, loess.smooth(x, curve1[,1])$y)


#other parameters


ts.plot(ss_fit@sim$samples[[1]]$"eta_sq[1]" )
ss_fit@sim$samples[[1]]$"eta_sq[2]" 
ss_fit@sim$samples[[1]]$"eta_sq[3]" 
ss_fit@sim$samples[[1]]$"eta_sq[4]" 
ss_fit@sim$samples[[1]]$"eta_sq[5]" 

ts.plot(ss_fit@sim$samples[[1]]$"rho_sq[1]" )
ts.plot(ss_fit@sim$samples[[1]]$"rho_sq[2]" )
ts.plot(ss_fit@sim$samples[[1]]$"rho_sq[3]" )
ts.plot(ss_fit@sim$samples[[1]]$"rho_sq[4]" )
ts.plot(ss_fit@sim$samples[[1]]$"rho_sq[5]" )



order_x=order(x)
plot(range(x), range(post_g2[sampling_idx,]),type="n",ylim=c(-0.5,0.5))
for(i in 1:50){
  lines(x[order_x], post_g2[(5001+i*100),1:N][order_x])
}



acf(post_g2[sampling_idx,1],lag.max=100)
acf(post_eta2[sampling_idx,1],lag.max=100)

acf(trace_mu2[,1],lag.max=40)












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

