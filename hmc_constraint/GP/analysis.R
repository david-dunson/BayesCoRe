setwd("~/git/empiricalTensor/hmc_constraint/GP")


#######################################################

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

load(file="gp1.RDa")

##### Data ###
N = 20
d = 5
p = 5

ss_fit = data$ss_fit
x = data$x
g = data$g
y = data$y

post_g = extractPosteriorMat("g",N,d,"ss_fit")


n_steps = nrow(post_g)
sampling_idx<- c((n_steps*3/4+1):n_steps)

par(mfrow=c(1,1))
acf(post_g[sampling_idx,1],lag.max=100)
ts.plot(post_g[sampling_idx,1])
ts.plot(post_g[sampling_idx,2])
ts.plot(post_g[sampling_idx,3])
ts.plot(post_g[sampling_idx,4])
ts.plot(post_g[sampling_idx,5])
ts.plot(post_g[sampling_idx,N+2])


# dev.off()

order_x=order(x)


plot_j= function(j){
	curve1 = matrix(0,N,100)

	for(i in 1:100){
	  curve1[,i]=post_g[((n_steps/2)+i*10),((j-1)*N+1):(j*N)]
	}

	plot(range(x), range(c(curve1)),type="n")
	for(i in 1:100){
		lines(x[order_x], curve1[order_x,i])
	}
}

par(mfrow=c(2,5))
plot_j(1)
plot_j(2)
plot_j(3)
plot_j(4)
plot_j(5)

plot(x[order_x],g[order_x,1],type="l")
plot(x[order_x],g[order_x,2],type="l")
plot(x[order_x],g[order_x,3],type="l")


curve1 = matrix(0,N,d)
for(j in 1:d){
  curve1[,j]=post_g[n_steps/2,((j-1)*N+1):(j*N)]
}
t(curve1)%*%curve1

post_eta = extractPosteriorMat("eta",d,p,"ss_fit")

acf(post_eta[sampling_idx,d+1],lag.max=100)
acf(post_eta[sampling_idx,1],lag.max=100)

trace_mu <- sapply(sampling_idx, function(i){

  g_est =post_g[i,]
  g_est = matrix(g_est, N,d)

  eta_est = post_eta[i,]
  eta_est = matrix(eta_est, d,p)

  mu_est = g_est %*% eta_est
  c(mu_est)
})

trace_mu= t(trace_mu)

dim(trace_mu)
acf(trace_mu[,2])

acf(trace_mu[,2])
acf(trace_mu[,3])
acf(trace_mu[,4])
acf(trace_mu[,N+1])
acf(trace_mu[,2*N+1])
acf(trace_mu[4000:5000,2*N+1])


mu = matrix(trace_mu[1,],N,p)

par(mfrow=c(1,p))
for(i in 1:p){
	plot(x[order_x], y[order_x,i])
	lines(x[order_x], mu[order_x,i])
}

post_tau = extractPosterior("tau",d,"ss_fit")

post_tau[n_steps,]


post_rho = extractPosterior("rho",d,"ss_fit")

post_rho[n_steps,]

ts.plot(post_rho[,5])