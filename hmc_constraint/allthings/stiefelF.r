setwd("~/git/empiricalTensor/hmc_constraint/")
# setwd("C:/Users/leo/git/empiricalTensor/hmc_constraint/")
# require("rstan")

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

##################################

require("rstiefel")

N= 40
d = 2
m = 100

# F0<-matrix(rnorm(N*d),N,d)
F0 <- matrix((rnorm(N*d))*3,N,d)

X_list = lapply(c(1:m), function(x) rmf.matrix(F0))
# X = do.call("rbind", lapply(c(1:m), function(x) rmf.matrix(F)))

Xsum= matrix(0, N, d)

for(i in 1:m){
	Xsum  = Xsum + X_list[[i]]
}

###

# ss_model = stan_model(file= "intractable_problem.stan")

# MC_n= 100
# lambda=1E3

# F0 = matrix(rnorm(N*d),N,d)
# Y0 = do.call("rbind", lapply(c(1:MC_n), function(x) rmf.matrix(F0)))


# input_dat <- list(N=N, d=d,m=m, MC_n= MC_n, Xsum= Xsum,lambda=lambda)
# init<- list(list(F= F0, Y = Y0))

# ss_fit <- sampling(ss_model, data = input_dat, init = init,iter = 20000, chains = 1, algorithm = "NUTS")

# listIP = list("ss_fit"= ss_fit, "F" = F)
# save(listIP,file="ip.RDa")

# load(file="ip2.RDa")
# ss_fit = listIP$ss_fit
# F = listIP$F

# sampling_idx<- c(10001:20000)

# post_sample = do.call("cbind",ss_fit@sim$samples[[1]])

# post_sample = post_sample[sampling_idx,] 

# post_F = post_sample[, substr( colnames(post_sample), 1, 1) =="F"]
# post_F_tensor = array(post_F, dim=c(10000,N,d))



# image(post_F_tensor[1,,])
# F

##exchange algorithm####

#loss function

langevin = function(F, X){
	sum(diag(t(F)%*%X))
}

curU = function(F, XprimeSum) { -langevin (F,Xsum - XprimeSum)}

computeDeriX = function(F, XprimeSum){
	 - ( Xsum - XprimeSum) 
}

updateXprime <- function(F){
	#generate latent Xprime
	Xprime_list = lapply(c(1:m), function(x) rmf.matrix(F))
	XprimeSum= matrix(0, N, d)
	for(i in 1:m){
		XprimeSum  = XprimeSum + Xprime_list[[i]]
	}
	XprimeSum
}


updateF <-function(F,eps = 0.1, L =10,steps = 1000,tuning = T, ideal_AR =0.234, burnin =T, microsteps = 100){
  #update x
  trace_F<- numeric()
  
  for(j in 1:(floor(steps/microsteps))){
  	accept = 0
  	for(i in 1:microsteps){

  		XprimeSum = updateXprime(F)

  		m_x = rep(1,N*d)
  		p_0 = matrix( rnorm(N*d) * sqrt(m_x), N, d)
  		q_x = F
  		p_x = p_0
  		p_x = p_x - eps/2 * computeDeriX(q_x, XprimeSum)
  		for(i in 1:L){
  			q_x = q_x + eps* p_x/m_x
  			p_x = p_x - eps * computeDeriX(q_x, XprimeSum)
  		}
  		q_x = q_x + eps* p_x/m_x
  		p_x = p_x - eps/2 * computeDeriX(q_x, XprimeSum)

  		cur_H = curU(F, XprimeSum) + sum(p_0^2/m_x/2)
  		new_H = curU(q_x, XprimeSum) + sum(p_x^2/m_x/2)


  		if(log(runif(1))< (cur_H-new_H)){
			F = q_x
  			if(tuning) { accept = accept+1}
  		}
  		if((!tuning) & (!burnin)) trace_F<- rbind(trace_F,c(F))
  	}
    #acceptance rate
    AR = accept / microsteps
    #reduce eps if AR too low, else increase eps
    if(tuning)  {
    	eps = eps* exp(AR-ideal_AR)
    	print(AR)
    }
    print(c("step: ",j, ", -log-L: ", curU(F, XprimeSum)))
}
return(list(F=F,trace_F=trace_F, eps=eps))
}




F1<- matrix(rnorm(N*d),N,d)

runX = updateF(F=F1,L = 1,eps= 0.1,steps = 100, tuning = FALSE,ideal_AR = 0.6, microsteps = 10)

runX = updateF(F=runX$F,L = 10,eps=0.1,steps = 100, tuning = FALSE,ideal_AR = 0.6, microsteps = 1, burnin = FALSE)


acf(runX$trace_F[,1:5])


ts.plot(runX$trace_F[,1])


plot(colMeans(runX$trace_F),c(F0))
