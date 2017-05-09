setwd("~/git/empiricalTensor/hmc_constraint/")
# setwd("C:/Users/leo/git/empiricalTensor/hmc_constraint/")
# require("rstan")

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

##################################

require("rstiefel")

N= 10
d = 5
m = 100

F<-matrix(1,N,d)

X_list = lapply(c(1:m), function(x) rmf.matrix(F))
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

loss_function = function(F, X){
	sum(diag(t(F)%*%X))
}


#generate latent Xprime
Xprime_list = lapply(c(1:m), function(x) rmf.matrix(F))
XprimeSum= matrix(0, N, d)
for(i in 1:m){
	XprimeSum  = XprimeSum + Xprime_list[[i]]
}


curU = function(F) { loss_function(F,Xsum) - loss_function(F,XprimeSum)}

computeDeriX = function(F){
	Xsum - XprimeSum
}



eps = 0.01
L = 1

m_x = rep(1,N*d)
p_0 = matrix( rnorm(N*d) * sqrt(m_x), N, d)
q_x = F
p_x = p_0
p_x = p_x - eps/2 * computeDeriX(q_x)
for(i in 1:L){
	q_x = q_x + eps* p_x/m_x
	p_x = p_x - eps * computeDeriX(x)
}
q_x = q_x + eps* p_x/m_x
p_x = p_x - eps/2 * computeDeriX(x)

cur_H = curU(F) + sum(p_0^2/m_x/2)
new_H = curU(q_x) + sum(p_x^2/m_x/2)
      
if(log(runif(1))< (cur_H-new_H)){
    F = q_x
    if(tuning)  accept = accept+1
}
      if((!tuning) & (!burnin)) trace_x<- rbind(trace_x,x)