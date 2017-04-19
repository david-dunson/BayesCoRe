options(max.print=100)

setwd("~/git/empiricalTensor/hmc_constraint/")


N= 3
d = 1
m = 10000

require("gtools")
require("rgl")
require("coda")

F0 = c(0.1,0.1,0.1)

y = rdirichlet(m,F0)


plot3d(y)


#loss function
loss = function(X, F=F0, lam =lambda){
	if(sum(X<0)>0 | sum(X>1)>0 ){
		l= Inf
	}else{
		l = sum( -(F-1)*log(X) ) +  lam * (sum(X)-1)^2
	}
	l
}

computeDeri = function(X, F=F0,lam =lambda){
  	-(F-1)/X + 2* lam * (sum(X)-1)*rep(1,length(X))
}


x0 = c(0.3,0.5,0.2)
(loss(x0) - loss(x0-c(1E-8,0,0)))/1E-8
(loss(x0) - loss(x0-c(0,1E-8,0)))/1E-8
(loss(x0) - loss(x0-c(0,0,1E-8)))/1E-8
computeDeri(x0)


updateX <-function(q_0,eps = 0.1, L =10,steps = 1000,tuning = T, ideal_AR =0.234, burnin =T, microsteps = 100, useSigmaP = FALSE,SigmaP = 0){
  #update x
  trace_x<- numeric()

  if(!useSigmaP){
  	SigmaP =diag( 1, length(q_0))
  }

  invSigmaP = solve(SigmaP)
  cholSigmaP = t(chol(SigmaP))
  
  for(j in 1:(floor(steps/microsteps))){
  	accept = 0
  	for(i in 1:microsteps){

  		p_0 = matrix( cholSigmaP %*% rnorm(N*d) , N, d)
  		q_x = q_0
  		p_x = p_0
  		p_x = p_x - eps/2 * computeDeri(q_x)
  		for(i in 1:L){
  			q_x = q_x + eps* invSigmaP %*% p_x
  			p_x = p_x - eps * computeDeri(q_x)
  		}
  		q_x = q_x + eps* invSigmaP %*% p_x
  		p_x = p_x - eps/2 * computeDeri(q_x)

  		cur_H = loss(q_0) + c( t(p_0) %*%invSigmaP  %*% p_0/2) #sum(p_0^2/m_x/2)
  		new_H = loss(q_x) + c( t(p_x) %*%invSigmaP  %*% p_x/2)#sum(p_x^2/m_x/2)

      diff_H = cur_H-new_H
      if(is.finite(diff_H)){
  		if(log(runif(1))< diff_H){
       q_0 = q_x
       { accept = accept+1}
      } 
      }
     if((!tuning) & (!burnin)) trace_x<- rbind(trace_x,c(q_0))
   }
    #acceptance rate
    AR = accept / microsteps
    #reduce eps if AR too low, else increase eps
    if(tuning)  {
    	eps = eps* exp(AR-ideal_AR)
    }
    print(c("step: ",j, ", -log-L: ", loss(q_0), " AR:",   AR))
  }
  return(list(x=q_0,trace_x=trace_x, eps=eps))
}

lambda = 1E3
x0 = rep(1/length(F0),length(F0))

# F0 = c(100,20,2)

L = 100


F0 = c(20,10,10)



x0 = c(0.1,0.1,0.1)
runX = updateX(q_0=x0, L = L,eps= 1E-1,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100)

runX = updateX(q_0=runX$x,L = L, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)


acf(runX$trace_x[,1:3])

plot3d(runX$trace_x,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1))





# ts.plot(runX$trace_x[,1])


# thin = 10 

# runX$trace_x = runX$trace_x[ c(1:(10000/thin))*thin,]


y = rdirichlet(m,F0)

plot3d(y,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1))


# shuffled = sample(1:10000,replace=T)
# runX$trace_x = runX$trace_x[shuffled,]

rowSums(runX$trace_x)
X_D=runX$trace_x /rowSums(runX$trace_x)

H_x = apply(runX$trace_x,1,loss)
H_xd = apply(X_D,1,loss)


trace_X_D=numeric()

X_D0= X_D[1,]
i = 1
cur_diff = (H_x[i] - H_xd[i])   # exp( - (H_xd - H_x) )
trace_X_D = rbind(trace_X_D, X_D0)

accept =0 
for(i in 2:nrow(X_D)){
  mh = (H_x[i] - H_xd[i])   - cur_diff
  if(log( runif(1))< mh){
    X_D0 = X_D[i,]
    accept = accept +1
  }
  trace_X_D = rbind(trace_X_D, X_D0)
}

acf(runX$trace_x)
acf(trace_X_D)
accept/10000

acf(trace_X_D)

plot3d(y,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1))
plot3d(trace_X_D,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1))

ts.plot(trace_X_D[,3])

# ## lam on mixing


# lambda = 1E3
# runX = updateX(q_0=x0, L = 100,eps= 0.05,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100,  burnin = F)
# runX1 = updateX(q_0=runX$x,L = 100, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)

# plot(runX1$trace_x[1:1000,],type="l")





# lambda = 1E4
# runX = updateX(q_0=x0, L = 100,eps= 0.05,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100,  burnin = F)
# runX2 = updateX(q_0=runX$x,L = 100, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)

# plot(runX2$trace_x[1:1000,],type="l")


# lambda = 1E5
# runX = updateX(q_0=x0, L = 100,eps= 0.05,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100,  burnin = F)
# runX3 = updateX(q_0=runX$x,L = 100, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)

# plot(runX3$trace_x[1:1000,],type="l")


# require("ggplot2")


# traceX = rbind( runX1$trace_x[1:1000,], runX2$trace_x[1:1000,], runX3$trace_x[1:1000,])

# df = data.frame( "x1" = traceX[,1], "x2"= traceX[,2], "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=1000 )))

# pdf("../draft/unit_circle_path.pdf",8,3)
# ggplot(data=df, aes(x=x1, y= x2))+ geom_path(size=.1)+ theme_bw()+facet_grid(~lambda)
# dev.off()


# acf1 = acf(runX1$trace_x[,1],lag.max=40,plot=F)$acf[,,1]
# acf2 = acf(runX2$trace_x[,1],lag.max=40,plot=F)$acf[,,1]
# acf3 = acf(runX3$trace_x[,1],lag.max=40,plot=F)$acf[,,1]


# df = data.frame( "ACF" = c(acf1,acf2,acf3),"Lag"=rep(c(0:40),3), "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=41 )))

# pdf("../draft/unit_circle_acf.pdf",8,3)
# ggplot(data=df, aes(x=Lag, y= ACF))+ geom_path(size=.75)+ theme_bw()+facet_grid(~lambda)
# dev.off()

# df = data.frame( "violation" = c(abs(rowSums(runX1$trace_x^2)-1), abs(rowSums(runX2$trace_x^2)-1),abs(rowSums(runX3$trace_x^2)-1)), "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=10000 )))

# pdf("../draft/unit_circle_violation.pdf",8,3)
# ggplot(data=df, aes(violation))+ geom_histogram(binwidth = 0.001)+ theme_bw()+facet_grid(~lambda)
# dev.off()


# hist(c(abs(rowSums(runX1$trace_x^2)-1)))

# hist(c(abs(rowSums(runX2$trace_x^2)-1)))

# hist(c(abs(rowSums(runX3$trace_x^2)-1)))


# max(c(abs(rowSums(runX3$trace_x^2)-1)))


# (10^3) = 1/ 2*x^2

# sqrt(1/(10^3)*2)/4

# runX1$eps
# runX2$eps
# runX3$eps