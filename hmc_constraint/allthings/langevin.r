setwd("~/git/empiricalTensor/hmc_constraint/")
require("rstiefel")

N= 2
d = 1
m = 10000

# F0<-matrix(rnorm(N*d),N,d)
F0 <- matrix(c(1,1),N,d)

y = t(sapply(c(1:m), function(x) rmf.matrix(F0)))

plot(y,xlim=c(-1,1),ylim=c(-1,1))



#loss function
loss = function(X, F=F0, lam =lambda){
	 - sum(diag(t(F)%*%X)) +  lam * sum((diag(t(X)%*%X)-1)^2) 
}

computeDeri = function(X, F=F0,lam =lambda){
  - F + 4 * lam * as.numeric(t(X)%*%X-1)* X
}


# (loss(y0) - loss(y0-c(1E-8,0)))/1E-8
# computeDeri(y0)


updateX <-function(q_0,eps = 0.1, L =10,steps = 1000,tuning = T, ideal_AR =0.234, burnin =T, microsteps = 100){
  #update x
  trace_x<- numeric()
  
  for(j in 1:(floor(steps/microsteps))){
  	accept = 0
  	for(i in 1:microsteps){

  		m_x = rep(1,N*d)
  		p_0 = matrix( rnorm(N*d) * sqrt(m_x), N, d)
  		q_x = q_0
  		p_x = p_0
  		p_x = p_x - eps/2 * computeDeri(q_x)
  		for(i in 1:L){
  			q_x = q_x + eps* p_x/m_x
  			p_x = p_x - eps * computeDeri(q_x)
  		}
  		q_x = q_x + eps* p_x/m_x
  		p_x = p_x - eps/2 * computeDeri(q_x)

  		cur_H = loss(q_0) + sum(p_0^2/m_x/2)
  		new_H = loss(q_x) + sum(p_x^2/m_x/2)

      diff_H = cur_H-new_H
      if(is.finite(diff_H)){
  		if(log(runif(1))< diff_H){
       q_0 = q_x
       if(tuning) { accept = accept+1}
      } 
      }
     if((!tuning) & (!burnin)) trace_x<- rbind(trace_x,c(q_0))
   }
    #acceptance rate
    AR = accept / microsteps
    #reduce eps if AR too low, else increase eps
    if(tuning)  {
    	eps = eps* exp(AR-ideal_AR)
    	print(AR)
    }
    print(c("step: ",j, ", -log-L: ", loss(q_0)))
  }
  return(list(x=q_0,trace_x=trace_x, eps=eps))
}

lambda = 1E3
x0 = matrix( c(0,1),N,d)

L= 100

runX = updateX(q_0=x0, L = L,eps= 0.1,steps = 10000, tuning = TRUE,ideal_AR = 0.6, microsteps = 100)

runX = updateX(q_0=runX$x,L = L, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = T)

runX = updateX(q_0=runX$x,L = L, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)



# ts.plot(runX$trace_x[,1])
acf(runX$trace_x)


plot(runX$trace_x[1:1000,],type="l")


# shuffled = sample(1:10000,replace=T)
# runX$trace_x = runX$trace_x[shuffled,]

# runX$trace_x= runX$trace_x[c(1:1000)*10,]
rowSums(runX1$trace_x^2)

X_D = runX1$trace_x / sqrt(rowSums(runX1$trace_x^2))

H_x = apply(runX1$trace_x,1,loss)
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

acf(runX1$trace_x)
acf(trace_X_D)
accept/10000

ts.plot(trace_X_D)

plot(y, xlim=c(-1,1),ylim=c(-1,1))
plot(trace_X_D, xlim=c(-1,1),ylim=c(-1,1))



## lam on mixing


lambda = 1E3
runX = updateX(q_0=x0, L = 100,eps= 0.05,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100,  burnin = F)
runX1 = updateX(q_0=runX$x,L = 100, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)

plot(runX1$trace_x[1:1000,],type="l")

mean(1E3* (rowSums(runX1$trace_x^2) -1 ))



lambda = 1E4
runX = updateX(q_0=x0, L = 100,eps= 0.05,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100,  burnin = F)
runX2 = updateX(q_0=runX$x,L = 100, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)

plot(runX2$trace_x[1:1000,],type="l")

mean(1E4* (rowSums(runX2$trace_x^2) -1 ))


lambda = 1E5
runX = updateX(q_0=x0, L = 100,eps= 0.05,steps = 10000, tuning = T,ideal_AR = 0.6, microsteps = 100,  burnin = F)
runX3 = updateX(q_0=runX$x,L = 100, eps=runX$eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)

plot(runX3$trace_x[1:1000,],type="l")

mean(1E5* (rowSums(runX2$trace_x^2) -1 ))

require("ggplot2")


traceX = rbind( runX1$trace_x[1:1000,], runX2$trace_x[1:1000,], runX3$trace_x[1:1000,])

df = data.frame( "x1" = traceX[,1], "x2"= traceX[,2], "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=1000 )))

pdf("../draft/unit_circle_path.pdf",8,3)
ggplot(data=df, aes(x=x1, y= x2))+ geom_path(size=.1)+ theme_bw()+facet_grid(~lambda)
dev.off()


acf1 = acf(runX1$trace_x[,1],lag.max=40,plot=F)$acf[,,1]
acf2 = acf(runX2$trace_x[,1],lag.max=40,plot=F)$acf[,,1]
acf3 = acf(runX3$trace_x[,1],lag.max=40,plot=F)$acf[,,1]


df = data.frame( "ACF" = c(acf1,acf2,acf3),"Lag"=rep(c(0:40),3), "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=41 )))

pdf("../draft/unit_circle_acf.pdf",8,3)
ggplot(data=df, aes(x=Lag, y= ACF))+ geom_path(size=.75)+ theme_bw()+facet_grid(~lambda)
dev.off()

df = data.frame( "relaxation" = c(abs(rowSums(runX1$trace_x^2)-1), abs(rowSums(runX2$trace_x^2)-1),abs(rowSums(runX3$trace_x^2)-1)), "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=10000 )))

pdf("../draft/unit_circle_violation.pdf",8,3)
ggplot(data=df, aes(relaxation))+ geom_histogram(binwidth = 0.001)+ theme_bw()+facet_grid(~lambda)
dev.off()


hist(c(abs(rowSums(runX1$trace_x^2)-1)))

hist(c(abs(rowSums(runX2$trace_x^2)-1)))

hist(c(abs(rowSums(runX3$trace_x^2)-1)))


max(c(abs(rowSums(runX3$trace_x^2)-1)))


(10^3) = 1/ 2*x^2

sqrt(1/(10^3)*2)/4

runX1$eps
runX2$eps
runX3$eps



contourX = seq(0,1,length.out=100)
contourY = sqrt(1- contourX^2)

eps = 0.02
lambda=1E3
runX1 = updateX(q_0=c(sqrt(2)/2,sqrt(2)/2),L = 1, eps=eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)



eps = 0.003
lambda=1E4
runX2 = updateX(q_0=c(sqrt(2)/2,sqrt(2)/2),L = 1, eps=eps,steps = 10000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)



eps = 0.0002
lambda=1E5
runX3 = updateX(q_0=c(sqrt(2)/2,sqrt(2)/2),L = 1, eps=eps,steps = 20000, tuning = FALSE,ideal_AR = 0.6, microsteps = 100, burnin = F)





plot(runX1$trace_x[1:1000,],type="l",xlim=c(0.4,0.9),ylim=c(0.4,0.9), col="red")
lines(contourX, sqrt(0.92- contourX^2),col="blue")
lines(contourX, sqrt(1.08- contourX^2),col="blue")

plot(runX2$trace_x[1:1000,],type="l",xlim=c(0.4,0.9),ylim=c(0.4,0.9), col="red")
lines(contourX, sqrt(0.97- contourX^2),col="blue")
lines(contourX, sqrt(1.03- contourX^2),col="blue")

plot(runX3$trace_x[1:1000,],type="l",xlim=c(0.4,0.9),ylim=c(0.4,0.9), col="red")
lines(contourX, sqrt(0.995- contourX^2),col="blue")
lines(contourX, sqrt(1.005- contourX^2),col="blue")


steps = 100

idx = c(1:steps)*10

traceX = rbind( runX1$trace_x[idx,], runX2$trace_x[idx,], runX3$trace_x[idx,])

d1 = data.frame( "x1" = traceX[,1], "x2"= traceX[,2], "lambda"= as.factor(rep(c("1,000","10,000","100,000"),each=steps )))


contourX = seq(0,1,length.out=steps)
d2 = data.frame("x"= contourX, "lb"= c(sqrt(0.92- contourX^2) ,sqrt(0.97- contourX^2), sqrt(0.995- contourX^2)) , "ub"=c(sqrt(1.08- contourX^2) ,sqrt(1.03- contourX^2), sqrt(1.005- contourX^2))  )

df = data.frame(d1,d2)



pdf("../draft/unit_circle_100steps.pdf",8,3)
ggplot(data=df)+ geom_path( aes(x=x1, y= x2),size=0.5,col="red") + 
 geom_line( aes(x=x, y= lb), size=0.5,linetype=2,alpha=0.2)+ 
 geom_line( aes(x=x, y= ub), size=0.5,linetype=2,alpha=0.2)+ theme_bw()+facet_grid(~lambda)
dev.off()




plot(rt(1000,df=1),rnorm(1000,0,10000))