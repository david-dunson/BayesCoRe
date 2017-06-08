# setwd("c:/Users/leo/git/empiricalTensor/gplvm/")
setwd("~/git/empiricalTensor/gplvm/")

n= 50
r = 50
p = 1

Sigma = rep(0.01, p)

w = seq(0,pi,length.out = r+1)
w = w[-r]
x<- runif(n)
# x<- c(1:n)
x[1]=0



spec<- function(rho, phi){
  # phi*exp(-1/4 * rho^2 * w^2)
  phi /(1+ rho^2 * w^2)
}

g = spec(1,10)

min(g)

g_theta = rep(g,2)

theta = matrix( rnorm(2*r*p)*sqrt(g_theta), 2*r)
# theta[,1] = abs(theta[,1])

xw = outer(x,w,"*")  
C = cbind( cos(xw), sin(xw))


Y = C%*%theta + matrix(rnorm(n*p),n)%*%diag(sqrt(Sigma),p)

plot(c(0,1),range(Y),type="n")
for(i in 1:p)lines(x[order(x)],Y[order(x),i])

# theta[,3]= -theta[,3]
# i=3
# plot(x[order(x)],Y[order(x),i])


# deriC0 = cbind( -sin(xw), cos(xw))
# deriC = t(t(deriC0)*w)
##test
# 
# theta_sum2=0
# theta_y_sum = 0
# for(i in 1:p){
#   theta_sum2 = theta_sum2 + theta[,i]%*%t(theta[,i])/Sigma[i]
#   theta_y_sum =  theta_y_sum+ theta[,i]%*%t(Y[,i])/Sigma[i]
# }
# 
# dim(theta_sum2)
# dim(theta_y_sum)
# 
# Chat = solve(theta_sum2 + diag(0.001,nrow = 2*r), theta_y_sum)
# 
# image(Chat,zlim = c(-1,1))
# image(C,zlim = c(-1,1))
# 
# Chat[2,]
# C[2,]
#not very good
######



source("./gplvm_hmc.r")
# computeDeriX(x) -deri_x
# deri_x = x
# for(i in 1:length(x)){
#   x1 = x
#   x1[i] = x[i]+1E-8
#   deri_x[i] = (computeU(x=x1)-computeU(x=x))/1E-8
# }

x1<- -runif(n)
x1[1]=0
# plot(x1,Y[,1])
# plot(x1,Y[,2])


runX = updateX(x=x1,L = 1,steps = 2000, tuning = T,ideal_AR = 0.6)

# runX = updateX(x=runX$x,L = 1,steps = 2000, tuning = T,ideal_AR = 0.5)

runX = updateX(x=runX$x,eps = runX$eps,L = 1,steps = 20000, tuning = F, ideal_AR=0.6)
runX = updateX(x=runX$x,eps = runX$eps,L = 1,steps = 2000, tuning = F,burnin = F)

computeXslope(runX$x)

plot(runX$x,Y[,1])
plot(runX$x,Y[,3])

#check on mixing
# runX = updateX(x=x,L = 1,steps = 2000, tuning = T,ideal_AR = 0.8)
# runX = updateX(x=x,eps = runX$eps,L = 1,steps = 2000, tuning = F,burnin = F)

ts.plot(runX$trace_x[,2])
acf(runX$trace_x[,2], lag.max = 100)
############

plot(c(0,1),range(Y),type="n")
for(i in 1:p)lines(x[order(x)],Y[order(x),i])
for(i in 1:p)lines(runX$x[order(runX$x)],Y[order(runX$x),i],col="red")

computeU(x,theta,Sigma)
# computeU(x1,theta,Sigma)
computeU(runX$x,theta,Sigma)
