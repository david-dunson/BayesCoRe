setwd("~/git/constrainedBayes/hmc_constraint/linear_inequal/")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#model
ss_model = stan_model(file= "inequality.stan")

sigma= 0.1


eps = 1E-3
lambda1=1E-10
n_steps = 2E4

##### Data 1 ###

N = 10
X<- matrix(rnorm(N*2),N)
theta = c(0.3,0.3)
y= rnorm(N, X%*%theta,sigma)

input_dat <- list(N=N, p=2, X = X, y=y,lambda1 = lambda1)
ss_fit <- sampling(ss_model, data = input_dat, 
                     iter = n_steps,
                     chains = 1, algorithm = "NUTS")

  
post_theta<- extract(ss_fit,"theta")
trace1 = post_theta[[1]]

##### Data 2 ###
N = 1E4
X<- matrix(rnorm(N*2),N)
theta = c(0.7,0.3)
y= rnorm(N, X%*%theta,sigma)

input_dat <- list(N=N, p=2, X = X, y=y,lambda1 = lambda1)
ss_fit <- sampling(ss_model, data = input_dat, 
                   iter = n_steps,
                   chains = 1, algorithm = "NUTS")

post_theta<- extract(ss_fit,"theta")
trace2 = post_theta[[1]]

#########
dat = list(trace1,trace2)
# save(dat,file="lin_inequ_fit.RDa")
load("lin_inequ_fit.RDa")
trace1 = dat[[1]]
trace2 = dat[[2]]

plot(trace1)
sum(rowSums(trace1)>1)

plot(trace2)
sum(rowSums(trace2)>1)

require("ggplot2")

#plot 1
df = data.frame(trace1)
colnames(df)<- c("theta1","theta2")


df2 = df

x1 = seq(min(df$theta1),max(df$theta1),length.out = 30)
x2 = seq(min(df$theta2),max(df$theta2),length.out = 30)
df2 = expand.grid(x1,x2)

z = apply(df2,MARGIN = 1, function(x){exp(- sum((x-c(0.3,0.3))^2)/2)})

df2 = cbind(df2,z)
colnames(df2)<- c("theta1","theta2","z")
df2$z[df2$theta1+df2$theta2>1]<-NA


pdf("../../draft/linear_inequal_1.pdf",5,3)
ggplot(df, aes(x=theta1, y=theta2))+
  geom_point(data=df,alpha=0.05, aes(x=theta1, y=theta2))+
  geom_contour(data=df2, aes(x=theta1, y=theta2,z=z))+
  geom_abline(slope = -1,intercept = 1,linetype=2)+
geom_vline(xintercept = 0,linetype=2)+ geom_hline(yintercept = 0,linetype=2)
dev.off()


#plot 2
df = data.frame(trace2)
colnames(df)<- c("theta1","theta2")

df2 = df

x1 = seq(min(df$theta1),max(df$theta1),length.out = 30)
x2 = seq(min(df$theta2),max(df$theta2),length.out = 30)
df2 = expand.grid(x1,x2)

z = apply(df2,MARGIN = 1, function(x){exp(- sum((x-c(0.7,0.3))^2)/0.1)})

df2 = cbind(df2,z)
colnames(df2)<- c("theta1","theta2","z")
df2$z[df2$theta1+df2$theta2>1]<-NA

pdf("../../draft/linear_inequal_2.pdf",5,3)
ggplot(df, aes(x=theta1, y=theta2))+
geom_point(data=df,alpha=0.05, aes(x=theta1, y=theta2))+
geom_contour(data=df2, aes(x=theta1, y=theta2,z=z))+
  geom_abline(slope = -1,intercept = 1,linetype=2)
  # geom_vline(xintercept = 0,linetype=2)+ geom_hline(yintercept = 0,linetype=2)
dev.off()

