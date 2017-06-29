setwd("~/git/constrainedBayes/hmc_constraint/unit_circle/")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#model
ss_model = stan_model(file= "unit_circle.stan")

##### Data ###
N = 2
sd= 1

theta = c(sqrt(3/4),sqrt(1/4))
y= cbind(rnorm(N,theta[1],sd),rnorm(N,theta[2],sd))

plot(y[,1],y[,2],xlim=c(-2,2),ylim=c(-2,2))
lines(theta[1],theta[2],col="red",type="p")
F=c(1,1)
n_steps = 2E4
eps = 1E-3

lambda1=1E-2


testUnitCircle<- function(lambda1,eps=0.001){
  input_dat <- list(N=N, F=F, y=y,lambda1 = lambda1)
  init = list(list("theta"=theta,"tau"=sd^2))
  ss_fit <- sampling(ss_model, data = input_dat, init=init,
                     iter = n_steps,
                     chains = 1, algorithm = "HMC",
                     control=list(
                       adapt_engaged=FALSE,
                       stepsize= eps,
                       int_time=100*eps
                       ))
  
  post_theta<- extract(ss_fit,"theta")
  post_theta[[1]]
}

trace1 = testUnitCircle(1E-3)
trace2 = testUnitCircle(1E-4, 1E-4)
trace3 = testUnitCircle(1E-5,1E-5)
trace4 = testUnitCircle(1E-8,1E-8)

plot(trace1[,1],trace1[,2])
plot(trace2[,1],trace2[,2])
plot(trace3[,1],trace3[,2])
plot(trace4[,1],trace4[,2])

acf(trace1[,1])
acf(trace2[,1])
acf(trace3[,1])
acf(trace4[,1])

require("coda")

effectiveSize(trace1[,1])
effectiveSize(trace2[,1])
effectiveSize(trace3[,1])
effectiveSize(trace4[,1])

abs(rowSums(trace1^2)-1)



df = data.frame( "relaxation" = c(abs(rowSums(trace1^2)-1),
                                  abs(rowSums(trace2^2)-1),
                                  abs(rowSums(trace3^2)-1)),
                 "lambda"= as.factor(rep(c("0.001","0.0001","0.00001"),each=10000 )))

df$lambda = factor(df$lambda,levels = c("0.001","0.0001","0.00001"))

pdf("../../draft/unit_circle_violation.pdf",8,3)
ggplot(data=df, aes(relaxation))+ geom_histogram(binwidth = 0.0001)+ theme_bw()+facet_grid(~lambda)
dev.off()





acf(post_theta[[1]][,1])
plot(post_theta[[1]][,1],post_theta[[1]][,2])
post_tau

# data = list("y"=y,"g"=g,"x"=x,"ss_fit"=ss_fit)
# save(data,file="gp1.RDa")

sampling_idx<- c((n_steps/2+1):n_steps)

post_U = extractPosteriorMat("U",N,d,"ss_fit")

acf(post_U[sampling_idx,1])
acf(post_U[sampling_idx,N+1])

ts.plot(post_U[sampling_idx,1])
ts.plot(post_U[sampling_idx,N+1])

post_eta = extractPosteriorMat("eta",d,d,"ss_fit")

ts.plot(post_eta[sampling_idx,1])
acf(post_eta[sampling_idx,2])


U1 = matrix(0,N,d)
for(j in 1:d){
  U1[,j]=post_U[n_steps/2,((j-1)*N+1):(j*N)]
}
t(U1)%*%U1



eta1 = matrix(0,d,d)
for(j in 1:d){
  eta1[,j]=post_eta[n_steps/2,((j-1)*d+1):(j*d)]
}

acf(post_eta[,1])
acf(post_eta[,2])
ts.plot(post_eta[,3])
acf(post_eta[,4])

