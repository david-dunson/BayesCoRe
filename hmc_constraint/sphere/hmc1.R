setwd("~/git/constrainedBayes/hmc_constraint/sphere/")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#model
ss_model = stan_model(file= "sphere.stan")

##### Data ###
N = 1
F = c(1,1,1)
F = F/ sqrt(sum(F^2))
lambda1=1E-2
sigma2 = 0.1
Sigma = diag(1,3)
Sigma[1,2]=Sigma[2,1]= 0.9
m=3

input_dat <- list(N=N, F=F, lambda1 = lambda1, sigma2=sigma2, Sigma=Sigma,m=m)

ss_fit <- sampling(ss_model, data = input_dat,
                     iter = 10000,
                     chains = 1, algorithm = "NUTS"
                     )
  
theta_vmf=extract(ss_fit,"theta_vmf",permuted=FALSE)
theta_fb=extract(ss_fit,"theta_fb",permuted=FALSE)
theta_t=extract(ss_fit,"theta_t",permuted=FALSE)

theta_vmf = theta_vmf[c(1:1000)*5,,]
theta_fb = theta_fb[c(1:1000)*5,,]
theta_t = theta_t[c(1:1000)*5,,]

proj<-function(x){
  t(apply(x,1,function(y){y/sqrt(sum(y^2))}))
}

vmf = proj(theta_vmf)
fb = proj(theta_fb)
t = proj(theta_t)


require("plot3D")

colnames(vmf)

scatter3D(x = vmf[,1],y=vmf[,2],z=vmf[,3],col="red", 
          xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1))

scatter3D(x = fb[,1],y=fb[,2],z=fb[,3],col="red", 
          xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1))

scatter3D(x = t[,1],y=t[,2],z=t[,3],col="red", 
          xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1))

vmf = data.frame(vmf)
fb = data.frame(fb)
t = data.frame(t)
colnames(vmf)<- c("x","y","z")
colnames(fb)<- c("x","y","z")
colnames(t)<- c("x","y","z")

pdf("sphere_vmf.pdf",6,6)
plot(vmf,xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1))
dev.off()
pdf("sphere_fb.pdf",6,6)
plot(fb,xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1))
dev.off()
pdf("sphere_t.pdf",6,6)
plot(t,xlim=c(-1,1),ylim=c(-1,1),zlim=c(-1,1))
dev.off()

scatterplot3d(vmf,angle = 60, box=FALSE)

effectiveSize(trace1[,1])
effectiveSize(trace2[,1])
effectiveSize(trace3[,1])
effectiveSize(trace4[,1])

mean(abs(rowSums(trace1^2)-1))
quantile(abs(rowSums(trace1^2)-1), c(0.025,0.975))
mean(abs(rowSums(trace2^2)-1))
quantile(abs(rowSums(trace2^2)-1), c(0.025,0.975))
mean(abs(rowSums(trace3^2)-1))
quantile(abs(rowSums(trace3^2)-1), c(0.025,0.975))







#compare with exact posterior
require("movMF")
require("transport")

exactY1 = movMF::rmovMF(1000,F+sum(y)/sd)
exactY2 = movMF::rmovMF(1000,F+sum(y)/sd)
tplan = transport(pp(exactY1[1:pick_n,]),pp(exactY2[1:pick_n,]),method = "shortsimplex")
wasserstein(pp(exactY1[1:pick_n,]),pp(exactY2[1:pick_n,]),p = 1,tplan = tplan)


plot(exactY1,xlim=c(0,1),ylim=c(-1,1))
plot(exactY2,xlim=c(0,1),ylim=c(-1,1))
plot(trace1,xlim=c(0,1),ylim=c(-1,1))

pick_n =1000
tplan = transport(pp(exactY1[1:pick_n,]),pp(exactY2[1:pick_n,]),method = "shortsimplex")
wasserstein(pp(exactY1[1:pick_n,]),pp(exactY2[1:pick_n,]),p = 1,tplan = tplan)

pick_n =1000
tplan = transport(pp(exactY1[1:pick_n,]),pp(trace1[1:pick_n,]),method = "shortsimplex")
wasserstein(pp(exactY1[1:pick_n,]),pp(trace1[1:pick_n,]),p = 1,tplan = tplan)

pick_n =1000
tplan = transport(pp(exactY1[1:pick_n,]),pp(trace2[1:pick_n,]),method = "shortsimplex")
wasserstein(pp(exactY1[1:pick_n,]),pp(trace2[1:pick_n,]),p = 1,tplan = tplan)

pick_n =1000
tplan = transport(pp(exactY1[1:pick_n,]),pp(trace3[1:pick_n,]),method = "shortsimplex")
wasserstein(pp(exactY1[1:pick_n,]),pp(trace3[1:pick_n,]),p = 1,tplan = tplan)

pick_n =1000
tplan = transport(pp(exactY1[1:pick_n,]),pp(trace4[1:pick_n,]),method = "shortsimplex")
wasserstein(pp(exactY1[1:pick_n,]),pp(trace4[1:pick_n,]),p = 1,tplan = tplan)



plot(trace2,xlim=c(0,1),ylim=c(-1,1))








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

