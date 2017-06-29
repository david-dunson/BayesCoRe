setwd("~/git/constrainedBayes/hmc_constraint/Network")
require("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


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


extractPosterior3D<-function(varname, d1,d2,d3, stan_fit){
  n = length( eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",1,",",1,",",1, "]`",sep=""))))
  L = matrix(0, n, d1*d2*d3)
  for(i in 1:d1){
    for(j in 1:d2){
      for(l in 1:d3){
          idx = (i-1)*d2*d3 + (j-1)*d3+ l;
          L[,idx] = eval(parse(text=paste(stan_fit,"@sim$samples[[1]]$`",varname,"[",i,",",j,",",l, "]`",sep="")))
      }
    }
  }
  L
}


#model
load("tensorA.RDa")
y=A
N = dim(A)[1]
p = dim(A)[3]

d1 = 5
d2 = 5
  
  
load(file="network_real_data_constrained.RDa")
trace1 = ss_fit

load(file="network_real_data_unconstrained.RDa")
trace2 = ss_fit

n_steps=1E4
sampling_idx<- c((n_steps/2+1):n_steps)

post_tau = extract(trace1,"tau")
ts.plot(post_tau[[1]])
# post_eta = extract(ss_fit,"eta")

post_U1 = extractPosteriorMat("U",N,d1,"trace1")
# acf(post_U[sampling_idx,1],lag.max = 40)
# acf(post_U[sampling_idx,N+1],lag.max = 40)
# acf(post_U[sampling_idx,3*N+1],lag.max = 40)
# ts.plot(post_U[sampling_idx,1])

post_U2 = extractPosteriorMat("U",N,d1,"trace2")
# acf(post_U[sampling_idx,1],lag.max = 40)
# acf(post_U[sampling_idx,N+1],lag.max = 40)
# ts.plot(post_U[sampling_idx,1])

df = data.frame( step =c(1:5000),
                 U=c(post_U1[sampling_idx,N+2], post_U2[sampling_idx,1]),
                 Constraint =as.factor(rep(c("Constrained","Unconstrained"), each=5000)))

pdf("../../draft/tucker_traceplot.pdf",8,3)
ggplot(data=df) + geom_path( aes(x=step, y=U),alpha =0.8,size=.8) +
  facet_wrap(~Constraint,scales = "free_y" )+ theme_bw()
dev.off()


ts.plot(post_U1[sampling_idx,1])



ts.plot(post_U1[,N-1])

acf_U1<- apply(post_U1[sampling_idx,], 2, function(x) acf(x,plot = F,lag.max = 40)$acf )
acf_U2<- apply(post_U2[sampling_idx,], 2, function(x)acf(x,plot = F,lag.max = 40)$acf)

df1 = data.frame("ACF"=c(acf_U1),"Lag"=as.factor(c(0:40)),Constraint="Constrained")
df2 = data.frame("ACF"=c(acf_U2),"Lag"=as.factor(c(0:40)),Constraint="Unconstrained")

df=rbind(df1,df2)

pdf("../../draft/tucker_acf.pdf",8,3)
ggplot(data = df, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
 scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)+
    facet_wrap(~Constraint)+ theme_bw()
dev.off()



ggplot(data = df, aes(x=Lag, y=ACF)) + geom_boxplot(outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.2,1))+ scale_x_discrete(breaks = c(0:8)*5)



# dev.off()

trace_eta = extract(trace1,"eta")
eta =  trace_eta[[1]][5000,,]

svd_eta = svd(eta)

plot(svd_eta$v[,1])
plot(svd_eta$v[,2])
plot(svd_eta$v[,3])
