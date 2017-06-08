options(max.print=100)

a = -Inf
b = 5

exactDensity = function(theta){
	if(theta<b &theta >a){ theta^2/2/5^2
	}else {1E9}
}

x= seq(4,b+5,length.out=1000)



approximateDensity = function(theta, lambda = lam, pow =1){
	loss = theta^2/2/5^2

	k2 = theta - b
	
	if(k2>0)   loss = loss +  (k2/lambda)^pow
	
	loss
}

pdf_e = sapply(x, exactDensity)

lam = 1
pdf_a1 = sapply(x,function(y) approximateDensity(y, lam=lam,1))
pdf_a2 = sapply(x,function(y) approximateDensity(y, lam=lam,2))

require("ggplot2")

theta = x
n= length(theta)

df = data.frame( "theta"= x, "unnormalized.density" = exp(-c(pdf_e,pdf_a1, pdf_a2)), "v"= rep(c("intrinsic","extrinsic.1","extrinsic.2"),each=n))

df$v = factor(df$v, levels =c("intrinsic","extrinsic.1","extrinsic.2"))


library(latex2exp)


setwd("~/git/constrainedBayes/draft")

pdf("./density_truncated_normal.pdf",5,3)
ggplot(data=df, aes(x=theta, y= unnormalized.density, linetype=v,color=v))+ geom_line(size=.75)+ theme_bw()+ xlab(TeX('$\\theta$')) + ylab(TeX('$\\pi_{R}(\\theta) K(\\theta)$'))
dev.off()



pdf_a1 = sapply(x,function(y) approximateDensity(y, lam=1,1))
pdf_a2 = sapply(x,function(y) approximateDensity(y, lam=0.1,1))
pdf_a3 = sapply(x,function(y) approximateDensity(y, lam=0.01,1))

require("ggplot2")

theta = x
n= length(theta)

df = data.frame( "theta"= x, "unnormalized.density" = exp(-c(pdf_e,pdf_a1, pdf_a2,pdf_a3)), "lambda"= rep(c("intrinsic","1","0.1","0.01"),each=n))

df$lambda = factor(df$lambda, levels =c("intrinsic","1","0.1","0.01"))

ggplot(data=df, aes(x=theta, y= unnormalized.density, linetype=lambda,color=lambda))+ geom_line(size=.75)+ theme_bw()+ xlab(TeX('$\\theta$')) + ylab(TeX('$\\pi_{R}(\\theta) K(\\theta)$'))

