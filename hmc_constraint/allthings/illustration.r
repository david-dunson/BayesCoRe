options(max.print=100)

a = 1
b = 2

exactDensity = function(theta){
	if(theta<b &theta >a){ theta^2/2
	}else {1E9}
}

x= seq(a-0.5,b+0.5,length.out=1000)





approximateDensity = function(theta, lambda = lam){
	loss = theta^2/2
	k1 = a - theta
	if(k1>0)   loss = loss + lambda * k1^2
	
	k2 = theta - b
	
	if(k2>0)   loss = loss + lambda * k2^2
	
	loss
}

pdf_e = sapply(x, exactDensity)

lam = 1E2
pdf_a1 = sapply(x, approximateDensity)
lam = 1E3
pdf_a2 = sapply(x, approximateDensity)
lam = 1E4
pdf_a3 = sapply(x, approximateDensity)
lam = 1E5
pdf_a4 = sapply(x, approximateDensity)



require("ggplot2")




theta = x
n= length(theta)

df = data.frame( "theta"= x, "unnormalized.density" = exp(-c(pdf_e, pdf_a1, pdf_a2, pdf_a3,pdf_a4)), "lambda"= rep(c("exact","100","1,000","10,000","100,000"),each=n))

setwd("~/git/empiricalTensor/draft")

pdf("./density_truncated_normal.pdf",5,3)
ggplot(data=df, aes(x=theta, y= unnormalized.density, color=lambda))+ geom_line(size=.75)+ theme_bw()
dev.off()


4/sqrt(2*100)