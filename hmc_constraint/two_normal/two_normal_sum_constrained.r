
l=2

Sigma <- function(l){
  matrix(c( (l+2)/(l+4), -2/(l+4), 
           -2/(l+4),(l+2)/(l+4)),2)
}

Mu <- function(l){
  2/(l+4)
}

sqrtMat<-function(Sigma){
  ev= eigen(Sigma)
  Q= ev$vectors
  V = ev$values
  V[V<0]=0
  Q%*%diag(sqrt(V))%*%t(Q)
  }


Wass<-function(l){
  C1= Sigma(l)
  C2= Sigma(0)
  C1half=sqrtMat(C1)
  M1= Mu(l)
  M2= Mu(0)
  
  W2=(M1-M2)^2*2+ sum(diag(C1+C2- 2* sqrtMat(C1half%*%C2%*%C1half)))
  sqrt(W2)
}

l= seq(0,10,length.out = 100)

result<- sapply(10^(-l),Wass)

require("ggplot2")

df=data.frame("log10.lambda"=l,"W.2"=result)

pdf("../../draft/two_normal_wass.pdf",5,3)
ggplot(df,aes(x=log10.lambda,y=W.2))+ geom_line()+theme_bw()
dev.off()