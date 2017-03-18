setwd("C:/Users/leo/git/empiricalTensor/network_data_gpu/")

U<- read.csv("./test/U.csv",header = F)$V1
D<- read.csv("./test/D.csv",header = F)$V1
avgA<- read.csv("./test/avgA.csv",header = F)$V1

n = 1105
r = 30
p = 44

U<- array(U,dim = c(n,r))
D<- array(D,dim = c(r,r,p))
avgA<- array(avgA, dim=c(n,n))

plotHugeMat<-function(img){
  img=(img -min(img))/(max(img)-min(img))
  plot(NA,xlim=c(0,nrow(img)),ylim=c(0,ncol(img)))
  rasterImage(img,0,0,nrow(img),ncol(img))
}


UDU = U%*%D[,,2]%*%t(U)

plotHugeMat(avgA)
plotHugeMat(pnorm(UDU))

UDU = U%*%D[,,2]%*%t(U)


plotHugeMat(D[,,1])

D[,,1]



trace_D<- read.csv("./test/trace_D.csv",header = F)$V1
trace_D = t(matrix(trace_D,ncol=1000))
acf(trace_D[400:1000,100],lag.max = 100)

hist(pnorm(UDU),breaks = 1000)

