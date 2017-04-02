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
  plot(NA,xlim=c(0,nrow(img)),ylim=c(0,ncol(img)),xlab="",ylab="")
  rasterImage(img,0,0,nrow(img),ncol(img))
}


UDU = U%*%D[,,2]%*%t(U)

setwd("c:/Users/leo/Desktop/ARI_Slides_032017/figures/")

png("large_network1.png",800,800)
plotHugeMat( U%*%D[,,3]%*%t(U))
dev.off()


png("large_network2.png",800,800)
plotHugeMat( U%*%D[,,4]%*%t(U))
dev.off()

png("large_network3.png",800,800)
plotHugeMat( U%*%D[,,5]%*%t(U))
dev.off()

png("large_network_core.png",800,800)
plotHugeMat(D[,,1])
dev.off()

D[1,1,1]=D[2,2,1]

png("large_network_core.png",800,800)
plotHugeMat(U[,1:10])
dev.off()



plotHugeMat(pnorm(UDU))

UDU = U%*%D[,,2]%*%t(U)






trace_D<- read.csv("./test/trace_D.csv",header = F)$V1
trace_D = t(matrix(trace_D,ncol=1000))
acf(trace_D[400:1000,100],lag.max = 100)

hist(pnorm(UDU),breaks = 1000)

