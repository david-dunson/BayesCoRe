setwd("C:/Users/leo/git/empiricalTensor/network_data_gpu/")

Z<- read.csv("./test/z.csv",header = F)$V1
D<- read.csv("./test/D.csv",header = F)$V1

A<- read.csv("./test/A.csv",header = F)$V1


i =1
d= matrix(D[((i*100)+1):((i+1)*100)],10)
range(d-t(d))
z = matrix(Z[((i*1E4)+1):((i+1)*1E4)],100)
range(z-t(z))


n=10
r=3
p=5

U = matrix(rnorm(n*r),r)

D_list=list()

for(i in 1:p){
D = matrix(rnorm(n^2),n)
D = D+t(D)
D_list[[i]]=D
}

UDU_list =list()
for(i in 1:p){
  UDU_list[[i]] = U%*% D_list[[i]]%*%t(U)
}

UDU_list


UD = U%*% matrix(unlist(D_list),n)

UD<-array(UD,dim = c(r,n,p))

for(i in 1:p){
  UD[,,i]<- t(UD[,,i])
}

UDU <- U%*% array(UD,dim=c(n,r*p))

UDU <- array(UDU,dim=c(r,r,p))


for(i in 1:p){
  print( max(abs(UDU[,,i] - UDU_list[[i]])))
}



z<- read.csv("./test/z.csv",header = F)$V1
p<- read.csv("./test/p.csv",header = F)$V1

set1 = (0.075 <= p) & (p <= 0.925)
hist(pnorm(z[set1])-(p[set1]))

set1 = (0.075 > p) | (p > 0.925)

  
hist(pnorm(z[set1])-p[set1])
