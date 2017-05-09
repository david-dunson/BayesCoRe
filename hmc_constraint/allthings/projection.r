setwd("~/git/empiricalTensor/hmc_constraint/")


require("rstiefel")

N= 10
d = 5
m = 100

F<-matrix(1,N,d)

X_list = lapply(c(1:m), function(x) rmf.matrix(F))
# X = do.call("rbind", lapply(c(1:m), function(x) rmf.matrix(F)))

Xsum= matrix(0, N, d)

for(i in 1:m){
	Xsum  = Xsum + X_list[[i]]
}



X = X_list[[1]]

svdX = svd(X)


t(X)%*%X


X = X+ rnorm( N*d) *0.01
A.qr=qr(X)

Y= qr.Q(A.qr)
qr.R(A.qr)


t(X)%*%X
t(Y)%*%Y


X

Y = Y%*%diag(round(diag(t(X)%*%Y)))

X
Y