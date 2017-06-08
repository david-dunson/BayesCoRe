n=10
p=5
Sigma =c(1:5)

A<- matrix(rnorm(n*p),n)

f=function(A,Sigma){
  sum(diag(A%*%diag(1/Sigma)%*%t(A)))
}

###derivative of A
deri_f<- function(a,b){
  (f(b,Sigma)-f(a,Sigma))/ sum(b-a)
}

deri_A =A

for(i in 1:length(deri_A)){
  A1=A
  A1[i]<- A1[i]+1E-8
  deri_A[i] = deri_f(A,A1)
}


A%*%(diag(2/Sigma))-deri_A

###derivative of sigma
deri_f<- function(a,b){
  (f(A,b)-f(A,a))/ sum(b-a)
}

deri_Sigma =Sigma

for(i in 1:length(Sigma)){
  Sigma1=Sigma
  Sigma1[i]<- Sigma1[i]+1E-8
  deri_Sigma[i] = deri_f(Sigma,Sigma1)
}

deri_Sigma

diag(-diag(1/Sigma)%*%t(A)%*%A%*%diag(1/Sigma))


########derivative of X and beta
r=2
n=10
p=5

y<- matrix(rnorm(n*p),n)
X<- matrix(rnorm(n*r),n)
beta<- matrix(rnorm(p*r),r)

loss<- function(X,beta,Sigma){
  A= y- X%*%beta
  sum(diag(A%*%diag(1/Sigma)%*%t(A)))
}

deri_X = X
for(i in 1:length(X)){
  X1 = X
  X1[i] = X[i]+1E-8
  deri_X[i] = (loss(X1,beta,Sigma) - loss(X,beta,Sigma))/1E-8
}

deri_X - A%*%(diag(2/Sigma))%*% t(-beta)

deri_beta = beta
for(i in 1:length(beta)){
  beta1 = beta
  beta1[i] = beta1[i]+1E-8
  deri_beta[i] = (loss(X,beta1,Sigma) - loss(X,beta,Sigma))/1E-8
}

deri_beta 
-t(X)%*%A%*%(diag(2/Sigma))

###################

r = 20
w = seq(0,pi,length.out = r/2+1)
w = w[-1]
x<- runif(n)
beta<- matrix(rnorm(r*p),r)

loss_x<- function(x){
  xw =outer(x,w,"*")
  X = cbind(cos(xw),sin(xw))
  loss(X,beta,Sigma)
}

deri_x = x
for(i in 1:length(x)){
  x1 = x
  x1[i] = x[i]+1E-8
  deri_x[i] = (loss_x(x1)-loss_x(x))/1E-8
}

xw =outer(x,w,"*")
X = cbind(cos(xw),sin(xw))
Xderi =  cbind(-sin(xw),cos(xw))
Xderi = t(t(Xderi)*w)

rowSums(((y-X%*%beta)%*%(diag(2/Sigma))%*% t(-beta)) * Xderi)
deri_x

