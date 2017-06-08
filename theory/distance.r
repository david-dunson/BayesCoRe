
K<- function(p,lambda){
  prod( exp(- abs(sum(p)-1)/lambda))
}


a= 1E-3

n<- 1E5
P = cbind(rnorm(n),rnorm(n))

mean(apply(P,1, function(x){K(x,a)}))

1/2/exp(1/4)/sqrt(pi)





2*a*(a*(exp(-1/a)-1)+1)

# a - a^3* exp(-2/a)* (exp(1/a)-1)^2

eigen(matrix(c(1,1,1,2),2))
