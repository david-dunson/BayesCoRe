HMC_step<- function(q0,lambda,t0=0.5){
  p0 = rnorm(1)
  a = p0
  c = q0
  margin = sqrt(a^2 - 2*lambda*c)
  u = (-a-margin)/lambda
  # l =  (-a+margin)/lambda
  l =  max( (-a+margin)/lambda, 0 )
  t = l*t0 +u*(1-t0)
  q1= lambda * t^2 /2 + a*t +c
  p1 = lambda * t + a

  # print( (- q1*lambda + p1^2/2) -(  - q0*lambda + p0^2/2))
  q1
}

# HMC_step<- function(q0,lambda,t0=0.5){
#   p0 = rexp(1)
#   a = p0
#   b = q0
#   u = a/(-lambda)
#   # l =  (-a+margin)/lambda
#   l =  0
#   t = l*t0 +u*(1-t0)
#   q1= lambda * t^2 /2 + a*t +c
#   p1 = lambda * t + a
# 
#   # print( (- q1*lambda + p1^2/2) -(  - q0*lambda + p0^2/2))
#   q1
# }


q=list()
q[1]<- 1
for(i in 1:5000){
  q[[i+1]]<- HMC_step(q[[i]],-1,t0 = .5)
}

q =  unlist(q)
q = q[-c(1:100)]

# q= q[c(1:490)*10]
mean(q)
sd(q)

hist((q),breaks = 100)
acf(q)

ts.plot(q)

acf(q)
