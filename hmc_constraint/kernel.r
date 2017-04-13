pcauchy(100) - pcauchy(-100)

hist(rcauchy(1000),breaks=1000, xlim=c(-1,1))


x= rcauchy(1000)

hist(x[x<1 & x>-1], breaks=100)