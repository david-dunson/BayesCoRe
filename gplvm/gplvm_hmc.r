lambda = 1E3

#monotone penalty
computeXslope<- function(x){
  xw =outer(x,w,"*")
  C = cbind(cos(xw),sin(xw))
  A= as.matrix(Y- C%*%theta)
  #likelihood function
  U0 = sum(diag(A%*%diag(1/Sigma,nrow = length(Sigma),ncol=length(Sigma))%*%t(A))) * 0.5
  #penalty for non-monotone (penalize positive slope)
  slopeC= t(t(cbind( -sin(xw), cos(xw)))*w)
  slopeX = slopeC%*%theta
  slopeX
}


# U: -log(L*prior)
computeU<- function(x, theta, Sigma){
  x[1]=0
  xw =outer(x,w,"*")
  C = cbind(cos(xw),sin(xw))
  A= as.matrix(Y- C%*%theta)
  #likelihood function
  U0 = sum(diag(A%*%diag(1/Sigma,nrow = length(Sigma),ncol=length(Sigma))%*%t(A))) * 0.5
  #penalty for non-monotone (penalize positive slope)
  
  slopeC= t(t(cbind( -sin(xw), cos(xw)))*w)
  slopeX = slopeC%*%theta
  penalty = sum(lambda*(slopeX* (slopeX>0)))

  #total
  U0 + penalty
}



computeDeriX<- function(x){
  x[1]=0
  xw = outer(x,w,"*")  
  C = cbind( cos(xw), sin(xw))
  deriC0 = cbind( -sin(xw), cos(xw))
  deriC = t(t(deriC0)*w)
  SigmaMat = diag(2/Sigma,nrow = length(Sigma),ncol=length(Sigma))
  
  #deri in U0
  deriX_U0 = rowSums(((Y-C%*%theta)%*%(SigmaMat)%*% t(-theta)) * deriC) * 0.5
  #deri in penalty
  slopeC= deriC
  slopeX = slopeC%*%theta
  slopeCderiX=   t(t(-C)*w^2)
  slopeXderiX = slopeCderiX%*%theta
  deriX_penalty = lambda* rowSums(slopeXderiX* (slopeX>0))
  #total
  deriX = deriX_U0 + deriX_penalty
  deriX[1] = 0
  deriX
}

updateX <-function(x,eps = 1E-5, L =10,steps = 1000,tuning = T, ideal_AR =0.234, burnin =T){
  #update x
  microsteps= 100
  trace_x<- numeric()
  
  for(j in 1:(floor(steps/microsteps))){
    accept = 0
    for(i in 1:microsteps){
      m_x = rep(1,n)
      p_0 = rnorm(n) * sqrt(m_x)
      p_0[1] = 0
      q_x = x
      p_x = p_0
      p_x = p_x - eps/2 * computeDeriX(x)
      for(i in 1:L){
        q_x = q_x + eps* p_x/m_x
        p_x = p_x - eps * computeDeriX(x)
      }
      q_x = q_x + eps* p_x/m_x
      p_x = p_x - eps/2 * computeDeriX(x)
      
      cur_H = computeU(x, theta=theta,Sigma = Sigma )+ sum(p_0^2/m_x/2)
      new_H = computeU(q_x, theta=theta,Sigma = Sigma) + sum(p_x^2/m_x/2)
      
      if(log(runif(1))< (cur_H-new_H)){
        x = q_x
        if(tuning)  accept = accept+1
      }
      if((!tuning) & (!burnin)) trace_x<- rbind(trace_x,x)
    }
    #acceptance rate
    AR = accept / microsteps
    #reduce eps if AR too low, else increase eps
    if(tuning)  {
      eps = eps* exp(AR-ideal_AR)
      print(AR)
    }
    print(c("step: ",j, ", -log-L: ", computeU(x,theta,Sigma)))
  }
  return(list(x=x,trace_x=trace_x, eps=eps))
}