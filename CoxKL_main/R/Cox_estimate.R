pl_cal_theta <- function(lp, delta, time){
  delta = delta[order(time)]
  lp = lp[order(time)]
  S0 <- rev(cumsum(rev(exp(lp))))
  pl <- sum(delta*(lp-log(S0)))
  return(pl)
}

cal_surv_prob <- function(theta, delta, time){
  delta = delta[order(time)]
  theta = theta[order(time)]
  time = time[order(time)]
  
  theta_mat <- as.matrix(theta)
  delta_mat <- as.matrix(delta)
  S0=cal_S0(theta_mat, delta_mat)$S0
  
  Lambda0=cumsum(delta/S0)
  tmax = length(Lambda0)
  n = length(delta)
  S <- matrix(rep(0, (n*tmax)), n, tmax)
  for (i in 1:tmax){
    S[,i]=exp(-Lambda0[i]*exp(theta_mat))
  }
  
  return_list <- list("S"=S)
  return(return_list)
}

###############
Cox_Estimate <- function(z, delta, time, tol=1.0e-7){
  
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]

  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  
  p <- ncol(z_mat)
  beta = as.matrix(rep(0,p))
  repeat{
    diff=ddloglik(z_mat,delta_mat,beta)
    G=diff$L1
    H=diff$L2
    S0=diff$S0
    
    Lambda=cumsum(delta/S0)
    S=exp(-Lambda)
    temp=solve(H)%*%t(G)
    beta=beta+temp
    
    if(max(abs(temp))<tol) break
  }
  
  return(beta)
}
