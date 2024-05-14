### Simulate continuous-time survival data
AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}

simu_z <- function(n, size.groups)
{
  Sigma_z1=diag(size.groups) # =p
  Corr1<-AR1(0.5,size.groups) #correlation structure 0.5 0.6
  diag(Corr1) <- 1
  Sigma_z1<- Corr1
  pre_z= rmvnorm(n, mean=rep(0,size.groups), sigma=Sigma_z1)
  return(pre_z)
}

### with interaction and group information (paper version: estimation)
sim.con_alt0 <- function(beta, p_h, N, upper_C, group){
  Z1 <- as.matrix(simu_z(N, p_h))
  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h)
  Z3 <- group*2 + rnorm(N)
  Z4 <- rnorm(N) - group*2
  Z <- cbind(Z1, Z2, Z3, Z4)
  z <- cbind(Z1, Z2, Z3, Z4)
  Z.char <- paste0('Z', 1:ncol(z))
  U=runif(N, 0, 1)
  #Weibull
  lambda=1
  nu=2
  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
  censoring=runif(N,0,upper_C) 
  tcens=(censoring<time) # censoring indicator
  delta=1-tcens
  time=time*(delta!=0)+censoring*(delta==0)
  
  ###order data; 
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  data <- as.data.frame(cbind(z, delta, time))
  colnames(data) <- c(Z.char, "status", "time")
  
  return(data)
}





