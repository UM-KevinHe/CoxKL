pl_cal_theta <- function(lp, delta, time){
  delta = delta[order(time)]
  lp = lp[order(time)]
  S0 <- rev(cumsum(rev(exp(lp))))
  pl <- sum(delta*(lp-log(S0)))
  return(pl)
}

cal_surv_prob <- function(z, delta, time, beta){
  delta = delta[order(time)]
  z = z[order(time),]
  time = time[order(time)]
  
  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  beta <- as.matrix(beta)
  diff=ddloglik(z_mat,delta_mat,beta)
  S0=diff$S0
  
  Lambda0=cumsum(delta/S0)
  tmax = length(Lambda0)
  n = nrow(z_mat)
  S <- matrix(rep(0, (n*tmax)), n, tmax)
  for (i in 1:tmax){
    S[,i]=exp(-Lambda0[i]*exp(z_mat%*%beta))
  }
  
  return_list <- list("S"=S)
  return(return_list)
}

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

### most-basic simulation for continuous data
sim.con <- function(beta, N, Z.char, upper_C){
  p_h <- 0.5*length(beta)
  Z1 <- as.matrix(simu_z(N, p_h))
  Z2 <- matrix(rbinom(N*p_1,1,0.5),N,p_h)
  Z <- cbind(Z1, Z2)
  U=runif(N, 0,1)
  #Exponential 
  lambda = 0.5
  time=-log(U)/(lambda*exp(Z%*%beta)) 
  #Weibull
  #lambda=0.5
  #nu=10
  #time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
  #censoring=runif(N,1,2) #0.9
  censoring=runif(N,0,upper_C) #0 or 0.5
  tcens=(censoring<time) # censoring indicator
  delta=1-tcens
  time=time*(delta!=0)+censoring*(delta==0)
  
  ###order data; 
  delta = delta[order(time)]
  Z = Z[order(time),]
  time = time[order(time)]
  data <- as.data.frame(cbind(Z, delta, time))
  colnames(data) <- c(Z.char, "status", "time")
  
  return(data)
}

### with interaction and group information (paper version: estimation)
sim.con_alt0 <- function(beta, p_h, N, upper_C, group){
  #p_h <- 0.5*(length(beta)-2)
  Z1 <- as.matrix(simu_z(N, p_h))
  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h)
  Z3 <- group*2 + rnorm(N)
  Z4 <- rnorm(N) - group*2
  #Z5 <- Z3 * Z1[,1]
  #Z6 <- Z4 * Z2[,1]
  #beta <- c(beta,-0.3, 0.3)
  #Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6)
  Z <- cbind(Z1, Z2, Z3, Z4)
  z <- cbind(Z1, Z2, Z3, Z4)
  Z.char <- paste0('Z', 1:ncol(z))
  U=runif(N, 0,1)
  #Exponential 
  #lambda = 5
  #time=-log(U)/(lambda*exp(Z%*%beta)) 
  #Weibull
  lambda=1
  nu=2
  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
  #censoring=runif(N,1,2) #0.9
  censoring=runif(N,0,upper_C) #0 or 0.5
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

## with interaction and group information (paper version: prediction)
sim.con_alt1 <- function(beta, p_h, N, upper_C, group){
  #p_h <- 0.5*(length(beta)-2)
  Z1 <- as.matrix(simu_z(N, p_h))
  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h)
  #Z3 <- group
  Z3 <- group*2 + rnorm(N)
  Z4 <- Z3 * Z1[,1]
  Z5 <- rbinom(N, 1, 0.5)
  Z6 <- rnorm(N, 0, 1)
  beta <- c(beta,1, 1)
  Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6)
  z <- cbind(Z1, Z2, Z3)
  Z.char <- paste0('Z', 1:ncol(z))
  U=runif(N, 0,1)
  #Exponential 
  #lambda = 5
  #time=-log(U)/(lambda*exp(Z%*%beta)) 
  #Weibull
  lambda=1
  nu=2
  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
  #censoring=runif(N,1,2) #0.9
  censoring=runif(N,0,upper_C) #0 or 0.5
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

## with interaction and group information (paper version: prediction - tree)
sim.con_alt <- function(beta, p_h, N, upper_C, group){
  #p_h <- 0.5*(length(beta)-2)
  Z1 <- as.matrix(simu_z(N, p_h))
  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h) #0.3, -0.3, 0.3, -0.3, -0.3, 0.3, 0.3, -0.3
  #Z3 <- group
  Z3 <- group*2 + rnorm(N)
  Z4 <- Z3 * Z1[,1]
  Z5 <- Z3 * Z2[,1]
  Z6 <- Z1[,2] * Z2[,2]
  Z7 <- rbinom(N, 1, 0.5)
  Z8 <- rnorm(N, 0, 1)
  beta <- c(beta, 0.5, -0.5, 1, 1)
  Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8)
  z <- cbind(Z1, Z2, Z3)
  Z.char <- paste0('Z', 1:ncol(z))
  U=runif(N, 0,1)
  #Exponential 
  #lambda = 5
  #time=-log(U)/(lambda*exp(Z%*%beta)) 
  #Weibull
  lambda=1
  nu=2
  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
  #censoring=runif(N,1,2) #0.9
  censoring=runif(N,0,upper_C) #0 or 0.5
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

### with interaction and group information (paper version: Nov1)
#sim.con_alt1 <- function(beta, p_h, N, upper_C, group){
#  #p_h <- 0.5*(length(beta)-2)
#  Z1 <- as.matrix(simu_z(N, p_h))
#  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h)
#  Z3 <- group*2 + rnorm(N)
#  Z4 <- rnorm(N) - group*2
#  Z5 <- Z3 * Z1[,1]
#  Z6 <- Z4 * Z2[,1]
#  #Z5 <- rbinom(N, 1, 0.5)
#  #Z6 <- rnorm(N, 0, 1)
#  beta <- c(beta,-0.3, 0.3)
#  Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6)
#  z <- cbind(Z1, Z2, Z3, Z4)
#  Z.char <- paste0('Z', 1:ncol(z))
#  U=runif(N, 0,1)
#  #Exponential 
#  #lambda = 5
#  #time=-log(U)/(lambda*exp(Z%*%beta)) 
#  #Weibull
#  lambda=1
#  nu=2
#  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
#  #censoring=runif(N,1,2) #0.9
#  censoring=runif(N,0,upper_C) #0 or 0.5
#  tcens=(censoring<time) # censoring indicator
#  delta=1-tcens
#  time=time*(delta!=0)+censoring*(delta==0)
#  
#  ###order data; 
#  delta = delta[order(time)]
#  z = z[order(time),]
#  time = time[order(time)]
#  data <- as.data.frame(cbind(z, delta, time))
#  colnames(data) <- c(Z.char, "status", "time")
#  
#  return(data)
#}

### with interaction and group information (basic version)
#sim.con_alt <- function(beta, p_h, N, upper_C, group){
#  #p_h <- 0.5*length(beta)
#  Z1 <- as.matrix(simu_z(N, p_h))
#  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h)
#  Z3 <- group*rbinom(N,1,0.5) #for different categrical covariate distribution of different group
#  Z4 <- group + rnorm(N) #for different continuous covariate distribution of different group
#  Z5 <- Z1[,p_h]*Z2[,p_h] #for interaction between last continuous variable and last categorical variable
#  Z6 <- Z3*Z4 # for interaction term between group-dependent continuous and categorical variable
#  Z7 <- rnorm(N)
#  Z8 <- rbinom(N,1,0.5)
#  #Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6)
#  Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8)
#  beta <- c(beta, c(1,-1))
#  Z.char <- paste0('Z', 1:ncol(Z))
#  U=runif(N, 0,1)
#  #Exponential 
#  #lambda = 0.5
#  #time=-log(U)/(lambda*exp(Z%*%beta)) 
#  #Weibull
#  lambda=1
#  nu=2
#  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
#  censoring=runif(N,0,upper_C) #0 or 0.5
#  tcens=(censoring<time) # censoring indicator
#  delta=1-tcens
#  time=time*(delta!=0)+censoring*(delta==0)
#  
#  ###order data; 
#  delta = delta[order(time)]
#  Z = Z[order(time),]
#  time = time[order(time)]
#  
#  data_all <- as.data.frame(cbind(Z, delta, time))
#  colnames(data_all) <- c(Z.char, "status", "time")
#  #data_sub <- as.data.frame(cbind(Z[,1:(ncol(Z)-2)], delta, time))
#  #colnames(data_sub) <- c(Z.char[1:(ncol(Z)-2)], "status", "time")
#  data_sub <- as.data.frame(cbind(Z[,1:(ncol(Z)-4)], delta, time))
#  colnames(data_sub) <- c(Z.char[1:(ncol(Z)-4)], "status", "time")
#  
#  #return(list("data_all"= data_all, "data_sub"=data_sub))
#  return(data_sub)
#}

### with more interactions and group information
#sim.con_alt0 <- function(beta, p_h, N, upper_C, group){
#  #p_h <- 0.5*(length(beta)-2)
#  Z1 <- as.matrix(simu_z(N, p_h)) #0.3, -0.3
#  Z2 <- matrix(rbinom(N*p_h,1,0.5),N,p_h) #0.3, -0.3
#  #Z3 <- group
#  Z3 <- group*2 + rnorm(N) #-0.3
#  Z4 <- Z3 * Z1[,1] #0.3
#  Z5 <- Z3 * Z2[,1] #0.3
#  Z6 <- Z1[,2] * Z2[,2] #-0.3
#  Z7 <- rbinom(N, 1, 0.5) #1
#  Z8 <- rnorm(N, 0, 1) #1
#  #beta <- c(beta,1, 1)
#  beta <- c(beta, 0.5, -0.5, 1, 1)
#  Z <- cbind(Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8)
#  z <- cbind(Z1, Z2, Z3)
#  Z.char <- paste0('Z', 1:ncol(z))
#  U=runif(N, 0,1)
#  #Exponential 
#  #lambda = 5
#  #time=-log(U)/(lambda*exp(Z%*%beta)) 
#  #Weibull
#  lambda=1
#  nu=2
#  time=(-log(U)/(lambda*exp(Z%*%beta)))^(1/nu) 
#  #censoring=runif(N,1,2) #0.9
#  censoring=runif(N,0,upper_C) #0 or 0.5
#  tcens=(censoring<time) # censoring indicator
#  delta=1-tcens
#  time=time*(delta!=0)+censoring*(delta==0)
#  
#  ###order data; 
#  delta = delta[order(time)]
#  z = z[order(time),]
#  time = time[order(time)]
#  data <- as.data.frame(cbind(z, delta, time))
#  colnames(data) <- c(Z.char, "status", "time")
#  
#  return(data)
#}


