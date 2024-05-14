kl_cox <- function(RiskScore, eta_min, eta_max, length.out, eta_minby, df_internal, criteria, fold)
{
  ######
  df_internal <- as.data.frame(df_internal)
  n_internal <- nrow(df_internal)
  X_internal <- dplyr::select(df_internal, -c("time", "status"))
  t_internal <- df_internal$time
  delta_internal <- df_internal$status
  y_internal <- Surv(t_internal, delta_internal)
  p <- ncol(X_internal)
  group <- as.factor(1:p)
  
  if(fold=="LOOCV"){
    fold=n_internal
  }
  
  eta_min_wkg <- eta_min
  eta_max_wkg <- eta_max
  
  while(((eta_max_wkg-eta_min_wkg)/(length.out-1))>eta_minby){
    eta_vec <- seq(from = eta_min_wkg, to = eta_max_wkg, length.out = length.out)
    likelihood <- rep(NA, length.out)
    likelihood <- cv_tuning_CoxKL(eta_vec, fold, df_internal, RiskScore, criteria, likelihood)
    max_loc <- which(likelihood==max(likelihood))
    eta_min_wkg <- eta_vec[max(1,max_loc-1)][1]
    eta_max_wkg <- eta_vec[min(length.out, max_loc+1)][1]
  }
  eta_where_max <- eta_vec[max_loc][1]
  
  return_list <- list("eta"= eta_where_max, "likelihood"=likelihood)
  
  return(return_list)
}

cv_tuning_CoxKL <- function(eta_vec, fold, df_internal, RiskScore, criteria, likelihood){
  
  X_internal <- dplyr::select(df_internal, -c("time", "status"))
  t_internal <- df_internal$time
  delta_internal <- df_internal$status
  folds <- get_fold(nfolds = fold, delta_internal)
  p <- ncol(X_internal)
  group <- as.factor(1:p)
  k=0
  for (eta in eta_vec)
  {
    k=k+1
    likelihood_cv = rep(0, fold)
    LP_internal <- rep(0, nrow(df_internal))
    for (cv in 1:fold)
    {
      testIndexes <- which(folds==cv,arr.ind=TRUE)
      df_test <- df_internal[testIndexes, ]
      df_train <- df_internal[-testIndexes, ]
      
      X_train <- X_internal[-testIndexes, ]
      t_train <- df_train$time
      delta_train <- df_train$status
      y_train <- Surv(t_train, delta_train)
      
      X_test <- X_internal[testIndexes, ]
      t_test<- df_test$time
      delta_test <- df_test$status
      y_test <- Surv(t_test, delta_test)   
      
      beta_train <- KL_Cox_Estimate(X_train, delta_train, t_train, RiskScore[-testIndexes], eta=eta)
      
      #V&VH criteria
      if (criteria == "V&VH")
      {
        LP_train <- as.matrix(X_train)%*%as.matrix(beta_train)
        LP_internal <- as.matrix(X_internal)%*%as.matrix(beta_train)
        
        likelihood_cv[cv] <- pl_cal_theta(LP_internal, delta_internal, t_internal) - pl_cal_theta(LP_train, delta_train, t_train)        
      }
      
      #Linear Predictors & C-Index
      if (criteria == "LP: log-partial likelihood")
      {
        LP_internal[testIndexes] <- as.matrix(X_test)%*%as.matrix(beta_train)
      }
      
      if (criteria == "LP: C-Index")
      {
        LP_internal[testIndexes] <- as.matrix(X_test)%*%as.matrix(beta_train)
      }
      
      #C-Index
      if (criteria == "C-Index"){
        LP_test <- as.matrix(X_test)%*%as.matrix(beta_train)
        likelihood_cv[cv] <- Cindex(LP_test, Surv(t_test, delta_test))        
      }
    }
    
    #V&VH
    if (criteria == "V&VH"){
      likelihood[k] <- mean(likelihood_cv)
    }
    
    #Linear Predictors
    if (criteria == "LP: log-partial likelihood"){
      likelihood[k] <- pl_cal_theta(LP_internal, delta_internal, t_internal)
    }
    
    #Cross-validated C-Index
    if (criteria == "LP: C-Index"){
      likelihood[k] <- glmnet::Cindex(LP_internal, Surv(t_internal, delta_internal))
    }
    
    #Classical C-Index
    if (criteria == "C-Index"){
      likelihood[k] <- mean(likelihood_cv)
    }
  }
  
  return(likelihood)
}

###############
KL_Cox_Estimate <- function(z, delta, time, RS_internal, eta, tol=1.0e-7){
  
  delta = delta[order(time)]
  z = z[order(time),]
  RS_internal = RS_internal[order(time)]
  time = time[order(time)]
  
  z_mat <- as.matrix(z)
  delta_mat <- as.matrix(delta)
  
  p <- ncol(z_mat)
  beta = as.matrix(rep(0,p))
  RS_internal = as.matrix(RS_internal)
  repeat{
    diff=ddloglik_KL_RS(z_mat,delta_mat,beta, RS_internal, eta)
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

#### get folds
#for dividing data to training and testing
get_fold <- function(nfolds = 5, delta){
  n <- length(delta)
  ind1 <- which(delta==1)
  ind0 <- which(delta==0)
  n1 <- length(ind1)
  n0 <- length(ind0)
  fold1 <- 1:n1 %% nfolds
  fold0 <- (n1 + 1:n0) %% nfolds
  fold1[fold1==0] <- nfolds
  fold0[fold0==0] <- nfolds
  fold <- integer(n)
  fold[delta==1] <- sample(fold1)
  fold[delta==0] <- sample(fold0)
  return(fold)
} 


