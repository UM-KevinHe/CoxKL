library(mvtnorm)
library(survival)
library(Rcpp)
library(Matrix)
library(devtools)
devtools::install_github("dwangi/CoxKL")
library(CoxKL)
library(parallel)
###### this simulation code can be applied with parallel computing
##### detect max number of cores on this machine
max.cores = parallel::detectCores()
rep = 500
length.out = 5000

#eta_seq
KL2 <- function(x){
  Hcindex_KL <- rep(0, rep)
  beta_KL <- as.data.frame(matrix(rep(0, 6*rep), rep, 6))
  for (i in 1:rep){
    data_internal <- data_internal_list[[i]]$data
    Z_internal <- data_internal[,1:p]
    delta_internal <- data_internal$status
    time_internal <- data_internal$time
    data_test <- data_internal_list[[i]]$data_test
    Z_test <- data_test[,1:p]
    delta_test <- data_test$status
    time_test <- data_test$time
    rs_external2 <- as.matrix(Z_internal)%*%as.matrix(beta_external2) #calculate predicted risk score
    beta_KL[i,] <- KL_Cox_Estimate(Z_internal, delta_internal, time_internal, rs_external2, eta=x)
    lp_KL <- as.matrix(Z_test)%*%t(as.matrix(beta_KL[i,]))
    Hcindex_KL[i] <- glmnet::Cindex(lp_KL, Surv(time_test, delta_test))
  }
  
  mean_Hcindex_KL <- mean(replace(Hcindex_KL, Hcindex_KL == 0, NA), na.rm = TRUE)
  mean_beta_KL <- colMeans(replace(beta_KL, beta_KL == 0, NA), na.rm = TRUE)
  mean_Bias_KL <- mean(abs(mean_beta_KL-beta))
  ESE_KL <- mean(apply(replace(beta_KL, beta_KL == 0, NA), MARGIN = 2, function(x){sd(x, na.rm = TRUE)}))
  MSE_KL <- mean_Bias_KL^2+ESE_KL^2
  
  return_list <- list("Hcindex" = mean_Hcindex_KL, "Bias" = mean_Bias_KL, 
                      "SE" = ESE_KL, "MSE" = MSE_KL)
  return(return_list)
}

#########General setting for ture model
p_h<-2
p <- 2+2*p_h
#length of beta should be 2*p_h+4
beta <- c(0.3,-0.3,0.3,-0.3,0.3,-0.3)
###############external data
N_external=2000
upper_C_external <- 1.5
group1=rbinom(2000, 1, 0.5) 
set.seed(41)
#local/external
data_external <- sim.con_alt0(beta, p_h, N_external, upper_C=upper_C_external, group=group1)
delta_external <- data_external$status
time_external <- data_external$time

#external scenario 2
external_var_index2 <- c(1,3,5,6)
stacked_var_index2 <- c(2,4)
Z_external2 <- data_external[,external_var_index2]
#external2 beta
beta_external2 <- rep(0, p)
var_external2 <- names(Z_external2)
formula_external2 <- as.formula(paste('Surv(time, status)~',paste(var_external2, collapse= "+")))
cox_external2=coxph(formula_external2, data = data_external)
beta_external2[external_var_index2] <- cox_external2$coef

############### internal data/test data
N=100
upper_C <- 1.5 #Censoring rate ~ 40%
###############
Hcindex_internal <- rep(0, rep)
Hcindex_KL <- rep(0, (length.out+1))
Bias_KL <- rep(0, (length.out+1))
SE_KL <- rep(0, (length.out+1))
MSE_KL <- rep(0, (length.out+1))
beta_internal <- as.data.frame(matrix(rep(0, 6*rep), rep, 6))

old_time_GD=proc.time()

data_internal_list <- vector(mode="list", length=rep)
eta_list <- vector(mode="list", length=(length.out+1))

for (i in 0:length.out){
  eta_list[[i+1]] <- 0.001*i
}

for (i in 1:rep){
  set.seed(i)
  data_internal_list[[i]]$data <- sim.con_alt0(beta, p_h, N, upper_C=upper_C, group=rep(1,N))
  data_internal_list[[i]]$data_test <- sim.con_alt0(beta, p_h, 1000, upper_C=upper_C, group=rep(1,1000))
}

models_KL2 <- mclapply(eta_list, KL2, mc.cores = floor(0.5*max.cores))

for (i in 1:rep){
  set.seed(i)
  ############internal data########################
  data_internal <-  data_internal_list[[i]]$data
  delta_internal <- data_internal$status
  Z_internal <- data_internal[,1:p]
  time_internal <- data_internal$time
  #censoring rate
  1-mean(delta_internal)
  ############test data########################
  data_test <- data_internal_list[[i]]$data_test
  delta_test <- data_test$status
  Z_test <- data_test[,1:p]
  time_test <- data_test$time
  #############estimation################################
  var_internal <- names(Z_internal)
  formula_internal <- as.formula(paste('Surv(time, status)~',paste(var_internal, collapse= "+")))
  
  #Internal model
  cox_int <- coxph(formula_internal, data = data_internal)
  beta_internal[i,] <- cox_int$coef
  lp_internal <- predict(cox_int, data_test, type="lp")
  Hcindex_internal[i] <- glmnet::Cindex(lp_internal, Surv(time_test, delta_test))
}

mean_Hcindex_internal <- mean(replace(Hcindex_internal, Hcindex_internal == 0, NA), na.rm = TRUE)
mean_beta_internal <- colMeans(replace(beta_internal, beta_internal == 0, NA), na.rm = TRUE)
mean_Bias_internal <- mean(abs(mean_beta_internal-beta))
ESE_internal <- mean(apply(replace(beta_internal, beta_internal == 0, NA), MARGIN = 2, function(x){sd(x, na.rm = TRUE)}))
MSE_internal <- mean_Bias_internal^2+ESE_internal^2

for (i in 1:(length.out+1)){
  {
    results_KL2 <- models_KL2[[i]]
    Bias_KL[i] <- results_KL2$Bias
    SE_KL[i] <- results_KL2$SE
    MSE_KL[i] <- results_KL2$MSE
    Hcindex_KL[i] <- results_KL2$Hcindex
  }
}

time_GD=proc.time()-old_time_GD
time_GD

mean_Hcindex2 <- mean_Hcindex_internal
mean_Bias2 <- mean_Bias_internal
mean_SE2 <- ESE_internal
mean_MSE2 <- MSE_internal
mean_Hcindex_KL2 <- Hcindex_KL
mean_Bias_KL2 <- Bias_KL
mean_SE_KL2 <- SE_KL
mean_MSE_KL2 <- MSE_KL

length.out = 5000
eta_1 <- seq(0,5,length.out=(length.out+1))

#Cindex
library(colorspace)
cols <- c("#E0B0A9", "#5E81B5")
library(ggplot2)

mr <- data.frame(Method = c(
  rep("KL", (length.out+1)),
  rep("Internal", (length.out+1))),
  value = c(
    mean_Hcindex_KL2,
    rep(mean_Hcindex2, (length.out+1))),
  eta = c(eta_1, eta_1)
)
mr$Method <- factor(mr$Method,
                    levels = c('Internal', 'KL'),ordered = TRUE)
ggplot(mr, aes(x=eta, y=value, group=Method)) + 
  geom_line(aes(linetype=Method, color=Method, size = Method)) + 
  geom_point(x=0, y=mean_Hcindex2, color="#E0B0A9", size=5)+
  labs(y = "C-Index", x = bquote(eta))+  scale_color_manual(values=cols)+ theme_bw() + 
  scale_linetype_manual(values=c("dotted", "solid"))+scale_size_manual( values = c(1,2) ) +
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 14))+theme(legend.position = "none")

#MSE
library(colorspace)
cols <- c("#E0B0A9", "#5E81B5")
library(ggplot2)

mr <- data.frame(Method = c(
  rep("KL", (length.out+1)),
  rep("Internal", (length.out+1))),
  value = c(
    mean_MSE_KL2,
    rep(mean_MSE2, (length.out+1))),
  eta = c(eta_1, eta_1)
)
mr$Method <- factor(mr$Method,
                    levels = c('Internal', 'KL'),ordered = TRUE)
ggplot(mr, aes(x=eta, y=value, group=Method)) + 
  geom_line(aes(linetype=Method, color=Method, size = Method)) + 
  geom_point(x=0, y=mean_MSE2, color="#E0B0A9", size=5)+
  labs(y = "MSE", x = bquote(eta))+  scale_color_manual(values=cols)+ theme_bw() + scale_linetype_manual(values=c("dotted", "solid"))+scale_size_manual( values = c(1,2) )+
  theme(panel.border = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 14))+theme(legend.position = "none")












