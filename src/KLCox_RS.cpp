#include <Rcpp.h>
#include <RcppEigen.h>
#include <math.h>
//[[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
//using Eigen::MatrixXd; 
  
// [[Rcpp::export]]
List rev_cumsum(Eigen::MatrixXd X){
  int n = X.rows(); // X should be n*1 matrix (vector)
  double sum_to_i = 0;
  Eigen::MatrixXd rev_X(n,1);
  rev_X.setZero(n,1); 
  Eigen::MatrixXd cumsum_rev_X(n,1);
  cumsum_rev_X.setZero(n,1); 
  Eigen::MatrixXd cumsum_X(n,1);
  cumsum_X.setZero(n,1);  

  rev_X = X.reverse();
  for (int i = 0; i < n; i++){
    sum_to_i += rev_X(i,0);
    cumsum_rev_X(i,0) = sum_to_i;
  }

  cumsum_X = cumsum_rev_X.reverse();

  List result;
  result["cumsum"]=cumsum_X;
   
  return result;
}

// [[Rcpp::export]]
List ddloglik_2KL_RS(Eigen::MatrixXd Z, Eigen::MatrixXd delta, Eigen::MatrixXd beta, Eigen::MatrixXd theta_tilde1, Eigen::MatrixXd theta_tilde2, double eta1, double eta2){
  int p = beta.rows();
  int n = delta.rows();
  double loglik=0;
  Eigen::MatrixXd L1(p,1);
  L1.setZero(p,1);
  Eigen::MatrixXd L2(p,p);
  L2.setZero(p,p);
  Eigen::MatrixXd temp_2(n,1);
  temp_2.setZero(n,1); 

  //Calcualte S0 and S1 for beta
  Eigen::MatrixXd theta(n,1);
  theta.setZero(n,1);
  Eigen::MatrixXd exp_theta(n,1);
  exp_theta.setZero(n,1);
  Eigen::MatrixXd S0(n,1);
  S0.setZero(n,1);
  Eigen::MatrixXd S1_pre(n,p);
  S1_pre.setZero(n,p);
  Eigen::MatrixXd S1(n,p);
  S1.setZero(n,p);
  Eigen::MatrixXd S1_coli(n,1);
  S1_coli.setZero(n,1);
  Eigen::MatrixXd temp_1(n,p);
  temp_1.setZero(n,p); 
  Eigen::MatrixXd S2_i1(n,p);
  S2_i1.setZero(n,p);
  Eigen::MatrixXd S2_i(n,p);
  S2_i.setZero(n,p); 
  Eigen::MatrixXd S2_i_colj(n,p);
  S2_i_colj.setZero(n,p);   
  Eigen::MatrixXd V_i(n,p);
  V_i.setZero(n,p);
  Eigen::MatrixXd V1_i(n,p);
  V1_i.setZero(n,p);

  //Calcualte S0 and S1 for beta_tilde1
  //Eigen::MatrixXd theta_tilde1(n,1);
  //theta_tilde1.setZero(n,1);
  Eigen::MatrixXd exp_theta_tilde1(n,1);
  exp_theta_tilde1.setZero(n,1);
  Eigen::MatrixXd S0_tilde1(n,1);
  S0_tilde1.setZero(n,1);
  Eigen::MatrixXd S1_pre_tilde1(n,p);
  S1_pre_tilde1.setZero(n,p);
  Eigen::MatrixXd S1_tilde1(n,p);
  S1_tilde1.setZero(n,p);
  Eigen::MatrixXd S1_coli_tilde1(n,1);
  S1_coli_tilde1.setZero(n,1);

  //Calcualte S0 and S1 for beta_tilde2
  //Eigen::MatrixXd theta_tilde2(n,1);
  //theta_tilde2.setZero(n,1);
  Eigen::MatrixXd exp_theta_tilde2(n,1);
  exp_theta_tilde2.setZero(n,1);
  Eigen::MatrixXd S0_tilde2(n,1);
  S0_tilde2.setZero(n,1);
  Eigen::MatrixXd S1_pre_tilde2(n,p);
  S1_pre_tilde2.setZero(n,p);
  Eigen::MatrixXd S1_tilde2(n,p);
  S1_tilde2.setZero(n,p);
  Eigen::MatrixXd S1_coli_tilde2(n,1);
  S1_coli_tilde2.setZero(n,1);  

  //Calculate S0 with beta
  theta = Z*beta;
  exp_theta = (theta.array().exp()).matrix();
  List S0_rev_cumsum = rev_cumsum(exp_theta);
  S0 = S0_rev_cumsum["cumsum"];

  //Calcualte S0_tilde1 with beta_tilde1
  //theta_tilde1= Z*beta_tilde1;
  exp_theta_tilde1 = (theta_tilde1.array().exp()).matrix();
  List S0_rev_cumsum_tilde1 = rev_cumsum(exp_theta_tilde1);
  S0_tilde1 = S0_rev_cumsum_tilde1["cumsum"]; 

  //Calcualte S0_tilde2 with beta_tilde2
  //theta_tilde2 = Z*beta_tilde2;
  exp_theta_tilde2 = (theta_tilde2.array().exp()).matrix();
  List S0_rev_cumsum_tilde2 = rev_cumsum(exp_theta_tilde2);
  S0_tilde2 = S0_rev_cumsum_tilde2["cumsum"];  
  
  for (int i = 0; i < p; i++){
   //Calculate S1 with beta
   S1_pre.col(i) = Z.col(i).cwiseProduct(exp_theta);
   List S1_col_rev_cumsum = rev_cumsum(S1_pre.col(i));
   S1_coli = S1_col_rev_cumsum["cumsum"];
   S1.col(i) = S1_coli;

   //Calcualte S1_tilde1 with beta_tilde1
   S1_pre_tilde1.col(i) = Z.col(i).cwiseProduct(exp_theta_tilde1);
   List S1_col_rev_cumsum_tilde1 = rev_cumsum(S1_pre_tilde1.col(i));
   S1_coli_tilde1 = S1_col_rev_cumsum_tilde1["cumsum"];
   S1_tilde1.col(i) = S1_coli_tilde1;

   //Calcualte S1_tilde2 with beta_tilde2
   S1_pre_tilde2.col(i) = Z.col(i).cwiseProduct(exp_theta_tilde2);
   List S1_col_rev_cumsum_tilde2 = rev_cumsum(S1_pre_tilde2.col(i));
   S1_coli_tilde2 = S1_col_rev_cumsum_tilde2["cumsum"];
   S1_tilde2.col(i) = S1_coli_tilde2;   

   temp_1.col(i) = delta.cwiseProduct(Z.col(i)+eta1*(S1_tilde1.col(i).cwiseQuotient(S0_tilde1))+eta2*(S1_tilde2.col(i).cwiseQuotient(S0_tilde2))-(1+eta1+eta2)*(S1.col(i).cwiseQuotient(S0)));
  } 

  L1 = temp_1.colwise().sum();
  temp_2 = delta.cwiseProduct(theta-S0.array().log().matrix());
  loglik = temp_2.sum();
  
  for (int i=0; i < p; i++){
    for (int j = 0; j < p; j++){
      S2_i1.col(j) = (Z.col(j).cwiseProduct(exp_theta)).cwiseProduct(Z.col(i));
      List S2i_col_rev_cumsum = rev_cumsum(S2_i1.col(j));
      //S1_coli = S1_col_rev_cumsum["cumsum"];
     // S1.col(i) = S1_coli;
      S2_i_colj = S2i_col_rev_cumsum["cumsum"];
      S2_i.col(j) = S2_i_colj;
      V1_i.col(j) = ((S2_i.col(j).cwiseQuotient(S0)))-((S1.col(j).cwiseProduct(S1.col(i))).cwiseQuotient(S0.array().square().matrix()));
      V_i.col(j) = V1_i.col(j).cwiseProduct(delta);
    } 
    L2.col(i) = (V_i.colwise().sum())*(1+eta1+eta2);
  }

  List result;
  result["loglik"]=loglik;
  result["L1"]=L1;
  result["L2"]=L2;
  result["S0"]=S0;
  
  return result;
}

// [[Rcpp::export]]
List ddloglik_KL_RS(Eigen::MatrixXd Z, Eigen::MatrixXd delta, Eigen::MatrixXd beta, Eigen::MatrixXd theta_tilde, double eta){
  int p = beta.rows();
  int n = delta.rows();
  double loglik=0;
  Eigen::MatrixXd L1(p,1);
  L1.setZero(p,1);
  Eigen::MatrixXd L2(p,p);
  L2.setZero(p,p);
  Eigen::MatrixXd temp_2(n,1);
  temp_2.setZero(n,1); 
  //Calcualte S0 and S1 for beta
  Eigen::MatrixXd theta(n,1);
  theta.setZero(n,1);
  Eigen::MatrixXd exp_theta(n,1);
  exp_theta.setZero(n,1);
  Eigen::MatrixXd S0(n,1);
  S0.setZero(n,1);
  Eigen::MatrixXd S1_pre(n,p);
  S1_pre.setZero(n,p);
  Eigen::MatrixXd S1(n,p);
  S1.setZero(n,p);
  Eigen::MatrixXd S1_coli(n,1);
  S1_coli.setZero(n,1);
  Eigen::MatrixXd temp_1(n,p);
  temp_1.setZero(n,p); 
  Eigen::MatrixXd S2_i1(n,p);
  S2_i1.setZero(n,p);
  Eigen::MatrixXd S2_i(n,p);
  S2_i.setZero(n,p); 
  Eigen::MatrixXd S2_i_colj(n,p);
  S2_i_colj.setZero(n,p);   
  Eigen::MatrixXd V_i(n,p);
  V_i.setZero(n,p);
  Eigen::MatrixXd V1_i(n,p);
  V1_i.setZero(n,p);
  //Calcualte S0 and S1 for beta_tilde
  //Eigen::MatrixXd theta_tilde(n,1);
  //theta_tilde.setZero(n,1);
  Eigen::MatrixXd exp_theta_tilde(n,1);
  exp_theta_tilde.setZero(n,1);
  Eigen::MatrixXd S0_tilde(n,1);
  S0_tilde.setZero(n,1);
  Eigen::MatrixXd S1_pre_tilde(n,p);
  S1_pre_tilde.setZero(n,p);
  Eigen::MatrixXd S1_tilde(n,p);
  S1_tilde.setZero(n,p);
  Eigen::MatrixXd S1_coli_tilde(n,1);
  S1_coli_tilde.setZero(n,1);

  //Calculate S0 with beta
  theta = Z*beta;
  exp_theta = (theta.array().exp()).matrix();
  List S0_rev_cumsum = rev_cumsum(exp_theta);
  S0 = S0_rev_cumsum["cumsum"];
  //Calcualte S0_tilde with beta_tilde
  //theta_tilde = Z*beta_tilde;
  exp_theta_tilde = (theta_tilde.array().exp()).matrix();
  List S0_rev_cumsum_tilde = rev_cumsum(exp_theta_tilde);
  S0_tilde = S0_rev_cumsum_tilde["cumsum"]; 
  
  for (int i = 0; i < p; i++){
   //Calculate S1 with beta
   S1_pre.col(i) = Z.col(i).cwiseProduct(exp_theta);
   List S1_col_rev_cumsum = rev_cumsum(S1_pre.col(i));
   S1_coli = S1_col_rev_cumsum["cumsum"];
   S1.col(i) = S1_coli;
   //Calcualte S1_tilde with beta_tilde
   S1_pre_tilde.col(i) = Z.col(i).cwiseProduct(exp_theta_tilde);
   List S1_col_rev_cumsum_tilde = rev_cumsum(S1_pre_tilde.col(i));
   S1_coli_tilde = S1_col_rev_cumsum_tilde["cumsum"];
   S1_tilde.col(i) = S1_coli_tilde;

   temp_1.col(i) = delta.cwiseProduct(Z.col(i)+eta*(S1_tilde.col(i).cwiseQuotient(S0_tilde))-(1+eta)*(S1.col(i).cwiseQuotient(S0)));
  } 

  L1 = temp_1.colwise().sum();
  temp_2 = delta.cwiseProduct(theta-S0.array().log().matrix());
  loglik = temp_2.sum();
  
  for (int i=0; i < p; i++){
    for (int j = 0; j < p; j++){
      S2_i1.col(j) = (Z.col(j).cwiseProduct(exp_theta)).cwiseProduct(Z.col(i));
      List S2i_col_rev_cumsum = rev_cumsum(S2_i1.col(j));
      //S1_coli = S1_col_rev_cumsum["cumsum"];
     // S1.col(i) = S1_coli;
      S2_i_colj = S2i_col_rev_cumsum["cumsum"];
      S2_i.col(j) = S2_i_colj;
      V1_i.col(j) = ((S2_i.col(j).cwiseQuotient(S0)))-((S1.col(j).cwiseProduct(S1.col(i))).cwiseQuotient(S0.array().square().matrix()));
      V_i.col(j) = V1_i.col(j).cwiseProduct(delta);
    } 
    L2.col(i) = (V_i.colwise().sum())*(1+eta);
  }

  List result;
  result["loglik"]=loglik;
  result["L1"]=L1;
  result["L2"]=L2;
  result["S0"]=S0;
  
  return result;
}

// [[Rcpp::export]]
List ddloglik(Eigen::MatrixXd Z, Eigen::MatrixXd delta, Eigen::MatrixXd beta){
  int p = beta.rows();
  int n = delta.rows();
  double loglik=0;
  Eigen::MatrixXd L1(p,1);
  L1.setZero(p,1);
  Eigen::MatrixXd L2(p,p);
  L2.setZero(p,p);
  Eigen::MatrixXd temp_2(n,1);
  temp_2.setZero(n,1); 
  Eigen::MatrixXd theta(n,1);
  theta.setZero(n,1);
  Eigen::MatrixXd exp_theta(n,1);
  exp_theta.setZero(n,1);
  Eigen::MatrixXd S0(n,1);
  S0.setZero(n,1);
  Eigen::MatrixXd S1_pre(n,p);
  S1_pre.setZero(n,p);
  Eigen::MatrixXd S1(n,p);
  S1.setZero(n,p);
  Eigen::MatrixXd S1_coli(n,1);
  S1_coli.setZero(n,1);
  Eigen::MatrixXd temp_1(n,p);
  temp_1.setZero(n,p); 
  Eigen::MatrixXd S2_i1(n,p);
  S2_i1.setZero(n,p);
  Eigen::MatrixXd S2_i(n,p);
  S2_i.setZero(n,p); 
  Eigen::MatrixXd S2_i_colj(n,p);
  S2_i_colj.setZero(n,p);   
  Eigen::MatrixXd V_i(n,p);
  V_i.setZero(n,p);
  Eigen::MatrixXd V1_i(n,p);
  V1_i.setZero(n,p);
  
  theta = Z*beta;
  exp_theta = (theta.array().exp()).matrix();
  List S0_rev_cumsum = rev_cumsum(exp_theta);
  S0 = S0_rev_cumsum["cumsum"];
  
  for (int i = 0; i < p; i++){
   S1_pre.col(i) = Z.col(i).cwiseProduct(exp_theta);
   List S1_col_rev_cumsum = rev_cumsum(S1_pre.col(i));
   S1_coli = S1_col_rev_cumsum["cumsum"];
   S1.col(i) = S1_coli;
   temp_1.col(i) = delta.cwiseProduct(Z.col(i)-(S1.col(i).cwiseQuotient(S0)));
  } 
  
  L1 = temp_1.colwise().sum();
  temp_2 = delta.cwiseProduct(theta-S0.array().log().matrix());
  loglik = temp_2.sum();
  
  for (int i=0; i < p; i++){
    for (int j = 0; j < p; j++){
      S2_i1.col(j) = (Z.col(j).cwiseProduct(exp_theta)).cwiseProduct(Z.col(i));
      List S2i_col_rev_cumsum = rev_cumsum(S2_i1.col(j));
      S2_i_colj = S2i_col_rev_cumsum["cumsum"];
      S2_i.col(j) = S2_i_colj;
      V1_i.col(j) = ((S2_i.col(j).cwiseQuotient(S0)))-((S1.col(j).cwiseProduct(S1.col(i))).cwiseQuotient(S0.array().square().matrix()));
      V_i.col(j) = V1_i.col(j).cwiseProduct(delta);
    } 
    L2.col(i) = V_i.colwise().sum();
  }

  List result;
  result["loglik"]=loglik;
  result["L1"]=L1;
  result["L2"]=L2;
  result["S0"]=S0;
  
  return result;
}

// [[Rcpp::export]]
List cal_S0(Eigen::MatrixXd theta, Eigen::MatrixXd delta){
  int n = theta.rows();
  Eigen::MatrixXd exp_theta(n,1);
  exp_theta.setZero(n,1);
  Eigen::MatrixXd S0(n,1);
  S0.setZero(n,1);
  
  exp_theta = (theta.array().exp()).matrix();
  List S0_rev_cumsum = rev_cumsum(exp_theta);
  S0 = S0_rev_cumsum["cumsum"];
  
  List result;
  result["S0"]=S0;
  
  return result;
}

