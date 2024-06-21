/* load: library(RcppEigen); library(Rcpp)  */ 
/* to compile, run: sourceCpp("funcs2.cpp") */ 

//#include <boost/math/special_functions/digamma.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


/*RPS function*/
// [[Rcpp::export]]
double rps_cpp(int m, Eigen::VectorXd pred, Eigen::VectorXd obs){
  double res = 0.0;
  for(int k = 0; k < m; k++){
    double rps = 0.0;
    for(int i = 0; i <= k; i++){
      rps += pred[i] - obs[i];
    }
    res += rps * rps;
  }
  return res;
}
// // [[Rcpp::export]]
// double rps_cpp(int m, Eigen::VectorXd pred, Eigen::VectorXd obs) {
//   double res = 0.0;
//   for (int k = 0; k < m; k++) {
//     double cum_pred = 0.0;
//     double cum_obs = 0.0;
//     for (int i = 0; i <= k; i++) {
//       cum_pred += pred[i];
//       cum_obs += obs[i];
//     }
//     double rps = cum_pred - cum_obs;
//     res += rps * rps;
//   }
//   return res;
// }


// #construct A matrix as generalized Erlang
Eigen::MatrixXd make_A(int m, Eigen::VectorXd lambda){
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ //iterate through rows
    A(i, i + 1) = lambda(i); //assign elements
  }
  A.diagonal() = -A.rowwise().sum(); //diag elements are equal to negative total rate of transmission
  return A;
}


// [[Rcpp::export]]
double discrete_rps_cpp_cov(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  int k = z.cols();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp(); 
  Eigen::VectorXd beta_covs = pars.segment(m-1,k); 
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::RowVectorXd Ptu_obs(m);
  
  double rps = 0;
  int start_idx = 0;
  int end_idx = 0;
  
  /* Compute the log-likelihood */
  for(int i = 0; i < n; i++){
    lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    At = make_A(m, lambda )*u(i);
    tpm = At.exp();
    
    start_idx = int(s1(i)-1);
    end_idx = int(s2(i)-1);
    
    /*Prediction*/
    Pt.setZero();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    /*Observation*/
    Ptu_obs.setZero();
    Ptu_obs(end_idx) = int(1);
    
    
    rps += rps_cpp(m, Ptu, Ptu_obs);
  }
  
  return rps/n;
}


// [[Rcpp::export]]
double discrete_rps_cpp(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  Eigen::VectorXd lambda = pars.segment(0,m-1).array().exp(); 
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::RowVectorXd Ptu_obs(m);
  
  double rps = 0;
  int start_idx = 0;
  int end_idx = 0;
  
  /* Compute the log-likelihood */
  for(int i = 0; i < n; i++){
    //lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    At = make_A(m, lambda )*u(i);
    tpm = At.exp();
    
    start_idx = int(s1(i)-1);
    end_idx = int(s2(i)-1);
    
    /*Prediction*/
    Pt.setZero();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    /*Observation*/
    Ptu_obs.setZero();
    Ptu_obs(end_idx) = int(1);
    
    
    rps += rps_cpp(m, Ptu, Ptu_obs);
  }
  
  return rps/n;
}


// [[Rcpp::export]]
double discrete_rps_loglik_cpp_cov(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  int k = z.cols();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp(); 
  Eigen::VectorXd beta_covs = pars.segment(m-1,k); 
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::RowVectorXd Ptu_obs(m);
  
  double rps = 0;
  double loglik = 0;
  int start_idx = 0;
  int end_idx = 0;
  
  /* Compute the log-likelihood */
  for(int i = 0; i < n; i++){
    lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    At = make_A(m, lambda )*u(i);
    tpm = At.exp();
    
    start_idx = int(s1(i)-1);
    end_idx = int(s2(i)-1);
    
    /*Prediction*/
    Pt.setZero();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    /*Observation*/
    Ptu_obs.setZero();
    Ptu_obs(end_idx) = int(1);
    
    loglik += log( Ptu( end_idx  )  );
    rps += rps_cpp(m, Ptu, Ptu_obs);
  }
  
  return (rps-loglik)/n;
}

// [[Rcpp::export]]
double discrete_rps_skill_cpp_cov(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars, Eigen::MatrixXd base){
  
  /* Initialization */
  int n = u.size();
  int k = z.cols();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp(); 
  Eigen::VectorXd beta_covs = pars.segment(m-1,k); 
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::RowVectorXd Ptu_obs(m);
  
  double rps = 0;
  double rps_base = 0; 
  double rpss = 0;
  int start_idx = 0;
  int end_idx = 0;
  
  /* Compute the log-likelihood */
  for(int i = 0; i < n; i++){
    lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    At = make_A(m, lambda )*u(i);
    tpm = At.exp();
    
    start_idx = int(s1(i)-1);
    end_idx = int(s2(i)-1);
    
    /*Prediction*/
    Pt.setZero();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    /*Observation*/
    Ptu_obs.setZero();
    Ptu_obs(end_idx) = int(1);
    
    //loglik += log( Ptu( end_idx  )  );
    rps = rps_cpp(m, Ptu, Ptu_obs);
    rps_base = rps_cpp(m, base.row(i), Ptu_obs);
    
    rpss += rps/(rps_base + 1);
    //rpss += rps_base;
    
  }
  
  return rpss;
}

