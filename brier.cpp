//#include <boost/math/special_functions/digamma.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]


/*Brier score function */
// [[Rcpp::export]]
double brier_cpp(const Eigen::VectorXd& pred, const Eigen::VectorXd& obs) {
  Eigen::VectorXd diff = pred - obs;
  Eigen::VectorXd squared_diff = diff.array().square();
  double res = squared_diff.sum(); 
  return res;
}


// #construct A matrix as generalized Erlang
// [[Rcpp::export]]
Eigen::MatrixXd make_A(int m, Eigen::VectorXd lambda){
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ //iterate through rows
    A(i, i + 1) = lambda(i); //assign elements
  }
  A.diagonal() = -A.rowwise().sum(); 
  return A;
}



// [[Rcpp::export]]
double discrete_brier_cpp_cov(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
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
  
  double brier = 0;
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
    
    
    brier += brier_cpp(Ptu, Ptu_obs);
  }
  
  return brier/n;
}


// [[Rcpp::export]]
double discrete_brier_cpp(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  Eigen::VectorXd lambda = pars.segment(0,m-1).array().exp(); 
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::RowVectorXd Ptu_obs(m);
  
  double brier = 0;
  int start_idx = 0;
  int end_idx = 0;
  
  /* Compute the log-likelihood */
  for(int i = 0; i < n; i++){
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
    
    
    brier += brier_cpp(Ptu, Ptu_obs);
  }
  
  return brier/n;
}


// [[Rcpp::export]]
double discrete_brier_loglik_cpp_cov(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
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
  
  double brier = 0;
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
    brier += brier_cpp(Ptu, Ptu_obs);
  }
  
  return (brier-loglik)/n;
}


