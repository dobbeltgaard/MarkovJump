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




// #construct A matrix as generalized Erlang
Eigen::MatrixXd make_A(int m, Eigen::VectorXd lambda){
  int count = 0; //init counter
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ //iterate through rows
        A(i, i + 1) = lambda(i); //assign elements
  }
  A.diagonal() = -A.rowwise().sum(); //diag elements are equal to negative total rate of transmission
  return A;
}

/*Log-likelihood*/
// [[Rcpp::export]]
double discrete_loglik_cpp(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  int k = z.cols();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp(); //(Eigen::seq(0, m-2));
  Eigen::VectorXd beta_covs = pars.segment(m-1,k); //beta(Eigen::seq(m-1, beta.size()-1));
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
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
    
    Pt.setZero(); //(s1(i) == states).cast<int>();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    loglik += log( Ptu( end_idx  )  );
  }
  
  return -loglik;
}




