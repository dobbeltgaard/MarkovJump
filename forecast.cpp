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

/* Prediction function */
// // [[Rcpp::export]]
// Eigen::VectorXd jump_prediction(Eigen::VectorXd Pt, double u, Eigen::MatrixXd A){
//   Eigen::MatrixXd At = A*u;
//   Eigen::MatrixXd tpm = At.exp();
//   Eigen::VectorXd Ptu = Pt * tpm;
//   return Ptu;
// }

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


// Predict probability vector by Markov Jump process with covariates
// [[Rcpp::export]]
Eigen::MatrixXd jump_prediction_cov(int m, Eigen::VectorXd s1, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){

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
  Eigen::MatrixXd solution(n,m);
  
  int start_idx = 0;

  for(int i = 0; i < n; i++){
    lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    At = make_A(m, lambda )*u(i);
    tpm = At.exp();
    
    start_idx = int(s1(i)-1);

    /*Prediction*/
    Pt.setZero();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    solution.row(i) = Ptu;
    
  }
  
  return solution;
  
}


// Predict probability vector by Markov Jump process with covariates
// [[Rcpp::export]]
Eigen::MatrixXd jump_prediction(int m, Eigen::VectorXd s1, Eigen::VectorXd u, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  Eigen::VectorXd lambda = pars.segment(0,m-1).array().exp(); 
  //Eigen::VectorXd beta_covs = pars.segment(m-1,k); 
  Eigen::MatrixXd At = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  //Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::MatrixXd solution(n,m);
  
  int start_idx = 0;
  
  for(int i = 0; i < n; i++){
    //lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    At = make_A(m, lambda )*u(i);
    tpm = At.exp();
    
    start_idx = int(s1(i)-1);
    
    /*Prediction*/
    Pt.setZero();
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    solution.row(i) = Ptu;
    
  }
  
  return solution;
  
}


// Predict probability vector with uniform distribution
// [[Rcpp::export]]
Eigen::MatrixXd uniform_prediction(int m, Eigen::VectorXd s1){
  int n = s1.size();
  Eigen::RowVectorXd Pt(m);
  Eigen::MatrixXd solution(n,m);
  int start_idx = 0;
  for(int i = 0; i < n; i++){
    start_idx = int(s1(i)-1);
    Pt.setZero();
    for(int j = 0; j < m; j++){
      if(j >= start_idx){
        Pt(j) = int(1);
      }
    }
    solution.row(i) = Pt / Pt.sum();
  }
  return solution;
}


// Create observation vectors
// [[Rcpp::export]]
Eigen::MatrixXd make_Ptu_obs(int m, Eigen::VectorXd s2){
  
  int n = s2.size();
  Eigen::RowVectorXd Ptu_obs(m);
  Eigen::MatrixXd solution(n,m);
  int idx = 0;
  
  for(int i = 0; i < n; i++){
    idx = int(s2(i)-1);
    Ptu_obs.setZero();
    Ptu_obs( idx) = int(1);
    solution.row(i) = Ptu_obs;
  }
  return solution;
}



/*RPS score function*/
// // [[Rcpp::export]]
// double rps_cpp(int m, Eigen::VectorXd pred, Eigen::VectorXd obs){
//   double res = 0.0; 
//   for(int k = 0; k < m; k++){
//     double rps = 0.0;
//     for(int i = 0; i <= k; i++){
//       rps += pred[i] - obs[i];
//     }
//     res += rps * rps; 
//   }
//   return res;
// }
// [[Rcpp::export]]
double rps_cpp(int m, Eigen::VectorXd pred, Eigen::VectorXd obs) {
  double res = 0.0;
  for (int k = 0; k < m; k++) {
    double cum_pred = 0.0;
    double cum_obs = 0.0;
    for (int i = 0; i <= k; i++) {
      cum_pred += pred[i];
      cum_obs += obs[i];
    }
    double rps = cum_pred - cum_obs;
    res += rps * rps;
  }
  return res;
}


/*RPS score function with vector input*/
// [[Rcpp::export]]
Eigen::VectorXd rps_vectors(int m, Eigen::MatrixXd pred, Eigen::MatrixXd obs){
  int n = pred.rows();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++){
    res(i) = rps_cpp(m, pred.row(i), obs.row(i)); 
  }
  return res;
}

/*log score function with vector input*/
// [[Rcpp::export]]
Eigen::VectorXd logscore_vectors(int m, Eigen::MatrixXd pred, Eigen::MatrixXd obs){
  int n = pred.rows();
  Eigen::VectorXd res(n);
  bool cond = true; 
  
  for(int i = 0; i < n; i++){
    
    // Find 1-index
    int j = 0;
    while(cond){
      if(obs(i,j) == 1){
        cond = false;
      }
      j++;
    }
    
    res(i) = -log(pred(i,j-1));
    cond = true;
  }
  return res;
}

/*Brier score function with vector input*/
// [[Rcpp::export]]
Eigen::VectorXd BrierScore_vectors(const Eigen::MatrixXd& pred, const Eigen::MatrixXd& obs) {
  Eigen::MatrixXd diff = pred - obs;
  Eigen::MatrixXd squared_diff = diff.array().square();
  Eigen::VectorXd res = squared_diff.rowwise().mean(); 
  return res;
}






