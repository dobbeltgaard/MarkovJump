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
    
    Pt.setZero(); 
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    loglik += log( Ptu( end_idx  )  );
  }
  
  return -loglik;
}


/* eigenspace matrix */
Eigen::MatrixXd eigenspace_U(int m, Eigen::VectorXd lambda){
  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(m, m);
  U.triangularView<Eigen::StrictlyLower>().setZero();
  
  for(int i = 0; i < (m-1); i++){
    for(int j = 0; j < (m-1); j++){
      if(j > i){
        for(int k = i; k < j; k++){
          U(i, j) *= (lambda(k) / (lambda(k) - lambda(j)));
        }
      }
    }
  }
  return U;
}


/* Inverse eigenspace matrix */
Eigen::MatrixXd eigenspace_U_inv(int m, Eigen::VectorXd lambda){
  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(m, m);
  U.triangularView<Eigen::StrictlyLower>().setZero();
  
  for(int i = 0; i < (m-1); i++){
    for(int j = 0; j < (m-1); j++){
      if(j > i){
        for(int k = i; k < j; k++){
          U(i, j) *= (lambda(k) / (lambda(i) - lambda(k + 1)));
        }
        U(i, j) *= pow(-1, (j - i));
      }
    }
    if( i < (m-1) ){ //this condition is unnecessary
      for(int k = i+1; k < m-1; k++){
        U(i, m-1) *= (lambda(k) / (lambda(i) - lambda(k)));
      }
      U(i, m-1) *= pow(-1, (m - i + 1));
    }
  }
  return U;
}


/*Log-likelihood with eigen decomposition*/
// [[Rcpp::export]]
double discrete_loglik_eigen_cpp(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  int k = z.cols();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp(); //(Eigen::seq(0, m-2));
  Eigen::VectorXd beta_covs = pars.segment(m-1,k); //beta(Eigen::seq(m-1, beta.size()-1));
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  double loglik = 0;
  int start_idx = 0;
  int end_idx = 0;
  Eigen::MatrixXd U = eigenspace_U(m, lambda_base);
  Eigen::MatrixXd U_inv = eigenspace_U_inv(m, lambda_base);
  Eigen::MatrixXd Delta = Eigen::MatrixXd::Zero(m, m); 
  Delta(m-1,m-1) = 1;
  
  /* Compute the log-likelihood */
  for(int i = 0; i < n; i++){
    lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
    for(int j = 0; j < m-1; ++j){ Delta(j, j) = exp(-lambda(j)*u(i)); }
    tpm = U * Delta * U_inv;
    
    start_idx = int(s1(i)-1);
    end_idx = int(s2(i)-1);
    
    Pt.setZero(); 
    Pt( start_idx) = int(1);
    Ptu = Pt * tpm;
    
    loglik += log( Ptu( end_idx  )  );
  }
  
  return -loglik;
}
 

/* Gradient of eigenspace wrt. the h'th lambda */
// [[Rcpp::export]]
Eigen::MatrixXd eigenspace_U_grad(int m, Eigen::VectorXd lambda, int h){
  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(m, m);
  Eigen::MatrixXd Us = Eigen::MatrixXd::Zero(m, m);
  
  for(int i = 0; i < (m-1); i++){
    for(int j = 0; j < (m-1); j++){
      if(j > i){
        for(int k = i; k < j; k++){
          if(h == j){
            Us(i, j) += 1 / (lambda(k) - lambda(j)) ;
          }
          if(h == k){
            Us(i, j) += lambda(j) / (lambda(k) * (-lambda(k) + lambda(j))) ;
          }
          U(i, j) *=  lambda(k) / (lambda(k) - lambda(j)) ;
        }
      }
    }
  }
  Eigen::MatrixXd result = U.array() * Us.array(); 
  return result;
}
 

/* Gradient of inverse eigenspace*/
// [[Rcpp::export]]
Eigen::MatrixXd eigenspace_U_inv_grad(int m, Eigen::VectorXd lambda, int h){
  Eigen::MatrixXd U = Eigen::MatrixXd::Ones(m, m);
  Eigen::MatrixXd Us = Eigen::MatrixXd::Zero(m, m);
  
  for(int i = 0; i < (m-1); i++){
    for(int j = 0; j < (m-1); j++){
      if(j > i){
        for(int k = i; k < j; k++){
          if(h == i){
            Us(i, j) -= 1 / (lambda(i) - lambda(k + 1)) ;
          }
          if(h == k){
            Us(i, j) += 1 / lambda(k) ;
          }
          if(h == (k + 1)){
            Us(i, j) += 1 / (lambda(i) - lambda(k + 1)) ;
          }
          U(i, j) *=  lambda(k) / (lambda(i) - lambda(k + 1)) ;
        }
        U(i, j) *= pow(-1, (j - i));
      }
    }
    if(i < (m-1)){ //this condition is unnecessary
      for(int k = i+1; k < m-1; k++){
        if(h == i){
          Us(i, m-1) -= 1 / (lambda(i) - lambda(k)) ;
        }
        if(h == k){
          Us(i, m-1) += lambda(i) / (lambda(k) * (lambda(i) - lambda(k))) ;
        }
        U(i, m-1) *=  lambda(k) / (lambda(i) - lambda(k)) ;
      }
      U(i, m-1) *= pow(-1, (m - i + 1));
    }
  }
  Eigen::MatrixXd result = U.array() * Us.array(); 
  return result; 
  
}


  
/* Gradient of log-likelihood with eigen decomposition*/
// [[Rcpp::export]]
Eigen::VectorXd discrete_loglik_eigen_grad_cpp(int m, Eigen::VectorXd s1, Eigen::VectorXd s2, Eigen::VectorXd u, Eigen::MatrixXd z, Eigen::VectorXd pars){
  
  /* Initialization */
  int n = u.size();
  int k = z.cols();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp(); //(Eigen::seq(0, m-2));
  Eigen::VectorXd beta_covs = pars.segment(m-1,k); //beta(Eigen::seq(m-1, beta.size()-1));
  Eigen::MatrixXd tpm = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd tpm_grad = Eigen::MatrixXd::Zero(m, m);
  
  Eigen::VectorXd lambda(m);
  double loglik = 0;
  int start_idx = 0;
  int end_idx = 0;
  Eigen::MatrixXd U = eigenspace_U(m, lambda_base);
  Eigen::MatrixXd U_inv = eigenspace_U_inv(m, lambda_base);
  Eigen::MatrixXd Delta = Eigen::MatrixXd::Zero(m, m); 
  Delta(m-1,m-1) = 1;
  Eigen::MatrixXd U_grad = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd U_inv_grad = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd lambda_grad = Eigen::MatrixXd::Zero(m, m);
  //Eigen::MatrixXd beta_grad = Eigen::MatrixXd::Zero(m, k);
  //Eigen::MatrixXd z_grad = Eigen::MatrixXd::Zero(n, k);
  Eigen::VectorXd grad(m-1+k);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu1(m);
  Eigen::RowVectorXd Ptu2(m);
  
  
  /* Loop over base parameters */
  for(int h = 0; h < m-1; ++h){
    double foo = 0;
    
    Eigen::MatrixXd Delta_gradient = Eigen::MatrixXd::Zero(m, m);
    U_grad = eigenspace_U_grad(m, lambda_base, h);
    U_inv_grad = -U_inv * U_grad * U_inv;//eigenspace_U_inv_grad(m, lambda_base, h);
    
    for(int i = 0; i < n; i++){
      lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
      for(int j = 0; j < m-1; ++j){ 
        Delta(j, j) = exp(-lambda(j)*u(i)); 
        }
      Delta_gradient(h,h) = -u(i) * exp( beta_covs.dot(z.row(i)) ) * exp(-lambda(h)*u(i)); // why Delta(h,h)?
      
      tpm_grad = (U_grad * Delta * U_inv) + (U * Delta_gradient * U_inv) + (U * Delta * U_inv_grad);
      tpm = U * Delta * U_inv;
      
      start_idx = int(s1(i)-1);
      end_idx = int(s2(i)-1);
      Pt.setZero(); 
      Pt( start_idx) = int(1);
      Ptu1 = Pt * tpm_grad;
      Ptu2 = Pt * tpm;
      
      foo +=  Ptu1( end_idx  )   / Ptu2( end_idx ) ;
    }
    grad(h) = foo;
    
  }
  
  /* Loop over covariate parameters */
  for(int h = 0; h < k; ++h){
    double foo = 0;
    
    Eigen::MatrixXd Delta_gradient = Eigen::MatrixXd::Zero(m, m);
    Delta_gradient(m-1,m-1) = 0;

    for(int i = 0; i < n; i++){
      lambda = lambda_base * exp( beta_covs.dot(z.row(i)) );
      for(int j = 0; j < m-1; ++j){ 
        Delta(j, j) = exp(-lambda(j)*u(i)); 
        Delta_gradient(j,j) = -lambda(j)*u(i)*z(i,h) * Delta(j,j); // = Delta_grad(m, lambda, u(i), Delta, h);
        }
      
      tpm_grad = U * Delta_gradient * U_inv;
      tpm = U * Delta * U_inv;
      
      start_idx = int(s1(i)-1);
      end_idx = int(s2(i)-1);
      Pt.setZero(); 
      Pt( start_idx) = int(1);
      Ptu1 = Pt * tpm_grad;
      Ptu2 = Pt * tpm;
      
      foo +=  Ptu1( end_idx  )   / Ptu2( end_idx ) ;
    }
    grad(h + m-1) = foo;
    
  }
  
  
  return -grad;
}


 
 
 
 
 
 
 
 
 
 
 
 
 
 