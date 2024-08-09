
//#include <boost/math/special_functions/digamma.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>

using namespace std;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]



/* --------------------------------------- */
/* FUNCTIONS to parameterize the generator */
/* --------------------------------------- */
// function to construct generalized Erlang distribution
// [[Rcpp::export]]
// #construct A matrix as generalized Erlang
Eigen::MatrixXd make_A1(int m, const Eigen::VectorXd& lambda){
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ //iterate through rows
    A(i, i + 1) = lambda(i); //assign elements
  }
  A.diagonal() = -A.rowwise().sum(); //diag elements are equal to negative total rate of transmission
  return A;
}

Eigen::MatrixXd make_A2(int m, const Eigen::VectorXd& lambda) {
  Eigen::MatrixXd A = make_A1(m, lambda);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      if (j - i > 1) {
        A(i, j) = (A(i, j - 1) * A(j - 1, j)) / (A(i, j - 1) + A(j - 1, j));
      }
    }
  }
  A.diagonal().setZero(); // Set diagonal elements to 0
  A.diagonal() = -A.rowwise().sum(); // Set diagonal elements to negative row sums
  
  return A;
}
// function to parameterize generator freely, filling all columns before changing row. m*(m-1)/2 free pars
// [[Rcpp::export]]
Eigen::MatrixXd make_A3(int m, const Eigen::VectorXd& lambda) {
  int count = 0;  // Initialize counter
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      if (i < j) {  
        A(i, j) = lambda(count);  
        count++;  
      }
    }
  }
  A.diagonal() = -A.rowwise().sum();
  return A;
}


/* ----------------------------------------- */
/* FUNCTIONS to compute probabilistic scores */
/* ----------------------------------------- */
// // [[Rcpp::export]]
// double rps_cpp(int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) {
//   double res = 0.0;
//   for(int k = 0; k < m; k++) {
//     double rps = 0.0;
//     for(int i = 0; i <= k; i++) {
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


double brier_cpp(int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) {
  Eigen::RowVectorXd diff = pred - obs;
  Eigen::RowVectorXd squared_diff = diff.array().square();
  double res = squared_diff.sum(); 
  return res;
}
// [[Rcpp::export]]
double log_cpp(int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) {
  Eigen::Index maxIndex;
  obs.maxCoeff(&maxIndex);
  return -log(pred(maxIndex));
}

/* ----------------------------------------- */
/* FUNCTIONS: Transient distribution methods */
/* ----------------------------------------- */

Eigen::RowVectorXd transient_dist_Pade(int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time){
  //function does not use eps, U and U_inv
  Eigen::MatrixXd At = A*cov_time;
  return Pt * At.exp();
}
// arma::rowvec transient_dist_Uni(int m, const arma::rowvec& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time){
//   //function does not use U and U_inv
//   Eigen::MatrixXd At = A*cov_time;
//   return (Pt,At,eps);
// }

/* eigenspace matrix */
Eigen::MatrixXd eigenspace_U(int m, const Eigen::VectorXd& lambda){
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
Eigen::MatrixXd eigenspace_U_inv(int m, const Eigen::VectorXd& lambda){
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



Eigen::RowVectorXd transient_dist_Eig(int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time){
  //function does not use eps
  Eigen::MatrixXd Delta = Eigen::MatrixXd::Zero(m, m);
  for(int i = 0; i < m; ++i){ Delta(i,i) = exp( D(i,i)*cov_time ) ; }
  Delta(m-1,m-1) = 1;
  Eigen::MatrixXd tpm = U * Delta * U_inv;
  return Pt * tpm;
}

bool are_rates_distinct(int m, const Eigen::VectorXd& x, double eps) {
  for (int i = 0; i < (m-1); ++i) {
    for (int j = i + 1; j < (m-1); ++j) {
      if (std::abs(x(i) - x(j)) <= eps) {
        return false;
      }
    }
  }
  return true;
}

std::string to_lowercase(const std::string& str) {
  std::string result = str;
  std::transform(result.begin(), result.end(), result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}


// [[Rcpp::export]]
double MJP_score(int m,
                 const Eigen::VectorXd& s1,
                 const Eigen::VectorXd& s2,
                 const Eigen::VectorXd& u,
                 const Eigen::VectorXd& pars,
                 const Eigen::MatrixXd& z,
                 const string& generator = "gerlang",
                 bool covs_bin = true,
                 bool likelihood_bin = true,
                 bool rps_bin = false,
                 bool brier_bin = false,
                 const string& transient_dist_method = "pade",
                 double eps = 2^(-52)){

  //Rcpp::Rcout << pars << std::endl;

  /* Initialization */
  int n = u.size();
  Eigen::VectorXd lambda_base = pars.segment(0,m-1).array().exp();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd U_inv = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::RowVectorXd obs(m);
  double res = 0;
  double cov_time = 0;
  int start_idx = 0;
  int end_idx = 0;
  bool eigen_solver_good = true;

  /* Determine which parameterization to use */
  std::function<Eigen::MatrixXd(int, const Eigen::VectorXd&)> make_A;
  if (to_lowercase(generator) == "gerlang") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A1(m, lambda); };
    U = eigenspace_U(m, lambda_base);
    U_inv = eigenspace_U_inv(m, lambda_base);
    A = make_A(m, lambda_base);
    for(int i; i < (m-1); i++){D(i,i) = -lambda_base(i); }
    bool distinct_rates = are_rates_distinct(m, lambda_base, 0.00000001); 
    if (!distinct_rates) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "gerlang_relax") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A2(m, lambda); };
    A = make_A(m, lambda_base);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
    Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
    Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
    D = D_complex.real();
    U = U_complex.real(); 
    U_inv = U.inverse();
    if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "free_upper_tri") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A3(m, lambda); };
    A = make_A(m, lambda_base);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
    Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
    Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
    D = D_complex.real();
    U = U_complex.real(); 
    U_inv = U.inverse();
    if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
  }

  /* Determine which score to use */
  std::function<double(int, const Eigen::RowVectorXd&, const Eigen::RowVectorXd&)> score_function;
  if (likelihood_bin && !rps_bin && !brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return log_cpp(m, pred, obs); };
  } else if (!likelihood_bin && rps_bin && !brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return rps_cpp(m, pred, obs); };
  } else if (!likelihood_bin && !rps_bin && brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return brier_cpp(m, pred, obs); };
  } else if (likelihood_bin && rps_bin && !brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return log_cpp(m, pred, obs) + rps_cpp(m, pred, obs); };
  } else if (likelihood_bin && !rps_bin && brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return log_cpp(m, pred, obs) + brier_cpp(m, pred, obs); };
  } else if (!likelihood_bin && rps_bin && brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return rps_cpp(m, pred, obs) + brier_cpp(m, pred, obs); };
  } else if (likelihood_bin && rps_bin && brier_bin) {
    score_function = [](int m, const Eigen::RowVectorXd& pred, const Eigen::RowVectorXd& obs) { return log_cpp(m, pred, obs) + rps_cpp(m, pred, obs) + brier_cpp(m, pred, obs); };
  } else {
    Rcpp::warning("A score metric needs to be specified.");
    return R_NaReal;
  }

  /* Determine how to calculate transient distribution */
  std::function<Eigen::MatrixXd(int, const Eigen::RowVectorXd&, const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double)> transient_dist;
  // if (to_lowercase(transient_dist_method) == "uniformization") {
  //   transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Uni(m, Pt, A, eps, U, U_inv, D, cov_time); };
  // } else if 
  if(to_lowercase(transient_dist_method) == "pade") {
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else if (to_lowercase(transient_dist_method) == "eigenvalue_decomp" && eigen_solver_good){
        //Rcpp::Rcout << "Using eigenvalue decomp." << std::endl;
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Eig(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else { 
      transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
  }
  // } else {
  //   Rcpp::warning("An appropiate calculation method for the transient distribution needs to be specified.");
  //   return R_NaReal;
  // }

  /* Check if covariates should be used */
  if (!covs_bin) { // Case: Covariates excluded
    for(int i = 0; i < n; i++){
      cov_time = u(i);
      //At = A;
      start_idx = int(s1(i)-1);
      end_idx = int(s2(i)-1);
      Pt.fill(0);
      obs.fill(0);
      Pt(start_idx) = 1 ;
      obs(end_idx) = 1 ;
      Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
      res += score_function(m, Ptu, obs);
    }
  } else {  // Case: Covariates included
    int k = z.cols();
    Eigen::VectorXd beta_covs = pars.segment(m-1,k);
    for(int i = 0; i < n; i++){
      cov_time = exp( beta_covs.dot(z.row(i)) ) * u(i);
      start_idx = int(s1(i)-1);
      end_idx = int(s2(i)-1);
      Pt.fill(0);
      obs.fill(0);
      Pt( start_idx) = 1 ;
      obs(end_idx) = 1 ;
      Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
      res += score_function(m, Ptu, obs);
    }
  }
  return res/n;
}
