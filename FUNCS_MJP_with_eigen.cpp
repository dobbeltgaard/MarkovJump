
//#include <boost/math/special_functions/digamma.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>


#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>


using namespace std;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]



/* --------------------------------------- */
/* FUNCTIONS to parameterize the generator */
/* --------------------------------------- */
// function to construct generalized Erlang distribution
// #construct A matrix as generalized Erlang
// [[Rcpp::export]]
Eigen::MatrixXd make_A1(int m, const Eigen::VectorXd& lambda){
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ //iterate through rows
    A(i, i + 1) = lambda(i); //assign elements
  }
  A.diagonal() = -A.rowwise().sum(); //diag elements are equal to negative total rate of transmission
  return A;
}

// [[Rcpp::export]]
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

/* ************** */
/* Uniformization */
/* ************** */

inline double h(double x) {
  return 1.0-x+x*log(x);
}

inline double hifunc(double rho, double B,double loge) {
  return rho +(B-loge)/3.0*(1.0+sqrt(1.0+18.0*rho/(B-loge)));
}
inline double lofunc(double rho, double A, double loge) {
  const double logroot2pi=0.5*log(2*3.14159265);
  return rho+sqrt(2*rho)*sqrt(-logroot2pi-loge-1.5*log(A)+log(A-1));
}
unsigned int get_mlo(unsigned int mhi, double rho) {
  unsigned int mlo;
  // Since using unsigned int, need to be careful of negative numbers
  double dmlo=double(2*(int)(rho-0.5))-(double)mhi;
  if (dmlo>0) {
    mlo=(unsigned int) dmlo;
  }
  else {
    mlo=0;
  }
  
  return mlo;
}

unsigned int get_m(double rho, double prec, unsigned int mlo=0) {
  if (rho > 4.2e9) {
    return 0;
  }
  
  const double logprec=log(prec), pi=3.14159265;
  double dmhi, dmlo;
  unsigned int mhi;
  
  dmhi= hifunc(rho,0.0,logprec)-1;
  dmlo=lofunc(rho,2*rho*h(dmhi/rho),logprec);
  if ((unsigned int)dmlo > mlo) {
    mlo=(unsigned int)dmlo;
  }
  
  if (log(boost::math::gamma_p((double)(mlo+1),rho))<logprec) {
    return mlo; // lower bound is the value we want - no binary search needed
  }
  else {
    const double B=-0.5*log(4*pi*rho*h(dmlo/rho));
    if (B>logprec) {
      dmhi=hifunc(rho,B,logprec);
    }
    mhi=(unsigned int)(dmhi+1);
    
    //    cout<<mhi<<"-"<<mlo<<"=width: "<< mhi-mlo<<"\n";
    
    while (mhi-mlo>1) {
      unsigned int mmid=mlo+(mhi-mlo)/2; // =(mlo+mhi)/2, rounds down
      double dm=(double)mmid;
      double loginv;
      
      loginv=log(boost::math::gamma_p(dm+1,rho));
      //    cout <<mlo<<", "<<mmid <<", "<<mhi<<", "<<loginv<<"\n";
      
      if (loginv<logprec) {
        mhi=mmid;
      }
      else {
        mlo=mmid;
      }
    }
  }
  return mhi;
}

Eigen::RowVectorXd transient_dist_Uni(int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time){
  const double big = 1e100;  
  Eigen::MatrixXd Q = A*cov_time;
  Eigen::VectorXd diagonal = -Q.diagonal();
  const double rho = diagonal.maxCoeff();
  Eigen::MatrixXd M = Q; for(int i; i<m; i++){M(i,i) += rho; } // M = Q + rho*I
  const bool t2 = true; 
  const unsigned int mhi=get_m(rho,eps/(1.0+(double)t2)),mlo=t2? get_mlo(mhi,rho): 0;
  
  
  double b = Pt.lpNorm<1>(); 
  double c = 0.0;

  Eigen::RowVectorXd v_sum = Pt;
  Eigen::RowVectorXd v_pro = v_sum;
  if(b > big){
    v_pro /= b;
    v_sum /= b;
    c += log(b); 
    b = 1.0;
  }
  
  int f = int(1);
  for(int i; i < mhi; i++){
    v_pro = v_pro * M;
    v_pro /= f;
    b *= rho/f;
    v_sum += v_pro;
    if(b > big){
      v_pro /= b;
      v_sum /= b;
      c += log(b); 
      b = 1.0; 
    }
    f += int(1);
  }
  return exp(c-rho)*v_sum;
}


Eigen::ArrayXd exp_vec(const Eigen::ArrayXd& x) {
  return x.exp();
}

Eigen::ArrayXd soft_plus_vec(const Eigen::ArrayXd& x) {
  return (x.exp() + 1).log();
}

Eigen::ArrayXd square_vec(const Eigen::ArrayXd& x) {
  return x.square();
}

double exp_double(double x){
  return exp(x);
}
double soft_plus_double(double x){
  return log(exp(x)+1);
}
double square_double(double x){
  return pow(x,2);
}


/* ********************************************************** */
/* FUNCTION: SCORE for Markov jump process given partial data */
/* ********************************************************** */

// [[Rcpp::export]]
double MJP_score(int m,
                 const Eigen::VectorXd& s1,
                 const Eigen::VectorXd& s2,
                 const Eigen::VectorXd& u,
                 const Eigen::VectorXd& pars,
                 const Eigen::MatrixXd& z,
                 const string& generator = "gerlang",
                 const string& link_type_base = "exp",
                 const string& link_type_covs = "exp",
                 bool covs_bin = true,
                 bool likelihood_bin = false,
                 bool rps_bin = false,
                 bool brier_bin = false,
                 const string& transient_dist_method = "pade",
                 double eps = 2^(-52)){

  //Rcpp::Rcout << pars << std::endl;

  /* Initialization */
  int n = u.size();
  Eigen::VectorXd lambda_base;
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
  
  /* Determine which link to use */
  std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)> link_function_base;
  if (link_type_base == "exp") {
    link_function_base = exp_vec;
  } else if (link_type_base == "softplus") {
    link_function_base = soft_plus_vec;
  } else if (link_type_base == "square") {
    link_function_base = square_vec;
  }
  std::function<double(const double&)> link_function_covs;
  if (link_type_covs == "exp") {
    link_function_covs = exp_double;
  } else if (link_type_covs == "softplus") {
    link_function_covs = soft_plus_double;
  } else if (link_type_covs == "square") {
    link_function_covs = square_double;
  }

  /* Determine which parameterization to use */
  std::function<Eigen::MatrixXd(int, const Eigen::VectorXd&)> make_A;
  if (to_lowercase(generator) == "gerlang") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A1(m, lambda); };
    lambda_base.resize(int(m-1));
    lambda_base = link_function_base(pars.segment(0,m-1).array()); 
    U = eigenspace_U(m, lambda_base);
    U_inv = eigenspace_U_inv(m, lambda_base);
    A = make_A(m, lambda_base);
    for(int i; i < (m-1); i++){D(i,i) = -lambda_base(i); }
    bool distinct_rates = are_rates_distinct(m, lambda_base, 0.00000001); 
    if (!distinct_rates) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "gerlang_relax") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A2(m, lambda); };
    lambda_base.resize(int(m-1));
    lambda_base = link_function_base(pars.segment(0,m-1).array()); 
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
    lambda_base.resize(int(m*(m-1)/2));
    lambda_base = link_function_base(pars.segment(0,int(m*(m-1)/2)).array()); 
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
  if (to_lowercase(transient_dist_method) == "uniformization") {
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Uni(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else if (to_lowercase(transient_dist_method) == "pade") {
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else if (to_lowercase(transient_dist_method) == "eigenvalue_decomp" && eigen_solver_good){
        //Rcpp::Rcout << "Using eigenvalue decomp." << std::endl;
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Eig(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else { 
      transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
  }
  
  /* Compute score */
  if (!covs_bin) { // Case: Covariates excluded
    for(int i = 0; i < n; i++){
      cov_time = u(i);
      start_idx = int(s1(i)-1);
      end_idx = int(s2(i)-1);
      Pt.setZero();
      obs.setZero();
      Pt(start_idx) = 1 ;
      obs(end_idx) = 1 ;
      Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
      res += score_function(m, Ptu, obs);
    }
  } else {  // Case: Covariates included
    int k = z.cols();
    Eigen::VectorXd beta_covs = pars.tail(k);
    for(int i = 0; i < n; i++){
      cov_time = link_function_covs( beta_covs.dot(z.row(i)) ) * u(i);
      start_idx = int(s1(i)-1);
      end_idx = int(s2(i)-1);
      Pt.setZero();
      obs.setZero();
      Pt( start_idx) = 1 ;
      obs(end_idx) = 1 ;
      Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
      res += score_function(m, Ptu, obs);
    }
  }
  return res/n;
}



/* ******************************************* */
/* FUNCTION: Forecast for Markov jump process  */
/* ******************************************* */

// [[Rcpp::export]]
Eigen::MatrixXd MJP_predict(int m,
                            const Eigen::VectorXd& s1,
                            const Eigen::VectorXd& u,
                            const Eigen::VectorXd& pars,
                            const Eigen::MatrixXd& z,
                            const string& generator = "gerlang",
                            const string& link_type_base = "exp",
                            const string& link_type_covs = "exp",
                            bool covs_bin = true,
                            const string& transient_dist_method = "pade",
                            double eps = 2^(-52)){
  
  //Rcpp::Rcout << pars << std::endl;
  
  /* Initialization */
  int n = u.size();
  Eigen::VectorXd lambda_base;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd U_inv = Eigen::MatrixXd::Zero(m, m);
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m, m);
  Eigen::VectorXd lambda(m);
  Eigen::RowVectorXd Pt(m);
  Eigen::RowVectorXd Ptu(m);
  Eigen::MatrixXd solution(n,m);
  double cov_time = 0;
  int start_idx = 0;
  bool eigen_solver_good = true;
  
  /* Determine which link to use */
  std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)> link_function_base;
  if (link_type_base == "exp") {
    link_function_base = exp_vec;
  } else if (link_type_base == "softplus") {
    link_function_base = soft_plus_vec;
  } else if (link_type_base == "square") {
    link_function_base = square_vec;
  }
  std::function<double(const double&)> link_function_covs;
  if (link_type_covs == "exp") {
    link_function_covs = exp_double;
  } else if (link_type_covs == "softplus") {
    link_function_covs = soft_plus_double;
  } else if (link_type_covs == "square") {
    link_function_covs = square_double;
  }
  
  /* Determine which parameterization to use */
  std::function<Eigen::MatrixXd(int, const Eigen::VectorXd&)> make_A;
  if (to_lowercase(generator) == "gerlang") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A1(m, lambda); };
    lambda_base.resize(int(m-1));
    lambda_base = link_function_base(pars.segment(0,m-1).array()); 
    U = eigenspace_U(m, lambda_base);
    U_inv = eigenspace_U_inv(m, lambda_base);
    A = make_A(m, lambda_base);
    for(int i; i < (m-1); i++){D(i,i) = -lambda_base(i); }
    bool distinct_rates = are_rates_distinct(m, lambda_base, 0.00000001); 
    if (!distinct_rates) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "gerlang_relax") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A2(m, lambda); };
    lambda_base.resize(int(m-1));
    lambda_base = link_function_base(pars.segment(0,m-1).array()); 
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
    lambda_base.resize(int(m*(m-1)/2));
    lambda_base = link_function_base(pars.segment(0,int(m*(m-1)/2)).array()); 
    A = make_A(m, lambda_base);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
    Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
    Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
    D = D_complex.real();
    U = U_complex.real(); 
    U_inv = U.inverse();
    if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
  }
  
  /* Determine how to calculate transient distribution */
  std::function<Eigen::MatrixXd(int, const Eigen::RowVectorXd&, const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double)> transient_dist;
  if (to_lowercase(transient_dist_method) == "uniformization") {
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Uni(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else if (to_lowercase(transient_dist_method) == "pade") {
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else if (to_lowercase(transient_dist_method) == "eigenvalue_decomp" && eigen_solver_good){
    //Rcpp::Rcout << "Using eigenvalue decomp." << std::endl;
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Eig(m, Pt, A, eps, U, U_inv, D, cov_time); };
  } else { 
    transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
  }
  
  /* Compute score */
  if (!covs_bin) { // Case: Covariates excluded
    for(int i = 0; i < n; i++){
      cov_time = u(i);
      start_idx = int(s1(i)-1);
      Pt.setZero();;
      Pt(start_idx) = 1 ;
      Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
      solution.row(i) = Ptu;
    }
  } else {  // Case: Covariates included
    int k = z.cols();
    Eigen::VectorXd beta_covs = pars.tail(k);
    for(int i = 0; i < n; i++){
      cov_time = link_function_covs( beta_covs.dot(z.row(i)) ) * u(i);
      start_idx = int(s1(i)-1);
      Pt.setZero();
      Pt( start_idx) = 1 ;
      Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
      solution.row(i) = Ptu;
    }
  }
  return solution;
}



/* ******************************* */
/* FUNCTIONS: Evaluation functions */
/* ******************************* */

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
Eigen::VectorXd BrierScore_vectors(int m, const Eigen::MatrixXd& pred, const Eigen::MatrixXd& obs) {
  int n = pred.rows();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++){
    res(i) = brier_cpp(m, pred.row(i), obs.row(i));
  }
  return res;
}




































































