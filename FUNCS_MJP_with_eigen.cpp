// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]   // or: cpp14

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
  for(int i = 0; i < (m-1); i++){ A(i, i + 1) = lambda(i); }
  A.diagonal() = -A.rowwise().sum();
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

// [[Rcpp::export]]
Eigen::MatrixXd make_A4(int m, const Eigen::VectorXd& lambda){
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ A(i, i + 1) = lambda(i); }
  for(int i = 0; i < (m-2); i++){ A(i, i + 2) = lambda(i+m-1); }
  A.diagonal() = -A.rowwise().sum();
  return A;
}

// [[Rcpp::export]]
Eigen::MatrixXd make_A5(int m, const Eigen::VectorXd& lambda){
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m); //init A
  for(int i = 0; i < (m-1); i++){ A(i, i + 1) = lambda(i); }
  for(int i = 0; i < (m-2); i++){ A(i, i + 2) = lambda(i+m-1); }
  for(int i = 0; i < (m-3); i++){ A(i, i + 3) = lambda(i+2*m-3); }
  A.diagonal() = -A.rowwise().sum();
  return A;
}

// // [[Rcpp::export]]
// Eigen::MatrixXd expand_A(int m, Eigen::MatrixXd A, int k){
//   int kp1 = k + 1;
//   int M = (m - 1) * k + m - k;
//   Eigen::MatrixXd B = Eigen::MatrixXd::Zero(M, M);
//   for(int i = 0; i < m; ++i){
//     int rowstart = i * kp1;
//     int rowend = rowstart + k;
//     for(int j = 0; j < m; ++j){
//       double aij = A(i, j);
//       if(j - i == 1 && i < m - 2){
//         for(int t = 0; t <= k; ++t){
//           int r = rowstart + t;
//           int c = rowstart + 1 + t;
//           B(r, c) = aij;
//         }
//       }
//       if(j - i == 1 && i == m - 2){
//         int r = rowstart;
//         int c = rowstart + 1;
//         B(r, c) = aij;
//       }
//       else if(j - i > 1 && j < m - 1){
//         int colstart = 1 + kp1 * (j - 1);
//         for(int t = 0; t <= k; ++t){
//           int r = rowstart + t;
//           int c = colstart;
//           B(r, c) = aij;
//         }
//       }
//       else if(j - i > 1 && j == m - 1){
//         int col = 1 + kp1 * (j - 1);
//         for(int t = 0; t <= k; ++t){
//           int r = rowstart + t;
//           B(r, col) = aij;
//         }
//       }
//     }
//   }
//   Eigen::VectorXd rs = B.rowwise().sum();
//   for(int r = 0; r < M; ++r) B(r, r) = -rs(r);
//   return B;
// }

// [[Rcpp::export]]
Eigen::MatrixXd expand_A1(int m, const Eigen::MatrixXd& A, int k){
  const int kp1 = k + 1;
  const int M   = (m - 1) * kp1 + 1;        // last state is single absorbing
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(M, M);

  for(int i = 0; i < m - 1; ++i){
    const int rowstart = i * kp1;

    const double a = A(i, i + 1);
    for(int t = 0; t <= k; ++t){
      const int r = rowstart + t;
      const int c = rowstart + t + 1;
      if(r < M && c < M) B(r, c) = a;
    }
  }

  B.diagonal() = -B.rowwise().sum();
  return B;
}
// 
// 
// // [[Rcpp::export]]
// Eigen::MatrixXd expand_A1(int m, Eigen::MatrixXd A, int k){
//   int kp1 = k + 1;
//   int M = (m - 1) * (k + 1) + 1; //(m - 1) * k + m - k+k;
//   Eigen::MatrixXd B = Eigen::MatrixXd::Zero(M, M);
//   for(int i = 0; i < m; ++i){
//     int rowstart = i * kp1;
//     int rowend = rowstart + k;
//     for(int j = 0; j < m; ++j){
//       double aij = A(i, j);
//       if(j - i == 1 ){//}&& i < m - 2){
//         for(int t = 0; t <= k; ++t){
//           int r = rowstart + t;
//           int c = rowstart + 1 + t;
//           B(r, c) = aij;
//         }
//       }
//     }
//   }
//   Eigen::VectorXd rs = B.rowwise().sum();
//   for(int r = 0; r < M; ++r) B(r, r) = -rs(r);
//   return B;
// }
inline int expanded_size(int m, int k){
  if(m <= 1) return 1;
  return (m - 1) * (k + 1) + 1;
}

// [[Rcpp::export]]
Eigen::RowVectorXd expand_dist(int m, Eigen::RowVectorXd Pt, int k){
  int M = expanded_size(m, k);
  Eigen::RowVectorXd y = Eigen::RowVectorXd::Zero(M);
  int pos = 0;
  for(int i = 0; i < m - 1; ++i){
    int b = k + 1;
    double v = Pt(i) / static_cast<double>(b);
    for(int r = 0; r < b; ++r) y(pos + r) = v;
    pos += b;
  }
  y(pos) = Pt(m - 1);
  return y;
}

// [[Rcpp::export]]
Eigen::RowVectorXd collapse_dist(int m, Eigen::RowVectorXd Pt, int k){
  Eigen::RowVectorXd z = Eigen::RowVectorXd::Zero(m);
  int pos = 0;
  for(int i = 0; i < m - 1; ++i){
    int b = k + 1;
    z(i) = Pt.segment(pos, b).sum();
    pos += b;
  }
  z(m - 1) = Pt(pos);
  return z;
}

Eigen::RowVectorXd trans_dist_append_states(int m, Eigen::RowVectorXd Pt, int k, Eigen::MatrixXd A, double u){
  Eigen::MatrixXd B = expand_A1(m, A, k);
  Eigen::RowVectorXd init = expand_dist(m, Pt, k);
  Eigen::MatrixXd Bt = B*u;
  Eigen::RowVectorXd Ptu = init * Bt.exp();
  return collapse_dist(m, Ptu, k); 
}



// // [[Rcpp::export]]
// Eigen::RowVectorXd expand_dist_weighted(
//     int m, 
//     const Eigen::RowVectorXd& Pt, 
//     int k, 
//     const Eigen::MatrixXd& B
// ) {
//   int M = expanded_size(m, k);
//   Eigen::RowVectorXd y = Eigen::RowVectorXd::Zero(M);
//   int pos = 0;
//   for (int i = 0; i < m; ++i) {
//     if (i < m - 2) {
//       int b = k + 1;
//       Eigen::MatrixXd Si = B.block(pos, pos, b, b);
//       Eigen::MatrixXd Minv = (-Si).inverse();
//       Eigen::RowVectorXd e1 = Eigen::RowVectorXd::Zero(b);
//       e1(0) = 1.0;
//       Eigen::RowVectorXd w = e1 * Minv;
//       double s = w.sum();
//       if (s > 0) w /= s;
//       y.segment(pos, b) = Pt(i) * w;
//       pos += b;
//     } else {
//       y(pos) = Pt(i);
//       pos += 1;
//     }
//   }
//   return y;
// }






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
Eigen::RowVectorXd expected_sojourn(int m, double u, Eigen::RowVectorXd p, Eigen::MatrixXd A, double dt = 0.005) {
  Eigen::RowVectorXd mu = Eigen::RowVectorXd::Zero(m);
  int K = static_cast<int>(u / dt);
  for (int k = 0; k < K; k++) {
    mu += p;
    p += dt * (p * A).eval();
  }
  mu *= dt;
  return mu;
}

// [[Rcpp::export]]
Eigen::RowVectorXd expected_sojourn_exact(int m, double u, Eigen::RowVectorXd p, Eigen::MatrixXd A) {
  Eigen::MatrixXd S = A.topLeftCorner(m - 1, m - 1);
  Eigen::MatrixXd expSu = (S * u).exp();
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m - 1, m - 1);
  Eigen::MatrixXd rhs = I - expSu;
  Eigen::MatrixXd F = (-S).fullPivLu().solve(rhs);
  Eigen::RowVectorXd alpha = p.head(m - 1);
  Eigen::RowVectorXd mu_trans = alpha * F;
  double mu_abs = u - mu_trans.sum();
  Eigen::RowVectorXd mu(m);
  mu.head(m - 1) = mu_trans;
  mu(m - 1) = mu_abs;
  return mu;
}

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
                 double eps = 2.220446049250313e-16,
                 bool warping = false){

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
    //lambda_base.resize(int(m-1));
    lambda_base = link_function_base(pars.segment(0,m-1).array()); 
    U = eigenspace_U(m, lambda_base);
    U_inv = eigenspace_U_inv(m, lambda_base);
    A = make_A(m, lambda_base);
    for(int i = 0; i < (m-1); i++){D(i,i) = -lambda_base(i); }
    bool distinct_rates = are_rates_distinct(m, lambda_base, 0.00000001); 
    if (!distinct_rates) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "gerlang_relax") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A2(m, lambda); };
    //lambda_base.resize(int(m-1));
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
    //lambda_base.resize(int(m*(m-1)/2));
    lambda_base = link_function_base(pars.segment(0,int(m*(m-1)/2)).array()); 
    A = make_A(m, lambda_base);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
    Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
    Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
    D = D_complex.real();
    U = U_complex.real(); 
    U_inv = U.inverse();
    if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "bidiagonal") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A4(m, lambda); };
    //lambda_base.resize(int(2*m-3));
    lambda_base = link_function_base(pars.segment(0,int(2*m-3)).array()); 
    A = make_A(m, lambda_base);
    Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
    Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
    Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
    D = D_complex.real();
    U = U_complex.real(); 
    U_inv = U.inverse();
    if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
  } else if (to_lowercase(generator) == "tridiagonal") {
    make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A5(m, lambda); };
    //lambda_base.resize(int(3*m-6));
    lambda_base = link_function_base(pars.segment(0,int(3*m-6)).array()); 
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
    if(!warping){
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
    } else{ // warping
      Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
      Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);  // initialize with 1s
      xi.head(m - 1) = exp_vec(xii); //soft_plus_vec(xii);
      for(int i = 0; i < n; i++){
        start_idx = int(s1(i)-1);
        end_idx = int(s2(i)-1);
        Pt.setZero();
        obs.setZero();
        Pt(start_idx) = 1 ;
        obs(end_idx) = 1 ;
        Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
        double tau_eff = 0;
        for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);
        cov_time = tau_eff; 
        Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
        res += score_function(m, Ptu, obs);
      }
    }
  } else {  // Case: Covariates included
    int k = z.cols();
    Eigen::VectorXd beta_covs = pars.tail(k);
    if(!warping){
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
    } else { //warping
      Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
      Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);  // initialize with 1s
      xi.head(m - 1) = exp_vec(xii); //soft_plus_vec(xii);            // overwrite first m-1 entries
      for(int i = 0; i < n; i++){
        start_idx = int(s1(i)-1);
        end_idx = int(s2(i)-1);
        Pt.setZero();
        Pt(start_idx) = 1.0;
        obs.setZero();
        obs(end_idx) = 1 ;
        Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
        double tau_eff = 0;
        for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);
        cov_time = link_function_covs( beta_covs.dot(z.row(i)) ) * tau_eff; 
        Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
        res += score_function(m, Ptu, obs);
      }
    }
  }
  return res/n;
}



/* ******************************************* */
/* FUNCTION: Forecast for Markov jump process  */
/* ******************************************* */

inline int base_len_from_generator(int m, const std::string& gen) {
  std::string g = to_lowercase(gen);
  if (g == "gerlang" || g == "gerlang_relax") return m - 1;
  if (g == "free_upper_tri") return int(m * (m - 1) / 2);
  if (g == "bidiagonal")     return int(2 * m - 3);
  if (g == "tridiagonal")    return int(3 * m - 6);
  return -1;
}

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
                            double eps = 2.220446049250313e-16,
                            bool warping = false,
                            bool state_covs = false,
                            int k = 1, bool append = false, bool mixture = false, int K = 1) {
  
  const int n = (int)u.size();
  
  auto link_base = (link_type_base == "softplus") ? std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)>(soft_plus_vec)
    : (link_type_base == "square")   ? std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)>(square_vec)
      : std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)>(exp_vec);
  
  auto link_cov  = (link_type_covs == "softplus") ? std::function<double(const double&)>(soft_plus_double)
    : (link_type_covs == "square")    ? std::function<double(const double&)>(square_double)
      : std::function<double(const double&)>(exp_double);
  
  auto dist_fun  = [&](const Eigen::RowVectorXd& Pt,
                       const Eigen::MatrixXd& A,
                       const Eigen::MatrixXd& U,
                       const Eigen::MatrixXd& U_inv,
                       const Eigen::MatrixXd& D,
                       double t) -> Eigen::RowVectorXd {
                         if (append) return trans_dist_append_states(m, Pt, k, A, t);
                         if (to_lowercase(transient_dist_method) == "uniformization") return transient_dist_Uni(m, Pt, A, eps, U, U_inv, D, t);
                         if (to_lowercase(transient_dist_method) == "eigenvalue_decomp") return transient_dist_Eig(m, Pt, A, eps, U, U_inv, D, t);
                         return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, t);
                       };
  
  auto build_single_AUD = [&](const Eigen::VectorXd& par_base,
                              Eigen::MatrixXd& A,
                              Eigen::MatrixXd& U,
                              Eigen::MatrixXd& U_inv,
                              Eigen::MatrixXd& D) {
    A.setZero(m, m); U.setZero(m, m); U_inv.setZero(m, m); D.setZero(m, m);
    const std::string g = to_lowercase(generator);
    
    if (g == "gerlang") {
      Eigen::VectorXd lam = link_base(par_base.array()).matrix();
      A = make_A1(m, lam);
      U = eigenspace_U(m, lam);
      U_inv = eigenspace_U_inv(m, lam);
      for (int i = 0; i < m - 1; ++i) D(i, i) = -lam(i);
      return;
    }
    
    auto eig_pack = [&](const Eigen::MatrixXd& A_in) {
      A = A_in;
      Eigen::EigenSolver<Eigen::MatrixXd> es(A);
      Eigen::MatrixXcd Dc = es.eigenvalues().asDiagonal();
      Eigen::MatrixXcd Uc = es.eigenvectors();
      D = Dc.real();
      U = Uc.real();
      U_inv = U.inverse();
    };
    
    if (g == "gerlang_relax") { eig_pack(make_A2(m, link_base(par_base.array()).matrix())); return; }
    if (g == "free_upper_tri") { eig_pack(make_A3(m, link_base(par_base.array()).matrix())); return; }
    if (g == "bidiagonal")     { eig_pack(make_A4(m, link_base(par_base.array()).matrix())); return; }
    if (g == "tridiagonal")    { eig_pack(make_A5(m, link_base(par_base.array()).matrix())); return; }
    
    Rcpp::stop("Unknown generator.");
  };
  
  if (mixture) {
    if (K < 2) Rcpp::stop("mixture=TRUE requires K >= 2");
    if (warping) Rcpp::stop("warping=TRUE not implemented for mixture branch in MJP_predict yet.");
    
    const int p = (int)z.cols();
    const int base_len = base_len_from_generator(m, generator);
    if (base_len < 0) Rcpp::stop("Unknown generator in mixture branch.");
    
    const int mix_len  = K - 1;
    const int beta_len = covs_bin ? (state_covs ? (m - 1) * p : p) : 0;
    const int needed   = mix_len + K * base_len + beta_len;
    if ((int)pars.size() != needed) {
      Rcpp::stop("pars length mismatch for mixture. Got %d, expected %d", (int)pars.size(), needed);
    }
    
    Eigen::VectorXd eta(K); eta.setZero();
    for (int kk = 0; kk < K - 1; ++kk) eta(kk) = pars(kk);
    eta(K - 1) = 0.0;
    Eigen::ArrayXd wexp = (eta.array() - eta.maxCoeff()).exp();
    Eigen::VectorXd pi = (wexp / wexp.sum()).matrix();
    
    Eigen::VectorXd beta_covs;
    if (beta_len > 0) beta_covs = pars.tail(beta_len);
    
    std::vector<Eigen::MatrixXd> As; As.reserve(K);
    int off = mix_len;
    for (int kk = 0; kk < K; ++kk) {
      Eigen::VectorXd th = pars.segment(off + kk * base_len, base_len);
      Eigen::VectorXd lam = link_base(th.array()).matrix();
      std::string g = to_lowercase(generator);
      if      (g == "gerlang")        As.push_back(make_A1(m, lam));
      else if (g == "gerlang_relax")  As.push_back(make_A2(m, lam));
      else if (g == "free_upper_tri") As.push_back(make_A3(m, lam));
      else if (g == "bidiagonal")     As.push_back(make_A4(m, lam));
      else if (g == "tridiagonal")    As.push_back(make_A5(m, lam));
    }
    
    Eigen::MatrixXd Ud = Eigen::MatrixXd::Zero(m, m), Uid = Eigen::MatrixXd::Zero(m, m), Dd = Eigen::MatrixXd::Zero(m, m);
    Eigen::MatrixXd sol(n, m); sol.setZero();
    Eigen::RowVectorXd Pt(m);
    
    for (int i = 0; i < n; ++i) {
      int s = (int)(s1(i) - 1);
      Pt.setZero(); Pt(s) = 1.0;
      
      double t_eff = (double)u(i);
      
      if (covs_bin) {
        if (!state_covs) {
          t_eff = link_cov(beta_covs.dot(z.row(i))) * (double)u(i);
        } else {
          Eigen::RowVectorXd mu_bar = Eigen::RowVectorXd::Zero(m);
          for (int kk = 0; kk < K; ++kk) mu_bar += pi(kk) * expected_sojourn_exact(m, u(i), Pt, As[kk]);
          double denom = mu_bar.sum(); if (denom <= 0) denom = 1.0;
          Eigen::RowVectorXd w = mu_bar / denom;
          double lp = 0.0;
          for (int r = 0; r < m - 1; ++r) lp += w(r) * beta_covs.segment(r * p, p).dot(z.row(i));
          t_eff = link_cov(lp) * (double)u(i);
        }
      }
      
      Eigen::RowVectorXd pred = Eigen::RowVectorXd::Zero(m);
      for (int kk = 0; kk < K; ++kk) pred += pi(kk) * dist_fun(Pt, As[kk], Ud, Uid, Dd, t_eff);
      sol.row(i) = pred;
    }
    return sol;
  }
  
  Eigen::MatrixXd A(m, m), U(m, m), U_inv(m, m), D(m, m);
  {
    const int base_len = base_len_from_generator(m, generator);
    Eigen::VectorXd par_base = pars.segment(0, base_len);
    build_single_AUD(par_base, A, U, U_inv, D);
  }
  
  Eigen::MatrixXd sol(n, m); sol.setZero();
  Eigen::RowVectorXd Pt(m);
  
  auto warp_time = [&](int i, int start_idx) -> double {
    Eigen::VectorXd xii = pars.segment(base_len_from_generator(m, generator), m - 1);
    Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);
    xi.head(m - 1) = exp_vec(xii);
    double meanlog = xi.head(m - 1).array().log().mean();
    xi.head(m - 1) *= std::exp(-meanlog);
    xi(m - 1) = 1.0;
    Pt.setZero(); Pt(start_idx) = 1.0;
    Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
    double tau = 0.0;
    for (int j = 0; j < m; ++j) tau += mu(j) * xi(j);
    return tau;
  };
  
  if (!covs_bin) {
    for (int i = 0; i < n; ++i) {
      int s = (int)(s1(i) - 1);
      Pt.setZero(); Pt(s) = 1.0;
      double t_eff = warping ? warp_time(i, s) : (double)u(i);
      sol.row(i) = dist_fun(Pt, A, U, U_inv, D, t_eff);
    }
    return sol;
  }
  
  const int p = (int)z.cols();
  const int beta_len = state_covs ? (m - 1) * p : p;
  Eigen::VectorXd beta_covs = pars.tail(beta_len);
  
  for (int i = 0; i < n; ++i) {
    int s = (int)(s1(i) - 1);
    Pt.setZero(); Pt(s) = 1.0;
    
    double base_t = warping ? warp_time(i, s) : (double)u(i);
    
    double lp = 0.0;
    if (!state_covs) {
      lp = beta_covs.dot(z.row(i));
    } else {
      Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
      Eigen::RowVectorXd w = mu / mu.sum();
      for (int r = 0; r < m - 1; ++r) lp += w(r) * beta_covs.segment(r * p, p).dot(z.row(i));
    }
    
    double t_eff = link_cov(lp) * (warping ? base_t : (double)u(i));
    if (warping) t_eff = link_cov(lp) * base_t;
    
    sol.row(i) = dist_fun(Pt, A, U, U_inv, D, t_eff);
  }
  
  return sol;
}

// [[Rcpp::export]]
Rcpp::List MJP_effective_time(int m,
                                    const Eigen::VectorXd& s1,
                                    const Eigen::VectorXd& u,
                                    const Eigen::VectorXd& pars,
                                    const Eigen::MatrixXd& z,
                                    const std::string& generator = "gerlang",
                                    const std::string& link_type_base = "exp",
                                    const std::string& link_type_covs = "exp",
                                    bool covs_bin = true,
                                    const std::string& transient_dist_method = "pade",
                                    double eps = 2.220446049250313e-16,
                                    bool warping = false,
                                    bool state_covs = true,
                                    int k = 1, bool append = false,
                                    bool mixture = false, int K = 1)
{
  const int n = (int)u.size();
  Eigen::VectorXd t_eff_vec(n);
  Eigen::MatrixXd W(n, m);
  W.setZero();
  
  auto link_base = (link_type_base == "softplus") ? std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)>(soft_plus_vec)
    : (link_type_base == "square")   ? std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)>(square_vec)
      : std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)>(exp_vec);
  
  auto link_cov  = (link_type_covs == "softplus") ? std::function<double(const double&)>(soft_plus_double)
    : (link_type_covs == "square")    ? std::function<double(const double&)>(square_double)
      : std::function<double(const double&)>(exp_double);
  
  auto dist_fun  = [&](const Eigen::RowVectorXd& Pt,
                       const Eigen::MatrixXd& A,
                       const Eigen::MatrixXd& U,
                       const Eigen::MatrixXd& U_inv,
                       const Eigen::MatrixXd& D,
                       double t) -> Eigen::RowVectorXd {
                         if (append) return trans_dist_append_states(m, Pt, k, A, t);
                         if (to_lowercase(transient_dist_method) == "uniformization") return transient_dist_Uni(m, Pt, A, eps, U, U_inv, D, t);
                         if (to_lowercase(transient_dist_method) == "eigenvalue_decomp") return transient_dist_Eig(m, Pt, A, eps, U, U_inv, D, t);
                         return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, t);
                       };
  
  auto build_single_AUD = [&](const Eigen::VectorXd& par_base,
                              Eigen::MatrixXd& A,
                              Eigen::MatrixXd& U,
                              Eigen::MatrixXd& U_inv,
                              Eigen::MatrixXd& D) {
    A.setZero(m, m); U.setZero(m, m); U_inv.setZero(m, m); D.setZero(m, m);
    const std::string g = to_lowercase(generator);
    
    if (g == "gerlang") {
      Eigen::VectorXd lam = link_base(par_base.array()).matrix();
      A = make_A1(m, lam);
      U = eigenspace_U(m, lam);
      U_inv = eigenspace_U_inv(m, lam);
      for (int i = 0; i < m - 1; ++i) D(i, i) = -lam(i);
      return;
    }
    
    auto eig_pack = [&](const Eigen::MatrixXd& A_in) {
      A = A_in;
      Eigen::EigenSolver<Eigen::MatrixXd> es(A);
      Eigen::MatrixXcd Dc = es.eigenvalues().asDiagonal();
      Eigen::MatrixXcd Uc = es.eigenvectors();
      D = Dc.real();
      U = Uc.real();
      U_inv = U.inverse();
    };
    
    if (g == "gerlang_relax")  { eig_pack(make_A2(m, link_base(par_base.array()).matrix())); return; }
    if (g == "free_upper_tri") { eig_pack(make_A3(m, link_base(par_base.array()).matrix())); return; }
    if (g == "bidiagonal")     { eig_pack(make_A4(m, link_base(par_base.array()).matrix())); return; }
    if (g == "tridiagonal")    { eig_pack(make_A5(m, link_base(par_base.array()).matrix())); return; }
    
    Rcpp::stop("Unknown generator.");
  };
  
  if (mixture) {
    Rcpp::stop("mixture=TRUE not implemented in MJP_effective_time_and_w (copy your mixture branch if needed).");
  }
  
  Eigen::MatrixXd A(m, m), Umat(m, m), U_inv(m, m), Dmat(m, m);
  const int base_len = base_len_from_generator(m, generator);
  if (base_len < 0) Rcpp::stop("Unknown generator.");
  build_single_AUD(pars.segment(0, base_len), A, Umat, U_inv, Dmat);
  
  Eigen::RowVectorXd Pt(m);
  
  auto warp_time = [&](int i, int start_idx) -> double {
    Eigen::VectorXd xii = pars.segment(base_len, m - 1);
    Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);
    xi.head(m - 1) = exp_vec(xii);
    double meanlog = xi.head(m - 1).array().log().mean();
    xi.head(m - 1) *= std::exp(-meanlog);
    xi(m - 1) = 1.0;
    
    Pt.setZero(); Pt(start_idx) = 1.0;
    Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
    
    double tau = 0.0;
    for (int j = 0; j < m; ++j) tau += mu(j) * xi(j);
    return tau;
  };
  
  if (!covs_bin) {
    for (int i = 0; i < n; ++i) {
      int s = (int)(s1(i) - 1);
      Pt.setZero(); Pt(s) = 1.0;
      
      double t_eff = warping ? warp_time(i, s) : (double)u(i);
      t_eff_vec(i) = t_eff;
      
      Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
      double denom = mu.sum(); if (denom <= 0.0) denom = 1.0;
      W.row(i) = (mu / denom);
    }
    return Rcpp::List::create(
      Rcpp::Named("t_eff") = t_eff_vec,
      Rcpp::Named("w")     = W
    );
  }
  
  const int p = (int)z.cols();
  const int beta_len = state_covs ? (m - 1) * p : p;
  Eigen::VectorXd beta_covs = pars.tail(beta_len);
  
  for (int i = 0; i < n; ++i) {
    int s = (int)(s1(i) - 1);
    Pt.setZero(); Pt(s) = 1.0;
    
    double base_t = warping ? warp_time(i, s) : (double)u(i);
    
    double lp = 0.0;
    if (!state_covs) {
      lp = beta_covs.dot(z.row(i));
      // w not used in this setting; still return something reasonable:
      Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
      double denom = mu.sum(); if (denom <= 0.0) denom = 1.0;
      W.row(i) = (mu / denom);
    } else {
      Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
      double denom = mu.sum(); if (denom <= 0.0) denom = 1.0;
      Eigen::RowVectorXd w = mu / denom;
      W.row(i) = w;
      
      for (int r = 0; r < m - 1; ++r) {
        lp += w(r) * beta_covs.segment(r * p, p).dot(z.row(i));
      }
    }
    
    double t_eff = link_cov(lp) * (warping ? base_t : (double)u(i));
    if (warping) t_eff = link_cov(lp) * base_t;
    
    t_eff_vec(i) = t_eff;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("t_eff") = t_eff_vec,
    Rcpp::Named("w")     = W,
    Rcpp::Named("A")     = A
  );
}

// 
// // [[Rcpp::export]]
// Eigen::MatrixXd MJP_predict(int m,
//                             const Eigen::VectorXd& s1,
//                             const Eigen::VectorXd& u,
//                             const Eigen::VectorXd& pars,
//                             const Eigen::MatrixXd& z,
//                             const string& generator = "gerlang",
//                             const string& link_type_base = "exp",
//                             const string& link_type_covs = "exp",
//                             bool covs_bin = true,
//                             const string& transient_dist_method = "pade",
//                             double eps = 2.220446049250313e-16,
//                             bool warping = false, 
//                             bool state_covs = false,
//                             int k = 1, bool append = false, bool mixture = false, int K = 1){
// 
//   //Rcpp::Rcout << pars << std::endl;
//   //Eigen::VectorXd random_effects,
//   //Eigen::VectorXi groups,
// 
//   /* Initialization */
//   int n = u.size();
//   Eigen::VectorXd lambda_base;
//   Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m, m);
//   Eigen::MatrixXd U = Eigen::MatrixXd::Zero(m, m);
//   Eigen::MatrixXd U_inv = Eigen::MatrixXd::Zero(m, m);
//   Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m, m);
//   Eigen::VectorXd lambda(m);
//   Eigen::RowVectorXd Pt(m);
//   Eigen::RowVectorXd Ptu(m);
//   Eigen::MatrixXd solution(n,m);
//   double cov_time = 0;
//   int start_idx = 0;
//   bool eigen_solver_good = true;
// 
//   /* Determine which link to use */
//   std::function<Eigen::ArrayXd(const Eigen::ArrayXd&)> link_function_base;
//   if (link_type_base == "exp") {
//     link_function_base = exp_vec;
//   } else if (link_type_base == "softplus") {
//     link_function_base = soft_plus_vec;
//   } else if (link_type_base == "square") {
//     link_function_base = square_vec;
//   }
//   std::function<double(const double&)> link_function_covs;
//   if (link_type_covs == "exp") {
//     link_function_covs = exp_double;
//   } else if (link_type_covs == "softplus") {
//     link_function_covs = soft_plus_double;
//   } else if (link_type_covs == "square") {
//     link_function_covs = square_double;
//   }
//   
//   /* Determine how to calculate transient distribution */
//   std::function<Eigen::MatrixXd(int, const Eigen::RowVectorXd&, const Eigen::MatrixXd&, double, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, double)> transient_dist;
//   if (to_lowercase(transient_dist_method) == "uniformization") {
//     transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Uni(m, Pt, A, eps, U, U_inv, D, cov_time); };
//   } else if (to_lowercase(transient_dist_method) == "pade") {
//     transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
//   } else if (to_lowercase(transient_dist_method) == "eigenvalue_decomp" && eigen_solver_good){
//     //Rcpp::Rcout << "Using eigenvalue decomp." << std::endl;
//     transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Eig(m, Pt, A, eps, U, U_inv, D, cov_time); };
//   } else {
//     transient_dist = [](int m, const Eigen::RowVectorXd& Pt, const Eigen::MatrixXd& A, double eps, const Eigen::MatrixXd& U, const Eigen::MatrixXd& U_inv, const Eigen::MatrixXd& D, double cov_time) { return transient_dist_Pade(m, Pt, A, eps, U, U_inv, D, cov_time); };
//   }
//   
//   // ---------------------------
//   // Mixture-of-MJPs branch
//   // ---------------------------
//   if (mixture) {
//     if (K < 2) Rcpp::stop("mixture=TRUE requires K >= 2");
//     
//     const int n = u.size();
//     const int p = z.cols();
//     Eigen::MatrixXd solution(n, m);
//     solution.setZero();
//     
//     // lengths
//     const int base_len = base_len_from_generator(m, generator);
//     if (base_len < 0) Rcpp::stop("Unknown generator in mixture branch.");
//     
//     const int mix_len  = K - 1;
//     const int beta_len = covs_bin ? (state_covs ? (m - 1) * p : p) : 0;
//     const int needed   = mix_len + K * base_len + beta_len;
//     
//     if (pars.size() != needed) {
//       Rcpp::stop("pars length mismatch for mixture. Got %d, expected %d",
//                  int(pars.size()), needed);
//     }
//     
//     // mixture weights pi via softmax of eta (last fixed to 0)
//     Eigen::VectorXd eta(K);
//     eta.setZero();
//     for (int kk = 0; kk < K - 1; ++kk) eta(kk) = pars(kk);
//     eta(K - 1) = 0.0;
//     
//     Eigen::ArrayXd wexp = (eta.array() - eta.maxCoeff()).exp();
//     Eigen::VectorXd pi  = (wexp / wexp.sum()).matrix();  // sums to 1
//     
//     // unpack beta (shared across components)
//     Eigen::VectorXd beta_covs;
//     if (beta_len > 0) beta_covs = pars.tail(beta_len);
//     
//     // build A_k
//     std::vector<Eigen::MatrixXd> As;
//     As.reserve(K);
//     
//     int offset = mix_len;
//     for (int kk = 0; kk < K; ++kk) {
//       Eigen::VectorXd theta_base = pars.segment(offset + kk * base_len, base_len);
//       Eigen::VectorXd lambda_k   = link_function_base(theta_base.array()).matrix();
//       
//       Eigen::MatrixXd Ak;
//       std::string g = to_lowercase(generator);
//       if (g == "gerlang")           Ak = make_A1(m, lambda_k);
//       else if (g == "gerlang_relax")Ak = make_A2(m, lambda_k);
//       else if (g == "free_upper_tri")Ak = make_A3(m, lambda_k);
//       else if (g == "bidiagonal")   Ak = make_A4(m, lambda_k);
//       else if (g == "tridiagonal")  Ak = make_A5(m, lambda_k);
//       else Rcpp::stop("Unknown generator in mixture branch.");
//       
//       As.push_back(Ak);
//     }
//     offset += K * base_len;
//     
//     // Choose transient_dist (reuse your existing selector).
//     // NOTE: transient_dist_Pade ignores U/U_inv/D internally, so we can pass dummies safely.
//     Eigen::MatrixXd Udummy = Eigen::MatrixXd::Zero(m, m);
//     Eigen::MatrixXd Uidummy = Eigen::MatrixXd::Zero(m, m);
//     Eigen::MatrixXd Ddummy = Eigen::MatrixXd::Zero(m, m);
//     
//     auto transient_dist_local = transient_dist; // uses your already-selected function object
//     
//     // Mixture prediction loop
//     Eigen::RowVectorXd Pt(m);
//     for (int i = 0; i < n; ++i) {
//       int start_idx = int(s1(i) - 1);
//       Pt.setZero();
//       Pt(start_idx) = 1.0;
//       
//       // effective time scaling
//       double t_eff = double(u(i));
//       
//       if (covs_bin) {
//         if (!state_covs) {
//           // global beta
//           double linpred = beta_covs.dot(z.row(i));
//           t_eff = link_function_covs(linpred) * double(u(i));
//         } else {
//           // state-weighted linpred using mu_bar = sum_k pi_k mu_k
//           Eigen::RowVectorXd mu_bar = Eigen::RowVectorXd::Zero(m);
//           for (int kk = 0; kk < K; ++kk) {
//             Eigen::RowVectorXd mu_k = expected_sojourn_exact(m, u(i), Pt, As[kk]);
//             mu_bar += pi(kk) * mu_k;
//           }
//           double denom = mu_bar.sum();
//           if (denom <= 0) denom = 1.0; // safety
//           
//           Eigen::RowVectorXd w = mu_bar / denom;
//           
//           double linpred = 0.0;
//           for (int r = 0; r < m - 1; ++r) {
//             linpred += w(r) * beta_covs.segment(r * p, p).dot(z.row(i));
//           }
//           t_eff = link_function_covs(linpred) * double(u(i));
//         }
//       }
//       
//       // (Optional) warping not supported in this mixture branch as written
//       if (warping) {
//         Rcpp::stop("warping=TRUE not implemented for mixture branch in MJP_predict yet.");
//       }
//       
//       // mixture predictive distribution
//       Eigen::RowVectorXd pred = Eigen::RowVectorXd::Zero(m);
//       
//       for (int kk = 0; kk < K; ++kk) {
//         Eigen::RowVectorXd pk;
//         if (!append) {
//           pk = transient_dist_local(m, Pt, As[kk], eps, Udummy, Uidummy, Ddummy, t_eff);
//         } else {
//           pk = trans_dist_append_states(m, Pt, k, As[kk], t_eff);
//         }
//         pred += pi(kk) * pk;
//       }
//       
//       solution.row(i) = pred;
//     }
//     
//     return solution;
//   }
// 
//   /* Determine which parameterization to use */
//   std::function<Eigen::MatrixXd(int, const Eigen::VectorXd&)> make_A;
//   if (to_lowercase(generator) == "gerlang") {
//     make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A1(m, lambda); };
//     //lambda_base.resize(int(m-1));
//     lambda_base = link_function_base(pars.segment(0,m-1).array()); 
//     U = eigenspace_U(m, lambda_base);
//     U_inv = eigenspace_U_inv(m, lambda_base);
//     A = make_A(m, lambda_base);
//     for(int i = 0; i < (m-1); i++){D(i,i) = -lambda_base(i); }
//     bool distinct_rates = are_rates_distinct(m, lambda_base, 0.00000001); 
//     if (!distinct_rates) {eigen_solver_good = false;}
//   } else if (to_lowercase(generator) == "gerlang_relax") {
//     make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A2(m, lambda); };
//     //lambda_base.resize(int(m-1));
//     lambda_base = link_function_base(pars.segment(0,m-1).array()); 
//     A = make_A(m, lambda_base);
//     Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
//     Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
//     Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
//     D = D_complex.real();
//     U = U_complex.real(); 
//     U_inv = U.inverse();
//     if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
//   } else if (to_lowercase(generator) == "free_upper_tri") {
//     make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A3(m, lambda); };
//     //lambda_base.resize(int(m*(m-1)/2));
//     lambda_base = link_function_base(pars.segment(0,int(m*(m-1)/2)).array()); 
//     A = make_A(m, lambda_base);
//     Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
//     Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
//     Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
//     D = D_complex.real();
//     U = U_complex.real(); 
//     U_inv = U.inverse();
//     if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
//   } else if (to_lowercase(generator) == "bidiagonal") {
//     make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A4(m, lambda); };
//     //lambda_base.resize(int(2*m-3));
//     lambda_base = link_function_base(pars.segment(0,int(2*m-3)).array()); 
//     A = make_A(m, lambda_base);
//     Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
//     Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
//     Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
//     D = D_complex.real();
//     U = U_complex.real(); 
//     U_inv = U.inverse();
//     if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
//   } else if (to_lowercase(generator) == "tridiagonal") {
//     make_A = [](int m, const Eigen::VectorXd& lambda) { return make_A5(m, lambda); };
//     //lambda_base.resize(int(3*m-6));
//     lambda_base = link_function_base(pars.segment(0,int(3*m-6)).array()); 
//     A = make_A(m, lambda_base);
//     Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(A);
//     Eigen::MatrixXcd D_complex = eigensolver.eigenvalues().asDiagonal(); 
//     Eigen::MatrixXcd U_complex = eigensolver.eigenvectors(); 
//     D = D_complex.real();
//     U = U_complex.real(); 
//     U_inv = U.inverse();
//     if (eigensolver.info() != Eigen::Success) {eigen_solver_good = false;}
//   }
//   
//   //bool use_re = (random_effects.size() > 0);
// 
//   /* Compute score */
//   if (!covs_bin) { // Case: Covariates excluded
//     if(!warping){
//       for(int i = 0; i < n; i++){
//         cov_time = u(i);
//         //if(use_re){cov_time *= std::exp(random_effects(groups(i)));}
//         start_idx = int(s1(i)-1);
//         Pt.setZero();;
//         Pt(start_idx) = 1 ;
//         Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
//         if(append){ Ptu = trans_dist_append_states(m,Pt, k, A, cov_time); } //append states
//         solution.row(i) = Ptu;
//       } 
//     } else { // using warping
//       //Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
//       //Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);  // initialize with 1s
//       //xi.head(m - 1) = exp_vec(xii); //soft_plus_vec(xii);            // overwrite first m-1 entries
//       Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
//       Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);
//       xi.head(m - 1) = exp_vec(xii);
//       double meanlog = xi.head(m - 1).array().log().mean();
//       xi.head(m - 1) *= std::exp(-meanlog);
//       xi(m - 1) = 1.0;
//       for(int i = 0; i < n; i++){
//         start_idx = int(s1(i)-1);
//         Pt.setZero();
//         Pt(start_idx) = 1.0;
//         Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
//         double tau_eff = 0;
//         for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);
//         cov_time = tau_eff;
//         //if(use_re){cov_time *= std::exp(random_effects(groups(i)));}
//         Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
//         if(append){ Ptu = trans_dist_append_states(m,Pt, k, A, cov_time); } //append states
//         solution.row(i) = Ptu;
//       }
//     }
//   } else {  // Case: Covariates included
//     //int k = z.cols();
//     //Eigen::VectorXd beta_covs = pars.tail(pars.size() - lambda_base.size());
//     
//     int p = z.cols();
//     int beta_len = state_covs ? (m - 1) * p : p;
//     Eigen::VectorXd beta_covs = pars.tail(beta_len);
// 
//     if(!warping){
//       for(int i = 0; i < n; i++){
//         start_idx = int(s1(i)-1);
//         Pt.setZero();
//         Pt( start_idx) = 1;
// 
//         Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
//         // Eigen::MatrixXd Az = A;
//         // if(!state_covs){
//         //   double linpred = beta_covs.dot(z.row(i));
//         //   double g = link_function_covs(linpred);
//         //   Az.topRows(m - 1) *= g;
//         // } else {
//         //   for(int r = 0; r < m - 1; ++r){
//         //     double linpred_r = beta_covs.segment(r * p, p).dot(z.row(i));
//         //     double g_r = link_function_covs(linpred_r);
//         //     Az.row(r) *= g_r;
//         //   }
//         // }
//         // cov_time = u(i);
//         double linpred = 0.0;
//         if(!state_covs){
//           linpred = beta_covs.dot(z.row(i));
//         } else {
//           Eigen::RowVectorXd w = mu / mu.sum();
//           for(int r = 0; r < m - 1; ++r){
//             linpred += w(r) * beta_covs.segment(r * p, p).dot(z.row(i));
//           }
//         }
//         cov_time = link_function_covs(linpred) * u(i);
// 
//         //if(use_re){cov_time *= std::exp(random_effects(groups(i)));}
//         Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
//         if(append){ Ptu = trans_dist_append_states(m,Pt, k, A, cov_time); } //append states
//         solution.row(i) = Ptu;
//       }
//     // int p = z.cols();
//     // int beta_len = state_covs ? (m - 1) * p : p;
//     // Eigen::VectorXd beta_covs = pars.tail(beta_len);
//     // 
//     // if(!warping){
//     //   for(int i = 0; i < n; i++){
//     //     int start_idx = int(s1(i)-1);
//     //     Pt.setZero();
//     //     Pt(start_idx) = 1.0;
//     //     
//     //     Eigen::MatrixXd Az = A;
//     //     
//     //     if(!state_covs){
//     //       double linpred = beta_covs.dot(z.row(i));
//     //       double g = link_function_covs(linpred);
//     //       Az.topRows(m - 1) *= g;
//     //     } else {
//     //       for(int r = 0; r < m - 1; ++r){
//     //         double linpred_r = beta_covs.segment(r * p, p).dot(z.row(i));
//     //         double g_r = link_function_covs(linpred_r);
//     //         Az.row(r) *= g_r;
//     //       }
//     //     }
//     //     
//     //     Ptu = transient_dist(m, Pt, Az, eps, U, U_inv, D, u(i));
//     //     if(append){ Ptu = trans_dist_append_states(m, Pt, k, Az, u(i)); }
//     //     solution.row(i) = Ptu;
//     //   }
//     } else { // using warping
//     //   Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
//     //   Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);
//     //   xi.head(m - 1) = exp_vec(xii);
//     //   double meanlog = xi.head(m - 1).array().log().mean();
//     //   xi.head(m - 1) *= std::exp(-meanlog);
//     //   xi(m - 1) = 1.0;
//     //   
//     //   for(int i = 0; i < n; i++){
//     //     start_idx = int(s1(i)-1);
//     //     Pt.setZero();
//     //     Pt(start_idx) = 1.0;
//     //     
//     //     Eigen::MatrixXd Az = A;
//     //     if(!state_covs){
//     //       double linpred = beta_covs.dot(z.row(i));
//     //       double g = link_function_covs(linpred);
//     //       Az.topRows(m - 1) *= g;
//     //     } else {
//     //       for(int r = 0; r < m - 1; ++r){
//     //         double linpred_r = beta_covs.segment(r * p, p).dot(z.row(i));
//     //         double g_r = link_function_covs(linpred_r);
//     //         Az.row(r) *= g_r;
//     //       }
//     //     }
//     //     
//     //     Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
//     //     
//     //     double tau_eff = 0.0;
//     //     for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);
//     //     
//     //     Ptu = transient_dist(m, Pt, Az, eps, U, U_inv, D, tau_eff);
//     //     if(append){ Ptu = trans_dist_append_states(m, Pt, k, Az, tau_eff); }
//     //     solution.row(i) = Ptu;
//     //   }
//     // }
//       //Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
//       //Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);  // initialize with 1s
//       //xi.head(m - 1) = exp_vec(xii); //soft_plus_vec(xii);            // overwrite first m-1 entries
//       Eigen::VectorXd xii = pars.segment(lambda_base.size(), m - 1);
//       Eigen::VectorXd xi = Eigen::VectorXd::Ones(m);
//       xi.head(m - 1) = exp_vec(xii);
//       double meanlog = xi.head(m - 1).array().log().mean();
//       xi.head(m - 1) *= std::exp(-meanlog);
//       xi(m - 1) = 1.0;
//       for(int i = 0; i < n; i++){
//         start_idx = int(s1(i)-1);
//         Pt.setZero();
//         Pt(start_idx) = 1.0;
// 
//         // Eigen::MatrixXd Az = A;
//         // if(!state_covs){
//         //   double linpred = beta_covs.dot(z.row(i));
//         //   double g = link_function_covs(linpred);
//         //   Az.topRows(m - 1) *= g;
//         // } else {
//         //   for(int r = 0; r < m - 1; ++r){
//         //     double linpred_r = beta_covs.segment(r * p, p).dot(z.row(i));
//         //     double g_r = link_function_covs(linpred_r);
//         //     Az.row(r) *= g_r;
//         //   }
//         // }
//         
// 
//         Eigen::RowVectorXd mu = expected_sojourn_exact(m, u(i), Pt, A);
//         double tau_eff = 0;
//         for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);
// 
//         //cov_time *= tau_eff;
//         
//         double linpred = 0.0;
//         if(!state_covs){
//           linpred = beta_covs.dot(z.row(i));
//         } else {
//           Eigen::RowVectorXd w = mu / mu.sum(); //w = mu / mu.sum();
//           for(int r = 0; r < m - 1; ++r){
//             linpred += w(r) * beta_covs.segment(r * p, p).dot(z.row(i));
//           }
//         }
//         cov_time = link_function_covs(linpred) * tau_eff;
// 
// 
//         // double linpred = 0.0;
//         // if(state_covs){
//         //   if(start_idx < m - 1){linpred = beta_covs.segment(start_idx * p, p).dot(z.row(i));} else {linpred = 0.0;}
//         // } else { linpred = beta_covs.dot(z.row(i));}
//         // cov_time = link_function_covs(linpred) * tau_eff
// 
// 
//         //if(use_re){cov_time *= std::exp(random_effects(groups(i)));}
//         Ptu = transient_dist(m, Pt, A, eps, U, U_inv, D, cov_time);
//         if(append){ Ptu = trans_dist_append_states(m,Pt, k, A, cov_time); } //append states
//         solution.row(i) = Ptu;
//       }
//     }
//   }
//   return solution;
// }



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



// ---------- reachability on the graph induced by Q ----------
static void reach(const Eigen::MatrixXd& Q,
                  int s,                      // 0-based
                  std::vector<int>& ok,       // output mask (0/1)
                  bool forward = true)        // true: edges a->b if q_ab>0; false: reverse
{
  const int m = Q.rows();
  ok.assign(m, 0);
  ok[s] = 1;
  bool changed = true;
  while (changed) {
    changed = false;
    for (int a = 0; a < m; ++a) if (ok[a]) {
      for (int b = 0; b < m; ++b) {
        if (a == b) continue;
        const bool edge = forward ? (Q(a,b) > 0.0) : (Q(b,a) > 0.0);
        if (edge && !ok[b]) { ok[b] = 1; changed = true; }
      }
    }
  }
}

// ---------- Van Loan + pruning (double) ----------
// [[Rcpp::export]]
static Eigen::VectorXd expected_sojourn_bridge_impl(const Eigen::MatrixXd& A,
                                                    double u,
                                                    int s1,  // 0-based
                                                    int s2)  // 0-based
{
  const int m = A.rows();
  Eigen::VectorXd mu = Eigen::VectorXd::Zero(m);
  if (u <= 0.0) return mu;
  
  // Denominator: P_ij(u) = [exp(A u)]_{s1,s2}
  Eigen::MatrixXd P = (A * u).exp();
  const double denom = P(s1, s2);
  if (!(denom > 1e-14)) {
    // numerically impossible bridge → zeros
    return mu;
  }
  
  // Prune states: must be reachable from s1 and must reach s2
  std::vector<int> from_s1(m), to_s2(m);
  reach(A, s1, from_s1, /*forward=*/true);
  reach(A, s2, to_s2,   /*forward=*/false);
  
  // Prebuild block matrix B = [A  0; 0  A]; then set TR = Δ_k per k
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(2*m, 2*m);
  B.block(0,   0,   m, m) = A;
  B.block(m,   m,   m, m) = A;
  
  for (int k = 0; k < m; ++k) {
    if (!(from_s1[k] && to_s2[k])) {
      mu(k) = 0.0;
      continue;
    }
    // Top-right block = Δ_k (all zeros except (k,k)=1)
    B.block(0, m, m, m).setZero();
    B(k, m + k) = 1.0;
    
    // Van Loan block exponential
    Eigen::MatrixXd E = (B * u).exp();
    
    // Upper-right block entry (s1, s2) gives ∫ e^{A t} Δ_k e^{A(u-t)} dt, entry (s1,s2)
    mu(k) = E(s1, m + s2) / denom;
  }
  return mu;
}

Rcpp::NumericVector expected_sojourn_bridge_rcpp(const Eigen::Map<Eigen::MatrixXd>& Q,
                                                 double u,
                                                 int s1,
                                                 int s2,
                                                 bool one_based = true){
  if (Q.rows() != Q.cols())
    Rcpp::stop("Q must be square");
  int m = Q.rows();
  int i = one_based ? (s1 - 1) : s1;
  int j = one_based ? (s2 - 1) : s2;
  if (i < 0 || i >= m || j < 0 || j >= m)
    Rcpp::stop("s1/s2 out of range after indexing adjustment");
  
  Eigen::VectorXd mu = expected_sojourn_bridge_impl(Q, u, i, j);
  return Rcpp::wrap(mu);
}

