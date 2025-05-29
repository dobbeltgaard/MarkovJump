#include <TMB.hpp>
using namespace atomic;

// template<class Type>
// vector<Type> softplus(const vector<Type> &x) {
//   vector<Type> res(x.size());
//   for (int i = 0; i < x.size(); ++i) {
//     res(i) = log(exp(x(i)) + 1); 
//   }
//   return res;
// }

template<class Type>
vector<Type> softplus(const vector<Type>& v) {
  return log(exp(v) + 1);
}



// Generalized Erlang (make_A1)
template<class Type>
matrix<Type> make_A1(int m, const vector<Type> &lambda) {
  matrix<Type> A(m, m);
  A.setZero();
  for (int i = 0; i < m - 1; ++i) {
    A(i, i + 1) = lambda(i);
  }
  for (int i = 0; i < m; ++i) {
    A(i, i) = -A.row(i).sum();
  }
  return A;
}


// Relaxed Erlang (make_A2)
template<class Type>
matrix<Type> make_A2(int m, const vector<Type> &lambda) {
  matrix<Type> A = make_A1(m, lambda);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      if (j - i > 1) {
        A(i, j) = (A(i, j - 1) * A(j - 1, j)) / (A(i, j - 1) + A(j - 1, j));
      }
    }
  }
  for (int i = 0; i < m; ++i) {
    A(i, i) = -A.row(i).sum();
  }
  return A;
}


// Free upper triangular (make_A3)
template<class Type>
matrix<Type> make_A3(int m, const vector<Type> &lambda) {
  matrix<Type> A(m, m);
  A.setZero();
  int count = 0;
  for (int i = 0; i < m; ++i) {
    for (int j = i + 1; j < m; ++j) {
      A(i, j) = lambda(count++);
    }
  }
  for (int i = 0; i < m; ++i) {
    A(i, i) = -A.row(i).sum();
  }
  return A;
}


// Log-score
template<class Type>
Type log_score(const vector<Type> &pred, const vector<Type> &obs) {
  for (int i = 0; i < obs.size(); ++i) {
    if (obs(i) == Type(1.0)) return -log(pred(i));
  }
  return Type(0.0);
}

// Brier score
template<class Type>
Type brier_score(const vector<Type> &pred, const vector<Type> &obs) {
  vector<Type> diff = pred - obs;
  return (diff * diff).sum();
}

// Ranked Probability Score
template<class Type>
Type rps_score(const vector<Type> &pred, const vector<Type> &obs) {
  Type res = 0.0;
  Type cum_pred = 0.0;
  Type cum_obs = 0.0;
  for (int i = 0; i < pred.size(); ++i) {
    cum_pred += pred(i);
    cum_obs += obs(i);
    Type diff = cum_pred - cum_obs;
    res += diff * diff;
  }
  return res;
}


template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(s1);
  DATA_VECTOR(s2);
  DATA_VECTOR(u);
  DATA_MATRIX(z);
  DATA_INTEGER(m);
  DATA_INTEGER(generator_type); // 0=A1, 1=A2, 2=A3
  DATA_INTEGER(cov_type);       // 0=no covs, 1=covs
  PARAMETER_VECTOR(theta);
  int n = s1.size();
  
  Type total_score = 0.0;
  matrix<Type> A(m, m);
  if (generator_type == 0) {
    vector<Type> theta_base = theta.segment(0, m - 1); 
    vector<Type> lambda = softplus(theta_base);
    A = make_A1(m, lambda);
  } else if (generator_type == 1) {
    vector<Type> theta_base = theta.segment(0, m - 1); 
    vector<Type> lambda = softplus(theta_base);
    A = make_A2(m, lambda);
  } else if (generator_type == 2) {
    vector<Type> theta_base = theta.segment(0, m * (m - 1) / 2); 
    vector<Type> lambda = softplus(theta_base);
    A = make_A3(m, lambda);
  } else {
    error("Invalid generator_type");
  }
  
  if (cov_type == 0) {
    for (int i = 0; i < n; ++i) {
      vector<Type> obs(m); obs.setZero();
      int start = CppAD::Integer(s1(i)) - 1;
      int end   = CppAD::Integer(s2(i)) - 1;
      obs(end)  = Type(1.0);
      matrix<Type> tpm = atomic::expm( matrix<Type>(A*u(i)) ); 
      vector<Type> pred = tpm.row(start).transpose();
      total_score += log_score(pred, obs);
      total_score += brier_score(pred, obs);
      total_score += rps_score(pred, obs);
    }
  } else if (cov_type == 1) {
    int theta_cov_start = (generator_type == 2) ? m * (m - 1) / 2 : m - 1;
    vector<Type> theta_cov = theta.segment(theta_cov_start, theta.size() - theta_cov_start);
    for (int i = 0; i < n; ++i) {
      vector<Type> obs(m); obs.setZero();
      int start = CppAD::Integer(s1(i)) - 1;
      int end   = CppAD::Integer(s2(i)) - 1;
      Type cov_linpred = 0.0;
      for (int j = 0; j < z.cols(); ++j) {cov_linpred += z(i, j) * theta_cov(j);}
      Type cov_time = log(1 + exp(cov_linpred)) * u(i);
      obs(end)  = Type(1.0);
      matrix<Type> tpm = atomic::expm( matrix<Type>(A*cov_time) ); 
      vector<Type> pred = tpm.row(start).transpose();  
      total_score += log_score(pred, obs);
      total_score += brier_score(pred, obs);
      total_score += rps_score(pred, obs);
    }
  }
  return total_score/ Type(n);
}

