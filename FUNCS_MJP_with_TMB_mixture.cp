#include <TMB.hpp>


template<class Type>
vector<Type> softplus(const vector<Type>& v) {
  return exp(v); 
  //return log(exp(v) + 1);
}



// Generalized Erlang (make_A1)
template<class Type>
matrix<Type> make_A1(int m, const vector<Type> &lambda) {
  matrix<Type> A(m, m);
  A.setZero();
  for (int i = 0; i < m - 1; ++i) {A(i, i + 1) = lambda(i);  }
  for (int i = 0; i < m; ++i) {A(i, i) = -A.row(i).sum(); }
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
  for (int i = 0; i < m; ++i) A(i,i) = Type(0);
  for (int i = 0; i < m; ++i) { A(i, i) = -A.row(i).sum();}
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
  for (int i = 0; i < m; ++i) {A(i, i) = -A.row(i).sum(); }
  return A;
}

// Upper bidiagonal (make_A4)
template<class Type>
matrix<Type> make_A4(int m, const vector<Type> &lambda) {
  matrix<Type> A(m, m);
  A.setZero();
  for (int i = 0; i < m - 1; ++i) {A(i, i + 1) = lambda(i);}
  for (int i = 0; i < m-2; ++i) {A(i, i + 2) = lambda(m-1+i);}
  for (int i = 0; i < m; ++i) {A(i, i) = -A.row(i).sum();}
  return A;
}

// Upper bidiagonal (make_A5)
template<class Type>
matrix<Type> make_A5(int m, const vector<Type> &lambda) {
  matrix<Type> A(m, m);
  A.setZero();
  for (int i = 0; i < m - 1; ++i) {A(i, i + 1) = lambda(i);}
  for (int i = 0; i < m-2; ++i) {A(i, i + 2) = lambda(m-1+i);}
  for (int i = 0; i < m-3; ++i) {A(i, i + 3) = lambda(2*m-3+i);}
  for (int i = 0; i < m; ++i) {A(i, i) = -A.row(i).sum();}
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
vector<Type> expected_sojourn_exact(int m, Type u, int s1, matrix<Type> A, Type dt = Type(0.005)) {
  vector<Type> mu(m);
  mu.setZero();
  if (s1 >= m - 1) {
    mu(m - 1) = u;
    return mu;
  }
  matrix<Type> S = A.block(0, 0, m - 1, m - 1);
  matrix<Type> expSu = atomic::expm(matrix<Type>(S * u));
  matrix<Type> I(m - 1, m - 1);
  I.setZero();
  for (int i = 0; i < m - 1; ++i) I(i, i) = Type(1);
  matrix<Type> rhs = I - expSu;
  matrix<Type> Minv = matrix<Type>(-S).inverse();
  matrix<Type> F = Minv * rhs;
  Type sum_trans = Type(0);
  for (int j = 0; j < m - 1; ++j) {
    mu(j) = F(s1, j);
    sum_trans += mu(j);
  }
  mu(m - 1) = u - sum_trans;
  
  return mu;
}

template<class Type>
vector<Type> softmax_vec(const vector<Type>& eta) {
  Type m = eta.maxCoeff();
  vector<Type> ex = exp(eta.array() - m);
  Type s = ex.sum();
  return ex / s;
}

inline int base_len_from_generator(int m, int generator_type) {
  if (generator_type == 0 || generator_type == 1) return m - 1;           // A1/A2
  if (generator_type == 2) return m * (m - 1) / 2;                        // A3
  if (generator_type == 3) return 2 * m - 3;                              // A4
  if (generator_type == 4) return 3 * m - 6;                              // A5
  return -1;
}

template<class Type>
matrix<Type> make_A_by_type(int m, int generator_type, const vector<Type>& lambda) {
  if (generator_type == 0) return make_A1(m, lambda);
  if (generator_type == 1) return make_A2(m, lambda);
  if (generator_type == 2) return make_A3(m, lambda);
  if (generator_type == 3) return make_A4(m, lambda);
  if (generator_type == 4) return make_A5(m, lambda);
  error("invalid generator_type");
  matrix<Type> dummy(1,1); dummy.setZero();
  return dummy;
}
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_IVECTOR(s1);
  DATA_IVECTOR(s2);
  DATA_VECTOR(u);
  DATA_MATRIX(z);
  DATA_INTEGER(m);
  DATA_INTEGER(generator_type);
  DATA_INTEGER(cov_type);
  DATA_INTEGER(use_log_score);
  DATA_INTEGER(use_rps_score);
  DATA_INTEGER(use_brier_score);
  DATA_INTEGER(K);
  
  PARAMETER_VECTOR(theta); // theta = < <eta>, <lambda>_1, ... , <lambda>_k, <beta> >
  
  const int n = s1.size();
  const int p = z.cols();
  Type total_score = 0.0;
  
  int base_len = base_len_from_generator(m, generator_type);
  if (base_len < 0) error("invalid generator");
  int mix_len = K - 1;
  int cov_len = 0;
  if (cov_type == 0) cov_len = 0;
  else if (cov_type == 1) cov_len = p;
  else if (cov_type == 2) cov_len = (m - 1) * p;
  else error("invalid cov");
  
  int needed = mix_len + K * base_len + cov_len;
  if ((int)theta.size() != needed) {error("theta has wrong length", (int)theta.size(), needed);}
  
  vector<Type> eta(K);
  eta.setZero();
  for (int k = 0; k < K - 1; ++k) eta(k) = theta(k);
  eta(K - 1) = Type(0);
  vector<Type> pi = softmax_vec(eta);
  
  // store A_k in a vector of matrices
  std::vector< matrix<Type> > As;
  As.reserve(K);
  
  int offset = mix_len;
  for (int k = 0; k < K; ++k) {
    vector<Type> theta_base = theta.segment(offset + k * base_len, base_len);
    vector<Type> lambda = softplus(theta_base);
    matrix<Type> A = make_A_by_type(m, generator_type, lambda);
    As.push_back(A);
  }
  offset += K * base_len;
  
  vector<Type> theta_cov;
  if (cov_len > 0) theta_cov = theta.segment(offset, cov_len);
  
  for (int i = 0; i < n; ++i) {
    vector<Type> obs(m); obs.setZero();
    int start = s1(i) - 1;
    int end   = s2(i) - 1;
    obs(end)  = Type(1);
    
    Type t_eff = u(i);
    if (cov_type == 1) {
      Type lp = 0.0;
      for (int j = 0; j < p; ++j) lp += z(i,j) * theta_cov(j);
      t_eff = exp(lp) * u(i);
    }
    
    if (cov_type == 2) {
      vector<Type> mu_bar(m); mu_bar.setZero();
      
      for (int k = 0; k < K; ++k) {
        vector<Type> mu_k = expected_sojourn_exact(m, u(i), start, As[k]);
        mu_bar += pi(k) * mu_k;
      }
      
      vector<Type> w = mu_bar;
      w /= w.sum();
      
      Type lp = 0.0;
      for (int ii = 0; ii < m - 1; ++ii) {
        Type lp_ii = 0.0;
        for (int j = 0; j < p; ++j) lp_ii += z(i,j) * theta_cov(ii * p + j);
        lp += w(ii) * lp_ii;
      }
      t_eff = exp(lp) * u(i);
    }
    
    vector<Type> pred(m); pred.setZero();
    
    for (int k = 0; k < K; ++k) {
      matrix<Type> tpm = atomic::expm(matrix<Type>(As[k] * t_eff));
      vector<Type> pred_k = tpm.row(start).transpose();
      pred += pi(k) * pred_k;
    }
    
    if (use_log_score)   total_score += log_score(pred, obs);
    if (use_brier_score) total_score += brier_score(pred, obs);
    if (use_rps_score)   total_score += rps_score(pred, obs);
  }
  
  SIMULATE {
    // pred_mat: n x m matrix of predictive distributions
    matrix<Type> pred_mat(n, m);
    pred_mat.setZero();
    
    for (int i = 0; i < n; ++i) {
      int start = s1(i) - 1;
      
      Type t_eff = u(i);
      
      if (cov_type == 1) {
        Type lp = Type(0);
        for (int j = 0; j < p; ++j) lp += z(i, j) * theta_cov(j);
        t_eff = exp(lp) * u(i);
      }
      
      if (cov_type == 2) {
        vector<Type> mu_bar(m); mu_bar.setZero();
        for (int kk = 0; kk < K; ++kk) {
          vector<Type> mu_k = expected_sojourn_exact(m, u(i), start, As[kk]);
          for (int j = 0; j < m; ++j) mu_bar(j) += pi(kk) * mu_k(j);
        }
        Type denom = Type(0);
        for (int j = 0; j < m; ++j) denom += mu_bar(j);
        vector<Type> w(m);
        for (int j = 0; j < m; ++j) w(j) = mu_bar(j) / denom;
        
        Type lp = Type(0);
        for (int ii = 0; ii < m - 1; ++ii) {
          Type lp_ii = Type(0);
          for (int j = 0; j < p; ++j) lp_ii += z(i, j) * theta_cov(ii * p + j);
          lp += w(ii) * lp_ii;
        }
        t_eff = exp(lp) * u(i);
      }
      
      vector<Type> pred(m); pred.setZero();
      
      for (int kk = 0; kk < K; ++kk) {
        matrix<Type> tpm = atomic::expm(matrix<Type>(As[kk] * t_eff));
        
        // pred_k = tpm.row(start).transpose()
        for (int j = 0; j < m; ++j) {
          pred(j) += pi(kk) * tpm(start, j);
        }
      }
      
      // store row i
      for (int j = 0; j < m; ++j) pred_mat(i, j) = pred(j);
    }
    
    REPORT(pred_mat);
  }
  
  return total_score;
}

