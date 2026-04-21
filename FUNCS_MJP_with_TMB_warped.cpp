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
  for (int i = 0; i < m; ++i) {
    A(i, i) = -A.row(i).sum();
  }
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

// Upper tridiagonal (make_A5)
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


// template<class Type>
// vector<Type> compute_conditional_sojourn(
//     const matrix<Type>& A,
//     int s1, int s2,
//     Type tau,
//     int K = 30)
// {
//   int m = A.rows();
//   Type delta = tau / Type(K);
//   matrix<Type> Ptau = atomic::expm(matrix<Type>(A * tau));
//   Type denom = Ptau(s1, s2);
//   vector<Type> mu(m);
//   mu.setZero();
// 
//   Type safe_denom = denom + Type(1e-12);
// 
//   for (int k = 0; k < K; ++k) {
//     Type t = delta * Type(k);
//     matrix<Type> At = atomic::expm(matrix<Type>(A * t));
//     matrix<Type> Atau_t = atomic::expm( matrix<Type> (A * (tau - t)));
// 
//     for (int j = 0; j < m; ++j) {
//       matrix<Type> Ej(m, m);
//       Ej.setZero();
//       Ej(j, j) = Type(1.0);
//       matrix<Type> W = At * Ej * Atau_t;
//       mu(j) += W(s1, s2) * delta;
//     }
//   }
// 
//   mu /= safe_denom;
//   return mu;
// }


template<class Type>
vector<Type> expected_sojourn(int m, Type u, int s1, matrix<Type> A, Type dt = Type(0.005)) {
  vector<Type> mu(m);
  mu.setZero();

  vector<Type> p(m);
  p.setZero();
  p(s1) = Type(1.0);
  matrix<Type> At = A.transpose();

  int K = CppAD::Integer(u / dt);
  for (int k = 0; k < K; ++k) {
    mu += p;
    p += dt * (At*p);
  }
  mu *= dt;
  return mu;
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
Type objective_function<Type>::operator() () {
  DATA_IVECTOR(s1);
  DATA_IVECTOR(s2);
  DATA_VECTOR(u);
  DATA_MATRIX(z);
  DATA_INTEGER(m);
  DATA_INTEGER(generator_type); // 0=A1, 1=A2, 2=A3
  DATA_INTEGER(cov_type);       // 0=no covs, 1=covs
  DATA_INTEGER(use_log_score);
  DATA_INTEGER(use_rps_score);
  DATA_INTEGER(use_brier_score);
  PARAMETER_VECTOR(theta);
  int n = s1.size();
  int p = z.cols();
  
  Type total_score = 0.0;
  matrix<Type> A(m, m);
  int lambda_len;
  if (generator_type == 0) {
    lambda_len = m-1;
    vector<Type> theta_base = theta.segment(0, lambda_len); 
    vector<Type> lambda = softplus(theta_base);
    A = make_A1(m, lambda);
  } else if (generator_type == 1) {
    lambda_len = m-1;
    vector<Type> theta_base = theta.segment(0, lambda_len);
    vector<Type> lambda = softplus(theta_base);
    A = make_A2(m, lambda);
  } else if (generator_type == 2) {
    lambda_len = int(m * (m - 1) / 2);
    vector<Type> theta_base = theta.segment(0, lambda_len);
    vector<Type> lambda = softplus(theta_base);
    A = make_A3(m, lambda);
  } else if (generator_type == 3) {
    lambda_len = int(2*m-3);
    vector<Type> theta_base = theta.segment(0, lambda_len); 
    vector<Type> lambda = softplus(theta_base);
    A = make_A4(m, lambda);
  } else if (generator_type == 4) {
    lambda_len = int(3*m-6);
    vector<Type> theta_base = theta.segment(0, lambda_len);
    vector<Type> lambda = softplus(theta_base);
    A = make_A5(m, lambda);
  } else {
    error("Invalid generator_type");
  }
  
  vector<Type> xii = theta.segment(lambda_len, m - 1);
  vector<Type> xi(m);
  xi.setOnes();                           // initialize all to 1
  xi.head(m - 1) = softplus(xii);         // overwrite first m - 1 entries
  xi(m - 1) = Type(1.0);
  
  Type meanlog = 0.0;
  for (int j = 0; j < m - 1; ++j) meanlog += log(xi(j));
  meanlog /= Type(m - 1);
  xi.head(m - 1) = xi.head(m - 1) * exp(-meanlog);
  
  if (cov_type == 0) {
    for (int i = 0; i < n; ++i) {
      vector<Type> obs(m); obs.setZero();
      int start = s1(i) - int(1);
      int end   = s2(i) - int(1);
      obs(end)  = Type(1.0);
      
      vector<Type> mu = expected_sojourn_exact(m, u(i), start, A);
      Type tau_eff = 0.0; 
      for (int j = 0; j < xi.size(); ++j) {tau_eff += mu(j) * xi(j);} // dot product
      
      matrix<Type> tpm = atomic::expm( matrix<Type>(A*tau_eff) ); 
      vector<Type> pred = tpm.row(start).transpose();
      if (use_log_score) { total_score += log_score(pred, obs);}
      if (use_brier_score) {total_score += brier_score(pred, obs);}
      if (use_rps_score) {total_score += rps_score(pred, obs);}
    }
  } else if (cov_type == 1) {
    vector<Type> theta_cov = theta.segment(lambda_len + xii.size(), z.cols());
    for (int i = 0; i < n; ++i) {
      vector<Type> obs(m); obs.setZero();
      int start = s1(i) - int(1);
      int end   = s2(i) - int(1);
      Type cov_linpred = 0.0;
      for (int j = 0; j < z.cols(); ++j) {cov_linpred += z(i, j) * theta_cov(j);}
      obs(end)  = Type(1.0);
      vector<Type> mu = expected_sojourn_exact(m, u(i), start, A);
      Type tau_eff = 0.0; 
      for (int j = 0; j < xi.size(); ++j) {tau_eff += mu(j) * xi(j);} // dot product
      Type cov_scaling = exp(cov_linpred); //log(1 + exp(cov_linpred));  // softplus
      tau_eff *= cov_scaling;  // final time warp
      matrix<Type> tpm = atomic::expm( matrix<Type>(A*tau_eff) );
      vector<Type> pred = tpm.row(start).transpose();
      if (use_log_score) { total_score += log_score(pred, obs);}
      if (use_brier_score) {total_score += brier_score(pred, obs);}
      if (use_rps_score) {total_score += rps_score(pred, obs);}
    } 
    // } else if (cov_type == 2) {
    //   vector<Type> theta_cov = theta.segment(lambda_len + xii.size(), p * (m - 1));
    //   for (int i = 0; i < n; ++i) {
    //     vector<Type> obs(m); obs.setZero();
    //     int start = s1(i) - int(1);
    //     int end   = s2(i) - int(1);
    //     obs(end)  = Type(1.0);
    //     
    //     matrix<Type> Az = A;
    //     for (int r = 0; r < m - 1; ++r) {
    //       Type lp = Type(0.0);
    //       for (int j = 0; j < p; ++j) lp += z(i,j) * theta_cov(r * p + j);
    //       Type g = exp(lp);
    //       Az.row(r) *= g;
    //     }
    //     
    //     vector<Type> mu = expected_sojourn_exact(m, u(i), start, A);
    //     
    //     Type tau_eff = Type(0.0);
    //     for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);
    //     
    //     matrix<Type> tpm = atomic::expm(matrix<Type>(Az * tau_eff));
    //     vector<Type> pred = tpm.row(start).transpose();
    //     
    //     if (use_log_score)   total_score += log_score(pred, obs);
    //     if (use_brier_score) total_score += brier_score(pred, obs);
    //     if (use_rps_score)   total_score += rps_score(pred, obs);
    //   }
    // }
  } else if (cov_type == 2) {
    vector<Type> theta_cov = theta.segment(lambda_len + xii.size(), p * (m - 1));
    for (int i = 0; i < n; ++i) {
      vector<Type> obs(m); obs.setZero();
      int start = s1(i) - int(1);
      int end   = s2(i) - int(1);
      obs(end)  = Type(1.0);

      vector<Type> mu = expected_sojourn_exact(m, u(i), start, A);
      Type tau_eff = 0.0;
      for (int j = 0; j < m; ++j) tau_eff += mu(j) * xi(j);

      //vector<Type> mu_trans = mu.head(m-1);
      //vector<Type> w = mu_trans / mu_trans.sum();
      vector<Type> w = mu;
      w /= mu.sum();

      Type cov_linpred = 0.0;
      for (int ii = 0; ii < m - 1; ++ii) {
        Type lp_ii = 0.0;
        for (int j = 0; j < p; ++j) lp_ii += z(i,j) * theta_cov(ii * p + j);
        cov_linpred += w(ii) * lp_ii;
      }

      Type cov_scaling = exp(cov_linpred);
      tau_eff *= cov_scaling;

      matrix<Type> tpm = atomic::expm(matrix<Type>(A * tau_eff));
      vector<Type> pred = tpm.row(start).transpose();

      if (use_log_score)   total_score += log_score(pred, obs);
      if (use_brier_score) total_score += brier_score(pred, obs);
      if (use_rps_score)   total_score += rps_score(pred, obs);
    }
  }
  return total_score;
}

















