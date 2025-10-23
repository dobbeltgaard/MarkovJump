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
// void reach(const matrix<Type>& Q, int s, vector<int>& ok, bool forward=true){
//   int m = Q.rows(); ok.setZero(); ok(s)=1; bool ch=true;
//   while(ch){
//     ch=false;
//     for(int a=0;a<m;a++) if(ok(a)){
//       for(int b=0;b<m;b++){
//         if(a==b) continue;
//         bool edge = forward ? (Q(a,b)>Type(0)) : (Q(b,a)>Type(0));
//         if(edge && !ok(b)){ ok(b)=1; ch=true; }
//       }
//     }
//   }
// }
// 
// // Van Loan + pruning
// template<class Type>
// vector<Type> expected_sojourn_bridge(int m, Type u, int s1, int s2, matrix<Type> A){
//   vector<Type> mu(m); mu.setZero();
//   
//   matrix<Type> P = atomic::expm(matrix<Type>(A * u));
//   Type denom = P(s1, s2);
//   if(denom < Type(1e-14)) return mu; // impossible bridge 
//   
//   vector<int> from_s1(m), to_s2(m);
//   reach(A, s1, from_s1, true);
//   reach(A, s2, to_s2, false);
//   
//   matrix<Type> B(2*m,2*m); B.setZero();
//   for(int i=0;i<m;i++) for(int j=0;j<m;j++){
//     B(i,j)     = A(i,j);
//     B(m+i,m+j) = A(i,j);
//   }
//   
//   for(int k=0;k<m;k++){
//     if(!(from_s1(k) && to_s2(k))){ mu(k)=Type(0); continue; }
//     
//     for(int i=0;i<m;i++) for(int j=0;j<m;j++) B(i,m+j)=Type(0);
//     B(k, m+k) = Type(1);
//     
//     matrix<Type> E = atomic::expm(matrix<Type>(B * u));
//     mu(k) = E(s1, m + s2) / denom;
//   }
//   return mu;
// }



template<class Type>
vector<Type> expm_row_uniformized(const matrix<Type>& A, Type t, int s,  Type eps = Type(1e-12), int  maxN = 2000)
{
  const int m = A.rows();
  
  Type q = Type(0);
  for (int i = 0; i < m; ++i) {
    Type cand = -A(i,i);
    q = CppAD::CondExpGt(cand, q, cand, q);
  }
  q = q + Type(1e-12);
  
  matrix<Type> P(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      P(i,j) = (i == j ? Type(1) : Type(0)) + A(i,j) / q;
    }
  }
  
  vector<Type> r(m); r.setZero(); r(s) = Type(1);
  vector<Type> sum_row = r;
  
  Type qt = q * t;
  Type coeff = Type(1);                 // (q t)^n / n!, starting at n=0
  
  for (int n = 1; n <= maxN; ++n) {
    vector<Type> r_next(m); r_next.setZero();
    for (int j = 0; j < m; ++j) {
      Type acc = Type(0);
      for (int k = 0; k < m; ++k) acc += r(k) * P(k,j);
      r_next(j) = acc;
    }
    r = r_next;
    
    coeff *= qt / Type(n);
    
    vector<Type> incr(m);
    for (int j = 0; j < m; ++j) incr(j) = r(j) * coeff;
    
    sum_row += incr;
    
    double l1 = 0.0;
    for (int j = 0; j < m; ++j) l1 += std::fabs(asDouble(incr(j)));
    if (l1 < asDouble(eps)) break;
  }
  
  sum_row *= exp(-asDouble(qt));
  return sum_row;
}


template<class Type>
vector<Type> expected_sojourn_bridge(int m, Type u, int s1, int s2, matrix<Type> A){ 
  vector<Type> mu(m); mu.setZero();
  
  matrix<Type> P = atomic::expm(matrix<Type>(A * u));
  Type denom = P(s1, s2);
  Type eps = Type(1e-12);
  denom = (denom > eps) ? denom : eps;
  
  matrix<Type> B(2*m, 2*m); B.setZero();
  for (int i=0;i<m;i++) for (int j=0;j<m;j++){
    B(i,   j)    = A(i,j);
    B(m+i, m+j)  = A(i,j);
  }
  
  for (int k=0;k<m;k++){
    for (int i=0;i<m;i++) for (int j=0;j<m;j++) B(i, m+j) = Type(0);
    B(k, m+k) = Type(1);
    
    matrix<Type> E = atomic::expm(matrix<Type>(B * u));
    mu(k) = E(s1, m + s2) / denom;   // small numerators just stay small
  }
  return mu;
}

template<class Type>
vector<Type> expected_sojourn_bridge_uniformized(int m,Type T, int a,  int b,const matrix<Type>& Q,Type tail_eps = Type(1e-12),int  maxN     = 20000)
{
  vector<Type> mu(m); mu.setZero();
  
  Type q = Type(0);
  for (int i=0;i<m;++i) {
    Type cand = -Q(i,i);
    q = CppAD::CondExpGt(cand, q, cand, q);
  }
  q = q + Type(1e-12); // bump to avoid q=0
  matrix<Type> R(m,m);
  for (int i=0;i<m;++i)
    for (int j=0;j<m;++j)
      R(i,j) = (i==j ? Type(1) : Type(0)) + Q(i,j)/q;
  
  double lam = asDouble(q*T);
  int N = std::min(maxN, std::max(1, (int)(lam + 10.0*std::sqrt(std::max(0.0,lam)) + 50.0)));
  
  std::vector< vector<Type> > alpha; alpha.reserve(N+1);
  std::vector< vector<Type> > beta;  beta.reserve(N+1);
  
  vector<Type> alpha0(m); alpha0.setZero(); alpha0(a) = Type(1);
  vector<Type> beta0(m);  beta0.setZero();  beta0(b)  = Type(1);
  alpha.push_back(alpha0);
  beta.push_back(beta0);
  
  for (int n=1; n<=N; ++n) {
    // alpha_{n} = alpha_{n-1} * R
    vector<Type> a_next(m); a_next.setZero();
    for (int j=0;j<m;++j){
      Type acc = Type(0);
      for (int k=0;k<m;++k) acc += alpha.back()(k) * R(k,j);
      a_next(j) = acc;
    }
    alpha.push_back(a_next);
    
    vector<Type> b_next(m); b_next.setZero();
    for (int i=0;i<m;++i){
      Type acc = Type(0);
      for (int k=0;k<m;++k) acc += R(i,k) * beta.back()(k);
      b_next(i) = acc;
    }
    beta.push_back(b_next);
  }
  
  vector<Type> w(N+1); w.setZero();
  w(0) = exp( - Type(lam) );
  for (int n=0; n<N; ++n) {
    w(n+1) = w(n) * Type(lam) / Type(n+1);
  }
  
  Type denom = Type(0);
  for (int n=0; n<=N; ++n) denom += w(n) * alpha[n](b);
  Type denom_floor = Type(1e-12);
  denom = CppAD::CondExpGt(denom, denom_floor, denom, denom_floor);
  
  for (int k=0; k<m; ++k) {
    Type num_k = Type(0);
    for (int n=0; n<=N; ++n) {
      Type conv_nk = Type(0);
      for (int l=0; l<=n; ++l) conv_nk += alpha[l](k) * beta[n-l](k);
      num_k += (T / Type(n+1)) * w(n) * conv_nk;
    }
    mu(k) = num_k / denom;
  }
  return mu;
}


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
  
  Type alpha = Type(0.5);
  
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
  
  if (cov_type == 0) {
    for (int i = 0; i < n; ++i) {
      vector<Type> obs(m); obs.setZero();
      int start = s1(i) - int(1);
      int end   = s2(i) - int(1);
      obs(end)  = Type(1.0);
      
      //vector<Type> mu = expected_sojourn(m, u(i), start, A);
      //vector<Type> mu1 = expected_sojourn_bridge(m, u(i), start, end, A);
      //vector<Type> mu = (1-alpha)*mu0 + alpha*mu1;
      //vector<Type> mu=compute_conditional_sojourn(A, start, end, u(i));  // correct
      //vector<Type> mu = expected_sojourn_bridge_uniformized(m, u(i), start, end, A);
      vector<Type> mu = expected_sojourn_bridge(m, u(i), start, end, A);
      //vector<Type> mu0 = expected_sojourn(m, u(i), start, A);
      //vector<Type> mu1 = expected_sojourn_bridge(m, u(i), start, end, A);
      //vector<Type> mu = (1-alpha)*mu0 + alpha*mu1;
      
      Type tau_eff = 0.0; 
      for (int j = 0; j < xi.size(); ++j) {tau_eff += mu(j) * xi(j);} // dot product
      
      matrix<Type> tpm = atomic::expm( matrix<Type>(A*tau_eff) ); 
      vector<Type> pred = tpm.row(start).transpose();
      //vector<Type> pred = expm_row_uniformized(A, tau_eff, start);
      
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
      //vector<Type> mu=compute_conditional_sojourn(A, start, end, u(i));  // correct
      vector<Type> mu = expected_sojourn_bridge(m, u(i), start, end, A);
      //vector<Type> mu = expected_sojourn_bridge_uniformized(m, u(i), start, end, A);
      //vector<Type> mu0 = expected_sojourn(m, u(i), start, A);
      //vector<Type> mu1 = expected_sojourn_bridge(m, u(i), start, end, A);
      //vector<Type> mu = (1-alpha)*mu0 + alpha*mu1;
      Type tau_eff = 0.0; 
      for (int j = 0; j < xi.size(); ++j) {tau_eff += mu(j) * xi(j);} // dot product
      Type cov_scaling = exp(cov_linpred); //log(1 + exp(cov_linpred));  // softplus
      tau_eff *= cov_scaling;  // final time warp
      
      matrix<Type> tpm = atomic::expm( matrix<Type>(A*tau_eff) );
      vector<Type> pred = tpm.row(start).transpose();
      //vector<Type> pred = expm_row_uniformized(A, tau_eff, start);
      
      if (use_log_score) { total_score += log_score(pred, obs);}
      if (use_brier_score) {total_score += brier_score(pred, obs);}
      if (use_rps_score) {total_score += rps_score(pred, obs);}
    }
  }
  
  return total_score; /// Type(n);
}



