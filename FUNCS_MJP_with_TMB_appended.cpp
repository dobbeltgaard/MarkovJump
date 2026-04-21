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
  for (int i = 0; i < m; ++i) {A(i, i) = -A.row(i).sum(); }
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
matrix<Type> expand_A1(int m, matrix<Type> A, int k){
  int kp1 = k + 1;
  int M = (m - 1) * (k + 1) + 1;//(m - 1) * k + m - k + k;
  matrix<Type> B(M, M);
  B.setZero();
  for(int i = 0; i < m; ++i){
    int rowstart = i * kp1;
    int rowend = rowstart + k;
    for(int j = 0; j < m; ++j){
      Type aij = A(i, j);
      if(j - i == 1){
        for(int t = 0; t <= k; ++t){
          int r = rowstart + t;
          int c = rowstart + 1 + t;
          if(r < M && c < M) B(r, c) = aij;
        }
      }
    }
  }
  vector<Type> rs(M);
  rs.setZero();
  for(int r = 0; r < M; ++r){
    for(int c = 0; c < M; ++c) rs(r) += B(r, c);
  }
  for(int r = 0; r < M; ++r) B(r, r) = -rs(r);
  return B;
}

template<class Type>
int expanded_size(int m, int k){
  if(m <= 1) return 1;
  return (m - 1) * (k + 1) + 1;
}

template<class Type>
vector<Type> expand_dist(int m, vector<Type> Pt, int k){
  int M = expanded_size<Type>(m, k);
  vector<Type> y(M);
  y.setZero();
  int pos = 0;
  for(int i = 0; i < m - 1; ++i){
    int b = k + 1;
    Type v = Pt(i) / Type(b);
    for(int r = 0; r < b; ++r) y(pos + r) = v;
    pos += b;
  }
  y(pos) = Pt(m - 1);
  return y;
}

template<class Type>
vector<Type> collapse_dist(int m, vector<Type> Pt, int k){
  vector<Type> z(m);
  z.setZero();
  int pos = 0;
  for(int i = 0; i < m - 1; ++i){
    int b = k + 1;
    Type s = 0;
    for(int r = 0; r < b; ++r) s += Pt(pos + r);
    z(i) = s;
    pos += b;
  }
  z(m - 1) = Pt(pos);
  return z;
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
  DATA_INTEGER(k);
  PARAMETER_VECTOR(theta);
  
  const int n = s1.size();
  Type total_score = 0.0;
  
  int base_len = 0;
  matrix<Type> A(m, m);
  {
    switch (generator_type) {
    case 0: case 1: {
    base_len = m - 1;
    vector<Type> theta_base = theta.segment(0, base_len);
    vector<Type> lambda = softplus(theta_base);
    A = (generator_type == 0) ? make_A1(m, lambda) : make_A2(m, lambda);
  } break;
    case 2: {
      base_len = m * (m - 1) / 2;
      vector<Type> theta_base = theta.segment(0, base_len);
      vector<Type> lambda = softplus(theta_base);
      A = make_A3(m, lambda);
    } break;
    case 3: {
      base_len = 2 * m - 3;
      vector<Type> theta_base = theta.segment(0, base_len);
      vector<Type> lambda = softplus(theta_base);
      A = make_A4(m, lambda);
    } break;
    case 4: {
      base_len = 3 * m - 6;
      vector<Type> theta_base = theta.segment(0, base_len);
      vector<Type> lambda = softplus(theta_base);
      A = make_A5(m, lambda);
    } break;
    default: error("Invalid generator_type");
    }
  }
  
  if (cov_type == 0) {
    for (int i = 0; i < n; ++i) {
      int start = s1(i) - 1;
      int end   = s2(i) - 1;
      
      vector<Type> obs(m); obs.setZero();
      obs(end) = Type(1);
      
      matrix<Type> B = expand_A1(m, A, k);
      matrix<Type> tpm = atomic::expm( matrix<Type>(B * u(i)) );
      
      vector<Type> init_obs(m); init_obs.setZero(); init_obs(start) = Type(1);
      vector<Type> init = expand_dist(m, init_obs, k);
      
      matrix<Type> tpmT = tpm.transpose();
      vector<Type> Ptu = vector<Type>(tpmT * init);//vector<Type> Ptu = tpm * init;
      vector<Type> pred = collapse_dist(m, Ptu, k);
      
      if (use_log_score)   total_score += log_score(pred, obs);
      if (use_brier_score) total_score += brier_score(pred, obs);
      if (use_rps_score)   total_score += rps_score(pred, obs);
    }
  } else {
    vector<Type> theta_cov = theta.segment(base_len, theta.size() - base_len);
    for (int i = 0; i < n; ++i) {
      int start = s1(i) - 1;
      int end   = s2(i) - 1;
      
      vector<Type> obs(m); obs.setZero();
      obs(end) = Type(1);
      
      Type cov_linpred = 0.0;
      for (int j = 0; j < z.cols(); ++j) cov_linpred += z(i,j) * theta_cov(j);
      Type cov_time = exp(cov_linpred) * u(i);
      
      matrix<Type> B = expand_A1(m, A, k);
      matrix<Type> tpm = atomic::expm( matrix<Type>(B * cov_time) );
      
      vector<Type> init_obs(m); init_obs.setZero(); init_obs(start) = Type(1);
      vector<Type> init = expand_dist(m, init_obs, k);
      
      matrix<Type> tpmT = tpm.transpose();
      vector<Type> Ptu = vector<Type>(tpmT * init); //vector<Type> Ptu = tpm * init;
      vector<Type> pred = collapse_dist(m, Ptu, k);
      
      if (use_log_score)   total_score += log_score(pred, obs);
      if (use_brier_score) total_score += brier_score(pred, obs);
      if (use_rps_score)   total_score += rps_score(pred, obs);
    }
  }
  
  SIMULATE {
    matrix<Type> pred_mat(n, m);
    pred_mat.setZero();
    
    matrix<Type> B = expand_A1(m, A, k); // generator in expanded space
    
    for (int i = 0; i < n; ++i) {
      int start = s1(i) - 1;
      
      Type time_eff = u(i);
      if (cov_type != 0) {
        vector<Type> theta_cov = theta.segment(base_len, theta.size() - base_len);
        Type cov_linpred = Type(0);
        for (int j = 0; j < z.cols(); ++j) cov_linpred += z(i, j) * theta_cov(j);
        time_eff = exp(cov_linpred) * u(i);
      }
      
      matrix<Type> tpm = atomic::expm(matrix<Type>(B * time_eff));
      matrix<Type> tpmT = tpm.transpose();
      
      vector<Type> init_obs(m);
      init_obs.setZero();
      init_obs(start) = Type(1);
      
      vector<Type> init = expand_dist(m, init_obs, k);
      vector<Type> Ptu  = vector<Type>(tpmT * init);
      vector<Type> pred = collapse_dist(m, Ptu, k);
      
      for (int j = 0; j < m; ++j) pred_mat(i, j) = pred(j);
    }
    
    REPORT(pred_mat);
  }
  
  return total_score;
}
