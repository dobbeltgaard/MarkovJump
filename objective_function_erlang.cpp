#include <TMB.hpp>

template<class Type>
vector<Type> transform(const vector<Type>& v) {
  return exp(v);
}


template<class Type>
matrix<Type> construct(const vector<Type>& Av) {
  const int m = Av.size()+1;
  
  matrix<Type> A(m, m); 
  A.fill(0.0);
  
  Type lambda;
  for (int i = 0; i < (m-1); i++) {
    lambda = Av(i);
    A(i, i  ) = -lambda;
    A(i, i+1) =  lambda;
  }
  
  return A;
}


template<class Type>
Type objective_function<Type>::operator() () {
  DATA_IVECTOR(s1);
  DATA_IVECTOR(s2);
  DATA_VECTOR(t);
  
  const vector<int> s1_ind((s1-1));
  const vector<int> s2_ind((s2-1));
  
  PARAMETER_VECTOR(Av);

  Type l((0));
  
  const vector<Type> Av_t((transform(Av)));
  
  const matrix<Type> A((construct(Av_t)));

  const int k = s1.size();
  for (int i = 0; i < k; i++) {
    l += log( atomic::expm( matrix<Type>(A*t(i)) )(s1_ind(i), s2_ind(i)) );
  }
  
  return -l;
}

