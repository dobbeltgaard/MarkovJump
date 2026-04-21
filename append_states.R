#ex

rm(list = ls()); gc(); 
library(MASS); library(Rcpp); library(RcppEigen); library(TMB); #library(trust)
sourceCpp("FUNCS_MJP_with_eigen.cpp")

m <- 4
A <- make_A1(m, rep(1:5^3, 5^3))
B = expand_A1(m, A, 1)
A
B

foo = expand_dist(m, t(c(0,1,0,0)), 1)
foo

collapse_dist(m, foo, 2)

B
expand_dist_weighted(m, t(c(1,0,0,0,0)), 1, B)
expand_dist_weighted(m, t(c(0,1,0,0,0)), 1, B)
expand_dist_weighted(m, t(c(0,0,1,0,0)), 1, B)
expand_dist_weighted(m, t(c(0,0,0,1,0)), 1, B)
expand_dist_weighted(m, t(c(0,0,0,0,1)), 1, B)



tpm = as.matrix(Matrix::expm(B*0.01))
collapse_dist(m, foo %*% tpm, 1)


c(0,0,1,0,0) %*% as.matrix(Matrix::expm(A*0.01))


expand_A = function(m, A, k){
  M = (m-1)*k + m - k + 1
  B = matrix(0, nrow = M, ncol = M)

  for(i in 1:m){
    for(j in 1:m){
      #diagonal blocks
      if(j - i == 1 ){#& i < (m-1)){
        start = i * (k+1) - k
        end = start + k
        B[start:end, (start:end)+1] = diag(A[i,j], nrow = k+1)
      }
      if(j - i == 1 & i == (m-1)){
        start = i * (k+1) - k
        end = start
        B[start:end, (start:end)+1] = A[i,j]
      }
      #enroute_ij
      else if(j-i > 1 & j < m) {
        enroute = matrix(0,nrow = k+1, ncol = k+1)
        enroute[, 1]= A[i,j]
        rowstart = i * (k+1) - k
        rowend = rowstart + k
        colstart = 2+(k+1)*(j - 2)
        colend = colstart + k
        B[rowstart:rowend, colstart:colend] = enroute
      }
      else if(j-i > 1 & j == m){
        enroute = matrix(0,nrow = k+1, ncol = 1)
        enroute[,1] = A[i,j]
        rowstart = i * (k+1) - k
        rowend = rowstart + k
        colstart = 2+(k+1)*(j - 2)
        colend = colstart
        B[rowstart:rowend, colstart:colend] = enroute
      }
    }
  }
  diag(B) <- -rowSums(B)
  return(B)
}
expand_A(4, make_A1(4, rep(1:5^3, 5^3)), k = 1)

