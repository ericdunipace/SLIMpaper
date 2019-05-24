#include "CoarsePosteriorSummary_types.h"

void distance_c(const refMatConst & A, const refMatConst & B, Rcpp::NumericMatrix & dist_matrix) {
  for (int i = 0; i < A.cols(); i++) {
    for (int j = 0; j < B.cols(); j++) {
      vector bvec = B.col(j);
      dist_matrix(i,j) = (A.col(i)-bvec).norm();
    }
  }
}

//[[Rcpp::export]]
Rcpp::NumericMatrix distance(const matrix & A_, const matrix & B_) {
  int N = A_.rows();
  int M = B_.rows();

  matrix A = A_.transpose();
  matrix B = B_.transpose();

  Rcpp::NumericMatrix dist_matrix(N,M);

  distance_c(A,B,dist_matrix);

  return dist_matrix;
}
