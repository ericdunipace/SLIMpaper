#include "CoarsePosteriorSummary_types.h"

using namespace Rcpp;

template <typename T, typename U>
matrix diag_pre_multiply_C(T & V, U & X){
  return X * V.asDiagonal();
}

// template <typename E>
// double mseM(matrix & mu, E & est) {
//   return (mu-est).array().square().mean();
// }

// double mseT(vector & mu, matrix & est) {
//   return ( est.colwise() - mu ).array().square().mean();
// }

template <typename E>
double mse(const vector & mu, E & est) {
  return ( est.colwise() - mu ).array().square().mean();
}

// [[Rcpp::export]]
double mse_C(vector & mu, matrix & est) {
  return ( est.colwise() - mu ).array().square().mean();
}

// [[Rcpp::export]]
NumericVector mse_idx_dist(NumericMatrix & V_, NumericMatrix & X_,
                  NumericMatrix & theta_, NumericVector & mu_,
                  IntegerVector & idx_, int max_p) {

  //maps to Eigen
  const Eigen::Map<matrix> V(as<Eigen::Map<matrix> >(V_));
  const Eigen::Map<matrix> X(as<Eigen::Map<matrix> >(X_));
  const Eigen::Map<matrix> theta(as<Eigen::Map<matrix> >(theta_));
  const Eigen::Map<vector> mu(as<Eigen::Map<vector> >(mu_));

  IntegerVector idx = idx_ - 1;
  int n = idx.size();
  int p = idx_(n-1);

  NumericVector out(max_p);
  out.fill(NumericVector::get_na());

  for(int i = 0; i < n; i++) {
    if(idx_(i) > p) continue;
    out(idx(i)) = mse(mu, X * diag_pre_multiply_C(V.col(i),theta).transpose());
    // out(idx(i)) = mseT(mu, X);

  }
  return(out);
}

// // [[Rcpp::export]]
// NumericVector mse_mat_dist(NumericMatrix & V_, NumericMatrix & X_,
//                            NumericMatrix & theta_, NumericMatrix & mu_,
//                            IntegerVector & idx_, int max_p) {
//
//   const Eigen::Map<matrix> V(as<Eigen::Map<matrix> >(V_));
//   const Eigen::Map<matrix> X(as<Eigen::Map<matrix> >(X_));
//   const Eigen::Map<matrix> theta(as<Eigen::Map<matrix> >(theta_));
//   const Eigen::Map<matrix> mu(as<Eigen::Map<matrix> >(mu_));
//   IntegerVector idx = idx_ - 1;
//   int n = idx.size();
//   int p = idx_(n-1);
//
//   NumericVector out(max_p);
//   out.fill(NumericVector::get_na());
//
//   for(int i = 0; i < n; i++) {
//     if(idx_(i) > p) continue;
//     out(idx(i)) = mseM(mu, X * diag_pre_multiply_C(V.col(i),theta).transpose());
//   }
//   return(out);
// }

// [[Rcpp::export]]
NumericVector mse_idx_expect(NumericMatrix & V_, NumericMatrix & X_,
                           NumericMatrix & theta_, NumericVector & mu_,
                           IntegerVector & idx_, int max_p) {

  const Eigen::Map<matrix> V(as<Eigen::Map<matrix> >(V_));
  const Eigen::Map<matrix> X(as<Eigen::Map<matrix> >(X_));
  const Eigen::Map<matrix> theta(as<Eigen::Map<matrix> >(theta_));
  const Eigen::Map<vector> mu(as<Eigen::Map<vector> >(mu_));
  IntegerVector idx = idx_ - 1;
  int n = idx.size();
  int p = idx_(n-1);

  NumericVector out(max_p);
  out.fill(NumericVector::get_na());

  for(int i = 0; i < n; i++) {
    if(idx_(i) > p) continue;
    out(idx(i)) = mse(mu, X * (diag_pre_multiply_C(V.col(i),theta).colwise().mean()).transpose());
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector mse_idx(NumericMatrix & V_, NumericMatrix & X_,
                      NumericVector & mu_,
                      IntegerVector & idx_, int max_p) {

  const Eigen::Map<matrix> V(as<Eigen::Map<matrix> >(V_));
  const Eigen::Map<matrix> X(as<Eigen::Map<matrix> >(X_));
  const Eigen::Map<vector> mu(as<Eigen::Map<vector> >(mu_));
  IntegerVector idx = idx_ - 1;
  int n = idx.size();
  int p = idx_(n-1);

  NumericVector out(max_p);
  out.fill(NumericVector::get_na());

  for(int i = 0; i < n; i++) {
    if(idx_(i) > p) continue;
    out(idx(i)) = mse(mu, X * V.col(i));
  }
  return(out);
}

// NumericVector mse_mat(NumericMatrix & V_, NumericMatrix & X_,
//                       NumericMatrix & mu_,
//                       IntegerVector & idx_, int max_p) {
//
//   const Eigen::Map<matrix> V(as<Eigen::Map<matrix> >(V_));
//   const Eigen::Map<matrix> X(as<Eigen::Map<matrix> >(X_));
//   const Eigen::Map<matrix> mu(as<Eigen::Map<matrix> >(mu_));
//   IntegerVector idx = idx_ - 1;
//   int n = idx.size();
//   int p = idx_(n-1);
//
//   NumericVector out(max_p);
//   out.fill(NumericVector::get_na());
//
//   for(int i = 0; i < n; i++) {
//     if(idx_(i) > p) continue;
//     out(idx(i)) = mse(mu, X * V.col(i));
//   }
//   return(out);
// }

/*** R
set.seed(11)
n <- 1000
p <- 10
samps <- 100

X_test <- matrix(rnorm(n*p), nrow=n, ncol=p)
theta <- matrix(rnorm(p * samps), nrow=samps, ncol = p)
gamma <- diag(1,p,p)
idx <- c(1L,3L,10L)
gamma <- gamma[,idx]
E_theta <- colMeans(theta)
mu <- X_test %*% E_theta
test_out <- sapply(1:length(idx), function(i) mean((matrix(mu,nrow=n, ncol=samps) - X_test %*%t(theta %*% diag(gamma[,i],p,p)))^2))

out <- mse_C(gamma, X_test, theta, mu, idx)
print(test_out)
print(out)

*/
