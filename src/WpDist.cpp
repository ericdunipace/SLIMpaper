#include "CoarsePosteriorSummary_types.h"
#include "sort.h"

using namespace Rcpp;

double W2dist_iid(refMat X, refMat Y){

  if(X.cols() != Y.cols()) {
    Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
  }

  if(X.rows() != Y.rows()) {
    Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
  }

  sort_matrix(X);
  sort_matrix(Y);

  return (X-Y).cwiseAbs2().mean();
}

double Wpdist_iid(refMat X, refMat Y, double p){

  if(X.cols() != Y.cols()){
    Rcpp::stop("Number of columns of first matrix don't match number of columns of second matrix");
  }

  if(X.rows() != Y.rows()){
    Rcpp::stop("Number of rows of first matrix don't match number of rows of second matrix");
  }

  sort_matrix(X);
  sort_matrix(Y);

  return (X-Y).array().pow(p).mean();
}

//[[Rcpp::export]]
double W2dist_normal(const SEXP & X_, const SEXP & Y_, double p){

  const matMap X(as<matMap >(X_));
  const matMap Y(as<matMap >(Y_));

  if(X.cols() != Y.cols()){
    Rcpp::stop("Number of columns of first matrix doesn't match number of columns of second matrix");
  }

  if(X.rows() != Y.rows()){
    Rcpp::stop("Number of rows of first matrix doesn't match number of rows of second matrix");
  }
  double df = double(X.cols()) - 1.0;
  vector mu_x = X.rowwise().mean();
  vector mu_y = Y.rowwise().mean();

  matrix cent_x = X.array().colwise() - mu_x.array();
  matrix cent_y = Y.array().colwise() - mu_y.array();

  matrix cov_x = cent_x * cent_x.transpose();
  matrix cov_y = cent_y * cent_y.transpose();
  cov_x.array() /= df;
  cov_y.array() /= df;

  Eigen::SelfAdjointEigenSolver<matrix> root1(cov_x);
  matrix sqrt_x = root1.operatorSqrt();
  // matrix check =  sqrt_x * sqrt_x;
  // Rcout << check(0,0) << " " << cov_x(0,0);
  matrix temp = sqrt_x * cov_y * sqrt_x;
  Eigen::SelfAdjointEigenSolver<matrix> root2(temp);
  matrix root_temp = root2.operatorSqrt();
  double mean_diff = (mu_x - mu_y).squaredNorm();
  double cov_diff = (0.5 * (cov_x + cov_y) -  root_temp).trace();

  return mean_diff + cov_diff;
}

//[[Rcpp::export]]
double Wasserstein_p_iid(const SEXP & X_, const SEXP & Y_, double p) {

  const matMap Xtemp(as<matMap >(X_));
  const matMap Ytemp(as<matMap >(Y_));

  matrix X = Xtemp;
  matrix Y = Ytemp;

  if(p == 2.0){
    // Rcout << "w2";
    return W2dist_iid(X,Y);
  } else {
    // Rcout << "wp";
    return Wpdist_iid(X,Y,p);
  }
}


//[[Rcpp::export]]
NumericVector W2_idx(const NumericMatrix & gamma_,
                    const NumericMatrix & X_, const NumericMatrix & theta_,
                    const NumericMatrix & mu_,
                    const IntegerVector & idx_,
                    int max_p){

  const matMap X(as<matMap >(X_));
  const matMap mu(as<matMap >(mu_));
  const matMap theta(as<matMap >(theta_));
  const matMap gamma(as<matMap >(gamma_));

  matrix muT(X.rows(),theta.rows());

  if(mu.rows() != theta.rows()) {
    // Rcout << "transpose" << std::endl;
    muT = mu.transpose();
  } else {
    muT = mu;
  }


  IntegerVector idx = idx_ - 1;

  int N = idx.size();
  int P = idx_(N-1);

  NumericVector out(max_p);
  out.fill(NumericVector::get_na());

  for(int i = 0; i < N; i++) {
    if(idx_(i) > P) continue;
    matrix product = theta * gamma.col(i).asDiagonal() * X.transpose();
    // Rcout << product.rows() << std::endl;
    // Rcout << product.cols() <<std::endl;
    // Rcout << muT.rows() << std::endl;
    // Rcout << muT.cols() << std::endl;
    out(idx(i)) = W2dist_iid(product, muT);
  }
  return(out);
}
