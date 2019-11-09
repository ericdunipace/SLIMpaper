#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rmvnorm(int nsamples, const NumericVector & mean, const NumericMatrix & covariance){
  RNGScope scope;
  int ncols = covariance.cols();
  const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
  Eigen::MatrixXd cholesky_covariance(covariance_.llt().matrixU());
  Eigen::MatrixXd Y(nsamples, ncols);
  for(int i = 0; i < ncols; i++){
    Y.col(i) = as<Eigen::ArrayXd>(rnorm(nsamples));
  }
  Y = Y * cholesky_covariance;
  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nsamples; i++){
      Y(i,j) = Y(i,j) + mean(j);
    }
  }
  return wrap(Y);
}

// [[Rcpp::export]]
NumericVector dmvnorm(const NumericMatrix & x, const NumericVector & mean, const NumericMatrix & covariance){
  const Eigen::Map<Eigen::MatrixXd> covariance_(as<Eigen::Map<Eigen::MatrixXd> >(covariance));
  const Eigen::Map<Eigen::MatrixXd> x_(as<Eigen::Map<Eigen::MatrixXd> >(x));
  Eigen::LLT<Eigen::MatrixXd> lltofcov(covariance_);
  Eigen::MatrixXd lower = lltofcov.matrixL();
  Eigen::MatrixXd xcentered(x_);
  double halflogdeterminant = lower.diagonal().array().log().sum();;
  double cst = - (halflogdeterminant) - (x.cols() * 0.9189385);
  for(int j = 0; j < x.cols(); j++){
    for(int i = 0; i < x.rows(); i++){
      xcentered(i,j) = xcentered(i,j) - mean(j);
    }
  }
  Eigen::VectorXd results = -0.5 * lower.triangularView<Eigen::Lower>().solve(xcentered.transpose()).colwise().squaredNorm();
  for (int i = 0; i < results.size(); i++){
    results(i) = results(i) + cst;
  }
  return wrap(results);
}

