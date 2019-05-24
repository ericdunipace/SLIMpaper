#include "CoarsePosteriorSummary_types.h"

using namespace Rcpp;

// [[Rcpp::export]]
List bayesConjRegNormal(int n_samp, const List & hyperparameters,
                        const NumericVector & Y_, const NumericMatrix X_){
  double a = double(hyperparameters["alpha"]);
  double b = double(hyperparameters["beta"]);
  const Eigen::Map<Eigen::VectorXd> m(as<Eigen::Map<Eigen::VectorXd> >(hyperparameters["mu"]));
  int n = Y_.size();
  int p = m.size();
  const Eigen::Map<Eigen::MatrixXd> Y(as<Eigen::Map<Eigen::MatrixXd> >(Y_));
  const Eigen::Map<Eigen::MatrixXd> X(as<Eigen::Map<Eigen::MatrixXd> >(X_));
  const Eigen::Map<Eigen::MatrixXd> Lambda(as<Eigen::Map<Eigen::MatrixXd> >(hyperparameters["Lambda"]));

  matrix post_prec_beta = X.transpose() * X + Lambda;
  // post_prec_beta.selfadjointView<Eigen::Lower>().rankUpdate(X.transpose() , 1.0);
  // post_prec_beta += Lambda.selfadjointView<Eigen::Lower>();

  Eigen::LLT<matrix> post_L_beta(post_prec_beta.selfadjointView<Eigen::Lower>().llt().solve( matrix::Identity(p,p)));
  vector post_mu = post_prec_beta.selfadjointView<Eigen::Lower>().ldlt().solve(X.transpose() * Y + Lambda * m);

  double post_a = a + double(n) * 0.5;
  double post_b = b + 0.5 * (Y.squaredNorm() - post_mu.transpose() * post_prec_beta * post_mu );

  vector post_sigma(n_samp);
  for(int i = 0; i < n_samp; i++) post_sigma(i) = 1/R::rgamma(post_a, 1/post_b);

  matrix post_theta(p, n_samp);
  vector Z(p);

  for(int i = 0; i < n_samp; i++) {
    for(int j = 0; j < p; j++) Z(j) = R::rnorm(0.0,1.0);
    post_theta.col(i) = post_mu + post_L_beta.matrixL() * Z * std::sqrt(post_sigma(i));
  }


  return Rcpp::List::create(Rcpp::Named("theta") = Rcpp::wrap(post_theta),
                            Rcpp::Named("sigma") = Rcpp::wrap(post_sigma));
}


/*** R
if(!("rmvnorm" %in% ls())) Rcpp::sourceCpp("dmvnorm.cpp")
set.seed(121)
p <- 10
n <- 1
n.samp <- 1E5

theta <- seq(-2,2, length.out = p)

x <- matrix(rnorm(p*n), nrow=n, ncol=p)
y <- rnorm(n, x %*% theta, sd=1)
a <- 10
b <- 10
m <- rep(0,p)
s <- diag(1,p,p)

Lambda <- solve(s)

post_prec_beta <- crossprod(x) + Lambda
post_cross_beta <- crossprod(x, y) + Lambda %*% m

post_mu <- solve(post_prec_beta, post_cross_beta)

post_a <- a + n * 0.5
post_b <- b + 0.5 * (sum(y^2) -  t(post_mu) %*% post_prec_beta %*% post_mu )

post_sigma <- 1/rgamma(n.samp, post_a, post_b)

post_theta <- sapply(1:n.samp, function(i) rmvnorm(1, post_mu, solve(post_prec_beta)*post_sigma[i]))

hyperparameters <- list()
hyperparameters$alpha <- a
hyperparameters$beta <- b
hyperparameters$Lambda <- Lambda
hyperparameters$mu <- m
hyperparameters$sigma <- s

check <- bayesConjRegNormal(n.samp, hyperparameters,
                            y, x)

# print(chol(solve(post_prec_beta)))

print(rbind(rowMeans(check$theta), rowMeans(post_theta)))

print(c(mean(check$sigma), mean(post_sigma)))

*/
