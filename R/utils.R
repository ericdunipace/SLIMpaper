augPseudo <- function(X, mu, theta, theta_norm, pseudo.obs, n, same=TRUE) {
  augDat <- augmentMatrix(X, mu, theta, same=same)
  augDat$XtX <- n/(n + pseudo.obs) * augDat$XtX + diag(theta_norm) * pseudo.obs/(n + pseudo.obs)
  augDat$XtY <- n/(n + pseudo.obs) * augDat$XtY + theta_norm * pseudo.obs/(n + pseudo.obs)
  return(augDat)
}

set_penalty_factor <- function(x, method){
  distances <- switch(method,
                      expectation = abs(colMeans(x)),
                      distance = sqrt(apply(x,2,crossprod)),
                      none = rep(1, ncol(x)),
                      wt_expect = abs(colMeans(x)/colSD(x))
  )
  penalties <- 1/distances
  penalties[1] <- 0
  return(penalties)
}

calc.lambdas <- function(aug, lambda.min.ratio, penfact, n.lambda) {
  xty <- aug$XtY
  xtx <- aug$XtX

  temp.xty <- abs(xty)/penfact
  temp.xty[penfact==0] <- 0
  max.lamb <- max(temp.xty)
  lambdas <- c(exp(seq(log(max.lamb), log(max.lamb) + log(lambda.min.ratio), length.out=n.lambda)),0)
  return(lambdas)
}
