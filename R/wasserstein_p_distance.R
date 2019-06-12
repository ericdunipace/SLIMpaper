WpDist_iid <- function(X,Y, p, observations = c("colwise","rowwise")) {

  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  obs <- match.arg(observations)
  if(obs == "rowwise"){
    X <- t(X)
    Y <- t(Y)
  }

  loss <- Wasserstein_p_iid(X,Y,p)
  return(loss^(1/p))
}

wass_trajectory <- function(traj, compare, p=2, idx, max) {

  dist <- rep(NA, max)
  B <- transport::pp(compare)

  val <- NA

  for (i in seq_along(idx) ) {
    val <- traj[[i]]
    if(nrow(val) != nrow(compare)) val <- t(val)
    A <- transport::pp(val)
    dist[idx[i]] <- transport::wasserstein(A, B,p=p)
  }
  return(dist)
}

WpDist_individual <- function(X,Y, p, observations = c("colwise","rowwise")) {
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(Y)) Y <- as.matrix(Y)
  obs <- match.arg(observations)
  if(obs == "rowwise"){
    X <- t(X)
    Y <- t(Y)
  }

  Xs <- apply(X,2,sort)
  Ys <- apply(Y,2,sort)

  loss <- colMeans((Xs - Ys)^p)

  return(loss^(1/p))

}

# WpDist_norm <- function(X,Y, p) {
#
#   if(!is.matrix(X)) X <- as.matrix(X)
#   if(!is.matrix(Y)) Y <- as.matrix(Y)
#   if(p != 2) stop("p must equal 2. Distances other than p=2 not yet implemented.")
#   # loss <- W2dist_normal(X,Y,p)
#
#   mu_x <- rowMeans(X)
#   mu_y <- rowMeans(Y)
#
#   cov_x <- cov(t(X))
#   cov_y <- cov(t(Y))
#
#   svd_x <- svd(cov_x)
#   sqrt_x <- svd_x$u %*% diag(sqrt(svd_x$d)) %*% t(svd_x$v)
#   temp <- sqrt_x %*% cov_y %*% sqrt_x
#   svd_t <- svd(temp)
#   sqrt_t <- svd_t$u %*% diag(sqrt(svd_t$d)) %*% t(svd_t$v)
#
#   loss <- sum((mu_x - mu_y)^2) + sum(diag(0.5 *(cov_x + cov_y - sqrt_t)))
#
#   return(loss)
# }
