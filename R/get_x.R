#### X variable ####
gen_x <- function() {
  rX <- function(n, corr, p, ...) {
    p <- p-1 # leave room for intercept
    if(is.matrix(corr)) {
      x_cov <- corr
    } else {
      x_cov <- matrix(corr, nrow=p, ncol=p)
    }
    diag(x_cov) <- 1
    if(p > 5) {
      p_lim <- 5
      p_lim1 <- p_lim + 1
      x_cov[p_lim1:p,1:p_lim] <- x_cov[1:p_lim, p_lim1:p] <- 0
    }
    if(p > 10) {
      p_lim <- 10
      p_lim1 <- p_lim + 1
      x_cov[p_lim1:p,1:p_lim] <- x_cov[1:p_lim, p_lim1:p] <- 0
    }
    if(p > 20) {
      p_lim <- 20
      p_lim1 <- p_lim + 1
      x_cov[p_lim1:p,1:p_lim] <- x_cov[1:p_lim, p_lim1:p] <- 0
    }
    if(p > 50) {
      p_lim <- 50
      p_lim1 <- p_lim + 1
      x_cov[p_lim1:p,1:p_lim] <- x_cov[1:p_lim, p_lim1:p] <- 0
    }
    X <- cbind(1,mvtnorm::rmvnorm(n, rep(0,p), x_cov))
    return(X)
  }

  scaleX <- function(x, sds,...) {
    if(all(x[,1]==1)) {
      stopifnot(length(sds) == (ncol(x)-1))
      x[,-1] <- x[,-1] %*% diag(sds)
    } else {
      stopifnot(length(sds) == ncol(x))
      x <- x %*% diag(sds)
    }
    return(x)
  }
  return(list(rX=rX, sX = scaleX))
}
