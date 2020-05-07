construct_interp_survival_matrix <- function(x, nT, intercept = FALSE) {

  n <- nrow(x)
  d <- ncol(x)

  xmat <- matrix(0, nrow=n * nT, ncol = d * nT)

  for(t in 1:nT) {
    rows <- seq(1, n * nT, nT) + (t-1)
    cols <- 1:d + (t-1) * d
    xmat[rows, cols] <- x
  }
  if(intercept) xmat <- cbind(1,xmat)
  return(xmat)
}
