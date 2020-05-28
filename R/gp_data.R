gp_data <- function(x,d,nT, times) {
  nS <- ncol(x)
  y <- matrix(NA, nrow = nT*nS, ncol=d)
  time_idx <- rep(as.numeric(factor(times)), each = nS)
  for(i in 1:nT) {
    y_rows <- 1:nS + (i-1) * nS
    x_rows <- 1:d + (i-1) * d
    y[y_rows, ] <- t(x[x_rows,])
  }
  return(list(y = y, time = time_idx))
}
