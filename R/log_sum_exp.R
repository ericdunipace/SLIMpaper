#' safely sums values on the exp scale from the log scale
#'
#' @param x vector to be summed
#'
#' @return real number
#' @export
log_sum_exp <- function(x) {
  # if(is.vector(x)) {
    if(all(is.infinite(x))) return(x[1])
    mx <- max(x)
    x_temp <- x - mx
    return(log(sum(exp(x_temp)))+ mx)
  # } else if (is.matrix(x)) {
  #   mx <- apply(x, 1, max)
  #   x_temp <- x - mx
  #   return(log(rowSums(exp(x_temp)))+ mx)
  # }
}

#' logged sum of exponentiated variables
#'
#' @param x a number
#' @param y another number
#'
#' @return \eqn{\log(\exp(x) + \exp(y))}
#' @export
log_sum_exp2 <- function(x,y) {
  mx <- pmax(x,y)
  # if(is.infinite(mx)) return(mx)

  temp <- cbind(x,y) - mx
  temp[mx == -Inf,] <- -Inf
  return(log(rowSums(exp(temp))) + mx)
}
