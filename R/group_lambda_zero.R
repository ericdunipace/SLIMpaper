group_lambda_zero <- function(x, y, groups, pf) {
  xty <- Matrix::crossprod(Matrix::Matrix(x), y)
  l2 <- function(x){sqrt(sum(x))}
  norms <- tapply(Matrix::rowSums(xty^2), INDEX = factor(groups), FUN = l2)
  return(max(norms/pf))
}
