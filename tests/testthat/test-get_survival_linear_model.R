testthat::test_that("surv pred function works", {
  set.seed(130824)
  library(SLIMpaper)
  target <- get_survival_linear_model()

  n <- 128
  s <- 100
  p <- 11
  bsurv <- apply(matrix(runif(n*s),n,s),2,sort)
  theta <- matrix(rnorm(p*s), p, s)
  x <- matrix(rnorm(n*p), n, p)

  pred.theta <- list(baseline = bsurv, param = theta)

  pred.fun <- target$sel.pred.fun("cox")

  check <- pred.fun(x, pred.theta)
  surv <- do.call("rbind", lapply(1:n, function(i) matrix(bsurv[i,],n,s, byrow=TRUE)^exp(x %*% theta)))

  testthat::expect_equivalent(surv, check)

  #### check linpred
  pred.fun <- target$sel.pred.fun("linpred")
  testthat::expect_equivalent(x %*% theta, pred.fun(x,theta))

  #### check bart
  #not implemented yet

})
