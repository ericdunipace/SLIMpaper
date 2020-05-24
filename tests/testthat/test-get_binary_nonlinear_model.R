testthat::test_that("modified friedman function works", {
  set.seed(092384)
  target <- get_binary_nonlinear_model()
  target$X$corr <- 0.5
  stan_dir <- "exec/Stan/logistic_horseshoe_noQR.stan"

  n <- 128
  s <- 100
  p <- 21
  p_star <- min(21, p)

  #param
  param <- target$rparam()


  X <- target$X$rX(n, target$X$corr, p)
  testthat::expect_silent(data <- target$rdata(n, X[,1:p_star,drop=FALSE], c(param$theta)[1:p_star],
                       param$sigma2, method = "modified.friedman", corr = target$X$corr))
  Y <- data$Y

  testthat::expect_warning(post_sample <- target$rpost(s, cbind(1, data$data_gen_functions( X ) ),
                              Y,
                              hyperparameters=list(),
                              method = "logistic",
                              stan_dir = stan_dir,
                              X.test = cbind(1, data$data_gen_functions(X ) ),
                              chains = 1, m0 = 20,
                              is.exponential = TRUE))

})

testthat::test_that("modified friedman check n of 1", {
  set.seed(092384)
  target <- get_binary_nonlinear_model()
  target$X$corr <- 0.5
  stan_dir <- "exec/Stan/logistic_horseshoe_noQR.stan"

  n <- 1
  s <- 100
  p <- 21
  p_star <- min(21, p)

  #param
  param <- target$rparam()


  X <- target$X$rX(n, target$X$corr, p)
  testthat::expect_silent( target$rdata(n, X[,1:p_star,drop=FALSE], c(param$theta)[1:p_star],
                       param$sigma2, method = "modified.friedman", corr = target$X$corr))

})


testthat::test_that("nn model works", {
  set.seed(9080808)
  target <- get_binary_nonlinear_model()
  target$X$corr <- 0.5

  n <- 128
  s <- 100
  p <- 21
  p_star <- min(21, p)

  param <- target$rparam()


  X <- target$X$rX(n, target$X$corr, p)
  dat <- target$rdata( n=n, x=X[,1:p_star,drop=FALSE], theta = c(param$theta)[1:p_star],
                                        method = "glm",
               corr = target$X$corr)

  hyperparam <- list(learning.rate = 1e-2,
                     lambda = 0.1,
                     batch.size = 32,
                     first.layer.width = 10,
                     hidden.layer.width = 10)

  # debugonce(target$rpost)
  # debugonce(nn_train)
  testthat::expect_warning(target$rpost(n.samp = s, x = X, y = dat$Y, hyperparameters = hyperparam,
               method = "nn",
               niter = 10, test.portion = 0, verbose = FALSE))


  X.test <- target$X$rX(32, target$X$corr, p)
  testthat::expect_warning(target$rpost(n.samp = s, x = X, y = dat$Y, hyperparameters = hyperparam,
               X.test = X.test,
               method = "nn",
               niter = 10, test.portion = 0, verbose = FALSE))


})
