skip_if_no_functions <- function() {
  have_scipy <- reticulate::py_module_available("scipy")
  have_numpy <- reticulate::py_module_available("numpy")
  have_torch <- reticulate::py_module_available("torch")
  if (!have_scipy)
    testthat::skip("scipy not available for testing")
  if (!have_numpy)
    testthat::skip("numpy not available for testing")
  if (!have_torch)
    testthat::skip("torch not available for testing")
}


testthat::test_that("nn python works", {
  skip_if_no_functions()
  set.seed(23408234)
  n <- (2L^10L)
  d <- 10L
  batch.size <- as.integer(2^7)
  first.layer.width <- (d*10L)
  hidden.layer.width <- 10L
  num.epoch <- 10L
  lr <- 1e-2
  lambda <- 0.01

  x <- matrix(rnorm(n*d), n, d)
  prob <- plogis(x %*% rnorm(d)/d*10 + rnorm(1))
  y <- as.matrix(rbinom(n, 1, prob))

  # debugonce(nn_train)
  testthat::expect_silent(res <- nn_train(x,y, niter = 10,
                  learning.rate = lr, lambda = lambda,
                  first.layer.width = first.layer.width,
                  hidden.layer.width = hidden.layer.width,
                  batch.size = 32L,
                  test.portion = 0.1, python.path = "/usr/local/bin/python3"))

})
