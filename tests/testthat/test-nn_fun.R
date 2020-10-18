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

  #check names
  testthat::expect_equal(names(res), c("yhat","model"))

  #check model
  testthat::expect_true(!is.null(res$model))

  #check yhat
  testthat::expect_true(!is.null(res$yhat))

  #check can get derivatives
  xt <- torch$FloatTensor(matrix(rnorm(n*d), n, d)[1,,drop=FALSE])
  xtv <- torch$autograd$Variable(xt, requires_grad = TRUE)
  xtv$retain_grad()
  holder <- res$model$predict(xtv)
  holder$backward()
  testthat::expect_true(all(xtv$grad$data$numpy()!=0))
  # xtv$grad$zero_()
  # xtv$grad

  # check derivative iteration works
  xt <- torch$FloatTensor(matrix(rnorm(50*d), 50, d))
  derivative <- matrix(NA, 50, d)
  for(i in 0:(xt$shape[0] - 1)) {
    xtv <- torch$autograd$Variable(xt[i], requires_grad = TRUE)
    xtv$retain_grad()
    holder <- res$model$predict(xtv)
    holder$backward()
    testthat::expect_true(all(xtv$grad$data$numpy()!=0))
    derivative[i + 1,] <- xtv$grad$data$numpy()
    xtv$grad$zero_()
    testthat::expect_true(all(xtv$grad$data$numpy()==0))
  }
  testthat::expect_true(all(derivative != 0))
})
