nn_train <- function(x, y, niter = 10L, learning.rate = 0.01, lambda,
                     first.layer.width = ncol(x) * 10L,
                     hidden.layer.width = 100L,
                     batch.size = 128L,
                     test.portion = .1, python.path,
                     model = NULL,
                     verbose = FALSE) {

  if(is.null(python.path)) python.path <- "/usr/local/bin/python3"
  reticulate::use_python(python.path)

  check_python_modules()

  # reticulate::source_python("Python/NN.py")
  # reticulate::source_python(nn_fun_path)
  NN <- reticulate::import_from_path(module = "NN", path = nn_fun_path)

  stopifnot(test.portion >= 0 & test.portion < 1)
  stopifnot(niter > 0)
  run.test <- test.portion > 0

  batch.size <- as.integer(batch.size)

  n <- nrow(x)
  d <- ncol(x)

  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(y)) y <- as.matrix(y)

  x_orig <- x

  if(run.test) {
    test.idx <- sample.int(n,floor(n*test.portion))
    xt <- x[test.idx,,drop=FALSE]
    x <- x[-test.idx,,drop=FALSE]

    yt <- y[test.idx,, drop=FALSE]
    y <- y[-test.idx,,drop=FALSE]
  }

  # make features as tensors
  train_features = torch$FloatTensor(x)
  train_target = torch$FloatTensor(y)

  # change to data laoders
  train_data = torch$utils$data$TensorDataset(train_features, train_target)
  train_loader = torch$utils$data$DataLoader(train_data, batch_size=batch.size, shuffle=TRUE)

  test_loader <- NULL

  if(run.test) {
    test_features = torch$FloatTensor(xt)
    test_target = torch$FloatTensor(yt)

    test_data = torch$utils$data$TensorDataset(test_features, test_target)
    test_loader = torch$utils$data$DataLoader(test_data, batch_size=batch.size, shuffle=FALSE)
  }


  if(is.null(model)) {
    model = NN$jit$jit_fun(as.integer(d), as.integer(first.layer.width), as.integer(hidden.layer.width))
  }

  #args for train
  #def train(model, train_iter, test_iter, D, nLayer, num_epoch, lr, test=False):
  # module is neural

  result <- NN$model$train(model = model, train_iter = train_loader,
                  test_iter = test_loader,
                  num_epoch = as.integer(niter),
                  lr = as.double(learning.rate),
                  lam = as.double(lambda), test = run.test,
                  verbose = verbose)

  # result[[3]]$eval()
  # yhat <- result[[3]]$predict(train_features)$data$numpy()

  yhat <- plogis(model$predict(torch$FloatTensor(x_orig))$data$numpy())

  # all.equal(yhat, yhat2)

  output <- list(yhat = yhat, model = model)
  return(output)

}
