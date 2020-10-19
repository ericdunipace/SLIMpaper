
get_binary_nonlinear_model <- function() {
  #### Prior Function ####
  rprior_coef <- function(n, mu, sigma){
    return(rmvnorm(n, mu, sigma))
  }

  rprior_sigma <- function(n, alpha, beta){
    return(NULL)
  }

  rprior<- function(n, hyperparameters) {
    stopifnot(is.list(hyperparameters))
    stopifnot(length(hyperparameters) >= 2)
    stopifnot(all(c("mu","sigma") %in% names(hyperparameters) ))
    stopifnot(length(hyperparameters$mu) == ncol(as.matrix(hyperparameters$sigma)))
    return(list(theta = rprior_coef(n, hyperparameters$mu, hyperparameters$sigma),
                sigma2 = NULL))
  }

  #### X Data ####
  rX <- gen_x()$rX

  #### Parameters ####
  rparam <- get_param()$binary
  #### Y Data ####
  brownian <- function(x) {
    n <- length(x)
    original.orders <- rank(x)
    s_x <- sort(c(0,x))
    zero_idx <- which(s_x == 0)
    sqrt_diffs <- sqrt(diff(s_x))
    increments <- rnorm(n, 0, sqrt_diffs)
    values <- c(cumsum(increments[(zero_idx-1):1]), cumsum(increments[zero_idx:n]))

    return(values[original.orders])
  }
  gp <- function(x) {
    sq_exp <- function(x) {
      L <- diag(0.5, ncol(x), ncol(x))
      sigma <- 1
      x_scale <- x %*% chol(L)
      distances <- as.matrix(dist(x_scale, diag = TRUE))
      Sigma <- sigma^2 * exp(-0.5*distances)
      return(Sigma)
    }
    n <- nrow(x)
    d <- ncol(x)
    z <- rnorm(n)
    cov <- 1*sq_exp(x)
    svds <- svd(cov)
    U <- svds$u
    V <- svds$v
    d_half <- sqrt(svds$d)
    sq_cov <- U %*% diag(d_half) %*% t(V)
    return(sq_cov %*% z)
  }
  mf_data <- function(x){
    if(all(x[,1]==0)) x <- x[,-1]
    if(all(x[,1]==1)) x <- x[,-1]

    x6sq <- x[,6]^2
    costerm <- cos(2 * x6sq + pi/2)
    x6term <- ifelse(2 * x6sq > 3/4 * pi, 1,  costerm)
    x1_2 <- sin(pi/8 * x[,1] * x[,2] + pi/2)
    x11_15 <- tanh(-pi/8 * x[,11] * x[,15]^2)
    x7_exp <- exp(-x[,7]^2/5)
    x7_3 <- (x[,7])^3 * x7_exp
    x7_2 <- (x[,7])^2 * x7_exp
    x7_1 <- x[,7]  * x7_exp

    return(cbind(x1_2, x6term, x11_15, x7_3, x7_2, x7_1))

  }
  modified_friedman <- function(x){
    # f(x) = 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3]-0.5)^2  10 * x[,4] + 5 * x[,5]
    x6sq <- x[,6, drop = FALSE]^2
    costerm <- cos(2 * x6sq + pi/2)
    x6term <- ifelse(2 * x6sq > 3/4 * pi, 1,  costerm)
    # x6term <- costerm
    return(1 * sin(pi/8 * x[,1, drop = FALSE] * x[,2, drop = FALSE] + pi/2) +
             2 * x6term  + tanh(-pi/8 * x[,11, drop = FALSE] * x[,15, drop = FALSE]^2) +
             (2/5 * (x[,7, drop = FALSE])^3 - 1/8 * (x[,7, drop = FALSE])^2 - 2 * x[,7, drop = FALSE])*exp(-x[,7, drop = FALSE]^2/5) )
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) + 2 * x6sq  + sin(-pi/8 * x[,11] * x[,15]^2) + 1/20 * x[,7]^3 - 1/24 * x[,7]^2 - 1 * x[,7])
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) + 2 * x6sq  + tanh(-x[,11] * x[,15]^2) + 0.25 * x[,7])
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) - 2 * x[,3]^2 + x[,4] + 0.5 * x[,5])
  }
  cart <- function(x) {
    require(rpart)
    # require(dbarts)
    temp_eta <- modified_friedman(x)
    # theta <- rnorm(5)
    # temp_eta <- apply(x[,1:5],2,gp) %*% abs(theta)/sum(abs(theta))
    # maxcol <- min(ncol(x), 15)
    temp_prob <- plogis(temp_eta )
    temp_y <- rbinom(nrow(x), 1, temp_prob)
    df <- data.frame(y = temp_y, x1=x[,1], x2=x[,2], x3 = x[,3], x4 = x[,4], x5=x[,5])
    # rpcntrl <- rpart.control(minsplit=10, minbucket=5, cp=0.008)
    # cartFit <- rpart(y ~ x1 + x2 + x3 + x4 + x5, method="class", data=df, control=rpcntrl)
    cartFit <- dbarts::bart(x,temp_y)
    phi <- cartFit$yhat.train
    prob <- apply(pnorm(phi),2,median)
    eta <- qlogis(prob)
    return(list(model = cartFit, eta=eta))
    plot(cartFit)
    sum(is.infinite(eta))
    hist(eta)
    hist(plogis(eta))
  }

  model_matrix <- function(X, ...) {
    # if(all(X[,1] == 1)) X <- X[,-1,drop=FALSE]
    # p <- ncol(X)
    # form <- formula(paste0("~. + ", paste0("I(X",1:p,"^2)", collapse=" + ")))
    #
    # X <- model.matrix(object=form, data = data.frame(X))
    return(X)
  }


  multiply.cols <- function(x) {
    return(apply(x,1,prod) )
  }
  data_gen_functions <- list(brownian=brownian,
                             modified_friedman = modified_friedman,
                             gp = gp)
  rdata <- function(n, x, param, ...) {

    dots <- list(...)
    method <- dots$method
    if(is.null(method)) method <- "cart"
    method <- match.arg(method, c("brownian","modified.friedman","cart","gp","glm", "polynomial","random.interaction"))


    is.intercept <- all(x[,1] == 1)
    if(is.intercept) {
      x <- x[,-1, drop=FALSE]
    }

    p <- ncol(x)
    logit.p <- NA
    extra <- NULL
    if (method == "brownian") {
      data_gen_functions <- brownian
      num.comb <- length(param) - ncol(x)
      if (num.comb > 0){
        potential.combinations <- combn(rep(1:p,2), 2)
        combinations <- potential.combinations[,sample(1:ncol(potential.combinations), num.comb)]
        design.matrix <- cbind(x,
                               sapply(1:ncol(combinations),
                                      function(i) multiply.cols(x[,combinations[,i]])))
      } else {
        design.matrix <- x[,1:(length(param)-1)]
      }
      browned <- apply(design.matrix,2,brownian)
      logit.p <-  cbind(1, browned) %*% param
      theta <- param
    }
    else if (method == "modified.friedman") {
      # data_gen_functions <- modified_friedman
      data_gen_functions <- mf_data
      stopifnot(ncol(x)>=6)
      logit.p <- modified_friedman(x)
      theta <- c(1, 2, 1, 2/5, -1/8, -2)
    }
    else if (method == "cart") {
      stopifnot(ncol(x)>=6)
      output <- cart(x)
      data_gen_functions <- output$model
      logit.p <- output$eta
    }
    else if (method == "gp") {
      logit.p <- gp(x)
      # theta <- param[1:(p+1)]
      # theta <- abs(theta)/sqrt(sum(theta^2))
      # logit.p <- cbind(1,output) %*% theta
      theta <- list(L = diag(0.5,p, p), sigma = 1)
      data_gen_functions <- gp
    }
    else if (method == "glm") {
      corr.x <- dots$corr
      scale <- dots$scale
      if(is.null(corr.x)) corr.x <- 0
      if(is.null(scale)) scale <- TRUE
      theta <- dots$theta

      if(is.intercept) {
        intercept <- theta[1]
        theta <- theta[-1]
      } else {
        intercept <- 0
      }

      if(p > length(theta)) {
        x <- x[, 1: length(theta), drop=FALSE]
        warning("Ncol X > length(theta). Only using first length(theta) columns of X.")
      }

      corr.mat <- corr_mat_construct(corr.x, p)
      diag(corr.mat) <- 1
      theta_norm <- c(t(theta) %*% corr.mat %*% theta)
      theta_scaled <- if(scale) {
        theta/sqrt(theta_norm)
      } else {
        theta
      }#* sqrt(sigma2)
      logit.p <- x %*% theta_scaled + intercept
    }
    else if (method == "random.interaction") {
      # param <- dots$param
      if(is.null(param) | length(param) != (20 +1)) stop("must specifiy 21 parameters")

      if(all(x[,1] == 1)) x <- x[,-1]
      p_star <- min(p, 20)
      x <- data.frame(x[,1:p_star])
      # pwrs <- c(rep(1, p_star + 1), 2*rpois(p_star, 0.5)+1)

      # formula <- formula(paste0("~ .^3 "))#, paste0("I(",colnames(x[,1:p_star]),"^2)", collapse=" + ")))
      # combinations <- model.matrix(formula, data=x)
      # select <- sample.int(ncol(combinations) - p_star, p_star) + p_star + 1
      # coef.names <- strsplit(colnames(combinations)[select], ":")
      # cmb <- matrix(NA, ncol=2*p_star + 1 , nrow = nrow(x))
      # cmb[,1:(p_star+1)] <- as.matrix(combinations[,1:(p_star + 1)])
      # colnames(cmb) <- c(colnames(combinations)[1:(p_star+1)], rep(NA, p_star))
      # sample.odd <- function(){
      #   rpois(1, 0.5) * 2 + 1
      # }
      # sample.even <- function() {
      #   rpois(1, 0.5) * 2 + 2
      # }
      # for(j in 1:p_star) {
      #   jj <- j + p_star + 1
      #   cur.coef <- sample(coef.names[[j]],length(coef.names[[j]]))
      #   ncoef <- length(cur.coef)
      #   pwrs <- rep(NA, ncoef)
      #   tmp <- rep(1, nrow(x))
      #   nm <- NULL
      #   for(i in 1:ncoef) {
      #     pwrs <- ifelse((i %% 2)==0, sample.odd(), sample.even())
      #     nm <- c(nm, paste0(cur.coef[i],"^",pwrs))
      #     tmp <- tmp*x[,cur.coef[i]] ^ pwrs
      #   }
      #   cmb[,jj] <- tmp
      #   colnames(cmb)[jj] <- paste0("(",nm,")", collapse=":")
      # }
      # cmb <- combinations[,c(1:(p_star+1),  select)] * pwrs
      # colnames(cmb) <- paste0("(", colnames(cmb), ")", "^",pwrs)
      # param <- c(param, rep(NA, p_star))
      # for(i in (p_star + 1):ncol(cmb)){
      #   param[i] <- runif(1, -0.05, 0.05) * 2/max(abs(cmb[,i]))
      # }
      formula <- formula(paste0("~ .^3 "))#, paste0("I(",colnames(x[,1:p_star]),"^2)", collapse=" + ")))
      combinations <- model.matrix(formula, data=x)
      prob.sel <- ((ncol(combinations)-1):1) / sum( (ncol(combinations) - 1):1 )
      select <- c(1,sample.int(ncol(combinations)-1, p_star, prob = prob.sel) + 1)
      cmb <- combinations[,select]
      colnames(cmb) <- colnames(combinations)[select]
      theta <- param
      extra <- list(Xgen = cmb, param = param)
      logit.p <- cmb %*% param

    }
    else if (method == "polynomial") {
      # param <- dots$param
      if (is.null(param) | length(param) < 6) stop("must specifiy at least  6parameters")

      if(all(x[,1] == 1)) x <- x[,-1]
      p_star <- min(p, 20)
      x <- data.frame(x[,1:p_star])

      formula <- formula("~ I(tanh(X1^3 * X20)) + I( pnorm(X8 * X15 * X20)*2 - 1) + X2  + I(log(1+X2^2*X12^2)-2) +  I( exp(X6) - 1.64) ")
      mm <- model.matrix(formula, data=x)
      intercept_adjust <- c(1,-1, 0, -2, -1.64)
      theta <- param[1:min(p_star,6)] * c(1,2,1,1,1,1)
      # -0.7 - 1 * 0.14  - 2 * 0.12 - 1.64 * 0.18
      theta[1] <- round(param[1:min(p_star,6)] %*% c(1,0, -1, 0, -2, -1.64)[1:min(p_star,6)], digits=3)
      extra <- list(Xgen = mm, param = theta)
      logit.p <- mm %*% theta

    }
    else {
      stop("Method not found!")
    }

    probs <- plogis(logit.p)
    return(list(Y = rbinom(n, 1, probs), mu = probs, eta = logit.p,
                link = binomial()$linkfun, invlink = binomial()$linkinv, param = theta,
                data_gen_functions=data_gen_functions,
                model_matrix = model_matrix,
                method=method, extra = extra))
  }

  #### Posterior on Coefficients ####
  rpost_coef <- function(x){NULL}
  rpost_sigma <- function(x){NULL}
  rnn <- function(n.samp, x, y, hyperparameters,...) {

    n <- nrow(x)

    x <- as.matrix(x)
    y <- as.matrix(y)

    dots <- list(...)

    if(is.null(dots$test.portion)) {
      test.portion <- 0.1
    } else {
      test.portion <- dots$test.portion
    }

    if(is.null(dots$niter)) {
      niter <- 20
    } else {
      niter <- dots$niter
    }

    if(is.null(hyperparameters$learning.rate)) {
      learning.rate <- 1e-2
    } else {
      learning.rate <- hyperparameters$learning.rate
    }

    if(is.null(hyperparameters$lambda)) {
      lambda <- 1
    } else {
      lambda <- hyperparameters$lambda
    }

    if(is.null(hyperparameters$batch.size)) {
      batch.size <- 128
    } else {
      batch.size <- hyperparameters$batch.size
    }

    if(is.null(hyperparameters$first.layer.width)) {
      first.layer.width <- ncol(x) * 10
    } else {
      first.layer.width <- hyperparameters$first.layer.width
    }
    if(is.null(hyperparameters$hidden.layer.width)) {
      hidden.layer.width <- as.integer(2/3 * ncol(x) + 1)
    } else {
      hidden.layer.width <- hyperparameters$hidden.layer.width
    }

    python.path <- dots$python.path

    if(is.null(python.path) | python.path == "") {
      python.path <- NULL
    }

    if(is.null(dots$verbose)) {
      verbose <- FALSE
    } else {
      verbose <- dots$verbose
    }

    res <- nn_train(x=x, y=y, niter = niter, learning.rate = learning.rate,
                    lambda = lambda,
                    test.portion = test.portion,
                    first.layer.width = first.layer.width,
                    hidden.layer.width = hidden.layer.width,
                    batch.size = batch.size,
                    python.path = python.path,
                    model = NULL,
                    verbose = verbose)
    xt <- dots[["X.test"]]
    run.test <- !is.null(xt)
    if (run.test) {
      if (all(xt[,1] == 1)) xt <- xt[,-1,drop=FALSE]
      xtt <- torch$FloatTensor(xt)
    }
    yhat.test <- plogis(temp$model$predict(xtt)$data$numpy())
    boots <- lapply(1:n.samp, function(i) {
      boot.idx <- sample.int(n,n,replace=TRUE)
      temp <- nn_train(x=x[boot.idx, , drop = FALSE], y=y[boot.idx, , drop=FALSE],
                       niter = niter, learning.rate = learning.rate,
                       lambda = lambda,
                       test.portion = test.portion,
                       first.layer.width = first.layer.width,
                       hidden.layer.width = hidden.layer.width,
                       batch.size = batch.size,
                       python.path = python.path,model = NULL)
      yhat <- temp$yhat
      yhat.test <- NULL
      derivative.x <- NULL
      if (run.test) {
        yhat.test <- plogis(temp$model$predict(xtt)$data$numpy())
        derivative.x <- matrix(NA, nrow(xt), ncol(xt))
        for(i in 1:nrow(xt)) {
          xtv <- torch$autograd$Variable(xtt[i-1], requires_grad = TRUE)
          xtv$retain_grad()
          temp.pred <- temp$model$predict(xtv)
          temp.pred$backward()
          derivative.x[i,] <- xtv$grad$data$numpy()
        }
      }
      return(list(mu = yhat, mu.test = yhat.test, derivative.x = derivative.x))
    })

    mu <- sapply(boots, function(b) b$mu)
    mu.test <- sapply(boots, function(b) b$mu.test)
    derivatives <- lapply(boots, function(b) b$derivative.x)

    return(list(mu = mu, mu.test = mu.test, derivatives = derivatives, model = res$model,
                yhat.model = list(train = res$yhat,
                                  test = yhat.test))
           )


  }

  rpost <- function(n.samp, x, y, hyperparameters,...) { # implements either BART or BNN
    dots <- list(...)
    theta <- NULL
    test <- list(eta  = NULL, mu = NULL)
    if (dots$method == "bnn"){
      require(rstan)

      stan_dir <- dots$stan_dir
      tol <- dots$tol
      iter <- dots$iter
      grad_samples <- dots$grad_samples
      eta <- dots$eta
      init <- dots$init
      eval_elbo <- dots$eval_elbo
      algorithm <- dots$algorithm

      if (is.null(tol)) tol <- 0.1
      if (is.null(iter)) iter <- 1e4
      if (is.null(grad_samples)) grad_samples <- 1
      if (is.null(eval_elbo)) eval_elbo <- 1e2
      if (is.null(eta))  {
        adapt_engaged <- TRUE
      } else {
        adapt_engaged <- FALSE
      }
      if (is.null(init)) init <- 0
      if (all(x[,1]==1)) x <- x[,-1]
      if (is.null(algorithm)) algorithm <- "meanfield"


      p <- ncol(x)
      n <- nrow(x)
      m0 <- hyperparameters$m0
      nodes <- hyperparameters$nodes
      L <- hyperparameters$L
      if(m0 >= p) {
        m0 <- p - 1
        warning("Adjusting m0 value. Must be less than number of predictors")
      }

      rstan::rstan_options(auto_write = TRUE)
      dnn <- rstan::stan_model(stan_dir)
      stan_dat <- list(N = n, P=p, L = L, X=x, Y=y, m0 = m0, nodes=nodes)
      vb.out <- rstan::vb(dnn,stan_dat, algorithm=algorithm, iter=iter,
                   grad_samples=grad_samples, eval_elbo=eval_elbo, eta=eta,
                   adapt_engaged=adapt_engaged, init=init, tol_rel_obj = tol,
                   output_samples = n.samp)
      pars <- extract(vb.out, pars=c("eta","prob"))
      eta <- pars$eta
      prob <- pars$prob

    }
    else if (dots$method == "bart") {
      # require(dbarts)
      if(all(x[,1] == 1)) x <- x[,-1]

      nskip <- dots[["nskip"]]
      ntree <- dots[["ntree"]]
      keepevery <- dots[["keepevery"]]
      X.test <- dots[["X.test"]]
      numcut <- dots[["numcut"]]
      chains <- dots[["chains"]]


      k <- hyperparameters[["k"]]
      power <- hyperparameters[["power"]]
      base <- hyperparameters[["base"]]

      if ( is.null(ntree) ) ntree <- 200
      if ( is.null(keepevery) ) keepevery <- 1
      # if ( is.null(test) ) test <- NULL
      if (is.null(chains)) chains <- 1
      if ( is.null(nskip) ) nskip <- n.samp * keepevery
      if ( is.null(numcut) ) numcut <- 100

      # if ( !("k" %in% names(hyperparameters)) ) k <- 2.0
      if ( is.null(k) ) k <- "chi(degreesOfFreedom = 1.25, scale = Inf)"
      if ( !is.numeric(k) ) if ( !is.character(k) ) stop( "k must be a number > 0 or a string of the form 'chi(degreesOfFreedom, scale )'")
      if ( is.null(power) ) power <- 2.0 # function default
      if ( is.null(base) ) base <- 0.95 # function default

      ndpost <- n.samp * keepevery / chains

      bartFit <- dbarts::bart(x.train = x, y.train = y, ndpost = ndpost, nskip = nskip,
                      k = k, base=base, power=power,
                      numcut = numcut, ntree=ntree, keepevery = keepevery,
                      printevery = n.samp*keepevery/10, keeptrees=FALSE, x.test = X.test,
                      nchain = chains, nthread = min(chains,parallel::detectCores()-1))
      phi <- bartFit$yhat.train
      mu <- pnorm(phi)
      eta <- qlogis(mu)
      model <- bartFit
      if(!is.null(X.test)){
        test$mu  <- pnorm(bartFit$yhat.test)
        test$eta <- qlogis(test$mu )
      }


    }
    else if (dots$method == "logistic"){
      # require(rstan)

      if (all(x[,1]==1)) x <- x[,-1]


      p <- ncol(x)
      n <- nrow(x)

      stan_dir <- dots$stan_dir
      m0 <- hyperparameters$m0
      scale_intercept <- hyperparameters$scale_intercept
      chains <- dots$chains
      X.test <- dots[["X.test"]]

      if(is.null(m0)) m0 <- round(0.1 * p)
      if(is.null(scale_intercept)) scale_intercept <- 2.5
      if(is.null(chains)) chains <- 4

      if(m0 >= p) {
        m0 <- p - 1
        warning("Adjusting m0 value. Must be less than number of predictors")
      }

      x_c <- scale(x, scale=FALSE)
      # {
      #   qrdecomp <- qr(x_c)
      #   Q_x <- qr.Q(qrdecomp)[,1:p] * sqrt(n-1)
      #   R_inv <- solve(qr.R(qrdecomp)/sqrt(n-1))
      # }

      stan_dat <- list(N = n,
                       P = p,
                       Y = y,
                       X = x,
                       mean_x = attr(x_c, "scaled:center"),
                       m0 = m0,
                       scale_intercept = scale_intercept)

      warmup <- max(n.samp*(2-1/chains), 1000)
      iter <- warmup + ceiling(n.samp/chains)

      rstan::rstan_options(auto_write = TRUE)
      stanModel <- rstan::stan_model(stan_dir)
      stanFit <- rstan::sampling(stanModel, data = stan_dat, iter = iter,
                          warmup = warmup, chains = chains, pars = c("eta","theta","prob", "beta", "intercept"))
      samples <- rstan::extract(stanFit, pars= c("eta","theta","prob"))
      eta <- t(samples$eta)
      mu <- t(samples$prob)
      theta <- t(samples$theta)
      model <- stanFit
      if(!is.null(X.test)){
        test$eta <- tcrossprod(X.test, samples$theta)
        test$mu <- plogis(test$eta )
      }


    }
    else if (dots$method == "vb.logistic"){
      require(rstan)

      if (all(x[,1]==1)) x <- x[,-1]

      stan_dir <- dots$stan_dir
      tol <- dots$tol
      iter <- dots$iter
      grad_samples <- dots$grad_samples
      eta <- dots$eta
      init <- dots$init
      eval_elbo <- dots$eval_elbo
      algorithm <- dots$algorithm

      if (is.null(tol)) tol <- 0.1
      if (is.null(iter)) iter <- 1e4
      if (is.null(grad_samples)) grad_samples <- 1
      if (is.null(eval_elbo)) eval_elbo <- 1e2
      if (is.null(eta))  {
        adapt_engaged <- TRUE
      } else {
        adapt_engaged <- FALSE
      }
      if (is.null(init)) init <- 0
      if (is.null(algorithm)) algorithm <- "meanfield"

      p <- ncol(x)
      n <- nrow(x)

      stan_dir <- dots$stan_dir
      m0 <- dots$m0
      scale_intercept <- dots$scale_intercept
      X.test <- dots[["X.test"]]

      if(is.null(m0)) m0 <- round(0.1 * p)
      if(is.null(scale_intercept)) scale_intercept <- 2.5

      x_c <- scale(x, scale=FALSE)
      qrdecomp <- qr(x_c)
      Q_x <- qr.Q(qrdecomp) * sqrt(n-1)
      R_inv <- solve(qr.R(qrdecomp)/sqrt(n-1))

      if(m0 >= p) {
        m0 <- p - 1
        warning("Adjusting m0 value. Must be less than number of predictors")
      }

      stan_dat <- list(N = n,
                       P = p,
                       Y = y,
                       X = x,
                       mean_x = attr(x_c, "scaled:center"),
                       Q_x = Q_x,
                       R_inv_x = R_inv,
                       m0 = m0,
                       scale_intercept = scale_intercept)
      rstan::rstan_options(auto_write = TRUE)
      stanModel <- rstan::stan_model(stan_dir)
      stanFit <- rstan::vb(stanModel, data=stan_dat, algorithm=algorithm, iter=iter,
                    grad_samples=grad_samples, eval_elbo=eval_elbo, eta=eta,
                    adapt_engaged=adapt_engaged, init=init, tol_rel_obj = tol,
                    output_samples = n.samp, pars=c("eta","prob","theta"))
      samples <- rstan::extract(stanFit, pars=c("eta","prob","theta"))
      eta <- samples$eta
      mu <- samples$prob
      theta <- t(samples$theta)
      model <- stanFit
      # test$eta <- tcrossprod(X.test, samples$theta)
      # test$mu <- plogis(test$eta)

    }
    else if (dots$method == "gamm") {
      # require(rstanarm)
      if(all(x[,1] == 1)) x <- x[,-1]
      p <- ncol(x)
      n <- nrow(x)

      # stan_dir <- dots$stan_dir
      m0 <- dots$m0
      scale_intercept <- dots$scale_intercept
      L <- dots$L
      chains <- dots$chains
      X.test <- dots[["X.test"]]
      prior <- hs_plus()
      warmup <- max(n.samp*(2-1/chains), 1000)
      iter <- warmup + ceiling(n.samp/chains)

      if(m0 >= p) {
        m0 <- p - 1
        warning("Adjusting m0 value. Must be less than number of predictors")
      }

      stan_dat <- data.frame(Y = Y)
      namesX <- paste0("X",1:p)
      for(i in 1:p) stan_dat[[namesX[i]]] <- x[,i]
      form <- as.formula(paste("Y ~ ", paste0("s(",namesX,")" , collapse= " + ")))

      stanFit <- rstanarm::stan_gamm4(form, data = stan_dat, family = binomial(),
                 chains = chains, iter = iter, algorithm="sampling",
                 warmup = warmup, prior = prior )
      eta <- rstanarm::posterior_linpred(stanFit)
      if(nrow(eta) > n.samp) eta <- eta[1:n.samp,]
      mu <- plogis(eta)
      params <- as.matrix(stanFit$stanfit)
      theta <- t(params[,seq_len(ncol(stanFit$x))])
      model <- stanFit
      if(!is.null(X.test)){
        if(all(X.test[,1]==1)) X.test <- X.test[,-1]
        X.test.lp <- mgcv::predict.gam(stanFit$jam, newdata=data.frame(X.test), type="lpmatrix")
        test$eta <- X.test.lp %*% theta
        test$mu <- plogis(test$eta)
      }
    }
    else if (dots$method == "gp"){
      # require(rstan)

      if (all(x[,1]==1)) x <- x[,-1]


      p <- ncol(x)
      n <- nrow(x)

      stan_dir <- dots$stan_dir
      # m0 <- hyperparameters$m0
      # scale_intercept <- hyperparameters$scale_intercept
      chains <- dots$chains
      X.test <- dots[["X.test"]]

      # if(is.null(m0)) m0 <- round(0.1 * p)
      # if(is.null(scale_intercept)) scale_intercept <- 2.5
      if(is.null(chains)) chains <- 4

      x_c <- scale(x, scale=FALSE)
      n_test <- 0
      x_test <- array(0, dim=c(0,p))
      if(!is.null(X.test)){
        if (all(X.test[,1]==1)) X.test <- X.test[,-1]
        x_test <- scale(X.test, center = attr(x_c, "scaled:center"),
                        scale = attr(x_c, "scaled:scale"))
      }
      # {
      #   qrdecomp <- qr(x_c)
      #   Q_x <- qr.Q(qrdecomp)[,1:p] * sqrt(n-1)
      #   R_inv <- solve(qr.R(qrdecomp)/sqrt(n-1))
      # }

      stan_dat <- list(N = n,
                       N_test = n_test,
                       P = p,
                       y = y,
                       x = x_c,
                       x_test = x_test)

      warmup <- max(n.samp*(2-1/chains), 1000)
      iter <- warmup + ceiling(n.samp/chains)

      rstan::rstan_options(auto_write = TRUE)
      stanModel <- rstan::stan_model(stan_dir)
      stanFit <- rstan::sampling(stanModel, data = stan_dat, iter = iter,
                                 warmup = warmup, chains = chains, pars = c("pred_eta"))
      samples <- rstan::extract(stanFit, pars= c("pred_eta"))
      eta <- t(samples$pred_eta[,1:n])
      mu <- plogis(eta)
      theta <- lm.fit(x=x, y = eta)$coefficients
      model <- stanFit
      if(!is.null(X.test)){
        test$eta <- t(samples$pred_eta[,(n+1):(n+n_test),drop=FALSE])
        test$mu <- plogis(test$eta )
      }


    }
    else if (dots$method == "nn") {
      warning("Not a Bayesian method! Does a NN with bootstrapped data")

      if (all(x[,1]==1)) x <- x[,-1, drop = FALSE]

      output <- rnn(n.samp, x, y, hyperparameters,...)
      theta <- NULL
      mu <- as.matrix(output$mu)
      eta <- qlogis(mu)

      model <- function(x) {
        return(as.matrix(output$model$predict(torch$FloatTensor(x))$data$numpy()))
      }

      if (!is.null(dots[["X.test"]])) {
        test$mu <- output$mu.test
        test$eta <- qlogis(test$mu)
        test$derivatives <- output$derivatives
      }
    }

    return(list(theta=theta, mu=mu, eta=eta, model=model, test = test))
  }


  #### mean functions for variable importance ####
  mf.friedman <- function(x,theta) {
    dat <- mf_data(x)
    eta <- x %*% theta
    return(binomial$linkinv(eta))
  }

  mf.linpred <- function(x,theta) {
    return(x %*% theta)
  }

  mf.logit <- function(x, theta) {
    eta <- x %*% theta
    return(binomial$linkinv(eta))
  }

  sel.pred.fun <- function(method = "linpred") {

    mf <- switch(method, "modified.friedman" = mf.friedman,
                 "glm" = mf.logit,
                 "logistic" = mf.logit,
                 "linpred" = mf.linpred)

    return(mf)
  }

  #### Output ####
  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = data_gen_functions,
              rparam = rparam,
              sel.pred.fun = sel.pred.fun,
              link = binomial()$linkfun,
              invlink = binomial()$linkinv,
              mf.data = mf_data))
}
