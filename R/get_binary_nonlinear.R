
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
    orders <- order(x)
    s_x <- sort(c(0,x))
    zero_idx <- which(s_x == 0)
    sqrt_diffs <- sqrt(diff(s_x))
    increments <- rnorm(n, 0, sqrt_diffs)
    values <- c(cumsum(increments[(zero_idx-1):1]), cumsum(increments[zero_idx:n]))

    return(values[orders])
  }
  gp <- function(x) {
    sq_exp <- function(x) {
      return(exp(-0.5*as.matrix(dist(x)^2)))
    }
    n <- length(x)
    z <- rnorm(n)
    cov <- 1*sq_exp(x)
    svds <- svd(cov)
    U <- svds$u
    V <- svds$v
    d_half <- sqrt(svds$d)
    sq_cov <- U %*% diag(d_half) %*% t(V)
    return(sq_cov %*% z)
  }
  modified_friedman <- function(x){
    # f(x) = 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3]-0.5)^2  10 * x[,4] + 5 * x[,5]
    x3sq <- x[,3]^2
    costerm <- cos(2 * x3sq + pi/2)
    x3term <- ifelse(2 * x3sq > 6/4 * pi, 1,  costerm)
    return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) + 2 * x3term  + tanh(-x[,4]) + 0.25 * x[,5])
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) - 2 * x[,3]^2 + x[,4] + 0.5 * x[,5])
    #sin(pi/8* X[,2]*X[,3] + pi/2)* cos(2*X[,3]^2 + pi/2)
  }
  cart <- function(x) {
    require(rpart)
    require(dbarts)
    temp_eta <- modified_friedman(x)
    # theta <- rnorm(5)
    # temp_eta <- apply(x[,1:5],2,gp) %*% abs(theta)/sum(abs(theta))
    # maxcol <- min(ncol(x), 15)
    temp_prob <- plogis(temp_eta )
    temp_y <- rbinom(nrow(x), 1, temp_prob)
    df <- data.frame(y = temp_y, x1=x[,1], x2=x[,2], x3 = x[,3], x4 = x[,4], x5=x[,5])
    # rpcntrl <- rpart.control(minsplit=10, minbucket=5, cp=0.008)
    # cartFit <- rpart(y ~ x1 + x2 + x3 + x4 + x5, method="class", data=df, control=rpcntrl)
    cartFit <- bart(x,temp_y)
    phi <- cartFit$yhat.train
    prob <- apply(pnorm(phi),2,median)
    eta <- qlogis(prob)
    return(list(model = cartFit, eta=eta))
    plot(cartFit)
    sum(is.infinite(eta))
    hist(eta)
    hist(plogis(eta))
  }
  multiply.cols <- function(x) {
    return(apply(x,1,prod) )
  }
  data_gen_functions <- list(brownian=brownian, modified_friedman = modified_friedman)
  rdata <- function(n, x, theta, ...) {

    dots <- list(...)
    method <- dots$method
    if(is.null(method)) method <- "cart"
    method <- match.arg(method, c("brownian","modified.friedman","cart","gp"))

    p <- ncol(x)-1
    logit.p <- NA
    if (method == "brownian") {
      data_gen_functions <- brownian
      num.comb <- length(theta) - ncol(x)
      potential.combinations <- combn(rep(1:p,2), 2)
      combinations <- potential.combinations[,sample(1:ncol(potential.combinations), num.comb)]
      if (num.comb > 0){
        design.matrix <- cbind(x[,-1],
                               sapply(1:ncol(combinations),
                                      function(i) multiply.cols(x[,1+combinations[,i]])))
      } else {
        design.matrix <- x[,-1]
      }
      browned <- apply(design.matrix,2,brownian)
      logit.p <-  cbind(1, browned) %*% theta
    } else if (method == "modified.friedman") {
      data_gen_functions <- modified_friedman
      stopifnot(ncol(x)>=6)
      logit.p <- modified_friedman(x[,-1])
      theta <- c(1,2,1,0.5)
    } else if (method == "cart") {
      stopifnot(ncol(x)>=6)
      output <- cart(x[,-1])
      data_gen_functions <- output$model
      logit.p <- output$eta
    } else if (method == "gp") {
      output <- apply(x[,-1],2,gp)
      theta <- theta[1:(p+1)]
      theta <- abs(theta)/sqrt(sum(theta^2))
      logit.p <- cbind(1,output) %*% theta
      data_gen_functions <- gp
    }

    probs <- plogis(logit.p)
    return(list(Y = rbinom(n, 1, probs), probs = probs, theta = theta, data_gen_functions=data_gen_functions, method=method))
  }

  #### Posterior on Coefficients ####
  rpost_coef <- function(x){NULL}
  rpost_sigma <- function(x){NULL}

  rpost <- function(n.samp, x, y, hyperparameters,...) { # implements either BART or BNN
    dots <- list(...)
    theta <- NULL
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


      dnn <- stan_model(stan_dir)
      stan_dat <- list(N = n, P=p, L = L, X=x, Y=y, m0 = m0, nodes=nodes)
      vb.out <- vb(dnn,stan_dat, algorithm=algorithm, iter=iter,
                   grad_samples=grad_samples, eval_elbo=eval_elbo, eta=eta,
                   adapt_engaged=adapt_engaged, init=init, tol_rel_obj = tol,
                   output_samples = n.samp)
      pars <- extract(vb.out, pars=c("eta","prob"))
      eta <- pars$eta
      prob <- pars$prob
    } else if (dots$method == "bart") {
      require(dbarts)
      nskip <- dots$nskip
      if(is.null(nskip)) nskip <- n.samp
      k <- dots$k
      power <- dots$power
      base <- dots$base
      ntree <- dots$ntree
      keepevery <- dots$keepevery
      test <- dots$test
      if(is.null(k)) k <- 2.0
      if(is.null(power)) power <- 2.0
      if(is.null(base)) base <- 0.95
      if(is.null(ntree)) ntree <- 200
      if(is.null(keepevery)) keepevery <- 1
      if (is.null(test)) test <- NULL

      bartFit <- bart(x, y, ndpost=n.samp*keepevery, nskip=nskip,
                      k=k, base=base, power=power, ntree=ntree, keepevery = keepevery,
                      printevery = n.samp*keepevery/10, keeptrees=FALSE, x.test = test)
      phi <- bartFit$yhat.train
      prob <- pnorm(phi)
      eta <- qlogis(prob)
      model <- bartFit

    } else if (dots$method == "logistic"){
      require(rstan)

      if (all(x[,1]==1)) x <- x[,-1]


      p <- ncol(x)
      n <- nrow(x)

      stan_dir <- dots$stan_dir
      m0 <- dots$m0
      scale_intercept <- dots$scale_intercept
      chains <- dots$chains

      if(is.null(m0)) m0 <- round(0.1 * p)
      if(is.null(scale_intercept)) scale_intercept <- 2.5
      if(is.null(chains)) chains <-4

      x_c <- scale(x, scale=FALSE)
      qrdecomp <- qr(x_c)
      Q_x <- qr.Q(qrdecomp) * sqrt(n-1)
      R_inv <- solve(qr.R(qrdecomp)/sqrt(n-1))

      stan_dat <- list(N = n,
                       P = p,
                       Y = y,
                       X = x,
                       mean_x = attr(x_c, "scaled:center"),
                       Q_x = Q_x,
                       R_inv_x = R_inv,
                       m0 = m0,
                       scale_intercept = scale_intercept)

      stanModel <- stan_model(stan_dir)
      stanFit <- sampling(stanModel, data=stan_dat, iter=n.samp*2,
                          warmup =n.samp* (2-1/chains) , chains=chains, pars = c("eta","theta","prob", "beta", "intercept"))
      samples <- extract(stanFit, pars= c("eta","theta","prob"))
      eta <- samples$eta
      prob <- samples$prob
      theta <- samples$theta
      model <- stanFit

    } else if (dots$method == "vb.logistic"){
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

      if(is.null(m0)) m0 <- round(0.1 * p)
      if(is.null(scale_intercept)) scale_intercept <- 2.5

      x_c <- scale(x, scale=FALSE)
      qrdecomp <- qr(x_c)
      Q_x <- qr.Q(qrdecomp) * sqrt(n-1)
      R_inv <- solve(qr.R(qrdecomp)/sqrt(n-1))

      stan_dat <- list(N = n,
                       P = p,
                       Y = y,
                       X = x,
                       mean_x = attr(x_c, "scaled:center"),
                       Q_x = Q_x,
                       R_inv_x = R_inv,
                       m0 = m0,
                       scale_intercept = scale_intercept)

      stanModel <- stan_model(stan_dir)
      stanFit <- vb(stanModel, data=stan_dat, algorithm=algorithm, iter=iter,
                    grad_samples=grad_samples, eval_elbo=eval_elbo, eta=eta,
                    adapt_engaged=adapt_engaged, init=init, tol_rel_obj = tol,
                    output_samples = n.samp, pars=c("eta","prob","theta"))
      samples <- extract(stanFit, pars=c("eta","prob","theta"))
      eta <- samples$eta
      prob <- samples$prob
      theta <- samples$theta
      model <- stanFit

    } else if (dots$method == "gamm") {
      require(rstanarm)
      if(all(x[,1] == 1)) x <- x[,-1]
      p <- ncol(x)
      n <- nrow(x)

      # stan_dir <- dots$stan_dir
      # m0 <- dots$m0
      # scale_intercept <- dots$scale_intercept
      chains <- dots$chains

      stan_dat <- data.frame(Y = Y)
      namesX <- paste0("X",1:p)
      for(i in 1:p) stan_dat[[namesX[i]]] <- x[,i]
      form <- as.formula(paste("Y ~ ", paste0("s(",namesX,")" , collapse= " + ")))

      stanFit <- stan_gamm4(form, data = stan_dat, family = binomial(),
                 chains = chains, iter = 2*n.samp, algorithm="sampling",
                 warmup = n.samp* (2-1/chains) )
      eta <- posterior_linpred(stanFit)
      if(nrow(eta) > n.samp) eta <- eta[1:n.samp,]
      prob <- plogis(eta)
      params <- as.matrix(stanFit$stanfit)
      theta <- params[,seq_len(ncol(stanFit$x))]
      model <- stanFit
    }

    return(list(theta=theta,prob=prob, eta=eta, model=model))

  }

  #### Output ####
  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = data_gen_functions,
              rparam = rparam))
}
