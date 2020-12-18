#' Get the Gaussian linear data simulation
#'
#' @return list with slots
#' `rprior`
#' `rdata`
#' `rpost`
#' `X = list(rX = rX, corr = NULL)`, rX being a function to simulate the x data = function(n, corr, p)
#'  `data_gen_functions`
#'  `rparam`
#'  `sel.pred.fun `
#'  `link` with the logistic link function
#'  `invlink` with the logistic inverse link function
#' @export
get_normal_linear_model <- function() {
  #### Prior Function ####
  rprior_coef <- function(n, mu, sigma){
    eigs <- eigen(sigma)
    sigmahalf <- eigs$vectors %*% diag(sqrt(eigs$values)) %*% t(eigs$vectors)
    p <- length(mu)
    stopifnot(ncol(sigma) == length(mu))

    mu <- matrix(mu, nrow=p, ncol=n)
    Z <- matrix(rnorm(n*p),p,n)
    return(mu + sigmahalf %*% Z )
  }

  rprior_sigma <- function(n, alpha, beta){
    return(1/rgamma(n,alpha,beta))
  }

  rprior<- function(n, hyperparameters) {
    stopifnot(is.list(hyperparameters))
    stopifnot(length(hyperparameters) >= 4)
    stopifnot(all(c("mu","sigma","alpha","beta") %in% names(hyperparameters) ))
    stopifnot(length(hyperparameters$mu) == ncol(as.matrix(hyperparameters$sigma)))
    return(list(theta = rprior_coef(n, hyperparameters$mu, hyperparameters$sigma),
                sigma2 = rprior_sigma(n, hyperparameters$alpha, hyperparameters$beta)))
  }

  #### X Data ####
  rX <- gen_x()$rX

  #### Y Data ####
  slim_model_mat <- function(x, n.theta) {

    stopifnot(is.matrix(x))
    add.int <- FALSE
    if(all(x[,1] == 1)) {
      x <- x[,-1, drop = FALSE]
      n.theta <- n.theta - 1
      add.int <- TRUE
    }
    p <- n.theta
    p.x <- ncol(x)

    if(p > p.x) {
      n.sq <- min(4,p - p.x)
    } else {
      warning("No square terms added. Length of theta less than ncol x")
      if(add.int) x <- cbind(1,x)
      return(x)
    }
    sq.vars <- c(1,3,7)
    sq.mat <- matrix(NA_real_, nrow=nrow(x), ncol = n.sq)
    for( i in 1:n.sq){
      if(i < 4) {
        sq.mat[,i] <- x[,sq.vars[i]]^2
      } else if(i == 4 & p.x >=15) {
        sq.mat[,i] <- x[,13] * x[,15]
      }

    }
    if(add.int) x <- cbind(1,x)
    return(cbind(x, sq.mat))
  }

  slim_sq_corr <- function(corr, n.x, n.theta) {

    p <- n.theta
    p.x <- n.x

    if(p > p.x) {
      n.sq <- min(4,p - p.x)
    } else {
      return(NULL)
    }
    sq.vars <- c(1,3,7)
    cross.vars <- c(13,15)
    cor.mat <- matrix(0, nrow=n.sq, ncol = n.sq)
    grps <- lapply(seq(0,20,5), function(nn) sq.vars <= (nn + 5) &
                     sq.vars > nn)
    for( i in grps){
      if(all(isFALSE(i))) next
      idx <- which(i)
      n.row <- length(idx)
      cor.mat[idx,idx] <- matrix(2 * corr^2, n.row, n.row)
      diag(cor.mat)[idx] <- 2
    }
    if (n.sq == 4) {

      if(any(sq.vars %in% cross.vars)) {
        idx <- which(sq.vars %in% cross.vars)
        cor.mat[idx,4] <- cor.mat[4,idx] <- 2 * corr #all vars are 1
      }
      cor.mat[4,4] <- 1 + corr * corr
    }

    return(cor.mat)
  }
  rdata <- function(n, x, theta, sigma2, ...) {
    dots <- list(...)
    corr.x <- dots$corr
    scale <- dots$scale
    if(is.null(corr.x)) corr.x <- 0
    if(is.null(scale)) scale <- TRUE
    # if(is.null(dots$method)) dots$method <- "linear"
    #
    # method <- match.arg(dots$method, c("linear","nonlinear"))

    if(ncol(x) > length(theta)) {
      x <- x[, 1: length(theta), drop=FALSE]
      warning("Ncol X > length(theta). Only using first length(theta) columns of X.")
    }
    is.intercept <- all(x[,1] == 1)
    if(is.intercept) {
      intercept <- theta[1]
      theta <- theta[-1]
      x <- x[,-1, drop=FALSE]
    } else {
      intercept <- 0
    }

    if(ncol(x) + 4 > length(theta)) {
      x <- x[,1:(length(theta) - 4), drop = FALSE]
    }
    p <- ncol(x)
    corr.mat <- corr_mat_construct(corr.x, p)
    diag(corr.mat) <- 1
    corr.mat <- as.matrix(
      Matrix::bdiag(corr.mat, slim_sq_corr(corr.x,p, length(theta)) ))
    x <- slim_model_mat(x, length(theta))



    theta_norm <- c(t(theta) %*% corr.mat %*% theta)
    theta_scaled <- if(scale) {
                      theta/sqrt(theta_norm)
                    } else {
                      theta
                    }#* sqrt(sigma2)
    mu <- x %*% theta_scaled + intercept
    Y <- rnorm(n, mu, sqrt(sigma2))
    theta <- c(intercept, theta_scaled)

    return(list(Y = Y,
                mu = mu, eta = mu,
                link = gaussian()$linkfun,
                invlink = gaussian()$linkinv,
                theta = theta,
                data_gen_function = slim_model_mat,
                model_matrix = slim_model_mat))
  }

  #### Parameters ####
  rparam <- get_param()$gaussian

  #### Posterior on Coefficients ####
  rpost_coef <- function(x){}
  rpost_sigma <- function(x){}

  rpost <- function(n.samp, x, y, hyperparameters = NULL ,...) {
    # a <- hyperparameters$alpha
    # b <- hyperparameters$beta
    # m <- hyperparameters$mu
    # s <- hyperparameters$sigma
    # n <- length(y)
    dots <- list(...)
    model <- NULL
    test <- list(eta = NULL, mu = NULL)


    method <- dots$method

    # post_prec_beta <- crossprod(x) + Lambda
    # post_cross_beta <- crossprod(x, y) + Lambda %*% m
    #
    # post_mu <- solve(post_prec_beta, post_cross_beta)
    #
    # post_a <- a + n * 0.5
    # post_b <- b + 0.5 * (sum(y^2) -  t(post_mu) %*% post_prec_beta %*% post_mu )
    #
    # post_sigma <- 1/rgamma(n.samp, post_a, post_b)
    #
    # post_theta <- sapply(1:n.samp, function(i) rmvnorm(1, post_mu, solve(post_prec_beta)*post_sigma[i]))
    if (method == "conjugate") {
      X.test <- dots$X.test
      if ( !("Lambda" %in% names(hyperparameters) ) ) {
        hyperparameters$Lambda <- solve(hyperparameters$sigma)
      }
      if(dim(hyperparameters$Lambda)[2] != ncol(x)) stop("dimensions of priors must equal ncol of x")
      conjFit <- bayesConjRegNormal(n.samp, hyperparameters,
                              y, x)
      theta <- conjFit$theta
      sigma <- conjFit$sigma
      model <- "conjugate"
      eta <- mu <- x %*% theta
      if(!is.null(X.test)){
        testEta <- testMu <- X.test %*% theta
        test <- list(eta = testEta, mu = testMu)
      }

    } else if (method == "stan") {
      require(rstan)

      if (all(x[,1]==1)) x <- x[,-1]


      p <- ncol(x)
      n <- nrow(x)

      stan_dir <- dots$stan_dir
      m0 <- dots$m0
      scale_intercept <- dots$scale_intercept
      chains <- dots$chains
      X.test <- dots$X.test

      if(is.null(m0)) m0 <- round(0.1 * p)
      if(is.null(scale_intercept)) scale_intercept <- 2.5
      if(is.null(chains)) chains <- 4

      x_c <- scale(x, scale=FALSE)
      qrdecomp <- qr(x_c)
      Q_x <- qr.Q(qrdecomp) * sqrt(n-1)
      R_inv <- solve(qr.R(qrdecomp)/sqrt(n-1))

      stan_dat <- list(N = n,
                       P = p,
                       Y = y,
                       X = x,
                       mean_x = attr(x_c, "scaled:center"),
                       Q = Q_x,
                       R_inv_x = R_inv,
                       m0 = m0,
                       scale_intercept = scale_intercept)

      warmup <- max(n.samp*(2-1/chains), 1000)
      iter <- warmup + ceiling(n.samp/chains)

      rstan::rstan_options(auto_write = TRUE)
      stanModel <- stan_model(stan_dir)
      stanFit <- sampling(stanModel, data=stan_dat, iter=iter,
                          warmup = warmup, chains=chains, pars = c("y_hat","theta","sigma"))
      samples <- extract(stanFit, pars= c("y_hat","theta","sigma"))
      eta <- mu <- t(samples$y_hat)
      theta <- t(samples$theta)
      sigma <- (samples$sigma^2)
      model <- stanFit
      if(!is.null(X.test)){
        testEta <- testMu <- X.test %*% theta
        test <- list(eta = testEta, mu = testMu)
      }
    } else {
      stopifnot(match.arg(method, c("conjugate","stan")))
    }
    return(list(theta=theta,
                sigma=sigma, eta = eta,
                mu = mu, model=model, test = test))

  }

  mf.linpred <- function(x, theta) {
    return(x %*% theta)
  }
  sel.pred.fun <- function(method="linpred") {
    return(mf.linpred)
  }

  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = slim_model_mat,
              rparam = rparam,
              link = gaussian()$linkfun,
              sel.pred.fun = sel.pred.fun,
              invlink = gaussian()$linkinv))
}
