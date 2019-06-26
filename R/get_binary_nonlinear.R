
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
    x6sq <- x[,6]^2
    costerm <- cos(2 * x6sq + pi/2)
    x6term <- ifelse(2 * x6sq > 3/4 * pi, 1,  costerm)
    # x6term <- costerm
    return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) + 2 * x6term  + tanh(-pi/8 * x[,11] * x[,15]^2) + (2/5 * (x[,7])^3 - 1/8 * (x[,7])^2 - 2 * x[,7])*exp(-x[,7]^2/5) )
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) + 2 * x6sq  + sin(-pi/8 * x[,11] * x[,15]^2) + 1/20 * x[,7]^3 - 1/24 * x[,7]^2 - 1 * x[,7])
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) + 2 * x6sq  + tanh(-x[,11] * x[,15]^2) + 0.25 * x[,7])
    # return(1 * sin(pi/8 * x[,1] * x[,2] + pi/2) - 2 * x[,3]^2 + x[,4] + 0.5 * x[,5])
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
  rdata <- function(n, x, param, ...) {

    dots <- list(...)
    method <- dots$method
    if(is.null(method)) method <- "cart"
    method <- match.arg(method, c("brownian","modified.friedman","cart","gp","glm", "polynomial","random.interaction"))

    if(all(x[,1] == 1)) {
      p <- ncol(x)-1
    } else {
      p <- ncol(x)
    }
    logit.p <- NA
    extra <- NULL
    if (method == "brownian") {
      data_gen_functions <- brownian
      num.comb <- length(param) - ncol(x)
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
      logit.p <-  cbind(1, browned) %*% param
    }
    else if (method == "modified.friedman") {
      data_gen_functions <- modified_friedman
      stopifnot(ncol(x)>=6)
      logit.p <- modified_friedman(x[,-1])
      theta <- c(1,2,1,1)
    }
    else if (method == "cart") {
      stopifnot(ncol(x)>=6)
      output <- cart(x[,-1])
      data_gen_functions <- output$model
      logit.p <- output$eta
    }
    else if (method == "gp") {
      output <- apply(x[,-1],2,gp)
      theta <- param[1:(p+1)]
      theta <- abs(theta)/sqrt(sum(theta^2))
      logit.p <- cbind(1,output) %*% theta
      data_gen_functions <- gp
    }
    else if (method == "glm") {
      logit.p <- x %*% dots$param
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
                data_gen_functions=data_gen_functions, method=method, extra = extra))
  }

  #### Posterior on Coefficients ####
  rpost_coef <- function(x){NULL}
  rpost_sigma <- function(x){NULL}

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


      dnn <- stan_model(stan_dir)
      stan_dat <- list(N = n, P=p, L = L, X=x, Y=y, m0 = m0, nodes=nodes)
      vb.out <- vb(dnn,stan_dat, algorithm=algorithm, iter=iter,
                   grad_samples=grad_samples, eval_elbo=eval_elbo, eta=eta,
                   adapt_engaged=adapt_engaged, init=init, tol_rel_obj = tol,
                   output_samples = n.samp)
      pars <- extract(vb.out, pars=c("eta","prob"))
      eta <- pars$eta
      prob <- pars$prob

    }
    else if (dots$method == "bart") {
      require(dbarts)
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

      bartFit <- bart(x.train = x, y.train = y, ndpost = ndpost, nskip = nskip,
                      k = k, base=base, power=power,
                      numcut = numcut, ntree=ntree, keepevery = keepevery,
                      printevery = n.samp*keepevery/10, keeptrees=FALSE, x.test = X.test,
                      nchain = chains, nthread = min(chains,parallel::detectCores()-1))
      phi <- bartFit$yhat.train
      mu <- pnorm(phi)
      eta <- qlogis(prob)
      model <- bartFit
      if(!is.null(X.test)){
        test$mu  <- pnorm(bartFit$yhat.test)
        test$eta <- qlogis(test$mu )
      }


    } else if (dots$method == "logistic"){
      require(rstan)

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

      warmup <- nsamp
      iter <- ceiling(n.samp/chains) + warmup

      stanModel <- stan_model(stan_dir)
      stanFit <- sampling(stanModel, data = stan_dat, iter = iter,
                          warmup = warmup, chains = chains, pars = c("eta","theta","prob", "beta", "intercept"))
      samples <- extract(stanFit, pars= c("eta","theta","prob"))
      eta <- samples$eta
      mu <- samples$prob
      theta <- samples$theta
      model <- stanFit
      if(!is.null(X.test)){
        test$eta <- tcrossprod(X.test, samples$theta)
        test$mu <- plogis(test$eta )
      }


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
      X.test <- dots[["X.test"]]

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
      mu <- samples$prob
      theta <- samples$theta
      model <- stanFit
      # test$eta <- tcrossprod(X.test, samples$theta)
      # test$mu <- plogis(test$eta)

    } else if (dots$method == "gamm") {
      require(rstanarm)
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
      warmup <- n.samp
      iter <- ceiling(n.samp/chains) + warmup

      stan_dat <- data.frame(Y = Y)
      namesX <- paste0("X",1:p)
      for(i in 1:p) stan_dat[[namesX[i]]] <- x[,i]
      form <- as.formula(paste("Y ~ ", paste0("s(",namesX,")" , collapse= " + ")))

      stanFit <- stan_gamm4(form, data = stan_dat, family = binomial(),
                 chains = chains, iter = iter, algorithm="sampling",
                 warmup = warmup, prior = prior )
      eta <- posterior_linpred(stanFit)
      if(nrow(eta) > n.samp) eta <- eta[1:n.samp,]
      mu <- plogis(eta)
      params <- as.matrix(stanFit$stanfit)
      theta <- params[,seq_len(ncol(stanFit$x))]
      model <- stanFit
      if(!is.null(X.test)){
        if(all(X.test[,1]==1)) X.test <- X.test[,-1]
        X.test.lp <- mgcv::predict.gam(stanFit$jam, newdata=data.frame(X.test), type="lpmatrix")
        test$eta <- tcrossprod(X.test.lp, theta)
        test$mu <- plogis(test$eta)
      }
    }

    return(list(theta=theta, mu=mu, eta=eta, model=model, test = test))

  }

  #### Output ####
  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = data_gen_functions,
              rparam = rparam,
              link = binomial()$linkfun,
              invlink = binomial()$linkinv))
}
