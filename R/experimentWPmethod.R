experimentWPMethod <- function(target, hyperparameters, conditions) {
  n <- conditions$n
  p <- conditions$p

  n.samps <- conditions$n.samps
  penalty <- conditions$penalty
  lambda.min.ratio <- conditions$lambda.min.ratio
  n.lambda <- conditions$n.lambda
  penalty_method <- conditions$penalty.factor
  family <- conditions$family
  epsilon <- 0.05
  otmaxit <- 100
  # pseudo.obs <- conditions$pseudo.obs
  stan_dir <- conditions$stan_dir
  calc_w2_post <- conditions$calc_w2_post
  if(is.null(calc_w2_post)) calc_w2_post <- FALSE
  # if(is.null(pseudo.obs)) pseudo.obs <- 0
  pseudo.obs <- 0
  posterior.method <- conditions$posterior.method
  transport.method <- conditions$transport.method
  if(is.null(transport.method)) transport.method <- "hilbert" #"univariate.approximation.pwr"
  not.only.timing <- conditions$not.only.timing
  if(is.null(not.only.timing)) not.only.timing <- FALSE
  L0 <- conditions$L0
  if(is.null(L0)) L0 <- FALSE

  sa_seq <- sort(unique(c(2,5,floor(seq(ceiling(p/5),p,floor(p/5))))))
  sa_max_time <- 64800
  sa_prop <- "random"
  # FSAiter <- 10*1:ceiling(p/2)
  # RSAiter <- if( p %% 2) { rev(FSAiter)[-1] } else { rev(FSAiter) }
  # SAiter <- c(FSAiter, RSAiter)
  SAiter <- 10 * ceiling(p/2)
  SAtemps <- 50
  # SAiter <- 1
  # SAtemps <- 1

  #IP sequence
  if( p < 200 ) {
    ip_seq <- 1:p
  } else {
    ip_seq <- sa_seq
  }

  #w2 dist param
  wp_alg <- conditions$wp_alg
  if(is.null(wp_alg)) wp_alg <- "greenkhorn"

  # MEAN FUNCTION FOR IMPORTANCE
  pred.fun <- switch(family,
                     "gaussian" = target$sel.pred.fun("linpred"),
                     "binomial" = target$sel.pred.fun("linpred"),
                     "exponential" = function(x,theta) {
                       return(x[,-1,drop=FALSE] %*% theta[-1,,drop=FALSE])
                     })

  # SETUP PARAMETERS
  param <- target$rparam()
  p_star <- min(length(param$theta),p)
  param$theta <- param$theta[1:p_star]
  full_param <- c(param$theta, rep(0, p-p_star))

  #set up data
  X <- target$X$rX(n, target$X$corr, p)
  X_sing <- matrix(c(target$X$rX(1, target$X$corr, p)), nrow=1, ncol=p)
  X_new <- target$X$rX(n, target$X$corr, p)
  X_neighborhood <- cbind(1, CoarsePosteriorSummary::rmvnorm(nsamples = p*3,
                                                    mean = X_sing[,-1,drop=FALSE],
                                                    covariance = cov(X[,-1,drop=FALSE])/n))


  data <- target$rdata(n, X[,1:p_star,drop=FALSE], c(param$theta),
                       param$sigma2, method = "modified.friedman", corr = target$X$corr)
  single_data <- target$rdata(1, X_sing[,1:p_star,drop=FALSE], c(param$theta),
                              param$sigma2, method = "modified.friedman", corr = target$X$corr)
  new_data <- target$rdata(n, X_new[,1:p_star,drop=FALSE], c(param$theta), param$sigma2, method = "modified.friedman", corr = target$X$corr)

  Y <- data$Y
  single_Y <- single_data$Y
  new_Y_new_X <- new_data$Y

  #True means
  true_mu <- data$mu
  new_mu <- new_data$mu
  new_mu_sing <- single_data$mu



  # add non-linearities if binomial
  if(family == "binomial") {
    cov_X <- cov(X[,-1,drop = FALSE])
    # form <- formula("~.**2")
    form <- formula(paste0("~",paste0("I(X",1:(p-1),"^2)", collapse=" + ")))

    X <- model.matrix(object=form, data = data.frame(X[,-1,drop = FALSE]))
    p <- ncol(X)
    X_new <- model.matrix(form, data = data.frame(X_new[,-1,drop = FALSE]))
    #resamp x_neighb for larger space
    X_neighborhood <- CoarsePosteriorSummary::rmvnorm(nsamples = p*3,
                                                      mean = X_sing[,-1,drop = FALSE],
                                                      covariance = cov_X/n)
    X_sing <- model.matrix(form, data = data.frame(X_sing[,-1,drop = FALSE]))
    X_neighborhood <- model.matrix(form, data = data.frame(X_neighborhood))

  }

  # indices for different data
  single_idx <- 1
  new_idx <- 2:(n+1)
  neighb_idx <- (n+2):(n+1 + p*3)

  #sample theta
  if(family != "binomial") {
    post_sample <- target$rpost(n.samps, X, Y, hyperparameters,
                              method = posterior.method, stan_dir = stan_dir,
                              X.test = rbind(X_sing, X_new, X_neighborhood),
                              chains = 1, m0 = 20,
                              is.exponential = TRUE)

    post_interp <- post_sample
    theta <- theta_new <- theta_sing <- post_interp$theta #regression coef
    if(nrow(theta) != ncol(X)) theta <- t(theta)
    t_theta <- t(theta)
    sigma <- post_interp$sigma
  } else {
    post_sample <- target$rpost(n.samps, cbind(1, data$data_gen_functions( X ) ),
                                Y,
                                hyperparameters,
                                method = posterior.method,
                                stan_dir = stan_dir,
                                X.test = cbind(1, data$data_gen_functions( rbind(X_sing, X_new, X_neighborhood) ) ),
                                chains = 1, m0 = 20,
                                is.exponential = TRUE)
    # post_interp <- target$rpost(n.samps, X, Y, NULL, method = "logistic", stan_dir = "exec/Stan/logistic_horseshoe_noQR.stan",
    #                             X.test = rbind(X_sing, X_new, X_neighborhood),
    #                             chains = 1, m0 = 20)
    theta <- lm.fit(X, post_sample$eta)$coefficients
    theta_new <-  lm.fit(X_new, post_sample$test$eta[new_idx,,drop=FALSE])$coefficients
    theta_sing <-  lm.fit(X_neighborhood, post_sample$test$eta[neighb_idx,,drop=FALSE])$coefficients
    post_interp <- list(theta = theta,
                        eta = post_sample$eta, #X %*% theta,
                        mu = post_sample$mu, #data$invlink(X %*% theta),
                        test = list(eta = post_sample$test$eta,#rbind(X_sing %*% theta_sing, X_new %*% theta_new, X_neighborhood %*% theta_sing)
                                    mu = post_sample$test$mu)
                        )
    t_theta <- t(theta)
    calc_w2_post <- FALSE
  }


  # if(family == "exponential" & posterior.method == "stan-cox") {
  #   log_dL0 <- rstan::extract(post_sample$model, pars= c("log_dL0"))$log_dL0
  #   intercept <- t(rowMeans(log_dL0))
  #   print(dim(post_sample$theta))
  #   post_sample$theta <-  rbind(intercept, post_sample$theta)
  #   print(dim(post_sample$theta))
  #   post_sample$eta <-  post_sample$eta + matrix(intercept, nrow = n, ncol=n.samps, byrow = TRUE)
  #   post_sample$test$eta <- post_sample$test$eta +
  #     matrix(intercept, nrow = n+1, ncol=n.samps, byrow = TRUE)
  # }
   #variance (if it exists for model)

  #functions of theta
  E_theta <- colMeans(theta)
  # theta_norm <- rowSums(theta^2)

  #conditional and marginal natural parameter means
  cond_eta <- post_sample$eta
  marg_eta <- rowMeans(cond_eta)

  cond_eta_new <- post_sample$test$eta[new_idx,, drop=FALSE]
  marg_eta_new <- rowMeans(cond_eta_new)

  cond_eta_sing <- post_sample$test$eta[single_idx,,drop=FALSE]
  marg_eta_sing <- rowMeans(cond_eta_sing)

  cond_eta_neighb <- post_sample$test$eta[neighb_idx,,drop=FALSE]
  marg_eta_neighb <- rowMeans(cond_eta_neighb)

  #conditional and marginal means
  # family != "exponential" &
  if(!grepl("cox", posterior.method)) {
    cond_mu <- post_sample$mu
    marg_mu <- rowMeans(cond_mu)

    cond_mu_new <- post_sample$test$mu[new_idx,, drop=FALSE]
    marg_mu_new <- rowMeans(cond_mu_new)

    cond_mu_sing <- post_sample$test$mu[single_idx,,drop=FALSE]
    marg_mu_sing <- rowMeans(cond_mu_sing)

    cond_mu_neighb <- post_sample$test$mu[neighb_idx,,drop=FALSE]
    marg_mu_neighb <- rowMeans(cond_mu_neighb)
  } else {
    # cond_mu <- post_sample$mu$S$surv
    # marg_mu <- apply(cond_mu, 2:3, mean)
    #
    # cond_mu_new <- post_sample$test$mu$S[,, -1,drop=FALSE]
    # marg_mu_new <- apply(cond_mu_new, 2:3, mean)
    #
    # cond_mu_sing <- post_sample$test$mu$S[,,1,drop=FALSE]
    # marg_mu_sing <- apply(cond_mu_sing, 2:3, mean)
    # X <- X[,-1]
    # X_new <- X_new[,-1]
    # X_sing <- X_sing[,-1]
    cond_mu <- data$invlink(post_sample$eta)
    marg_mu <- rowMeans(cond_mu)

    cond_mu_new <- data$invlink(post_sample$test$eta[new_idx,, drop=FALSE])
    marg_mu_new <- rowMeans(cond_mu_new)

    cond_mu_sing <- data$invlink(post_sample$test$eta[single_idx,,drop=FALSE])
    marg_mu_sing <- rowMeans(cond_mu_sing)

    cond_mu_neighb <- data$invlink(post_sample$test$eta[neighb_idx,,drop=FALSE])
    marg_mu_neighb <- rowMeans(cond_mu_neighb)
  }


  #penalty terms
  penalty_fact <- set_penalty_factor(theta = theta, method = penalty_method, intercept = TRUE,
                                     x = X, y = cond_eta, transport.method = transport.method)
  penalty_factN <- set_penalty_factor(theta = theta, method = penalty_method, intercept = TRUE,
                                      x = X_new,  y = cond_eta_new, transport.method = transport.method)
  penalty_factO <- set_penalty_factor(theta = theta, method = penalty_method, intercept = TRUE,
                                      x = X_sing, y = cond_eta_sing, transport.method = transport.method)
  proj_penalty_fact <- set_penalty_factor(theta = theta, "distance", intercept = TRUE)
  HC_penalty_fact <- set_penalty_factor(theta = theta, "expectation", intercept = TRUE)

  #optional L0
  if(L0){
    cat(paste0("Running L0 methods only, same data: ", date(), "\n"))
    L0list <- list(Selection = NULL,
                   Projection = NULL)
    L0list$Selection <- WPL0(X = X, Y = cond_eta, theta = theta,
                             p = 2, ground_p = 2, method = "selection.variable",
                             transport.method = transport.method, epsilon = epsilon,
                             maxit = otmaxit)
    L0list$Projection <- WPL0(X = X, Y = cond_eta, theta = theta,
                              p = 2, ground_p = 2, method = "projection",
                              transport.method = transport.method, epsilon = epsilon,
                              maxit = otmaxit)

    cat("L0 Distance Calculations\n")
    W2L0 <- distCompare(L0list,
                        target = list(posterior = NULL,
                                      mean = cond_mu),
                        method = wp_alg,
                        quantity=c("mean"),
                        parallel=NULL,
                        transform = data$invlink,
                        epsilon = epsilon,
                        niter = otmaxit)
    mseL0 <- distCompare(L0list,
                         target = list(posterior = NULL,
                                       mean = true_mu),
                         method = "mse",
                         quantity="mean",
                         parallel=NULL,
                         transform = data$invlink,
                         epsilon = epsilon,
                         niter = otmaxit)

    outList <- list (
      W2_dist = W2L0,
      mse = mseL0,
      time = NULL,
      order = list(selection = L0list$Selection$minCombPerActive,
                   projection = L0list$Projection$minCombPerActive)
    )
    return(outList)
  }

  # augDat <- augPseudo(X, cond_eta, theta, theta_norm, pseudo.obs, n, same=TRUE)
  # lambdas <- calc.lambdas(augDat, lambda.min.ratio, penalty_fact, n.lambda)
  cat(paste0("Running methods, same data: ", date(),"\n"))

  #### In sample ####
  #IP
  cat("   Selection: IP,\n")
  time <- proc.time()
  ip <- W2IP(X = X, Y = cond_eta, theta = theta,
             display.progress=FALSE,
             transport.method = transport.method,
             model.size = ip_seq,
             infimum.maxit = 100, solution.method = "cone",
             parallel = NULL)
  ipTime <- proc.time() - time
  # trajSel <- selDist$theta

  cat("   Selection: IP, Lasso, ")
  #selection variable
  time <- proc.time()
  lassoSel <- W2L1(X, cond_eta, theta, family="gaussian", penalty=penalty,
                   penalty.factor = penalty_fact, nlambda = n.lambda, alpha = 0.99,
                   gamma = 1.1,
                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=1e2,
                   maxit = 1e6,
                   display.progress=FALSE,
                   transport.method = transport.method,
                   method = "selection.variable")
  selTime <- proc.time() - time
  # trajSel <- selDist$theta

  cat(" HC, ")
  #carvalho method
  time <- proc.time()
  lassoHC <- HC(X, cond_eta, theta = theta,
                family="gaussian", penalty=penalty, method = "selection.variable",
                penalty.factor=HC_penalty_fact, nlambda = n.lambda, alpha = 0.99, gamma = 1.1,
                lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  hcTime <- proc.time() - time

  cat(" SW, ")
  #stepwise
  time <- proc.time()
  step <- WPSW(X, Y = cond_eta, theta, force=1, p=2,
               direction = "backward", method = "selection.variable",
               transport.method = transport.method,
               display.progress = FALSE)
  stepTime <- proc.time() - time
  # trajStep <- step$theta

  cat(" SA\n")
  #simulated annealing
  annealTime <- NULL
  # if(n > 512 | p > 11){
  #   anneal <- WPSA(X=X, Y=cond_eta, theta=theta,
  #                  force = 1, p=2, model.size = 5, iter = SAiter, temps = SAtemps,
  #                  options = list(method = "selection.variable",
  #                                 energy.distribution = "boltzman",
  #                                 transport.method = transport.method,
  #                                 cooling.schedule="exponential"),
  #                  display.progress=TRUE)
  #   annealTime <- NULL
  # }
  # else {
  time <- proc.time()
  anneal <- WPSA(X=X, Y=cond_eta, theta=theta,
                 force = 1, p=2, model.size = sa_seq, iter = SAiter, temps = SAtemps,
                 options = list(method = "selection.variable",
                                energy.distribution = "boltzman",
                                transport.method = transport.method,
                                cooling.schedule="exponential",
                                proposal.method = sa_prop),
                 display.progress=FALSE, max.time = sa_max_time)
  annealTime <- proc.time() - time
  cat(paste0(anneal$message,"\n"))
  if(anneal$message != "completed") {
    annealTime <- paste0(">", annealTime)
  }

  cat("   Projection: Lasso, ")
  #projection
  time <- proc.time()
  lassoProj <- W2L1(X, cond_eta, theta, family="gaussian", penalty=penalty,
                    penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                    maxit = 1e6, alpha = 0.99, gamma = 1.1,
                    display.progress=FALSE,
                    transport.method = transport.method,
                    method = "projection")
  projTime <- proc.time() - time
  # trajProj <- projDist$theta

  cat(" HC, ")
  time <- proc.time()
  PlassoHC <- HC(X, cond_eta, theta = theta,
                 family="gaussian", penalty=penalty, method = "projection",
                 alpha = 0.99, gamma = 1.1,
                 penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                 lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  PhcTime <- proc.time() - time

  cat(" SW, ")
  #stepwise
  time <- proc.time()
  Pstep <- WPSW(X, Y = cond_eta, theta, force=1, p=2,
                direction = "backward", method = "projection",
                transport.method = transport.method,
                display.progress = FALSE)
  PstepTime <- proc.time() - time
  # trajStep <- step$theta

  cat(" SA\n")
  #simulated annealing
  annealTime <- NULL

  time <- proc.time()
  Panneal <- WPSA(X=X, Y=cond_eta, theta=theta,
                  force = 1, p=2, model.size = sa_seq, iter = SAiter, temps = SAtemps,
                  options = list(method = "projection",
                                 energy.distribution = "boltzman",
                                 transport.method = transport.method,
                                 cooling.schedule="exponential",
                                 proposal.method = sa_prop),
                  display.progress=FALSE, max.time = sa_max_time)
  PannealTime <- proc.time() - time

  # }
  # trajAnneal <- anneal$theta
  cat(paste0(Panneal$message,"\n"))
  if(Panneal$message != "completed") {
    PannealTime <- paste0(">", PannealTime)
  }

  if (not.only.timing) {
    inSampModels <- list("Binary Programming" = ip,
                         "Lasso" = lassoSel,
                         "Simulated Annealing" = anneal,
                         "Stepwise" = step,
                         "Hahn-Carvalho" = lassoHC)
    inSampModelsProj <- list("Lasso" = lassoProj,
                             "Simulated Annealing" = Panneal,
                             "Stepwise" = Pstep,
                             "Hahn-Carvalho" = PlassoHC)
    rm("ip","lassoSel", "anneal","step","lassoHC")
    rm("lassoProj", "Panneal","Pstep","PlassoHC")

    cat("Calculating distances\n")
    if( calc_w2_post){
      W2_insamp <- distCompare(inSampModels, target = list(posterior = theta,
                                                           mean = cond_mu),
                               method = wp_alg,
                               quantity=c("posterior","mean"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit)
      mse_insamp <- distCompare(inSampModels, target = list(posterior = full_param,
                                                            mean = true_mu),
                                method = "mse",
                                quantity=c("posterior","mean"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      PW2_insamp <- distCompare(inSampModelsProj, target = list(posterior = theta,
                                                                mean = cond_mu),
                                method = wp_alg,
                                quantity=c("posterior","mean"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      Pmse_insamp <- distCompare(inSampModelsProj, target = list(posterior = full_param,
                                                                 mean = true_mu),
                                 method = "mse",
                                 quantity=c("posterior","mean"),
                                 parallel=NULL,
                                 transform = data$invlink,
                                 epsilon = epsilon,
                                 niter = otmaxit)
    }
    else {
      W2_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
                                                           mean = cond_mu),
                               method = wp_alg,
                               quantity=c("mean"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit)
      mse_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
                                                            mean = true_mu),
                                method = "mse",
                                quantity="mean",
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      PW2_insamp <- distCompare(inSampModelsProj, target = list(posterior = NULL,
                                                                mean = cond_mu),
                                method = wp_alg,
                                quantity=c("mean"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      Pmse_insamp <- distCompare(inSampModelsProj, target = list(posterior = NULL,
                                                                 mean = true_mu),
                                 method = "mse",
                                 quantity="mean",
                                 parallel=NULL,
                                 transform = data$invlink,
                                 epsilon = epsilon,
                                 niter = otmaxit)
    }

    rm(inSampModels)
    rm(inSampModelsProj)

    #### new X variable ####
    #mse on new outcome data from same paramters and different X
    #new method
    # augDatN <- augPseudo(X_new, cond_mu_new, theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatN, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("Running methods, new X variable: ", date(),"\n"))
    cat("  Selection: IP")
    ipN <- W2IP(X = X_new, Y = cond_eta_new, theta = theta_new,
                display.progress=TRUE,
                transport.method = transport.method,
                model.size = ip_seq,
                infimum.maxit = 100, solution.method = "cone",
                parallel = NULL)

    cat(" Lasso")
    lassoSelN <- W2L1(X_new, cond_eta_new, theta_new, family="gaussian",
                      penalty=penalty, alpha = 0.99, gamma = 1.1,
                      penalty.factor=penalty_factN, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=100,
                      maxit = 1e5,
                      transport.method = transport.method,
                      display.progress=TRUE, method = "selection.variable")
    # trajSelN <- extractCoef(lassoSelN)

    cat(" HC, ")
    #carvalho method
    lassoHCN <- HC(X_new, cond_eta_new, theta = theta_new,
                   family="gaussian", penalty=penalty, alpha = 0.99, gamma = 1.1,
                   penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, maxit = 1e5)

    cat(" SW")
    #stepwise
    stepN <- WPSW(X_new, cond_eta_new, theta_new, force=1, p=2,
                  direction = "backward", method = "selection.variable",
                  transport.method = transport.method,
                  display.progress = TRUE)

    cat(" SA")
    annealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=theta_new,
                     force = 1, p=2, model.size = sa_seq, iter = SAiter,
                     temps = SAtemps,
                     options = list(method = "selection.variable",
                                    energy.distribution = "boltzman",
                                    transport.method = transport.method,
                                    cooling.schedule="exponential",
                                    proposal.method = sa_prop),
                     display.progress = TRUE, max.time = sa_max_time)
    cat(paste0(annealN$message,"\n"))
    cat("\n")

    cat("  Projection: Lasso")
    #projection
    lassoProjN <- W2L1(X_new, cond_eta_new, theta_new, family="gaussian", penalty=penalty,
                       penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                       lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                       maxit=1e5, alpha = 0.99, gamma = 1.1,
                       transport.method = transport.method,
                       display.progress=TRUE, method = "projection")

    cat(" SW")
    #stepwise
    PstepN <- WPSW(X_new, cond_eta_new, theta_new, force=1, p=2,
                   direction = "backward", method = "projection",
                   transport.method = transport.method,
                   display.progress = TRUE)
    cat(" HC, ")
    #HC
    PlassoHCN <- HC(X_new, cond_eta_new, theta = theta_new, alpha = 0.99, gamma = 1.1,
                    family="gaussian", penalty=penalty, method = "projection",
                    penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
    cat(" SA")
    # anneal
    PannealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=theta_new,
                      force = 1, p=2, model.size = sa_seq, iter = SAiter,
                      temps = SAtemps,
                      options = list(method = "projection",
                                     energy.distribution = "boltzman",
                                     transport.method = transport.method,
                                     cooling.schedule="exponential",
                                     proposal.method = sa_prop),
                      display.progress = TRUE, max.time = sa_max_time)
    cat(paste0(PannealN$message,"\n"))
    cat("\n")
    # }
    # trajAnnealN <- annealCoef(annealN, theta)
    newXModels <- list("Binary Programming" = ipN,
                       "Lasso" = lassoSelN,
                       "Simulated Annealing" = annealN,
                       "Stepwise" = stepN,
                       "Hahn-Carvalho" = lassoHCN)
    rm("ipN","lassoSelN", "annealN","stepN","lassoHCN")
    newXModelsP <- list("Lasso" = lassoProjN,
                        "Simulated Annealing" = PannealN,
                        "Stepwise" = PstepN,
                        "Hahn-Carvalho" = PlassoHCN)
    rm("lassoProjN",
       "PannealN", "PstepN","PlassoHCN")

    cat("Calculating distances\n")
    if( calc_w2_post){
      W2_newX <- distCompare(newXModels, target = list(posterior = theta_new,
                                                       mean = cond_mu_new),
                             method = wp_alg,
                             quantity=c("posterior","mean"),
                             parallel=NULL,
                             transform = data$invlink,
                             epsilon = epsilon,
                             niter = otmaxit)
      mse_newX <- distCompare(newXModels, target = list(posterior = full_param,
                                                        mean = new_mu),
                              method = "mse",
                              quantity=c("posterior","mean"),
                              parallel=NULL,
                              transform = data$invlink,
                              epsilon = epsilon,
                              niter = otmaxit)
      PW2_newX <- distCompare(newXModelsP, target = list(posterior = theta_new,
                                                         mean = cond_mu_new),
                              method = wp_alg,
                              quantity=c("posterior","mean"),
                              parallel=NULL,
                              transform = data$invlink,
                              epsilon = epsilon,
                              niter = otmaxit)
      Pmse_newX <- distCompare(newXModelsP, target = list(posterior = full_param,
                                                          mean = new_mu),
                               method = "mse",
                               quantity=c("posterior","mean"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit)
    }
    else {
      W2_newX <- distCompare(newXModels, target = list(posterior = NULL,
                                                       mean = cond_mu_new),
                             method = wp_alg,
                             quantity=c("mean"),
                             parallel=NULL,
                             transform = data$invlink,
                             epsilon = epsilon,
                             niter = otmaxit)
      mse_newX <- distCompare(newXModels, target = list(posterior = NULL,
                                                        mean = new_mu),
                              method = "mse",
                              quantity="mean",
                              parallel=NULL,
                              transform = data$invlink,
                              epsilon = epsilon,
                              niter = otmaxit)
      PW2_newX <- distCompare(newXModelsP, target = list(posterior = NULL,
                                                         mean = cond_mu_new),
                              method = wp_alg,
                              quantity=c("mean"),
                              parallel=NULL,
                              transform = data$invlink,
                              epsilon = epsilon,
                              niter = otmaxit)
      Pmse_newX <- distCompare(newXModelsP, target = list(posterior = NULL,
                                                          mean = new_mu),
                               method = "mse",
                               quantity="mean",
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit)
    }


    rm(newXModels, newXModelsP)



    #### new method, single datapoint ####
    # augDatO <- augPseudo(X_sing, cond_mu_sing, theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatO, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("Running methods, single data point: ", date(),"\n"))
    cat("  Selection: IP")
    ipO <- W2IP(X = X_sing, Y = cond_eta_sing, theta = theta_sing,
                display.progress=TRUE,
                transport.method = transport.method,
                model.size = ip_seq,
                infimum.maxit = 100, solution.method = "cone",
                parallel = NULL)

    cat(" Lasso")
    lassoSelO <- W2L1(X_sing, cond_eta_sing, theta_sing, family="gaussian", penalty=penalty,
                      penalty.factor=penalty_factO, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=100,
                      maxit=1e5, alpha = 0.99, gamma = 1.1,
                      transport.method = transport.method,
                      display.progress=TRUE, method = "selection.variable")

    cat(" SW")
    #stepwise
    stepO <- WPSW(X_sing, cond_eta_sing, theta_sing, force=1, p=2,
                  direction = "backward",
                  method = "selection.variable",
                  transport.method = transport.method,
                  display.progress = TRUE)

    cat(" SA")
    #simulated annealing
    annealO <- WPSA( X = X_sing, Y = cond_eta_sing, theta = theta_sing,
                     force = 1, p=2, model.size = sa_seq, iter = SAiter, temps = SAtemps,
                     options = list(method = "selection.variable",
                                    energy.distribution = "boltzman",
                                    transport.method = transport.method,
                                    cooling.schedule="exponential",
                                    proposal.method = sa_prop),
                     display.progress = TRUE , max.time = sa_max_time)
    cat(paste0(annealO$message,"\n"))
    cat("\n")

    cat("  Projection: Lasso")
    #projection
    lassoProjO <- W2L1(X_neighborhood, cond_eta_neighb, theta_sing, family="gaussian", penalty=penalty,
                       penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                       lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                       maxit=1e5, alpha = 0.99, gamma = 1.1,
                       transport.method = transport.method,
                       display.progress=TRUE, method = "projection")

    cat(" SW")
    #stepwise
    PstepO <- WPSW(X_neighborhood, cond_eta_neighb, theta_sing, force=1, p=2,
                   direction = "backward", method = "projection",
                   transport.method = transport.method,
                   display.progress = TRUE)
    cat(" HC, ")
    #HC
    PlassoHCO <- HC(X_neighborhood, cond_eta_neighb, theta = theta_sing, alpha = 0.99, gamma = 1.1,
                    family="gaussian", penalty=penalty, method = "projection",
                    penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
    cat(" SA")
    # anneal
    PannealO <-  WPSA(X_neighborhood, cond_eta_neighb, theta=theta_sing,
                      force = 1, p=2, model.size = sa_seq, iter = SAiter,
                      temps = SAtemps,
                      options = list(method = "projection",
                                     energy.distribution = "boltzman",
                                     transport.method = transport.method,
                                     cooling.schedule="exponential",
                                     proposal.method = sa_prop),
                      display.progress = TRUE, max.time = sa_max_time)
    cat(paste0(PannealO$message,"\n"))
    cat("\n")
    # }
    # }
    # trajAnnealN <- annealCoef(annealN, theta)
    singleModels <- list("Binary Programming" = ipO,
                         "Lasso" = lassoSelO,
                         "Simulated Annealing" = annealO,
                         "Stepwise" = stepO#,
    )
    rm("ipO", "lassoSelO", "annealO","stepO")
    singleModelsP <- list("Lasso" = lassoProjO,
                          "Simulated Annealing" = PannealO,
                          "Stepwise" = PstepO,
                          "Hahn-Carvalho" = PlassoHCO)
    rm("lassoProjO",
       "PannealO", "PstepO","PlassoHCO")
    #recalculate values for single obs
    singleModelsP <- lapply(singleModelsP, function(x) {
      x$eta <- lapply(x$theta, function(tt) X_sing %*% tt)
      return(x)
    })

    cat("Calculating distances\n")
    if( calc_w2_post){
      W2_single <- distCompare(singleModels, target = list(posterior = theta_sing,
                                                           mean = cond_mu_sing),
                               method = wp_alg,
                               quantity=c("posterior","mean"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit)
      mse_single <- distCompare(singleModels, target = list(posterior = full_param,
                                                            mean = new_mu_sing),
                                method = "mse",
                                quantity=c("posterior","mean"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      PW2_single <- distCompare(singleModelsP, target = list(posterior = theta_sing,
                                                             mean = cond_mu_sing),
                                method = wp_alg,
                                quantity=c("posterior","mean"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      Pmse_single <- distCompare(singleModelsP, target = list(posterior = full_param,
                                                              mean = new_mu_sing),
                                 method = "mse",
                                 quantity=c("posterior","mean"),
                                 parallel=NULL,
                                 transform = data$invlink,
                                 epsilon = epsilon,
                                 niter = otmaxit)
    }
    else {

      W2_single <- distCompare(singleModels, target = list(posterior = NULL,
                                                           mean = cond_mu_sing),
                               method = wp_alg,
                               quantity=c("mean"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit)
      mse_single <- distCompare(singleModels, target = list(posterior = NULL,
                                                            mean = new_mu_sing),
                                method = "mse",
                                quantity="mean",
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      PW2_single <- distCompare(singleModelsP, target = list(posterior = NULL,
                                                             mean = cond_mu_sing),
                                method = wp_alg,
                                quantity=c("mean"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      Pmse_single <- distCompare(singleModelsP, target = list(posterior = NULL,
                                                              mean = new_mu_sing),
                                 method = "mse",
                                 quantity="mean",
                                 parallel=NULL,
                                 transform = data$invlink,
                                 epsilon = epsilon,
                                 niter = otmaxit)
    }

    rm(singleModels)


    #### Variable importance ####
    cat(paste0("Running variable importance: ", date(),"\n"))
    # pred.fun.theta <- switch(family,
    #                          "gaussian" = theta,
    #                          "binomial" = theta,
    #                          "exponential" = list(baseline = post_sample$mu$S$base,
    #                                               param = theta))

    imp_insamp <- WPVI(X = X, Y = NULL, theta = theta, pred.fun = pred.fun,
                       p = 2, ground_p = 2, transport.method = "exact")
    imp_newX <- WPVI(X = X_new, Y = NULL, theta = theta_new, pred.fun = pred.fun,
                     p = 2, ground_p = 2, transport.method = "exact")
    imp_single <- WPVI(X = X_sing, Y = NULL, theta = theta_sing, pred.fun = pred.fun,
                     p = 2, ground_p = 2, transport.method = "exact")

  } else {
    W2_insamp <- W2_newX <- W2_single <- mse_insamp <-  mse_newX <-  mse_single <- NULL
    imp_insamp <- imp_newX <- imp_single <- NULL
  }




  #### generate output list ####

  outList <- list (
    W2_dist = list(inSamp = list(selection = W2_insamp,
                                 projection = PW2_insamp),
                   newX = list(selection = W2_newX,
                               projection = PW2_newX),
                   single = list(selection = W2_single,
                                 projection = PW2_single) ),
    mse = list(inSamp = list(selection = mse_insamp,
                             projection = Pmse_insamp),
               newX = list(selection = mse_newX,
                           projection = Pmse_newX),
               single = list(selection = mse_single,
                             projection = Pmse_single) ),
    importance = list(inSamp = imp_insamp,
                      newX = imp_newX,
                      single = imp_single),
    time = list(selection = list(ip = ipTime[3], lasso = selTime[3],
                                 HC = hcTime[3],
                                 step = stepTime[3],
                                 anneal = annealTime[3]),
                projection = list(lasso = projTime[3],
                                  HC = PhcTime[3], step = PstepTime[3],
                                  anneal = PannealTime[3]) )#,
  )

  return(outList)
}
