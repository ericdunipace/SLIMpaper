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
  if(is.null(transport.method)) transport.method <- "greenkhorn" #"univariate.approximation.pwr"
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
  # if( p < 100 | (p < 200 & method == "exact")) {
  #   ip_seq <- 1:p
  # } else {
    ip_seq <- sa_seq
  # }

  #w2 dist param
  wp_alg <- conditions$wp_alg
  if(is.null(wp_alg)) wp_alg <- "gandkhorn"

  # SETUP PARAMETERS
  param <- target$rparam()
  p_star <- min(length(param$theta),p)
  param$theta <- param$theta[1:p_star]
  full_param <- c(param$theta, rep(0, p-p_star))

  #set up data
  X <- target$X$rX(n, target$X$corr, p)
  X_sing <- matrix(c(target$X$rX(1, target$X$corr, p)), nrow=1, ncol=p)
  X_new <- target$X$rX(n, target$X$corr, p)

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

  #sample theta
  post_sample <- target$rpost(n.samps, X, Y, hyperparameters, method = posterior.method, stan_dir = stan_dir, X.test = rbind(X_sing, X_new))
  if(family != "binomial") {
    post_interp <- post_sample
  } else {
    post_interp <- target$rpost(n.samps, X, Y, NULL, method = "logistic", stan_dir = stan_dir,
                                X.test = rbind(X_sing, X_new))
    calc_w2_post <- FALSE
  }

  theta <- post_interp$theta #regression coef
  if(nrow(theta) != ncol(X)) theta <- t(theta)
  t_theta <- t(theta)
  sigma <- post_interp$sigma #variance (if it exists for model)

  #functions of theta
  E_theta <- colMeans(theta)
  # theta_norm <- rowSums(theta^2)

  #conditional and marginal natural parameter means
  cond_eta <- post_sample$eta
  marg_eta <- rowMeans(cond_eta)

  cond_eta_new <- post_sample$test$eta[-1,, drop=FALSE]
  marg_eta_new <- rowMeans(cond_eta_new)

  cond_eta_sing <- post_sample$test$eta[1,,drop=FALSE]
  marg_eta_sing <- rowMeans(cond_eta_sing)

  #conditional and marginal means
  cond_mu <- post_sample$mu
  marg_mu <- rowMeans(cond_mu)

  cond_mu_new <- post_sample$test$mu[-1,, drop=FALSE]
  marg_mu_new <- rowMeans(cond_mu_new)

  cond_mu_sing <- post_sample$test$mu[1,,drop=FALSE]
  marg_mu_sing <- rowMeans(cond_mu_sing)

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
  time <- proc.time()
  ip <- W2IP(X = X, Y = cond_eta, theta = theta,
                   display.progress=FALSE,
                   transport.method = transport.method,
                   model.size = ip_seq,
                   infimum.maxit = 10, solution.method = "lp",
                   parallel = NULL)
  ipTime <- proc.time() - time
  # trajSel <- selDist$theta

  #selection variable
  time <- proc.time()
  lassoSel <- W2L1(X, cond_eta, theta, family="gaussian", penalty="selection.lasso",
                   penalty.factor = penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=1e2,
                   maxit = 1e6,
                   display.progress=FALSE,
                   transport.method = transport.method,
                   gamma = 1, method = "selection.variable")
  selTime <- proc.time() - time
  # trajSel <- selDist$theta

  #carvalho method
  time <- proc.time()
  lassoHC <- HC(X, cond_eta, theta = theta,
                family=family, penalty=penalty, method = "selection.variable",
                penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  hcTime <- proc.time() - time

  #stepwise
  time <- proc.time()
  step <- WPSW(X, Y = cond_eta, theta, force=1, p=2,
               direction = "backward", method = "selection.variable",
               transport.method = transport.method,
               display.progress = FALSE)
  stepTime <- proc.time() - time
  # trajStep <- step$theta

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


  #projection
  time <- proc.time()
  lassoProj <- W2L1(X, cond_eta, theta, family="gaussian", penalty=penalty,
                    penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                    maxit = 1e6,
                    display.progress=FALSE,
                    transport.method = transport.method,
                    gamma = 1, method = "projection")
  projTime <- proc.time() - time
  # trajProj <- projDist$theta

  time <- proc.time()
  PlassoHC <- HC(X, cond_eta, theta = theta,
                family=family, penalty=penalty, method = "projection",
                penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  PhcTime <- proc.time() - time

  #stepwise
  time <- proc.time()
  Pstep <- WPSW(X, Y = cond_eta, theta, force=1, p=2,
               direction = "backward", method = "projection",
               transport.method = transport.method,
               display.progress = FALSE)
  PstepTime <- proc.time() - time
  # trajStep <- step$theta

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
  cat(paste0(anneal$message,"\n"))
  if(anneal$message != "completed") {
    annealTime <- paste0(">", annealTime)
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
                               method = ,
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
                               method = ,
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
    }

    rm(inSampModels)
    rm(inSampModelsProj)

    #### new X variable ####
    #mse on new outcome data from same paramters and different X
    #new method
    # augDatN <- augPseudo(X_new, cond_mu_new, theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatN, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("Running methods, new X variable: ", date(),"\n"))
    ipN <- W2IP(X = X_new, Y = cond_eta_new, theta = theta,
                display.progress=FALSE,
                transport.method = transport.method,
                model.size = ip_seq,
                infimum.maxit = 100, solution.method = "lp",
                parallel = NULL)

    lassoSelN <- W2L1(X_new, cond_eta_new, theta, family="gaussian",
                      penalty="selection.lasso",
                      penalty.factor=penalty_factN, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=100,
                      maxit = 1e5,
                      transport.method = transport.method,
                      display.progress=TRUE, gamma = 1, method = "selection.variable")
    # trajSelN <- extractCoef(lassoSelN)

    #carvalho method
    lassoHCN <- HC(X_new, cond_eta_new, theta = theta,
                   family=family, penalty=penalty,
                   penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, maxit = 1e5)

    #stepwise
    stepN <- WPSW(X_new, cond_eta_new, theta, force=1, p=2,
                  direction = "backward", method = "selection.variable",
                  transport.method = transport.method,
                  display.progress = TRUE)

    annealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=theta,
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

    #projection
    lassoProjN <- W2L1(X_new, cond_eta_new, theta, family="gaussian", penalty=penalty,
                       penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                       lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                       maxit=1e5,
                       transport.method = transport.method,
                       display.progress=TRUE, method = "projection")

    #stepwise
    PstepN <- WPSW(X_new, cond_eta_new, theta, force=1, p=2,
                  direction = "backward", method = "projection",
                  transport.method = transport.method,
                  display.progress = TRUE)
    #HC
    PlassoHCN <- HC(X_new, cond_eta_new, theta = theta,
                   family=family, penalty=penalty, method = "projection",
                   penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
    # anneal
    PannealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=theta,
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
      W2_newX <- distCompare(newXModels, target = list(posterior = theta,
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
      PW2_newX <- distCompare(newXModelsP, target = list(posterior = theta,
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



    #new method, single datapoint
    # augDatO <- augPseudo(X_sing, cond_mu_sing, theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatO, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("Running methods, single data point: ", date(),"\n"))
    ipO <- W2IP(X = X_sing, Y = cond_eta_sing, theta = theta,
               display.progress=FALSE,
               transport.method = transport.method,
               model.size = ip_seq,
               infimum.maxit = 10, solution.method = "lp",
               parallel = NULL)

    lassoSelO <- W2L1(X_sing, cond_eta_sing, theta, family="gaussian", penalty="selection.lasso",
                      penalty.factor=penalty_factO, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=100,
                      maxit=1e5,
                      transport.method = transport.method,
                      display.progress=TRUE, gamma = 1, method = "selection.variable")


    #stepwise
    stepO <- WPSW(X_sing, cond_eta_sing, theta, force=1, p=2,
                  direction = "backward",
                  method = "selection.variable",
                  transport.method = transport.method,
                  display.progress = TRUE)

    #simulated annealing
    annealO <- WPSA( X = X_sing, Y = cond_eta_sing, theta = theta,
                     force = 1, p=2, model.size = sa_seq, iter = SAiter, temps = SAtemps,
                     options = list(method = "selection.variable",
                                    energy.distribution = "boltzman",
                                    transport.method = transport.method,
                                    cooling.schedule="exponential",
                                    proposal.method = sa_prop),
                     display.progress = TRUE , max.time = sa_max_time)
    cat(paste0(annealO$message,"\n"))
    cat("\n")

    singleModels <- list("Binary Programming" = ipO,
                         "Lasso" = lassoSelO,
                         "Simulated Annealing" = annealO,
                         "Stepwise" = stepO#,
    )
    rm("ipO", "lassoSelO", "annealO","stepO")

    cat("Calculating distances\n")
    if( calc_w2_post){
      W2_single <- distCompare(singleModels, target = list(posterior = theta,
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
    }

    rm(singleModels)

    # trajAnnealO <- annealCoef(annealO, theta)

    # list of models
    # inSampModels <- list("Selection" = lassoSel,
    #                      "Simulated Annealing" = anneal,
    #                      "Stepwise" = step,
    #                      "Projection" = lassoProj,
    #                      "Hahn-Carvalho" = lassoHC)
    # newXModels <- list("Selection" = lassoSelN,
    #                    "Simulated Annealing" = annealN,
    #                    "Stepwise" = stepN,
    #                    "Projection" = lassoProjN,
    #                    "Hahn-Carvalho" = lassoHCN)
    # singleModels <- list("Selection" = lassoSelO,
    #                      "Simulated Annealing" = annealO,
    #                      "Stepwise" = stepO#,
    #                      # "Projection" = lassoProjO,
    #                      # "Hahn-Carvalho" = lassoHCO
    # )
    #
    # #w2 from "true" posterior mean
    # cat("\nCalculating W2 distances")
    # if( calc_w2_post){
    #   W2_insamp <- distCompare(inSampModels, target = list(posterior = theta,
    #                                                        mean = cond_mu),
    #                            method = "sinkhorn",
    #                            quantity=c("posterior","mean"),
    #                            parallel=NULL,
    #                            transform = data$invlink)
    #   W2_newX <- distCompare(newXModels, target = list(posterior = theta,
    #                                                    mean = cond_mu_new),
    #                          method = "sinkhorn",
    #                          quantity=c("posterior","mean"),
    #                          parallel=NULL,
    #                          transform = data$invlink)
    #   W2_single <- distCompare(singleModels, target = list(posterior = theta,
    #                                                        mean = cond_mu_sing),
    #                            method = "sinkhorn",
    #                            quantity=c("posterior","mean"),
    #                            parallel=NULL,
    #                            transform = data$invlink)
    # } else {
    #   W2_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
    #                                                        mean = cond_mu),
    #                            method = "sinkhorn",
    #                            quantity=c("mean"),
    #                            parallel=NULL,
    #                            transform = data$invlink)
    #   W2_newX <- distCompare(newXModels, target = list(posterior = NULL,
    #                                                    mean = cond_mu_new),
    #                          method = "sinkhorn",
    #                          quantity=c("mean"),
    #                          parallel=NULL,
    #                          transform = data$invlink)
    #   W2_single <- distCompare(singleModels, target = list(posterior = NULL,
    #                                                        mean = cond_mu_sing),
    #                            method = "sinkhorn",
    #                            quantity=c("mean"),
    #                            parallel=NULL,
    #                            transform = data$invlink)
    # }
    #
    # #mse from true mean
    # cat("\nCalculating MSEs")
    # mse_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
    #                                                       mean = true_mu),
    #                           method = "mse",
    #                           quantity="mean",
    #                           parallel=NULL,
    #                           transform = data$invlink)
    # mse_newX <- distCompare(newXModels, target = list(posterior = NULL,
    #                                                   mean = new_mu),
    #                         method = "mse",
    #                         quantity="mean",
    #                         parallel=NULL,
    #                         transform = data$invlink)
    # mse_single <- distCompare(singleModels, target = list(posterior = NULL,
    #                                                       mean = new_mu_sing),
    #                           method = "mse",
    #                           quantity="mean",
    #                           parallel=NULL,
    #                           transform = data$invlink)

    #mse on means of original data
    # mean_mseSel <- mse_idx_dist(trajSel$coefs, X, theta, true_mu, trajSel$nzero,p)
    # mean_mseSelE <- mse_idx_expect(trajSel$coefs, X, theta, true_mu, trajSel$nzero, p)
    # mean_msePerm <- mse_idx_dist(trajPerm$coefs, X, theta, true_mu, trajPerm$nzero,p)
    # mean_msePermE <- mse_idx_expect(trajPerm$coefs, X, theta, true_mu, trajPerm$nzero,p)
    # mean_mseHC <- mse_idx_dist(trajHCdist$coefs, X, theta, true_mu, trajHCdist$nzero,p)
    # mean_mseHCE <- mse_idx(trajHC$coefs, X, true_mu, trajHC$nzero, p)
    # mean_mseStep <- mse_idx_dist(trajStep$coefs, X, theta, true_mu, trajStep$nzero,p)
    # mean_mseAnneal <- mse_idx_dist(trajAnneal$coefs, X, theta, true_mu, trajAnneal$nzero,p)


    #mse on new outcome data from same paramters and same X
    # new_Y <- target$rdata(n, X, c(param$theta),
    #                          param$sigma2)
    # newY_mseSel <- mse_idx_dist(trajSel$coefs, X, theta, new_Y, trajSel$nzero, p)
    # newY_mseSelE <- mse_idx_expect(trajSel$coefs, X, theta, new_Y, trajSel$nzero, p)
    # newY_msePerm <- mse_idx_dist(trajPerm$coefs, X, theta, new_Y, trajPerm$nzero,p)
    # newY_msePermE <- mse_idx_expect(trajPerm$coefs, X, theta, new_Y, trajPerm$nzero,p)
    # newY_mseHC <- mse_idx_dist(trajHCdist$coefs, X, theta, new_Y, trajHCdist$nzero,p)
    # newY_mseHCE <- mse_idx(trajHC$coefs, X, new_Y, trajHC$nzero, p)
    # newY_mseStep <- mse_idx_dist( trajStep$coefs, X, theta, new_Y, trajStep$nzero, p )
    # newY_mseAnneal <- mse_idx_dist( trajAnneal$coefs, X, theta, new_Y, trajAnneal$nzero, p )




    #new mean
    # newMu_mseSel <- mse_idx_dist(trajSelN$coefs, X_new, theta, new_mu, trajSelN$nzero, p)
    # newMu_mseSelE <- mse_idx_expect(trajSelN$coefs, X_new, theta, new_mu, trajSelN$nzero, p)
    # newMu_msePerm <- mse_idx_dist(trajPermN$coefs, X_new, theta, new_mu, trajPermN$nzero,p)
    # newMu_msePermE <- mse_idx_expect(trajPermN$coefs, X_new, theta, new_mu, trajPermN$nzero,p)
    # newMu_mseHC <- mse_idx_dist(trajHCNdist$coefs, X_new, theta, new_mu, trajHCNdist$nzero,p)
    # newMu_mseHCE <- mse_idx(trajHCN$coefs, X_new, new_mu, trajHCN$nzero, p)
    # newMu_mseStep <- mse_idx_dist(trajStepN$coefs, X_new, theta, new_mu, trajStepN$nzero, p)
    # newMu_mseAnneal <- mse_idx_dist(trajAnnealN$coefs, X_new, theta, new_mu, trajAnnealN$nzero, p)

    #new data
    # newX_mseSel <- mse_idx_dist(trajSelN$coefs, X_new, theta, new_Y_new_X, trajSelN$nzero, p)
    # newX_mseSelE <- mse_idx_expect(trajSelN$coefs, X_new, theta, new_Y_new_X, trajSelN$nzero, p)
    # newX_msePerm <- mse_idx_dist(trajPermN$coefs, X_new, theta, new_Y_new_X, trajPermN$nzero,p)
    # newX_msePermE <- mse_idx_expect(trajPermN$coefs, X_new, theta, new_Y_new_X, trajPermN$nzero,p)
    # newX_mseHC <-  mse_idx(trajHCNdist$coefs, X_new, new_Y_new_X, trajHCNdist$nzero, p)
    # newX_mseHCE <-  mse_idx(trajHCN$coefs, X_new, new_Y_new_X, trajHCN$nzero, p)
    # newX_mseStep <- mse_idx_dist(trajStepN$coefs, X_new, theta, new_Y_new_X, trajStepN$nzero, p)
    # newX_mseAnneal <- mse_idx_dist(trajAnnealN$coefs, X_new, theta, new_Y_new_X, trajAnnealN$nzero, p)


    #single new mean
    # singleMu_mseSel <- mse_idx_dist(trajSelO$coefs, X_sing, theta, new_mu_sing, trajSelO$nzero, p)
    # singleMu_mseSelE <- mse_idx_expect(trajSelO$coefs, X_sing, theta, new_mu_sing, trajSelO$nzero, p)
    # singleMu_msePerm <- mse_idx_dist(trajPermO$coefs, X_sing, theta, new_mu_sing, trajPermO$nzero,p)
    # singleMu_msePermE <- mse_idx_expect(trajPermO$coefs, X_sing, theta, new_mu_sing, trajPermO$nzero,p)
    # singleMu_mseHC <- mse_idx_dist(trajHCOdist$coefs, X_sing, theta, new_mu_sing, trajHCOdist$nzero,p)
    # singleMu_mseHCE <- mse_idx(trajHCO$coefs, X_sing, new_mu_sing, trajHCO$nzero, p)
    # singleMu_mseStep <- mse_idx_dist(trajStepO$coefs, X_sing, theta, new_mu_sing, trajStepO$nzero, p)
    # singleMu_mseAnneal <- mse_idx_dist(trajAnnealO$coefs, X_sing, theta, new_mu_sing, trajAnnealO$nzero, p)

    #single new data
    # singleY_mseSel <- mse_idx_dist(trajSelO$coefs, X_sing, theta, single_Y, trajSelO$nzero, p)
    # singleY_mseSelE <- mse_idx_expect(trajSelO$coefs, X_sing, theta, single_Y, trajSelO$nzero, p)
    # singleY_msePerm <- mse_idx_dist(trajPermO$coefs, X_sing, theta, single_Y, trajPermO$nzero,p)
    # singleY_msePermE <- mse_idx_expect(trajPermO$coefs, X_sing, theta, single_Y, trajPermO$nzero,p)
    # singleY_mseHC <- mse_idx_dist(trajHCOdist$coefs, X_sing, theta, single_Y, trajHCOdist$nzero,p)
    # singleY_mseHCE <-  mse_idx(trajHCO$coefs, X_sing, single_Y, trajHCO$nzero, p)
    # singleY_mseStep <- mse_idx_dist(trajStepO$coefs, X_sing, theta, single_Y, trajStepO$nzero, p)
    # singleY_mseAnneal <- mse_idx_dist(trajAnnealO$coefs, X_sing, theta, single_Y, trajAnnealO$nzero, p)

    #### Using W2 to measure closeness to full posterior ####
    # X_ident <- diag(1,p,p)
    # postW2_Sel <- W2_idx(trajSel$coefs, X_ident, theta, post_sample$theta, trajSel$nzero,p)
    # postW2_Perm <- W2_idx(trajPerm$coefs, X_ident, theta, post_sample$theta, trajPerm$nzero,p)
    # postW2_Step <- W2_idx(trajStep$coefs, X_ident, theta, post_sample$theta, trajStep$nzero,p)
    # postW2_Anneal <- W2_idx(trajAnneal$coefs, X_ident, theta, post_sample$theta, trajAnneal$nzero,p)
    # postW2_HC <- W2_idx(trajHCdist$coefs, X_ident, theta, post_sample$theta, trajHCdist$nzero,p)

  } else {
    W2_insamp <- W2_newX <- W2_single <- mse_insamp <-  mse_newX <-  mse_single <- NULL
  }

  #### generate output list ####

  outList <- list (
    W2_dist = list(inSamp = list(selection = W2_insamp,
                                 projection = PW2_insamp),
                   newX = list(selection = W2_newX,
                               projection = PW2_newX),
                   single = list(selection = W2_single,
                                 projection = NULL) ),
    mse = list(inSamp = list(selection = mse_insamp,
                             projection = Pmse_insamp),
               newX = list(selection = mse_newX,
                           projection = Pmse_newX),
               single = list(selection = mse_single,
                             projection = NULL) ),
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
