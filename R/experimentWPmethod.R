experimentWPMethod <- function(target, hyperparameters, conditions, w2=FALSE) {
  n <- conditions$n
  p <- conditions$p

  n.samps <- conditions$n.samps
  penalty <- conditions$penalty
  lambda.min.ratio <- conditions$lambda.min.ratio
  n.lambda <- conditions$n.lambda
  penalty_method <- conditions$penalty.factor
  family <- conditions$family
  # pseudo.obs <- conditions$pseudo.obs
  stan_dir <- conditions$stan_dir
  calc_w2_post <- conditions$calc_w2_post
  if(is.null(calc_w2_post)) calc_w2_post <- FALSE
  # if(is.null(pseudo.obs)) pseudo.obs <- 0
  pseudo.obs <- 0
  posterior.method <- conditions$posterior.method
  transport.method <- conditions$transport.method
  if(is.null(transport.method)) transport.method <- "sinkhorn" #"univariate.approximation.pwr"
  not.only.timing <- conditions$not.only.timing
  if(is.null(not.only.timing)) not.only.timing <- FALSE

  sa_seq <- sort(unique(c(2,5,floor(seq(ceiling(p/5),p,floor(p/5))))))
  sa_max_time <- 64800
  sa_prop <- "random"
  # FSAiter <- 10*1:ceiling(p/2)
  # RSAiter <- if( p %% 2) { rev(FSAiter)[-1] } else { rev(FSAiter) }
  # SAiter <- c(FSAiter, RSAiter)
  SAiter <- 10*ceiling(p/2)
  SAtemps <- 50
  # SAiter <- 1
  # SAtemps <- 1

  #IP sequence
  if( p < 200 ) {
    ip_seq <- 1:p
  } else {
    ip_seq <- sa_seq
  }

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
  t_theta <- t(theta)
  sigma <- post_interp$sigma #variance (if it exists for model)

  #functions of theta
  E_theta <- rowMeans(theta)
  # theta_norm <- rowSums(theta^2)

  #penalty terms
  penalty_fact <- set_penalty_factor(theta, penalty_method)
  proj_penalty_fact <- set_penalty_factor(theta, "distance")
  HC_penalty_fact <- set_penalty_factor(theta, "expectation")

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

  #optional Wasserstein
  # if(w2){
  #   w2s <- w2l0(X,theta)
  # }

  # augDat <- augPseudo(X, cond_eta, t_theta, theta_norm, pseudo.obs, n, same=TRUE)
  # lambdas <- calc.lambdas(augDat, lambda.min.ratio, penalty_fact, n.lambda)
  cat(paste0("\nRunning methods, same data: ", date()))

  #### In sample ####
  #IP
  time <- proc.time()
  ip <- W2IP(X = X, Y = cond_eta, theta = t_theta,
                   display.progress=FALSE,
                   transport.method = transport.method,
                   model.size = ip_seq,
                   infimum.maxit = 100, solution.method = "cone",
                   parallel = NULL)
  ipTime <- proc.time() - time
  # trajSel <- selDist$theta

  #selection variable
  time <- proc.time()
  lassoSel <- W2L1(X, cond_eta, t_theta, family="gaussian", penalty="selection.lasso",
                   penalty.factor=penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                   maxit = 1e5,
                   display.progress=FALSE,
                   transport.method = transport.method,
                   gamma = 1, method = "selection.variable")
  selTime <- proc.time() - time
  # trajSel <- selDist$theta

  #projection
  time <- proc.time()
  lassoProj <- W2L1(X, cond_eta, t_theta, family="gaussian", penalty=penalty,
                    penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                    maxit = 1e5,
                    display.progress=FALSE,
                    transport.method = transport.method,
                    gamma = 1, method = "projection")
  projTime <- proc.time() - time
  # trajProj <- projDist$theta


  #carvalho method
  time <- proc.time()
  lassoHC <- HC(X, marg_eta, theta = t_theta,
                family=family, penalty=penalty,
                penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  hcTime <- proc.time() - time
  # trajHC <- lassoHC$theta

  # trajHCdist <- list(coefs = NULL, nzero = trajHC$nzero)
  # trajHCdist$coefs<- matrix(0, nrow=p, ncol = length(trajHC$nzero))
  # HC_non_zero_idx <- which(trajHC$coefs != 0)
  # trajHCdist$coefs[HC_non_zero_idx] <- 1

  #stepwise
  time <- proc.time()
  step <- WPSW(X, Y = cond_eta, t_theta, force=1, p=2,
               direction = "backward", method = "selection.variable",
               transport.method = transport.method,
               display.progress = FALSE)
  stepTime <- proc.time() - time
  # trajStep <- step$theta

  #simulated annealing
  annealTime <- NULL
  # if(n > 512 | p > 11){
  #   anneal <- WPSA(X=X, Y=cond_eta, theta=t_theta,
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
  anneal <- WPSA(X=X, Y=cond_eta, theta=t_theta,
                 force = 1, p=2, model.size = sa_seq, iter = SAiter, temps = SAtemps,
                 options = list(method = "selection.variable",
                                energy.distribution = "boltzman",
                                transport.method = transport.method,
                                cooling.schedule="exponential",
                                proposal.method = sa_prop),
                 display.progress=FALSE, max.time = sa_max_time)
  annealTime <- proc.time() - time
  # }
  # trajAnneal <- anneal$theta
  cat(anneal$message)
  if(anneal$message != "completed") {
    annealTime <- paste0(">", annealTime)
  }

  if (not.only.timing) {
    inSampModels <- list("I.P." = ip,
                         "Selection" = lassoSel,
                         "Simulated Annealing" = anneal,
                         "Stepwise" = step,
                         "Projection" = lassoProj,
                         "H.C." = lassoHC)
    cat("\nCalculating distances")
    if( calc_w2_post){
      W2_insamp <- distCompare(inSampModels, target = list(posterior = theta,
                                                           mean = cond_mu),
                               method = "exact",
                               quantity=c("posterior","mean"),
                               parallel=FALSE,
                               transform = data$invlink)
      mse_insamp <- distCompare(inSampModels, target = list(posterior = full_param,
                                                            mean = true_mu),
                                method = "mse",
                                quantity=c("posterior","mean"),
                                parallel=FALSE,
                                transform = data$invlink)
    }
    else {
      W2_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
                                                           mean = cond_mu),
                               method = "exact",
                               quantity=c("mean"),
                               parallel=FALSE,
                               transform = data$invlink)
      mse_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
                                                            mean = true_mu),
                                method = "mse",
                                quantity="mean",
                                parallel=FALSE,
                                transform = data$invlink)
    }

    rm(inSampModels)
    rm("ip","lassoSel", "anneal","step","lassoProj","lassoHC")

    #### new X variable ####
    #mse on new outcome data from same paramters and different X
    #new method
    # augDatN <- augPseudo(X_new, cond_mu_new, t_theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatN, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("\nRunning methods, new X variable: ", date()))
    ipN <- W2IP(X = X_new, Y = cond_eta_new, theta = t_theta,
                display.progress=FALSE,
                transport.method = transport.method,
                model.size = ip_seq,
                infimum.maxit = 100, solution.method = "cone",
                parallel = NULL)

    lassoSelN <- W2L1(X_new, cond_eta_new, t_theta, family="gaussian",
                      penalty="selection.lasso",
                      penalty.factor=penalty_fact, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                      maxit = 1e5,
                      transport.method = transport.method,
                      display.progress=TRUE, gamma = 1, method = "selection.variable")
    # trajSelN <- extractCoef(lassoSelN)

    #carvalho method
    lassoHCN <- HC(X_new, marg_eta_new, theta = t_theta,
                   family=family, penalty=penalty,
                   penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
    # trajHCN <- extractCoef(lassoHCN)
    #
    # trajHCNdist <- list(coefs = NULL, nzero = trajHCN$nzero)
    # trajHCNdist$coefs<- matrix(0, nrow=p, ncol = length(trajHCN$nzero))
    # HCN_non_zero_idx <- which(trajHCN$coefs != 0)
    # trajHCNdist$coefs[HCN_non_zero_idx] <- 1

    #permutation
    lassoProjN <- W2L1(X_new, cond_eta_new, t_theta, family="gaussian", penalty=penalty,
                       penalty.factor=proj_penalty_fact, nlambda = n.lambda,
                       lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
                       maxit=1e5,
                       transport.method = transport.method,
                       display.progress=TRUE, method = "projection")
    # trajPermN <- extractCoef(projDistN)

    #stepwise
    stepN <- WPSW(X_new, cond_eta_new, t_theta, force=1, p=2,
                  direction = "backward", method = "selection.variable",
                  transport.method = transport.method,
                  display.progress = TRUE)
    # trajStepN <- stepCoef(stepN, t_theta)

    #simulated annealing
    # if(n > 512 | p > 11){
    #   annealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=t_theta,
    #                    force = 1, p=2, model.size = 5, iter = SAiter,
    #                    temps = SAtemps,
    #                    options = list(method = "selection.variable",
    #                                   energy.distribution = "boltzman",
    #                                   transport.method = transport.method,
    #                                   cooling.schedule="exponential"),
    #                    display.progress = TRUE)
    # } else {
    annealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=t_theta,
                     force = 1, p=2, model.size = sa_seq, iter = SAiter,
                     temps = SAtemps,
                     options = list(method = "selection.variable",
                                    energy.distribution = "boltzman",
                                    transport.method = transport.method,
                                    cooling.schedule="exponential",
                                    proposal.method = sa_prop),
                     display.progress = TRUE, max.time = sa_max_time)
    cat(annealN$message)
    cat("\n")
    # }
    # trajAnnealN <- annealCoef(annealN, t_theta)
    newXModels <- list("I.P." = ipN,
                       "Selection" = lassoSelN,
                       "Simulated Annealing" = annealN,
                       "Stepwise" = stepN,
                       "Projection" = lassoProjN,
                       "H.C." = lassoHCN)

    cat("\nCalculating distances")
    if( calc_w2_post){
      W2_newX <- distCompare(newXModels, target = list(posterior = theta,
                                                       mean = cond_mu_new),
                             method = "exact",
                             quantity=c("posterior","mean"),
                             parallel=FALSE,
                             transform = data$invlink)
      mse_newX <- distCompare(newXModels, target = list(posterior = full_param,
                                                        mean = new_mu),
                              method = "mse",
                              quantity=c("posterior","mean"),
                              parallel=FALSE,
                              transform = data$invlink)
    }
    else {
      W2_newX <- distCompare(newXModels, target = list(posterior = NULL,
                                                       mean = cond_mu_new),
                             method = "exact",
                             quantity=c("mean"),
                             parallel=FALSE,
                             transform = data$invlink)
      mse_newX <- distCompare(newXModels, target = list(posterior = NULL,
                                                        mean = new_mu),
                              method = "mse",
                              quantity="mean",
                              parallel=FALSE,
                              transform = data$invlink)
    }


    rm(newXModels)
    rm("ipN","lassoSelN", "annealN","stepN","lassoProjN","lassoHCN")


    #new method, single datapoint
    # augDatO <- augPseudo(X_sing, cond_mu_sing, t_theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatO, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("\nRunning methods, single data point: ", date()))
    ipO <- W2IP(X = X_sing, Y = cond_eta_sing, theta = t_theta,
               display.progress=FALSE,
               transport.method = transport.method,
               model.size = ip_seq,
               infimum.maxit = 100, solution.method = "cone",
               parallel = NULL)

    lassoSelO <- W2L1(X_sing, cond_eta_sing, t_theta, family="gaussian", penalty="selection.lasso",
                      penalty.factor=penalty_fact, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                      maxit=1e5,
                      transport.method = transport.method,
                      display.progress=TRUE, gamma = 1, method = "selection.variable")
    # trajSelO <- extractCoef(lassoSelO)

    #carvalho method, single datapoint
    # lassoHCO <- HC(X, NULL, theta,
    #                     family=family, penalty=penalty,
    #                     penalty.factor=penalty_fact, nlambda = n.lambda,
    #                     lambda.min.ratio = lambda.min.ratio, lambda=lambdas)
    # trajHCO <- extractCoef(lassoHCO)

    # trajHCOdist <- list(coefs = NULL, nzero = trajHCO$nzero)
    # trajHCOdist$coefs<- matrix(0, nrow=p, ncol = length(trajHCO$nzero))
    # HCO_non_zero_idx <- which(trajHCO$coefs != 0)
    # trajHCOdist$coefs[HCO_non_zero_idx] <- 1

    #permutation, single datapoint
    # lassoProjO <- W2L1(X_sing, NULL, theta, family="gaussian", penalty=penalty,
    #                   penalty.factor=penalty_fact, nlambda = n.lambda,
    #                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
    #                   maxit=1000,
    #                   pseudo_observations = pseudo.obs,
    #                   display.progress=TRUE, method = "projection")
    # trajPermO <- extractCoef(permDistO)

    #stepwise
    stepO <- WPSW(X_sing, cond_eta_sing, t_theta, force=1, p=2,
                  direction = "backward",
                  method = "selection.variable",
                  transport.method = transport.method,
                  display.progress = TRUE)
    # trajStepO <- stepCoef(stepO, t_theta)

    #simulated annealing
    # if(n > 512 | p > 11){
    #   annealO <- WPSA( X = X_sing, Y = cond_eta_sing, theta = t_theta,
    #                    force = 1, p=2, model.size = 5, iter = SAiter, temps = SAtemps,
    #                    options = list(method = "selection.variable",
    #                                   energy.distribution = "boltzman",
    #                                   transport.method = transport.method,
    #                                   cooling.schedule="exponential"),
    #                    display.progress = TRUE )
    # } else {
    annealO <- WPSA( X = X_sing, Y = cond_eta_sing, theta = t_theta,
                     force = 1, p=2, model.size = sa_seq, iter = SAiter, temps = SAtemps,
                     options = list(method = "selection.variable",
                                    energy.distribution = "boltzman",
                                    transport.method = transport.method,
                                    cooling.schedule="exponential",
                                    proposal.method = sa_prop),
                     display.progress = TRUE , max.time = sa_max_time)
    cat(annealO$message)
    cat("\n")
    # }
    singleModels <- list("I.P." = ipO,
                         "Selection" = lassoSelO,
                         "Simulated Annealing" = annealO,
                         "Stepwise" = stepO#,
                         # "Projection" = lassoProjO,
                         # "H.C." = lassoHCO
    )

    cat("\nCalculating distances")
    if( calc_w2_post){
      W2_single <- distCompare(singleModels, target = list(posterior = theta,
                                                           mean = cond_mu_sing),
                               method = "exact",
                               quantity=c("posterior","mean"),
                               parallel=FALSE,
                               transform = data$invlink)
      mse_single <- distCompare(singleModels, target = list(posterior = full_param,
                                                            mean = new_mu_sing),
                                method = "mse",
                                quantity=c("posterior","mean"),
                                parallel=FALSE,
                                transform = data$invlink)
    }
    else {

      W2_single <- distCompare(singleModels, target = list(posterior = NULL,
                                                           mean = cond_mu_sing),
                               method = "exact",
                               quantity=c("mean"),
                               parallel=FALSE,
                               transform = data$invlink)
      mse_single <- distCompare(singleModels, target = list(posterior = NULL,
                                                            mean = new_mu_sing),
                                method = "mse",
                                quantity="mean",
                                parallel=FALSE,
                                transform = data$invlink)
    }

    rm(singleModels)
    rm("ipO", "lassoSelO", "annealO","stepO")
    # trajAnnealO <- annealCoef(annealO, t_theta)

    # list of models
    # inSampModels <- list("Selection" = lassoSel,
    #                      "Simulated Annealing" = anneal,
    #                      "Stepwise" = step,
    #                      "Projection" = lassoProj,
    #                      "H.C." = lassoHC)
    # newXModels <- list("Selection" = lassoSelN,
    #                    "Simulated Annealing" = annealN,
    #                    "Stepwise" = stepN,
    #                    "Projection" = lassoProjN,
    #                    "H.C." = lassoHCN)
    # singleModels <- list("Selection" = lassoSelO,
    #                      "Simulated Annealing" = annealO,
    #                      "Stepwise" = stepO#,
    #                      # "Projection" = lassoProjO,
    #                      # "H.C." = lassoHCO
    # )
    #
    # #w2 from "true" posterior mean
    # cat("\nCalculating W2 distances")
    # if( calc_w2_post){
    #   W2_insamp <- distCompare(inSampModels, target = list(posterior = theta,
    #                                                        mean = cond_mu),
    #                            method = "sinkhorn",
    #                            quantity=c("posterior","mean"),
    #                            parallel=FALSE,
    #                            transform = data$invlink)
    #   W2_newX <- distCompare(newXModels, target = list(posterior = theta,
    #                                                    mean = cond_mu_new),
    #                          method = "sinkhorn",
    #                          quantity=c("posterior","mean"),
    #                          parallel=FALSE,
    #                          transform = data$invlink)
    #   W2_single <- distCompare(singleModels, target = list(posterior = theta,
    #                                                        mean = cond_mu_sing),
    #                            method = "sinkhorn",
    #                            quantity=c("posterior","mean"),
    #                            parallel=FALSE,
    #                            transform = data$invlink)
    # } else {
    #   W2_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
    #                                                        mean = cond_mu),
    #                            method = "sinkhorn",
    #                            quantity=c("mean"),
    #                            parallel=FALSE,
    #                            transform = data$invlink)
    #   W2_newX <- distCompare(newXModels, target = list(posterior = NULL,
    #                                                    mean = cond_mu_new),
    #                          method = "sinkhorn",
    #                          quantity=c("mean"),
    #                          parallel=FALSE,
    #                          transform = data$invlink)
    #   W2_single <- distCompare(singleModels, target = list(posterior = NULL,
    #                                                        mean = cond_mu_sing),
    #                            method = "sinkhorn",
    #                            quantity=c("mean"),
    #                            parallel=FALSE,
    #                            transform = data$invlink)
    # }
    #
    # #mse from true mean
    # cat("\nCalculating MSEs")
    # mse_insamp <- distCompare(inSampModels, target = list(posterior = NULL,
    #                                                       mean = true_mu),
    #                           method = "mse",
    #                           quantity="mean",
    #                           parallel=FALSE,
    #                           transform = data$invlink)
    # mse_newX <- distCompare(newXModels, target = list(posterior = NULL,
    #                                                   mean = new_mu),
    #                         method = "mse",
    #                         quantity="mean",
    #                         parallel=FALSE,
    #                         transform = data$invlink)
    # mse_single <- distCompare(singleModels, target = list(posterior = NULL,
    #                                                       mean = new_mu_sing),
    #                           method = "mse",
    #                           quantity="mean",
    #                           parallel=FALSE,
    #                           transform = data$invlink)

    #mse on means of original data
    # mean_mseSel <- mse_idx_dist(trajSel$coefs, X, t_theta, true_mu, trajSel$nzero,p)
    # mean_mseSelE <- mse_idx_expect(trajSel$coefs, X, t_theta, true_mu, trajSel$nzero, p)
    # mean_msePerm <- mse_idx_dist(trajPerm$coefs, X, t_theta, true_mu, trajPerm$nzero,p)
    # mean_msePermE <- mse_idx_expect(trajPerm$coefs, X, t_theta, true_mu, trajPerm$nzero,p)
    # mean_mseHC <- mse_idx_dist(trajHCdist$coefs, X, t_theta, true_mu, trajHCdist$nzero,p)
    # mean_mseHCE <- mse_idx(trajHC$coefs, X, true_mu, trajHC$nzero, p)
    # mean_mseStep <- mse_idx_dist(trajStep$coefs, X, t_theta, true_mu, trajStep$nzero,p)
    # mean_mseAnneal <- mse_idx_dist(trajAnneal$coefs, X, t_theta, true_mu, trajAnneal$nzero,p)


    #mse on new outcome data from same paramters and same X
    # new_Y <- target$rdata(n, X, c(param$theta),
    #                          param$sigma2)
    # newY_mseSel <- mse_idx_dist(trajSel$coefs, X, t_theta, new_Y, trajSel$nzero, p)
    # newY_mseSelE <- mse_idx_expect(trajSel$coefs, X, t_theta, new_Y, trajSel$nzero, p)
    # newY_msePerm <- mse_idx_dist(trajPerm$coefs, X, t_theta, new_Y, trajPerm$nzero,p)
    # newY_msePermE <- mse_idx_expect(trajPerm$coefs, X, t_theta, new_Y, trajPerm$nzero,p)
    # newY_mseHC <- mse_idx_dist(trajHCdist$coefs, X, t_theta, new_Y, trajHCdist$nzero,p)
    # newY_mseHCE <- mse_idx(trajHC$coefs, X, new_Y, trajHC$nzero, p)
    # newY_mseStep <- mse_idx_dist( trajStep$coefs, X, t_theta, new_Y, trajStep$nzero, p )
    # newY_mseAnneal <- mse_idx_dist( trajAnneal$coefs, X, t_theta, new_Y, trajAnneal$nzero, p )




    #new mean
    # newMu_mseSel <- mse_idx_dist(trajSelN$coefs, X_new, t_theta, new_mu, trajSelN$nzero, p)
    # newMu_mseSelE <- mse_idx_expect(trajSelN$coefs, X_new, t_theta, new_mu, trajSelN$nzero, p)
    # newMu_msePerm <- mse_idx_dist(trajPermN$coefs, X_new, t_theta, new_mu, trajPermN$nzero,p)
    # newMu_msePermE <- mse_idx_expect(trajPermN$coefs, X_new, t_theta, new_mu, trajPermN$nzero,p)
    # newMu_mseHC <- mse_idx_dist(trajHCNdist$coefs, X_new, t_theta, new_mu, trajHCNdist$nzero,p)
    # newMu_mseHCE <- mse_idx(trajHCN$coefs, X_new, new_mu, trajHCN$nzero, p)
    # newMu_mseStep <- mse_idx_dist(trajStepN$coefs, X_new, t_theta, new_mu, trajStepN$nzero, p)
    # newMu_mseAnneal <- mse_idx_dist(trajAnnealN$coefs, X_new, t_theta, new_mu, trajAnnealN$nzero, p)

    #new data
    # newX_mseSel <- mse_idx_dist(trajSelN$coefs, X_new, t_theta, new_Y_new_X, trajSelN$nzero, p)
    # newX_mseSelE <- mse_idx_expect(trajSelN$coefs, X_new, t_theta, new_Y_new_X, trajSelN$nzero, p)
    # newX_msePerm <- mse_idx_dist(trajPermN$coefs, X_new, t_theta, new_Y_new_X, trajPermN$nzero,p)
    # newX_msePermE <- mse_idx_expect(trajPermN$coefs, X_new, t_theta, new_Y_new_X, trajPermN$nzero,p)
    # newX_mseHC <-  mse_idx(trajHCNdist$coefs, X_new, new_Y_new_X, trajHCNdist$nzero, p)
    # newX_mseHCE <-  mse_idx(trajHCN$coefs, X_new, new_Y_new_X, trajHCN$nzero, p)
    # newX_mseStep <- mse_idx_dist(trajStepN$coefs, X_new, t_theta, new_Y_new_X, trajStepN$nzero, p)
    # newX_mseAnneal <- mse_idx_dist(trajAnnealN$coefs, X_new, t_theta, new_Y_new_X, trajAnnealN$nzero, p)


    #single new mean
    # singleMu_mseSel <- mse_idx_dist(trajSelO$coefs, X_sing, t_theta, new_mu_sing, trajSelO$nzero, p)
    # singleMu_mseSelE <- mse_idx_expect(trajSelO$coefs, X_sing, t_theta, new_mu_sing, trajSelO$nzero, p)
    # singleMu_msePerm <- mse_idx_dist(trajPermO$coefs, X_sing, t_theta, new_mu_sing, trajPermO$nzero,p)
    # singleMu_msePermE <- mse_idx_expect(trajPermO$coefs, X_sing, t_theta, new_mu_sing, trajPermO$nzero,p)
    # singleMu_mseHC <- mse_idx_dist(trajHCOdist$coefs, X_sing, t_theta, new_mu_sing, trajHCOdist$nzero,p)
    # singleMu_mseHCE <- mse_idx(trajHCO$coefs, X_sing, new_mu_sing, trajHCO$nzero, p)
    # singleMu_mseStep <- mse_idx_dist(trajStepO$coefs, X_sing, t_theta, new_mu_sing, trajStepO$nzero, p)
    # singleMu_mseAnneal <- mse_idx_dist(trajAnnealO$coefs, X_sing, t_theta, new_mu_sing, trajAnnealO$nzero, p)

    #single new data
    # singleY_mseSel <- mse_idx_dist(trajSelO$coefs, X_sing, t_theta, single_Y, trajSelO$nzero, p)
    # singleY_mseSelE <- mse_idx_expect(trajSelO$coefs, X_sing, t_theta, single_Y, trajSelO$nzero, p)
    # singleY_msePerm <- mse_idx_dist(trajPermO$coefs, X_sing, t_theta, single_Y, trajPermO$nzero,p)
    # singleY_msePermE <- mse_idx_expect(trajPermO$coefs, X_sing, t_theta, single_Y, trajPermO$nzero,p)
    # singleY_mseHC <- mse_idx_dist(trajHCOdist$coefs, X_sing, t_theta, single_Y, trajHCOdist$nzero,p)
    # singleY_mseHCE <-  mse_idx(trajHCO$coefs, X_sing, single_Y, trajHCO$nzero, p)
    # singleY_mseStep <- mse_idx_dist(trajStepO$coefs, X_sing, t_theta, single_Y, trajStepO$nzero, p)
    # singleY_mseAnneal <- mse_idx_dist(trajAnnealO$coefs, X_sing, t_theta, single_Y, trajAnnealO$nzero, p)

    #### Using W2 to measure closeness to full posterior ####
    # X_ident <- diag(1,p,p)
    # postW2_Sel <- W2_idx(trajSel$coefs, X_ident, t_theta, post_sample$theta, trajSel$nzero,p)
    # postW2_Perm <- W2_idx(trajPerm$coefs, X_ident, t_theta, post_sample$theta, trajPerm$nzero,p)
    # postW2_Step <- W2_idx(trajStep$coefs, X_ident, t_theta, post_sample$theta, trajStep$nzero,p)
    # postW2_Anneal <- W2_idx(trajAnneal$coefs, X_ident, t_theta, post_sample$theta, trajAnneal$nzero,p)
    # postW2_HC <- W2_idx(trajHCdist$coefs, X_ident, t_theta, post_sample$theta, trajHCdist$nzero,p)

  } else {
    W2_insamp <- W2_newX <- W2_single <- mse_insamp <-  mse_newX <-  mse_single <- NULL
  }

  #### generate output list ####

  outList <- list (
    W2_dist = list(inSamp = W2_insamp, newX = W2_newX, single = W2_single),
    mse = list(inSamp = mse_insamp, newX = mse_newX, single = mse_single),
    # mseNewY = list(sel=newY_mseDSel, sel_E = newY_mseSelE, perm = newY_msePerm,
    #                perm_E = newY_msePermE, hc=newY_mseHC, step = newY_mseStep),
    # mseNewMu = list(sel = newMu_mseSel, sel_E = newMu_mseSelE, perm = newMu_msePerm,
    # perm_E = newMu_msePermE, hc=newMu_mseHC, step = newMu_mseStep),
    # mseNewX = list(sel=newX_mseSel, sel_E = newX_mseSelE, perm = newX_msePerm,
    #                perm_E = newX_msePermE, hc=newX_mseHC, step = newX_mseStep),
    # singleNewMu = list(sel = singleMu_mseSel, dist_E = singleMu_mseSelE, perm = singleMu_msePerm,
    # perm_E = singleMu_msePermE, hc=singleMu_mseHC, step = singleMu_mseStep),
    time = list(ip = ipTime[3], selection = selTime[3],
                projection = projTime[3], HC = hcTime[3],
                step = stepTime[3], anneal = annealTime[3])#,
    # W2_dist_post = list(sel = postW2_Sel, perm = postW2_Perm,
    #                     step = postW2_Step, anneal = postW2_Anneal,
    #                     HC = postW2_HC)
  )

  #add w2 data if necessary
  # if( w2 ) {
  #   W2_w2 <- sapply(w2s$minCombPerActive, function(i)
  #     mse(apply(cond_mu,1,sort), apply(X[,i] %*% t(theta[,i]),1,sort)))
  #   mean_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(true_mu, X[,i] %*% t(theta[,i])))
  #   newY_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(new_Y, X[,i] %*% t(theta[,i])))
  #   newMu_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(new_mu, X_new[,i] %*% t(theta[,i])))
  #   newX_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(new_Y_new_X, X_new[,i] %*% t(theta[,i])))
  #
  #   #save in output
  #   outList$W2_dist$w2  <- W2_w2
  #   outList$mseMean$w2  <- mean_w2
  #   outList$mseNewY$w2  <- newY_w2
  #   outList$mseNewMu$w2 <- newMu_w2
  #   outList$mseNewX$w2  <- newX_w2
  # }
  # ranks <- lapply(outList, function(o) apply(sapply(o, function(cc) cc), 1, rank))
  # for(ii in seq_along(names(outList))) {
  #   ranks[[ii]][which(is.na(t(sapply(outList[[ii]], function(cc) cc))),arr.ind=TRUE)] <- NA
  # }
  # outList$rank <- lapply(ranks, rowMeans, na.rm=TRUE)
  return(outList)
}
