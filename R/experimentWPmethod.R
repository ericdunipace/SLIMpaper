#' Simulation study function
#'
#' @param target function that is the output of the get_binary_nonlinear, get_normal_linear, etc
#' @param hyperparameters hyperparameters for the Bayesian methodology
#' @param conditions conditions for the simulation
#'
#' @return
#' @export
experimentWPMethod <- function(target, hyperparameters, conditions) {
  n <- conditions$n
  p <- conditions$p

  n.samps <- conditions$n.samps
  penalty <- conditions$penalty
  lambda.min.ratio <- conditions$lambda.min.ratio
  n.lambda <- conditions$n.lambda
  penalty_method <- conditions$penalty.factor
  family  <- conditions$family
  solver <- conditions$solver
  python.path <- conditions$python.path
  recalc <- conditions$recalculate
  methods.to.run <- conditions$methods.to.run

  if (is.null(solver) ) {
    solver <- "mosek"
  }
  if (solver == "") {
    solver <- "mosek"
  }
  solver <- match.arg(solver, choices = c("mosek", "gurobi",
                                          "cplex", "cone", "lp"))
  if(is.null(recalc)) recalc <- FALSE
  if(recalc == "") recalc <- FALSE

  if (is.null(methods.to.run)) methods.to.run <- c("approximate binary program", "binary program", "L1", "simulated annealing", "miscellaneous")

  # pseudo.obs <- conditions$pseudo.obs
  stan_dir <- conditions$stan_dir
  calc_w2_post <- conditions$calc_w2_post
  if(is.null(calc_w2_post)) calc_w2_post <- FALSE
  # if(is.null(pseudo.obs)) pseudo.obs <- 0
  pseudo.obs <- 0
  posterior.method <- conditions$posterior.method
  transport.method <- conditions$transport.method
  if(is.null(transport.method)) transport.method <- "hilbert"
  epsilon <- 0.05
  otmaxit <- switch(transport.method, "exact" = 0, 100)
  #"univariate.approximation.pwr"
  not.only.timing <- conditions$not.only.timing
  if(is.null(not.only.timing)) not.only.timing <- FALSE
  L0 <- conditions$L0
  if(is.null(L0)) L0 <- FALSE

  sa_max_time <- 64800
  sa_prop <- "random"
  # FSAiter <- 10*1:ceiling(p/2)
  # RSAiter <- if( p %% 2) { rev(FSAiter)[-1] } else { rev(FSAiter) }
  # SAiter <- c(FSAiter, RSAiter)
  # SAiter <- 10 * ceiling(p_star/2)
  # SAtemps <- 50
  # SAiter <- 1
  # SAtemps <- 1

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
  p_star <- if(family != "gaussian") {
    min(length(param$theta),p)
  } else {
    min(length(param$theta),p + 4)
  }

  if (p > p_star) {
    p_temp <- p
  } else {
    p_temp <- p_star
  }

  sa_seq <- sort(unique(c(2,5,floor(seq(ceiling(p_temp/5),p_temp,floor(p_temp/5))))))
  SAiter <- 10 * ceiling(p_temp/2)
  SAtemps <- 50

  #IP sequence
  if( p_temp < 30 ) {
    ip_seq <- 1:p_temp
  } else {
    ip_seq <- sort(unique(c(1,floor(seq(ceiling(p_temp/20),41,length.out = 20)))))
  }

  param$theta <- param$theta[1:p_star]
  full_param <- c(param$theta, rep(0, max(p-p_star,0)))

  #set up data
  X <- target$X$rX(n, target$X$corr, p)
  X_sing <- matrix(c(target$X$rX(1, target$X$corr, p)), nrow=1, ncol=p)
  X_new <- target$X$rX(n, target$X$corr, p)
  X_neighborhood <- cbind(1, SLIMpaper::rmvnorm(nsamples = p*10,
                                                    mean = X_sing[,-1,drop=FALSE],
                                                    covariance = cov(X[,-1,drop=FALSE])/n))


  data <- target$rdata(n, X, c(param$theta),
                       param$sigma2, method = "modified.friedman", corr = target$X$corr)
  single_data <- target$rdata(1, X_sing, c(param$theta),
                              param$sigma2, method = "modified.friedman", corr = target$X$corr)
  new_data <- target$rdata(n, X_new, c(param$theta), param$sigma2, method = "modified.friedman", corr = target$X$corr)

  neighb_data <- target$rdata(n, X_neighborhood, c(param$theta), param$sigma2, method = "modified.friedman", corr = target$X$corr)

  Y <- data$Y
  single_Y <- single_data$Y
  new_Y_new_X <- new_data$Y

  #True means
  true_mu <- data$mu
  new_mu <- new_data$mu
  new_mu_sing <- single_data$mu
  new_mu_neighb <- neighb_data$mu

  # indices for different data
  single_idx <- 1
  new_idx <- 2:(n+1)
  neighb_idx <- (n+2):(n+1 + p*10)

  #derivative store
  derivative.test <- NA

  #sample theta
  if(family != "binomial") {
    X <- data$model_matrix(X,p_star)
    X_sing <- data$model_matrix(X_sing,p_star)
    X_new <- data$model_matrix(X_new, p_star)
    X_neighborhood <- data$model_matrix(X_neighborhood, p_star)

    if(is.null(hyperparameters$Lambda) || dim(hyperparameters$Lambda)[1] < p_star) {
      hyperparameters$Lambda <- matrix(1, p_star, p_star)
    }
    if (is.null(hyperparameters$mu) || length(hyperparameters$mu) < p_star) {
      hyperparameters$mu <- rep(0, p_star)
    }

    post_sample <- target$rpost(n.samps, X,
                                Y, hyperparameters,
                                method = posterior.method, stan_dir = stan_dir,
                                X.test = rbind(X_sing, X_new, X_neighborhood),
                                chains = 1, m0 = 20,
                                is.exponential = TRUE,
                                python.path = python.path)

    post_interp <- post_sample
    theta <- theta_new <- theta_sing <- post_interp$theta #regression coef
    if(nrow(theta) != ncol(X)) theta <- t(theta)
    if(!is.null(theta)) {
      t_theta <- t(theta)
    } else {
      t_theta <- theta
    }
    sigma <- post_interp$sigma

    derivative.test <- theta

    pred.orig <- deriv.orig <- NULL

    # X <- target$model_matrix(X)
    # X_sing <- target$model_matrix(X_sing)

  } else {
    hyperparameters <- list(learning.rate = 1e-3,
                            lambda = 0.001,
                            batch.size = 128,
                            first.layer.width = as.integer(p*10),
                            hidden.layer.width = as.integer(2/3 * p + 1))
    post_sample <- target$rpost(n.samps, X,
                                Y,
                                hyperparameters,
                                method = posterior.method,
                                stan_dir = stan_dir,
                                X.test = rbind(X_sing, X_new, X_neighborhood),
                                chains = 1, m0 = 20,
                                is.exponential = TRUE,
                                test.portion = 0,
                                learning.rate = 1e-3,
                                niter = 5,
                                python.path = python.path,
                                verbose = FALSE)
    # post_interp <- target$rpost(n.samps, X, Y, NULL, method = "logistic", stan_dir = "exec/Stan/logistic_horseshoe_noQR.stan",
    #                             X.test = rbind(X_sing, X_new, X_neighborhood),
    #                             chains = 1, m0 = 20)
    # add non-linearities if binomial
    # if(family == "binomial") {
      # cov_X <- cov(X[,-1,drop = FALSE])
      # form <- formula("~.**2")
    X <- data$model_matrix(X)
    X_sing <- data$model_matrix(X_sing)
    X_new <- data$model_matrix(X_new)
    X_neighborhood <- data$model_matrix(X_neighborhood)

    # }
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
    derivative.test <- sapply(post_sample$test$derivatives, function(d) d[1,])
    pred.orig <- c(post_sample$yhat$test)[c(single_idx,neighb_idx)]
    deriv.orig <- post_sample$yhat$derivatives[c(single_idx,neighb_idx),]
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
    L0list$Selection <- WpProj:::WPL0(X = X, Y = cond_eta, theta = theta,
                             power = 2,  method = "selection.variable",
                             transport.method = transport.method, epsilon = epsilon,
                             maxit = otmaxit)
    L0list$Projection <- WpProj:::WPL0(X = X, Y = cond_eta, theta = theta,
                              power = 2,  method = "projection",
                              transport.method = transport.method, epsilon = epsilon,
                              maxit = otmaxit)

    cat("L0 Distance Calculations\n")
    W2L0 <- WpProj:::distCompare(L0list,
                        target = list(parameters = NULL,
                                      predictions = cond_mu),
                        method = wp_alg,
                        quantity=c("predictions"),
                        parallel=NULL,
                        transform = data$invlink,
                        epsilon = epsilon,
                        niter = otmaxit)
    mseL0 <- WpProj:::distCompare(L0list,
                         target = list(parameters = NULL,
                                       predictions = true_mu),
                         method = "mse",
                         quantity="predictions",
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

  time <- list(NULL)

  PW2_insamp    <- W2_insamp  <-
    PW2_newX    <- W2_newX    <-
    PW2_single  <- W2_single  <-
    mse_insamp  <-  mse_newX  <-  mse_single  <-
    Pmse_insamp <-  Pmse_newX <-  Pmse_single <- NULL
  imp_insamp    <- imp_newX   <- imp_single   <- NULL
  w2_r2_single  <- Pw2_r2_single              <- NULL
  w2_r2_single_null <- Pw2_r2_single_null <- w1_r2_single <- NULL
  Pw1_r2_single <- w1_r2_single_null <- w1_r2_single <-Pw1_r2_single_null <- NULL
  singleModels <- W2_insamp <- singleModelsP <- W2_insamp <- NULL

  if (not.only.timing) {
    #### new method, single datapoint ####
    # augDatO <- augPseudo(X_sing, cond_mu_sing, theta, theta_norm, pseudo.obs, n, same=TRUE)
    # lambdas <- calc.lambdas(augDatO, lambda.min.ratio, penalty_fact, n.lambda)
    cat(paste0("Running methods, single data point: ", date(),"\n"))
    if (!(family == "binomial" && posterior.method == "nn") ) {
      cat("  Selection: IP\n")
      # ipO <- WpProj:::W2IP(X = X_sing, Y = cond_eta_sing,
      #                      theta = theta_sing,
      #             display.progress=TRUE,
      #             transport.method = transport.method,
      #             nvars = ip_seq, maxit = 1e4,
      #             infimum.maxit = 100, solver = solver,
      #             parallel = NULL)
      ipO <- WpProj(X = X_sing, eta = cond_eta_sing,
                    theta = theta_sing, power = 2,
                    method = "binary program",
                    solver = solver,
                    options = binary_program_method_options(
                      maxit = 1e4, infimum.maxit = 100,
                      nvars = ip_seq,
                      display.progress=TRUE
                    ))

      cat(" Lagrange BP\n")
      ipLO <- LBP(X = X_sing, eta = cond_eta_sing, theta = theta_sing,
                 solver = solver,
                 nvars = ip_seq
      )

      cat(" Lasso\n")
      # lassoSelO <- WpProj:::W2L1(X_sing, cond_eta_sing, theta_sing, penalty=penalty,
      #                   penalty.factor=penalty_factO, nlambda = n.lambda,
      #                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=100,
      #                   maxit=1e5, alpha = 0.99, gamma = 1.1,
      #                   transport.method = transport.method,
      #                   display.progress=TRUE, method = "selection.variable")
      lassoSelO <- WpProj(X = X_sing, eta = cond_eta_sing,
                          theta = theta_sing, power = 2,
                          method = "binary program",
                          solver = "lasso",
                          options = binary_program_method_options(
                            maxit=1e5,
                            infimum.maxit=100,
                            transport.method = transport.method,
                            display.progress=TRUE,
                            solver.options = L1_method_options(penalty = penalty,
                                                               nlambda = n.lambda,
                                                               lambda.min.ratio = lambda.min.ratio,
                                                               alpha = 0.99, gamma = 1.1,
                                                               solver.options = list(penalty.factor=penalty_factO))
                          ))
      cat(" SW\n")
      #stepwise
      # stepO <- WpProj:::WPSW(X_sing, cond_eta_sing, theta_sing, force=1, power=2,
      #               direction = "backward",
      #               method = "selection.variable",
      #               transport.method = transport.method,
      #               display.progress = TRUE)
      stepO <- WpProj(X = X_sing, eta = cond_eta_sing,
                      theta = theta_sing, power=2,
                      method = "stepwise",
                      options = stepwise_method_options(
                        force=1,
                        direction = "backward",
                        method = "binary program",
                        transport.method = transport.method,
                        display.progress = TRUE))

      cat(" SA\n")
      #simulated annealing
      # annealO <- WpProj:::WPSA( X = X_sing, Y = cond_eta_sing, theta = theta_sing,
      #                  force = 1, power=2, nvars = sa_seq, maxit = SAiter, temps = SAtemps,
      #                  options = list(method = "selection.variable",
      #                                 energy.distribution = "boltzman",
      #                                 transport.method = transport.method,
      #                                 cooling.schedule="exponential",
      #                                 proposal.method = sa_prop),
      #                  display.progress = TRUE , max.time = sa_max_time)
      annealO <- WpProj( X = X_sing, eta = cond_eta_sing,
                         theta = theta_sing,power=2,
                         method = "simulated annealing",
                         solver = "lasso",
                         options = simulated_annealing_method_options(
                           method = "binary program",
                           force = 1, nvars = sa_seq,
                           maxit = SAiter, temps = SAtemps,
                           energy.distribution = "boltzman",
                           transport.method = transport.method,
                           cooling.schedule = "exponential",
                           proposal.method = switch(sa_prop, "random" = "uniform", sa_prop),
                           display.progress = TRUE , max.time = sa_max_time)
      )
      cat(paste0(annealO$message,"\n"))
      cat("\n")
    }


    cat("  Projection: Lasso, ",date(),"\n")
    #projection
    # lassoProjO <- WpProj:::W2L1(X_neighborhood, cond_eta_neighb, theta_sing, penalty=penalty,
    #                    penalty.factor=proj_penalty_fact, nlambda = n.lambda,
    #                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
    #                    maxit=1e6, alpha = 0.99, gamma = 1.1,
    #                    transport.method = transport.method,
    #                    display.progress=TRUE, method = "projection")
    lassoProjO <- WpProj(X = X_neighborhood, eta = cond_eta_neighb,
                         theta = theta_sing,
                         method = "L1",
                         solver = "lasso",
                         options = L1_method_options(
                           nlambda = n.lambda,
                           lambda.min.ratio = lambda.min.ratio,
                           maxit=1e6, alpha = 0.99, gamma = 1.1,
                           display.progress = TRUE,
                           solver.options = list(penalty.factor = proj_penalty_fact
                                                 )
                         ))
    cat(" SW, ",date(),"\n")
    #stepwise
    # PstepO <- WpProj:::WPSW(X_neighborhood, cond_eta_neighb, theta_sing, force=1, power=2,
    #                direction = "backward", method = "projection",
    #                transport.method = transport.method,
    #                display.progress = TRUE)
    PstepO <- WpProj(X = X_neighborhood, eta = cond_eta_neighb,
                     theta = theta_sing, power=2,
                     method = "stepwise",
                     options = stepwise_method_options(
                       force=1,
                       direction = "backward",
                       method = "projection",
                       display.progress = TRUE))
    cat(" HC, ",date(),"\n")
    #HC
    PlassoHCO <- HC(X_neighborhood, cond_eta_neighb, theta = theta_sing, alpha = 0.99, gamma = 1.1,
                    family="gaussian", penalty=penalty, method = "projection",
                    penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
    cat(" SA, ",date(),"\n")
    # anneal
    # PannealO <-  WpProj:::WPSA(X_neighborhood, cond_eta_neighb, theta=theta_sing,
    #                   force = 1, power=2, nvars = sa_seq, maxit = SAiter,
    #                   temps = SAtemps,
    #                   options = list(method = "projection",
    #                                  energy.distribution = "boltzman",
    #                                  transport.method = transport.method,
    #                                  cooling.schedule="exponential",
    #                                  proposal.method = sa_prop),
    #                   display.progress = TRUE, max.time = sa_max_time)
    PannealO <-  WpProj( X = X_neighborhood, eta = cond_eta_neighb,
                         theta = theta_sing, power=2,
                         method = "simulated annealing",
                         solver = "lasso",
                         options = simulated_annealing_method_options(
                           method = "projection",
                           force = 1, nvars = sa_seq,
                           maxit = SAiter, temps = SAtemps,
                           energy.distribution = "boltzman",
                           transport.method = transport.method,
                           cooling.schedule = "exponential",
                           proposal.method = switch(sa_prop, "random" = "uniform", sa_prop),
                           display.progress = TRUE , max.time = sa_max_time)
    )
    # PL1O <- NULL
    cat("\n L1, ",date(),"\n")
    # PL1O <- WpProj:::W1L1(X=X_neighborhood, Y=cond_eta_neighb,
    #              solver = solver,
    #              nlambda = n.lambda, penalty = penalty,
    #              lambda.min.ratio = lambda.min.ratio, maxit = 1e5,
    #              gamma = 1.1, display.progress = TRUE)
    PL1O <- WpProj(X = X_neighborhood, eta = cond_eta_neighb,
                   power = 1,
                   solver = solver,
                   method = "L1",
                   options = L1_method_options(
                     nlambda = n.lambda, penalty = penalty,
                     lambda.min.ratio = lambda.min.ratio, maxit = 1e5,
                     gamma = 1.1, display.progress = TRUE
                   ))
    # cat("\n L3\n")
    # PL3O <- WPL1(X=X_neighborhood, Y=cond_eta_neighb, p = 3,
    #                  nlambda = n.lambda, penalty = penalty,
    #                  lambda.min.ratio = lambda.min.ratio,
    #                  gamma = 1.1,
    #                  display.progress = TRUE)
    # PLInfO <- NULL
    cat("\n LInfinity, ",date(),"\n")
    # PLInfO <- WpProj:::WInfL1(X=X_neighborhood, Y=cond_eta_neighb,
    #                  nlambda = n.lambda, penalty = penalty,
    #                  lambda.min.ratio = lambda.min.ratio,
    #                  gamma = 1.1, solver = solver,
    #                  display.progress = TRUE)
    PLInfO <- WpProj(X = X_neighborhood, eta = cond_eta_neighb,
           power = Inf,
           solver = solver,
           method = "L1",
           options = L1_method_options(
             nlambda = n.lambda, penalty = penalty,
             lambda.min.ratio = lambda.min.ratio, maxit = 1e5,
             gamma = 1.1, display.progress = TRUE
           ))

    cat("\n")
    cat(paste0(PannealO$message))
    cat("\n")
    # }
    # }
    # trajAnnealN <- annealCoef(annealN, theta)
    if (!(family == "binomial" && posterior.method == "nn") ) {
    singleModels <- list("Binary Programming" = ipO,
                         "Lagrange Binary Programming" = ipLO,
                         "Lasso" = lassoSelO,
                         "Simulated Annealing" = annealO,
                         "Stepwise" = stepO#,
    )
    }
    else {
      singleModels <- NULL
    }
    # rm("ipO", "lassoSelO", "annealO","stepO")
    singleModelsP <- list("Lasso" = lassoProjO,
                          "Simulated Annealing" = PannealO,
                          "Stepwise" = PstepO,
                          "Hahn-Carvalho" = PlassoHCO
                          , "L1" = PL1O
                          # "L3" = PL3O,
                          , "LInf" = PLInfO
                          )
    # rm("lassoProjO",
    #    "PannealO", "PstepO","PlassoHCO")
    #recalculate values for single obs
    if( recalc ) {
      singleModelsP <- lapply(singleModelsP, function(x) {
        x$eta <- lapply(x$theta, function(tt) X_sing %*% tt)
        return(x)
      })
      cond_mu_calc <- cond_mu_sing
      new_mu_calc <- new_mu_sing
    } else {
      cond_mu_calc <- cond_mu_neighb
      new_mu_calc <- new_mu_neighb
    }

    cat("Calculating distances\n")
    if( calc_w2_post){
      W2_single <- distCompare(singleModels, target = list(parameters = theta_sing,
                                                           predictions = cond_mu_sing),
                               method = wp_alg,
                               quantity=c("parameters","predictions"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit,
                               power = 2)
      W1_single <-  distCompare(singleModels, target = list(parameters = theta_sing,
                                                           predictions = cond_mu_sing),
                               method = wp_alg,
                               quantity=c("parameters","predictions"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit,
                               power = 1)
      # cat("W2 selection\n")
      w2_r2_single <-  WPR2(predictions = cond_mu_sing, projected_model = W2_single, p = 2, method = wp_alg)
      w1_r2_single <-  WPR2(predictions = cond_mu_sing, projected_model = W1_single, p = 1, method = wp_alg)
      w2_r2_single_null <-  WPR2(predictions = NULL, projected_model = W2_single, p = 2, method = wp_alg)
      w1_r2_single_null <-  WPR2(predictions = NULL, projected_model = W1_single, p = 1, method = wp_alg)

      mse_single <-  distCompare(singleModels, target = list(parameters = full_param,
                                                            predictions = new_mu_sing),
                                method = "mse",
                                quantity=c("parameters", "predictions"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      mse_single$predictions$dist <- mse_single$predictions$dist/mean((c(cond_mu_sing) - c(new_mu_sing))^2)
      mse_single$parameters$dist <- mse_single$parameters$dist/mean(((full_param) - c(theta))^2)
      # cat("W2 projection\n")
      PW2_single <-  distCompare(singleModelsP, target = list(parameters = theta_sing,
                                                             predictions = cond_mu_calc),
                                method = wp_alg,
                                quantity=c("parameters", "predictions"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      PW1_single <-  distCompare(singleModelsP, target = list(parameters = theta_sing,
                                                             predictions = cond_mu_calc),
                                method = wp_alg,
                                quantity=c("parameters","predictions"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit,
                                power = 1.0)
      Pw2_r2_single <-  WPR2(predictions = cond_mu_calc, projected_model = PW2_single, p = 2, method = wp_alg)
      Pw1_r2_single <-  WPR2(predictions = cond_mu_calc, projected_model = PW1_single, p = 1, method = wp_alg)
      if(recalc) {
        Pw2_r2_single_null <-  WPR2(predictions = NULL, projected_model = PW2_single, p = 2, method = wp_alg)
        Pw1_r2_single_null <-  WPR2(predictions = NULL, projected_model = PW1_single, p = 1, method = wp_alg)
      } else {
        Pw2_r2_single_null <-  WPR2(predictions = cond_mu_calc, projected_model = singleModelsP, p = 2, method = wp_alg)
        Pw1_r2_single_null <-  WPR2(predictions = cond_mu_calc, projected_model = singleModelsP, p = 1, method = wp_alg)
      }

      Pmse_single <-  distCompare(singleModelsP, target = list(parameters = full_param,
                                                              predictions = new_mu_calc),
                                 method = "mse",
                                 quantity=c("parameters","predictions"),
                                 parallel=NULL,
                                 transform = data$invlink,
                                 epsilon = epsilon,
                                 niter = otmaxit)
      Pmse_single$predictions$dist <- Pmse_single$predictions$dist/mean((cond_mu_calc - c(new_mu_calc))^2)
      Pmse_single$parameters$dist <- Pmse_single$parameters$dist/mean(((full_param) - c(theta))^2)

    }
    else {

      W2_single <-  W1_single <- w2_r2_single <-
        w1_r2_single <- w2_r2_single_null <-
        w1_r2_single_null <-mse_single <- NULL
      if (!(family == "binomial" && posterior.method == "nn") ) {
      W2_single <-  distCompare(singleModels, target = list(parameters = NULL,
                                                           predictions = cond_mu_sing),
                               method = wp_alg,
                               quantity=c("predictions"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit,
                               power = 2)
      W1_single <-  distCompare(singleModels, target = list(parameters = NULL,
                                                           predictions = cond_mu_sing),
                               method = wp_alg,
                               quantity=c("predictions"),
                               parallel=NULL,
                               transform = data$invlink,
                               epsilon = epsilon,
                               niter = otmaxit,
                               power = 1)
      # cat("W2 selection\n")
      w2_r2_single <-  WPR2(predictions = cond_mu_sing, projected_model = W2_single, p = 2, method = wp_alg)
      w1_r2_single <-  WPR2(predictions = cond_mu_sing, projected_model = W1_single, p = 1, method = wp_alg)
      w2_r2_single_null <-  WPR2(predictions = cond_mu_sing, projected_model = singleModels, p = 2, method = wp_alg, base = data$invlink(colMeans(cond_eta)))
      w1_r2_single_null <-  WPR2(predictions = cond_mu_sing, projected_model = singleModels, p = 1, method = wp_alg, base = data$invlink(colMeans(cond_eta)))

      mse_single <-  distCompare(singleModels,
                                 target = list(parameters = NULL,
                                               predictions = new_mu_sing),
                                method = "mse",
                                quantity="predictions",
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      mse_single$predictions$dist <- mse_single$predictions$dist/mean((c(cond_mu_sing) - c(new_mu_sing))^2)
      # mse_single$parameters$dist <- mse_single$parameters$dist/mean((as.matrix(full_param) - c(theta))^2)
      }
      # cat("W2 projection\n")
      PW2_single <- distCompare(singleModelsP, target = list(parameters = NULL,
                                                             predictions = cond_mu_calc),
                                method = wp_alg,
                                quantity=c("predictions"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit)
      PW1_single <- distCompare(singleModelsP, target = list(parameters = NULL,
                                                             predictions = cond_mu_calc),
                                method = wp_alg,
                                quantity=c("predictions"),
                                parallel=NULL,
                                transform = data$invlink,
                                epsilon = epsilon,
                                niter = otmaxit,
                                power = 1)
      Pw2_r2_single <- WPR2(predictions = cond_mu_calc, projected_model = PW2_single, p = 2, method = wp_alg)
      Pw1_r2_single <- WPR2(predictions = cond_mu_calc, projected_model = PW1_single, p = 1, method = wp_alg)
      if(recalc) {
        Pw2_r2_single_null <- WPR2(predictions = cond_mu_sing, projected_model = PW2_single, p = 2, method = wp_alg)
        Pw1_r2_single_null <- WPR2(predictions = cond_mu_sing, projected_model = PW1_single, p = 1, method = wp_alg)
      } else {
        Pw2_r2_single_null <- WPR2(predictions = cond_mu_calc, projected_model = singleModelsP, p = 2, method = wp_alg)
        Pw1_r2_single_null <- WPR2(predictions = cond_mu_calc, projected_model = singleModelsP, p = 1, method = wp_alg)
      }


      Pmse_single <- distCompare(singleModelsP, target = list(parameters = NULL,
                                                              predictions = new_mu_calc),
                                 method = "mse",
                                 quantity="predictions",
                                 parallel=NULL,
                                 transform = data$invlink,
                                 epsilon = epsilon,
                                 niter = otmaxit)
      Pmse_single$predictions$dist <- Pmse_single$predictions$dist/mean((cond_mu_calc - c(new_mu_calc))^2)
      # Pmse_single$parameters$dist <- Pmse_single$parameters$dist/mean((as.matrix(full_param) - c(theta))^2)

    }

    # rm(singleModels)

    #### Variable importance ####
    cat(paste0("Running variable importance: ", date(),"\n"))
    # pred.fun.theta <- switch(family,
    #                          "gaussian" = theta,
    #                          "binomial" = theta,
    #                          "exponential" = list(baseline = post_sample$mu$S$base,
    #                                               param = theta))
    # imp_newX <- WPVI(X = X_new, Y = NULL, theta = theta_new, pred.fun = pred.fun,
    #                  p = 2, ground_p = 2, transport.method = "exact")
    if (!(family == "binomial" && posterior.method == "nn") ) {
    imp_single <- WPVI(X = X_sing, eta = NULL, theta = theta_sing, pred.fun = pred.fun,
                       p = 2, ground_p = 2, transport.method = "exact")
    }

  }
  else {
    cat(paste0("Running methods, same data: ", date(),"\n"))

    #### In sample ####
    ipTime <- selTime <- ipLTime<- NULL
    if("binary program" %in% methods.to.run) {
      #IP
      cat("   Selection: IP,\n")
      time <- proc.time()
      ip <- WpProj(X = X, eta = cond_eta, theta = theta, power = 2,
                   method = "binary program", solver = solver,
                   options = binary_program_method_options(
                     infimum.maxit = 100, nvars = ip_seq,
                     transport.method = transport.method,
                     parallel = NULL, display.progress = FALSE
                   )
      )
      ipTime <- proc.time() - time
      # trajSel <- selDist$theta
    }
    if ("lagrange binary program" %in% methods.to.run) {
      cat("   Selection: Lagrange IP, ")
      #IP
      time <- proc.time()
      ipL <- LBP(X = X, eta = cond_eta, theta = theta,
                   solver = solver,
                 nvars = ip_seq
      )
      ipLTime <- proc.time() - time
    }
    if ("approximate binary program" %in% methods.to.run) {
      cat("   Selection: IP, Lasso, ")
      #selection variable
      time <- proc.time()
      lassoSel <- WpProj(X = X, eta = cond_eta ,
                 theta = theta , power = 2,
                 method = "binary program",
                 solver = "lasso",
                 options = binary_program_method_options(
                   maxit=1e5,
                   infimum.maxit=100,
                   transport.method = transport.method,
                   display.progress=FALSE,
                   solver.options = L1_method_options(
                    penalty = penalty,
                    nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio,
                    alpha = 0.99, gamma = 1.1,
                    solver.options = list(penalty.factor=penalty_factO))
  ))
      selTime <- proc.time() - time
    }

    # trajSel <- selDist$theta
    PstepTime <- stepTime <- hcTime <- PhcTime <- NULL
    if("miscellaneous" %in% methods.to.run) {
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
      step <- WpProj(X = X, eta = cond_eta, theta = theta,
                     method = "stepwise",
                     options = stepwise_method_options(
                       force = 1, direction = "backward", method = "binary program", transport.method = transport.method,
                       display.progress = FALSE
                     ))
      stepTime <- proc.time() - time

      cat(" HC, ")
      time <- proc.time()
      PlassoHC <-  HC(X, cond_eta, theta = theta, alpha = 0.99, gamma = 1.1,
                      family="gaussian", penalty=penalty, method = "projection",
                      penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
      PhcTime <- proc.time() - time

      cat(" SW, ")
      #stepwise
      time <- proc.time()
      Pstep <- WpProj(X = X, eta = cond_eta,
                      theta = theta, power=2,
                      method = "stepwise",
                      options = stepwise_method_options(
                        force=1,
                        direction = "backward",
                        method = "projection",
                        display.progress = FALSE))
      PstepTime <- proc.time() - time
    }


    # trajStep <- step$theta
    annealTime <- PannealTime <- NULL
    if("simulated annealing" %in% methods.to.run) {
      cat(" SA\n")
      #simulated annealing
      # if(n > 512 | p > 11){
      #   anneal <- WPSA(X=X, Y=cond_eta, theta=theta,
      #                  force = 1, p=2, nvars = 5, maxit = SAiter, temps = SAtemps,
      #                  options = list(method = "selection.variable",
      #                                 energy.distribution = "boltzman",
      #                                 transport.method = transport.method,
      #                                 cooling.schedule="exponential"),
      #                  display.progress=TRUE)
      #   annealTime <- NULL
      # }
      # else {
      time <- proc.time()
      anneal <- WpProj( X = X , eta = cond_eta ,
                        theta = theta, power=2,
                        method = "simulated annealing",
                        solver = "lasso",
                        options = simulated_annealing_method_options(
                          method = "binary program",
                          force = 1, nvars = sa_seq,
                          maxit = SAiter, temps = SAtemps,
                          energy.distribution = "boltzman",
                          transport.method = transport.method,
                          cooling.schedule = "exponential",
                          proposal.method = switch(sa_prop, "random" = "uniform", sa_prop),
                          display.progress = FALSE , max.time = sa_max_time))
      annealTime <- proc.time() - time

      cat(paste0(anneal$message,"\n"))
      if(anneal$message != "completed") {
        annealTime <- paste0(">", annealTime)
      }



      time <- proc.time()
      Panneal <- WpProj( X = X, eta = cond_eta ,
                         theta = theta , power=2,
                         method = "simulated annealing",
                         solver = "lasso",
                         options = simulated_annealing_method_options(
                           method = "projection",
                           force = 1, nvars = sa_seq,
                           maxit = SAiter, temps = SAtemps,
                           energy.distribution = "boltzman",
                           transport.method = transport.method,
                           cooling.schedule = "exponential",
                           proposal.method = switch(sa_prop, "random" = "uniform", sa_prop),
                           display.progress = FALSE , max.time = sa_max_time)
      )
      PannealTime <- proc.time() - time

      cat(paste0(Panneal$message,"\n"))
      if(Panneal$message != "completed") {
        PannealTime <- paste0(">", PannealTime)
      }

    }

    # anneal <- WpProj:::WPSA(X=X, Y=cond_eta, theta=theta,
    #                force = 1, power=2, nvars = sa_seq, maxit = SAiter, temps = SAtemps,
    #                options = list(method = "selection.variable",
    #                               energy.distribution = "boltzman",
    #                               transport.method = transport.method,
    #                               cooling.schedule="exponential",
    #                               proposal.method = sa_prop),
    #                display.progress=FALSE, max.time = sa_max_time)

    projTime <- PL1Time <- PLInfTime <- NULL
    if("L1" %in% methods.to.run) {
      cat("   Projection: Lasso, ")
      #projection
      time <- proc.time()
      lassoProj <- WpProj(X = X_neighborhood, eta = cond_eta_neighb,
             theta = theta_sing,
             method = "L1",
             solver = "lasso",
             options = L1_method_options(
               nlambda = n.lambda,
               lambda.min.ratio = lambda.min.ratio,
               maxit=1e6, alpha = 0.99, gamma = 1.1,
               display.progress = FALSE,
               solver.options = list(penalty.factor = proj_penalty_fact
               )
             ))
      projTime <- proc.time() - time
      # lassoProj <- WpProj:::W2L1(X, cond_eta, theta, penalty=penalty,
      #                   penalty.factor=proj_penalty_fact, nlambda = n.lambda,
      #                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
      #                   maxit = 1e6, alpha = 0.99, gamma = 1.1,
      #                   display.progress=FALSE,
      #                   transport.method = transport.method,
      #                   method = "projection")



      time <- proc.time()
      PL1 <- WpProj(X = X , eta = cond_eta ,
                    power = 1,
                    solver = solver,
                    method = "L1",
                    options = L1_method_options(
                      nlambda = n.lambda, penalty = penalty,
                      lambda.min.ratio = lambda.min.ratio, maxit = 1e5,
                      gamma = 1.1, display.progress = FALSE
                    ))
      PL1Time <- proc.time() - time

        time <- proc.time()
        PLInf <-WpProj(X = X, eta = cond_eta ,
                       power = Inf,
                       solver = solver,
                       method = "L1",
                       options = L1_method_options(
                         nlambda = n.lambda, penalty = penalty,
                         lambda.min.ratio = lambda.min.ratio, maxit = 1e5,
                         gamma = 1.1, display.progress = FALSE
                       ))
        PLInfTime <- proc.time() - time
    }

    time <- list(selection = list(ip = as.numeric(ipTime[3]),
                                  lbp = as.numeric(ipLTime[3]),
                                  lasso = as.numeric(selTime[3]),
                                  HC = as.numeric(hcTime[3]),
                                  step = as.numeric(stepTime[3]),
                                  anneal = as.numeric(annealTime[3])),
                 projection = list(lasso = as.numeric(projTime[3]),
                                   L1 = as.numeric(PL1Time[3]),
                                   # L3 = as.numeric(PL3Time[3]),
                                   LInf = as.numeric(PLInfTime[3]),
                                   HC = as.numeric(PhcTime[3]), step = as.numeric(PstepTime[3]),
                                   anneal = as.numeric(PannealTime[3])))


  }

    imp_insamp <- WpProj:::WPVI(X = X, eta = NULL, theta = theta, pred.fun = pred.fun,
                                p = 2, ground_p = 2, transport.method = "exact")
  # augDat <- augPseudo(X, cond_eta, theta, theta_norm, pseudo.obs, n, same=TRUE)
  # lambdas <- calc.lambdas(augDat, lambda.min.ratio, penalty_fact, n.lambda)

  # }
  # trajAnneal <- anneal$theta


  # if (not.only.timing) {
  #   inSampModels <- list("Binary Programming" = ip,
  #                        "L2" = lassoSel,
  #                        "Simulated Annealing" = anneal,
  #                        "Stepwise" = step,
  #                        "Hahn-Carvalho" = lassoHC)
  #   inSampModelsProj <- list("L2" = lassoProj,
  #                            "L1" = PL1,
  #                            "LInf" = PLInf,
  #                            "Simulated Annealing" = Panneal,
  #                            "Stepwise" = Pstep,
  #                            "Hahn-Carvalho" = PlassoHC)
  #   rm("ip","lassoSel", "anneal","step","lassoHC")
  #   rm("lassoProj","PL1", "PLInf", "Panneal","Pstep","PlassoHC")
  #
  #   cat("Calculating distances\n")
  #   if( calc_w2_post){
  #     W2_insamp <- distCompare(inSampModels, target = list(parameters = theta,
  #                                                          predictions = cond_mu),
  #                              method = wp_alg,
  #                              quantity=c("","predictions"),
  #                              parallel=NULL,
  #                              transform = data$invlink,
  #                              epsilon = epsilon,
  #                              niter = otmaxit)
  #     mse_insamp <- distCompare(inSampModels, target = list(parameters = full_param,
  #                                                           predictions = true_mu),
  #                               method = "mse",
  #                               quantity=c("parameters","predictions"),
  #                               parallel=NULL,
  #                               transform = data$invlink,
  #                               epsilon = epsilon,
  #                               niter = otmaxit)
  #     PW2_insamp <- distCompare(inSampModelsProj, target = list(parameters = theta,
  #                                                               predictions = cond_mu),
  #                               method = wp_alg,
  #                               quantity=c("parameters","predictions"),
  #                               parallel=NULL,
  #                               transform = data$invlink,
  #                               epsilon = epsilon,
  #                               niter = otmaxit)
  #     Pmse_insamp <- distCompare(inSampModelsProj, target = list(parameters = full_param,
  #                                                                predictions = true_mu),
  #                                method = "mse",
  #                                quantity=c("parameters","predictions"),
  #                                parallel=NULL,
  #                                transform = data$invlink,
  #                                epsilon = epsilon,
  #                                niter = otmaxit)
  #   }
  #   else {
  #     W2_insamp <- distCompare(inSampModels, target = list(parameters = NULL,
  #                                                          predictions = cond_mu),
  #                              method = wp_alg,
  #                              quantity=c("predictions"),
  #                              parallel=NULL,
  #                              transform = data$invlink,
  #                              epsilon = epsilon,
  #                              niter = otmaxit)
  #     mse_insamp <- distCompare(inSampModels, target = list(parameters = NULL,
  #                                                           predictions = true_mu),
  #                               method = "mse",
  #                               quantity="predictions",
  #                               parallel=NULL,
  #                               transform = data$invlink,
  #                               epsilon = epsilon,
  #                               niter = otmaxit)
  #     PW2_insamp <- distCompare(inSampModelsProj, target = list(parameters = NULL,
  #                                                               predictions = cond_mu),
  #                               method = wp_alg,
  #                               quantity=c("predictions"),
  #                               parallel=NULL,
  #                               transform = data$invlink,
  #                               epsilon = epsilon,
  #                               niter = otmaxit)
  #     Pmse_insamp <- distCompare(inSampModelsProj, target = list(parameters = NULL,
  #                                                                predictions = true_mu),
  #                                method = "mse",
  #                                quantity="predictions",
  #                                parallel=NULL,
  #                                transform = data$invlink,
  #                                epsilon = epsilon,
  #                                niter = otmaxit)
  #   }
  #
  #   rm(inSampModels)
  #   rm(inSampModelsProj)
  #
  #   #### new X variable ####
  #   #mse on new outcome data from same paramters and different X
  #   #new method
  #   # augDatN <- augPseudo(X_new, cond_mu_new, theta, theta_norm, pseudo.obs, n, same=TRUE)
  #   # lambdas <- calc.lambdas(augDatN, lambda.min.ratio, penalty_fact, n.lambda)
  #   cat(paste0("Running methods, new X variable: ", date(),"\n"))
  #   cat("  Selection: IP")
  #   ipN <- W2IP(X = X_new, Y = cond_eta_new, theta = theta_new,
  #               display.progress=TRUE,
  #               transport.method = transport.method,
  #               nvars = ip_seq,
  #               infimum.maxit = 100, solution.method = "cone",
  #               parallel = NULL)
  #
  #   cat(" Lasso")
  #   lassoSelN <- W2L1(X_new, cond_eta_new, theta_new, family="gaussian",
  #                     penalty=penalty, alpha = 0.99, gamma = 1.1,
  #                     penalty.factor=penalty_factN, nlambda = n.lambda,
  #                     lambda.min.ratio = lambda.min.ratio, infimum.maxit=100,
  #                     maxit = 1e5,
  #                     transport.method = transport.method,
  #                     display.progress=TRUE, method = "selection.variable")
  #   # trajSelN <- extractCoef(lassoSelN)
  #
  #   cat(" HC, ")
  #   #carvalho method
  #   lassoHCN <- HC(X_new, cond_eta_new, theta = theta_new,
  #                  family="gaussian", penalty=penalty, alpha = 0.99, gamma = 1.1,
  #                  penalty.factor=HC_penalty_fact, nlambda = n.lambda,
  #                  lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  #
  #   cat(" SW")
  #   #stepwise
  #   stepN <- WPSW(X_new, cond_eta_new, theta_new, force=1, p=2,
  #                 direction = "backward", method = "selection.variable",
  #                 transport.method = transport.method,
  #                 display.progress = TRUE)
  #
  #   cat(" SA")
  #   annealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=theta_new,
  #                    force = 1, p=2, nvars = sa_seq, maxit = SAiter,
  #                    temps = SAtemps,
  #                    options = list(method = "selection.variable",
  #                                   energy.distribution = "boltzman",
  #                                   transport.method = transport.method,
  #                                   cooling.schedule="exponential",
  #                                   proposal.method = sa_prop),
  #                    display.progress = TRUE, max.time = sa_max_time)
  #   cat(paste0(annealN$message,"\n"))
  #   cat("\n")
  #
  #   cat("  Projection: Lasso")
  #   #projection
  #   lassoProjN <- W2L1(X_new, cond_eta_new, theta_new, family="gaussian", penalty=penalty,
  #                      penalty.factor=proj_penalty_fact, nlambda = n.lambda,
  #                      lambda.min.ratio = lambda.min.ratio, infimum.maxit=1,
  #                      maxit=1e5, alpha = 0.99, gamma = 1.1,
  #                      transport.method = transport.method,
  #                      display.progress=TRUE, method = "projection")
  #
  #   cat(" SW")
  #   #stepwise
  #   PstepN <- WPSW(X_new, cond_eta_new, theta_new, force=1, p=2,
  #                  direction = "backward", method = "projection",
  #                  transport.method = transport.method,
  #                  display.progress = TRUE)
  #   cat(" HC, ")
  #   #HC
  #   PlassoHCN <- HC(X_new, cond_eta_new, theta = theta_new, alpha = 0.99, gamma = 1.1,
  #                   family="gaussian", penalty=penalty, method = "projection",
  #                   penalty.factor=HC_penalty_fact, nlambda = n.lambda,
  #                   lambda.min.ratio = lambda.min.ratio, maxit = 1e5)
  #   cat(" SA")
  #   # anneal
  #   PannealN <-  WPSA(X=X_new, Y=cond_eta_new, theta=theta_new,
  #                     force = 1, p=2, nvars = sa_seq, maxit = SAiter,
  #                     temps = SAtemps,
  #                     options = list(method = "projection",
  #                                    energy.distribution = "boltzman",
  #                                    transport.method = transport.method,
  #                                    cooling.schedule="exponential",
  #                                    proposal.method = sa_prop),
  #                     display.progress = TRUE, max.time = sa_max_time)
  #   cat(paste0(PannealN$message,"\n"))
  #   cat("\n")
  #   # }
  #   # trajAnnealN <- annealCoef(annealN, theta)
  #   newXModels <- list("Binary Programming" = ipN,
  #                      "Lasso" = lassoSelN,
  #                      "Simulated Annealing" = annealN,
  #                      "Stepwise" = stepN,
  #                      "Hahn-Carvalho" = lassoHCN)
  #   rm("ipN","lassoSelN", "annealN","stepN","lassoHCN")
  #   newXModelsP <- list("Lasso" = lassoProjN,
  #                       "Simulated Annealing" = PannealN,
  #                       "Stepwise" = PstepN,
  #                       "Hahn-Carvalho" = PlassoHCN)
  #   rm("lassoProjN",
  #      "PannealN", "PstepN","PlassoHCN")
  #
  #   cat("Calculating distances\n")
  #   if( calc_w2_post){
  #     W2_newX <- distCompare(newXModels, target = list(parameters = theta_new,
  #                                                      predictions = cond_mu_new),
  #                            method = wp_alg,
  #                            quantity=c("parameters","predictions"),
  #                            parallel=NULL,
  #                            transform = data$invlink,
  #                            epsilon = epsilon,
  #                            niter = otmaxit)
  #     mse_newX <- distCompare(newXModels, target = list(parameters = full_param,
  #                                                       predictions = new_mu),
  #                             method = "mse",
  #                             quantity=c("parameters","predictions"),
  #                             parallel=NULL,
  #                             transform = data$invlink,
  #                             epsilon = epsilon,
  #                             niter = otmaxit)
  #     PW2_newX <- distCompare(newXModelsP, target = list(parameters = theta_new,
  #                                                        predictions = cond_mu_new),
  #                             method = wp_alg,
  #                             quantity=c("parameters","predictions"),
  #                             parallel=NULL,
  #                             transform = data$invlink,
  #                             epsilon = epsilon,
  #                             niter = otmaxit)
  #     Pmse_newX <- distCompare(newXModelsP, target = list(parameters = full_param,
  #                                                         predictions = new_mu),
  #                              method = "mse",
  #                              quantity=c("parameters","predictions"),
  #                              parallel=NULL,
  #                              transform = data$invlink,
  #                              epsilon = epsilon,
  #                              niter = otmaxit)
  #   }
  #   else {
  #     W2_newX <- distCompare(newXModels, target = list(parameters = NULL,
  #                                                      predictions = cond_mu_new),
  #                            method = wp_alg,
  #                            quantity=c("predictions"),
  #                            parallel=NULL,
  #                            transform = data$invlink,
  #                            epsilon = epsilon,
  #                            niter = otmaxit)
  #     mse_newX <- distCompare(newXModels, target = list(parameters = NULL,
  #                                                       predictions = new_mu),
  #                             method = "mse",
  #                             quantity="predictions",
  #                             parallel=NULL,
  #                             transform = data$invlink,
  #                             epsilon = epsilon,
  #                             niter = otmaxit)
  #     PW2_newX <- distCompare(newXModelsP, target = list(parameters = NULL,
  #                                                        predictions = cond_mu_new),
  #                             method = wp_alg,
  #                             quantity=c("predictions"),
  #                             parallel=NULL,
  #                             transform = data$invlink,
  #                             epsilon = epsilon,
  #                             niter = otmaxit)
  #     Pmse_newX <- distCompare(newXModelsP, target = list(parameters = NULL,
  #                                                         predictions = new_mu),
  #                              method = "mse",
  #                              quantity="predictions",
  #                              parallel=NULL,
  #                              transform = data$invlink,
  #                              epsilon = epsilon,
  #                              niter = otmaxit)
  #   }
  #
  #
  #   rm(newXModels, newXModelsP)






  #### generate output list ####

  outList <- list (
    W2_dist = list(inSamp = list(selection = W2_insamp,
                                 projection = PW2_insamp),
                   # newX = list(selection = W2_newX,
                   #             projection = PW2_newX),
                   single = list(selection = W2_single,
                                 projection = PW2_single) ),
    mse = list(inSamp = list(selection = mse_insamp,
                             projection = Pmse_insamp),
               # newX = list(selection = mse_newX,
               #             projection = Pmse_newX),
               single = list(selection = mse_single,
                             projection = Pmse_single) ),
    importance = list(inSamp = imp_insamp,
                      # newX = imp_newX,
                      single = imp_single),
    W2_r2 = list(expectation = list(
                                   inSamp = list(selection = NULL,
                                                 projection = NULL),
                                   single = list(selection = w2_r2_single,
                                               projection = Pw2_r2_single) ),
                 null        = list(
                                   inSamp = list(selection = NULL,
                                                 projection = NULL),
                                   single = list(selection = w2_r2_single_null,
                                                 projection = Pw2_r2_single_null) )
                                 ),
    W1_r2 = list(expectation = list(
                                    inSamp = list(selection = NULL,
                                                  projection = NULL),
                                    single = list(selection = w1_r2_single,
                                                  projection = Pw1_r2_single) ),
                 null        = list(
                                    inSamp = list(selection = NULL,
                                                  projection = NULL),
                                    single = list(selection = w1_r2_single_null,
                                                  projection = Pw1_r2_single_null) )
              ),
    time = time, #,
    data = list(test = list(
                            X = X_sing,
                            mean = new_mu_sing),
                neighborhood = list(X = X_neighborhood,
                                    mean = new_mu_neighb)
    ),
    models = list(
      selection = lapply(singleModels, function(s) s$theta),
      projection = lapply(singleModelsP, function(s) s$theta),
      estimation = list(boot = post_sample$test$eta[c(single_idx,neighb_idx),],
                        original = list(prediction = pred.orig,
                                        derivatives = deriv.orig)
                        ),
    derivatives = derivative.test
  ))

  return(outList)
}
