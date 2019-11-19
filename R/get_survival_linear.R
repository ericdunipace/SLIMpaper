get_survival_linear_model <- function() {
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

  # trajectory <- function(nT, cov, p, ARM) {
  #   sds <- sqrt(diag(cov))
  #   corr <- diag(1/sds) %*% cov %*% diag(1/sds)
  #   x <- rmvnorm(nT, rep(0, p), corr)
  #   for(i in 2:nT) x[i,] <- x[i-1,] * ARM + x[i,]
  #
  #   return(x)
  # }

  # rX <- function(n, corr, p,...)
  # {
  #   dots <- list(...)
  #   nT <- dots$nT
  #   if(is.null(nT)) nT <- 1
  #   ARM <- dots$ARM
  #   if(is.null(ARM)) ARM <- rnorm(p-1)
  #   if (is.matrix(corr)) {
  #     x_cov <- corr
  #   }
  #   else {
  #     x_cov <- matrix(corr, nrow = p - 1, ncol = p - 1)
  #   }
  #   diag(x_cov) <- 1
  #   x_list <- lapply(1:n, function(i) trajectory(nT, x_cov, p-1, ARM))
  #   x <- matrix(NA, nrow=n*nT, ncol=p-1)
  #   for(i in 1:n) x[(nT*(i-1)+1):(nT*i),] <- x_list[[i]]
  #   X <- cbind(1,x)
  #   return(X)
  # }

  rX <- gen_x()$rX

  #### Y  Data ####
  rdata <- function(n, x, theta, ...) {
    dots <- list(...)
    # id <- dots$id

    log.rate <- x %*% c(theta)
    mu <- exp(log.rate)
    Y <- rexp(n, mu)
    return(list(Y= Y, mu = mu, eta = log.rate, link = Gamma(link = log)$linkfun, invlink = Gamma(link = log)$linkinv, param = theta))
  }

  #### Parameters ####
  rparam <- get_param()$exponential

  #### Posterior on Coefficients
  rpost_coef <- function(x){NULL}
  rpost_sigma <- function(x){NULL}

  rpost <- function(n.samp, x, y, hyperparameters,...) { #uses RJAGS
    # require(rjags)
    dots <- list(...)
    id <- dots$id
    follow.up <- y
    fail <- dots$fail
    if (is.matrix(follow.up)) {
      if(ncol(follow.up) == 2 & is.null(fail)) {
        fail <- follow.up[,2]
        follow.up <- follow.up[,1]
      } else if (ncol(follow.up) == 2 & !is.null(fail) ){
        stop("if fail is null, then follow.up must be a two column matrix with the first column the failure times and the second column the event indicator. Do not specify both fail and make follow.up a matrix")
      }
    }
    jags_dir <- dots$jags_dir
    thin <- dots$thin
    model <- dots$model
    nchain <- dots$nchain
    X.test <- dots$X.test
    method <- dots$method
    seed <- dots$seed
    parallel <- dots$parallel
    cutpoints <- dots$cutpoints
    n.intervals <- dots$n.intervals
    if(is.null(method)) method <- "ss-jags"

    test <- list(eta = NULL, mu = NULL)

    if(is.null(thin)) thin <- 1
    if(all(x[,1]==1) || all(x[,1]==0)) x <- x[,-1]
    if(is.null(nchain)) nchain <- 1
    if(is.null(id)) id <- 1:length(y)


    # a <- hyperparameters$alpha
    # b <- hyperparameters$beta
    m <- hyperparameters$mu
    s <- hyperparameters$sigma
    prior_frac <- hyperparameters$prior_frac
    spike_a <- hyperparameters$spike_a
    spike_b <- hyperparameters$spike_b

    if(is.null(spike_a)) spike_a <- 1
    if(is.null(spike_b)) spike_b <- 1

    lambda <- 1/s
    # n <- length(y)

    # get_x <- function(x, id, follow.up, time) {
    #   X <- matrix(NA, nrow=length(id) * length(follow.up), ncol=ncol(x))
    #   id_unique <- unique(id)
    #   for(i in seq_along(id_unique)){
    #     fup <- follow.up[i]
    #     current_id <- id == id_unique[i]
    #     curr_x <- x[current_id ,]
    #
    #   }
    # }

    # if(any(follow.up > dots$time)) follow.up[follow.up > dots$time] <- max(dots$time)
    times <- sort(c(0,unique(follow.up)))
    mu <- eta <- alpha <- theta <- NULL


    # obs.time <- tapply(dots$time, id, max)
    obs.time <- follow.up
    if (is.null(fail )) {
      fail <- rep(1, length(unique(id))) # ifelse(obs.time < 5, 1, 0)
    }

    surv.calc <- function(eta, fit) {
      n <- nrow(eta)
      max.t <- max(fit$Surv[, 1], na.rm = T)
      last.d <- nrow(fit$d.scaled)
      d.scaled <- fit$d.scaled
      d.scaled[last.d,] <- max.t
      length.int <- apply(d.scaled, 2, diff)
      H <- apply(length.int * fit$h.scaled[-1,], 2, cumsum)
      mu <- lapply(1:n, function(i) exp(-H * exp(eta[i,])))
      for(i in 1:n) rownames(mu[[i]]) <- paste0("(",format(d.scaled[1:(last.d-1),1]),", ", format(d.scaled[2:last.d,1]),"]")
      return(mu)
    }

    if(!is.null(parallel)) {
      if( parallel ) {
        parallel.MPI <- TRUE
        ncpu <- parallel::detectCores()-1
      } else {
        parallel.MPI <- FALSE
        ncpu <- 1
      }
    } else {
      parallel.MPI <- FALSE
      ncpu <- 1
    }

    if(method == "ss-jags") {
      mu_vec <- rep(mu_0, ncol(x))
      prec_matrix <- matrix(0, ncol(x), ncol(x))
      diag(prec_matrix) <- lambda
      jags_data <-
        list(
          N = length(unique(id)),
          T = length(times)-1,
          P = ncol(x),
          obs.time = obs.time,
          time = times,
          eps = 0.0001,
          X = x,
          fail = fail,
          mu_0 = as.double(m),
          prec_0 = as.double(lambda),
          prior_frac = prior_frac,
          spike_a = spike_a,
          spike_b = spike_b,
          mu_vec = mu_vec,
          prec_matrix = prec_matrix
        )
      if (is.null(model)) {
        model <- rjags::jags.model(file=jags_dir,data = jags_data, n.chains = nchain, n.adapt = n.samp*thin*.1)
      } else {
        model$recompile()
      }
      total_iter <- n.samp*thin + ceiling(n.samp*thin/nchain)
      samples <- rjags::jags.samples(model,
                                     variable.names = c("beta","alpha","prob","baseHaz"),
                                     n.iter = total_iter, thin = thin)
      remove_burnin.idx <- round(seq(nsamp+1, total_iter, length.out = nsamp))
      final_nsamp <- length(remove_burnin.idx)
      # nsamp_portion <- floor(final_nsamp/nchain)
      # adjust_remove.idx <- remove_burnin.idx[((nsamp - nsamp_portion + 1):nsamp)]
      adjust_remove.idx <- remove_burnin.idx[((nsamp - final_nsamp + 1):nsamp)]

      theta_samp <- samples$beta[ , adjust_remove.idx, ]
      alpha_samp <- samples$alpha[, adjust_remove.idx, ]

      if(nchain > 1) {
        theta <- matrix(NA, ncol=final_nsamp, nrow = dim(theta_samp)[1])
        alpha <- matrix(NA, ncol=final_nsamp, nrow = dim(alpha_samp)[1])
        for(i in 1:nchain){
          get_rows <- (i-1)*nsamp_portion + (1:nsamp_portion)
          theta[,get_rows] <- theta_samp[,,i]
          alpha[,get_rows] <- alpha_samp[,,i]
        }
      } else {
        theta <- theta_samp
        alpha <- alpha_samp
      }

      eta <- cbind(1,x) %*% theta
      if(!is.null(X.test)){
        if(all(X.test[,1]==1) || all(X.test[,1]==0)) X.test <- cbind(1,X.test)
        test$eta <- X.test %*% theta
      }
    } else if(method == "jags") {
      mu_vec <- rep(mu_0, ncol(x))
      prec_matrix <- matrix(0, ncol(x), ncol(x))
      diag(prec_matrix) <- lambda
      jags_data <-
        list(
          N = length(unique(id)),
          T = length(times)-1,
          P = ncol(x),
          obs.time = obs.time,
          time = times,
          eps = 0.0001,
          X = x,
          fail = fail,
          mu_0 = as.double(m),
          prec_0 = as.double(lambda),
          prior_frac = prior_frac,
          spike_a = spike_a,
          spike_b = spike_b,
          mu_vec = mu_vec,
          prec_matrix = prec_matrix
        )
      if (is.null(model)) {
        model <- rjags::jags.model(file=jags_dir,data = jags_data, n.chains = nchain, n.adapt = n.samp*thin*.1)
      } else {
        model$recompile()
      }
      total_iter <- n.samp*thin + ceiling(n.samp*thin/nchain)
      samples <- rjags::jags.samples(model,
                                     variable.names = c("beta", "baseHaz"),
                                     n.iter = total_iter, thin = thin)
      saved_iter <- n.samp + ceiling(n.samp / nchain)
      nsamp_portion <- ceiling(n.samp/nchain)
      remove_burnin.idx <- round(seq(nsamp+1, saved_iter, length.out = nsamp_portion))
      # final_nsamp <- length(remove_burnin.idx)
      # adjust_remove.idx <- remove_burnin.idx[((nsamp - nsamp_portion + 1):nsamp)]

      theta_samp <- samples$beta[ , remove_burnin.idx, ]
      # alpha_samp <- samples$alpha[, adjust_remove.idx, ]

      if(nchain > 1) {
        theta <- matrix(NA, ncol=nsamp, nrow = dim(theta_samp)[1])
        # alpha <- matrix(NA, nrow=final_nsamp, ncol = dim(alpha_samp)[1])
        for(i in 1:nchain){
          get_rows <- (i-1)*nsamp_portion + (1:nsamp_portion)
          theta[,get_rows] <- theta_samp[,,i]
          # alpha[get_rows,] <- t(alpha_samp[,,i])
        }
      } else {
        theta <- theta_samp
        # alpha <- t(alpha_samp)
      }
      if(nrow(theta) == (ncol(x) + 1)){
        X.prod <- cbind(1,x)
      } else{
        X.prod <- x
      }
      eta <- X.prod %*% theta
      if(!is.null(X.test)){
        if(ncol(theta) == (ncol(X.test) + 1)) X.test <- cbind(1,X.test)
        test$eta <- X.test %*% theta
      }
    } else if (method == "bvs-cox") {
      # require(BVSNLP)
      # select median prob model first so don't have to run all 13000+ covariates in bayesian regression
      cat("Selecting median probability model (MPM)")
      resp <- cbind(follow.up, fail)
      xdf <- as.data.frame(log(x))
      sel <- BVSNLP::bvs(X = xdf, resp = resp, family = "survival",
                         niter = n.samp*10, prep = TRUE, mod_prior = "unif", logT = FALSE,
                         inseed = seed, ncpu = ncpu, parallel.MPI = parallel.MPI)
      mpm <- sel$MPM
      if(length(mpm)==0) {
        mpm <- sort(unique(unlist(sel$max_models[sel$max_prob_vec > (log(0.01) + sel$max_prob)])))
      }
      # eta <- BVSNLP::predBMA(fit, X=df, resp = resp, prep=TRUE, logT=FALSE, family = "survival")
      # eta <- fit$des_mat %*% fit$beta_hat
      # alpha <- fit$inc_probs
      #
      # test$eta <- X.test[,colnames(fit$des_mat)] %*% fit$beta_hat
      survform <- formula(survival::Surv(time = follow.up, event = fail) ~ .)
      # x_sc <- scale(log(x))
      # xt_sc <- NULL
      # if(!is.null(X.test)) {
      #   xt_sc <- scale(log(X.test), center = attr(x_sc,"scaled:center"), scale = attr(x_sc, "scaled:scale"))
      # }
      x_sc <- xdf[,mpm, drop=FALSE]
      df <- as.data.frame(cbind(follow.up, fail, x_sc))
      colnames(df) <- c("follow.up","fail", colnames(x_sc))
      pred <- list(xpred = x_sc)
      if(!is.null(X.test)) {
        xt_sc <- scale(log(X.test[,mpm,drop=FALSE]), center = colMeans(x_sc), scale = colSD(x_sc))
        pred <- list(xpred = rbind(scale(x_sc), xt_sc))
      }
      cat("Running Bayesian cox on MPM")
      fit <- spBayesSurv::indeptCoxph(survform, data = df, prediction = pred,
                                      mcmc = list(nburn = n.samp, nsave = n.samp, nskip = 0, ndisplay = 100),
                                      prior = NULL, state = NULL, scale.designX = TRUE)

      eta <- fit$X.scaled %*% fit$beta.scaled
      mu <- list(time = NULL, S = NULL)
      mu$time <- fit$Tpred[1:n,,drop=FALSE]
      mu$S <- surv.calc(eta, fit)
      if(!is.null(X.test)) {
        test$eta <- xt_sc %*% fit$beta.scaled
        test$mu <- list(time = NULL, S = NULL)
        test$mu$time <- fit$Tpred[(n+1):nrow(fit$Tpred),,drop=FALSE]
        test$mu$S <- surv.calc(test$eta, fit)
      }
      model <- fit
    } else if (method == "cox") {

      x_sc <- scale(x) #log(x)
      df <- as.data.frame(cbind(follow.up, fail, x_sc))
      colnames(df) <- c("follow.up","fail", colnames(x_sc))
      pred <- list(xpred = x_sc)
      if(!is.null(X.test)) {
        xt_sc <- scale(X.test[,mpm,drop=FALSE], center = colMeans(x_sc), scale = colSD(x_sc))
        pred <- list(xpred = rbind(scale(x_sc), xt_sc))
      }
      survform <- formula(survival::Surv(time = follow.up, event = fail) ~ .)
      fit <- spBayesSurv::indeptCoxph(survform, data = df, prediction = pred,
                                      mcmc = list(nburn = n.samp*thin, nsave = n.samp, nskip = thin,
                                                  ndisplay = n.samp/10),
                                      prior = NULL, state = NULL, scale.designX = TRUE)

      eta <- fit$X.scaled %*% fit$beta.scaled
      theta <- fit$beta.scaled
      mu <- list(time = NULL, S = NULL)
      mu$time <- fit$Tpred[1:n,,drop=FALSE]
      mu$S <- surv.calc(eta, fit)
      if(!is.null(X.test)) {
        test$eta <- xt_sc %*% fit$beta.scaled
        test$mu <- list(time = NULL, S = NULL)
        test$mu$time <- fit$Tpred[(n+1):nrow(fit$Tpred),,drop=FALSE]
        test$mu$S <- surv.calc(test$eta, fit)
      }
      model <- fit
    } else if (method == "bvs") {
      theta <- NULL
      alpha <- NULL
      mu <- NULL
      eta <- NULL
      # require(BVSNLP)
      # require(doParallel)
      # require(foreach)

      resp <- cbind(follow.up, fail)
      xdf <- as.data.frame(log(x))
      model <- bvsmod(X = xdf, resp = resp, family = "survival",
                         niter = n.samp, prep = TRUE, mod_prior = "unif", logT = FALSE,
                         inseed = seed, ncpu = ncpu, parallel.MPI = parallel.MPI)
    } else if (method == "inla") {
      # if(all(x[,1]==1) || all(x[,1]==0)) x <- x[,-1]
      sx <- scale(x)
      xdf <- as.data.frame(sx)
      pred.names <- colnames(xdf)
      colnames(xdf) <- paste0("x",1:ncol(xdf))
      intercept1 <- rep(1,n)
      # index_df <- matrix(1,nrow=n, ncol=ncol(xdf))
      # index_names <- paste0("index_", colnames(xdf))
      # colnames(index_df) <- index_names
      # df <- cbind(data.frame(time = follow.up, event = fail), xdf, index_df)
      #lambda = 0.4712777
      # hc <- "expression:
      #         lambda = 0.01;
      #         precision = exp(log_precision);
      #         logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
      #         log_jacobian = log_precision;
      #         return(logdens+log_jacobian);"
      # hcprior <- list(prec = list(prior = hc))

      # fs <- paste0("f(", index_names,",", colnames(xdf), ", model = 'iid', hyper = hcprior)")
      # survform <- formula(paste(c("inla.surv(time, event) ~ 1 ", fs), collapse = " + "))
      df <- cbind(data.frame(time = follow.up, event = fail), intercept1, xdf)
      survform <- formula(paste(c("inla.surv(time, event) ~ 0 + intercept1", colnames(xdf)), collapse = " + "))
      # if(is.null(cutpoints)) cutpoints <- NULL
      if(is.null(n.intervals)) n.intervals <- 15
      cox.call <- inla.coxph(survform, df, control.hazard = list(model = "rw1",
                                                                 n.intervals = n.intervals,
                                                                 cutpoints = cutpoints,
                                                                 constr = TRUE))
      # timings <- proc.time()
      # model <- inla(survform, family = "coxph",
      #               data = df, control.fixed=list(mean=m[1], prec=1/(s[1])),
      #               debug = TRUE, verbose=FALSE)
      model.gauss <- inla(cox.call$formula, family = cox.call$family,
                    data = c(as.list(cox.call$data), cox.call$data.list),
                    E = cox.call$E,
                    control.fixed=list(mean=m[1], prec=lambda[1]),
                    control.inla = list(strategy = "gaussian")
                    )
      model <- inla(cox.call$formula, family = cox.call$family,
                    data = c(as.list(cox.call$data), cox.call$data.list),
                    E = cox.call$E,
                    control.fixed=list(mean=m[1], prec=lambda[1]),
                    control.inla = list(strategy = "laplace", fast = FALSE),
                    control.mode = list(result = model.gauss, restart = TRUE),
                    control.compute=list(config=TRUE))
      # print(proc.time() - timings)
      samples <- inla.posterior.sample(n = n.samp, result = model, intern = FALSE,
                                     use.improved.mean = TRUE,
                                     add.names = FALSE, seed = -1L, num.threads = 1L)
      fe.names <- model$names.fixed
      st.names <- rownames(samples[[1]]$latent)
      pred.exclude <- !grepl("Predictor", st.names)
      bh.exclude <- !grepl("baseline.hazard", st.names)
      which.fe <- which(pred.exclude & bh.exclude)
      which.not.pred <- which(pred.exclude)
      theta <- matrix(sapply(samples, function(ss) ss$latent[which.fe,]), ncol=n.samp)
      rownames(theta) <- fe.names
      save.samples <- sapply(samples, function(ss) ss$latent[which.not.pred,])

      model$samples <- save.samples
      surv.calc <- function(model, theta, x) {
        n <- nrow(x)
        intercept <- theta[1,]
        theta_reg <- theta[-1,]
        eta <- x %*% theta_reg
        haz.times <- model$.args$data$baseline.hazard.values
        # if(is.null(times)) times <- haz.times

        log_BH <- model$samples[grep("baseline.hazard", rownames(model$samples)),]
        nT <- nrow(log_BH)
        nS <- ncol(theta)
        log_BH <- log_BH + matrix(intercept, nT, nS, byrow=TRUE)
        # cutTimes <- cut(times, haz.times, include.lowest = TRUE)
        BH <- exp(log_BH) * diff(c(haz.times, Inf))
        # BH_times <- BH
        cumHaz <- apply(BH,2,cumsum)
        baseSurv <- exp(-cumHaz)
        rownames(baseSurv) <- haz.times

        Surv <- simplify2array(lapply(1:n, function(i) baseSurv^matrix(exp(eta[i,]), nT, nS, byrow=TRUE)))

        return(list(surv = Surv, base = baseSurv))
      }

      survlist <- surv.calc(model, theta, sx)
      mu$S <- survlist
      # mu$time
      eta <- sx %*% theta[-1,]

      if(!is.null(X.test)) {
        mm <- attr(sx, "scaled:center")
        ss <- attr(sx, "scaled:scale")
        sxt <- scale(X.test, center =mm, scale = ss)
        test$eta <- sxt %*% theta[-1,]
        test$mu <- list()
        test$mu$S <- surv.calc(model, theta, sxt)
      }



    } else if (method == "inla.cens") {
      # if(all(x[,1]==1) || all(x[,1]==0)) x <- x[,-1]
      # xdf <- as.data.frame(scale(x))
      # pred.names <- colnames(xdf)
      # colnames(xdf) <- paste0("x",1:ncol(xdf))
      intercept1 <- rep(1,n)
      # index_df <- matrix(1,nrow=n, ncol=ncol(xdf))
      # index_names <- paste0("index_", colnames(xdf))
      # colnames(index_df) <- index_names
      # df <- cbind(data.frame(time = follow.up, event = fail), xdf, index_df)
      #lambda = 0.4712777
      # hc <- "expression:
      #         lambda = 0.01;
      #         precision = exp(log_precision);
      #         logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
      #         log_jacobian = log_precision;
      #         return(logdens+log_jacobian);"
      # hcprior <- list(prec = list(prior = hc))

      # fs <- paste0("f(", index_names,",", colnames(xdf), ", model = 'iid', hyper = hcprior)")
      # survform <- formula(paste(c("inla.surv(time, event) ~ 1 ", fs), collapse = " + "))
      df <- cbind(data.frame(time = follow.up, event = fail), intercept1)
      survform <- formula(inla.surv(time, event) ~ 0 + intercept1)
      # if(is.null(cutpoints)) cutpoints <- NULL
      if(is.null(n.intervals)) n.intervals <- 15
      cox.call <- inla.coxph(survform, df, control.hazard = list(model = "rw1",
                                                                 n.intervals = n.intervals,
                                                                 cutpoints = cutpoints,
                                                                 constr = TRUE))
      # timings <- proc.time()
      # model <- inla(survform, family = "coxph",
      #               data = df, control.fixed=list(mean=m[1], prec=1/(s[1])),
      #               debug = TRUE, verbose=FALSE)
      model.gauss <- inla(cox.call$formula, family = cox.call$family,
                          data = c(as.list(cox.call$data), cox.call$data.list),
                          E = cox.call$E,
                          control.fixed=list(mean=m[1], prec=lambda[1]),
                          control.inla = list(strategy = "gaussian")
      )
      model <- inla(cox.call$formula, family = cox.call$family,
                    data = c(as.list(cox.call$data), cox.call$data.list),
                    E = cox.call$E,
                    control.fixed=list(mean=m[1], prec=lambda[1]),
                    control.inla = list(strategy = "laplace", fast = FALSE),
                    control.mode = list(result = model.gauss, restart = TRUE),
                    control.compute=list(config=TRUE))
      # print(proc.time() - timings)
      samples <- inla.posterior.sample(n = n.samp, result = model, intern = FALSE,
                                       use.improved.mean = TRUE,
                                       add.names = FALSE, seed = -1L, num.threads = 1L)
      fe.names <- model$names.fixed
      st.names <- rownames(samples[[1]]$latent)
      pred.exclude <- !grepl("Predictor", st.names)
      bh.exclude <- !grepl("baseline.hazard", st.names)
      which.fe <- which(pred.exclude & bh.exclude)
      which.not.pred <- which(pred.exclude)
      theta <- matrix(sapply(samples, function(ss) ss$latent[which.fe,]), ncol=n.samp)
      rownames(theta) <- fe.names
      save.samples <- sapply(samples, function(ss) ss$latent[which.not.pred,])

      model$samples <- save.samples
      surv.calc <- function(model, theta) {
        n <- nrow(x)
        intercept <- theta[1,]
        # theta_reg <- theta[-1,]
        # eta <- x %*% theta_reg
        haz.times <- model$.args$data$baseline.hazard.values
        # if(is.null(times)) times <- haz.times

        log_BH <- model$samples[grep("baseline.hazard", rownames(model$samples)),]
        nT <- nrow(log_BH)
        nS <- ncol(theta)
        log_BH <- log_BH + matrix(intercept, nT, nS, byrow=TRUE)
        # cutTimes <- cut(times, haz.times, include.lowest = TRUE)
        BH <- exp(log_BH) * diff(c(haz.times, Inf))
        # BH_times <- BH
        cumHaz <- apply(BH,2,cumsum)
        baseSurv <- exp(-cumHaz)
        rownames(baseSurv) <- haz.times

        # Surv <- simplify2array(lapply(1:n, function(i) baseSurv^matrix(exp(eta[i,]), nT, nS, byrow=TRUE)))
        Surv <- baseSurv
        return(list(surv = Surv, base = baseSurv))
      }

      survlist <- surv.calc(model, theta)
      mu$S <- survlist
      # mu$time

      # eta <- x %*% theta[-1,]
    }


    return(list(theta=theta, alpha = alpha, mu = mu, eta = eta, test = test, model = model, surv.calc = surv.calc))

  }

  cindex <- function(times, event, risk = NULL, surv = NULL, surv.times=NULL, cens = NULL,...) {


    # ntimes <- length(surv.times)
    # if (is.null(risk) & is.null(surv) ) {
    #   stop("Must provide either survival data or risk data")
    # } else if (is.null(surv)) {
    #
    #   pred <- array(risk, dim = c(ntimes, dim(risk)[2:1]))
    #
    # } else if (is.null(risk)) {
    #   pred <- 1 - surv
    # } else {
    #   pred <- array(risk, dim = c(ntimes, dim(risk)))
    # }
    #
    # if (is.null(surv.times)) {
    #   surv.times <- as.numeric(sapply(strsplit(dimnames(surv)[[1]],":"), function(x) x[2]))
    #   ntimes <- length(surv.times) - 1
    #   surv.times <- surv.times[-1]
    # }
    # cut.times <- cut(times, c(0,surv.times), include.lowest = TRUE)
    #
    # n <- dim(pred)[3]
    if (is.null(risk) & is.null(surv) ) {
      stop("Must provide either survival data or risk data")
    } else if (is.null(surv)) {

      pred <- -risk

    } else if (is.null(risk)) {
      pred <- surv
    } else {
      pred <- -risk
    }
    nsamp <- dim(pred)[2]
    #
    # stopifnot(n == length(event))
    #
    # # sort times if using rank method
    # # ordtime <- order(times)
    # # event <- event[ordtime]
    # # times <- times[ordtime]
    #
    # # indicator of event times
    # eventbycol <- matrix(event, ncol=n, nrow=n)
    # eventbyrow <- matrix(event, ncol=n,nrow=n, byrow = TRUE)
    # timeless <- sapply(times, function(tt) as.integer(tt < times))
    # timeequal <- sapply(times, function(tt) as.integer(tt == times))
    # K <- (timeless  + timeequal * eventbycol) * eventbyrow
    # # diag(K) <- 0 #don't care for same person
    #
    # # ranks
    # # timgings <- proc.time()
    #   rank_fun <- function(i, n, surv) {
    #     ranks <- apply(surv[,i:n, drop=FALSE], 1, rank)
    #     avg_rank <- colMeans(ranks)
    #     return(2*(ranks[1,] - avg_rank))
    #   }
    #   # U <- D <- matrix(NA, nrow= ntimes , ncol = nsamp)
    #   # for(tt in seq_along(surv.times)) {
    #   #   timetrue <- times <= surv.times[tt]
    #   #   idx.dead <- which(timetrue & event == 1)
    #   #   idx.time <- which(timetrue)
    #   #
    #   #   U[tt, ] <- rowSums(sapply(idx.dead, rank_fun, n = n, surv = pred[tt,,]))
    #   #   D[tt, ] <- Reduce("+", sapply(idx.time, function(i) if(i < n) {sum(K[(i+1):n,i])} else {0}))
    #   # }
    #   # cstat <- 0.5*(U/D + 1)
    #   # rowMeans(cstat)
    #   # print(proc.time() - timgings)
    # #risk of predictors over time and over posterior samples
    # # risk_array <- sapply(1:n, function(i) K[i,-i] * (as.numeric(pred[,,-i] < pred[,,i]) +
    #                            # as.numeric(pred[,,-i] == pred[,,i])/2))
    #   # timgings <- proc.time()
    # risk_mat <- matrix(NA, nrow=ntimes, ncol=nsamp)
    # compare <- lapply(1:n, function(i) matrix(K[-i,i], nrow=nsamp, ncol=n-1, byrow=TRUE) *
    #                     (pred[cut.times[i],,i] > pred[cut.times[i],,-i] + (pred[cut.times[i],,i] == pred[cut.times[i],,-i])/2))
    # oldidx <- idx <- newidx <- prev <- NULL
    #
    # for (tt in seq_along(surv.times)) {
    #   idx <- which(times <= surv.times[tt] & event == 1)
    #   newidx <- idx[!(idx %in% oldidx)]
    #   # not.idx <- which(times > surv.times[tt])
    #   if(tt>1) {
    #     prev <- risk_mat[tt-1,]
    #   } else {
    #     prev <- 0
    #   }
    #   risk_mat[tt,] <- prev + rowSums(Reduce("+", compare[newidx]))
    #   oldidx <- idx
    # }
    #
    # denoms <- sapply(surv.times, function(ss) sum(K[,times <= ss]))
    # cstat <- risk_mat/matrix(denoms, nrow=ntimes, ncol=nsamp)
    survobj <- survival::Surv(time = times, event = event)
    cstat <- sapply(1:nsamp, function(i) survival::concordancefit(y = survobj, x = pred[,i],
                                                                  timefix = TRUE)$concordance)
    # print(proc.time() - timgings)
    output <- list(mean = mean(cstat), low = quantile(cstat, 0.025),
                   high = quantile(cstat, 0.975), cindex = cstat)
    return(output)
  }

  brier.score <- function(times, event, surv, surv.times, cens_prob, ...) {
    # stopifnot(all(dim(event) == dim(probs) ))
    readTime <- sort(unique(times))
    nst <- length(surv.times)
    ntimes <- length(readTime)
    # nct <- length(cens.times)

    if (nst != ntimes) {
      if(dim(surv)[1] != nst) stop("surv.times must match first dimension of surv probs")
        # surv <- expandPred(readTime, surv, surv.times)
    }
    if(dim(cens_prob)[1] != ntimes) {
      if(dim(cens_prob)[1] != nst) stop("first dimension of censoring probabilities must match length of surv.times or number of unique times in data")
      # cens_prob <- expandPred(readTime, cens_prob, surv.times)
      cens.times <- surv.times
    } else {
      cens.times <- sort(unique(c(0,readTime)))
    }
    # if (nct != ntimes) {
    #   cens_prob <- expandPred(readTime, cens, surv.times)
    # }
    # ot <- order(times)
    # times <- times[ot]
    # event <- event[ot]
    n <- dim(surv)[3]
    ncens <- dim(cens_prob)[3]
    nsamp <- dim(surv)[2]

    eventHappen <- t(sapply(readTime, function(curTime)
      as.integer( (times <= curTime) & (event == 1) )))
    eventNotHappen <- t(sapply(readTime, function(curTime)
      as.integer( times > curTime )))
    # cens_mat <- vector("list", ncens)
    bs <- matrix(0, nrow=ntimes, ncol=nsamp)
    time.idx <- as.numeric(cut(times, readTime, include.lowest = TRUE))
    cens_prob_at_event <- matrix(NA, nrow = nsamp, ncol=n)
    if(any(is.na(cens_prob))) {
      repl_idx <- which(is.na(cens_prob),arr.ind=TRUE)
      cens_prob[repl_idx] <- min(cens_prob, na.rm=TRUE)
    }
    if(any(cens_prob == 0)) {
      repl_idx <- which(cens_prob == 0, arr.ind=TRUE)
      cens_prob[repl_idx] <- Inf
    }

    for(i in 1:n) cens_prob_at_event[,i] <- cens_prob[time.idx[i],,i]

    # for ( curTime in readTime ) {
    #   idx_time <- which(readTime == curTime)
    #  for(s in 1:nsamp){
    #    for(i in 1:n) {
    #      bs[idx_time, s] <- bs[idx_time, s] + (0 - surv[idx_time,s,i])^2 *eventHappen[idx_time, i]/cens_prob[time.idx[i],s,i]/n +
    #        (1 - surv[idx_time,s,i])^2 * eventNotHappen[idx_time,i]/cens_prob[idx_time,s,i]/n
    #    }
    #  }
    #   # bs[idx_time, ] <- bs[idx_time, ]
    # }
    # lcpe <- log(cens_prob_at_event)
    # lcp <- log(cens_prob)
    # l_eventHappen <- log(eventHappen)
    # l_eventNotHappen <- log(eventNotHappen)
    l_normalize <- log_uwt1 <- log_uwt2 <- eh <- enh <- lterm1 <- lterm2 <- lsum1 <- lsum2 <-
      surv_at_time <- cnes_at_time <- normalize <- NULL
    # ttt <- proc.time()
    for ( curTime in readTime ) {
      idx_time <- which(readTime == curTime)
      eh <- matrix(eventHappen[idx_time,], nrow = nsamp, ncol = n, byrow = TRUE)
      enh <- matrix(eventNotHappen[idx_time,], nrow = nsamp, ncol = n, byrow=TRUE)
      # l_eh <- matrix(l_eventHappen[idx_time,], nrow = nsamp, ncol = n, byrow = TRUE)
      # l_enh <- matrix(l_eventNotHappen[idx_time,], nrow = nsamp, ncol = n, byrow = TRUE)
      # log_uwt1 <- apply(l_eh - lcpe, 1, log_sum_exp)
      # log_uwt2 <- apply(l_enh - lcp[idx_time,,], 1, log_sum_exp)
      # l_normalize <- log_sum_exp2(log_uwt1, log_uwt2)
      # lterm1 <- 2 * log(abs(0 - surv[idx_time,,])) + l_eh - lcpe - l_normalize
      # lterm2 <- 2 * log(1 - surv[idx_time,,]) + l_enh - lcp[idx_time,,] - l_normalize
      # lsum1 <- apply(lterm1, 1, log_sum_exp)
      # lsum2 <- apply(lterm2, 1, log_sum_exp)
      # bs[idx_time, ] <- exp(log_sum_exp2(lsum1, lsum2))
      surv_at_time <- predAtTime(curTime, surv, surv.times)
      cens_at_time <- predAtTime(curTime, cens_prob, cens.times)
      normalize <- rowSums(eh / cens_prob_at_event + enh / cens_at_time)
      bs[idx_time, ] <- rowSums((0 -  surv_at_time)^2 * eh / cens_prob_at_event / normalize +
      ( 1 - surv_at_time  )^2 * enh / cens_at_time/ normalize )
      # for (samp in nsamp) {
      #   for(cc in 1:ncens) {
      #     cens_mat[[cc]] <-  ((0 - surv[,samp,] )^2 /cens_prob[,cc,] * eventHappen +
      #     ( 1 - surv[,samp,] )^2 /cens_prob[idx_time,cc,] * eventNotHappen)/ncens
      #   }
      #   int_cens <- Reduce("+", cens_mat)
      #   bs[idx_time, samp] <- mean( int_cens )
      # }
    }
    # print(proc.time() - ttt)
    # idx <- 2:ntimes
    # intbs <- diff(readTime) %*% ( ( bs[idx -1,] + bs[idx,]) / 2 )
    # intbs <- intbs/diff(range(readTime))
    intbs <- (diff(c(0,readTime)) %*% bs)/max(readTime)


    output <- list(
        brier.score = list(
          mean = rowMeans(bs),
          low = apply(bs, 1, quantile, 0.025),
          high = apply(bs,1, quantile, 0.975),
          median = apply(bs,1, median),
          bscore = bs
          ),
        int.BS = list(
         mean = mean(intbs),
         low = quantile(intbs, 0.025),
         high = quantile(intbs, 0.975),
         median = median(intbs),
                      intBS = c(intbs))
     )

    return(output)

  }

  evalfit <- function(times, event, risk = NULL, surv = NULL, cens=NULL, surv.times= NULL, method = c("c-index","brier")) {
    meth <- match.arg(method)
    efun <- switch(meth, "c-index" = cindex,
           "brier" = brier.score)
    stat <- efun(times = times, event = event, risk=risk, surv=surv, surv.times=surv.times, cens = cens)

    return(stat)
  }

  predAtTime <- function(time, pred, predTimes) {
    dict <- cut(time, predTimes, include.lowest = TRUE)
    idx <- as.numeric(dict)
    predTime <- pred[idx,,]
    return(predTime)
  }

  expandPred <- function(times, pred, predTimes) {
    dims <- if (length(dim(pred)) > 1) {
      dim(pred)
    } else {
      length(pred)
    }
    stopifnot(dims[1] == length(predTimes))
    st <- sort(unique(times))
    dict <- cut(st, predTimes, include.lowest = TRUE)
    idx <- as.numeric(dict)
    predExp <- pred[idx,,]
    return(predExp)
  }

  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = NULL,
              rparam = rparam,
              link = Gamma(link = log)$linkfun,
              invlink = Gamma(link = log)$linkinv,
              evalfit = evalfit,
              expandPred = expandPred))
}
