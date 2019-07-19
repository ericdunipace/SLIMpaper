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
    if(is.null(method)) method <- "ss-jags"

    test <- list(eta = NULL, mu = NULL)

    if(is.null(thin)) thin <- 1
    if(all(x[,1]==1)) x <- x[,-1]
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
          spike_b = spike_b
        )
      if (is.null(model)) {
        model <- rjags::jags.model(file=jags_dir,data = jags_data, n.chains = nchain, n.adapt = n.samp*thin*.1)
      } else {
        model$recompile()
      }
      samples <- rjags::jags.samples(model,
                                     variable.names = c("beta","alpha","prob"),
                                     n.iter = n.samp*2*thin, thin = thin)
      remove_burnin.idx <- round(seq(nsamp+1, nsamp*2, length.out = nsamp))
      final_nsamp <- length(remove_burnin.idx)
      nsamp_portion <- floor(final_nsamp/nchain)
      adjust_remove.idx <- remove_burnin.idx[((nsamp - nsamp_portion + 1):nsamp)]


      theta_samp <- samples$beta[ , adjust_remove.idx, ]
      alpha_samp <- samples$alpha[, adjust_remove.idx, ]

      if(nchain > 1) {
        theta <- matrix(NA, nrow=final_nsamp, ncol = dim(theta_samp)[1])
        alpha <- matrix(NA, nrow=final_nsamp, ncol = dim(alpha_samp)[1])
        for(i in 1:nchain){
          get_rows <- (i-1)*nsamp_portion + (1:nsamp_portion)
          theta[get_rows,] <- t(theta_samp[,,i])
          alpha[get_rows,] <- t(alpha_samp[,,i])
        }
      } else {
        theta <- t(theta_samp)
        alpha <- t(alpha_samp)
      }

      eta <- tcrossprod(cbind(1,x), theta)
      if(!is.null(X.test)){
        if(!(all(X.test[,1]==1))) X.test <- cbind(1,X.test)
        test$eta <- tcrossprod(X.test, theta)
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

      x_sc <- log(x)
      df <- as.data.frame(cbind(follow.up, fail, x_sc))
      colnames(df) <- c("follow.up","fail", colnames(x_sc))
      pred <- list(xpred = x_sc)
      if(!is.null(X.test)) {
        xt_sc <- scale(log(X.test[,mpm,drop=FALSE]), center = colMeans(x_sc), scale = colSD(x_sc))
        pred <- list(xpred = rbind(scale(x_sc), xt_sc))
      }
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
    } else if (method == "bvs") {
      theta <- NULL
      alpha <- NULL
      mu <- NULL
      eta <- NULL

      resp <- cbind(follow.up, fail)
      xdf <- as.data.frame(log(x))
      model <- BVSNLP::bvs(X = xdf, resp = resp, family = "survival",
                         niter = n.samp*10, prep = TRUE, mod_prior = "unif", logT = FALSE,
                         inseed = seed, ncpu = ncpu, parallel.MPI = parallel.MPI)
    }


    return(list(theta=theta, alpha = alpha, mu = mu, eta = eta, test = test, model = model, surv.calc = surv.calc))


  }

  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = NULL,
              rparam = rparam,
              link = Gamma(link = log)$linkfun,
              invlink = Gamma(link = log)$linkinv))
}
