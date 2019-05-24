
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
  
  rX <- get_x()$rX
  
  #### Y  Data ####
  rdata <- function(n, x, theta, ...) {
    dots <- list(...)
    id <- dots$id
    
    log.rate <- x %*% c(theta)
    return(rexp(n, exp(log.rate)))
  }
  
  #### Parameters ####
  rparam <- get_param()$exponential
  
  #### Posterior on Coefficients
  rpost_coef <- function(x){NULL}
  rpost_sigma <- function(x){NULL}
  
  rpost <- function(n.samp, x, y, hyperparameters,...) { #uses RJAGS
    require(rjags)
    dots <- list(...)
    id <- dots$id
    follow.up <- y
    jags_dir <- dots$jags_dir
    thin <- dots$thin
    model <- dots$model
    
    if(is.null(thin)) thin <- 1
    if(all(x[,1]==1)) x <- x[,-1]
    
    # a <- hyperparameters$alpha
    # b <- hyperparameters$beta
    m <- hyperparameters$mu
    s <- hyperparameters$sigma
    prior_frac <- hyperparameters$prior_frac
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
    
    
    # obs.time <- tapply(dots$time, id, max)
    obs.time <- follow.up
    fail <- rep(1, length(unique(id))) # ifelse(obs.time < 5, 1, 0)
    
    jags_data <- list(
      N = length(unique(id)),
      T = length(times)-1,
      P = ncol(x),
      obs.time = obs.time,
      time = times,
      eps = 0.0001,
      X = x,
      fail = fail,
      mu_0 = m,
      prec_0 = lambda,
      prior_frac = prior_frac
    )
    if (is.null(model)) {
      model <- jags.model(file=jags_dir,data = jags_data, n.chains = 4, n.adapt = n.samp*.1)
    } else {
      model$recompile()
    }
    samples <- jags.samples(model, 
                            variable.names = c("beta","alpha"), n.iter = n.samp*2, thin = thin)
    remove_burnin.idx <- round(n.samp/thin+1): round(n.samp/thin*2)
    final_nsamp <- length(remove_burnin.idx)
    theta_samp <- samples$beta[ , remove_burnin.idx, ]
    alpha_samp <- samples$alpha[, remove_burnin.idx, ]
    
    theta <- alpha <- matrix(NA, nrow=final_nsamp, ncol = dim(theta_samp)[1])
    for(i in 1:dim(theta_samp)[3]){
      get_rows <- (i-1)*final_nsamp + (1:final_nsamp)
      theta[get_rows,] <- theta_samp[,,i]
      alpha[get_rows,] <- alph_samp[,,i]
    }
    
    return(list(theta=theta, alpha = alpha, model=model))
    
    
  }
  
  return(list(rprior=rprior,
              rdata = rdata,
              rpost = rpost,
              X = list(rX = rX, corr = NULL),
              data_gen_function = NULL,
              rparam = rparam))
}
