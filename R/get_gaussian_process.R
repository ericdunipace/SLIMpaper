get_gp <- function() {
  kernelFn <- function(sigma, L, x,n) {
    xmat <- matrix(x, ncol=n,nrow=n)
    xmatT <- matrix(x, ncol=n, nrow=n, byrow = TRUE)
    Sigma <- sigma^2 * exp(-0.5/L^2 * (xmat - xmatT)^2)
    return(Sigma)
  }

  GP_prior <- function(sigma_f, L,m,x, iter){
    require(mvtnorm)
    n <- length(x)

    Sigma <- kernelFn(sigma_f, L, x, n)
    f <- t(rmvnorm(iter, mean=rep(m,n), sigma=Sigma, method="svd"))
    return(f)
  }

  post_f <- function(x, y, L, sigma_f, sigma_y, grid, iter) {
    # require(mvtnorm)
    #make united X vector
    X <- c(grid,x)
    nX <- length(X)

    #no which ones are observed data
    obs <- seq(length(grid)+1,nX,1)
    not_obs <- 1:length(grid)

    #store length of observed data
    n <- length(y)

    #calculate covariance
    K <- kernelFn(sigma_f, L, X, nX)
    K_c <- K[obs,]
    K_tilde <- K[obs,obs]
    diag(K_tilde) <- diag(K_tilde) + sigma_y^2

    # SVD because matrix close to singular
    svd_Ktilde <- svd(K_tilde)
    Kdot <- solve(svd_Ktilde$u,
                  diag(1/svd_Ktilde$d)%*%solve(t(svd_Ktilde$v), K_c))

    # posterior moments
    post_var <- K - crossprod(K_c, solve(K_tilde, K_c))
    post_mean <- crossprod(K_c, solve(K_tilde,y))


    #posterior samples
    z <- matrix(rnorm(iter * length(post_mean)), nrow = length(post_mean), ncol=iter)
    f_star <-  post_mean + t(chol(post_var)) %*% z

    # posterior predictive
    y_star <- f_star + matrix(rnorm(iter*(nX),
                                    sd=sigma_y), nrow=iter, ncol=nX)

    # distribution summary
    quant_ystar <- t(apply(y_star, 2, function(x)
      c(mean(x), quantile(x, prob=c(0.02,0.975)))))

    return(list(posterior=t(f_star), post_pred=t(y_star), obs=obs, not_obs=not_obs, summary=quant_ystar, X=X))
  }



  f_dens <- function(f, K) {
    dmvnorm(f, rep(0,length(f)), K, log=TRUE)
  }

  likelihood <- function(Y, f, sigma_y) {
    dnorm(Y, mean=f, sd=sigma_y, log=TRUE)
  }

  L_prior <- function(L) {
    dgamma(L, 10, 10, log=TRUE) + L
  }

  sigma_f_prior <- function(sig) {
    dgamma(sig, 10, 10, log=TRUE)
  }

  L_sample <- function(L_current, sd_adapt){
    L_new <- exp(rnorm(1, log(L_current), sd_adapt))
  }

  kern_dens <- function(f, K, L, sigma) {
    f_dens(f,K) + L_prior(L) + sigma_f_prior(sigma)
  }

  L_update <- function(f,exp_cost2, K, L2_old, R, sigma2, sd_adapt) {
    dens_old <- kern_dens(f, K, L2_old, sigma2)
    L2_new <- L_sample(L2_old, sd_adapt)
    KR_new <- kernelinvLFn(sigma2, L2_new, cost2, nrow(K))
    dens_new <- kern_dens(f, KR_new$K, L2_new, sigma2)
    if(log(runif(1)) <= dens_new - dens_old) {
      return(list(K=KR_new$K, L=L2_new, R=KR_new$R))
    } else {
      return(list(K=K, L=L2_old, R=R))
    }
  }

  kernelinvLFn <- function(sigma2, L2, exp_cost2, n) {
    R <- Rupdate(L2, exp_cost2, n)
    Sigma <- (R)/sigma2
    if(!isSymmetric(Sigma)){
      upper.tri(Sigma) <- lower.tri(Sigma)
    }
    return(list(R= R, K=Sigma))
  }

  Rupdate <- function(L2, exp_cost2 ,n) {
    # xmat <- matrix(x, ncol=n,nrow=n)
    # xmatT <- matrix(x, ncol=n, nrow=n, byrow = TRUE)
    R <- exp_cost2 *exp(L2) + diag(0.001, nrow(R)) # exp(-0.5*L2 * cost2)
    diag(R) <- diag(R) + 0.001
    return(R)
  }

  Kupdate <- function(sigma2, R) {
    Sigma <- (R )/sigma2
    return(Sigma)
  }


  invGammaSamp <- function(n, alpha, beta){
    1/rgamma(n, alpha, beta)
  }

  varUpdate <- function(alpha, beta, y, mu,n){
    alpha_star <- alpha + n * 0.5
    beta_star <- beta + sum((y-mu)^2)/2
    return(invGammaSamp(1, alpha_star, beta_star))
  }

  sigmafUpdate <- function(alpha, beta, f, K,n,sigma2){
    alpha_star <- alpha + n * 0.5
    beta_star <- beta + (crossprod(f, solve(K*sigma2, f)))/2
    return(rgamma(1, alpha_star, beta_star))
  }

  f_sample <- function(K,obs, sigma_y, y) {
    #calculate covariance
    K_c <- K[obs,]
    K_tilde <- K[obs,obs]
    diag(K_tilde) <- diag(K_tilde) + sigma_y^2

    # posterior moments
    post_var <- K - crossprod(K_c, solve(K_tilde, K_c))
    post_mean <- c(crossprod(K_c, solve(K_tilde,y)))

    #samples
    f_star <- c(rmvnorm(1, post_mean, post_var))
    return(f_star)
  }

  mse <- function(f, obs, y,burnin=NULL, iter=NULL){
    if(is.null(burnin)) burnin <- 0
    if(is.null(iter)) iter <- ncol(f)

    f_obs <- f[obs,(burnin+1):iter]

    mean(apply(f_obs, 2, function(x) (x-y)^2))
  }

  GP <- function(x, y, grid = seq(max(x), min(x), length.out = 100),
                 hyperparameters = list(), iter=2000, burnin=NULL,
                 display.progress = TRUE)
  {
    if(is.null(burnin)) burnin = iter*0.5
    adapt <- burnin * 0.1
    nsamps <- iter - burnin

    # pull out hyperparam
    alpha_f <- hyperparameters$alpha_f
    beta_f <- hyperparameters$beta_f
    alpha_y <- hyperparameters$alpha_y
    beta_y <- hyperparameters$beta_y

    if(is.null(alpha_f)) alpha_f <- 1
    if(is.null(beta_f)) beta_f <- 1
    if(is.null(alpha_y)) alpha_y <- 1
    if(is.null(beta_y)) beta_y <- 1

    #make united X vector
    X <- c(grid,x)
    nX <- length(X)
    exp_cost2 <- exp(-0.5* cost_calc(t(X),t(X),2)^2 )

    #no which ones are observed data
    obs <- seq(length(grid)+1,nX,1)
    not_obs <- 1:length(grid)

    #store length of observed data
    n <- length(y)

    # generate initial variables
    param <- list(
      L = rgamma(1,alpha_f,beta_f),
      sigma_f = rgamma(1,alpha_f,beta_f),
      sigma_y = rgamma(1,alpha_y,beta_y),
      R = NULL,
      K = NULL,
      f = NULL
    )
    genKR   <- kernelinvLFn(param$sigma_f, param$L, exp_cost2, nX)
    param$K <- genKR$K
    param$R <- genKR$R
    param$f <- f_sample(param$K,obs, param$sigma_y,y)
    sd_adapt <- 2.38


    #setup posterior sample storage
    store <- list(
      L = rep(NA, iter),
      sigma_f = rep(NA, iter),
      sigma_y = rep(NA, iter),
      f = matrix(NA, nrow=nX, ncol=iter)
    )

    # temporary values
    L_temp <- NULL
    L_adapt <- rep(NA, adapt)

    # progress bar
    if(display.progress){
      pb <- txtProgressBar(min = 0, max = iter, style = 3)
    }

    for(i in 1:iter) {
      L_temp  <- L_update(param$f, exp_cost2, param$K,
                          param$L, param$R, param$sigma_f,
                          sd_adapt)
      param$L <- L_temp$L
      param$K <- L_temp$K
      param$R <- L_temp$R
      param$sigma_f <- sigmafUpdate(alpha_f, beta_f, param$f,
                                    param$K,nX, param$sigma_f)
      param$K <- Kupdate(param$sigma_f, param$R)
      param$f <- f_sample(param$K,obs, param$sigma_y, y)
      param$sigma_y <- varUpdate(alpha_y, beta_y, y, param$f[obs], n)

      store$L[i] <- param$L
      store$sigma_f[i] <- param$sigma_f
      store$sigma_y[i] <- param$sigma_y
      store$f[,i] <- param$f

      if(i < adapt){
        L_adapt[i] <- param$L
        if(i %% 20 == 0){
          sd_adapt <- 0.75*sd(log(L_adapt[(i-19):i])) + 0.25*sd_adapt
        }
      }
      if(display.progress) setTxtProgressBar(pb, i)
    }
    if(display.progress) close(pb)
    # posterior predictive
    store$y_star <- store$f + matrix(rnorm(iter*(nX),
                                           sd=store$sigma_y), ncol=iter, nrow=nX)

    # distribution summary
    quant_ystar <- t(apply(store$y_star[,(burnin+1):iter],
                           1, function(x) c(mean(x),
                                            quantile(x, prob=c(0.025,0.975)))))

    return(list(param=store, obs=obs, not_obs=not_obs,
                summary=quant_ystar, X=X, MHsd=sd_adapt,
                burnin=burnin, iter=iter))
  }

  return(list(rpost = GP, rprior = GP_prior))
}
