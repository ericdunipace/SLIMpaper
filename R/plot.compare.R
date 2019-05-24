plot.compare <- function(models, target = NULL, X, theta, method = c("w2", "mse"), quantity=c("posterior","mean"), parallel=FALSE, transform = function(x){return(x)}) {

  method <- match.arg(method)
  quantity <- match.arg(quantity,several.ok = TRUE)
  if(length(quantity)>1 & method == "mse") stop("Can only do one quantity with mse")
  if(parallel) {
    require(doParallel)
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    doParallel::registerDoParallel(cl)
  }
  mu_fun <- function(tclist, X){
    if(parallel) {
      mu_coarse <- foreach::foreach(tc=tclist$theta) %dopar% {X %*% tc}
    } else {
      mu_coarse <- lapply(tclist$theta, function(tc) X %*% tc)
    }

    return(list(mu = mu_coarse, nzero=tclist$nzero))
  }

  dist_fun <- function(mulist, method, mu) {
    dist <-  if(method == "w2") {
      if(parallel) {
         foreach::foreach(m=mulist, .combine = c, .export='transform') %dopar% {
           mm <- as.matrix(transform(t(m)))
           costm <- distance(mm, as.matrix(mu))
           mass <- rep(1,nrow(mu))
           if(any(is.infinite(costm^2) | is.nan(costm))) {
             return(NA)
           } else {
             return(transport::wasserstein(mass, mass, p=2, tplan=NULL, costm=costm, method="shortsimplex"))
           }
         }
      } else {
        sapply(mulist, function(m) {
          mm <- as.matrix(transform(t(m)))
          costm <- distance(mm, as.matrix(mu))
          mass <- rep(1,nrow(mu))
          if(any(is.infinite(costm^2) | is.nan(costm))) {
            return(NA)
          } else {
            return(transport::wasserstein(mass, mass, p=2, tplan=NULL, costm=costm, method="shortsimplex"))
          }
        })
      }

    } else if (method == "mse"){
      if(parallel) {
        foreach::foreach(m=mulist, .combine = c, .export='transform') %dopar% {

          mm <- as.matrix(transform(t(m)))
          return(mean((mm - as.matrix(mu))^2))
          }
      } else {
        sapply(mulist, function(m) {
          mm <- as.matrix(transform(t(m)))
          return(mean((mm - as.matrix(mu))^2))
          }
        )
      }
    }
    return(dist)
  }

  method <- match.arg(method)
  if (!is.list(models)) models <- list(models)

  # theta <- models[[1]]$call$theta
  # X <- models[[1]]$call$X
  n <- nrow(X)
  p <- ncol(X)


  if ( p != nrow(theta) ) theta <- t(theta)
  if ( is.null(target) ) target <- models[[1]]$Y
  if ( is.null(target) ) target <- X %*% theta
  s <- ncol(theta)

  if (method == "mse" & (is.vector(target) | any(dim(target)[2]==1))) {
    target <- if(quantity == "mean") {
      matrix(target, n, s)
    } else {
      matrix(target, p, s)
    }

  }

  dist_df <- dist_mu_df <- nactive <- groups <- plot <- plot_mu <- NULL

  theta_coarse <- lapply(models, extractTheta, theta=theta)

  if (method == "w2") {
    ylab <- "2-Wasserstein Distance"
  } else if (method == "mse") {
    ylab <- "MSE"
  }

  if("posterior" %in% quantity){
    dist_list <- if ( method == "w2" ){
       lapply(theta_coarse, function(mc) dist_fun(mc$theta, method=method, mu=t(theta)))
    } else if ( method == "mse" ) {
       lapply(theta_coarse, function(mc) dist_fun(mc$theta, method=method, mu=t(target)))
    }

    dist <- unlist(dist_list)
    nactive <- sapply(theta_coarse, function(d) d$nzero)
    groups <- mapply(function(x,z){return(rep(x, each=z))}, x=names(models), z=sapply(dist_list, length))

    dist_df <- data.frame(dist = dist,
                          nactive = unlist(nactive),
                          groups=factor(unlist(groups)))
    if(all(is.na(dist_df$dist))) {
      plot <- NULL
    } else {
      plot <- ggplot( dist_df, aes(x=nactive, y=dist, color = groups, group=groups )) +
        geom_line() + scale_color_jama() + labs(color ="Method") +
        xlab("Number of active coefficients") + ylab(ylab) + theme_bw() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    }

  }

  if("mean" %in% quantity){
    mu_coarse <- lapply(theta_coarse, function(tc) mu_fun(tc, X=X))
    if(any(dim(target) == dim(mu_coarse[[1]][[1]][[1]]))) {
      if(nrow(target) == nrow(mu_coarse[[1]][[1]][[1]])) target <- t(target)
    }
    dist_list_mu <- lapply(mu_coarse, function(mc) dist_fun(mc$mu, method=method, mu=target))

    dist_mu <- unlist(dist_list_mu)
    if(is.null(nactive)) nactive <- sapply(theta_coarse, function(d) d$nzero)

    if(is.null(groups)) groups <- mapply(function(x,z){ return(rep(x, each=z)) }, x=names(models), z=sapply(dist_list_mu, length))

    dist_mu_df <- data.frame(dist = dist_mu,
                          nactive = unlist(nactive),
                          groups=factor(unlist(groups)))
    if ( all(is.na(dist_mu_df$dist))){
      plot_mu <- NULL
    } else {
      plot_mu <- ggplot( dist_mu_df, aes(x=nactive, y=dist, color = groups, group=groups )) +
        geom_line() + scale_color_jama() + labs(color ="Method") +
        xlab("Number of active coefficients") + ylab(ylab) + theme_bw() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    }

  }

  if (parallel) parallel::stopCluster(cl)
  plots <- list(posterior = plot, mean = plot_mu)
  data <- list(posterior = dist_df, mean = dist_mu_df)

  return(list(plot = plots, data = data))
}
