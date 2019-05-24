extractCoef <- function(fit, drop=FALSE) {

  lambda <- unlist(fit$lambda)
  nlambda <- length(lambda)

  if(any(class(fit) == "glmnet")){
    nvar <- fit$df+1
    extractCoefVar <- tapply(lambda, nvar, min)
    coefs <- coef(fit, s = extractCoefVar)
    idx <- which(lambda %in% extractCoefVar)
  } else {
    if(any(class(fit)=="oem")){
      nvar <- colSums(fit$beta[[1]] != 0)
      extractCoefVar <- tapply(lambda, nvar, min)
      idx <- which(lambda %in% extractCoefVar)
      coefs <- fit$beta[[1]][,idx,drop=FALSE]
    }
    if(any(class(fit)=="sparse-posterior")){
      nvar <- colSums(fit$beta != 0)
      extractCoefVar <- tapply(lambda, nvar, min)
      iters <- fit$niter
      maxit <- eval(fit$call$maxit)
      lastHitMax <- (iters[nrow(iters),nlambda]>=maxit)
      if(lastHitMax == TRUE & sum(nvar==max(nvar))>1) {
        last.idx <- which(nvar==max(nvar))
        hitMax <- ifelse(iters[nrow(iters),last.idx]>=maxit,1,0)
        cummulative <- cumsum(hitMax)
        diffs <- diff(c(0,cummulative))
        rep.idx <- max(last.idx[diffs == 0])
        extractCoefVar[length(extractCoefVar)] <- lambda[rep.idx]
      }
      # hitInf <- is.infinite(colSums(fit$beta))
      # if(lastHitMax == TRUE & sum(nvar==max(nvar))>1) {
      #   last.idx <- which(nvar==max(nvar))
      #   hitMax <- ifelse(iters[nrow(iters),last.idx]>=maxit,1,0)
      #   cummulative <- cumsum(hitMax)
      #   diffs <- diff(c(0,cummulative))
      #   rep.idx <- max(last.idx[diffs == 0])
      #   extractCoefVar[length(extractCoefVar)] <- lambda[rep.idx]
      # }
      idx <- which(lambda %in% extractCoefVar)

      coefs <- fit$beta[,idx,drop=FALSE]

      orders <- order(nvar[idx])
      coefs <- coefs[,orders,drop=FALSE]
      idx <- idx[orders]
    }
  }
  selLamb <- lambda[idx]
  selVar <- nvar[idx]
  if(drop!=FALSE) coefs <- coefs[-drop,,drop=FALSE]
  return(list(coefs = coefs, lambda = selLamb, nzero = selVar))
}

annealCoef <- function(fit, theta) {
  idx <- lapply(fit, function(f) f$optimal$index)
  numActive <- sapply(idx, length)
  gamma <- lapply(fit, function(f) f$optimal$gamma)
  p <- ncol(theta)
  force <- fit[[1]]$force

  coefs <- matrix(0, nrow=p, ncol=p)
  coefs[force,1:length(force)] <- 1
  coefs[,p] <- 1

  for(j in seq_along(numActive)){
    jj <- numActive[j]
    coefs[idx[[j]], jj] <- gamma[[j]]
  }
  numActive <- c(length(force),numActive,p)

  return(list(coefs = coefs, nzero = numActive))
}

stepCoef <- function(fit, theta) {
  idx <- fit$index
  p <- ncol(theta)

  coefs <- sapply(idx, function(i) {
    ones <- rep(0, p)
    ones[i] <- 1
    return(ones)
  })
  return(list(coefs = coefs, nzero = fit$numActive))
}

extractTheta <- function(fit, theta) {
  if(any(class(fit)=="optimization")) {
    coefficient_trajectory <- extractCoef(fit)
    method <- fit$method
  } else if ( any(class(fit)=="stepwise" ) ) {
    coefficient_trajectory <- annealCoef(fit)
    method <- "stepwise"
  } else if ( any(class(fit)=="annealing" )) {
    coefficient_trajectory <- stepCoef(fit)
    method <- "annealing"
  }

  nzero <- coefficient_trajectory$nzero
  n.coefsets <- length(nzero)
  p <- nrow(theta)
  s <- ncol(theta)

  theta_list <- vector("list", n.coefsets)

  if ( method == "projection" ) {
    for ( i in seq_along(theta_list) ) {
      theta_list[[i]] <- matrix(coefficient_trajectory$coef[,i],p,s)
    }
    nzero <- nzero/s
  } else {
    is_loc_scale <- (method == "location.scale")

    if ( is_loc_scale ) {
      mean_theta <- matrix(rowMeans(theta),p, s)
      c_theta    <- theta - mean_theta
      theta <- rbind(c_theta, mean_theta)
      nzero <- nzero/2
    }

    for ( i in seq_along(theta_list) ) {
      theta_list[[i]] <- diag(coefficient_trajectory$coef[,i] ) %*% theta
      if ( is_loc_scale ) {
        theta_list[[i]] <- theta_list[[i]][1:p,] + theta_list[[i]][( p + 1):(2 * p),]
      }
    }
  }

  return(list(theta=theta_list, nzero=nzero))

}
