LBP <- function(X, eta, theta,
    solver,
    nvars,
    infm_maxit = 100L,
    tol = 1e-5,
    solver.options = list(verbose = 0L),
    parallel = NULL
) {
  require(doRNG)
  power <- 2
  solver <- match.arg(solver, c("mosek"))

  X <- as.matrix(X)
  Y <- as.matrix(eta)
  theta <- as.matrix(theta)

  n <- nrow(X)
  p <- ncol(X)

  if (missing(nvars) || is.null(nvars)) {
    nvars <- 1:p
  }
  nvars <- as.integer(nvars)
  p_star <- length(nvars)

  # rotate vars accordingly
  if(ncol(theta) == ncol(X)){
    theta_ <- t(theta)
  } else {
    theta_ <- theta
  }
  if(nrow(theta_) != p) stop("dimensions of theta must match X")
  theta_save <- theta_

  #transpose X
  X_ <- t(X)

  same <- FALSE
  if(is.null(Y)) {
    same <- TRUE
    Y_ <- crossprod(X_,theta_)
  } else{
    if(!any(dim(Y) %in% dim(X_))) stop("dimensions of Y must match X")
    if(!is.matrix(Y)) Y <- as.matrix(Y)
    if(nrow(Y) == ncol(X_)){
      # print("Transpose")
      Y_ <- Y
    } else{
      Y_ <- t(Y)
    }
    if(all(Y_==crossprod(X_, theta_))) same <- TRUE
  }
  if(ncol(Y_) != ncol(theta_)) stop("ncol of Y should be same as ncols of theta")
  if(nrow(Y_) != ncol(X_)) stop("The number of observations in Y and X don't line up. Make sure X is input with observations in rows.")
  rmv.idx <- NULL
  if(any(apply(theta_,1, function(x) all(x == 0)))) {
    rmv.idx <- which(apply(theta_,1, function(x) all(x == 0)))

    X_ <- X_[-rmv.idx, ]
    theta_ <- theta_[-rmv.idx,]
    penalty.factor <- penalty.factor[-rmv.idx]
    warning("Some dimensions of theta have no variation. These have been removed")
  }

  # setup vars
  transport.method <- "exact"
  infm.maxit <- as.integer(100L)
  epsilon <- as.double(0.05)
  OTmaxit <- as.integer(0L)

  OToptions <- list(same = same,
                    method = "selection.variable",
                    transport.method = transport.method,
                    epsilon = epsilon,
                    niter = OTmaxit)


  ss <- WpProj:::sufficientStatistics(X, Y_, theta_, OToptions)
  xtx <- ss$XtX
  xty <- xty_init <- ss$XtY
  Ytemp <- Y_
  # lam_max <- log(max(abs(xty)))
  # lambda <- exp(seq(lam_max, log(1e-4) + lam_max, length.out = nlambda))


  # obj function
  obj <- function(x, xtx, xty) {
    c(0.5 * t(x) %*% (2 * xtx) %*% x - 2 * xty %*% x)
  }
  to_feasible <- function(x) {
    x[x < 0] <- 0
    x[x > 1] <- 1
    round(x)
  }
  # setup problem
  cc        <- c(-2 * xty)
  Upper <- chol(2 * xtx)
  A   <- as.matrix(t(c(rep(1,p))))
  b   <- 1L
  prob <- Rmosek::mosek_qptoprob(F = Upper, f = cc,
                                 Aeq = A, beq = b,
                         lb = rep(0,p),
                         ub = rep(1, p))
  contprob <- prob
  contprob$bx[1,1:p] <- -Inf
  contprob$bx[2,1:p] <- Inf

  prob$intsub <- 1:p
  prob$c <- c(prob$c, 0, 0)
  num.vars <- length(prob$c)
  prob$A  <- rbind(
    # add slack var `t` such that `t = b-Ax`, then in objective
    # we have `min_t lambda * t`
    cbind(prob$A, c(1, rep(0, nrow(prob$A) - 1L)),0)
    , c(rep(0, num.vars- 2L), 1, -1) # `t <= u`
    , c(rep(0, num.vars- 2L), 1,  1)) # `-u <= t`
  prob$bx <- cbind(prob$bx, c(-Inf, Inf), c(0, Inf))
  prob$bc <- cbind(prob$bc, c(-Inf, 0), c(0, Inf))

  # l2_cone <- list(
  #   type = "MSK_CT_QUAD",
  #   inds = c(num.vars, num.vars - 1L)
  # )
  #
  # prob$cones <- cbind(prob$cones, l2_cone)


  # optimization loop
  z_lp <- ret_cont <- ret <- k <- alpha <- rho <- lambda <- NULL
  alpha <- alpha_old <- test_feasible <- rep(0,p)

  output <- list(beta = list(), niter = list())
  output$beta <- vector("list", length(nvars))
  output$niter <- vector("list", length(nvars))

  if(!is.null(parallel)){
    if(!inherits(parallel, "cluster") && !is.numeric(parallel)) {
      stop("parallel must be a registered cluster backend or the number of cores desired")
    }
    doParallel::registerDoParallel(parallel)
    display.progress <- FALSE
  } else{
    foreach::registerDoSEQ()

  }

  comb <- function(x, ...) {
    # from https://stackoverflow.com/questions/19791609/saving-multiple-outputs-of-foreach-dopar-loop
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }

  prob_init <- prob
  contprob_init <- contprob

  output <- foreach::foreach(k=nvars, .combine='comb', .multicombine=TRUE,
                             .init=list(list(), list()),
                             .errorhandling = 'pass',
                             .inorder = FALSE) %dorng%
    {
  # for(k in nvars) {
      prob <- prob_init
      contprob <- contprob_init
      xty <- xty_init
      z_lp <- ret_cont <- ret <- alpha <- rho <- lambda <- NULL
      alpha <- alpha_old <- test_feasible <- rep(0,p)
      prob$bc[,1] <- c(k,k)
      test_feasible[1:k] <- 1
      if(k == p) {
        alpha <- rep(1,p)
        inf <- 0L
      } else {
        for(inf in 1:infm_maxit) {
          ret_cont <- Rmosek::mosek(contprob, opts = solver.options)
          alpha_cont <- ret_cont$sol$itr$xx[1:p]
          z_lp <- obj(alpha_cont, xtx, c(xty))

          rho  <- max(obj(test_feasible, xtx, c(xty)) - z_lp, 0) # delta = 1 in this case
          lambda <- max(ret_cont$sol$itr$suc[1], ret_cont$sol$itr$slc[1])

          prob$c[length(prob$c )] <- rho #max(rho * 10, lambda)
          prob$c[length(prob$c ) - 1] <- lambda
          ret <- Rmosek::mosek(prob, opts = solver.options)
          alpha <- round(ret$sol$int$xx[1:p])


          obj_val <- obj(alpha, xtx, c(xty))

          if(WpProj:::not.converged(alpha, alpha_old, tol) ||
             WpProj:::not.converged(obj_val, obj_val_old, tol)){
            alpha_old <- alpha
            obj_val_old <- obj_val

            Ytemp <- WpProj:::selVarMeanGen(X_, theta_, as.double(alpha))
            xty   <- WpProj:::xtyUpdate(X, Ytemp, theta_,
                                        result_ = as.double(alpha),
                                        OToptions)
            contprob$c[1:p] <- -2 * xty
            prob$c[1:p] <- -2 * xty

          } else {
            break
          }

        }
      }
      return(list(beta = alpha, niter = inf))
      # output$beta[[idx]] <- alpha
      # output$niter[[idx]] <- inf
  }
  names(output) <- c("beta","niter")
  output$beta <- do.call("cbind", output$beta)
  output$niter <- unlist(output$niter)
  output[c("xtx", "xty_init","xty_final")] <- list(xtx, xty_init, xty)

  output$nvars <- p
  output$varnames <- colnames(X)
  output$call <- formals(LBP)
  output$remove.idx <- rmv.idx
  output$nonzero_beta <- colSums(output$beta != 0)
  # output$nzero <- nz
  class(output) <- c("WpProj","IP")

  extract <- WpProj:::extractTheta(output, theta_)
  output$nzero <- extract$nzero
  output$eta <- lapply(extract$theta, function(tt) crossprod(X_, tt))
  output$theta <- extract$theta
  if(!is.null(rmv.idx)) {
    for(i in seq_along(output$theta)){
      output$theta[[i]] <- theta_save
      output$theta[[i]][-rmv.idx,] <- extract$theta[[i]]
    }
  }

  return_value <- list(
    call = output$call,
    theta = output$theta,
    fitted.values = output$eta,
    power = 2.0,
    method = "binary program",
    solver = solver,
    niter = output$niter,
    nzero = output$nzero
  )

  class(return_value) <- c("WpProj","LBP")
  return(return_value)

}
