abart.fix <- function (x.train, times, delta, x.test = matrix(0, 0, 0), K = 100,
          type = "abart", ntype = 1, sparse = FALSE, theta = 0, omega = 1,
          a = 0.5, b = 1, augment = FALSE, rho = NULL, xinfo = matrix(0,
                                                                      0, 0), usequants = FALSE, rm.const = TRUE, sigest = NA,
          sigdf = 3, sigquant = 0.9, k = 2, power = 2, base = 0.95,
          lambda = NA, tau.num = c(NA, 3, 6)[ntype], offset = NULL,
          w = rep(1, length(times)), ntree = c(200L, 50L, 50L)[ntype],
          numcut = 100L, ndpost = 1000L, nskip = 100L, keepevery = c(1L,
                                                                     10L, 10L)[ntype], printevery = 100L, transposed = FALSE,
          mc.cores = 1L, nice = 19L, seed = 99L)
{
  if (type != "abart")
    stop("type must be \"abart\"")
  if (ntype != 1)
    stop("ntype must be 1")
  y.train = log(times)
  n = length(y.train)
  if (n != length(delta))
    stop("length of times and delta must be equal")
  delta = as.integer(delta)
  if (!transposed) {
    temp = BART:::bartModelMatrix(x.train, numcut, usequants = usequants,
                           xinfo = xinfo, rm.const = rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if (length(x.test) > 0) {
      x.test = BART:::bartModelMatrix(x.test)
      x.test = t(x.test[, temp$rm.const])
    }
    rm.const <- temp$rm.const
    grp <- temp$grp
    rm(temp)
  }
  else {
    rm.const <- NULL
    grp <- NULL
  }
  if (n != ncol(x.train))
    stop("The length of times and the number of rows in x.train must be identical")
  p = nrow(x.train)
  np = ncol(x.test)
  if (length(rho) == 0)
    rho = p
  if (length(rm.const) == 0)
    rm.const <- 1:p
  if (length(grp) == 0)
    grp <- 1:p
  if (length(offset) == 0) {
    offset = mean(y.train)
  }
  if (type == "abart") {
    y.train = y.train - offset
    if (is.na(lambda)) {
      if (is.na(sigest)) {
        if (p < n)
          sigest = summary(lm(y.train ~ ., data.frame(t(x.train),
                                                      y.train)))$sigma
        else sigest = sd(y.train)
      }
      qchi = qchisq(1 - sigquant, sigdf)
      lambda = (sigest^2) * qchi/sigdf
    }
    else {
      sigest = sqrt(lambda)
    }
    if (is.na(tau.num)) {
      tau = (max(y.train) - min(y.train))/(2 * k * sqrt(ntree))
    }
    else {
      tau = tau.num/(k * sqrt(ntree))
    }
  }
  else {
    lambda = 1
    sigest = 1
    tau = tau.num/(k * sqrt(ntree))
  }
  ptm <- proc.time()
  res = .Call("cabart", ntype, n, p, np, x.train, y.train,
              delta, x.test, ntree, numcut, ndpost * keepevery, nskip,
              keepevery, power, base, offset, tau, sigdf, lambda,
              sigest, w, sparse, theta, omega, grp, a, b, rho, augment,
              printevery, xinfo, PACKAGE = "BART")
  res$proc.time <- proc.time() - ptm
  K <- min(n, K)
  events = unique(sort(times))
  if (length(events) > K) {
    events <- unique(quantile(times, probs = (1:K)/K))
    attr(events, "names") <- NULL
  }
  K <- length(events)
  if (type == "abart") {
    res$surv.train <- matrix(nrow = ndpost, ncol = n * K)
    for (i in 1:n) for (j in 1:K) {
      h <- (i - 1) * K + j
      res$surv.train[, h] <- pnorm(log(events[j]),
                                   mean = res$yhat.train[,i],
                                   sd = res$sigma[seq(nskip+keepevery,keepevery*ndpost, length.out = ndpost)],
                                   lower.tail = FALSE)
    }
    res$yhat.train.mean <- apply(res$yhat.train, 2, mean)
    res$surv.train.mean <- apply(res$surv.train, 2, mean)
  }
  else {
    if (type == "pbart")
      res$prob.train = pnorm(res$yhat.train)
    else if (type == "lbart")
      res$prob.train = plogis(res$yhat.train)
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
  }
  if (np > 0) {
    if (type == "abart") {
      res$surv.test <- matrix(nrow = ndpost, ncol = np *
                                K)
      for (i in 1:np) for (j in 1:K) {
        h <- (i - 1) * K + j
        res$surv.test[, h] <- pnorm(log(events[j]),
                                    mean = res$yhat.test[, i],
                                    sd = res$sigma[seq(nskip+keepevery,keepevery*ndpost, length.out = ndpost)],
                                    lower.tail = FALSE)
      }
      res$yhat.test.mean <- apply(res$yhat.test, 2, mean)
      res$surv.test.mean <- apply(res$surv.test, 2, mean)
    }
    else {
      if (type == "pbart")
        res$prob.test = pnorm(res$yhat.test)
      else if (type == "lbart")
        res$prob.test = plogis(res$yhat.test)
      res$prob.test.mean <- apply(res$prob.test, 2, mean)
    }
  }
  res$times = events
  res$K = K
  res$offset = offset
  names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
  dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const
  attr(res, "class") <- type
  return(res)
}

mc.abart.fix <- function (x.train, times, delta, x.test = matrix(0, 0, 0), K = 100,
                          type = "abart", ntype = 1, sparse = FALSE, theta = 0, omega = 1,
                          a = 0.5, b = 1, augment = FALSE, rho = NULL, xinfo = matrix(0,
                                                                                      0, 0), usequants = FALSE, rm.const = TRUE, sigest = NA,
                          sigdf = 3, sigquant = 0.9, k = 2, power = 2, base = 0.95,
                          lambda = NA, tau.num = c(NA, 3, 6)[ntype], offset = NULL,
                          w = rep(1, length(times)), ntree = c(200L, 50L, 50L)[ntype],
                          numcut = 100L, ndpost = 1000L, nskip = 100L, keepevery = c(1L,
                                                                                     10L, 10L)[ntype], printevery = 100L, transposed = FALSE,
                          mc.cores = 2L, nice = 19L, seed = 99L)
{
  if (type != "abart")
    stop("type must be \"abart\"")
  if (ntype != 1)
    stop("ntype must be 1")
  if (.Platform$OS.type != "unix")
    stop("parallel::mcparallel/mccollect do not exist on windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  if (!transposed) {
    temp = bartModelMatrix(x.train, numcut, usequants = usequants,
                           xinfo = xinfo, rm.const = rm.const)
    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    if (length(x.test) > 0) {
      x.test = bartModelMatrix(x.test)
      x.test = t(x.test[, temp$rm.const])
    }
    rm.const <- temp$rm.const
    rm(temp)
  }
  mc.cores.detected <- parallel::detectCores()
  if (mc.cores > mc.cores.detected)
    mc.cores <- mc.cores.detected
  mc.ndpost <- ceiling(ndpost/mc.cores)
  for (i in 1:mc.cores) {
    parallel::mcparallel({
      tools::psnice(value = nice)
      abart.fix(x.train = x.train, times = times, delta = delta,
            x.test = x.test, K = K, type = type, ntype = ntype,
            sparse = sparse, theta = theta, omega = omega,
            a = a, b = b, augment = augment, rho = rho, xinfo = xinfo,
            usequants = usequants, rm.const = rm.const, sigest = sigest,
            sigdf = sigdf, sigquant = sigquant, k = k, power = power,
            base = base, lambda = lambda, tau.num = tau.num,
            offset = offset, w = w, ntree = ntree, numcut = numcut,
            ndpost = mc.ndpost, nskip = nskip, keepevery = keepevery,
            printevery = printevery, transposed = TRUE)
    }, silent = (i != 1))
  }
  post.list <- parallel::mccollect()
  post <- post.list[[1]]
  if (mc.cores == 1 | attr(post, "class") != type)
    return(post)
  else {
    if (class(rm.const)[1] != "logical")
      post$rm.const <- rm.const
    post$ndpost <- mc.cores * mc.ndpost
    p <- nrow(x.train[post$rm.const, ])
    old.text <- paste0(as.character(mc.ndpost), " ", as.character(ntree),
                       " ", as.character(p))
    old.stop <- nchar(old.text)
    post$treedraws$trees <- sub(old.text, paste0(as.character(post$ndpost),
                                                 " ", as.character(ntree), " ", as.character(p)),
                                post$treedraws$trees)
    keeptest <- length(x.test) > 0
    for (i in 2:mc.cores) {
      post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)
      post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)
      post$surv.train <- rbind(post$surv.train, post.list[[i]]$surv.train)
      if (keeptest) {
        post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)
        post$surv.test <- rbind(post$surv.test, post.list[[i]]$surv.test)
      }
      post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
      post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
      post$treedraws$trees <- paste0(post$treedraws$trees,
                                     substr(post.list[[i]]$treedraws$trees, old.stop +
                                              2, nchar(post.list[[i]]$treedraws$trees)))
      post$proc.time["elapsed"] <- max(post$proc.time["elapsed"],
                                       post.list[[i]]$proc.time["elapsed"])
      for (j in 1:5) if (j != 3)
        post$proc.time[j] <- post$proc.time[j] + post.list[[i]]$proc.time[j]
    }
    if (type == "abart") {
      post$yhat.train.mean <- apply(post$yhat.train, 2,
                                    mean)
      post$surv.train.mean <- apply(post$surv.train, 2,
                                    mean)
      if (keeptest) {
        post$yhat.test.mean <- apply(post$yhat.test,
                                     2, mean)
        post$surv.test.mean <- apply(post$surv.test,
                                     2, mean)
      }
    }
    else {
      post$prob.train.mean <- apply(post$prob.train, 2,
                                    mean)
      if (keeptest)
        post$prob.test.mean <- apply(post$prob.test,
                                     2, mean)
    }
    post$varcount.mean <- apply(post$varcount, 2, mean)
    post$varprob.mean <- apply(post$varprob, 2, mean)
    attr(post, "class") <- type
    return(post)
  }
}


mc.surv.bart.fix <- function (x.train = matrix(0, 0, 0), y.train = NULL, times = NULL,
          delta = NULL, x.test = matrix(0, 0, 0), K = NULL, events = NULL,
          ztimes = NULL, zdelta = NULL, sparse = FALSE, theta = 0,
          omega = 1, a = 0.5, b = 1, augment = FALSE, rho = NULL, xinfo = matrix(0,
                                                                                 0, 0), usequants = FALSE, rm.const = TRUE, type = "pbart",
          ntype = as.integer(factor(type, levels = c("wbart", "pbart",
                                                     "lbart"))), k = 2, power = 2, base = 0.95, offset = NULL,
          tau.num = c(NA, 3, 6)[ntype], ntree = 50L, numcut = 100L,
          ndpost = 1000L, nskip = 250L, keepevery = 10L, printevery = 100L,
          id = NULL, seed = 99L, mc.cores = 2L, nice = 19L)
{
  if (.Platform$OS.type != "unix")
    stop("parallel::mcparallel/mccollect do not exist on windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  if (is.na(ntype) || ntype == 1)
    stop("type argument must be set to either 'pbart' or 'lbart'")
  x.train <- bartModelMatrix(x.train)
  x.test <- bartModelMatrix(x.test)
  if (length(y.train) == 0) {
    pre <- BART::surv.pre.bart(times, delta, x.train, x.test, K = K,
                         events = events, ztimes = ztimes, zdelta = zdelta)
    y.train <- pre$y.train
    x.train <- pre$tx.train
    x.test <- pre$tx.test
  }
  else {
    if (length(unique(sort(y.train))) > 2)
      stop("y.train has >2 values; make sure you specify times=times & delta=delta")
  }
  H <- 1
  Mx <- 2^31 - 1
  Nx <- max(nrow(x.train), nrow(x.test))
  if (Nx > Mx%/%ndpost) {
    H <- ceiling(ndpost/(Mx%/%Nx))
    ndpost <- ndpost%/%H
  }
  mc.cores.detected <- parallel::detectCores()
  if (mc.cores > mc.cores.detected) {
    message("The number of cores requested, ", mc.cores,
            ",\n exceeds the number of cores detected via detectCores() ",
            "reducing to ", mc.cores.detected)
    mc.cores <- mc.cores.detected
  }
  mc.ndpost <- ceiling(ndpost/mc.cores)
  post.list <- list()
  for (h in 1:H) {
    for (i in 1:mc.cores) {
      parallel::mcparallel({
        tools::psnice(value = nice)
        surv.bart(x.train = x.train, y.train = y.train,
                  x.test = x.test, sparse = sparse, theta = theta,
                  omega = omega, a = a, b = b, augment = augment,
                  rho = rho, xinfo = xinfo, usequants = usequants,
                  rm.const = rm.const, type = type, k = k, power = power,
                  base = base, offset = offset, tau.num = tau.num,
                  ntree = ntree, numcut = numcut, ndpost = mc.ndpost,
                  nskip = nskip, keepevery = keepevery, printevery = printevery,
                  id = id)
      }, silent = (i != 1))
    }
    post.list[[h]] <- parallel::mccollect()
  }
  if ((H == 1 & mc.cores == 1) | attr(post.list[[1]][[1]],
                                      "class") != "survbart")
    return(post.list[[1]][[1]])
  else {
    for (h in 1:H) for (i in mc.cores:1) {
      if (h == 1 & i == mc.cores) {
        post <- post.list[[1]][[mc.cores]]
        post$ndpost <- H * mc.cores * mc.ndpost
        p <- ncol(x.train[, post$rm.const])
        old.text <- paste0(as.character(mc.ndpost), " ",
                           as.character(ntree), " ", as.character(p))
        old.stop <- nchar(old.text)
        post$treedraws$trees <- sub(old.text, paste0(as.character(post$ndpost),
                                                     " ", as.character(ntree), " ", as.character(p)),
                                    post$treedraws$trees)
      }
      else {
        if (length(x.test) > 0) {
          post$yhat.test <- rbind(post$yhat.test, post.list[[h]][[i]]$yhat.test)
          post$prob.test <- rbind(post$prob.test, post.list[[h]][[i]]$prob.test)
          post$surv.test <- rbind(post$surv.test, post.list[[h]][[i]]$surv.test)
        }
        post$varcount <- rbind(post$varcount, post.list[[h]][[i]]$varcount)
        post$varprob <- rbind(post$varprob, post.list[[h]][[i]]$varprob)
        post$treedraws$trees <- paste0(post$treedraws$trees,
                                       substr(post.list[[h]][[i]]$treedraws$trees,
                                              old.stop + 2, nchar(post.list[[h]][[i]]$treedraws$trees)))
        post$proc.time["elapsed"] <- max(post$proc.time["elapsed"],
                                         post.list[[h]][[i]]$proc.time["elapsed"])
        for (j in 1:5) if (j != 3)
          post$proc.time[j] <- post$proc.time[j] + post.list[[h]][[i]]$proc.time[j]
      }
      post.list[[h]][[i]] <- NULL
    }
    if (length(x.test) > 0)
      post$surv.test.mean <- apply(post$surv.test, 2, mean)
    post$varcount.mean <- apply(post$varcount, 2, mean)
    post$varprob.mean <- apply(post$varprob, 2, mean)
    attr(post, "class") <- "survbart"
    return(post)
  }
}
