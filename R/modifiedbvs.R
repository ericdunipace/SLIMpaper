bvsmod <- function (X, resp, prep = TRUE, logT = FALSE, fixed_cols = NULL,
          eff_size = 0.5, family = c("logistic", "survival"), hselect = TRUE,
          nlptype = "piMOM", r = 1, tau = 0.25, niter = 30, mod_prior = c("unif",
                                                                          "beta"), inseed = NULL, cplng = FALSE, ncpu = 4, parallel.MPI = FALSE)
{
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  if (!class(X) == "data.frame")
    stop("input X should be a data frame!")
  ol <- BVSNLP:::matprep(X, fixed_cols, prep, logT)
  X <- ol$fulmat
  gnames <- ol$gnames
  nf <- ol$nf
  if (family == "logistic") {
    y <- as.numeric(resp)
    dx <- dim(X)
    n <- dx[1]
    p <- dx[2]
    X <- cbind(rep(1, n), X)
    gnames <- c("Intercept", gnames)
    cons <- 0
    prp <- p/n
    ar <- 2^n
    if (prp > 4 && ar < Inf) {
      ac <- 0
      cons <- 0
      while (ar > ac) {
        cons <- cons + 1
        ac <- choose(p, cons)
      }
    }
    else {
      cons <- ceiling(log(p))
    }
    cons <- min(cons, ceiling(log(p)))
    if (mod_prior == "beta") {
      a <- cons
      b <- p - a
    }
    if (mod_prior == "unif") {
      a <- 1
      b <- 1
    }
    if (hselect) {
      hyper <- BVSNLP:::HyperSelect(X, y, eff_size, nlptype, 20000,
                           mod_prior, family)
      tau <- hyper$tau
      r <- 1
    }
    initProb <- cons/p
    exmat <- cbind(y, X)
    if (nlptype == "piMOM")
      nlptype_int <- 0
    if (nlptype == "pMOM")
      nlptype_int <- 1
    if (!cplng) {
      schain <- p
      while (schain > cons || schain == 0) {
        chain1 <- rbinom(p - nf, 1, initProb)
        schain <- sum(chain1)
      }
      chain1 <- as.numeric(c(rep(1, nf + 1), chain1))
      chain2 <- chain1
      Lregout <- BVSNLP:::logreg_bvs(exmat, chain1, nf, tau, r,
                            nlptype_int, a, b, cons, niter, cplng, chain2)
      Hash_Key <- Lregout$hash_key
      all_probs <- Lregout$hash_prob
      VisCovs <- Lregout$vis_covs
      inds <- which(all_probs != 0)
      Hash_Key <- Hash_Key[inds]
      all_probs <- all_probs[inds]
      VisCovs <- VisCovs[inds]
      nvm <- length(unique(Hash_Key))
      uinds <- which(!duplicated(Hash_Key))
      all_probs <- all_probs[uinds]
      list_vis_covs <- VisCovs[uinds]
      outnum <- min(nvm, 1000)
      sout <- sort(all_probs, decreasing = T, index.return = T)
      MaxMargs <- sout$x[1:outnum]
      minds <- sout$ix[1:outnum]
      max_marg <- MaxMargs[1]
      indmax <- minds[1]
      sel_model <- list_vis_covs[[indmax]]
      gnames2 <- gnames[sel_model + 1]
      beta_hat <- Lregout$beta_hat
      names(beta_hat) <- gnames2
      gnames <- gnames[-1]
      sel_model <- sel_model[-1]
      MaxModels <- list(NULL)
      for (i in 1:outnum) {
        MaxModels[[i]] <- list_vis_covs[[minds[i]]][-1]
      }
      inc_probs <- inc_prob_calc(all_probs, list_vis_covs,
                                 p + 1)
      inc_probs <- inc_probs[-1]
      median_model <- which(inc_probs >= 0.5)
      return(list(max_prob = max_marg, HPM = sel_model,
                  beta_hat = beta_hat, MPM = median_model, inc_probs = inc_probs,
                  max_prob_vec = MaxMargs, max_models = MaxModels,
                  num_vis_models = nvm, nlptype = nlptype, des_mat = X,
                  gene_names = gnames, r = r, tau = tau))
    }
    else {
      comb <- function(x, ...) {
        lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...),
                                                          function(y) y[[i]])))
      }
      if (parallel.MPI) {
        if (!requireNamespace("doMPI", quietly = TRUE)) {
          stop("Package doMPI needed for this function to work. Please install it.",
               call. = FALSE)
        }
        else {
          cl <- doMPI::startMPIcluster(count = ncpu)
          doMPI::registerDoMPI(cl)
          parout <- foreach::foreach(j = 1:ncpu, .combine = "comb",
                            .multicombine = TRUE, .init = list(list(),
                                                               list(), list(), list()), .packages = "BVSNLP",
                            .options.mpi = list(seed = inseed)) %dopar%
            {
              schain <- p
              while (schain > cons || schain == 0) {
                chain1 <- rbinom(p - nf, 1, initProb)
                schain <- sum(chain1)
              }
              chain1 <- as.numeric(c(rep(1, nf + 1),
                                     chain1))
              schain <- p
              while (schain > cons || schain == 0) {
                chain2 <- rbinom(p - nf, 1, initProb)
                schain <- sum(chain2)
              }
              chain2 <- as.numeric(c(rep(1, nf + 1),
                                     chain2))
              Lregout <- BVSNLP:::logreg_bvs(exmat, chain1, nf,
                                    tau, r, nlptype_int, a, b, cons, niter,
                                    cplng, chain2)
              maxChain <- as.logical(Lregout$max_chain)
              maxMarg <- Lregout$max_prob
              cflag <- Lregout$cplng_flag
              bhat <- numeric(p + 1)
              bhat[maxChain] <- Lregout$beta_hat
              list(maxChain, maxMarg, cflag, bhat)
            }
          doMPI::closeCluster(cl)
        }
      }
      else {
        cl <- parallel::makeCluster(ncpu)
        doParallel::registerDoParallel(cl)
        opts <- list(preschedule = TRUE)
        if (!is.null(inseed)) {
          parallel::clusterSetRNGStream(cl, inseed)
        }
        ParOut <- foreach::foreach(j = 1:ncpu, .combine = "comb",
                          .multicombine = TRUE, .init = list(list(),
                                                             list(), list(), list()), .packages = "BVSNLP",
                          .options.snow = opts) %dopar% {
                            schain <- p
                            while (schain > cons || schain == 0) {
                              chain1 <- rbinom(p - nf, 1, initProb)
                              schain <- sum(chain1)
                            }
                            chain1 <- as.numeric(c(rep(1, nf + 1), chain1))
                            schain <- p
                            while (schain > cons || schain == 0) {
                              chain2 <- rbinom(p - nf, 1, initProb)
                              schain <- sum(chain2)
                            }
                            chain2 <- as.numeric(c(rep(1, nf + 1), chain2))
                            Lregout <- BVSNLP:::logreg_bvs(exmat, chain1, nf, tau,
                                                  r, nlptype_int, a, b, cons, niter, cplng,
                                                  chain2)
                            maxChain <- as.logical(Lregout$max_chain)
                            maxMarg <- Lregout$max_prob
                            cflag <- Lregout$cplng_flag
                            bhat <- numeric(p + 1)
                            bhat[maxChain] <- Lregout$beta_hat
                            list(maxChain, maxMarg, cflag, bhat)
                          }
        stopCluster(cl)
      }
      MaxChain <- matrix(unlist(ParOut[[1]]), ncol = (p +
                                                        1), byrow = T)
      MaxMarg <- unlist(ParOut[[2]])
      cpl_flag <- unlist(ParOut[[3]])
      bhat <- matrix(unlist(ParOut[[4]]), ncol = (p + 1),
                     byrow = T)
      cpl_percent <- sum(cpl_flag)/ncpu
      Final_Marg <- MaxMarg
      Final_Chains <- MaxChain
      D <- as.data.frame(cbind(Final_Chains, Final_Marg))
      Counts <- rep(1, length(Final_Marg))
      A <- aggregate(Counts, by = as.list(D), FUN = sum)
      Freq <- A[, p + 3]
      Probs <- A[, p + 2]
      UniqModels <- apply(A[, 1:(p + 1)], 1, function(x) which(x >
                                                                 0))
      return(list(cpl_percent = cpl_percent, margin_probs = Final_Marg,
                  chains = Final_Chains, cpl_flags = cpl_flag,
                  beta_hat = bhat, freq = Freq, probs = Probs,
                  uniq_models = UniqModels, nlptype = nlptype,
                  gene_names = gnames, r = r, tau = tau))
    }
  }
  if (family == "survival") {
    TS <- resp
    time <- TS[, 1]
    status <- TS[, 2]
    sfidx <- nf + 1
    dx <- dim(X)
    n <- dx[1]
    p <- dx[2]
    exmat <- cbind(time, status, X)
    if (nlptype == "piMOM")
      nlptype_int <- 0
    if (nlptype == "pMOM")
      nlptype_int <- 1
    cons <- 1 + nf
    if (mod_prior == "beta") {
      a <- cons
      b <- p - a
    }
    if (mod_prior == "unif") {
      a <- 1
      b <- 1
    }
    if (hselect) {
      hyper <- BVSNLP:::HyperSelect(X, TS, eff_size, nlptype, 5000,
                           mod_prior, family)
      tau <- hyper$tau
      r <- 1
    }
    ntimes <- 10
    d <- 2 * ceiling(log(p))
    temps <- seq(3, 1, length.out = ntimes)
    L <- ntimes
    J <- niter
    comb <- function(x, ...) {
      lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...),
                                                        function(y) y[[i]])))
    }
    if (parallel.MPI) {
      if (!requireNamespace("doMPI", quietly = TRUE)) {
        stop("Package doMPI needed for this function to work. Please install it.",
             call. = FALSE)
      }
      else {
        cl <- doMPI::startMPIcluster(count = ncpu)
        doMPI::registerDoMPI(cl)
        parout <- foreach(j = 1:ncpu, .combine = "comb",
                          .multicombine = TRUE, .init = list(list(),
                                                             list(), list(), list(), list(), list()),
                          .packages = "BVSNLP", .options.mpi = list(seed = inseed)) %dopar%
          {
            cur_model <- sample(sfidx:p, floor(p/100))
            if (nf > 0)
              cur_model <- c(1:nf, cur_model)
            coxout <- BVSNLP:::cox_bvs(exmat, cur_model, nf, tau,
                              r, nlptype_int, a, b, d, L, J, temps)
            maxmod <- coxout$max_model
            maxprob <- coxout$max_prob
            hashkey <- coxout$hash_key
            allprobs <- coxout$all_probs
            viscovs <- coxout$vis_covs_list
            list(maxmod, maxprob, hashkey, allprobs,
                 cur_model, viscovs)
          }
        doMPI::closeCluster(cl)
      }
    }
    else {
      cl <- parallel::makeCluster(ncpu)
      doParallel::registerDoParallel(cl)
      opts <- list(preschedule = TRUE)
      if (!is.null(inseed)) {
        parallel::clusterSetRNGStream(cl, inseed)
      }
      parout <- foreach::foreach(j = 1:ncpu, .combine = "comb",
                        .multicombine = TRUE, .init = list(list(), list(),
                                                           list(), list(), list(), list()), .packages = "BVSNLP",
                        .options.snow = opts) %dopar% {
                          cur_model <- sample(sfidx:p, floor(p/100))
                          if (nf > 0)
                            cur_model <- c(1:nf, cur_model)
                          coxout <- BVSNLP:::cox_bvs(exmat, cur_model, nf, tau,
                                            r, nlptype_int, a, b, d, L, J, temps)
                          maxmod <- coxout$max_model
                          maxprob <- coxout$max_prob
                          hashkey <- coxout$hash_key
                          allprobs <- coxout$all_probs
                          viscovs <- coxout$vis_covs_list
                          list(maxmod, maxprob, hashkey, allprobs, cur_model,
                               viscovs)
                        }
      parallel::stopCluster(cl)
    }
    Hash_Key <- unlist(parout[[3]])
    All_Probs <- unlist(parout[[4]])
    CurModel <- matrix(unlist(parout[[5]]), ncol = 3, byrow = T)
    VisCovs <- NULL
    for (i in 1:ncpu) {
      VisCovs <- c(VisCovs, parout[[6]][[i]])
    }
    num_vis_models <- length(unique(Hash_Key))
    uinds <- which(!duplicated(Hash_Key))
    all_probs <- All_Probs[uinds]
    list_vis_covs <- VisCovs[uinds]
    outnum <- min(num_vis_models, 1000)
    sout <- sort(all_probs, decreasing = T, index.return = T)
    MaxMargs <- sout$x[1:outnum]
    minds <- sout$ix[1:outnum]
    max_marg <- MaxMargs[1]
    indmax <- minds[1]
    sel_model <- list_vis_covs[[indmax]] + 1
    MaxModels <- list(NULL)
    for (i in 1:outnum) {
      MaxModels[[i]] <- list_vis_covs[[minds[i]]] + 1
    }
    inc_probs <- BVSNLP:::inc_prob_calc(all_probs, list_vis_covs,
                               p)
    median_model <- which(inc_probs >= 0.5)
    beta_hat <- BVSNLP:::CoefEst(X, TS, sel_model, nlptype, tau, r,
                        "survival")
    return(list(num_vis_models = num_vis_models, max_prob = max_marg,
                HPM = sel_model, MPM = median_model, beta_hat = beta_hat,
                max_prob_vec = MaxMargs, max_models = MaxModels,
                inc_probs = inc_probs, nlptype = nlptype, des_mat = X,
                start_models = CurModel, r = r, tau = tau, gene_names = gnames,
                all_probs = all_probs, vis_mod = list_vis_covs))
  }
}
