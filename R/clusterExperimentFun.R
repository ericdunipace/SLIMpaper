experimentW2Method <- function(target, hyperparameters, conditions, w2=FALSE) {
  n <- conditions$n
  p <- conditions$p

  n.samps <- conditions$n.samps
  penalty <- conditions$penalty
  lambda.min.ratio <- conditions$lambda.min.ratio
  n.lambda <- conditions$n.lambda
  penalty_method <- conditions$penalty.factor
  family <- conditions$family
  pseudo.obs <- conditions$pseudo.obs
  stan_dir <- conditions$stan_dir
  if(is.null(pseudo.obs)) pseudo.obs <- 0
  posterior.method <- conditions$posterior.method

  # SETUP PARAMETERS
  param <- target$rparam()
  p_star <- min(length(param$theta),p)
  param$theta <- param$theta[1:p_star]

  #set up data
  X <- target$X$rX(n, target$X$corr, p)
  Y <- target$rdata(n, X[,1:p_star,drop=FALSE], c(param$theta),
                    param$sigma2)

  X_sing <- matrix(c(target$X$rX(1, target$X$corr, p)), nrow=1, ncol=p)
  single_Y <- target$rdata(1, X_sing[,1:p_star,drop=FALSE], c(param$theta),
                           param$sigma2)

  X_new <- target$X$rX(n, target$X$corr, p)
  new_Y_new_X <- target$rdata(n, X_new[,1:p_star,drop=FALSE], c(param$theta), param$sigma2)

  #True means
  true_mu <- X[,1:p_star,drop=FALSE] %*% c(param$theta)
  new_mu <- X_new[,1:p_star,drop=FALSE] %*% c(param$theta)
  new_mu_sing <- X_sing[,1:p_star,drop=FALSE] %*% c(param$theta)

  #sample theta
  post_sample <- target$rpost(n.samps, X, Y, hyperparameters, method = posterior.method, stan_dir = stan_dir)
  theta <- post_sample$theta #regression coef
  t_theta <- t(theta)
  sigma <- post_sample$sigma #variance (if it exists for model)

  #functions of theta
  E_theta <- rowMeans(theta)
  theta_norm <- rowSums(theta^2)

  #penalty terms
  penalty_fact <- set_penalty_factor(t_theta, penalty_method)
  HC_penalty_fact <- set_penalty_factor(t_theta, "expectation")

  #conditional and marginal means
  cond_mu <- X %*% theta
  marg_mu <- X %*% E_theta

  cond_mu_new <- tcrossprod(X_new, t_theta)
  marg_mu_new <- X_new %*% E_theta

  cond_mu_sing <- tcrossprod(X_sing,t_theta)
  marg_mu_sing <- X_sing %*% c(E_theta)

  #optional Wasserstein
  if(w2){
    w2s <- w2l0(X,theta)
  }

  augDat <- augPseudo(X, cond_mu, t_theta, theta_norm, pseudo.obs, n, same=TRUE)
  lambdas <- calc.lambdas(augDat, lambda.min.ratio, penalty_fact, n.lambda)

  #New Method
  time <- proc.time()
  # augDat <- augmentMatrix(X, cond_mu, t_theta, same=TRUE)
  selDist <- W2L1(X, NULL, theta, family=family, penalty="selection.lasso",
                  penalty.factor=penalty_fact, nlambda = n.lambda,
                  lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                  maxit=1000,
                  pseudo_observations = pseudo.obs,
                  display.progress=TRUE, gamma = 1)
  selTime <- proc.time() - time
  trajSel <- extractCoef(selDist)

  #permutation
  time <- proc.time()
  permDist <- W2L1(X, NULL, theta, family=family, penalty=penalty,
                   penalty.factor=penalty_fact, nlambda = n.lambda,
                   lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                   maxit=1000,
                   pseudo_observations = pseudo.obs,
                   display.progress=TRUE, gamma = 1)
  trueTime <- proc.time() - time
  trajPerm <- extractCoef(permDist)


  #carvalho method
  lassoHC <- oem.xtx(crossprod(X)/n, crossprod(X,marg_mu)/n,
                     family=family, penalty=penalty,
                     penalty.factor=HC_penalty_fact, nlambda = n.lambda,
                     lambda.min.ratio = lambda.min.ratio, lambda=lambdas, gamma=1)
  trajHC <- extractCoef(lassoHC)

  trajHCdist <- list(coefs = NULL, nzero = trajHC$nzero)
  trajHCdist$coefs<- matrix(0, nrow=p, ncol = length(trajHC$nzero))
  HC_non_zero_idx <- which(trajHC$coefs != 0)
  trajHCdist$coefs[HC_non_zero_idx] <- 1

  #stepwise
  time <- proc.time()
  step <- w2.stepwise(X, theta, force=1, direction = "backward", gamma.method = "identity", pseudo.obs = pseudo.obs)
  trajStep <- stepCoef(step, t_theta)
  stepTime <- proc.time() - time

  #simulated annealing
  time <- proc.time()
  anneal <- lapply(2:(p-1), function(j)
    W2SA(X=X, Y=NULL, theta=t_theta,
         force = 1, model.size = j, iter = 10*j, temps = 100,  pseudo.obs = pseudo.obs,
         options = list(gamma.method = "identity", energy.distribution = "boltzman",
                        cooling.schedule="exponential"))
  )
  trajAnneal <- annealCoef(anneal, t_theta)
  annealTime <- proc.time() - time

  #w2 from "true" posterior mean
  sel_theta <- extractCoef(selDist, theta)
  perm_theta <- extractCoef(permDist, theta)
  step_theta <- extractCoef(step, theta)
  anneal_theta <- extractCoef(anneal, theta)
  hc_theta <- list()
  hc_theta$theta <- lapply(trajHC, function(h) diag(h) %*% theta)
  hc_theta$nzero <- trajHC$nzero

  sel_mu <- lapply(sel_theta$theta, function(tt) X %*% tt)
  perm_mu <- lapply(perm_theta$theta, function(tt) X %*% tt)
  step_mu <- lapply(step_theta$theta, function(tt) X %*% tt)
  anneal_mu <- lapply(anneal_theta$theta, function(tt) X %*% tt)
  hc_mu <- lapply(hc_theta$theta, function(tt) X %*% tt)

  wass_trajectory(sel_theta, t_theta, sel_theta$nzero, p)

  W2_sel <- wass_trajectory(sel_mu, t(cond_mu), sel_theta$nzero, p)
  W2_perm <- wass_trajectory(perm_mu, t(cond_mu), perm_theta$nzero, p)
  W2_step <- wass_trajectory(step_mu, t(cond_mu), step_theta$nzero, p)
  W2_anneal <- wass_trajectory(anneal_mu, t(cond_mu), anneal_theta$nzero, p)
  W2_HC <- wass_trajectory(hc_mu, t(cond_mu), hc_theta$nzero, p)

  # W2_sel <- W2_idx(trajSel$coefs, X, t_theta, cond_mu, trajSel$nzero, p)
  # W2_perm <- W2_idx(trajPerm$coefs, X, t_theta, cond_mu, trajPerm$nzero, p)
  # W2_step <- W2_idx(trajStep$coefs, X, t_theta, cond_mu, trajStep$nzero, p)
  # W2_anneal <- W2_idx(trajAnneal$coefs, X, t_theta, cond_mu, trajAnneal$nzero, p)
  # W2_HC <- W2_idx(trajHCdist$coefs, X, t_theta, cond_mu, trajHCdist$nzero, p)

  #mse on means of original data
  mean_mseSel <- mse_idx_dist(trajSel$coefs, X, t_theta, true_mu, trajSel$nzero,p)
  mean_mseSelE <- mse_idx_expect(trajSel$coefs, X, t_theta, true_mu, trajSel$nzero, p)
  mean_msePerm <- mse_idx_dist(trajPerm$coefs, X, t_theta, true_mu, trajPerm$nzero,p)
  mean_msePermE <- mse_idx_expect(trajPerm$coefs, X, t_theta, true_mu, trajPerm$nzero,p)
  mean_mseHC <- mse_idx_dist(trajHCdist$coefs, X, t_theta, true_mu, trajHCdist$nzero,p)
  mean_mseHCE <- mse_idx(trajHC$coefs, X, true_mu, trajHC$nzero, p)
  mean_mseStep <- mse_idx_dist(trajStep$coefs, X, t_theta, true_mu, trajStep$nzero,p)
  mean_mseAnneal <- mse_idx_dist(trajAnneal$coefs, X, t_theta, true_mu, trajAnneal$nzero,p)


  #mse on new outcome data from same paramters and same X
  new_Y <- target$rdata(n, X, c(param$theta),
                        param$sigma2)
  # newY_mseSel <- mse_idx_dist(trajSel$coefs, X, t_theta, new_Y, trajSel$nzero, p)
  # newY_mseSelE <- mse_idx_expect(trajSel$coefs, X, t_theta, new_Y, trajSel$nzero, p)
  # newY_msePerm <- mse_idx_dist(trajPerm$coefs, X, t_theta, new_Y, trajPerm$nzero,p)
  # newY_msePermE <- mse_idx_expect(trajPerm$coefs, X, t_theta, new_Y, trajPerm$nzero,p)
  # newY_mseHC <- mse_idx_dist(trajHCdist$coefs, X, t_theta, new_Y, trajHCdist$nzero,p)
  # newY_mseHCE <- mse_idx(trajHC$coefs, X, new_Y, trajHC$nzero, p)
  # newY_mseStep <- mse_idx_dist( trajStep$coefs, X, t_theta, new_Y, trajStep$nzero, p )
  # newY_mseAnneal <- mse_idx_dist( trajAnneal$coefs, X, t_theta, new_Y, trajAnneal$nzero, p )


  #mse on new outcome data from same paramters and different X
  #new method
  augDatN <- augPseudo(X_new, cond_mu_new, t_theta, theta_norm, pseudo.obs, n, same=TRUE)
  lambdas <- calc.lambdas(augDatN, lambda.min.ratio, penalty_fact, n.lambda)

  lassoSelN <- W2L1(X_new, NULL, theta, family=family, penalty="selection.lasso",
                    penalty.factor=penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                    maxit=1000,
                    pseudo_observations = pseudo.obs, lambda = lambdas,
                    display.progress=TRUE, gamma = 1)
  trajSelN <- extractCoef(lassoSelN)

  #carvalho method
  lassoHCN <- oem.xtx(crossprod(X_new)/n, crossprod(X_new,marg_mu_new)/n,
                      family=family, penalty=penalty,
                      penalty.factor=penalty_fact, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, lambda=lambdas)
  trajHCN <- extractCoef(lassoHCN)

  trajHCNdist <- list(coefs = NULL, nzero = trajHCN$nzero)
  trajHCNdist$coefs<- matrix(0, nrow=p, ncol = length(trajHCN$nzero))
  HCN_non_zero_idx <- which(trajHCN$coefs != 0)
  trajHCNdist$coefs[HCN_non_zero_idx] <- 1

  #permutation
  permDistN <- W2L1(X_new, NULL, theta, family=family, penalty=penalty,
                    penalty.factor=penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                    maxit=1000,
                    pseudo_observations = pseudo.obs, lambda = lambdas,
                    display.progress=TRUE)
  trajPermN <- extractCoef(permDistN)

  #stepwise
  stepN <- w2.stepwise(X_new, t_theta, force=1, direction = "backward", gamma.method = "identity", pseudo.obs = pseudo.obs)
  trajStepN <- stepCoef(stepN, t_theta)

  #simulated annealing
  annealN <- lapply(2:(p-1), function(j)
    W2SA(X=X_new, Y=NULL, theta=t_theta,
         force = 1, model.size = j, iter = 10*j, temps = 100,  pseudo.obs = pseudo.obs,
         options = list(gamma.method = "identity", energy.distribution = "boltzman",
                        cooling.schedule="exponential"))
  )
  trajAnnealN <- annealCoef(annealN, t_theta)

  #new method, single datapoint
  augDatO <- augPseudo(X_sing, cond_mu_sing, t_theta, theta_norm, pseudo.obs, n, same=TRUE)
  lambdas <- calc.lambdas(augDatO, lambda.min.ratio, penalty_fact, n.lambda)

  lassoSelO <- W2L1(X_sing, NULL, theta, family=family, penalty="selection.lasso",
                    penalty.factor=penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                    maxit=1000,
                    pseudo_observations = pseudo.obs, lambda = lambdas,
                    display.progress=TRUE, gamma = 1)
  trajSelO <- extractCoef(lassoSelO)

  #carvalho method, single datapoint
  lassoHCO <- oem.xtx(crossprod(X_sing)/n, crossprod(X_sing,marg_mu_sing)/n,
                      family=family, penalty=penalty,
                      penalty.factor=penalty_fact, nlambda = n.lambda,
                      lambda.min.ratio = lambda.min.ratio, lambda=lambdas)
  trajHCO <- extractCoef(lassoHCO)

  trajHCOdist <- list(coefs = NULL, nzero = trajHCO$nzero)
  trajHCOdist$coefs<- matrix(0, nrow=p, ncol = length(trajHCO$nzero))
  HCO_non_zero_idx <- which(trajHCO$coefs != 0)
  trajHCOdist$coefs[HCO_non_zero_idx] <- 1

  #permutation, single datapoint
  permDistO <- W2L1(X_sing, NULL, theta, family=family, penalty=penalty,
                    penalty.factor=penalty_fact, nlambda = n.lambda,
                    lambda.min.ratio = lambda.min.ratio, infimum.maxit=10000,
                    maxit=1000,
                    pseudo_observations = pseudo.obs,
                    display.progress=TRUE)
  trajPermO <- extractCoef(permDistO)

  #stepwise
  stepO <- w2.stepwise(X_sing, t_theta, force=1, direction = "backward", gamma.method = "identity", pseudo.obs = pseudo.obs)
  trajStepO <- stepCoef(stepO, t_theta)

  #simulated annealing
  annealO <- lapply(2:(p-1), function(j)
    W2SA(X=X_sing, Y=NULL, theta=t_theta,
         force = 1, model.size = j, iter = 10*j, temps = 100,  pseudo.obs = pseudo.obs,
         options = list(gamma.method = "identity", energy.distribution = "boltzman",
                        cooling.schedule="exponential"))
  )
  trajAnnealO <- annealCoef(annealO, t_theta)

  #new mean
  sel_thetaN <- extractCoef(selDistN, theta)
  perm_theta <- extractCoef(permDistN, theta)
  step_thetaN <- extractCoef(stepN, theta)
  anneal_thetaN <- extractCoef(annealN, theta)
  hc_thetaN <- list()
  hc_thetaN$theta <- lapply(trajHCN, function(h) diag(h) %*% theta)
  hc_thetaN$nzero <- trajHCN$nzero

  sel_muN <- lapply(sel_thetaN$theta, function(tt) X %*% tt)
  perm_muN <- lapply(perm_thetaN$theta, function(tt) X %*% tt)
  step_muN <- lapply(step_thetaN$theta, function(tt) X %*% tt)
  anneal_muN <- lapply(anneal_thetaN$theta, function(tt) X %*% tt)
  hc_muN <- lapply(hc_thetaN$theta, function(tt) X %*% tt)

  W2_selN <- wass_trajectory(sel_muN, t(cond_muN), sel_thetaN$nzero, p)
  W2_permN <- wass_trajectory(perm_muN, t(cond_mu), perm_thetaN$nzero, p)
  W2_stepN <- wass_trajectory(step_muN, t(cond_mu), step_thetaN$nzero, p)
  W2_annealN <- wass_trajectory(anneal_muN, t(cond_mu), anneal_thetaN$nzero, p)
  W2_HCN <- wass_trajectory(hc_muN, t(cond_mu), hc_thetaN$nzero, p)

  newMu_mseSel <- mse_idx_dist(trajSelN$coefs, X_new, t_theta, new_mu, trajSelN$nzero, p)
  newMu_mseSelE <- mse_idx_expect(trajSelN$coefs, X_new, t_theta, new_mu, trajSelN$nzero, p)
  newMu_msePerm <- mse_idx_dist(trajPermN$coefs, X_new, t_theta, new_mu, trajPermN$nzero,p)
  newMu_msePermE <- mse_idx_expect(trajPermN$coefs, X_new, t_theta, new_mu, trajPermN$nzero,p)
  newMu_mseHC <- mse_idx_dist(trajHCNdist$coefs, X_new, t_theta, new_mu, trajHCNdist$nzero,p)
  newMu_mseHCE <- mse_idx(trajHCN$coefs, X_new, new_mu, trajHCN$nzero, p)
  newMu_mseStep <- mse_idx_dist(trajStepN$coefs, X_new, t_theta, new_mu, trajStepN$nzero, p)
  newMu_mseAnneal <- mse_idx_dist(trajAnnealN$coefs, X_new, t_theta, new_mu, trajAnnealN$nzero, p)

  #new data
  # newX_mseSel <- mse_idx_dist(trajSelN$coefs, X_new, t_theta, new_Y_new_X, trajSelN$nzero, p)
  # newX_mseSelE <- mse_idx_expect(trajSelN$coefs, X_new, t_theta, new_Y_new_X, trajSelN$nzero, p)
  # newX_msePerm <- mse_idx_dist(trajPermN$coefs, X_new, t_theta, new_Y_new_X, trajPermN$nzero,p)
  # newX_msePermE <- mse_idx_expect(trajPermN$coefs, X_new, t_theta, new_Y_new_X, trajPermN$nzero,p)
  # newX_mseHC <-  mse_idx(trajHCNdist$coefs, X_new, new_Y_new_X, trajHCNdist$nzero, p)
  # newX_mseHCE <-  mse_idx(trajHCN$coefs, X_new, new_Y_new_X, trajHCN$nzero, p)
  # newX_mseStep <- mse_idx_dist(trajStepN$coefs, X_new, t_theta, new_Y_new_X, trajStepN$nzero, p)
  # newX_mseAnneal <- mse_idx_dist(trajAnnealN$coefs, X_new, t_theta, new_Y_new_X, trajAnnealN$nzero, p)


  #single new mean
  sel_thetaO <- extractCoef(selDistO, theta)
  perm_theta <- extractCoef(permDistO, theta)
  step_thetaO <- extractCoef(stepO, theta)
  anneal_thetaO <- extractCoef(annealO, theta)
  hc_thetaO <- list()
  hc_thetaO$theta <- lapply(trajHCO, function(h) diag(h) %*% theta)
  hc_thetaO$nzero <- trajHCO$nzero

  sel_muO <- lapply(sel_thetaO$theta, function(tt) X %*% tt)
  perm_muO <- lapply(perm_thetaO$theta, function(tt) X %*% tt)
  step_muO <- lapply(step_thetaO$theta, function(tt) X %*% tt)
  anneal_muO <- lapply(anneal_thetaO$theta, function(tt) X %*% tt)
  hc_muO <- lapply(hc_thetaO$theta, function(tt) X %*% tt)

  W2_selO <- wass_trajectory(sel_muO, t(cond_muO), sel_thetaO$nzero, p)
  W2_permO <- wass_trajectory(perm_muO, t(cond_mu), perm_thetaO$nzero, p)
  W2_stepO <- wass_trajectory(step_muO, t(cond_mu), step_thetaO$nzero, p)
  W2_annealO <- wass_trajectory(anneal_muO, t(cond_mu), anneal_thetaO$nzero, p)
  W2_HCO <- wass_trajectory(hc_muO, t(cond_mu), hc_thetaO$nzero, p)

  singleMu_mseSel <- mse_idx_dist(trajSelO$coefs, X_sing, t_theta, new_mu_sing, trajSelO$nzero, p)
  singleMu_mseSelE <- mse_idx_expect(trajSelO$coefs, X_sing, t_theta, new_mu_sing, trajSelO$nzero, p)
  singleMu_msePerm <- mse_idx_dist(trajPermO$coefs, X_sing, t_theta, new_mu_sing, trajPermO$nzero,p)
  singleMu_msePermE <- mse_idx_expect(trajPermO$coefs, X_sing, t_theta, new_mu_sing, trajPermO$nzero,p)
  singleMu_mseHC <- mse_idx_dist(trajHCOdist$coefs, X_sing, t_theta, new_mu_sing, trajHCOdist$nzero,p)
  singleMu_mseHCE <- mse_idx(trajHCO$coefs, X_sing, new_mu_sing, trajHCO$nzero, p)
  singleMu_mseStep <- mse_idx_dist(trajStepO$coefs, X_sing, t_theta, new_mu_sing, trajStepO$nzero, p)
  singleMu_mseAnneal <- mse_idx_dist(trajAnnealO$coefs, X_sing, t_theta, new_mu_sing, trajAnnealO$nzero, p)

  #single new data
  # singleY_mseSel <- mse_idx_dist(trajSelO$coefs, X_sing, t_theta, single_Y, trajSelO$nzero, p)
  # singleY_mseSelE <- mse_idx_expect(trajSelO$coefs, X_sing, t_theta, single_Y, trajSelO$nzero, p)
  # singleY_msePerm <- mse_idx_dist(trajPermO$coefs, X_sing, t_theta, single_Y, trajPermO$nzero,p)
  # singleY_msePermE <- mse_idx_expect(trajPermO$coefs, X_sing, t_theta, single_Y, trajPermO$nzero,p)
  # singleY_mseHC <- mse_idx_dist(trajHCOdist$coefs, X_sing, t_theta, single_Y, trajHCOdist$nzero,p)
  # singleY_mseHCE <-  mse_idx(trajHCO$coefs, X_sing, single_Y, trajHCO$nzero, p)
  # singleY_mseStep <- mse_idx_dist(trajStepO$coefs, X_sing, t_theta, single_Y, trajStepO$nzero, p)
  # singleY_mseAnneal <- mse_idx_dist(trajAnnealO$coefs, X_sing, t_theta, single_Y, trajAnnealO$nzero, p)

  #### Using W2 to measure closeness to full posterior ####
  X_ident <- diag(1,p,p)
  postW2_Sel <- wass_trajectory(sel_theta, t_theta, sel_theta$nzero, p)
  postW2_Perm <- wass_trajectory(perm_theta, t_theta, perm_theta$nzero, p)
  postW2_Step <- wass_trajectory(step_theta, t_theta, step_theta$nzero, p)
  postW2_Anneal <- wass_trajectory(anneal_theta, t_theta, anneal_theta$nzero, p)
  postW2_HC <- wass_trajectory(hc_theta, t_theta, hc$nzero, p)

  #### generate output list ####

  outList <- list (
    W2_dist = list(sel = W2_sel, perm = W2_perm, step = W2_step, anneal = W2_anneal,
                   HC = W2_HC),
    W2_dist_post = list(sel = postW2_Sel, perm = postW2_Perm,
                        step = postW2_Step, anneal = postW2_Anneal,
                        HC = postW2_HC),
    W2_new = list(sel = W2_selN, perm = W2_permN, step = W2_stepN, anneal = W2_annealN,
                   HC = W2_HCN),
    W2_single = list(sel = W2_selO, perm = W2_permO, step = W2_stepO, anneal = W2_annealO,
                   HC = W2_HCO),
    mseMean = list(sel = mean_mseSel, sel_E = mean_mseSelE, perm = mean_msePerm,
                   perm_E = mean_msePermE, hc=mean_mseHC, step = mean_mseStep,
                   anneal = mean_mseAnneal),
    # mseNewY = list(sel=newY_mseDSel, sel_E = newY_mseSelE, perm = newY_msePerm,
    #                perm_E = newY_msePermE, hc=newY_mseHC, step = newY_mseStep),
    mseNewMu = list(sel = newMu_mseSel, sel_E = newMu_mseSelE, perm = newMu_msePerm,
                    perm_E = newMu_msePermE, hc=newMu_mseHC, step = newMu_mseStep,
                    anneal = newMu_mseAnneal),
    # mseNewX = list(sel=newX_mseSel, sel_E = newX_mseSelE, perm = newX_msePerm,
    #                perm_E = newX_msePermE, hc=newX_mseHC, step = newX_mseStep),
    singleNewMu = list(sel = singleMu_mseSel, dist_E = singleMu_mseSelE,
                       perm = singleMu_msePerm,
                       perm_E = singleMu_msePermE, hc=singleMu_mseHC,
                       step = singleMu_mseStep, anneal=singleMu_mseAnneal),
    time = list(selTime, trueTime, stepTime, annealTime)
  )

  #add w2 data if necessary
  if( w2 ) {
    W2_w2 <- sapply(w2s$minCombPerActive, function(i)
      mse(apply(cond_mu,1,sort), apply(X[,i] %*% t(theta[,i]),1,sort)))
    mean_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(true_mu, X[,i] %*% t(theta[,i])))
    newY_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(new_Y, X[,i] %*% t(theta[,i])))
    newMu_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(new_mu, X_new[,i] %*% t(theta[,i])))
    newX_w2 <- sapply(w2s$minCombPerActive, function(i) mse_C(new_Y_new_X, X_new[,i] %*% t(theta[,i])))

    #save in output
    outList$W2_dist$w2  <- W2_w2
    outList$mseMean$w2  <- mean_w2
    outList$mseNewY$w2  <- newY_w2
    outList$mseNewMu$w2 <- newMu_w2
    outList$mseNewX$w2  <- newX_w2
  }
  # ranks <- lapply(outList, function(o) apply(sapply(o, function(cc) cc), 1, rank))
  # for(ii in seq_along(names(outList))) {
  #   ranks[[ii]][which(is.na(t(sapply(outList[[ii]], function(cc) cc))),arr.ind=TRUE)] <- NA
  # }
  # outList$rank <- lapply(ranks, rowMeans, na.rm=TRUE)
  return(outList)
}
