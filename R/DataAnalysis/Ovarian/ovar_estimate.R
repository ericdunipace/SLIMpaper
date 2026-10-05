#### Set seed for reproducibility ####
set.seed(601135235)
cens.seed <- sample.int(.Machine$integer.max,1)
cox.seed <- sample.int(.Machine$integer.max,1)

#### Require Packages ####
require(SLIMpaper)
require(survival)

#### Load Data ####
data(ovar, package = "SLIMpaper")

X <- ovar$train$X
time <- ovar$train$time
event <- ovar$train$event


#### Set Dir ####
direct <- "Output/Ovar"
elements <- strsplit(direct, "/")[[1]][1:2]
if(!dir.exists(file.path(elements[1],elements[2]))) {
  dir.create(file.path(elements[1],elements[2]))
}
if(!dir.exists(direct)){
  dir.create(direct)
}

#### Parallel ####
cores <- parallel::detectCores() - 1
options(mc.cores = cores)

#### Survival Functions ####
target <- get_survival_linear_model()

#### control variables ####
n <- nrow(X)
nsamp <- 100
method <- "survival.pkg-cox"

# concordance(fit)
# pred <- survfit(fit)
# class(pred) <- "survfit"
# ipred::sbrier(Surv(time = time, event = event), pred)

#### Linear Outcome Model ####
# debugonce(target$rpost)
cox_file <- file.path(direct, "recurrence_cox.rds")
cstatfn <- file.path(direct, "recur_cstat.rds")
if(!file.exists(cox_file)){
  recur_cox <- list()
  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  # varsel <- glmnet::cv.glmnet(X[,-1], cbind(time=time, status=event), family="cox",
  #                             parallel=TRUE)
  # parallel::stopCluster(cl)

  # mincvm <- which.min(varsel$cvm)
  # sel <- max( which ( varsel$cvm < (varsel$cvm[mincvm] + varsel$cvup[mincvm]) ) )
  # covar.id <- which( as.numeric(coef(varsel, varsel$lambda[sel]))!=0 )+1
  # covar.small <- which( as.numeric(coef(varsel, varsel$lambda.min))!=0 ) + 1

  # tempp <- 1:ncol(log_X)
  # tempn <- 1:n

  recur_cox <- target$rpost(n.samp = nsamp, x = X, y = time, fail = event,
                            hyperparameters = list(),
                            X.test = ovar$val$X,
                            id = 1:n, seed = cox.seed, method = method,
                            parallel = FALSE)

  theta_noint <- recur_cox$theta
  eta <- recur_cox$eta
  surv <- recur_cox$mu$S$surv
  baseSurv <- recur_cox$mu$S$base

  valid_data <- data.frame(ovar$val$X,
                           follow.up = ovar$val$time,
                           fail = ovar$val$event)
  # names(valid_data) <- c("x", "fail", "follow.up")
  # colnames(valid_data) <- c(names(recur_cox$model$coefficients), "fail","follow.up")
  cstat <- survival::concordance(recur_cox$model, newdata = valid_data)
  print(cstat$concordance)
  #0.7578797
  saveRDS(cstat, file = cstatfn)


  # original model
  print(survival::concordance(recur_cox$model)$concordance)
  #0.8892564
  # pred <- lapply(1:nrow(ovar$val$X), function(i) {
  #   x <- data.frame(ovar$val$X[i,,drop=FALSE])
  #   colnames(x) <- names(recur_cox$model$coefficients)
  #   survival::survfit(recur_cox$model, newdata = x)
  #   })
  # for(i in seq_along(pred)) class(pred[[i]]) <- "survfit"
  # ipred::sbrier(survival::Surv(time = ovar$val$time, event = ovar$val$event), pred)
  # # 0.2207157
  # mean.val.pred <- survival::survfit(recur_cox$model, newdata = valid_data)
  # class(mean.val.pred) <- "survfit"
  # ipred::sbrier(survival::Surv(time = ovar$val$time, event = ovar$val$event), mean.val.pred)
  # # 0.1967461
  saveRDS(recur_cox, file = cox_file)
  # rm(recur_cox)
} else {
  recur_cox <- readRDS(cox_file)
  # cox_non.inla <- recur_cox$non.inla
  # cox_inla  <- recur_cox$inla
  theta_noint <- recur_cox$theta
  eta <- recur_cox$eta
  surv <- recur_cox$mu$S$surv
  baseSurv <- recur_cox$mu$S$base
  # rm(recur_cox)
}

#### censoring model ####
cens_file <- file.path(direct, "recur_cens.rds")
brierfn <- file.path(direct, "recur_brier.rds")
if(!file.exists(cens_file)) {
  # recur_cens <- list(non.inla = NULL, inla = NULL, bart = NULL)
  cens <- (1-event)
  # if(method != "inla") {
  #   recur_cens$non.inla <- target$rpost(n.samp = nsamp, x = NULL,
  #                                       y = time, fail = cens,
  #                                       hyperparameters = hyperparameters,
  #                                       id = 1:n, nchain = nchain, jags_dir = jags_dir,
  #                                       d
  # }
  # cutpoints <- quantile(time, seq(0,1,length.out=nintervals))
  # cutpoints[1] <- 0
  # censpoints <- c(0, unique(sort(time)))
  # recur_cens$inla  <- target$rpost(n.samp = nsamp, x = recur_cox$log_lX,
  #                                  y = time, fail = cens,
  #                                  hyperparameters = hyperparameters,
  #                                  id = 1:n, nchain = nchain, jags_dir = jags_dir,
  #                                  method = "inla", seed = sample.int(.Machine$integer.max,1),
  #                                  parallel = FALSE, thin = thin,
  #                                  cutpoints = cutpoints)#sort(unique(time)))
  recur_cens  <- target$rpost(n.samp = nsamp, x = NULL,
                              y = time, fail = cens,
                              hyperparameters = list(),
                              id = 1:n,
                              method = "survival.pkg-km", seed = cens.seed)
  cens_array <- recur_cens$mu$S$surv


  # mean.val.pred <- survival::survfit(recur_cox$model)
  # class(mean.val.pred) <- "survfit"
  # ipred::sbrier(recur_cox$model$y, mean.val.pred)
  #0.1353387
  # debugonce(target$evalfit)
  brier <-  target$evalfit(time, event, surv.times = sort(unique(time)), surv = surv,
                           cens = cens_array, method = "brier")
  # debugonce(brier.score)
  print(brier$int.BS$mean)
  #0.03453203
  train_data <- data.frame(ovar$train$X,
                           follow.up = ovar$train$time,
                           fail = ovar$train$event)
  valid_data <- data.frame(ovar$val$X,
                           follow.up = ovar$val$time,
                           fail = ovar$val$event)

  # bs_calc_train <- pec(recur_cox$model, formula=Surv(follow.up,fail)~1, data = train_data)
  # bs_train <- crps(object = bs_calc_train, times = bs_calc_train$time,
  #                start = bs_calc_train$start,
  #                what = "AppErr")
  # print(bs_train[,ncol(bs_train)])
  # 0.04019197

  bs_calc_val <- pec(recur_cox$model, formula=Surv(follow.up,fail)~1, data = valid_data)
  bs_val <- crps(object = bs_calc_val, times = bs_calc_val$time,
                 start = bs_calc_val$start,
                 what = "AppErr")
  print(bs_val[,ncol(bs_val)])
  #0.093

  #brier score of samples

  brier <-  target$evalfit(ovar$val$time, ovar$val$event,
                           surv.times = sort(unique(time)),
                           surv = recur_cox$test$mu$S,
                           cens = cens_array[,,1:nrow(ovar$val$X)], method = "brier")
  # debugonce(brier.score)
  print(brier$int.BS$mean)
  # 0.06202829

  #brier score of mle
  mle.surv.val <- recur_cox$surv.calc(matrix(survival::survfit(recur_cox$model)$surv), ovar$val$X, matrix(recur_cox$model$coefficients))
  mle.cens.val <- array(recur_cens$model$surv, dim = dim(mle.surv.val))
  brier <-  target$evalfit(ovar$val$time, ovar$val$event,
                           surv.times = sort(unique(time)),
                           surv = mle.surv.val,
                           cens = mle.cens.val, method = "brier")
  # debugonce(brier.score)
  print(brier$int.BS$mean)
  # 0.1388079


  # saveRDS(brier, file = brierfn)
  # recur_cens$bart  <- survBart(x.train = scale(recur_cox$log_sX),
  #                              y.train = NULL,
  #                              times = time,
  #                              delta = cens, printevery = 10,
  #                              ndpost = nsamp,
  #                              nskip = nsamp * thin,
  #                              keepevery = thin)
  saveRDS(recur_cens, file = cens_file)
  # cens_inla <- recur_cens$inla
  # rm(recur_cens)

  # for(i in 1:n) cens_array[,,i] <- cens
  rm(recur_cens)
} else {
  recur_cens <- readRDS(cens_file)
  # cens_inla <- recur_cens$inla
  cens_array <- recur_cens$mu$S$surv
  censpoints <- c(0, unique(sort(time)))
  # cens_array <- array(NA, dim=c(nrow(cens), nsamp,n))
  # for(i in 1:n) cens_array[,,i] <- cens
  rm(recur_cens)
  # rm(cens_inla)
}

