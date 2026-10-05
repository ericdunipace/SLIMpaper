#### Set seed for reproducibility ####
set.seed(1848078034)

#### Require Packages ####
require(SLIMpaper)
require(mosek)
require(dplyr)
require(gridExtra)

arraynum <- Sys.getenv('SLURM_ARRAY_TASK_ID')

#### Load Data ####
data(ovar, package = "SLIMpaper")

X <- ovar$train$X
time <- ovar$train$time
event <- ovar$train$event
pos <- which(apply(X, 2, function(x) all(x>0)))

XS <- ovar$test$X
XV <- ovar$val$X

# rm(ovar)

#### Set Dir ####
direct <- "Output/Ovar"
elements <- strsplit(direct, "/")[[1]][1:2]
if(!dir.exists(file.path(elements[1],elements[2]))) {
  dir.create(file.path(elements[1],elements[2]))
}
if(!dir.exists(direct)){
  dir.create(direct)
}

#### Cox Hyperparameters ####
target <- get_survival_linear_model()

#### control variables ####
n <- nrow(X)
nsamp <- 100
# nintervals <- 20
model.size <- 10
# individual.idx <- sample.int(nrow(X),1)
# print(individual.idx)

#### Load Outcome Model ####
# debugonce(target$rpost)
cox_file <- file.path(direct, "recurrence_cox.rds")
{
  # cutpoints <- survcuts <- seq(0,max(time),length.out=nintervals+1)
  # cutpoints[1] <- survcuts[1] <- 0
  recur_cox <- readRDS(cox_file)
  # cox_non.inla <- recur_cox$non.inla
  # cox_inla  <- recur_cox$inla
  theta_noint <- recur_cox$theta
  eta <- recur_cox$eta
  surv <- recur_cox$mu$S$surv
  baseSurv <- recur_cox$mu$S$base
  eta_sing <- XS %*% theta_noint
  eta_valid <- XV %*% theta_noint
  # rm(recur_cox)
}

cens_file <- file.path(direct, "recur_cens.rds")
{
  recur_cens <- readRDS(cens_file)
  cens_array <- recur_cens$mu$S$surv
}


#### Select random observation for local model ####
neighb <- SLIMpaper::rmvnorm(ncol(X) * 2, mean = c(XS), covariance = cov(X)/nrow(X))
eta_neighb <- neighb %*% theta_noint

#### Global model ####
global_file <- file.path(direct, "ovar_global_model.RDS")
if ( !file.exists(global_file) ) {
  global <- list()

  global$`B.P.` <- W2IP(X = X, Y = eta, theta = theta_noint, model.size = 1:model.size,
                        transport.method = "exact", display.progress = TRUE,
                        solution.method = "gurobi")
  global$`Relaxed B.P.` <- W2L1(X = X, Y = eta, theta = theta_noint, model.size = model.size,
                                penalty = "mcp.net", method = "selection.variable",
                                transport.method = "exact", nlambda = 1e3, alpha = 0.99, gamma = 1.5,
                                display.progress = TRUE, infimum.maxit = 2e1, maxit = 1e6)
  # global$L1    <- WPL1(X = X, Y = eta, theta = theta_noint, power = 1, model.size = model.size,
  #                      penalty = "mcp", solver = "mosek",
  #                      nlambda = 1e3, gamma = 1.5, options = list(solver_opts = list(verbose = 10)),
  #                      display.progress = TRUE, maxit = 1e6, solver_opts = list(verbose = TRUE))
  global$L2    <- WPL1(X = X, Y = eta, power = 2, theta = theta_noint, model.size = model.size,
                       family = "gaussian", penalty = "mcp.net", method = "projection",
                       nlambda = 1e3, alpha = 0.99, gamma = 1.5,
                       display.progress = TRUE, infimum.maxit = 1, maxit = 1e3)
  # global$LInf  <- WPL1(X = X, Y = eta, theta = theta_noint, power = Inf, model.size = model.size,
  #                      penalty = "mcp", solver = "mosek",
  #                      nlambda = 1e3, gamma = 1.5, options = list(solver_opts = list(verbose = 10)),
  #                      display.progress = TRUE, maxit = 1e6)

  saveRDS(global, global_file)
} else {
  global <- readRDS(global_file)
}


#### Local model ####
local_file <- file.path(direct, "ovar_local_model.RDS")
if ( !file.exists(local_file) ) {
  local <- list()

  local$`B.P.` <- W2IP(X = neighb, Y = eta_neighb, theta = theta_noint, model.size = 1:model.size,
                       transport.method = "exact", display.progress = TRUE,
                       solution.method = "gurobi")
  local$`Relaxed B.P.` <- W2L1(X = neighb, Y = eta_neighb, theta = theta_noint, model.size = model.size,
                               penalty = "mcp.net", method = "selection.variable",
                               transport.method = "exact", nlambda = 1e3, alpha = 0.99, gamma = 1.5,
                               display.progress = TRUE, infimum.maxit = 2e1, maxit = 1e6)
  # local$L1    <- WPL1(X = neighb, Y = eta_neighb, theta = theta_noint, power = 1, model.size = model.size,
  #                     penalty = "mcp", display.progress = TRUE, solver = "mosek",
  #                     transport.method = "exact", nlambda = 1e3, gamma = 1.5,
  #                     display.progress = TRUE, maxit = 1e6)
  local$L2    <- WPL1(X = neighb, Y = eta_neighb, power = 2, theta = theta_noint, model.size = model.size,
                      family = "gaussian", penalty = "mcp.net", method = "projection",
                      nlambda = 1e3, alpha = 0.99, gamma = 1.5,
                      display.progress = TRUE, infimum.maxit = 1, maxit = 1e3)
  # local$LInf  <- WPL1(X = neighb, Y = eta_neighb, theta = theta_noint, power = Inf, model.size = model.size,
  #                     penalty = "mcp", display.progress = TRUE, solver = "mosek",
  #                     nlambda = 1e3, gamma = 1.5,
  #                     display.progress = TRUE, maxit = 1e6)

  saveRDS(local, local_file)
} else {
  local <- readRDS(local_file)
}

#### add in models from the cluster ####
global_full.fn <- "Output/Ovar/ovar_global_model_full.RDS"
local_full.fn <- "Output/Ovar/ovar_local_model_full.RDS"
if(!file.exists(global_full.fn) | !file.exists(local_full.fn)) {
  # cluster.fn <- "Output/Ovar/cluster"
  # all.files <- list.files(cluster.fn)
  # gfn.string <- sapply(strsplit(all.files[grepl("global", all.files)], "_"), function(i) i[[3]])
  # gfn <- all.files[grepl("global", all.files)][order(sapply(strsplit(gfn.string, "[.]"), function(i) as.numeric(i[[1]])))]
  # lfn.string <- sapply(strsplit(all.files[grepl("local", all.files)], "_"), function(i) i[[3]])
  # lfn <- all.files[grepl("local", all.files)][order(sapply(strsplit(lfn.string, "[.]"), function(i) as.numeric(i[[1]])))]
  #
  # global.cluster.files <- lapply(gfn, function(f) readRDS(file.path(cluster.fn, f)))
  # local.cluster.files <- lapply(lfn, function(f) readRDS(file.path(cluster.fn, f)))
  #
  # combine.cluster <- function(cluster.files, X, model.size) {
  #   cluster <- cluster.files[[1]]
  #   d <- nrow(cluster[[1]]$theta[[1]])
  #   s <- ncol(cluster[[1]]$theta[[1]])
  #   zeros <- matrix(0, d,s)
  #   for(j in names(cluster)){
  #     cluster[[j]]$beta <- sapply(cluster.files, function(bb) bb[[j]]$beta)
  #     cluster[[j]]$lambda <- sapply(cluster.files, function(bb) bb[[j]]$lambda)
  #     cluster[[j]]$nonzero_beta <- sapply(cluster.files, function(bb) bb[[j]]$nonzero_beta)
  #     # cluster[[j]]$nzero <- sapply(cluster.files, function(bb) bb[[j]]$nzero)
  #     # cluster[[j]]$eta <- lapply(cluster.files, function(bb) bb[[j]]$eta)
  #     # cluster[[j]]$theta <- lapply(cluster.files, function(bb) bb[[j]]$theta)
  #     extracted <- WpProj::extractTheta(cluster[[j]], zeros)
  #     keep <- which(extracted$nzero <= model.size)
  #     cluster[[j]]$nzero <- extracted$nzero[keep]
  #     cluster[[j]]$theta <- extracted$theta[keep]
  #
  #     cluster[[j]]$eta <- lapply(cluster[[j]]$theta, function(tt) X %*%tt )
  #   }
  #
  #   return(cluster)
  # }
  #
  # global.cluster <- combine.cluster(global.cluster.files, X, model.size = model.size)
  # local.cluster  <- combine.cluster(local.cluster.files, XS, model.size = model.size)

  global.cluster <- readRDS("Output/Ovar/global_model_cluster.RDS")
  local.cluster <- readRDS("Output/Ovar/local_model_cluster.RDS")

  global <- c(global, global.cluster)
  local <- c(local, local.cluster)

  global <- global[c("B.P." ,"Relaxed B.P.","L1", "L2",  "LInf" )]
  local <- local[c("B.P." ,"Relaxed B.P.","L1", "L2",  "LInf" )]

  saveRDS(global, file = global_full.fn)
  saveRDS(local, file = local_full.fn)
} else {
  global <- readRDS(global_full.fn)
  local  <- readRDS(local_full.fn)
}



#### Calculate evaluation metrics ####
distfn <- "Output/Ovar/distances.rds"
if(!file.exists(distfn)) {
  global_dist <- distCompare(global, target = list(posterior = theta_noint, mean = eta), p = 2, ground_p = 2, method = "exact")
  # local.models.adjust <- local
  neighb_dist <- distCompare(local, target = list(posterior = theta_noint, mean = eta_neighb), p = 2, ground_p = 2, method = "exact")

  local.models.adjust <- lapply(local, function(x) {
    x$eta <- lapply(x$theta, function(tt) XS %*% tt)
    return(x)
  })
  local_dist <- distCompare(local.models.adjust, target = list(posterior = theta_noint, mean = eta_sing), p = 2, ground_p = 2, method = "exact")

  validation.models <- lapply(global, function(x) {
    x$eta <- lapply(x$theta, function(tt) XV %*% tt)
    return(x)
  })
  validation_dist <- distCompare(validation.models, target = list(posterior = theta_noint, mean = eta_valid), p = 2, ground_p = 2, method = "exact")

  saveRDS(list(global = global_dist, validation = validation_dist, neighb = neighb_dist, local = local_dist), file = distfn)
} else {
  dist <- readRDS(distfn)
  global_dist <- dist$global
  validation_dist <- dist$validation
  local_dist <- dist$local
  neighb_dist <- dist$neighb
  local.models.adjust <- lapply(local, function(x) {
    x$eta <- lapply(x$theta, function(tt) XS %*% tt)
    return(x)
  })
  validation.models <- lapply(global, function(x) {
    x$eta <- lapply(x$theta, function(tt) XV %*% tt)
    return(x)
  })
  rm(dist)
}

#### Brier score for SLIM selection ####
brierspfn <- "Output/Ovar/recur_brier_sparse.rds"
if(!file.exists(brierspfn)){
  # baseSurv <- recur_cox$mu$S$base
  nT <- nrow(baseSurv)
  nS <- ncol(baseSurv)
  temp_model <- recur_cox$model
  brier.sparse.global <- parallel::mclapply(global, function(m){
    print(class(m))
    sapply(m$theta, function(tt) {
      baseSurv_temp <-  sapply(1:ncol(tt), function(t_iter) {
        temp_model$coefficients <- tt[,t_iter]
        return(survival::survfit(temp_model)$surv)
      })
      surv_temp <- recur_cox$surv.calc(baseSurv_temp, X, tt) #simplify2array(lapply(1:n, function(i) baseSurv^matrix(exp(ee[i,]), nT, nS, byrow=TRUE)))
      brier <- target$evalfit(time, event, surv.times = sort(unique(time)),
                              surv = surv_temp,
                              cens = cens_array, method = "brier")$int.BS$intBS
      return(brier)
    })
  })
  brier.sparse.validation <- parallel::mclapply(validation.models, function(m){
    print(class(m))
    sapply(m$theta, function(tt) {
      baseSurv_temp <-  sapply(1:ncol(tt), function(t_iter) {
        temp_model$coefficients <- tt[,t_iter]
        return(survival::survfit(temp_model)$surv)
      })
      surv_temp <- recur_cox$surv.calc(baseSurv_temp, XV, tt) #simplify2array(lapply(1:n, function(i) baseSurv^matrix(exp(ee[i,]), nT, nS, byrow=TRUE)))
      brier <- target$evalfit(ovar$val$time, ovar$val$event, surv.times = sort(unique(time)),
                              surv = surv_temp,
                              cens = cens_array[,,1:nrow(XV)], method = "brier")$int.BS$intBS
      return(brier)
    })
  })
  brier.sparse.local <- parallel::mclapply(local.models.adjust, function(m){
    print(class(m))
    sapply(m$theta, function(tt) {
      baseSurv_temp <-  sapply(1:ncol(tt), function(t_iter) {
        temp_model$coefficients <- tt[,t_iter]
        return(survival::survfit(temp_model)$surv)
      })
      surv_temp <- recur_cox$surv.calc(baseSurv_temp, XS, tt) #simplify2array(lapply(1:n, function(i) baseSurv^matrix(exp(ee[i,]), nT, nS, byrow=TRUE)))
      brier <- target$evalfit(ovar$test$time, ovar$test$event, surv.times = sort(unique(time)),
                              surv = surv_temp,
                              cens = cens_array[,,1,drop=FALSE], method = "brier")$int.BS[c("intBS")]
      return(brier)
    })
  })

  saveRDS(list(global = brier.sparse.global,
               validation = brier.sparse.validation,
               local = brier.sparse.local),
          brierspfn)
  rm(nT)
  rm(nS)
} else {
  brier.sparse <- readRDS(brierspfn)
  brier.sparse.global <- brier.sparse$global
  brier.sparse.validation <- brier.sparse$validation
  brier.sparse.local <- brier.sparse$local
}

brier.global <-  target$evalfit(time, event, surv.times = sort(unique(time)), surv = surv,
                                cens = cens_array, method = "brier")
brier.validation <-  target$evalfit(time = ovar$val$time, event = ovar$val$event,
                                    surv.times = sort(unique(time)),
                                    surv = recur_cox$surv.calc(baseSurv, XV, theta_noint),
                                    cens = cens_array[,,1:nrow(XV), drop = FALSE], method = "brier")
brier.local <-  target$evalfit(time = ovar$test$time, event = ovar$test$event,
                               surv.times = sort(unique(time)),
                               surv = recur_cox$surv.calc(baseSurv, XS, theta_noint),
                               cens = cens_array[,,1, drop = FALSE], method = "brier")


# bslist <- lapply(seq_along(brier.sparse), function(i) data.frame(brier.sparse[[i]], model = names(use.models)[i], nzero = use.models[[i]]$nzero))
bslist.global <- lapply(seq_along(brier.sparse.global), function(i) {
  nzero <- global[[i]]$nzero
  # intbs <- (brier.sparse.global[[i]] - matrix(brier.global$int.BS$intBS, nrow=nsamp,
  #                                                  ncol = length(nzero)))/sd(brier.global$int.BS$intBS)
  # intbs <- (brier.sparse.global[[i]]/ matrix(brier.global$int.BS$intBS, nrow=nsamp,
  #                                             ncol = length(nzero)))
  intbs <- brier.sparse.global[[i]]
  bs.mean <- colMeans(intbs)
  # bs.mean <- apply(intbs, 2, median)

  return(data.frame(mean = bs.mean, nzero = nzero,
                    lwr = apply(intbs,2,quantile, 0.025),
                    upr = apply(intbs,2,quantile, 0.975),
                    Method = names(brier.sparse.global)[i]))
})

bsdf.global <- do.call("rbind", bslist.global)
bplot <- ggplot2::ggplot(data = bsdf.global,#[bsdf$model %in% c("Stepwise", "Selection, W2", "HC", "Projection"),],
                         ggplot2::aes(x = nzero, y = mean, fill = Method)) +
  ggsci::scale_color_jama() + ggsci::scale_fill_jama() + ggplot2::theme_bw() +
  ggplot2::xlab("Number of active coefficients") +
  ggplot2::ylab("Brier Score") +
  ggplot2::scale_x_continuous(limits = c(0,10), minor_breaks = 1:10, breaks = seq(0,10,2))

# pdf("inst/figure/applied/Ovar/global_brier.pdf", width = 5, height = 7)
# print(bplot +
#         ggplot2::geom_hline(yintercept = mean(brier.global$int.BS$intBS)) +
#         # ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
#         ggplot2::geom_line(ggplot2::aes(color = Method), size = 1.5) +
#         ggplot2::geom_point(ggplot2::aes(color = Method),  position = ggplot2::position_dodge(width=0.25),
#                             size = 3)) +
#         ggplot2::scale_y_continuous(limits = c(0.075, 0.225))
#
# dev.off()
bplot.global <- bplot +
  ggplot2::geom_hline(yintercept = mean(brier.global$int.BS$intBS)) +
  # ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
  ggplot2::geom_line(ggplot2::aes(color = Method)) +
  ggplot2::geom_point(ggplot2::aes(color = Method),  position = ggplot2::position_dodge(width=0.25)) +
  ggplot2::scale_y_continuous(limits = c(0, 0.225)) +
  ggplot2::geom_text(ggplot2::aes(x=nzero, y = mean, label = lab),
                     data=data.frame(nzero=0, mean=0.02, lab=c("Training Data"),
                                     Method = "B.P."),
                     vjust=1,
                     hjust = 0)
# bplot + ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
#   ggplot2::scale_x_continuous(limits = c(0,5.5), breaks = 0:5) + ggplot2::geom_line(ggplot2::aes(color = Method), size = 1.5) +
#   ggplot2::geom_point(ggplot2::aes(color = Method), position = ggplot2::position_dodge(width=0.25),
#                       size = 3)

bslist.validation <- lapply(seq_along(brier.sparse.validation), function(i) {
  nzero <- validation.models[[i]]$nzero
  # intbs <- (brier.sparse.validation[[i]] - matrix(brier.validation$int.BS$intBS, nrow=nsamp,
  #                                             ncol = length(nzero)))/sd(brier.validation$int.BS$intBS)
  # intbs <- (brier.sparse.validation[[i]]/matrix(brier.validation$int.BS$intBS, nrow=nsamp,
  #                                                 ncol = length(nzero)))
  intbs <- brier.sparse.validation[[i]]
  bs.mean <- colMeans(intbs)

  return(data.frame(mean = bs.mean, nzero = nzero,
                    lwr = apply(intbs,2,quantile, 0.025),
                    upr = apply(intbs,2,quantile, 0.975),
                    Method = names(brier.sparse.validation)[i]))
})

bsdf.validation <- do.call("rbind", bslist.validation)
bplot.validation <- ggplot2::ggplot(data = bsdf.validation,#[bsdf$model %in% c("Stepwise", "Selection, W2", "HC", "Projection"),],
                                    ggplot2::aes(x = nzero, y = mean, fill = Method)) +
  ggsci::scale_color_jama() + ggsci::scale_fill_jama() + ggplot2::theme_bw() +
  ggplot2::xlab("Number of active coefficients") +
  ggplot2::ylab("Brier Score") +
  ggplot2::scale_x_continuous(limits = c(0,10), minor_breaks = 1:10, breaks = seq(0,10,2))

# pdf("inst/figure/applied/Ovar/validation_brier.pdf", width = 5, height = 7)
# print(bplot.validation + #ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
#         ggplot2::geom_hline(yintercept = mean(brier.validation$int.BS$intBS)) +
#         ggplot2::geom_line(ggplot2::aes(color = Method), size = 1.5) +
#         ggplot2::geom_point(ggplot2::aes(color = Method),  position = ggplot2::position_dodge(width=0.25),
#                             size = 3)) +
#         ggplot2::scale_y_continuous(limits = c(0.075, 0.225))
# dev.off()

bplot.validation <- bplot.validation + #ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
  ggplot2::geom_hline(yintercept = mean(brier.validation$int.BS$intBS)) +
  ggplot2::geom_line(ggplot2::aes(color = Method)) +
  ggplot2::geom_point(ggplot2::aes(color = Method),  position = ggplot2::position_dodge(width=0.25)) +
  ggplot2::scale_y_continuous(limits = c(0.0, 0.225)) +
  ggplot2::geom_text(ggplot2::aes(x=nzero, y = mean, label = lab),
                     data=data.frame(nzero=0, mean=0.02, lab=c("Validation Data"),
                                     Method = "B.P."),
                     vjust=1,
                     hjust = 0)


bslist.local <- lapply(seq_along(brier.sparse.local), function(i) {
  nzero <- local[[i]]$nzero
  # intbs <- (do.call("cbind", brier.sparse.local[[i]])-matrix(brier.local$int.BS$intBS, nrow=nsamp,
  #                                                 ncol = length(nzero)))/sd(brier.local$int.BS$intBS)
  intbs <- do.call("cbind", brier.sparse.local[[i]])
  bs.mean <- colMeans(intbs)

  return(data.frame(mean = bs.mean, nzero = nzero,
                    lwr = apply(intbs,2,quantile, 0.025),
                    upr = apply(intbs,2,quantile, 0.975),
                    Method = names(brier.sparse.local)[i]))
})

bsdf.local <- do.call("rbind", bslist.local)
bplot.local <- ggplot2::ggplot(data = bsdf.local,#[bsdf$model %in% c("Stepwise", "Selection, W2", "HC", "Projection"),],
                               ggplot2::aes(x = nzero, y = mean, fill = Method)) +
  ggsci::scale_color_jama() + ggsci::scale_fill_jama() + ggplot2::theme_bw() +
  ggplot2::xlab("Number of active coefficients") +
  ggplot2::ylab("Brier Score") +
  ggplot2::scale_x_continuous(limits = c(0,10), minor_breaks = 1:10, breaks = seq(0,10,2))

# pdf("inst/figure/applied/Ovar/local_brier.pdf", width = 5, height = 7)
# print(bplot.local + #ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
#         ggplot2::geom_hline(yintercept = mean(brier.local$int.BS$intBS)) +
#         ggplot2::geom_line(ggplot2::aes(color = Method), size = 1.5) +
#         ggplot2::geom_point(ggplot2::aes(color = Method),  position = ggplot2::position_dodge(width=0.25),
#                             size = 3))
# dev.off()
bplot.local <- bplot.local + #ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = Method), alpha = .2) +
  ggplot2::geom_hline(yintercept = mean(brier.local$int.BS$intBS)) +
  ggplot2::geom_line(ggplot2::aes(color = Method)) +
  ggplot2::geom_point(ggplot2::aes(color = Method),  position = ggplot2::position_dodge(width=0.25))

global.b.var <- lapply(brier.sparse.global, FUN=apply, MARGIN = 2, sd)
validation.b.var <- lapply(brier.sparse.validation, FUN=apply, MARGIN = 2, sd)
local.b.var <- lapply(brier.sparse.local, function(bb) sapply(bb, sd))



#### Distances for Selection ####
global.dist <- distCompare(global, target = list( posterior = NULL,
                                                  mean = eta),
                           p = 2, ground_p = 2, method = c("univariate.approximation.pwr"),
                           quantity = c("mean"))
valid.dist <- distCompare(validation.models, target = list( posterior = NULL,
                                                            mean = eta_valid),
                          p = 2, ground_p = 2, method = c("univariate.approximation.pwr"),
                          quantity = c("mean"))
local.dist <- distCompare(local.models.adjust, target = list( posterior = NULL,
                                                              mean = eta_sing),
                          p = 2, ground_p = 2, method = c("exact"),
                          quantity = c("mean"))

global.dist$mean$dist <- sqrt(global.dist$mean$dist)
valid.dist$mean$dist  <- sqrt(valid.dist$mean$dist)
# local.dist$mean$dist  <- sqrt(local.dist$mean$dist) #exact is the same for single obs

global.dist$mean$groups <- forcats::fct_relevel(global.dist$mean$groups,
                                                "Relaxed B.P.", after = 1)
valid.dist$mean$groups  <- forcats::fct_relevel(valid.dist$mean$groups,
                                                "Relaxed B.P.", after = 1)
local.dist$mean$groups  <- forcats::fct_relevel(local.dist$mean$groups,
                                                "Relaxed B.P.", after = 1)

pdist.global <- plot(global.dist,  ylab = "Average 2-Wasserstein",
                     xlim = c(0,10))$mean + ggplot2::scale_y_continuous(limits = c(2,3.1)) +
  ggplot2::scale_x_continuous(breaks = seq(0,10,2), minor_breaks = 0:10,
                              limits = c(0,10), expand = c(0,0))

pdist.valid <- plot(valid.dist,  ylab = "Average 2-Wasserstein",
                    xlim = c(0,10))$mean + ggplot2::scale_y_continuous(limits = c(2,3.1)) +
  ggplot2::scale_x_continuous(breaks = seq(0,10,2), minor_breaks = 0:10,
                              limits = c(0,10), expand = c(0,0))

pdist.local <- plot(local.dist,  ylab = "Average 2-Wasserstein",
                    xlim = c(0,10))$mean + ggplot2::scale_y_continuous(limits = c(0,4.5)) +
  ggplot2::scale_x_continuous(breaks = seq(0,10,2), minor_breaks = 0:10,
                              limits = c(0,10), expand = c(0,0))

#### WpR2 ####
w2r2.global <- WPR2(Y = eta, nu = global, p = 2, method = "exact")
w2r2.valid  <- WPR2(Y = eta_valid, nu = validation.models, p = 2, method = "exact", base = colMeans(eta))
w2r2.local  <- WPR2(nu = local.dist, p = 2, method = "exact")

w2r2.global <- w2r2.global %>%  mutate(groups = factor(groups, labels = c("B.P.","Relaxed B.P.", "L1", "L2", "LInf" )))
w2r2.valid  <- w2r2.valid  %>%  mutate(groups = factor(groups, labels = c("B.P.","Relaxed B.P.", "L1", "L2", "LInf" )))
w2r2.local  <- w2r2.local  %>%  mutate(groups = factor(groups, labels = c("B.P.","Relaxed B.P.", "L1", "L2", "LInf" )))

pw2.global  <- plot(w2r2.global, xlim = c(0,10), ylim = c(0,0.35)) +
  ggplot2::scale_x_continuous(breaks = seq(0,10,2), minor_breaks = 0:10,
                              limits = c(0,10), expand = c(0,0))
pw2.valid   <- plot(w2r2.valid,  xlim = c(0,10), ylim = c(0,0.35)) +
  ggplot2::scale_x_continuous(breaks = seq(0,10,2), minor_breaks = 0:10,
                              limits = c(0,10), expand = c(0,0))
pw2.local   <- plot(w2r2.local,  xlim = c(0,10)) +
  ggplot2::scale_x_continuous(breaks = seq(0,10,2), minor_breaks = 0:10,
                              limits = c(0,10), expand = c(0,0))


#### combine plots and save ####
plot.legend <- get_legend(bplot.global + ggplot2::theme(legend.position="bottom") +
                            ggsci::scale_colour_jama(name = "Method:",
                                                     labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity])) +
                            ggsci::scale_fill_jama(name = "Method:",
                                                   labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity])))
gv.brier.plots <- arrangeGrob(bplot.global + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
                              bplot.validation + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
                              nrow = 2, ncol = 1,
                              left = grid::textGrob("Brier Score",
                                                    vjust = 2.5, rot = 90,
                                                    hjust = 0.3))
gv.dist.plots <- arrangeGrob(pdist.global + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
                             pdist.valid + ggplot2::theme(legend.position="none")+ ggplot2::xlab("")+ ggplot2::ylab(""),
                             nrow = 2, ncol = 1,
                             left = grid::textGrob("Average 2-Wasserstein",
                                                   vjust = 2.5, rot = 90,
                                                   hjust = 0.5))
gv.r2.plots <- arrangeGrob(
  pw2.global + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
  pw2.valid + ggplot2::theme(legend.position="none") + ggplot2::xlab("")+ ggplot2::ylab(""),
  nrow = 2, ncol = 1,
  left = grid::textGrob(expression(W[2]~R^2),
                        vjust = 2.15, rot = 90,
                        hjust = 0.2))
grob.global.valid <- arrangeGrob(gv.brier.plots, gv.dist.plots, gv.r2.plots,
                                 nrow=1, ncol = 3,
                                 bottom = grid::textGrob("Number of active coefficients",
                                                         vjust = -2, hjust=0.4)
)

# grob.global.valid <- arrangeGrob(bplot.global + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
#                                  pdist.global + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab("Avg. 2-Wass."),
#                                  pw2.global + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
#                                  bplot.validation + ggplot2::theme(legend.position="none") + ggplot2::xlab("") + ggplot2::ylab(""),
#                                  pdist.valid + ggplot2::theme(legend.position="none")+ ggplot2::xlab("")+ ggplot2::ylab(""),
#                                  pw2.valid + ggplot2::theme(legend.position="none") + ggplot2::xlab("")+ ggplot2::ylab(""),
#                                  nrow=2, ncol = 3,
#                                  bottom = grid::textGrob("Number of active coefficients",
#                                                          vjust = -2, hjust=0.45),
#                                  )
# middle = "Avg. 2-Wass")
pdf(file = "inst/figure/applied/Ovar/ovar_train_validation.pdf", width = 7, height = 5)
grid.arrange(grob.global.valid,
             plot.legend, nrow = 2, heights = c(10,.5))
dev.off()


local.grob <- arrangeGrob(bplot.local + ggplot2::theme(legend.position="none") + ggplot2::xlab("") ,
                          pdist.local + ggplot2::theme(legend.position="none") + ggplot2::xlab(""),
                          pw2.local + ggplot2::theme(legend.position="none") + ggplot2::xlab(""),
                          nrow=1, ncol = 3,
                          bottom = grid::textGrob("Number of active coefficients",
                                                  vjust = -2, hjust=0.4))
pdf(file = "inst/figure/applied/Ovar/ovar_local.pdf", width = 7, height = 3)
grid.arrange(local.grob,
             plot.legend, nrow = 2, heights = c(10,.5))
dev.off()

#### Get Genes ####
indexes <- lapply(global, function(m) lapply(m$theta, function(T) which(rowSums(T) != 0)))
for(i in seq_along(indexes)) {
  rmv <- rep(NA,length(indexes[[i]]))
  for(j in seq_along(indexes[[i]])){
    if(length(indexes[[i]][[j]]) > 10) rmv[j] <- j
  }
  rmv <- rmv[!is.na(rmv)]
  indexes[[i]][rmv] <- NULL
}
cnx <- colnames(X)
indx10 <- lapply(indexes, function(ii) cnx[unlist(ii[length(ii)])])
# NFX1 most commonly chosen first reduces induction of MHC-II molecules by IFN-gamma
# PIK3C2A may be involved in mitosis among many other things, may participate in EGF signaling cascade
table(unlist(indx10))

#### orders ####

ords <- ranking(validation.models$B.P., full = eta_valid, p = 2 , minCoef = 5, maxCoef = 10, quantiles = c(0,0.5,1))
rp   <- ridgePlot(validation.models, index = ords$index, minCoef = 5, maxCoef = 10, scale = 1,
                  alpha = 0.5, full = eta_valid)

ridge.grob <- arrangeGrob(rp[[1]] + ggplot2::theme_bw(11) + ggplot2::theme(legend.position = "none") + ggplot2::xlab("") +
                            ggplot2::ggtitle("Best prediction"),
                          rp[[2]] + ggplot2::theme_bw(11) + ggplot2::theme(legend.position = "none") + ggplot2::xlab("") + ggplot2::ylab("") +
                            ggplot2::ggtitle("Median prediction"),
                          rp[[3]] + ggplot2::theme_bw(11) + ggplot2::theme(legend.position = "none") + ggplot2::xlab("") + ggplot2::ylab("") +
                            ggplot2::ggtitle("Worst prediction"),
                          nrow = 1,
                          bottom = grid::textGrob("Log-Hazard Ratio",
                                                  vjust = -2, hjust=0.4))
fill.legend <- get_legend(rp[[2]] + ggplot2::theme_bw(11) +
                            ggplot2::theme(legend.position="bottom") +
                            ggplot2::scale_fill_manual(name = "Method:",
                                                       breaks= c("B.P.", "Relaxed B.P.", "L1", "L2", "LInf"),
                                                       labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity]),
                                                       values = c(ggsci::pal_jama("default")(5), "#e41a1c")))
pdf("inst/figure/applied/Ovar/valid_ridgeplot.pdf", width = 7.5, height  = 3.5)
grid.arrange(ridge.grob, nrow = 2,
             fill.legend,
             heights = c(10,1)
)
dev.off()

# sel <- sapply(theta_noint[2:5], function(x) which(rowSums(x) !=0))
# print(sel)
# colnames(recur_cox$model$X.scaled)[sel[[4]]]

