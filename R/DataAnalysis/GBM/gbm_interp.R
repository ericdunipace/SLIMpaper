rm(list=ls())
set.seed(-903881470) #from random.org

arraynum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


#### Load Packages ####
require(SLIMpaper)
require(Matrix)
require(gamm4)
# require(qgam)
require(mgcv)
require(BART)
require(ggplot2)
require(dplyr)
require(gridExtra)
require(xtable)


#### Load Data ####
data.file   <- "../Data/GBM/GBM_DF.RData"
out.dir     <- "Output/GBM"

load(data.file)

n <- K <- nT <- nrow(DF_train$X)
n_valid <- nrow(DF_valid$X)
n_test <- nrow(DF_test$X)
n_neighb <- nrow(DF_test$neighb)

mm.df <- model.matrix(~ scale(Age) + Gender + KPS.C + MGMT + Resection+ 0, data= DF_train$X)
mm.df_test <- model.matrix(~ scale(Age, center = mean(DF_train$X$Age),
                                   scale = sd(DF_train$X$Age)) + Gender + KPS.C + MGMT + Resection+ 0,
                           data= rbind(DF_test$X, DF_valid$X))
colnames(mm.df_test) <- colnames(mm.df)

mm.df_valid <- mm.df_test[-1,,drop=FALSE]
mm.df_test <- mm.df_test[1,,drop=FALSE]

mm.df_neighb <- DF_test$neighb

preDF <- surv.pre.bart(x.train = mm.df,
                       x.test = mm.df_neighb,
                       times = DF_train$Y$T.OS, delta = DF_train$Y$C.OS)
df_train <- surv.pre.bart(x.train = mm.df,
                          x.test = mm.df,
                          times = DF_train$Y$T.OS, delta = DF_train$Y$C.OS)$tx.test[,-1]
df_valid <- surv.pre.bart(x.train = mm.df,
                          x.test = mm.df_valid,
                          times = DF_train$Y$T.OS, delta = DF_train$Y$C.OS)$tx.test[,-1]
df_test  <- surv.pre.bart(x.train = mm.df,
                          x.test = mm.df_test,
                          times = DF_train$Y$T.OS, delta = DF_train$Y$C.OS)$tx.test[,-1]

#age is scaled
#gender 1 = Female, 0 = male
#KPS 1 if KPS >= 90
#MGMT 1 if any methylation
#Resection biopsy, sub, gross

#### can load BART model, but not neccessary ####
# df.bart     <- readRDS(file = file.path(out.dir, "gbm_model.rds"))
# df.cens.bart<- readRDS(file = file.path(out.dir, "gbm_cens.rds"))
# predDF      <- readRDS(file = file.path(out.dir, "gbm_model_pred.rds"))
# censDF      <- readRDS(file = file.path(out.dir, "gbm_cens_pred.rds"))

#### set-up interp data ####
interpY_full     <- readRDS(file = file.path(out.dir, "gbm_model_Y.rds"))
interpY_neighb   <- interpY_full[seq(1,by=1, length.out = K * n_neighb),]
interpY_train    <- readRDS(file = file.path(out.dir, "gbm_model_Y_training.rds"))
interpY_valid    <- interpY_full[-seq(1,by=1, length.out = K * (n_neighb + 1)),]
interpY_test     <- interpY_full[((K * n_neighb)+1):((K * n_neighb)+K),, drop = FALSE]

# interpX     <- construct_interp_survival_matrix(mm.df_neighb, K, intercept=FALSE)
nS          <- ncol(interpY_neighb)
nT          <- preDF$K
dat_neighb  <- gamm_interp_data_gbm(preDF$tx.test[,-1], times = preDF$times)
dat_valid   <- gamm_interp_data_gbm(df_valid, times = preDF$times)
dat_train   <- gamm_interp_data_gbm(df_train, times = preDF$times)
dat_test   <- gamm_interp_data_gbm(df_test, times = preDF$times)



#### Interp hyperparam ####
# gamma       <- 3
# alpha       <- 0.5
# nlambda     <- 100
# lambda.min.ratio <- 1e-10
# p           <- 2
# groups      <- rep(c(1:5,5,5), nT)
# # edges       <- c(sapply(0:(d*(nT-1)-1), function(i) c(1,8) + i))
# # gr          <- graph(edges, directed = FALSE)
# sparseX     <- Matrix::Matrix(interpX)
# qrX         <- Matrix::qr(sparseX)
# # xtx         <- Matrix::crossprod(sparseX)
# # xty         <- Matrix::crossprod(sparseX, interpY)
# new_times   <- runif(15)
# std_times   <- c(preDF$times/max(preDF$times), new_times)
# cost        <- cost_calc(t(std_times), t(std_times), 2.0)^2
# x_smooth    <- cbind(1, matrix(preDF$times))
# f_groups    <- rep(groups,nS)
# pf          <- rep(1, nT * ncol(sparseX))
# lambda_f    <- 1
# lambda_g    <- group_lambda_zero(sparseX, interpY, groups, sqrt(tapply(rep(1,length(f_groups)), f_groups, sum)))
#

#### Formula ####
gam_form  <- formula(y ~ Resection +
                       # Age + Gender + KPS +  MGMT +#Age:Resection +
                       s(time, by = Age) +
                       s(time, by = Gender) +
                       s(time, by = KPS) +
                       s(time, by = MGMT) +
                       s(time, by = Resection, pc = 0) - 1)

#### Interpretable Models ####
df_imp_file <- file.path(out.dir, "varImportance_df.rds")

if( !file.exists(df_imp_file) ) {
  # df_grp    <- W2L1(X=interpX, Y=interpY, theta = NULL,
  #                   groups = groups,
  #                   penalty="grp.mcp.net", alpha=alpha, gamma = gamma,
  #                   nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
  #                   infimum.maxit=1, maxit = 1e3,
  #                   display.progress = TRUE,
  #                   penalty.factor = pf,
  #                   lambda = exp(log(lambda_g) + seq(0, log(lambda.min.ratio), length.out = nlambda)),
  #                   method="projection", tol = 1e-3)
  # gamm_df   <- data.frame(y = interpY[,1], gammXmod, time = gammT)
  # gamm_df   <- data.frame(y = plogis(interpY[,1]), gammXmod, time = gammT)
  # gamm_mm  <- model.matrix(~Age + Gender + KPS + MGMT + Resection + 1, data = gamm_df)[,-1]

  # gam_form  <- formula(y ~ Resection + Age:Resection +
  #                        s(time, by = Age) +
  #                        s(time, by = I(Age^2)) +
  #                        s(time, by = I(Age*Gender)) +
  #                        s(time, by = I(Age*KPS)) +
  #                        s(time, by = I(Age*MGMT)) +
  #                        s(time, by = I(Gender*MGMT)) +
  #                        s(time, by = I(Gender*KPS)) +
  #                        s(time, by = I(KPS*MGMT)) +
  #                        s(time, by = Age.ResectionSub) +
  #                        s(time, by = Age.ResectionGross) +
  #                        s(time, by = Gender) +
  #                        s(time, by = KPS) +
  #                        s(time, by = MGMT) +
  #                        s(time, by = Resection) - 1)
  # gam_form  <- formula(y ~ Resection +
  #                        Age + Gender + KPS +  MGMT +#Age:Resection +
  #                        s(time, by = Age, pc = 0) +
  #                        s(time, by = Gender, pc = 0) +
  #                        s(time, by = KPS, pc = 0) +
  #                        s(time, by = MGMT, pc = 0) +
  #                        s(time, by = Resection, pc = 0) - 1)
  # gaml_form  <- formula(y ~ X +#Age:Resection +
  #                        s(time, by = Age, pc = 0) +
  #                        s(time, by = Gender, pc = 0) +
  #                        s(time, by = KPS, pc = 0) +
  #                        s(time, by = MGMT, pc = 0) +
  #                        s(time, by = Resection, pc = 0) - 1)
  # gamm_gl  <- gamm_df
  # gamm_gl$X <- gamm_mm
  # df_gam   <- gamlasso(formula = gaml_form,
  #                      data = gamm_gl,
  #                      # linear.terms = c("Age","Gender",
  #                      #                  "KPS","MGMT",
  #                      #                  "ResectionBiopsy",
  #                      #                  "ResectionSub", "ResectionGross"),
  #                      linear.penalty = "l1",
  #                      smooth.penalty = "l1")
  # df_gam   <- gam_iterate(gam_form, y = interpY_train,
  #                         x = dat_train$gammX, extract = dat_train$extract_terms,
  #                         time = dat_train$times, nT=nT,
  #                         which.gam = "gam")
  # df_gam_L1<- gam_iterate(gam_form, y = interpY, x = gammX, extract = extract_terms, time = gammT, nT=nT,
  #                         which.gam = "qgam", lsig = -1.606527, err = 0.15)

  # df_ols   <- as.matrix(Matrix::qr.coef(qrX, interpY))
  # dat <- gp_data(df_ols, d, nT, preDF$times)
  # stan_dat    <- list(N = nS*nT,
  #                     nT = nT,
  #                     N_total = nT + 15,
  #                     P = d,
  #                     y = dat$y,
  #                     cost = cost,
  #                     time_idx = dat$time
  # )
  # df_smooth_fit <- lapply(1:d, function(j) lm.fit(x = x_smooth, y = df_ols[seq(j,d*nT,d),])$coef)
  # df_smooth_cut <- lapply(seq_along(df_smooth_fit), function(dd) {
  #   pred <- x_smooth %*% df_smooth_fit[[dd]]
  #   return(data.frame(mean = rowMeans(pred), lwr = apply(pred,1,quantile, 0.025),
  #               upr = apply(pred,1,quantile, 0.975), Coefficient = dd, Time = preDF$times))
  #   })
  # df_smooth_df  <- do.call("rbind", df_smooth_cut)
  # df_smooth     <- list(data = df_smooth_df, fit = df_smooth_fit)
  # gp_mod   <- rstan::stan_model(file = "exec/Stan/normal_gp.stan")
  # df_gp    <- rstan::sampling(gp_mod, iter = 2000, pars = c("sds","L_corr","pred_eta"), data = stan_dat)
  # df_fuse <- lapply(1:nS, function(i) fusedlasso(y = interpY[,i], X = interpX, gamma = 0, graph = gr))

  #combine from cluster
  all.files <- list.files("Output/GBM/cluster")
  gfn.string <- sapply(strsplit(all.files[grepl("global", all.files)], "_"), function(i) i[[3]])
  gfn <- all.files[grepl("global", all.files)][order(sapply(strsplit(gfn.string, "[.]"), function(i) as.numeric(i[[1]])))]
  lfn.string <- sapply(strsplit(all.files[grepl("local", all.files)], "_"), function(i) i[[3]])
  lfn <- all.files[grepl("local", all.files)][order(sapply(strsplit(lfn.string, "[.]"), function(i) as.numeric(i[[1]])))]

  global.cluster.files <- lapply(gfn, function(f) readRDS(file.path("Output/GBM/cluster", f)))
  local.cluster.files <- lapply(lfn, function(f) readRDS(file.path("Output/GBM/cluster", f)))

  combine.cluster <- function(cluster.files) {
    n <- length(cluster.files)
    x <- array(data = NA, dim = c(dim(cluster.files[[1]][[1]])[1:2], n))
    for(i in 1:n) {
      x[,,i] <- cluster.files[[i]][[1]]
    }

    return(x)
  }

  global.cluster <- combine.cluster(global.cluster.files)
  local.cluster  <- combine.cluster(local.cluster.files)
  dimnames(global.cluster) <- list(time = preDF$times,
                                   coefficient = c("Age", "Gender", "KPS ≥ 90", "Methylation", "Resection: Biopsy",
                                                   "Resection: Sub Total", "Resection: Gross Total"),
                                   iteration = 1:length(gfn.string)
  )
    dimnames(local.cluster) <- list(time = preDF$times,
                                   coefficient = c("Age", "Gender", "KPS ≥ 90", "Methylation", "Resection: Biopsy",
                                                   "Resection: Sub Total", "Resection: Gross Total"),
                                   iteration = 1:length(lfn.string)
                                   )
  df_imp        <- list(global = list(
                                # glasso = df_grp,
                                gam    = global.cluster
                                # , gam_l1 = df_gam_L1
                                # ols    = df_ols,
                                # smooth = df_smooth
                                ),
                        local = list(gam = local.cluster))
  saveRDS(df_imp, df_imp_file)
  local_gam <- df_imp$local$gam
  global_gam <- df_imp$global$gam
} else {
  df_imp    <- readRDS(df_imp_file)
  # df_smooth_df  <- df_imp$smooth$data
  # df_gam <- df_imp$gam
  local_gam <- df_imp$local$gam
  global_gam <- df_imp$global$gam
  # df_gam_L1 <- df_imp$gam_l1
}

#### calculate distances ####

# W2 bw the interp and full
if(!file.exists(file.path(out.dir, "wass_calc.rds"))) {
  coef.to.pred <- function(coef, formula, x) {
    nS <- dim(coef)[3]
    d  <- dim(coef)[2]
    nT <- dim(coef)[1]
    rep.time <- nrow(x)/nT
    rep.row  <- rep(1:nT, rep.time)

    # gam.no.fit <- gam(formula, data = data, fit = FALSE)
    x.pred <- model.matrix(~ . + 0, data = x) #predict(gam.no.fit, type = "lpmatrix")

    pred <- matrix(NA_real_, nrow= nrow(x), ncol = nS)
    for(i in 1:nS) {
      pred[,i] <- rowSums(x.pred * coef[rep.row,,i])
    }
    return(pred)
  }
  base.calc <- function(eta, nT, intercept.mat= NULL) {
    if(is.null(intercept.mat)) intercept.mat <- rep(1, nrow(eta))
    base <- matrix(NA_real_, nrow = nrow(eta), ncol = ncol(eta))
    idx  <- rep(1:nT, nrow(eta)/nT)
    df <- data.frame(idx = factor(idx), int = intercept.mat, eta = eta)
    x  <- Matrix::sparse.model.matrix(~ int * idx, data = df)
    lm.fit <- sapply(1:ncol(eta), function(e) MatrixModels:::lm.fit.sparse(x = x, y = c(eta[,e])))
    base <- x %*% lm.fit
    # means <- sapply(1:nT, function(i)
    #                 tapply(eta[idx==i,], intercept.mat[idx==i], colMeans))
    # for(i in 1:nT){
    #   base[idx == i,] <- matrix(means[,i], nrow = nT, ncol = ncol(eta), byrow = TRUE)
    # }
    return(list(base=as.matrix(base), coef = lm.fit))
  }
  # train_gam_data <- data.frame(y = interpY_train[,1], dat_train$gammX, time = preDF$times)
  train_pred <- coef.to.pred(global_gam, formula = gam_form, x = dat_train$gammX)

  base_calc <- base.calc(interpY_train[,1:ncol(train_pred)],
                          nT = nT, intercept.mat = dat_train$gammX$Resection)

  base_train <- base_calc$base

  wpr2.global <-
    WPR2(Y = plogis(interpY_train)[,1:ncol(train_pred)],
       nu =plogis(train_pred),
       p = 2,
       method = "exact",
       base = plogis(base_train),
       niter = 1e6)

  valid_pred <- coef.to.pred(global_gam, formula = gam_form, x = dat_valid$gammX)

  vdf <- data.frame(int = dat_valid$gammX$Resection,
                    idx = factor(rep(1:nT, nrow(interpY_valid)/nT)))
  valid.base <- as.matrix(sparse.model.matrix(~ int * idx, data = vdf) %*% base_calc$coef)
  wpr2.valid <-
    WPR2(Y = plogis(interpY_valid)[,1:ncol(valid_pred)],
         nu =plogis(valid_pred),
         p = 2,
         method = "exact",
         base = plogis(valid.base[,1:ncol(valid_pred)]),
         niter = 1e6)

  local_pred <- coef.to.pred(local_gam, formula = gam_form, x = dat_test$gammX)
  ldf <- data.frame(int = dat_test$gammX$Resection,
                    idx = factor(rep(1:nT, nrow(interpY_test)/nT)))
  local.base <- as.matrix(sparse.model.matrix(~ int * idx, data = ldf) %*% base_calc$coef)
  wpr2.local <-
    WPR2(Y = plogis(interpY_test)[,1:ncol(local_pred)],
         nu =plogis(local_pred),
         p = 2,
         method = "exact",
         base = plogis(local.base[,1:ncol(local_pred)]),
         niter = 1e6)
  wass_met <- list(global = wpr2.global,
                   validation = wpr2.valid,
                   local = wpr2.valid)

  saveRDS(wass_met, file.path(out.dir, "wass_calc.rds"))
} else {
  wass_met <- readRDS(file.path(out.dir, "wass_calc.rds"))
  wpr2.global <- wass_met$global
  wpr2.valid <- wass_met$validation
  wpr2.local <- wass_met$local
}

wpr2.mat <- matrix(sapply(wass_met, function(w) w$r2), nrow = 1)
colnames(wpr2.mat) <- c("Training data", "Validation data",
                        "Test data")
rownames(wpr2.mat) <- c("$W_2 R^2$")
xwpr2 <- xtable(wpr2.mat,
                caption = "$W_2 R^2$ statistics for the the training, validation, and test data.
             The null models in each case are an intercept only model using the resection status
             over time. The models are the same in the training and validation cases meaning that
             the validation data evaluates an out of sample fit. The test data model was
             fit on an interpretable neighborhood around the test point and then evaluated
             at the test data.")
print(xwpr2,
      sanitize.text.function = function(x) {x},
      sanitize.colnames.function = function(x){x},
      sanitize.rownames.function = function(x){x},
      file = "inst/table/applied/GBM/gbm_w2r2.tex"
      )

#### Plots ####
#https://stackoverflow.com/questions/11979017/changing-facet-label-to-math-formula-in-ggplot2
# facet_wrap_labeller <- function(gg.plot,labels=NULL) {
#   #works with R 3.0.1 and ggplot2 0.9.3.1
#   require(gridExtra)
#
#   g <- ggplot2::ggplotGrob(gg.plot)
#   gg <- g$grobs
#   strips <- grep("strip_t", names(gg))
#
#   for(ii in seq_along(labels))  {
#     modgrob <- grid::getGrob(gg[[strips[ii]]], "strip.text",
#                        grep=TRUE, global=TRUE)
#     gg[[strips[ii]]]$children[[modgrob$name]] <- grid::editGrob(modgrob,label=labels[ii])
#   }
#
#   g$grobs <- gg
#   class(g) = c("arrange", "ggplot",class(g))
#   g
# }

trajectories.plot <- function(gam_mod1, gam_mod2, plot.fn, facet.labels) {
  data.prep <- function(gam_mod) {
    dft <- as.data.frame.table(gam_mod)
    colnames(dft) <- c("time", "variable", "iter", "value")
    levels(dft$variable) <- c("Age", "Gender", "KPS >= 90", "Methylation", "Resection: Biopsy",
                              "Resection: Sub-total", "Resection: Total")
    dft$time <- as.numeric(as.character(dft$time))
    dft$iter <- as.numeric(dft$iter)
    dft$or <- exp(dft$value)

    mean_dft <- as.data.frame.table(apply(gam_mod, 1:2, mean))
    colnames(mean_dft) <- c("time", "variable", "lo")

    mean_dft$time <- as.numeric(as.character(mean_dft$time))
    levels(mean_dft$variable) <- c("Age", "Gender", "KPS >= 90", "Methylation", "Resection: Biopsy",
                                   "Resection: Sub-total", "Resection: Total")
    mean_dft$or <- as.data.frame.table(apply(exp(gam_mod), 1:2, mean))$Freq
    mean_dft$iter <- factor(rep(1, nrow(mean_dft)))
    return(list(full = dft, mean = mean_dft))
  }

  gm1 <- data.prep(gam_mod1)
  gm2 <- data.prep(gam_mod2)
  dft <- rbind(gm1$full %>% mutate(sample = "global"),
          gm2$full %>% mutate(sample = "local"))
  mean_dft <- rbind(gm1$mean %>% mutate(sample = "global"),
                           gm2$mean %>% mutate(sample = "local"))
  dft$sample <- factor(dft$sample)
  mean_dft$sample <- factor(mean_dft$sample)

  small_coef <- c("Age", "Gender", "KPS >= 90", "Methylation")
  p <- dft %>%
    filter(iter %in% as.integer(seq(1,100, length.out = 100))) %>%
    filter(variable %in% small_coef) %>%
    ggplot2::ggplot(aes(x = time, y = or, group = iter)) +
    ggplot2::geom_line(color = "grey", alpha = 0.1) +
    ggplot2::geom_line(data = mean_dft[mean_dft$variable %in% small_coef,], aes(x = time, y = or), color="blue") +
    ggplot2::facet_grid(sample ~ variable,
                        labeller = label_parsed) +
    ggplot2::theme_bw() +
    ggplot2::geom_abline(slope = 0, intercept = 1, color = "red", linetype = "dotted" ) +
    ggplot2::ylab("Coefficient value (Odds Ratio Scale)") +
    ggplot2::xlab("Time (days)") + ggplot2::coord_cartesian(ylim=c(0, 7.5), expand = FALSE) +
    ggplot2::expand_limits(y = c(0,NA)) +
    ggplot2::theme(
      # panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::theme(
                  strip.text.y = ggplot2::element_blank(),
                   # strip.text.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank()) +
    ggplot2::geom_text(aes(time, or, label=lab),
                       data=data.frame(time=200, or=7, lab=c("Training Data","","","",
                                                        "Test Data","","",""),
                                       variable=rep(small_coef,2),
                                       iter = "1",
                                       sample=rep(c("global", "local"), each = 4)),
                                       vjust=1,
                                       hjust = 0)
  pdf(plot.fn, height = 5, width = 7.5)
  print(p)
  dev.off()
}


trajectories.plot(global_gam, local_gam,
                  plot.fn = "inst/figure/applied/GBM/oddsratio_gbm.pdf",
                  facet.labels = c("Age", "Gender", expression(KPS >= 90),
                                   "Methylation"))

# df_gam_df <- as.data.frame.table(df_gam)
# colnames(df_gam_df) <- c("time", "variable", "iter", "value")
# levels(df_gam_df$time) <- sort(unique(preDF$times))
# df_gam_df$time <- as.numeric(as.character(df_gam_df$time))
# levels(df_gam_df$variable) <- c("Age (S.D.)", "Gender", "KPS >= 90", "Methylation", "Resection: Biopsy",
#                                 "Resection: Sub Total", "Resection: Gross Total")
# df_gam_df$iter <- as.numeric(df_gam_df$iter)
# df_gam_df$or <- exp(df_gam_df$value)

# df_gam_mean <- as.data.frame.table(apply(df_gam, 1:2, mean))
# colnames(df_gam_mean) <- c("time", "variable", "lo")
# levels(df_gam_mean$time) <- sort(unique(preDF$times))
# df_gam_mean$time <- as.numeric(as.character(preDF$times))
# levels(df_gam_mean$variable) <- c("Age (S.D.)", "Gender", "KPS >= 90", "Methylation", "Resection: Biopsy",
#                           "Resection: Sub Total", "Resection: Gross Total")
# df_gam_mean$or <- as.data.frame.table(apply(exp(df_gam), 1:2, mean))$Freq
# df_gam_mean$iter <- factor(rep(1, nrow(df_gam_mean)))
#
# # plot(x = preDF$times, y = Matrix::rowMeans(df_ols)[seq(4, d*nT, by = d)])
# # plot(x = preDF$times, y = rowMeans(df_gam[,4,]))
# # plot(x = preDF$times, y = Matrix::rowMeans(df_imp$glasso$theta[[2]])[seq(4, d*nT, by = d)], type = "l")
# # plot(x = rep(preDF$times,nS), y = c(x_smooth %*% df_smooth_fit[[1]]), pch = ".")
# # ggplot2::ggplot(df_smooth_df, aes(x = Time, y = mean)) +
# #   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) +
# #   geom_line() +
# #   facet_wrap(Coefficient ~.)
#
# pdf("inst/figure/applied/GBM/logodds_gbm.pdf", height = 5, width = 7.5)
# df_gam_df %>%
#   filter(iter %in% as.integer(seq(1,2000, length.out = 100))) %>%
#   ggplot2::ggplot(aes(x = time, y = value, group = iter)) +
#   ggplot2::geom_line(color = "grey", alpha = 0.1) +
#   ggplot2::geom_line(data = df_gam_mean, aes(x = time, y = lo), color="blue") +
#   ggplot2::facet_wrap(variable ~.) + ggplot2::theme_bw() +
#   ggplot2::geom_abline(slope = 0, intercept = 0, color = "red", linetype = "dotted" ) +
#   ggplot2::ylab("Coefficient value (Log-Odds Scale)") +
#   ggplot2::xlab("Time (Days)") +
#   ggplot2::theme(
#     # panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank())
# dev.off()
#
# pdf("inst/figure/applied/GBM/oddsratio_gbm.pdf", height = 5, width = 7.5)
# small_coef <- c("Age (S.D.)", "Gender", "KPS >= 90", "Methylation")
# df_gam_df %>%
#   filter(iter %in% as.integer(seq(1,2000, length.out = 100))) %>%
#   # filter(variable %in% small_coef) %>%
#   ggplot2::ggplot(aes(x = time, y = or, group = iter)) +
#   ggplot2::geom_line(color = "grey", alpha = 0.1) +
#   ggplot2::geom_line(data = df_gam_mean[df_gam_mean$variable %in% small_coef,], aes(x = time, y = or), color="blue") +
#   ggplot2::facet_wrap(variable ~.) + ggplot2::theme_bw() +
#   ggplot2::geom_abline(slope = 0, intercept = 1, color = "red", linetype = "dotted" ) +
#   ggplot2::ylab("Coefficient value (Odds Ratio Scale)") +
#   ggplot2::xlab("Time (Days)") + ggplot2::coord_cartesian(ylim=c(0, 10.5), expand = FALSE) +
#   ggplot2::expand_limits(y = c(0,NA)) +
#   ggplot2::theme(
#     # panel.grid.major = ggplot2::element_blank(),
#                  panel.grid.minor = ggplot2::element_blank())
# dev.off()
#
#
# pdf("inst/figure/applied/GBM/logodds_gbm.pdf", height = 5, width = 7.5)
# df_gam_df %>%
#   filter(iter %in% as.integer(seq(1,2000, length.out = 100))) %>%
#   ggplot2::ggplot(aes(x = time, y = value, group = iter)) +
#   ggplot2::geom_line(color = "grey", alpha = 0.1) +
#   ggplot2::geom_line(data = df_gam_mean, aes(x = time, y = lo), color="blue") +
#   ggplot2::facet_wrap(variable ~.) + ggplot2::theme_bw() +
#   ggplot2::geom_abline(slope = 0, intercept = 0, color = "red", linetype = "dotted" ) +
#   ggplot2::ylab("Coefficient value (Log-Odds Scale)") +
#   ggplot2::xlab("Time (Days)") +
#   ggplot2::theme(
#     # panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank())
# dev.off()
#
# pdf("inst/figure/applied/GBM/oddsratio_gbm.pdf", height = 5, width = 7.5)
# small_coef <- c("Age (S.D.)", "Gender", "KPS >= 90", "Methylation")
# df_gam_df %>%
#   filter(iter %in% as.integer(seq(1,2000, length.out = 100))) %>%
#   # filter(variable %in% small_coef) %>%
#   ggplot2::ggplot(aes(x = time, y = or, group = iter)) +
#   ggplot2::geom_line(color = "grey", alpha = 0.1) +
#   ggplot2::geom_line(data = df_gam_mean[df_gam_mean$variable %in% small_coef,], aes(x = time, y = or), color="blue") +
#   ggplot2::facet_wrap(variable ~.) + ggplot2::theme_bw() +
#   ggplot2::geom_abline(slope = 0, intercept = 1, color = "red", linetype = "dotted" ) +
#   ggplot2::ylab("Coefficient value (Odds Ratio Scale)") +
#   ggplot2::xlab("Time (Days)") + ggplot2::coord_cartesian(ylim=c(0, 10.5), expand = FALSE) +
#   ggplot2::expand_limits(y = c(0,NA)) +
#   ggplot2::theme(
#     # panel.grid.major = ggplot2::element_blank(),
#     panel.grid.minor = ggplot2::element_blank())
# dev.off()

# rm(list=ls())
# set.seed(-903881470) #from random.org
#
# #### Load Packages ####
# require(SLIMpaper)
# require(BART)
# data.file   <- "../Data/GBM/GBM_DF.RData"
# out.dir     <- "Output/GBM"
#
# load(data.file)
#
# mm.df       <- model.matrix(~  scale(Age) + Gender + KPS.C + MGMT + Resection+ 0, data= DF.X)
# rm(DF.X)
# rm(DF.Y)
# bart.file <- file.path(out.dir, "gbm_model.rds")
# df.bart <- readRDS(file = bart.file)
#
# cfDF <- mm.df
# cfDF[,"Gender"] <- 1
# cfDF2 <- mm.df
# cfDF2[,"Gender"] <- 0
# rm(mm.df)
#
# test.pred0  <- predict(df.bart, newdata = cfDF2, mc.cores = 3)
#
# test.pred1  <- predict(df.bart, newdata = cfDF, mc.cores = 3)
#
# mean0 <- tapply(test.pred0$surv.test.mean, rep(test.pred0$times, nrow(cfDF)), mean)
# mean1 <- tapply(test.pred1$surv.test.mean, rep(test.pred0$times, nrow(cfDF)), mean)
# or <- (mean1/(1-mean1))/(mean0/(1-mean0))
# plot(test.pred0$times, mean0)
# points(test.pred0$times, mean1, col = "gray")
#
# plot(test.pred0$times, or)

#### Heatmap OR plots ####
approx.coef <- function(coef, times) {
  approx_mat <- array(NA, dim = dim(coef), dimnames = dimnames(coef))
  int.time <- ceiling(min(times)):floor(max(times))
  out <- apply(coef, 2:3, function(y)
    approx(y = y, x = times, xout = int.time)$y)
  return(list(coef = out, time = int.time))
}

approx_local <- approx.coef(local_gam, preDF$times)
local_gam_mean <- apply(approx_local$coef, 1:2, mean)



local.df <- data.frame(value = exp(c(local_gam_mean)),
                       coef = factor(rep(dimnames(local_gam_mean)[[2]], each = dim(local_gam_mean)[1]),
                                                                    levels = dimnames(local_gam_mean)[[2]]),
                        time = approx_local$time)

heat.local <- ggplot(local.df,
                      aes(x = time,
                          y = coef,
                          fill = value)) +
  geom_tile() + theme_bw() +
  viridis::scale_color_viridis(name="Odds\nRatio",option="inferno",
                               trans = "log10", labels = scales::comma) +
  viridis::scale_fill_viridis(name="O.R. of\nSurvival",option="inferno",
                              trans = "log10", labels = scales::comma,
                              direction = -1) +
  xlab("Time (days)") + ylab("Variable")
# heat.local

get_mode <- function(x) {
  u   <- unique(x)
  tab <- table(x)
  return(u[which.max(tab)])
}


approx_global <- approx.coef(global_gam, preDF$times)
global_gam_mean <- apply(approx_global$coef, 1:2, mean)
global_gam_rank_list <- lapply(1:100, function(i) apply(approx_global$coef[,,i],1,rank))
global_gam_rank <- array(NA, dim = dim(approx_global$coef), dimnames = dimnames(approx_global$coef))
# for(i in 1:100) global_gam_rank[,,i] <- global_gam_rank_list[[i]]
# global_gam_rank_median <- apply(global_gam_rank, 1:2, median)

global.df <- data.frame(value = exp(c(global_gam_mean)), coef = factor(rep(dimnames(global_gam_mean)[[2]], each = dim(global_gam_mean)[1]),
                                                             levels = dimnames(global_gam_mean)[[2]]),
                       time = approx_global$time)

heat.global <- ggplot(global.df,
                      aes(x = time,
                         y = coef,
                         fill = value)) +
  geom_tile() + theme_bw() +
  viridis::scale_color_viridis(name="Odds\nRatio",option="inferno") +
  viridis::scale_fill_viridis(name="O.R. of\nSurvival",option="inferno",
                              trans = "log10", direction = -1) +
  xlab("Time (days)") + ylab("Variable")

heat.legend <- get_legend(heat.global, location = "right")
heat.grob <- arrangeGrob(heat.global + theme(legend.position = "none") +
                           ggtitle("Training Data") + xlab("") +
                           scale_y_discrete(labels = c("Age", "Gender", expression(KPS >= 90),
                                                       "Methylation", "Resection:\nBiopsy",
                                                       "Resection:\nSub-total", "Resection:\nTotal"),
                                            expand = c(0,0)),
                         heat.local + ggtitle("Test Data") +
                           theme(axis.text.y = element_blank()) +
                           ylab("") + xlab("") +
                           scale_y_discrete(expand=c(0,0)),
                         nrow = 1,
                         bottom = grid::textGrob("Time (days)",
                                                 vjust = -1.8, hjust=0.6),
                         widths = c(3.5,4))

pdf(file = "inst/figure/applied/GBM/heat_total.pdf", width = 7, height = 3.5)
grid.arrange(heat.grob)
dev.off()


# heat.rank.global <- data.frame(value = (c(global_gam_rank_median)),
#                                coef = factor(rep(dimnames(global_gam_rank_median)[[2]],
#                                                  each = dim(global_gam_rank_median)[1]),
#                                   levels = dimnames(global_gam_rank_median)[[2]]),
#                                time = approx_global$time) %>%
#   ggplot(
#           aes(x = time,
#               y = coef,
#               fill = value)) +
#   geom_tile() + theme_bw() +
#   viridis::scale_color_viridis(name="Coefficient Value",option="inferno") +
#   viridis::scale_fill_viridis(name="Coefficient Value",option="inferno") +
#   xlab("Time (days)") + ylab("Variable")
# heat.rank.global

#### Heatmap rank plots ####
approx.coef <- function(coef, times) {
  approx_mat <- array(NA, dim = dim(coef), dimnames = dimnames(coef))
  int.time <- ceiling(min(times)):floor(max(times))
  out <- apply(coef, 2:3, function(y)
    approx(y = y, x = times, xout = int.time)$y)
  return(list(coef = out, time = int.time))
}

approx_local <- approx.coef(local_gam, preDF$times)
local_gam_mean <- apply(approx_local$coef, 1:2, mean)
local_gam_rank <- t(apply(abs(local_gam_mean), 1, rank))


local.df <- data.frame(value = as.factor(c(local_gam_rank)),
                       coef = factor(rep(dimnames(local_gam_rank)[[2]], each = dim(local_gam_rank)[1]),
                                     levels = dimnames(local_gam_rank)[[2]]),
                       time = approx_local$time)

heat.local <- ggplot(local.df,
                     aes(x = time,
                         y = coef,
                         fill = value)) +
  geom_tile() + theme_bw() +
  scale_color_brewer(name="Rank",
                     type = "seq", palette = "YlGnBu") +
  scale_fill_brewer(name="Rank",
                     type = "seq", palette = "YlGnBu") +
  xlab("Time (days)") + ylab("Variable")
# heat.local

get_mode <- function(x) {
  u   <- unique(x)
  tab <- table(x)
  return(u[which.max(tab)])
}


approx_global <- approx.coef(global_gam, preDF$times)
global_gam_mean <- apply(approx_global$coef, 1:2, mean)
# global_gam_rank_list <- lapply(1:100, function(i) apply(approx_global$coef[,,i],1,rank))
# global_gam_rank <- array(NA, dim = dim(approx_global$coef), dimnames = dimnames(approx_global$coef))
# for(i in 1:100) global_gam_rank[,,i] <- global_gam_rank_list[[i]]
# global_gam_rank_median <- apply(global_gam_rank, 1:2, median)
global_gam_rank <- t(apply(abs(global_gam_mean),1,rank))

global.df <- data.frame(value = as.factor(c(global_gam_rank)), coef = factor(rep(dimnames(global_gam_rank)[[2]], each = dim(global_gam_rank)[1]),
                                                                       levels = dimnames(global_gam_rank)[[2]]),
                        time = approx_global$time)

heat.global <- ggplot(global.df,
                      aes(x = time,
                          y = coef,
                          fill = value)) +
  geom_tile() + theme_bw() +
  scale_color_brewer(name="Rank",
                     type = "seq", palette = "YlGnBu") +
  scale_fill_brewer(name="Rank",
                    type = "seq", palette = "YlGnBu") +
  xlab("Time (days)") + ylab("Variable")

heat.legend <- get_legend(heat.global, location = "right")
heat.grob <- arrangeGrob(heat.global + theme(legend.position = "none") +
                           ggtitle("Training Data") + xlab("") +
                           scale_y_discrete(labels = c("Age", "Gender", expression(KPS >= 90),
                                                       "Methylation", "Resection:\nBiopsy",
                                                       "Resection:\nSub-total", "Resection:\nTotal"),
                                            expand = c(0,0)),
                         heat.local + ggtitle("Test Data") +
                           theme(axis.text.y = element_blank()) +
                           ylab("") + xlab("") +
                           scale_y_discrete(expand=c(0,0)),
                         nrow = 1,
                         bottom = grid::textGrob("Time (days)",
                                                 vjust = -1.8, hjust=0.6),
                         widths = c(3.5,4))

pdf(file = "inst/figure/applied/GBM/rank_total.pdf", width = 7, height = 3.5)
grid.arrange(heat.grob)
dev.off()


# heat.rank.global <- data.frame(value = (c(global_gam_rank_median)),
#                                coef = factor(rep(dimnames(global_gam_rank_median)[[2]],
#                                                  each = dim(global_gam_rank_median)[1]),
#                                   levels = dimnames(global_gam_rank_median)[[2]]),
#                                time = approx_global$time) %>%
#   ggplot(
#           aes(x = time,
#               y = coef,
#               fill = value)) +
#   geom_tile() + theme_bw() +
#   viridis::scale_color_viridis(name="Coefficient Value",option="inferno") +
#   viridis::scale_fill_viridis(name="Coefficient Value",option="inferno") +
#   xlab("Time (days)") + ylab("Variable")
# heat.rank.global
