rm(list=ls())
set.seed(236161905)
require(SLIMpaper)
figure.path <- file.path("inst","figure","toy_eg")
corrs <- c(0,0.5,0.9)
pres <- FALSE # is this for presentation
interactive <- FALSE # is this an interactive run

#### function to get ranks and plot ####
rankplot <- function(mod, color, label, base_size,
                     xlab = NULL, ylab = NULL) {
  which.sel.orig <- sapply(mod$theta, function(tt) as.numeric(rowSums(abs(tt)) > 0))
  # nzero <- mod$nzero[-which(colSums(which.sel.orig)==0)]
  rmv <- -which(colSums(which.sel.orig)==0)
  if(length(rmv) == 0) {
    nzero <-mod$nzero
  } else {
    nzero <-mod$nzero[rmv]
    which.sel.orig <- which.sel.orig[,rmv]
  }
  # which.sel.orig <- which.sel.orig[,-which(colSums(which.sel.orig)==0)]
  which.sel <- matrix(0, nrow=p, ncol=p)
  for(i in seq_along(nzero)){
    which.sel[,nzero[i]] <- which.sel.orig[,i]
  }


  included.df_singcompare <- data.frame(Included = c(t(which.sel)),
                                        Variable = factor(rep(c(1:p),each=p)),
                                        Active = factor(rep(1:p, p)))

  if(is.null(xlab)) xlab <- "Variable Number"
  if(is.null(ylab)) ylab <- "Number Active Coefficients"
  p <- ggplot2::ggplot(included.df_singcompare, ggplot2::aes(y=Active, x=Variable)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
    # ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::scale_fill_gradient(low = "white", high = color) +
    ggplot2::theme_bw(base_size) + ggplot2::ylab(ylab) +
    ggplot2::xlab(xlab) +
    ggplot2::theme(legend.position="none") +
    # ggplot2::ggtitle(bquote(L[1] ~ "Selection"))
    ggplot2::ggtitle(label)

  return(p)
}

#### Setup original data for posterior and for interpretable  ####
target <- get_normal_linear_model()

n.samp <- 1000
n <- 1024
p_star <- 6
base_size <- if(pres) {
  20
} else {
  11
}

param <- target$rparam()
param$theta <- param$theta[1:p_star]
param$theta <-c(2, -0.1, -0.2, 1.3, 1.4, 1.5)
X_single <- cbind(1,100,90,0.01,0.01,0.01)
cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)

for( corr in corrs){
  corr_fn <- gsub("[.]", "_",corr)
  cat("Correlation:\n")
  print(corr)
  target$X$corr <- corr
  X <- target$X$rX(n, target$X$corr, p_star)
  # newCorr <- diag(p-1)
  # newCorr[1:2,1:2] <- newCorr[3:4,3:4] <- newCorr[5:6,5:6] <- 0.5
  # diag(newCorr) <- 1
  # # X[,-1] <- X[,-1] %*% chol(newCorr)


  data <- target$rdata(n, X, c(param$theta),
                       param$sigma2, corr = 0, scale = FALSE)
  single_data <- target$rdata(1, X_single, c(param$theta),
                              param$sigma2, corr = 0, scale = FALSE)
  Y <- data$Y - 2
  single_Y <- single_data$Y - 2

  X <- X[,-1]
  X_sing <- X_single[,-1,drop=FALSE]
  p <- p_star-1

  hyperparameters <- list(mu = rep(0,p), sigma=diag(p), alpha = 10, beta = 10)
  post <- target$rpost(n.samp, X, Y, hyperparameters,
                       method = "conjugate", X.test = X_sing)

  cat("Variable Importance Order\n")
  print(WPVI(X=X_sing, Y=post$test$eta, theta=post$theta, mean.fun = NULL,
             p = 2, ground_p = 2, transport.method = "exact",
             parallel = cl))
  suffstat <- sufficientStatistics(X_sing, post$test$eta, post$theta, list(same = TRUE,
                                                                           method = "selection.variable",
                                                                           transport.method = "exact",
                                                                           epsilon = 0.05,
                                                                           niter = 100))
  selection <- W2L1(X_sing, NULL, post$theta, family="gaussian",
                    method = "selection.variable",
                    penalty="mcp.net", nlambda = 1e2, alpha = .99,
                    infimum.maxit = 1e2, gamma = 1.5,
                    maxit = 1e4,
                    transport.method = "hilbert", penalty.factor = 1/abs(suffstat$XtY),
                    lambda.min.ratio = 1e-10, display.progress = TRUE)#,
  plot.limep_L1 <- rankplot(mod = selection, color = "steelblue",
                            label = bquote(L[1]),
                            base_size = base_size,
                            xlab = "", ylab = "")
  if(interactive) print(plot.limep_L1)
  # which.sel.orig <- sapply(selection$theta, function(tt) as.numeric(rowSums(abs(tt)) > 0))
  # nzero <- selection$nzero[-which(colSums(which.sel.orig)==0)]
  # which.sel.orig <- which.sel.orig[,-which(colSums(which.sel.orig)==0)]
  # which.sel <- matrix(0, nrow=p, ncol=p)
  # for(i in seq_along(nzero)){
  #   which.sel[,nzero[i]] <- which.sel.orig[,i]
  # }
  #
  # # included.df_sing <- data.frame(Included = c(t(which.sel.orig)),
  # #                                Variable = factor(rep(c(1:p),each=length(nzero))),
  # #                                Active = factor(rep(nzero, p)))
  # # ggplot2::ggplot(included.df_sing, ggplot2::aes(y=Active, x=Variable)) +
  # #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  # #   ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  # #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  # #   ggplot2::theme(legend.position="none")
  #
  # included.df_singcompare <- data.frame(Included = c(t(which.sel)),
  #                                       Variable = factor(rep(c(1:p),each=p)),
  #                                       Active = factor(rep(1:p, p)))
  # plot.limep_L1 <- ggplot2::ggplot(included.df_singcompare, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[1] ~ "Selection"))

  #### IP (LIME-p) ####
  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  BP <- W2IP(X_sing, post$test$eta, post$theta, transport.method = "univariate.approximation.pwr",
             solution.method = "cone", display.progress = TRUE,
             parallel = cl)
  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(cl)
  plot.limep_BP <- rankplot(mod = BP, color = "forestgreen",
                            label = "Binary Program",
                            base_size = base_size,
                            xlab = "", ylab = "")
  if(interactive) print(plot.limep_BP)
  # which.sel.orig.ip <- sapply(IP$theta, function(tt) as.numeric(rowSums(abs(tt)) > 0))
  # rmv <- -which(colSums(which.sel.orig.ip)==0)
  # if(length(rmv) == 0) {
  #   nzero.ip <-IP$nzero
  # } else {
  #   nzero.ip <-IP$nzero[rmv]
  #   which.sel.orig.ip <- which.sel.orig.ip[,rmv]
  # }
  #
  # which.sel.ip <- matrix(0, nrow=p, ncol=p)
  # for(i in seq_along(nzero.ip)){
  #   which.sel.ip[,nzero.ip[i]] <- which.sel.orig.ip[,i]
  # }
  #
  # ip.df_sing <- data.frame(Included = c(t(which.sel.ip)),
  #                          Variable = factor(rep(c(1:p),each=p)),
  #                          Active = factor(rep(1:p, p)))
  #
  # plot.limep_IP <-ggplot2::ggplot(ip.df_sing, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "forestgreen") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle("I.P. Selection")

  #### L0 (LIME-p) ####

  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  l0ideal <- WPL0(X_sing, NULL, post$theta,p=2,ground_p = 2,
                  transport.method = "hilbert",
                  method = "selection.variable",
                  parallel = cl)
  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(cl)
  plot.limep_L0 <- rankplot(mod = l0ideal, color = "firebrick3",
                            label = bquote(L[0]),
                            base_size = base_size)
  if(interactive) print(plot.limep_L0)

  # which.ideal.l0 <- matrix(0, nrow=p, ncol=p)
  # for (i in 1:p) {
  #   which.ideal.l0[l0ideal$minCombPerActive[[i]],i] <- 1
  # }
  #
  # ideal.df_sing <- data.frame(Included = c(t(which.ideal.l0)),
  #                             Variable = factor(rep(c(1:p),each=p)),
  #                             Active = factor(rep(1:p, p)))
  #
  # plot.limep_L0 <-ggplot2::ggplot(ideal.df_sing, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "firebrick3") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[0] ~ "Selection"))

  filename <- paste0("selection_order_limep_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3)
  gridExtra::grid.arrange(plot.limep_L0,
                          plot.limep_L1, plot.limep_BP, nrow=1)
  dev.off()

  #### Ridge Plots (LIME-p) ####
  # rr <- ridgePlot(list("L0" = l0ideal,
  #                      "L1" = selection,
  #                      "B.P." = BP),
  #                 minCoef=1, maxCoef=6, full = c(post$test$eta), xlab="")
  # rr <- rr + ggridges::theme_ridges(base_size) + ggplot2::scale_fill_manual(
  #   breaks = c("L0","L1", "B.P."),
  #   labels = c(bquote(L[0]),bquote(L[1]), "B.P."),
  #   values = c("forestgreen", "firebrick3", "steelblue","red")) +
  #   ggplot2::ggtitle("LIME-p") + ggplot2::ylab("")
  # # rr$data$Method <- factor(rr$data$Method, levels = c("L0","L1"), labels = c(expression(L[0]), expression(L[1])))
  # filename <- paste0("ridge_plot_limep_",corr,".pdf")
  # pdf(file.path(figure.path, filename),
  #     width = 4, height = 3)
  # print(rr)
  # dev.off()

  #### LIME-a estimation ####
  Sigma <- cov(X)
  n_neighborhood <- 100
  X_neighborhood <- SLIMpaper::rmvnorm(n_neighborhood,
                                                    mean = X_sing,
                                                    covariance = Sigma/n)
  proj <-W2L1(X_neighborhood, NULL, post$theta, family="gaussian",
              method = "projection",
              penalty="mcp.net", nlambda = 1e2, alpha = .99,
              infimum.maxit = 1, gamma = 1.01,
              maxit = 1e4, lambda.min.ratio = 1e-10, display.progress = TRUE)

  plot.limea_L1 <- rankplot(mod = proj, color = "steelblue",
                            label = bquote(L[1] ~ "Selection"),
                            base_size = base_size)
  if(interactive) print(plot.limea_L1)

  # which.sel.orig.proj <- sapply(proj$theta, function(tt) as.numeric(rowSums(abs(tt)) > 0))
  # rmv <- -which(colSums(which.sel.orig.proj)==0)
  # if(length(rmv) == 0) {
  #   nzero.proj <-proj$nzero
  # } else {
  #   nzero.proj <-proj$nzero[rmv]
  #   which.sel.orig.proj <- which.sel.orig.proj[,rmv]
  # }
  #
  # which.sel.proj <- matrix(0, nrow=p, ncol=p)
  # for(i in seq_along(nzero.proj)){
  #   which.sel.proj[,nzero.proj[i]] <- which.sel.orig.proj[,i]
  # }
  #
  # included.df_proj <- data.frame(Included = c(t(which.sel.proj)),
  #                                Variable = factor(rep(c(1:p),each=p)),
  #                                Active = factor(rep(1:p, p)))
  # plot.limea_L1 <- ggplot2::ggplot(included.df_proj, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[1] ~ "Selection"))
  #
  #
  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  l0_limea <- WPL0(X_neighborhood, NULL, post$theta,p=2, ground_p = 2,
                   transport.method = "hilbert",
                   method = "projection",
                   parallel = cl)
  # doParallel::stopImplicitCluster()
  plot.limea_l0 <- rankplot(mod = l0_limea, color = "firebrick3",
                            label = bquote(L[0] ~ "Selection"),
                            base_size = base_size)
  if(interactive) print(plot.limea_l0)
  #
  # which.sel.orig.l0_limea <- sapply(l0_limea$theta, function(tt) as.numeric(rowSums(abs(tt)) > 0))
  # rmv <- -which(colSums(which.sel.orig.l0_limea)==0)
  # if(length(rmv) == 0) {
  #   nzero.l0_limea <-l0_limea$nzero
  # } else {
  #   nzero.l0_limea <-l0_limea$nzero[rmv]
  #   which.sel.orig.l0_limea <- which.sel.orig.l0_limea[,rmv]
  # }
  #
  # which.sel.l0_limea <- matrix(0, nrow=p, ncol=p)
  # for(i in seq_along(nzero.l0_limea)){
  #   which.sel.l0_limea[,nzero.l0_limea[i]] <- which.sel.orig.l0_limea[,i]
  # }
  #
  # included.df_l0_limea <- data.frame(Included = c(t(which.sel.l0_limea)),
  #                                    Variable = factor(rep(c(1:p),each=p)),
  #                                    Active = factor(rep(1:p, p)))
  # plot.limea_l0 <- ggplot2::ggplot(included.df_l0_limea, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "firebrick3") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[0] ~ "Selection"))

  filename <- paste0("selection_order_limea_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3)
  gridExtra::grid.arrange(plot.limep_L0,
                          plot.limep_L1, nrow=1, ncol=3)
  dev.off()

  #### Ridge Plots (LIME-a) ####
  # rr2 <- ridgePlot(list("L0" = l0_limea,
  #                       "L1" = proj),
  #                  minCoef=1, maxCoef=6, full = c(post$test$eta),
  #                  xlab="Posterior Predictive Mean",
  #                  scale = 4)
  # rr2 <- rr2 + ggridges::theme_ridges(base_size) +
  #   ggplot2::scale_fill_manual(
  #     breaks = c("L0","L1"),
  #     labels = c(bquote(L[0]),bquote(L[1])),
  #     values = c( "firebrick3", "steelblue","red")) +
  #   ggplot2::ggtitle("LIME-a") +
  #   ggplot2::theme(legend.position = "none")
  # # rr$data$Method <- factor(rr$data$Method, levels = c("L0","L1"), labels = c(expression(L[0]), expression(L[1])))
  # filename <- paste0("ridge_plot_limea_",corr,".pdf")
  # pdf(file.path(figure.path, filename),
  #     width = 3.25, height = 3)
  # print(rr2)
  # dev.off()
  scale <- if(corr == 0.5) {
    1.3
  } else {
    4
  }
  rr <- ridgePlot(list("LIME-a" = list("L0" = l0_limea,
                                       "L1" = proj),
                       "LIME-p" = list("L0" = l0ideal,
                                       "L1" = selection,
                                       "B.P." = BP)),
                  minCoef=1, maxCoef=5, full = c(post$test$eta),
                  xlab="Posterior Predictive Mean",
                  scale = scale)
  rr <- rr + ggridges::theme_ridges(base_size) + ggplot2::scale_fill_manual(
    breaks = c("L0","L1", "B.P."),
    labels = c(bquote(L[0]),bquote(L[1]), "B.P."),
    values = c("firebrick3", "steelblue","forestgreen","red"))
  # rr$data$Method <- factor(rr$data$Method, levels = c("L0","L1"), labels = c(expression(L[0]), expression(L[1])))
  filename <- paste0("ridge_plot_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3.5)
  print(rr)
  dev.off()
}

doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
