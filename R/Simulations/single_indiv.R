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

  if(is.null(xlab)) xlab <- "Variable number"
  if(is.null(ylab)) ylab <- "Number of active coefficients"
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


  Y <- X %*%  c(param$theta) + rnorm(n, param$sigma2) - 2
  single_Y <- X_single %*%  c(param$theta) + rnorm(1, param$sigma2) - 2
  # Y <- data$Y - 2
  # single_Y <- single_data$Y - 2

  X <- X[,-1]
  X_sing <- X_single[,-1,drop=FALSE]
  p <- p_star-1

  hyperparameters <- list(mu = rep(0,p), sigma=diag(p), alpha = 10, beta = 10)
  post <- target$rpost(n.samp, X, Y, hyperparameters,
                       method = "conjugate", X.test = X_sing)

  cat("Variable Importance Order\n")
  print(WPVI(X=X_sing, Y=post$test$eta, theta=post$theta,
       p = 2, ground_p = 2, transport.method = "exact",
       parallel = cl))
  # suffstat <- sufficientStatistics(X_sing, post$test$eta, post$theta, list(same = TRUE,
  #                                                                          method = "selection.variable",
  #                                                                          transport.method = "exact",
  #                                                                          epsilon = 0.05,
  #                                                                          niter = 100))
  penalty_fact <- set_penalty_factor(theta = post$theta, method = "covar", intercept = FALSE,
                                     x = X_sing, y = post$test$eta,
                                     transport.method = "exact")

  selection <- W2L1(X_sing, NULL, post$theta, family="gaussian",
                    method = "selection.variable",
                    penalty="mcp.net", nlambda = 1e2, alpha = .99,
                    infimum.maxit = 1e2, gamma = 1.5,
                    maxit = 1e4,
                    transport.method = "exact", penalty.factor = penalty_fact,
                    lambda.min.ratio = 1e-10, display.progress = TRUE)#,
  plot.slimp_L1 <- rankplot(mod = selection, color = ggsci::pal_jama()(5)[2],
                            label = bquote(atop(W[2] ~ "selection,",L[1]~"penalty")),
                            base_size = base_size,
                            xlab = "", ylab = "")
  if(interactive) print(plot.slimp_L1)
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
  # plot.slimp_L1 <- ggplot2::ggplot(included.df_singcompare, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[1] ~ "Selection"))

  #### IP (SLIM-p) ####
  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  BP <- W2IP(X_sing, post$test$eta, post$theta, transport.method = "exact",
             solution.method = "gurobi", display.progress = TRUE)
  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(cl)
  plot.slimp_BP <- rankplot(mod = BP, color = ggsci::pal_jama()(5)[1],
                            label = bquote(atop(W[2]~"selection,","B.P.")),
                            base_size = base_size,
                            xlab = "", ylab = "")
  if(interactive) print(plot.slimp_BP)
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
  # plot.slimp_IP <-ggplot2::ggplot(ip.df_sing, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "forestgreen") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle("I.P. Selection")

  #### L0 (SLIM-p) ####

  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  l0ideal <- WPL0(X_sing, NULL, post$theta,p=2,ground_p = 2,
                  transport.method = "univariate.approximation.pwr", #same as exact for univariate
                  method = "selection.variable")
  # doParallel::stopImplicitCluster()
  # parallel::stopCluster(cl)
  plot.slimp_L0 <- rankplot(mod = l0ideal, color = ggsci::pal_jama()(6)[6],
                            label = bquote(atop(W[2] ~ "selection,",L[0]~"Penalty")),
                            base_size = base_size)
  if(interactive) print(plot.slimp_L0)

  # which.ideal.l0 <- matrix(0, nrow=p, ncol=p)
  # for (i in 1:p) {
  #   which.ideal.l0[l0ideal$minCombPerActive[[i]],i] <- 1
  # }
  #
  # ideal.df_sing <- data.frame(Included = c(t(which.ideal.l0)),
  #                             Variable = factor(rep(c(1:p),each=p)),
  #                             Active = factor(rep(1:p, p)))
  #
  # plot.slimp_L0 <-ggplot2::ggplot(ideal.df_sing, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "firebrick3") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[0] ~ "Selection"))

  filename <- paste0("selection_order_slimp_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3)
  gridExtra::grid.arrange(plot.slimp_L0,
                          plot.slimp_L1, plot.slimp_BP, nrow=1)
  dev.off()

  #### Ridge Plots (SLIM-p) ####
  # rr <- ridgePlot(list("L0" = l0ideal,
  #                      "L1" = selection,
  #                      "B.P." = BP),
  #                 minCoef=1, maxCoef=6, full = c(post$test$eta), xlab="")
  # rr <- rr + ggridges::theme_ridges(base_size) + ggplot2::scale_fill_manual(
  #   breaks = c("L0","L1", "B.P."),
  #   labels = c(bquote(L[0]),bquote(L[1]), "B.P."),
  #   values = c("forestgreen", "firebrick3", "steelblue","red")) +
  #   ggplot2::ggtitle("SLIM-p") + ggplot2::ylab("")
  # # rr$data$Method <- factor(rr$data$Method, levels = c("L0","L1"), labels = c(expression(L[0]), expression(L[1])))
  # filename <- paste0("ridge_plot_slimp_",corr,".pdf")
  # pdf(file.path(figure.path, filename),
  #     width = 4, height = 3)
  # print(rr)
  # dev.off()

  #### SLIM-a estimation ####
  Sigma <- cov(X)
  n_neighborhood <- 100
  X_neighborhood <- SLIMpaper::rmvnorm(n_neighborhood,
                                                    mean = X_sing,
                                                    covariance = Sigma/n)
  proj <-W2L1(X_neighborhood, NULL, post$theta, family="gaussian",
              method = "projection",
              penalty="mcp.net", nlambda = 1e3, alpha = .99,
              infimum.maxit = 1, gamma = 1.01,
              maxit = 1e4, lambda.min.ratio = 1e-10, display.progress = TRUE)

  plot.slima_L1 <- rankplot(mod = proj, color = ggsci::pal_jama()(5)[4],
                            label = bquote(atop(W[2] ~ "projection,",
                                                L[1] ~ "Penalty")),
                            base_size = base_size)
  if(interactive) print(plot.slima_L1)

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
  # plot.slima_L1 <- ggplot2::ggplot(included.df_proj, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[1] ~ "Selection"))
  #
  #
  # cl <- parallel::makeCluster(parallel::detectCores()-1)
  # doParallel::registerDoParallel(cl)
  l0_slima <- WPL0(X_neighborhood, NULL, post$theta,p=2, ground_p = 2,
                   transport.method = "exact",
                   method = "projection")
  # doParallel::stopImplicitCluster()
  plot.slima_L0 <- rankplot(mod = l0_slima, color = ggsci::pal_jama()(6)[6],
                            label = bquote(atop(W[2] ~ "projection,",
                                                L[0]~"penalty")),
                            base_size = base_size)

  if(interactive) print(plot.slima_L0)
  #
  # which.sel.orig.l0_slima <- sapply(l0_slima$theta, function(tt) as.numeric(rowSums(abs(tt)) > 0))
  # rmv <- -which(colSums(which.sel.orig.l0_slima)==0)
  # if(length(rmv) == 0) {
  #   nzero.l0_slima <-l0_slima$nzero
  # } else {
  #   nzero.l0_slima <-l0_slima$nzero[rmv]
  #   which.sel.orig.l0_slima <- which.sel.orig.l0_slima[,rmv]
  # }
  #
  # which.sel.l0_slima <- matrix(0, nrow=p, ncol=p)
  # for(i in seq_along(nzero.l0_slima)){
  #   which.sel.l0_slima[,nzero.l0_slima[i]] <- which.sel.orig.l0_slima[,i]
  # }
  #
  # included.df_l0_slima <- data.frame(Included = c(t(which.sel.l0_slima)),
  #                                    Variable = factor(rep(c(1:p),each=p)),
  #                                    Active = factor(rep(1:p, p)))
  # plot.slima_l0 <- ggplot2::ggplot(included.df_l0_slima, ggplot2::aes(y=Active, x=Variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = Included), colour = "white") +
  #   ggplot2::scale_fill_gradient(low = "white", high = "firebrick3") +
  #   ggplot2::theme_bw(base_size) + ggplot2::ylab("Number Active Coefficients") + ggplot2::xlab("Variable Number") +
  #   ggplot2::theme(legend.position="none") + ggplot2::ggtitle(bquote(L[0] ~ "Selection"))

  filename <- paste0("selection_order_slima_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3)
  gridExtra::grid.arrange(plot.slima_L0,
                          plot.slima_L1, nrow=1, ncol=3)
  dev.off()

  #### combine all plots ####
  filename <- paste0("selection_order_slim_total_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3)
  gridExtra::grid.arrange(
                          plot.slimp_L0,
                          plot.slimp_BP,
                          plot.slimp_L1 ,
                          plot.slima_L1 + ggplot2::xlab("")+ggplot2::ylab(""),
                          nrow=1, ncol=4)
  dev.off()

  #### Ridge Plots (SLIM-a) ####
  # rr2 <- ridgePlot(list("L0" = l0_slima,
  #                       "L1" = proj),
  #                  minCoef=1, maxCoef=6, full = c(post$test$eta),
  #                  xlab="Posterior Predictive Mean",
  #                  scale = 4)
  # rr2 <- rr2 + ggridges::theme_ridges(base_size) +
  #   ggplot2::scale_fill_manual(
  #     breaks = c("L0","L1"),
  #     labels = c(bquote(L[0]),bquote(L[1])),
  #     values = c( "firebrick3", "steelblue","red")) +
  #   ggplot2::ggtitle("SLIM-a") +
  #   ggplot2::theme(legend.position = "none")
  # # rr$data$Method <- factor(rr$data$Method, levels = c("L0","L1"), labels = c(expression(L[0]), expression(L[1])))
  # filename <- paste0("ridge_plot_slima_",corr,".pdf")
  # pdf(file.path(figure.path, filename),
  #     width = 3.25, height = 3)
  # print(rr2)
  # dev.off()
  scale <- if(corr == 0.9) {
    1
  } else {
    3.7
  }
  # rr <- ridgePlot(list("SLIM-a" = list("L0" = l0_slima,
  #                                      "Relaxed B.P." = proj),
  #                       "SLIM-p" = list("L0" = l0ideal,
  #                                        "Relaxed B.P" = selection,
  #                                        "B.P." = BP)),
  #                 minCoef=1, maxCoef=5, full = c(post$test$eta),
  #                 xlab="Posterior Predictive Mean",
  #                 scale = scale)
  rr <- ridgePlot(list("L0" = l0ideal,
                       "B.P." = BP,
                       "Relaxed B.P." = selection,
                       "Proj" = proj),
                  minCoef=1, maxCoef=5, full = c(post$test$eta),
                  xlab="Posterior Predictive Mean",
                  scale = scale)
  rr <- rr + ggridges::theme_ridges(base_size) + ggplot2::scale_fill_manual(
    breaks = c("L0", "B.P.","Relaxed B.P.","Proj"),
    labels = c(bquote(W[2]~"Selection,"~L[0]~" Penalty"),
               bquote(W[2]~"Selection, B.P."),
               bquote(W[2]~"Selection,"~L[1] ~" Penalty"),
               bquote(W[2]~"Projection,"~L[1] ~" Penalty")),
    values = c(ggsci::pal_jama()(6)[c(6,1,2,4)], "red")
      # c("firebrick3", "steelblue","forestgreen","red")
    ) +
    ggplot2::theme(strip.background = ggplot2::element_blank())
  # rr$data$Method <- factor(rr$data$Method, levels = c("L0","L1"), labels = c(expression(L[0]), expression(L[1])))
  filename <- paste0("ridge_plot_",corr_fn,".pdf")
  pdf(file.path(figure.path, filename),
      width = 7.5, height = 3.5)
  print(rr)
  dev.off()
}

doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
