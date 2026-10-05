require(SLIMpaper)
library(dplyr)
library(forcats)

family <- "binomial"
penalty <- "mcp.net"
pf <- "none"
method <- "hilbert"
corr <- "Corr_0" #Corr_0, Corr_0.5, Corr_0.9
date <- "2020-06-15 12:00:00"
n <- 131072 #1024, 131072
p <- 21

folder <- file.path("Output", "Single_point", family, penalty, pf, method, corr,n,p)
files <- list.files(folder, full.names = TRUE)
files <- files[date.fun(files, date.start = date)]
niter <- length(files)

outputs <- lapply(files, readRDS)

label.df_sel <- data.frame(iter = rep(1:niter, each = 5),
                           transport.method = method,
                           method = factor(rep(c("Binary Programming","W2","Hahn-Carvalho","Stepwise","Simulated Annealing"),
                                               niter))
)

label.df_proj <- data.frame(iter = rep(1:niter, each = 6),
                            transport.method = method,
                            method = factor(rep(c("W2","Hahn-Carvalho","Stepwise","Simulated Annealing", "W1", "WInfty"),
                                                niter))
)


dist.list <- list(posterior = list(inSamp = NULL,
                                   single = NULL),
                  mean = list(inSamp = NULL,
                              single = NULL))
mses <- dist.list
distance <- dist.list

for(i in c("inSamp", "newX","single")) {
  for(j in c("posterior","mean")){
    mses[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$mse[[i]][[j]]))
    distance[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$W2_dist[[i]][[j]]))
  }
}

#### Single Point ####
w2single <- msesingle <- list()
w2single$selection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$single$selection))
w2single$projection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$single$projection))
msesingle$selection <- combine.dist.compare(lapply(outputs, function(o) o$mse$single$selection))
msesingle$projection <- combine.dist.compare(lapply(outputs, function(o) o$mse$single$projection))

#
w2_opt <- mse_opt <- list()
# cols <- c("dist","nactive", "groups", "method")
levs <- c("Binary Programming", "Relaxed B.P.","L1","L2","LInf")

#Wass
w2_opt$selection <- lapply(w2single$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
w2_opt$selection$p <- w2single$selection$p
w2_opt$projection <- lapply(w2single$projection[1:2], function(res) res[res$groups %in% c("L1", "Lasso", "LInf"),])
w2_opt$projection$p <- w2single$projection$p

for(pp in names(w2_opt)) {
  for(nn in names(w2_opt$selection)[1:2]) {
    if(is.null(w2_opt[[pp]][[nn]])) next
    w2_opt[[pp]][[nn]]$groups <- droplevels(w2_opt[[pp]][[nn]]$groups)

    w2_opt[[pp]][[nn]]$groups <- droplevels(w2_opt[[pp]][[nn]]$groups)
    if (pp == "projection") {
      levels(w2_opt[[pp]][[nn]]$groups) <- c("L1","L2","LInf")
      w2_opt[[pp]][[nn]]$groups <- factor(as.character(w2_opt[[pp]][[nn]]$groups), levels = levs)
    } else {
      levels(w2_opt[[pp]][[nn]]$groups) <- levs
    }


  }
}
class(w2_opt$selection) <- class(w2_opt$projection) <- class(w2single$selection)

w2_opt$combined <- lapply(names(w2_opt$selection)[1:2], function(nn) do.call("rbind", list(w2_opt$selection[[nn]], w2_opt$projection[[nn]])))
w2_opt$combined[[3]] <- w2_opt$selection$p
names(w2_opt$combined) <- c("posterior", "mean","p")
class(w2_opt$combined) <- class(w2single$selection)

# mse
mse_opt$selection <- lapply(msesingle$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
mse_opt$selection$p <- msesingle$selection$p
mse_opt$projection <- lapply(msesingle$projection[1:2], function(res) res[res$groups %in% c("L1","Lasso","LInf"),])
mse_opt$projection$p <- msesingle$projection$p

for(pp in names(mse_opt)) {
  for(nn in names(msesingle$selection)[1:2]) {

    mse_opt[[pp]][[nn]]$groups <- droplevels(mse_opt[[pp]][[nn]]$groups)

    mse_opt[[pp]][[nn]]$groups <- droplevels(mse_opt[[pp]][[nn]]$groups)
    if (pp == "projection") {
      levels(mse_opt[[pp]][[nn]]$groups) <- c("L1","L2","LInf")
      mse_opt[[pp]][[nn]]$groups <- factor(as.character(mse_opt[[pp]][[nn]]$groups), levels = levs)
    } else {
      levels(mse_opt[[pp]][[nn]]$groups) <- levs
    }
  }
}

class(mse_opt$selection) <- class(mse_opt$projection) <- class(msesingle$selection)

mse_opt$combined <- lapply(names(mse_opt$selection)[1:2], function(nn) do.call("rbind", list(mse_opt$selection[[nn]], mse_opt$projection[[nn]])))
mse_opt$combined[[3]] <- mse_opt$selection$p
names(mse_opt$combined) <- c("posterior", "mean","p")
class(mse_opt$combined) <- class(msesingle$selection)

# debugonce(plot.combine.dist.compare)
plot(w2single$selection, alpha = 0.5)
plot(msesingle$selection, alpha = 0.5)
plot(w2single$projection, alpha = 0.5)
plot(msesingle$projection, alpha = 0.5)

if(family == "binomial") {
  mse_opt$combined$mean <- mse_opt$combined$mean %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
  w2_opt$combined$mean <- w2_opt$combined$mean %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
}

msep <- plot(mse_opt$combined, alpha = 0.5, base_size = 20, ylab = "MSE")
w2p <- plot(w2_opt$combined,  alpha = 0.5, base_size = 20, ylab = "2-Wasserstein Distance")

if(family == "binomial") {
  msep$mean <- msep$mean + ggplot2::scale_color_manual(values = ggsci::pal_jama()(4)[2:4])
  w2p$mean <-  w2p$mean + ggplot2::scale_color_manual(values = ggsci::pal_jama()(4)[2:4])
}

w2fn <- file.path("inst", "figure","rawplots", paste0(paste0(c("single",family,penalty , pf ,method, corr, "w2plot"), collapse="_"), ".rds"))
saveRDS(list(mse = msep, w2 = w2p), file = w2fn)


# pdf(file = "inst/figure/simulation/mse_normal_single.pdf", width = 7, height = 5)
# msep$mean
# dev.off()
#
# pdf(file = "inst/figure/simulation/w2_normal_single.pdf", width = 7, height = 5)
# w2p$mean
# dev.off()


#### W2 r2 stuff ####
w2r2single <- w1r2single <- list()
fix.fun <- function(oo) {
  class(oo) <- c("WPR2", class(oo))
  return(oo)
}

w2r2single$selection  <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W2_r2$null$single$selection)))
w2r2single$projection <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W2_r2$null$single$projection)))
w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W1_r2$null$single$selection)))
w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W1_r2$null$single$projection)))

plot(w2r2single$selection)
plot(w2r2single$projection)

w2r2 <- combine.WPR2(w2r2single$selection %>% filter(groups %in% c("Binary Programming",  "Lasso")) %>%
                       mutate(groups = fct_recode(groups, "Relaxed B.P." = "Lasso")) %>%
                       mutate(groups = fct_drop(groups)),
                     w2r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
                       mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
                       mutate(groups = fct_drop(groups)))

w1r2 <- combine.WPR2(w1r2single$selection %>% filter(groups %in% c("Binary Programming",  "Lasso")) %>%
                       mutate(groups = fct_recode(groups, "Relaxed B.P." = "Lasso")) %>%
                       mutate(groups = fct_drop(groups)),
                     w1r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
                       mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
                       mutate(groups = fct_drop(groups)))

# pdf(file = "inst/figure/simulation/w2r2_normal_single.pdf", width = 7, height = 5)
# plot(w2r2, alpha = 0.5, base_size = 20)
# dev.off()
#
# pdf(file = "inst/figure/simulation/w1r2_normal_single.pdf", width = 7, height = 5)
# plot(w1r2, alpha = 0.5, base_size = 20)
# dev.off()

if(family == "binomial") {
  w1r2 <- w1r2 %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
  w2r2 <- w2r2 %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
}

w1r2p <- plot(w1r2, alpha = 0.5, base_size = 20)
w2r2p <- plot(w2r2, alpha = 0.5, base_size = 20)


if(family == "binomial") {
  w1r2p <- w1r2p + ggplot2::scale_color_manual(values = ggsci::pal_jama()(4)[2:4])
  w2r2p <-  w2r2p + ggplot2::scale_color_manual(values = ggsci::pal_jama()(4)[2:4])
}

wpr2fn <- file.path("inst", "figure","rawplots", paste0(paste0(c("single", family,penalty , pf ,method, corr, "wpr2plot"), collapse="_"), ".rds"))
saveRDS(list(w1r2 = w1r2p, w2r2 = w2r2p), file = wpr2fn)
