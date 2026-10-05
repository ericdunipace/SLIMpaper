require(SLIMpaper)
require(SparsePosterior)
family <- "gaussian"
penalty <- "mcp.net"
pf <- "none"
method <- "hilbert"
corr <- "Corr_0" #Corr_0.5"
date <- "2020-01-09 12:00:00"
n <- 1024
p <- 101

folder <- file.path("Output",family, penalty, pf, method, corr,n,p)
files <- list.files(folder, full.names = TRUE)
files <- files[date.fun(files, date.start = date)]
niter <- length(files)

outputs <- lapply(files, readRDS)

label.df_sel <- data.frame(iter = rep(1:niter, each = 5),
                       transport.method = method,
                       method = factor(rep(c("Binary Programming","Lasso","Hahn-Carvalho","Stepwise","Simulated Annealing"),
                                    niter))
)

label.df_proj <- data.frame(iter = rep(1:niter, each = 4),
                           transport.method = method,
                           method = factor(rep(c("Lasso","Hahn-Carvalho","Stepwise","Simulated Annealing"),
                                               niter))
)
for(iter in 1:niter ) {
  for(name in c("selection", "projection")) {
    if(is.null(outputs[[iter]]$time[[name]]$anneal)) outputs[[iter]]$time[[name]]$anneal <- c("elapsed" = 64800)
  }
}
timings.sel <- cbind(label.df_sel, time = unlist(lapply(outputs, function(o) o$time$selection)),n,p)
timings.proj<- cbind(label.df_proj, time = unlist(lapply(outputs, function(o) o$time$projection)),n,p)
time.file.sel <- file.path("Output","Timings", paste0(c("timings_selection",family,corr,n,p,".rds"),collapse="_"))
time.file.proj <- file.path("Output","Timings", paste0(c("timings_projection",family,corr,n,p,".rds"),collapse="_"))
if(!dir.exists(dirname(time.file.sel))) dir.create(dirname(time.file.sel))
if(!dir.exists(dirname(time.file.proj))) dir.create(dirname(time.file.proj))
saveRDS(timings.proj, time.file.sel)
saveRDS(timings.proj, time.file.proj)

dist.list <- list(posterior = list(inSamp = NULL,
                                   newX = NULL,
                                   single = NULL),
                  mean = list(inSamp = NULL,
                              newX = NULL,
                              single = NULL))
mses <- dist.list
distance <- dist.list

for(i in c("inSamp", "newX","single")) {
  for(j in c("posterior","mean")){
    mses[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$mse[[i]][[j]]))
    distance[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$W2_dist[[i]][[j]]))
  }
}
w2insamp <- mseinsamp <- list()
w2insamp$selection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$inSamp$selection))
w2insamp$projection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$inSamp$projection))
mseinsamp$selection <- combine.dist.compare(lapply(outputs, function(o) o$mse$inSamp$selection))
mseinsamp$projection <- combine.dist.compare(lapply(outputs, function(o) o$mse$inSamp$projection))

# debugonce(plot.combine.dist.compare)
plot(w2insamp$selection, alpha = 0.5)
plot(mseinsamp$selection, alpha = 0.5)
plot(w2insamp$projection, alpha = 0.5)
plot(mseinsamp$projection, alpha = 0.5)

# filter for Lassos
w2_opt <- list()
w2_opt$selection <- lapply(w2insamp$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
w2_opt$projection <- lapply(w2insamp$projection[1:2], function(res) res[res$groups %in% c("Lasso"),])


# debugonce(plot_time)
plot_time(timings.sel)
plot_time(timings.proj)

#remove projection
rmout <- outputs
for(h in 1:niter){
  for(k in c("W2_dist","mse")){
    for(i in c("inSamp", "newX","single")) {
      for(j in c("posterior","mean")){
        idx <- which(rmout[[h]][[k]][[i]][[j]]$groups == "Projection")
        rmout[[h]][[k]][[i]][[j]] <- rmout[[h]][[k]][[i]][[j]][-idx,]
      }
    }
  }
}
rmtime <- timings[timings$method != "Projection",]

#plots for jsm
w2insamp <- combine.dist.compare(lapply(rmout, function(o) o$W2_dist$inSamp))
mseinsamp <- combine.dist.compare(lapply(rmout, function(o) o$mse$inSamp))

# debugonce(plot.combine.dist.compare)
plot(w2insamp, alpha = 0.3, base_size = 20, ylab = "2-Wasserstein Distance")
plot(mseinsamp, alpha = 0.3, base_size = 20, ylab = "MSE")

# debugonce(plot_time)
plot_time(rmtime, ylabs = "Minutes", alpha = 0.5, base_size = 20, scale.time=60)
plot_time(rmtime[rmtime$method %in% c("H.C.", "Selection"),], ylabs = "Seconds", alpha = 0.5, base_size = 20)


# debugonce(plot_ranks)
# plot_ranks(w2insamp, alpha = 0.1)
# plot_ranks(mseinsamp, alpha = 0.1, ylim = c(1,5))

#### Single Point ####
w2single <- msesingle <- list()
w2single$selection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$single$selection))
w2single$projection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$single$projection))
msesingle$selection <- combine.dist.compare(lapply(outputs, function(o) o$mse$single$selection))
msesingle$projection <- combine.dist.compare(lapply(outputs, function(o) o$mse$single$projection))

#
w2_opt <- mse_opt <- list()
# cols <- c("dist","nactive", "groups", "method")
levs <- c("Binary Programming", "Relaxed Lasso","Projection Lasso")

#Wass
w2_opt$selection <- lapply(w2single$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
w2_opt$selection$p <- w2single$selection$p
w2_opt$projection <- lapply(w2single$projection[1:2], function(res) res[res$groups %in% c("Lasso"),])
w2_opt$projection$p <- w2_opt$projection$p

for(pp in names(w2_opt)) {
  for(nn in names(w2single$selection)[1:2]) {

    w2_opt[[pp]][[nn]]$groups <- droplevels(w2_opt[[pp]][[nn]]$groups)

    w2_opt[[pp]][[nn]]$groups <- droplevels(w2_opt[[pp]][[nn]]$groups)
    if (pp == "projection") {
      levels(w2_opt[[pp]][[nn]]$groups) <- "Projection Lasso"
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
mse_opt$projection <- lapply(msesingle$projection[1:2], function(res) res[res$groups %in% c("Lasso"),])
mse_opt$projection$p <- msesingle$projection$p

for(pp in names(mse_opt)) {
  for(nn in names(msesingle$selection)[1:2]) {

    mse_opt[[pp]][[nn]]$groups <- droplevels(mse_opt[[pp]][[nn]]$groups)

    mse_opt[[pp]][[nn]]$groups <- droplevels(mse_opt[[pp]][[nn]]$groups)
    if (pp == "projection") {
      levels(mse_opt[[pp]][[nn]]$groups) <- "Projection Lasso"
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

msep <- plot(mse_opt$combined, alpha = 0.5, base_size = 20, ylab = "MSE")
w2p <- plot(w2_opt$combined,  alpha = 0.5, base_size = 20, ylab = "2-Wasserstein Distance")


pdf(file = "inst/figure/simulation/mse_normal.pdf", width = 7, height = 5)
msep$mean
dev.off()

pdf(file = "inst/figure/simulation/w2_normal.pdf", width = 7, height = 5)
w2p$mean
dev.off()
