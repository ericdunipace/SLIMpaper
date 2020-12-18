#### Load Packages ####
library(SLIMpaper)

#### File save ####
figure.path <- file.path("inst","figure","transport")
corrfn <- "0_0"
distfn <- "dist_transport"
# msefn <- "mse_transport"

#### Load Data ####
family <- "gaussian"
penalty <- "lasso"
pf <- "none"
transport.methods <- c("exact", "sinkhorn","greenkhorn","gandkhorn","hilbert")
# method <- "univariate.approximation.pwr"
corr <- "Corr_0"
date <- "2019-11-15 23:00:00"
n <- 1024
p <- 11
label.list <- outputs.list <- outputs.list_limea <- list()

for(method in transport.methods) {
  folder <- file.path("Output",family, penalty, pf, method, corr,n,p)
  files <- list.files(folder, full.names = TRUE)
  files <- files[date.fun(files, date.start = date)]
  niter <- length(files)

  outputs.list[[method]] <- lapply(files, readRDS)

  label.list[[method]] <- data.frame(iter = rep(1:niter, each = 5),
                                     transport.method = method,
                                     method = "L0")

}
outputs.list_limea <- outputs.list
for(method in transport.methods) {
  for(i in seq_along(outputs.list[[method]])) {
    for(j in c("W2_dist","mse")) {
      sel <- outputs.list[[method]][[i]][[j]]$mean$dist.Selection[1:p]
      proj <- outputs.list[[method]][[i]][[j]]$mean$dist.Projection[1:p]
      if(nrow(outputs.list[[method]][[i]][[j]]$mean) != 2*p){
        outputs.list[[method]][[i]][[j]]$mean <- rbind(outputs.list[[method]][[i]][[j]]$mean,
                                                          outputs.list[[method]][[i]][[j]]$mean)
      }
      outputs.list[[method]][[i]][[j]]$mean$dist <- c(sel,proj)
      outputs.list[[method]][[i]][[j]]$mean$method <- factor(rep(c("Selection","Projection"), each = p))
      outputs.list[[method]][[i]][[j]]$mean$groups <- factor(method)
      outputs.list_limea[[method]][[i]][[j]]$mean <- outputs.list[[method]][[i]][[j]]$mean[(p+1):(2*p),]
      outputs.list[[method]][[i]][[j]]$mean <- outputs.list[[method]][[i]][[j]]$mean[1:p,]
    }


    # sel <- outputs.list[[method]][[i]]$mse$mean$dist.Selection[1:11]
    # proj <- outputs.list[[method]][[i]]$mse$mean$dist.Projection[1:11]
    # if(nrow(outputs.list[[method]][[i]]$mse$mean) != 22){
    #   outputs.list[[method]][[i]]$mse$mean <- rbind(outputs.list[[method]][[i]]$mse$mean,
    #                                                     outputs.list[[method]][[i]]$mse$mean)
    # }
    # outputs.list[[method]][[i]]$mse$mean$dist <- c(sel,proj)
    # outputs.list[[method]][[i]]$mse$mean$groups <- factor(rep(c("Selection","Projection"), each = 11))
    # # outputs.list[[method]][[i]]$mse$mean$groups <- factor(method)
    # outputs.list[[method]][[i]]$mse$mean <- outputs.list[[method]][[i]]$mse$mean[1:11,]
  }
}
label.df <- do.call("rbind",label.list)
w2insamp <- combine.dist.compare(unlist(lapply(outputs.list, function(o) lapply(o, function(o2) o2$W2_dist)),recursive=FALSE))
mseinsamp <- combine.dist.compare(unlist(lapply(outputs.list, function(o) lapply(o, function(o2) o2$mse)),recursive=FALSE))

w2insamp_limea <- combine.dist.compare(unlist(lapply(outputs.list_limea, function(o) lapply(o, function(o2) o2$W2_dist)),recursive=FALSE))
mseinsamp_limea <- combine.dist.compare(unlist(lapply(outputs.list_limea, function(o) lapply(o, function(o2) o2$mse)),recursive=FALSE))

#### plots ####
pw2 <- plot(w2insamp, alpha = 0.8, ribbon = FALSE, ylab = "2-Wasserstein",xlab = "", base_size = 11)
pmse <- plot(mseinsamp, alpha = 0.8, ribbon = FALSE, ylab = "MSE",  xlab = "",base_size = 11)
pw2_la <- plot(w2insamp_limea, alpha = 0.8, ribbon = FALSE, ylab = "2-Wasserstein", base_size = 11)
pmse_la <- plot(mseinsamp_limea, alpha = 0.8, ribbon = FALSE, ylab = "MSE",  base_size = 11)

#### Save plots ####
w2mean <- pw2$mean + ggplot2::theme(legend.position = c(0.79,0.7)) +
  ggplot2::ggtitle("LIME-p")
msemean <- pmse$mean + ggplot2::theme(legend.position = c(0.79,0.7))+
  ggplot2::ggtitle("LIME-p")
w2mean_a <- pw2_la$mean + ggplot2::theme(legend.position = "none")+
  ggplot2::ggtitle("LIME-a")
msemean_a <- pmse_la$mean +  ggplot2::theme(legend.position = "none")+
  ggplot2::ggtitle("LIME-a")

filename <- paste0(distfn,"_", "wass","_",corrfn,".pdf" )
pdf(file.path(figure.path, filename),
    width = 7.5, height = 3.5)
gridExtra::grid.arrange(w2mean_a, w2mean,nrow=1)
dev.off()

filename <- paste0(distfn,"_", "mse","_",corrfn,".pdf" )
pdf(file.path(figure.path, filename),
    width = 7.5, height = 3.5)
gridExtra::grid.arrange(msemean_a, msemean,nrow=1)
dev.off()
