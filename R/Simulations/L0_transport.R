# Test transport methods

#### Load Packages ####
library(SLIMpaper)
set.seed(851090325) #from random.org
#### Data Generating Parameters
n <- 2^10
p <- 11
target <- get_normal_linear_model()
param <- target$rparam()
beta <- param$theta
sigma2 <- param$sigma2
corrs <- c(0,0.5,0.9)


#### Transportation methods ####
transp.methods <- c("exact","sinkhorn","greenkhorn",
                    "randkhorn", "gandkhorn",
                    "hilbert")

#### Generate Data ####
x <- target$X$rX(n, corrs[1], p)
data.list <- target$rdata(n,x,beta,sigma2, corrs[1])
y <- data.list$Y

#### posterior ####
n.samps <- 1000
hyperparameters <- list(mu = rep(0,p), sigma = diag(1,p,p), alpha = 1, beta = 1)
posterior.method <- "conjugate"
post <- target$rpost(n.samps, x, y, hyperparameters, method = posterior.method)

#### test transport methods ####
cl <- parallel::makeCluster(parallel::detectCores()-1)
doParallel::registerDoParallel(cl)
ex <- WPL0(x, post$eta, post$theta, p = 2, ground_p =2, method = "selection.variable", transport.method = "exact", parallel = cl)
un <-  WPL0(x, post$eta, post$theta, p = 2, ground_p =2, method = "selection.variable", transport.method = "univariate.approximation.pwr",
            parallel = cl)
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)

sapply(1:11, function(i) all.equal(ex$minCombPerActive[[i]], un$minCombPerActive[[i]]))

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
label.list <- outputs.list <- list()

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
for(method in transport.methods) {
  for(i in seq_along(outputs.list[[method]])) {
    sel <- outputs.list[[method]][[i]]$W2_dist$mean$dist.Selection[1:11]
    proj <- outputs.list[[method]][[i]]$W2_dist$mean$dist.Projection[1:11]
    outputs.list[[method]][[i]]$W2_dist$mean$dist <- c(sel,proj)
    outputs.list[[method]][[i]]$W2_dist$mean$groups <- factor(rep(c("Selection","Projection"), each = 11))
    outputs.list[[method]][[i]]$W2_dist$mean$groups <- factor(method)
    outputs.list[[method]][[i]]$W2_dist$mean <- outputs.list[[method]][[i]]$W2_dist$mean[1:11,]

    sel <- outputs.list[[method]][[i]]$mse$mean$dist.Selection[1:11]
    proj <- outputs.list[[method]][[i]]$mse$mean$dist.Projection[1:11]
    outputs.list[[method]][[i]]$mse$mean$dist <- c(sel,proj)
    outputs.list[[method]][[i]]$mse$mean$groups <- factor(rep(c("Selection","Projection"), each = 11))
    outputs.list[[method]][[i]]$mse$mean$groups <- factor(method)
    outputs.list[[method]][[i]]$mse$mean <- outputs.list[[method]][[i]]$mse$mean[1:11,]
  }
}
label.df <- do.call("rbind",label.list)
w2insamp <- combine.dist.compare(unlist(lapply(outputs.list, function(o) lapply(o, function(o2) o2$W2_dist)),recursive=FALSE))
mseinsamp <- combine.dist.compare(unlist(lapply(outputs.list, function(o) lapply(o, function(o2) o2$mse)),recursive=FALSE))

#### plots ####
pw2 <- plot(w2insamp, alpha = 0.8, ribbon = FALSE, ylab = "2-Wasserstein", base_size = 11)
pmse <- plot(mseinsamp, alpha = 0.8, ribbon = FALSE, ylab = "MSE", xlab = "", base_size = 11)

#### Save plots ####
w2mean <- pw2$mean + ggplot2::theme(legend.position = "none")
msemean <- pmse$mean
filename <- paste0(distfn,"_", corrfn,".pdf" )
pdf(file.path(figure.path, filename),
    width = 7.5, height = 3.5)
gridExtra::grid.arrange(w2mean, msemean,nrow=1)
dev.off()
