rm(list=ls())
set.seed(148021652) #from random.org

seed.bgam <- sample.int(.Machine$integer.max, size = 1)
seed.surv <- sample.int(.Machine$integer.max, size = 1)
seed.cens <- sample.int(.Machine$integer.max, size = 1)

print(seed.bgam)
# 359770598
print(seed.surv)
# 1772790390
print(seed.cens)
# 874712390

#### Load Packages ####
require(SLIMpaper)
require(survival)
require(rjags)
require(BART)
require(ggplot2)
require(gamm4)

#### Load Data ####
data.file <- "../Data/GBM/Merged.ndGBM.RCTS.RWE.Subset.OS.RData"
load(data.file)
Y$T.OS <- Y$T.OS + runif(nrow(Y), -0.001, 0.001) # add small epsilon to not worry about ties.
Y$Control <- droplevels(Y$Control)
X$Control <- droplevels(X$Control)

sourceDF <- X$Control=="rwe.c0"
DF.X <- X[sourceDF,]
DF.Y <- Y[sourceDF,]

PB.X <- X[!sourceDF,]
PB.Y <- Y[!sourceDF,]

DF.X$Gender <- abs(as.numeric(DF.X$Gender)-2)
DF.X$MGMT <- as.numeric(DF.X$MGMT)
DF.X$id <- 1:nrow(DF.X)

valid_test_idx <- sample(DF.X$id, ceiling(nrow(DF.X) * 0.1), replace = FALSE)
test_idx       <- valid_test_idx[1]
validate_idx   <- valid_test_idx[-1]

DF_train  <- list(X = DF.X[ !(DF.X$id %in% valid_test_idx), , drop = FALSE],
                  Y = DF.Y[!(DF.X$id %in% valid_test_idx),, drop = FALSE])
DF_valid  <- list(X = DF.X[validate_idx,, drop = FALSE],
                 Y = DF.Y[validate_idx,, drop = FALSE])
DF_test   <- list(X = DF.X[test_idx,, drop = FALSE],
                  Y = DF.Y[test_idx,, drop = FALSE])


mm.df <- model.matrix(~ scale(Age) + Gender + KPS.C + MGMT + Resection+ 0, data= DF_train$X)
n <- nrow (mm.df)
mm.df_test <- model.matrix(~ scale(Age, center = mean(DF_train$X$Age),
                                   scale = sd(DF_train$X$Age)) + Gender + KPS.C + MGMT + Resection+ 0,
                           data= rbind(DF_test$X, DF_valid$X))

test_binary <- interaction(data.frame(Gender = c(0,1), KPS.C = c(0,1),
                                      MGMT = c(0,1)))
bin.matrix <- do.call("rbind", strsplit(levels(test_binary), "[.]"))
storage.mode(bin.matrix) <- "integer"
cat.mat <- diag(1, 3, 3)
ham.check.mat <- cbind(bin.matrix, 0)
ham.check.mat <- ham.check.mat[c(1:8,1:8,1:8),]
ham.check.mat[1:8,4] <- 1
ham.check.mat[9:16,4] <- 2
ham.check.mat[17:24,4] <- 3
ham.test <- c(mm.df_test[1,2:4],2)
hamming.dist <- rep(NA, nrow(ham.check.mat))
for(i in 1:nrow(ham.check.mat)) hamming.dist[i] <- sum(ham.test != ham.check.mat[i,])


design.matrix <- cbind(bin.matrix, matrix(0, 8,3))
design.matrix <- design.matrix[c(1:8,1:8,1:8),]
design.matrix[1:8,4] <- design.matrix[9:16,5] <- design.matrix[17:24,6] <- 1
design.matrix <- design.matrix[hamming.dist <= 1,]

mm.df_neighb <- do.call("rbind", lapply(1:100, function(i) cbind(rnorm(nrow(design.matrix),
                                                                        mm.df_test[1,1], sd = 1/sqrt(n)),
                                                                  design.matrix)))



# mm.df_neighb <- rmvnorm(100, mean = c(mm.df_test[1,]), covariance = cov(mm.df)/n)
colnames(mm.df_neighb) <- colnames(mm.df_test) <- colnames(mm.df)

DF_test$neighb <- mm.df_neighb
save(DF_train, DF_valid, DF_test, file = "../Data/GBM/GBM_DF.RData")

#### Examine Data ####
KM <- survfit(Surv(DF_train$Y$T.OS, DF_train$Y$C.OS)~1)
plot(KM)

#### Check PH assumption ####
mlefit <- coxph(Surv(T.OS, C.OS) ~ Age + Gender + KPS.C + MGMT +
                  Resection, data=cbind(DF_train$X, DF_train$Y))
mlezph <- cox.zph(mlefit)
print(mlezph)
plot(mlezph) # looks like PH assumption might be violated

# try stratified
stratafit <- coxph(Surv(T.OS, C.OS) ~ Age + Gender +  MGMT +
                     Resection + strata(KPS.C), data=cbind(DF_train$X, DF_train$Y))
stratazph <- cox.zph(stratafit)
print(stratazph)
plot(stratazph) # looks like PH assumption might be ok but can still improve with non PH

# check brier score
pred <- survival::survfit(mlefit)
class(pred) <- "survfit"
surv.obj <- survival::Surv(DF_train$Y$T.OS, DF_train$Y$C.OS)
print(ipred::sbrier(surv.obj, pred)[1])
# 0.1211049

pred <- survival::survfit(mlefit, newdata = cbind(DF_valid$X, DF_valid$Y))
class(pred) <- "survfit"
surv.obj <- survival::Surv(DF_valid$Y$T.OS, DF_valid$Y$C.OS)
print(ipred::sbrier(surv.obj, pred)[1])
# 0.1690142

#### Set-up Parameters ####
out.dir <- "Output/GBM"
if(!dir.exists(out.dir)) dir.create(dir)
n.samp <- 2000
nchain <- 4
thin <- 100
gamma <- 1
nlambda <- 100
lambda.min.ratio <- 1e-10


#### GAM Model ####
target <- get_survival_linear_model()
gam.file <- file.path(out.dir, "gbm_model_gam_compare.rds")
cens_gam.file <- file.path(out.dir, "gbm_gam_cens.rds")
# debugonce(target$rpost)
# if (!file.exists(gam.file) | !file.exists(cens_gam.file)) {
#   gam_post <- target$rpost(n.samp = n.samp, x = mm.df[,-5],
#                            y = DF.Y$T.OS,
#                            fail = DF.Y$C.OS,
#                hyperparameters = list(mu = 0, sigma = 1),
#                id = (1:n), nchain = nchain, jags_dir = NULL,
#                stan_dir = NULL,
#                method = "inla-GAM", seed = sample.int(.Machine$integer.max,1),
#                parallel = FALSE, thin = thin,
#                X.test = mm.df[,-5])
#   gam_cens <- target$rpost(n.samp = n.samp, x = mm.df[,-5],
#                            y = DF.Y$T.OS,
#                            fail = (1-DF.Y$C.OS),
#                            hyperparameters = list(mu = 0, sigma = 1),
#                            id = (1:n), nchain = nchain, jags_dir = NULL,
#                            stan_dir = NULL,
#                            method = "inla-GAM", seed = sample.int(.Machine$integer.max,1),
#                            parallel = FALSE, thin = thin,
#                            X.test = mm.df[,-5])
#   saveRDS(gam_post, file = gam.file)
#   saveRDS(gam_cens, file = cens_gam.file)
# } else {
#   gam_post <- readRDS(gam.file)
#   gam_cens <- readRDS(cens_gam.file)
# }
# gam.df <- data.frame(y = DF.Y$T.OS, DF.X)
# gam.df$Age <- scale(gam.df$Age)
# gam.form <- formula(y ~ Resection +
#                       # Age + Gender + KPS +  MGMT +#Age:Resection +
#                       s(Age) +
#                       Gender +
#                       KPS.C +
#                       MGMT - 1)
# gam.mle <- mgcv::gam( gam.form, data = gam.df ,
#                          weight = DF.Y$C.OS, family = cox.ph(link = "identity"))
# survobj <- survival::Surv(gam.df$y, DF.Y$C.OS)
# pred <- predict(gam.mle, type="response")
# fit <- survival::survfit(survival::coxph(survobj ~ pred))
# class(fit) <- "survfit"
# ipred::sbrier(survobj, fit)
# #0.1389285
#
# df.gam.mle <- mgcv::gam( gam.form, data = gam.df ,
#                          weight = DF.Y$C.OS, family = cox.ph(link = "identity"),
#                          fit = FALSE)
# df.gam.post <- mgcv::ginla(G = df.gam.mle, A=NULL,nk=16,nb=100,J=1,interactive=FALSE,int=0,approx=0)

poisdf <- pois_dat(time = seq(0,max(DF_train$Y$T.OS), length.out=15),
                   event.times = DF_train$Y$T.OS,
                   event = DF_train$Y$C.OS,
                   x = DF_train$X)
gam.df <- data.frame(y = poisdf$y, time = poisdf$time,
                     offset = (poisdf$offset), poisdf$x,
                     event.time = factor(poisdf$time),
                     id = poisdf$id)
# gam.df <- gam_dat(times = DF_train$Y$T.OS,
#                   event.times = DF_train$Y$T.OS,
#                                      event = DF_train$Y$C.OS,
#                                      x = DF_train$X,
#                   n.time.groups = 15)
gam.df$Age <- scale(gam.df$Age)
gam.form <- formula(y ~ Resection  + #event.time +
                      # Age + Gender + KPS.C + MGMT +
                      s(time, by = Age) +
                      s(time, by = Gender) +
                      s(time, by = KPS.C) +
                      s(time, by = MGMT) +
                      s(time, by = Resection) - 1
                    )
gam.mle <- mgcv::bam( gam.form, data = gam.df ,
                         family = poisson(link = "log"),
                         discrete = TRUE,
                         fit = TRUE,
                      offset = log(gam.df$offset))
# gam.mle <- mgcv::bam( gam.form, data = gam.df ,
#                       family = poisson(link = "log"),
#                       discrete = TRUE,
#                       fit = TRUE)
predictor <- predict(gam.mle, type = "link")
lambda <- gam.mle$fitted.values

times.bam2 <- poisdf$time
times.bam1 <- rep(NA, length(times.bam2))
times.bam1.list <- (tapply(times.bam2, poisdf$id, function(x){
  if(length(x) > 1) {
    c(0, x[1:(length(x)-1)])
  } else {
    0
  }
}))
for(i in poisdf$id) {
  times.bam1[poisdf$id == i] <- times.bam1.list[[i]]
}


survobj.bam <- survival::Surv(times.bam1, times.bam2, poisdf$y)
fit <- survival::survfit(survival::coxph(survobj.bam ~ lambda))
class(fit) <- "survfit"
surv.obj <- survival::Surv(DF_train$Y$T.OS, DF_train$Y$C.OS)
print(ipred::sbrier(surv.obj, fit)[1])
#0.1220538

poisdf_val <- pois_dat(time = seq(0,max(DF_train$Y$T.OS), length.out=15),
                       event.times = DF_valid$Y$T.OS,
                       event = DF_valid$Y$C.OS,
                       x = DF_valid$X)
times.bam2 <- poisdf_val$time
times.bam1 <- rep(NA, length(poisdf_val))
times.bam1.list <- (tapply(times.bam2, poisdf_val$id, function(x){
  if(length(x) > 1) {
    c(0, x[1:(length(x)-1)])
  } else {
    0
  }
}))
for(i in poisdf_val$id) {
  times.bam1[poisdf_val$id == i] <- times.bam1.list[[i]]
}

gam.df_val <- data.frame(y = poisdf_val$y, time = poisdf_val$time,
                     offset = (poisdf_val$offset), poisdf_val$x,
                     event.time = factor(poisdf_val$time),
                     id = poisdf_val$id)
gam.df_val$Age <- scale(gam.df_val$Age,
                        center = attributes(gam.df$Age)$`scaled:center`,
                        scale = attributes(gam.df$Age)$`scaled:scale`)
lambda_valid <- predict(gam.mle, newdata = gam.df_val, type = "response") *gam.df_val$offset
survobj.bam <- survival::Surv(times.bam1, times.bam2, poisdf_val$y)
fit <- survival::survfit(survival::coxph(survobj.bam ~ lambda_valid))
class(fit) <- "survfit"
surv.obj <- survival::Surv(DF_valid$Y$T.OS, DF_valid$Y$C.OS)
print(ipred::sbrier(surv.obj, fit)[1])
#0.2862567

if(!file.exists(gam.file)) {
  df.gam.mle <- mgcv::bam( gam.form, data = gam.df ,
                           family = poisson(link = "log"),
                           discrete = TRUE,
                           fit = FALSE,
                           offset = log(gam.df$offset))
  warning("mgcv::ginla function is experimental. May change in future versions of mgcv.
          Used version 1.8-31 in R 3.6.3 to run analysis")
  set.seed(seed.bgam)
  gam.post <- mgcv::ginla(G = df.gam.mle, A=NULL,nk=100,nb=1000,J=1,
                          interactive=FALSE,int=0,approx=0) #note: integration fails for this data
  saveRDS(gam.post, file = gam.file)

} else {
  gam.post <- readRDS(file = gam.file)

}
dens.prod <- apply(gam.post$density,2,prod)
weights <- dens.prod/sum(dens.prod)

# see if coefficients seem reasonable
all.equal(gam.mle$coefficients, c(gam.post$beta %*% weights), check.attributes = FALSE)

Xp <- predict(gam.mle, newdata = gam.df_val, type = "lpmatrix")
bayes.link <- Xp %*% gam.post$beta
bayes.pred <- exp(bayes.link) * gam.df_val$offset


#unfortunately there are some low probability bad predictions that mess up the
#survfit calculation. using the mean prediction instead


# survobj <- survival::Surv(DF_valid$Y$T.OS, DF_valid$Y$C.OS)
# brier.vals <- sapply(1:ncol(bayes.pred), function(b) {
#   temp_mod <- survival::coxph(survobj.bam ~ bayes.pred[,b], iter.max = 1e9)
#   fit <- survival::survfit(temp_mod)
#   class(fit) <- "survfit"
#   ipred::sbrier(survobj, fit)
# })

# E_pred <- bayes.pred %*% weights
# cat("Check expectation close to MLE\n")
# c(MLE = lambda_valid[1], post.mean = E_pred[1])
# all.equal(c(lambda_valid), c(E_pred), check.attributes = FALSE)
#
# sum(brier.vals * weights)


approx_E_pred <- exp(Xp %*% c(gam.post$beta %*% weights)) * gam.df_val$offset
cat("Check expectation close to MLE use approximate prediction\n")
c(MLE = lambda_valid[1], post.mean = approx_E_pred[1])
all.equal(c(lambda_valid), c(approx_E_pred), check.attributes = FALSE)
fit <- survival::survfit(survival::coxph(survobj.bam ~ approx_E_pred))
class(fit) <- "survfit"
survobj <- survival::Surv(DF_valid$Y$T.OS, DF_valid$Y$C.OS)
ipred::sbrier(survobj, fit)
# 0.3686498

#### BART model DF.X ####
bart.file <- file.path(out.dir, "gbm_model.rds")
cens.file <- file.path(out.dir, "gbm_cens.rds")
if(!file.exists(bart.file) | !file.exists(cens.file)) {
  # df.bart <- mc.surv.bart(x.train = mm.df,
  #                      y.train = NULL,
  #                      times = DF.Y$T.OS,
  #                      delta = DF.Y$C.OS,
  #                      printevery = n.samp/10,
  #                      ndpost = n.samp,
  #                      nskip = n.samp * thin,
  #                      keepevery = thin,
  #                      k=100,
  #                      base = .95,
  #                      power = 0.25,
  #                      ntree = 50,
  #                      # id = 1:nrow(mm.df),
  #                      seed = seed.surv,
  #                      mc.cores = 3,
  #                      nice = 0L,
  #                      )
  df.bart <- mc.abart.fix(x.train = mm.df,
                          times = DF_train$Y$T.OS,
                          delta = DF_train$Y$C.OS,
                          x.test = rbind(mm.df_neighb, mm.df_test),
                          K = length(unique(DF_train$Y$T.OS)),
                          printevery = n.samp * thin/10,
                          ndpost = n.samp,
                          nskip = n.samp * thin,
                          keepevery = thin,
                          k=2,
                          base = .95,
                          power = .25,
                          ntree = 50,
                          seed = seed.surv,
                          mc.cores = 3,
                          nice = 0L,
  )
  # df.cens.bart <- mc.surv.bart(x.train = mm.df,
  #                           y.train = NULL,
  #                           times = DF.Y$T.OS,
  #                           delta = (1-DF.Y$C.OS),
  #                           printevery = n.samp/10,
  #                           ndpost = n.samp,
  #                           nskip = n.samp * thin,
  #                           keepevery = thin,
  #                           k=100,
  #                           base = .95,
  #                           power = 0.25,
  #                           ntree = 50,
  #                           id = 1:nrow(mm.df),
  #                           seed = seed.cens,
  #                           mc.cores = 3,
  #                           nice = 0L)
  df.cens.bart <- mc.abart.fix(x.train = mm.df,
                               times = DF_train$Y$T.OS,
                               delta = 1-DF_train$Y$C.OS,
                               x.test = rbind(mm.df_neighb,mm.df_test),
                               K = length(unique(DF_train$Y$T.OS)),
                               printevery = n.samp * thin/10,
                               ndpost = n.samp,
                               nskip = n.samp * thin,
                               keepevery = thin,
                               k=2,
                               base = .95,
                               power = 0.25,
                               ntree = 50,
                               seed = seed.cens,
                               mc.cores = 3,
                               nice = 0L)
  saveRDS(df.bart, file = bart.file)
  # saveRDS(df.bart, file = "Output/GBM/gbm_gen_bart_model.rds")
  saveRDS(df.cens.bart, file  = cens.file)
} else {
  df.bart <- readRDS(file = bart.file)
  df.cens.bart <- readRDS(file  = cens.file)
}

# preDF       <- surv.pre.bart(x.train = cbind(id=DF.X$id, mm.df),
#                              x.test = cbind(id=DF.X$id, mm.df),
#                              times = DF.Y$T.OS, delta = DF.Y$C.OS)
# testDF      <- preDF$tx.test[,-c(1:2)]
# predictions <- predict(df.bart, newdata = mm.df, mc.cores=3)
# debugonce(predict.abart) # currently right for yhat, wrong surv...
# predictions <- predict.abart(df.bart, newdata = mm.df, mc.cores=3)
interpY     <- t(qlogis(df.bart$surv.test))
saveRDS(interpY, file = file.path(out.dir, "gbm_model_Y.rds"))

interpY_train     <- t(qlogis(df.bart$surv.train))
saveRDS(interpY_train, file = file.path(out.dir, "gbm_model_Y_training.rds"))

#### Eval model performance ####
# cens.pred  <- predict(df.cens.bart, newdata = testDF, mc.cores = 3)
cens.pred  <- df.cens.bart$surv.test
#setup arrays
surv.array <- cens.array <- array(NA,
                                  dim = c(preDF$K,
                                          nrow(df.bart$yhat.test),
                                          n),
                                  dimnames = list(
                                    times = 1:preDF$K,
                                    iterations=1:nrow(df.bart$yhat.test),
                                    observations = 1:n)
                                  )

# for(i in 1:n) {
#     surv.array[,,i] <- t(predictions$surv.test[, (i-1) * preDF$K + 1:preDF$K])
#     cens.array[,,i] <- t(cens.pred$surv.test[, (i-1) * preDF$K + 1:preDF$K])
# }

for(i in 1:n) {
  surv.array[,,i] <- t(df.bart$surv.test[, (i-1) * preDF$K + 1:preDF$K])
  cens.array[,,i] <- t(df.cens.bart$surv.test[, (i-1) * preDF$K + 1:preDF$K])
}


#concordance
# concordance <- get_survival_linear_model()$evalfit(times = DF.Y$T.OS,
#                                                    event = DF.Y$C.OS,
#                                                    surv.times = preDF$times,
#                                                    surv = surv.array,
#                                                    method = "c-index")

#check brier score
# test.pred <- predict(df.bart, newdata = preDF$tx.train[,-2], mc.cores = 3)
# id <- preDF$tx.train[,"id"]
# time.bart.2 <- test.pred$tx.test[,"t"]
# time.bart.1 <- rep(NA, length(time.bart.2))
# time.bart.list <- (tapply(time.bart.2, id, function(x){
#   if(length(x) > 1) {
#     c(0, x[1:(length(x)-1)])
#   } else {
#     0
#   }
# }))
# for(i in id) {
#   time.bart.1[id == i] <- time.bart.list[[i]]
# }
# survobj.bart <- survival::Surv(time.bart.1, time.bart.2, preDF$y.train)
# surv.obj <- survival::Surv(DF.Y$T.OS, DF.Y$C.OS)
# yhat.rows <- rep(NA, n)
#
# for(i in 1:n) {
#   idx <- which(preDF$tx.test[1:preDF$K + (i-1)*preDF$K,"t"] == DF.Y$T.OS[i])
#   yhat.rows[i] <- (i-1)*preDF$K + idx
# }
#
# yhat.sub <- predictions$yhat.test[,yhat.rows]
# sbrier.bart.vals <- sapply(1:nrow(test.pred$yhat.test), function(b) {
#   fit <- survival::survfit(survival::coxph(survobj.bart ~ test.pred$yhat.test[b,]))
#   class(fit) <- "survfit"
#   ipred::sbrier(surv.obj, fit)
# })
# mean(sbrier.bart.vals)
# 0.1321027
# 0.1309973, k = 20, power = .25
# 0.1773833, k = 1, base = 1, power = .25
# 0.1726473, k = 2, base = 1, power = .25
# 0.1265414, k = 100, base = 0.95, power = .25
surv.obj <- survival::Surv(DF_valid$Y$T.OS, DF_valid$Y$C.OS)
sbrier.bart.vals <- sapply(1:nrow(df.bart$yhat.test), function(b) {
  fit <- survival::survfit(survival::coxph(surv.obj ~ df.bart$yhat.test[b,-c(1:(nrow(mm.df_neighb) + 1))]))#test.pred$yhat.test[b,]))
  class(fit) <- "survfit"
  ipred::sbrier(surv.obj, fit)
})
mean(sbrier.bart.vals)
# 0.1477806

#0.1262651 k = 2, base = 0.95, power =.25
# 0.1261279, k = 20, base = 0.95, power = .25
# 0.1262904, k = 1, base = 1, power = .25
# 0.1262436, k = 2, base = 1, power = .25
#  0.1261631, k = 2, base = 1, power = .01

#brier score
# BS <- get_survival_linear_model()$evalfit(times = DF.Y$T.OS,
#                                           event = DF.Y$C.OS,
#                                           surv.times = preDF$times,
#                                           surv = surv.array,
#                                           cens = cens.array,
#                                           method = "brier")
# print(BS$int.BS$mean)
# 0.1104898
# 0.1171551
# 0.1087906, k = 2, power = 0.25
# 0.1181439, k = 20, power = .25
# 0.1116987, k = 1, base = 1, power = .25
# 0.1087082, k = 2, base = 1, power = .25
# 0.1219672, k = 100, base = .95, power = .25


#brier aft
if(!file.exists(file.path(out.dir, "gbm_brier_score.rds"))){
  n_valid <- nrow(DF_valid$X)
  surv.bart <- cens.bart <- array(NA,
                                    dim = c(df.bart$K,
                                            nrow(df.bart$yhat.test[,-c(1:(nrow(mm.df_neighb) + 1))]),
                                            n_valid),
                                    dimnames = list(
                                      times = 1:df.bart$K,
                                      iterations=1:nrow(df.bart$yhat.test[,-c(1:(nrow(mm.df_neighb) + 1))]),
                                      observations = 1:n_valid)
  )
  for(i in 1:n_valid) {
    surv.bart[,,i] <- t(df.bart$surv.test[, (i-1 +(nrow(mm.df_neighb) + 1)) * df.bart$K + 1:df.bart$K])
    cens.bart[,,i] <- t(df.cens.bart$surv.test[, (i-1 + (nrow(mm.df_neighb) + 1)) * df.bart$K + 1:df.bart$K])
  }
  BS <- get_survival_linear_model()$evalfit(times = DF_valid$Y$T.OS,
                                            event = DF_valid$Y$C.OS,
                                            surv.times = df.bart$times,
                                            surv = surv.bart,
                                            cens = cens.bart,
                                            method = "brier")
  print(BS$int.BS$mean)
  # 0.1169188
  saveRDS(BS, file = file.path(out.dir, "gbm_brier_score.rds"))

} else {
  BS <- readRDS(file.path(out.dir, "gbm_brier_score.rds"))
}
# gam brier
# BS_gam <- get_survival_linear_model()$evalfit(times = DF.Y$T.OS,
#                                               event = DF.Y$C.OS,
#                                               surv.times = seq(0, max(preDF$times), length.out = 16)[-1],
#                                               surv = gam_post$test$mu$S$surv,
#                                               cens = gam_cens$test$mu$S$surv,
#                                               method = "brier")
# saveRDS(BS_gam, file = file.path(out.dir, "gbm_brier_score_gam.rds"))
#### Plots ####
# brier score over time
num.lines <- 100
# gam_times <- seq(0, max(df.bart$times), length.out = 16)[-1]
bs_lines <- data.frame(x = rep(df.bart$times, num.lines),
                       y = c(BS$brier.score$bscore[,seq(1,n.samp, length.out = num.lines)]),
                       group = 1:num.lines)
bs_ribbon <- data.frame(x = df.bart$times,
                        low = BS$brier.score$low,
                        high = BS$brier.score$high)

# bs_lines_gam <- data.frame(x = rep(gam_times, num.lines),
#                        y = c(BS_gam$brier.score$bscore[,seq(1,n.samp, length.out = num.lines)]),
#                        group = 1:num.lines)
# bs_ribbon_gam <- data.frame(x = gam_times,
#                         low = BS_gam$brier.score$low,
#                         high = BS_gam$brier.score$high)

pdf(file = "inst/figure/applied/GBM/bs_time_gbm.pdf", width = 4, height = 4)
ggplot() +
  # geom_line(data = bs_lines, mapping = aes(x = x, y = y, group = group) ,
  #           color = "gray", alpha = 0.5) +
  geom_ribbon(data = bs_ribbon,  aes(x = x, ymin = low, ymax = high),
              fill = "gray", alpha = 0.4) +
  geom_line(aes(x = df.bart$times, y = BS$brier.score$mean),
            color = "blue", size = 1) +
  theme_bw() +
  xlab("Time (days)") +
  ylab("Brier Score") +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.5)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,max(df.bart$times)*1.05))
dev.off()

# pdf(file = "inst/figure/applied/GBM/bs_time_gbm_gam.pdf", width = 4, height = 4)
# ggplot() +
#   # geom_line(data = bs_lines_gam, mapping = aes(x = x, y = y, group = group) ,
#   #           color = "gray", alpha = 0.5) +
#   geom_ribbon(data = bs_ribbon_gam,  aes(x = x, ymin = low, ymax = high),
#               fill = "gray", alpha = 0.4) +
#   geom_line(aes(x = gam_times, y = BS_gam$brier.score$mean),
#             color = "blue", size = 1) +
#   theme_bw() +
#   xlab("Time (days)") +
#   ylab("Brier Score") +
#   scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0,max(preDF$times)*1.05))
# dev.off()

# integrated brier score
pdf(file = "inst/figure/applied/GBM/intbs_gbm.pdf", width = 4, height = 4)
ggplot(data = data.frame(x=BS$int.BS$intBS), aes(x = x)) +
  geom_histogram(fill = "blue", alpha = 0.5, binwidth = 0.005) +
  geom_vline(xintercept=BS$int.BS$mean) +
  theme_bw() +
  ylab("Counts") +
  xlab("Integrated Brier Score") +
  # scale_y_continuous(expand = c(0,0), limits = c(0,350)) +
  # scale_x_continuous(expand = c(0,0), limits = c(0,0.79))
  scale_y_continuous(expand = c(0,0), limits = c(0,400)) +
  scale_x_continuous(expand = c(0,0), limits = c(0,.2))
dev.off()

# pdf(file = "inst/figure/applied/GBM/intbs_gbm_gam.pdf", width = 4, height = 4)
# ggplot(data = data.frame(x=BS_gam$int.BS$intBS), aes(x = x)) +
#   geom_histogram(fill = "blue", alpha = 0.5, binwidth = 0.005) +
#   theme_bw() +
#   ylab("Counts") +
#   xlab("Integrated Brier Score") +
#   scale_y_continuous(expand = c(0,0), limits = c(0,350)) +
#   scale_x_continuous(expand = c(0,0), limits = c(0,0.79))
# dev.off()




