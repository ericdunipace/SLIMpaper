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


#### Load Data ####
data.file   <- "GBM_DF.RData"
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
interpY_full     <- readRDS(file = file.path("gbm_model_Y.rds"))
interpY_neighb   <- interpY_full[seq(1,by=1, length.out = K * n_neighb),]
interpY_train    <- readRDS(file = file.path("gbm_model_Y_training.rds"))
interpY_valid     <- interpY_full[-seq(1,by=1, length.out = K * (n_neighb + 1)),]

# interpX     <- construct_interp_survival_matrix(mm.df_neighb, K, intercept=FALSE)
nS          <- ncol(interpY_neighb)
nT          <- preDF$K
dat_neighb  <- gamm_interp_data_gbm(preDF$tx.test[,-1], times = preDF$times)
# dat_valid   <- gamm_interp_data_gbm(df_valid, times = preDF$times)
dat_train   <- gamm_interp_data_gbm(df_train, times = preDF$times)




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
# # xty         <- Matrix::crossprod(sparseX, interpY_neighb)
# new_times   <- runif(15)
# std_times   <- c(preDF$times/max(preDF$times), new_times)
# cost        <- cost_calc(t(std_times), t(std_times), 2.0)^2
# x_smooth    <- cbind(1, matrix(preDF$times))
# f_groups    <- rep(groups,nS)
# pf          <- rep(1, nT * ncol(sparseX))
# lambda_f    <- 1
# lambda_g    <- group_lambda_zero(sparseX, interpY_neighb, groups, sqrt(tapply(rep(1,length(f_groups)), f_groups, sum)))

#### Formula ####
gam_form  <- formula(y ~ Resection +
                       # Age + Gender + KPS +  MGMT +#Age:Resection +
                       s(time, by = Age) +
                       s(time, by = Gender) +
                       s(time, by = KPS) +
                       s(time, by = MGMT) +
                       s(time, by = Resection, pc = 0) - 1)

#### Global model ####
global_file <- file.path(out.dir, paste0("global_model_",arraynum, ".RDS"))

global <- list()

# global$L1    <- gam_iterate(gam_form, y = interpY_train[,arraynum,drop=FALSE],
#                             x = dat_train$gammX, extract = dat_train$extract_terms,
#                             time = dat_train$times, nT=nT,
#                             which.gam = "qgam"
#                             , lsig = -1.606527, err = 0.15
# )
# saveRDS(global, global_file)

global$L2  <- gam_iterate(gam_form, y = interpY_train[,arraynum,drop=FALSE],
                          x = dat_train$gammX, extract = dat_train$extract_terms,
                          time = dat_train$times, nT=nT,
                          which.gam = "gam")

saveRDS(global, global_file)

#### validation model ####
# valid_file <- file.path(out.dir, paste0("validation_model_",arraynum, ".RDS"))
#
# valid <- list()
#
# valid$L1    <- gam_iterate(gam_form, y = interpY_valid[,arraynum,drop=FALSE],
#                             x = dat_valid$gammX, extract = dat_valid$extract_terms,
# time = dat_valid$times,, nT=nT,
#                             which.gam = "qgam", lsig = -1.606527, err = 0.15)
# saveRDS(valid, valid_file)
#
# valid$L2  <- gam_iterate(gam_form, y = interpY_valid[,arraynum,drop=FALSE],
#                           x = dat_valid$gammX,
#                           extract = dat_valid$extract_terms,
#                           time = dat_valid$time, nT=nT,
#                           which.gam = "gam")
#
# saveRDS(valid, valid_file)

#### Local model ####
local_file <- file.path(out.dir, paste0("local_model_",arraynum, ".RDS"))
local <- list()

# local$L1    <- gam_iterate(gam_form, y = interpY_neighb[,arraynum,drop=FALSE],
#                            x = dat_neighb$gammX,
#                            extract = dat_neighb$extract_terms,
#                            time = dat_neighb$time, nT=nT,
#                            which.gam = "qgam"
#                            , lsig = -1.606527, err = 0.15
# )
# saveRDS(local, local_file)

local$L2  <- gam_iterate(gam_form, y = interpY_neighb[,arraynum,drop=FALSE],
                         x = dat_neighb$gammX,
                         extract = dat_neighb$extract_terms,
                         time = dat_neighb$time, nT=nT,
                         which.gam = "gam")

saveRDS(local, local_file)

q("no")
