#### Set seed for reproducibility ####
set.seed(1848078034)

#### Require Packages ####
require(SLIMpaper)
require(mosek)

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

#### Generate Lambdas ####
global_l1_lambda <- max(sqrt(colSums(X^2)))
global_linf_lambda <- max(sqrt(rowSums(theta_noint^2)))
local_l1_lambda <- max(sqrt(colSums(neighb^2)))
local_linf_lambda <- max(sqrt(rowSums(theta_noint^2)))

#### Global model ####
global_file <- file.path(direct, paste0("global_model_",arraynum, ".RDS"))

  global <- list()

  global$L1    <- WPL1(X = X, Y = eta, theta = theta_noint, power = 1, model.size = model.size,
                       penalty = "mcp", solver = "mosek",
                       lambda = global_l1_lambda[arraynum], gamma = 1.5, options = list(solver_opts = list(verbose = 10)),
                       display.progress = TRUE, maxit = 1e6, solver_opts = list(verbose = TRUE))
  global$LInf  <- WPL1(X = X, Y = eta, theta = theta_noint, power = Inf, model.size = model.size,
                       penalty = "mcp", solver = "mosek",
                       lambda = global_linf_lambda[arraynum], gamma = 1.5, options = list(solver_opts = list(verbose = 10)),
                       display.progress = TRUE, maxit = 1e6)

  saveRDS(global, global_file)


#### Local model ####
local_file <- file.path(direct, paste0("local_model_",arraynum, ".RDS"))
  local <- list()

  local$L1    <- WPL1(X = neighb, Y = eta_neighb, theta = theta_noint, power = 1, model.size = model.size,
                      penalty = "mcp", display.progress = TRUE, solver = "mosek",
                      lambda = local_l1_lambda[arraynum], gamma = 1.5,
                      display.progress = TRUE, maxit = 1e6)
  local$LInf  <- WPL1(X = neighb, Y = eta_neighb, theta = theta_noint, power = Inf, model.size = model.size,
                      penalty = "mcp", display.progress = TRUE, solver = "mosek",
                      lambda = local_linf_lambda[arraynum], gamma = 1.5,
                      display.progress = TRUE, maxit = 1e6)

  saveRDS(local, local_file)

