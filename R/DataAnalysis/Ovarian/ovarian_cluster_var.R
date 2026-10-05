rm(list=ls())

#### Set seed for reproducibility ####
arraynum <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
get.seed <- readRDS("ovar_seeds.rds")
set.seed(get.seed[arraynum])


#### Require Packages ####
require(SLIMpaper)

#### Load Data ####
ovar <- readRDS(file = "../Data/Ovarian/tcga_ovar.rds")

X <- ovar$recurr$X
time <- ovar$recurr$time
event <- ovar$recurr$event

#### Cox Hyperparameters ####
target <- get_survival_linear_model()
mu_0 <- 0
sigma_0 <- 1
hyperparameters <- list(mu = mu_0,
                        sigma = sigma_0)

#### Sampling Hyperparameters ####
n <- nrow(X)
nsamp <- 1e5
method <- "bvs"
nchain <- 1

#### Var sel ####
sel <- target$rpost(n.samp = nsamp, x = X, y = time, fail = event,
                    hyperparameters = hyperparameters,
                    id = 1:n, nchain = nchain,
                    method = method, seed = sample.int(.Machine$integer.max,1),
                    parallel = FALSE)

date <- gsub(" ", "_", as.name(as.character(Sys.time())))
date <- gsub(":", "=", date)
term <- paste0(c(date, ".rds"), collapse="")
varfn <- paste0(c("ovarsel",arraynum,term),collapse="_")
filename <- file.path("..", "Output", "Ovarian","Selection",varfn)
saveRDS(sel, filename)
q("no")
