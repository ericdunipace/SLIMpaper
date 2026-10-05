## This file will rerun the experiments for the paper
## "Interpretable Summaries"


# These simulations were originally run in parallel on a cluster, hence
# the various seeds. The seeds can be generated via the "seed.R" file.

rm(list=ls())
# Load packages
library(SLIMpaper)
library(doParallel)
library(parallel)


#### Set Conditions ####
arraynumbers <- 1:100

n <- as.numeric(1024L)
p <- as.numeric(21L)
n.samps <- as.numeric(100L)

#priors
prior.sigma <- c(1,1)
alpha <- as.numeric(prior.sigma[1])
beta <- as.numeric(prior.sigma[2])
mu_prior <- rep(0, p)
sigma_prior <- diag(1, p, p)
posterior.method <- "conjugate"


penalty_factor <- "none"
transport.method <- "exact"
wp_dist_alg <- "exact"
only.timing <- FALSE
solver <- "mosek" # can also use "cone", which is free and uses the ECOS solver
recalculate <- FALSE
run_L0 <- FALSE

Sys.unsetenv("RETICULATE_PYTHON")
python.path <- NULL # put in your favorite python location, if NULL will install python 3.10.16 and required packages in virtual env
python_packages <- c("numpy==2.2.1","torch==2.5.1","scipy==1.15.1")
stan_dir <- NULL

#x structure
correlation_x <- as.numeric(c(0, 0.5, 0.9))

#model conditions
# penalties <- source("penalties.Rdmped")
penalty_type <- "mcp.net" # can change to preferred penalty
penalty <- match.arg(as.character(penalty_type),WpProj::L1_penalty_options())
lambda.min.ratio <- as.numeric(1e-4)

# families <- source("families.Rdmped")
families <- c("gaussian", "binomial")
families <- "gaussian"

L0 <- as.logical(FALSE)
calc_w2_post_pre <- as.logical(TRUE)
only.timing <- FALSE
not.only.timing <- !(as.logical(only.timing))

penalty.factor <- "none"

#### Set seeds ####
data("seed_array", package = "SLIMpaper")

#### Load Target ####


#### Setup Target ####
# source("gen_x.R")
# source("gen_param.R")
# rXvars <- generateX(x_method)
# target$X$rX <- rXvars$rX
# target$X$rXnew <- rXvars$rXnew
# target$rparam <- rparam(meth,target)

#### Setup condition list ####
# p <- 6
# mu_prior <- rep(0, p)
# sigma_prior <- diag(1, p, p)
conditions <- list(family = NULL)
conditions$penalty <- penalty
conditions$lambda.min.ratio <- lambda.min.ratio
conditions$n <- n
conditions$p <- p
conditions$L0 <- L0
conditions$calc_w2_post <- calc_w2_post_pre
conditions$wp_alg <- wp_dist_alg
conditions$not.only.timing <- not.only.timing
conditions$n.lambda <- 100L # max(p*10,100) # 100 #
conditions$penalty.factor <- penalty.factor
conditions$posterior.method <- posterior.method
conditions$transport.method <- transport.method
conditions$stan_dir <- stan_dir
conditions$python.path <- python.path
conditions$recalculate <- isTRUE(as.logical(recalculate))
# conditions$n.experiment <- n.experiment
conditions$n.samps <- n.samps
conditions$solver <- solver
hyperparameters <- list(mu = NULL, sigma = NULL,
                        alpha = NULL, beta = NULL,
                        Lambda = NULL)
hyperparameters$mu <- mu_prior
hyperparameters$sigma <- sigma_prior
hyperparameters$alpha <- alpha
hyperparameters$beta <- beta
hyperparameters$Lambda <- solve(sigma_prior)

# set up parallel clusters
cl <- parallel::makeCluster(min(parallel::detectCores()-1L, 8L))
doParallel::registerDoParallel(cl)

# initialize variables used in loop
target <- NULL
output <- NULL
seed   <- NULL
date   <- NULL
term   <- NULL
spfn   <- NULL
out    <- NULL

# run experiments
for(family in families) {
  print(family)
  conditions$family <- family
  if(family == "gaussian") {
    target <- get_normal_linear_model()
  } else if (family == "exponential" | family == "survival" | family == "cox") {
    target <- get_survival_linear_model()
  } else if (family == "binomial") {
    target <- get_binary_nonlinear_model()
  }
  if(family == "binomial") {
    conditions$calc_w2_post <- FALSE
    conditions$n <- 2^17
    conditions$posterior.method <- "nn"
  } else {
    conditions$calc_w2_post <- calc_w2_post_pre
    conditions$n <- n
    conditions$posterior.method <- posterior.method
  }

  for (corr.x in correlation_x) {
    print(corr.x)
    target$X$corr <- corr.x

    export <- c("seed_array", "target","conditions",
                "hyperparameters","n","p", "python.path",
                "python_packages")
    # for (arraynum in arraynumbers) {
    out <- foreach::foreach(arraynum = arraynumbers,
                            .export = export,
                            .packages = c("SLIMpaper",
                                          "WpProj")) %dopar%
      {
      # set seed
      seed <- seed_array[family, paste(corr.x), paste(n),paste(p),paste(arraynum)]
      set.seed(seed, kind = "default", normal.kind = "default")
      if (family == "binomial") {
        if (is.null(python.path)) {
          version <-"3.10.16"
          reticulate::install_python(version)
          if(!reticulate::virtualenv_exists("SLIM")) {
            reticulate::virtualenv_create("SLIM", python = version, packages = python_packages)
          }
          reticulate::use_virtualenv("SLIM")
          python.path <- reticulate::py_exe()
        } else {
          reticulate::use_python(python.path, require = TRUE)
          if(!reticulate::virtualenv_exists("SLIM")) {
            reticulate::virtualenv_create("SLIM", python = version, packages = python_packages)
          }
          reticulate::use_virtualenv("SLIM")
        }
        # torch$set_num_threads(1L)
        # torch$set_num_interop_threads(1L)
        conditions$python.path <- python.path
      }
      #### Run experiment ####
      # debugonce(experimentWPMethod)
      output <- SLIMpaper::experimentWPMethod(target, hyperparameters, conditions)

      #### Save File ####
      # p <- as.numeric(n.coef)
      date <- gsub(" ", "_", as.name(as.character(Sys.time())))
      date <- gsub(":", "=", date)
      term <- paste0(c(date, ".rds"), collapse="")
      spfn <- paste0(c("SP",family,conditions$transport.method,"Corr",target$X$corr,n,p,arraynum,term),collapse="_")
      spdir <- file.path("Output", conditions$family, conditions$penalty,conditions$penalty.factor,
                         conditions$transport.method, paste0("Corr_",corr.x),conditions$n,conditions$p)
      svfn <- file.path(spdir, spfn)
      if(!dir.exists(spdir)) dir.create(spdir,recursive = TRUE)
      if(!is.null(output)) saveRDS(output, file=svfn)
    }

  }


}
parallel::stopCluster(cl)

# print(proc.time() - tt)
warnings()
