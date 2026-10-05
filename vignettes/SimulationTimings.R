## This file will rerun the experiments for the paper
## "Interpretable Summaries"


rm(list=ls())
# Load packages
library(SLIMpaper)
library(doRNG)
library(doParallel)
library(parallel)


#### Set Conditions ####
arraynumbers <- 1:100

p.exp <- 1L:5L
Ps <- c(17, 31, 41, 51, 61, 71, 81, 91, 10L^2L * p.exp +1L)
Ts <- as.integer(c(10,100,1000,2000,4000))
Ns <- 2L^(seq(6L,16L, 2L))

penalty_factor <- "none"
transport.method <- "exact"
wp_dist_alg <- "exact"
only.timing <- FALSE
solver <- "mosek" # can also use "cone", which is free and uses the ECOS solver
recalculate <- FALSE
run_L0 <- FALSE


#x structure
corr.x <- 0.0

#model conditions
# penalties <- source("penalties.Rdmped")
penalty_type <- "mcp.net" # can change to preferred penalty
penalty <- match.arg(as.character(penalty_type),WpProj::L1_penalty_options())
lambda.min.ratio <- as.numeric(1e-4)

family <- "gaussian"

L0 <- as.logical(FALSE)
calc_w2_post_pre <- as.logical(TRUE)
only.timing <- TRUE
not.only.timing <- !(as.logical(only.timing))

penalty.factor <- "none"

# initialize variables used in loop
target <- NULL
output <- NULL
seed   <- NULL
date   <- NULL
term   <- NULL
spfn   <- NULL
out    <- NULL

#### Setup condition list ####
conditions <- list(family = family)
conditions$penalty <- penalty
conditions$lambda.min.ratio <- lambda.min.ratio
conditions$n <- list(NULL)
conditions$p <- list(NULL)
conditions$L0 <- L0
conditions$calc_w2_post <- calc_w2_post_pre
conditions$wp_alg <- wp_dist_alg
conditions$not.only.timing <- not.only.timing
conditions$n.lambda <- 21L # so has same number of evaluations as BP method
conditions$penalty.factor <- penalty.factor
conditions$posterior.method <- "conjugate"
conditions$transport.method <- transport.method
conditions$stan_dir <- list(NULL)
conditions$python.path <- list(NULL)
conditions$recalculate <- isTRUE(as.logical(recalculate))
conditions$n.samps <- list(NULL)
conditions$solver <- solver
conditions$methods.to.run <- c("approximate binary program","binary program","lagrange binary program")

target <- get_normal_linear_model()
target$X$corr <- corr.x

#priors
prior.sigma <- c(1,1)

conditions$p <- 17L # will be 21 after adding nonlinear terms
conditions$n.samps  <- 100L
conditions$n <- 1024

hyperparameters <- list(mu = NULL, sigma = NULL,
                        alpha = NULL, beta = NULL,
                        Lambda = NULL)
hyperparameters$mu <- rep(0, conditions$p )
hyperparameters$sigma <- diag(1, conditions$p , conditions$p )
hyperparameters$alpha <- as.numeric(prior.sigma[1])
hyperparameters$beta <- as.numeric(prior.sigma[2])
hyperparameters$Lambda <- solve( diag(1, conditions$p , conditions$p ))

# experiment function so don't have to repeat
exper_and_save <- function(target, hyperparameters, conditions,arraynum) {
  output <- SLIMpaper::experimentWPMethod(target, hyperparameters, conditions)
  names(output$time$selection) <- c("binary program", "L.R. B.P.","relaxed B.P.","sel HC", "sel stepwise",
                                    "sel anneal")
  names(output$time$projection) <- c("W2","W1","WInf", "HC", "stepwise",
                                     "anneal")
  ptem <- ifelse(conditions$p < 30, conditions$p+4,conditions$p)
  sel_timings <- unlist(output$time$selection)
  proj_timings <- unlist(output$time$projection)
  comb_timings <- c(sel_timings, proj_timings)
  timings <- data.frame(time = comb_timings,
                        n = conditions$n,
                        p = ptem,
                        n.samp = conditions$n.samps,
                        method =  names(comb_timings)
  )

  #### Save File ####
  # p <- as.numeric(n.coef)
  date <- gsub(" ", "_", as.name(as.character(Sys.time())))
  date <- gsub(":", "=", date)
  term <- paste0(c(date, ".rds"), collapse="")
  spfn <- paste0(c("SP",family,conditions$transport.method,"Corr",target$X$corr,conditions$n,ptem,arraynum,term),collapse="_")
  spdir <- file.path("Output_timing", conditions$family, conditions$penalty,conditions$penalty.factor,
                     conditions$transport.method, paste0("Corr_",corr.x),conditions$n,ptem)
  svfn <- file.path(spdir, spfn)
  if(!dir.exists(spdir)) dir.create(spdir,recursive = TRUE)
  if(!is.null(output)) saveRDS(timings, file=svfn)
}
export <- c("target","conditions",
            "hyperparameters", "exper_and_save")

# run experiments
# set up parallel clusters
#### Change N ####
conditions$methods.to.run <- "binary program"
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(796843192) # from Random.org
for(n in Ns) {
  print(sprintf("n = %i", n))
  conditions$n <- n

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj"),
                          .errorhandling = "pass") %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)

    }
}
parallel::stopCluster(cl)

conditions$methods.to.run <- "lagrange binary program"
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(796843192) # from Random.org
for(n in Ns) {
  print(sprintf("n = %i", n))
  conditions$n <- n

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj"),
                          .errorhandling = "pass") %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)

    }
}
parallel::stopCluster(cl)

conditions$methods.to.run <- "approximate binary program"
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,8L))
doParallel::registerDoParallel(cl)
set.seed(796843192) # from Random.org
for(n in Ns) {
  print(sprintf("n = %i", n))
  conditions$n <- n

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj"),
                          .errorhandling = "pass") %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)

    }
}
parallel::stopCluster(cl)

#### Change T ####
conditions$methods.to.run <- "binary program"
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(32580726) # from Random.org
conditions$n <- 2L^10L
conditions$p <- 17L
for(n.samps in Ts) {
  print(sprintf("T = %i", n.samps))
  conditions$n.samps <- n.samps

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj")) %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)
    }
}
parallel::stopCluster(cl)

conditions$methods.to.run <- "lagrange binary program"
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(32580726) # from Random.org
conditions$n <- 2L^10L
conditions$p <- 17L
for(n.samps in Ts) {
  print(sprintf("T = %i", n.samps))
  conditions$n.samps <- n.samps

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj")) %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)
    }
}
parallel::stopCluster(cl)

conditions$methods.to.run <- "approximate binary program"
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(32580726) # from Random.org
conditions$n <- 2L^10L
conditions$p <- 17L
for(n.samps in Ts) {
  print(sprintf("T = %i", n.samps))
  conditions$n.samps <- n.samps

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj")) %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)
    }
}
parallel::stopCluster(cl)

#### Change P ####
#full bp
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(943408348) # from Random.org
conditions$n <- 2L^10L
conditions$n.samps  <- 100L
conditions$methods.to.run <- c("binary program")
for(p in Ps) {
  print(sprintf("p = %i", p))
  conditions$p <- p

  if (p > 59) {
    break
  }

  #priors
  p_param <- p
  hyperparameters$mu <- rep(0, p_param)
  hyperparameters$sigma <- diag(1, p_param, p_param)
  hyperparameters$alpha <- as.numeric(prior.sigma[1])
  hyperparameters$beta <- as.numeric(prior.sigma[2])
  hyperparameters$Lambda <- solve( diag(1, p_param, p_param))

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj")) %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)
    }
}
parallel::stopCluster(cl)

#lagrange BP
cl <- parallel::makeCluster(min(parallel::detectCores()-1L,5L))
doParallel::registerDoParallel(cl)
set.seed(943408348) # from Random.org
conditions$n <- 2L^10L
conditions$n.samps  <- 100L
conditions$methods.to.run <- c("lagrange binary program")
for(p in Ps) {
  print(sprintf("p = %i", p))
  conditions$p <- p

  if (p > 89) {
    break
  }

  #priors
  p_param <- p
  hyperparameters$mu <- rep(0, p_param)
  hyperparameters$sigma <- diag(1, p_param, p_param)
  hyperparameters$alpha <- as.numeric(prior.sigma[1])
  hyperparameters$beta <- as.numeric(prior.sigma[2])
  hyperparameters$Lambda <- solve( diag(1, p_param, p_param))

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj")) %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)
    }
}
parallel::stopCluster(cl)

# approx bp
cl <- parallel::makeCluster(min(parallel::detectCores()-1,5L))
doParallel::registerDoParallel(cl)
set.seed(943408348) # from Random.org
conditions$n <- 2L^10L
conditions$n.samps  <- 100L
conditions$methods.to.run <- "approximate binary program"
for(p in Ps) {
  print(sprintf("p = %i", p))
  conditions$p <- p

  #priors
  p_param <- p
  hyperparameters$mu <- rep(0, p_param)
  hyperparameters$sigma <- diag(1, p_param, p_param)
  hyperparameters$alpha <- as.numeric(prior.sigma[1])
  hyperparameters$beta <- as.numeric(prior.sigma[2])
  hyperparameters$Lambda <- solve( diag(1, p_param, p_param))

  out <- foreach::foreach(arraynum = arraynumbers,
                          .export = export,
                          .packages = c("SLIMpaper",
                                        "WpProj")) %dorng%
    {
      exper_and_save(target, hyperparameters, conditions, arraynum)
    }
}
parallel::stopCluster(cl)

#### stop ####
warnings()
