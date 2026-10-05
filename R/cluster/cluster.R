rm(list=ls())
# tt <- proc.time()
#### Get Commands from environment ####
  arraynum <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  jobid <- Sys.getenv('SLURM_ARRAY_JOB_ID')
  method <-  Sys.getenv('METHOD')
  # x_method <- Sys.getenv("X_METHOD")
  n.obs <- Sys.getenv("NOBS")
  n.coef <- Sys.getenv("NCOEF")
  n.samps <- Sys.getenv("NSAMPS")
  prior.sigma <- Sys.getenv("PRIOR_SIGMA")
  correlation_x <- Sys.getenv("CORR_X")
  penalty_type <- Sys.getenv("PENALTY")
  distribution_family <- Sys.getenv("FAM")
  run_L0 <- Sys.getenv("L0")
  w2_post <- Sys.getenv("W2_POST")
  penalty_factor <- Sys.getenv("PENALTY_FACTOR")
  lambda_min_ratio <- Sys.getenv("LAMBDA_MIN_RATIO")
  # pseudo.obs <- Sys.getenv("PSEUDO_OBSERVATIONS")
  posterior.method <- Sys.getenv("POSTERIOR_METHOD")
  transport.method <- Sys.getenv("TRANSPORT_METHOD")
  wp_dist_alg <- Sys.getenv("WP_DIST_ALG")
  only.timing <- Sys.getenv("TIMING")
  solver <- Sys.getenv("SOLVER")
  python.path <- Sys.getenv("PYTHON_PATH")
  recalculate <- Sys.getenv("RECALC_SINGLE")
  commands <- commandArgs(trailingOnly=T)
  count <- 1
  transport.methods <- penalty_factors <- prior.sigma <- penalties <- families <- x_methods <-stan_dir<- NULL
  for(i in seq_along(commands)){
    if(commands[i] == "cc"){
      count <- count + 1
      next
    }
    if(count == 1) transport.methods <- c(transport.methods, commands[i])
    if(count == 2) penalty_factors <- c(penalty_factors, commands[i])
    if(count == 3) prior.sigma <- c(prior.sigma, commands[i])
    if(count == 4) penalties <- c(penalties, commands[i])
    if(count == 5) families <- c(families, commands[i])
    # if(count == 6) x_methods <- c(x_methods, commands[i])
    if(count == 6) stan_dir <- c(stan_dir, commands[i])
  }

#### Set error handler ####
  options(error = quote(dump.frames(paste0("error_dump_",arraynum), TRUE)))

#### Set Method ####
  print(transport.methods)
  transport.method <- match.arg(transport.method, transport.methods)
  wp_dist_alg <- match.arg(wp_dist_alg, transport.methods)
  print(paste0("Transport method: ", transport.method))
  print(paste0("Distance method: ", wp_dist_alg))


#### packages ####
  source("packages.R", echo=TRUE)

#### load functions ####
  # source("functions.R")

#### Load Conditions ####
  n <- as.numeric(n.obs)
  p <- as.numeric(n.coef)
  n.samps <- as.numeric(n.samps)

  #priors
  alpha <- as.numeric(prior.sigma[1])
  beta <- as.numeric(prior.sigma[2])
  mu_prior <- rep(0, p)
  sigma_prior <- diag(1, p, p)

  #x structure
  corr.x <- as.numeric(correlation_x)

  #model conditions
  # penalties <- source("penalties.Rdmped")
  print(penalties)
  penalty <- match.arg(as.character(penalty_type),penalties)
  lambda.min.ratio <- as.numeric(lambda_min_ratio)

  # families <- source("families.Rdmped")
  print(families)
  if(is.null(distribution_family)){
    family <- "gaussian"
  } else {
    family <- match.arg(as.character(distribution_family), families)
  }
  L0 <- as.logical(run_L0)
  calc_w2_post <- as.logical(w2_post)
  if(family == "binomial") calc_w2_post <- FALSE
  not.only.timing <- !(as.logical(only.timing))

  penalty.factor <- match.arg(as.character(penalty_factor), penalty_factors)

#### Set seeds ####
  seed.file <- file.path("seeds.Rdmped")
  source(seed.file)
  seed <- seed_array[family, paste(corr.x), paste(n),paste(p),paste(arraynum)]
  set.seed(seed)

#### Load Target ####
  if(family == "gaussian") {
    target <- get_normal_linear_model()
  } else if (family == "exponential" | family == "survival" | family == "cox") {
    target <- get_survival_linear_model()
  } else if (family == "binomial") {
    target <- get_binary_nonlinear_model()
  }

#### Setup Target ####
  # source("gen_x.R")
  # source("gen_param.R")
  # rXvars <- generateX(x_method)
  target$X$corr <- corr.x
  # target$X$rX <- rXvars$rX
  # target$X$rXnew <- rXvars$rXnew
  # target$rparam <- rparam(meth,target)

#### Setup condition list ####
  # p <- 6
  # mu_prior <- rep(0, p)
  # sigma_prior <- diag(1, p, p)
  conditions <- list()
  conditions$family <- family
  conditions$penalty <- penalty
  conditions$lambda.min.ratio <- lambda.min.ratio
  conditions$n <- n
  conditions$p <- p
  conditions$L0 <- L0
  conditions$calc_w2_post <- calc_w2_post
  conditions$wp_alg <- wp_dist_alg
  conditions$not.only.timing <- not.only.timing
  conditions$n.lambda <- 100 # max(p*10,100) # 100 #
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

#### Run experiment ####
  # debugonce(experimentWPMethod)
  output <- NULL
  output <- experimentWPMethod(target, hyperparameters, conditions)

#### Save File ####
  # p <- as.numeric(n.coef)
  date <- gsub(" ", "_", as.name(as.character(Sys.time())))
  date <- gsub(":", "=", date)
  term <- paste0(c(date, ".rds"), collapse="")
  spfn <- paste0(c("SP",family,transport.method,"Corr",corr.x,n,p,jobid,arraynum,term),collapse="_")
  svfn <- file.path("Output", family, penalty,penalty.factor,
                    transport.method, paste0("Corr_",corr.x),n,p, spfn)
  if(!is.null(output)) saveRDS(output, file=svfn)
# print(proc.time() - tt)
  warnings()
  q("no")

#interactive
  # R --no-save "--args ${methods[*]} cc ${penalty_factors[*]} cc ${PRIOR_SIGMA[*]} cc ${penalties[*]} cc ${families[*]} cc ${x_methods[*]}"
# bash interactive
# R CMD BATCH --no-save "--args ${methods[*]} cc ${penalty_factors[*]} cc ${PRIOR_SIGMA[*]} cc ${penalties[*]} cc ${families[*]} cc ${x_methods[*]}" cluster.R
