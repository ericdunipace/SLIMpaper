rm(list=ls())
set.seed(3696797) #from random.org

commands   <- commandArgs(trailingOnly = TRUE)

count      <- 1
families    <- corrs <- numobs <- numcoefs <- penalty_factors <- NULL
for(i in seq_along(commands)){
  if(commands[i] == "cc"){
    count  <- count + 1
    next
  }
  if(count == 1) families <- c(families, commands[i])
  if(count == 2) corrs <- c(corrs, commands[i])
  if(count == 3) numobs <- c(numobs, commands[i])
  if(count == 4) numcoefs <- c(numcoefs, commands[i])
  if(count == 5) penalty_factors <- c(penalty_factors, commands[i])
  # if(count == 5) x_methods <- c(x_methods, commands[i])
}


n.fam      <- length(families)
n.obs      <- length(numobs)
n.coef     <- length(numcoefs)
n.corr     <- length(corrs)
n.exper    <- as.numeric(Sys.getenv("NUMEXPERIMENT"))
n.pfs      <- length(penalty_factors)
# n.xmeth    <- length(x_methods)

n.total    <- prod(c(n.fam, n.obs, n.coef, n.exper, n.corr))
seeds      <- sample.int(.Machine$integer.max, n.total)

seed_array <- array(seeds, dim=c(n.fam, n.corr, n.obs, n.coef, n.exper),
                    dimnames = list(families = families,
                                    corr = corrs,
                                    n.obs = numobs,
                                    n.coefs = numcoefs,
                                    exper.num = 1:n.exper)
                    )
dump("seed_array", file="seeds.Rdmped")

q("no")

#  R --no-save --no-restore "--args ${methods[*]} cc ${numobs[*]} cc ${numcoef[*]} cc ${penalty_factors[*]} cc ${x_methods[*]}" seed.R Output/seed.txt
