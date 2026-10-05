# simple_adapt_subset <- readRDS("examples/binary/simple_adapt_subset.RDS")
# n <- 2^10
# xtx <- simple_adapt_subset$selection$xtx
# xty <- simple_adapt_subset$selection$xty
#
# max.lamb <- max(xty)
# beta_temp <- rep(0, ncol(xtx))
# beta_temp[which.max(xty)] <- 1
# candidates <- xty - xtx %*% beta_temp
#
# max(candidates)
# print(candidates)
require(SparsePosterior)
require(CoarsePosteriorSummary)
require(oem)
require(glmnet)

lambdas <- simple_adapt$selection$lambda
lambda1 <- lambdas[1]
suffstat <- sufficientStatistics(X = X, Y = t(post.bart$eta), post.stan$theta, FALSE, 0.0, "selection.variable")


beta_old <- beta <- rep(0, ncol(xtx))
beta_store <- beta_store_sel <- beta_store_sel_sort <- matrix(0, nrow=length(beta),ncol=length(lambdas))


xtx <- suffstat$XtX
xty <- suffstat$XtY
xty_update <- xty
xtx_spec <- xtx
diag(xtx_spec) <- 0
xtx_spec <- crossprod(X)/n
diag(xtx_spec) <- 0
Ytemp <- X %*% c(1:21)/100 + rnorm(n)
xty_update <- crossprod(X,Ytemp)/n
for(ll in seq_along(lambdas)){
  lambda <- lambdas[ll]
  for(i in 1:1000) {
    for(k in 1:ncol(xtx_spec)) {
      u <- xty_update[k,drop=FALSE] - xtx_spec[k,,drop=FALSE] %*% beta
      beta[k] <- pmax(0, abs(u) - lambda) * sign(u)
    }
    if(SparsePosterior::not.converged(beta, beta_old, 1e-7)){
      beta_old <- beta
    } else {
      break
    }
  }
  beta_store[,ll] <- beta
}

beta <- beta_old <- rep(0, ncol(X))
for(ll in seq_along(lambdas)){
  lambda <- lambdas[ll]
  for(i in 1:1000) {
    for(k in 1:ncol(X)) {
      u <- crossprod(X[,k,drop=FALSE], Ytemp - X[,-k,drop=FALSE] %*% beta[-k])/n
      beta[k] <- pmax(0, abs(u) - lambda) * sign(u)
    }
    if(SparsePosterior::not.converged(beta, beta_old, 1e-7)){
      beta_old <- beta
    } else {
      break
    }
  }
  beta_store[,ll] <- beta
}

beta <- rep(0, ncol(xtx))
xty_update <- xty
xtx_spec <- xtx
diag(xtx_spec) <- 0
for(ll in seq_along(lambdas)){
  lambda <- lambdas[ll]
  for(i in 1:1000) {
    for (k in 1:ncol(xtx)) {
      u <- xty_update[k] - xtx_spec[k,,drop=FALSE] %*% beta
      beta[k] <- pmax(0, abs(u) - lambda) * sign(u)
      beta[k] <- ifelse(beta[k] >0,1,0)
    }
  }
  beta_store_sel[,ll] <- beta
}

beta <- beta_old <- rep(0, ncol(xtx))
xty_update <- xty
xtx_spec <- xtx
diag(xtx_spec) <- 0
for(ll in seq_along(lambdas)){
  lambda <- lambdas[ll]
  for(i in 1:10000) {
    for (k in 1:ncol(xtx)) {
      u <- xty_update[k] - xtx_spec[k,,drop=FALSE] %*% beta
      beta[k] <- pmax(0, abs(u) - lambda) * sign(u)
      # beta[k] <- ifelse(beta[k] >0,1,0)

    }
    if(SparsePosterior::not.converged(beta, beta_old, 1e-7)){
      beta_old <- beta
    } else {
      break
    }
  }
  beta_store_sel[,ll] <- beta
}

check <- oem.xtx(xtx=xtx, xty=xty, family="gaussian",penalty="ols", lambda=0,maxit=10000)
checkscale <- oem.xtx(xtx=xtx, xty=xty, family="gaussian",penalty="ols", lambda=0,maxit=10000, scale.factor = sqrt(diag(xtx)))
mine <- W2L1(X=X, Y = post.bart$eta,
             theta=post.stan$theta, penalty="ols",
             nlambda = 1, lambda.min.ratio = lambda.min.ratio,
             infimum.maxit=1, maxit = 1e4, gamma = gamma, lambda=0,
             pseudo_observations = pseudo.observations, display.progress = FALSE,
             method="scale")
cbind(check$beta[[1]], checkscale$beta[[1]], mine$beta,
      solve(xtx,xty), beta)
# scaled version does perform better
c(check$d, checkscale$d, mine$d) #eigenvalues look good
all.equal(xtx, mine$xtx)
all.equal(xty, mine$xty)

beta <- beta_old <- beta_sort <- rep(0, ncol(xtx))
xty_update <- xty
xtx_spec <- xtx
diag(xtx_spec) <- 0
for(ll in seq_along(lambdas)){
  lambda <- lambdas[ll]
  for(i in 1:50){
    for(j in 1:100) {
      for (k in 1:ncol(xtx)) {
        u <- xty_update[k] - xtx_spec[k,,drop=FALSE] %*% beta
        beta[k] <- pmax(0, abs(u) - lambda) * sign(u)
      }
      beta <- ifelse(beta > 0,1,0)
      if(SparsePosterior::not.converged(beta, beta_old, 1e-7)){
        beta_old <- beta
      } else {
        break
      }
    }
    if(SparsePosterior::not.converged(beta, beta_sort, 1e-7)){
      xty_update <- xtyUpdate(X, t(post.bart$eta), post.stan$theta, beta, 0.0, "selection.variable")
      beta_sort <- beta
    } else {
      break
    }

  }
  beta_store_sel_sort[,ll] <- beta
}


dat <- list(temp=matrix(0, n, p), xtx = matrix(0,p,p), xty = rep_len(0, p),
            mu = rep(0, n), idx_mu = rep(0, n),
            sort_y = rep(0, n))
dat$xty <- rep(0,p)
for(i in 1:n) {
  dat$temp <- post.stan$theta * matrix(X[i,,drop=FALSE], nsamp,p, byrow = TRUE)
  dat$temp2 <- post.stan$theta * matrix(beta, nsamp,p, byrow=TRUE)  * matrix(X[i,,drop=FALSE], nsamp,p, byrow = TRUE)
  dat$mu <- rowSums(dat$temp2)
  dat$idx_mu <- order(dat$mu)
  dat$sort_y <- sort(post.bart$eta[,i])
  dat$xty = dat$xty + crossprod(dat$temp[dat$idx_mu,,drop=FALSE], dat$sort_y)/(n*nsamp)
}
cbind(dat$xty, xtyUpdate(X, t(post.bart$eta), post.stan$theta, beta, 0.0, "selection.variable"), xtyUpdate(X, t(post.bart$eta), post.stan$theta, beta, 0.0, "scale"), sufficientStatistics(X = X, Y = t(post.bart$eta), post.stan$theta, FALSE, 0.0, "selection.variable")$XtY)

oem.check <- oem.xtx(xtx=xtx, xty=xty, family="gaussian",penalty="lasso", lambda=lambdas,maxit=10000)
selection = W2L1(X=X, Y = post.bart$eta,
                 theta=post.stan$theta, penalty="selection.lasso",
                 nlambda = 10, lambda.min.ratio = lambda.min.ratio,
                 infimum.maxit=1e4, maxit = 1e3, gamma = gamma, lambda=lambdas[1],
                 pseudo_observations = pseudo.observations, display.progress = FALSE,
                 penalty.factor = penalty.factor.simple, method="selection.variable")
selection2 = W2L1(X=X, Y = post.bart$eta,
                  theta=post.stan$theta, penalty="selection.lasso",
                  nlambda = 10, lambda.min.ratio = lambda.min.ratio,
                  infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                  pseudo_observations = pseudo.observations, display.progress = TRUE,
                  penalty.factor = penalty.factor.simple, method="selection.variable")
selection3 = W2L1(X=X, Y = post.bart$eta,
                  theta=post.stan$theta, penalty="selection.lasso",
                  lambda.min.ratio = lambda.min.ratio, lambda=lambdas,
                  infimum.maxit=1e4, maxit = 1e3, gamma = gamma,
                  pseudo_observations = pseudo.observations, display.progress = FALSE,
                  penalty.factor = penalty.factor.simple, method="selection.variable")
selection$beta
selection2$beta
oem.check$beta[[1]]
selection3$beta

suffstat <- sufficientStatistics(X = X, Y = t(post.bart$eta), post.stan$theta, FALSE, 0.0, "projection")

xtx <- suffstat$XtX
xty <- suffstat$XtY

mcptest = W2L1(X = X, Y = post.bart$eta,
               theta = post.stan$theta, penalty = "mcp",
               lambda.min.ratio = 1e-3,
               nlambda = 5,
               infimum.maxit = 10, maxit = 1e3, gamma = gamma,
               pseudo_observations = pseudo.observations, display.progress = FALSE,
               penalty.factor = rep(1,p), method = "projection")
mcptest1 = W2L1(X = X, Y = post.stan$eta,
                theta = post.stan$theta, penalty = "mcp",
                lambda.min.ratio = 5e-1, nlambda = 3,
                infimum.maxit = 10, maxit = 1e3, gamma = gamma,
                pseudo_observations = pseudo.observations, display.progress = FALSE,
                penalty.factor = rep(1,p), method = "projection")
mcptest1A = W2L1(X = X, Y = post.stan$eta,
                 theta = post.stan$theta, penalty = "ols",
                 lambda.min.ratio = lambda.min.ratio, nlambda = 1,
                 infimum.maxit = 1, maxit = 1e3, gamma = gamma,
                 pseudo_observations = pseudo.observations, display.progress = FALSE,
                 penalty.factor = rep(1,p), method = "projection")
mcptest2 =  W2L1(X = X, Y = post.bart$eta,
                 theta = post.stan$theta, penalty = "mcp",
                 lambda.min.ratio = lambda.min.ratio,
                 nlambda=50,
                 infimum.maxit = 5, maxit = 1e3, gamma = gamma,
                 pseudo_observations = pseudo.observations, display.progress = FALSE,
                 penalty.factor = rep(1,p), method = "projection")
mcptest3 = W2L1(X = X, Y = post.stan$eta,
                theta = post.stan$theta, penalty = "mcp",
                lambda.min.ratio = lambda.min.ratio, nlambda = 50,
                infimum.maxit = 5, maxit = 1e3, gamma = gamma,
                pseudo_observations = pseudo.observations, display.progress = FALSE,
                penalty.factor = rep(1,p), method = "projection")

mcplist <- list(sort_bart = mcptest, sort_logis =mcptest1, bart = mcptest2, logis = mcptest3, ols_logis = mcptest1A)
mcpplot <- plot.compare(mcplist, post.bart$eta, X,
                        t(post.stan$theta), "w2", "mean", TRUE)
mcpplot$plot$mean + scale_y_continuous(limits=c(0,50)) + geom_point()

mcpplot_mse <- plot.compare(mcplist, prob, X,
                            t(post.stan$theta), "mse", "mean", TRUE, plogis)
mcpplot_mse$plot$mean + scale_y_continuous(limits=c(0,.1)) + geom_point() + geom_hline(yintercept = 0.05225924)

theta <- post.stan$theta
X_std <- scale(X)
X_std[,1] <- 1
eta <- X_std %*% t(cbind(extract(post.stan$model,"intercept")[[1]],
                         extract(post.stan$model,"beta")[[1]]))
all.equal(c(t(eta)), c(post.stan$eta))
eta_test <- X %*% t(post.stan$theta)
all.equal(c(t(eta)), c(t(eta_test)))
all.equal(c(post.stan$eta), c(t(eta_test)))


f <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
    10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

mad <- function(y.train, y.train.hat)
  mean(abs(y.train - apply(y.train.hat, 1L, mean)))


require(dbarts)
## low iteration numbers to to run quickly
xval <- xbart(x, y, n.samples = 15L, n.reps = 4L, n.burn = c(10L, 3L, 1L),
              n.trees = c(5L, 7L),
              k = c(1, 2, 4),
              power = c(1.5, 2),
              base = c(0.75, 0.8, 0.95), n.threads = 1L,
              loss = mad)

xval <- xbart(X, Y, n.samples = 15L, n.reps = 4L, n.burn = c(10L, 3L, 1L),
              n.trees = c(5L, 7L),
              k = c(1, 2, 4),
              power = c(1.5, 2),
              base = c(0.75, 0.8, 0.95), n.threads = 1L,
              loss = mad)

xval <- xbart(X, Y, n.samples = 15L, n.reps = 4L, n.burn = c(10L, 3L, 1L),
              n.trees = c(5L, 7L),
              k = c(1, 2, 4),
              power = c(1.5, 2),
              base = c(0.75, 0.8, 0.95), n.threads = 1L,
              loss = "log")

true.loss <- function(ytrain, yhat) {
  mean((prob[names(ytrain)] - yhat)^2)
}

# names(Y) <- names(prob) <- 1:n
test <- xbart(X,Y, verbose = TRUE,
              method = "k-fold", n.reps = 5, n.test = 5, n.burn = c(1000,750, 250),
              n.samples = 100, drop = FALSE,
              k = c(1), power = seq(0.1, 4, length.out = 5),
              base = seq(0.1, 0.95, length.out = 5), n.trees = c(200),
              n.threads = 3, loss = "log")
means <- apply(test,2:5,mean)
arrayInd(which.min(means), dim(means), dimnames(means), useNames = TRUE)
