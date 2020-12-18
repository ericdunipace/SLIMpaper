library(ggplot2)
library(gridExtra)
library(SLIMpaper)
set.seed(23987)

# set function to generate outcome and range of values
x_range <- c(-1,5)
fun <- function(x) {
  cos( pi * x ) + 0.5 * x
}
y_sd <- 0.25

# generate data used to estimate function
nfull <- 2^9
xfull <- runif(nfull, x_range[1], x_range[2])
yfull <- fun(xfull) + rnorm(nfull, 0, y_sd)

dffull <- data.frame(x = xfull, y =  yfull)

#design matrix
design.x <- function(x) {
  model.matrix(~ I(cos( pi * x)) + I(0.5 * x))
}
predict_conj <- function(obj, x) {
  eta <- x %*% obj$theta
  fit <- rowMeans(eta)
  lwr <- apply(eta, 1, quantile, 0.025)
  upr <- apply(eta, 1, quantile, 0.975)
  return(list(predict = data.frame(fit,lwr,upr), eta=eta))
}

#conjugate model
nsamp <- 1000
hyperparameters <- list(mu = rep(0,3), sigma = diag(1,3,3),
                        alpha = 1, beta = 1)
hyperparameters.linear <-  list(mu = rep(0,2), sigma = diag(1,2,2),
                                alpha = 1, beta = 1)

target <- get_normal_linear_model()
posterior <- target$rpost(n.samp = nsamp, x = design.x(xfull), y=yfull, hyperparameters = hyperparameters,
                          method = "conjugate")

predict.fit.nonlinear <- predict_conj(posterior, design.x(dffull$x))
dffull <- cbind(dffull, predict.fit.nonlinear$predict)



lmfit.full <- target$rpost(n.samp = nsamp, x = cbind(1,xfull), y=yfull,
                           hyperparameters = hyperparameters.linear,
                           method = "conjugate")
predict_global <- predict_conj(lmfit.full, cbind(1,xfull))
dffull_line <- data.frame(x=dffull$x, y=dffull$y, predict_global$predict)

# subset data
n <- 100
x <- runif(n, 2, 3)
y <- fun(x) + rnorm(n,0,y_sd)

df <- data.frame(x,y)
pred.lm <- predict_conj(posterior, design.x(df$x))

# get global linear line
predict_global_small <- predict_conj(lmfit.full, cbind(1,x))
df_global <- cbind(x=df$x, y=df$y, predict_global_small$predict)

# estimate local line from predictions
lmfit.small <- lm(pred.lm$eta~ df$x)
local.model <- list(theta = coef(lmfit.small))
predfull <- predict_conj(local.model, cbind(1,dffull$x))
dffull_smalllm <- data.frame(x = xfull, y = yfull, predfull$predict)

predict.smalllm <- predict_conj(local.model, cbind(1,df$x))
df <- cbind(df,predict.smalllm$predict)

# create graphs
filepath <- file.path("inst","figure","nonlinear_approx.pdf")

pdf(file = filepath, width = 7.5, height = 3)

p0 <- ggplot(data = dffull, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill = "gray40", alpha = 0.5) +
  geom_line(aes(y = fit), color = "blue") +
  theme_bw() +
  ylab("y") + xlab("") +
  coord_cartesian(xlim = x_range, ylim = c(min(yfull), max(yfull)))

p1 <- p0 + geom_point(alpha = 0.5, shape = 20, size = 0.5) + ggtitle("(a)")

p0gray <- ggplot(data = dffull, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill = "gray40", alpha = 0.4) +
  geom_line(aes(y = fit), color = "black", alpha = 0.2) +
  theme_bw() +
  ylab("") + xlab("x") +
  coord_cartesian(xlim = x_range, ylim = c(min(yfull), max(yfull)))

p2 <- p0gray +
  # geom_ribbon(data = df, aes(ymin = lwr, ymax = upr), fill = "gray40", alpha = 0.5) +
  # geom_line(data = df, aes(x = x, y = fit), color = "blue") +
  geom_point(data = df, aes(x=x, y = y), alpha = 0.5, shape = 20, size = 0.5)


p3 <- p2 +
  geom_ribbon(data = dffull_line, aes(ymin = lwr, ymax = upr), fill = "gray40", alpha = 0.25) +
  geom_line(data = dffull_line, aes(x = x, y = fit), color = "black", alpha = 0.5) +
  # geom_ribbon(data = df_global, aes(ymin = lwr, ymax = upr), fill = "gray40", alpha = 0.5) +
  geom_line(data = df_global, aes(x = x, y = fit), color = "red") + ggtitle("(b)")


p4 <- p2 +
  geom_line(data = dffull_smalllm, aes(x = x, y = fit),  color = "black", alpha = 0.5) +
  geom_ribbon(data = dffull_smalllm, aes(ymin = lwr, ymax = upr), fill = "gray40", alpha = 0.25) +
  geom_line(data = df, aes(x = x, y = fit), color = "red") +
  # geom_ribbon(data = df, aes(ymin = lwr, ymax = upr), fill = "gray40", alpha = 0.5) +
  ggtitle("(c)") + xlab("") + ylab("")

grid.arrange(p1, p3, p4, nrow=1)

dev.off()

