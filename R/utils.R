augPseudo <- function(X, mu, theta, theta_norm, pseudo.obs, n, same=TRUE) {
  augDat <- augmentMatrix(X, mu, theta, same=same)
  augDat$XtX <- n/(n + pseudo.obs) * augDat$XtX + diag(theta_norm) * pseudo.obs/(n + pseudo.obs)
  augDat$XtY <- n/(n + pseudo.obs) * augDat$XtY + theta_norm * pseudo.obs/(n + pseudo.obs)
  return(augDat)
}

set_penalty_factor <- function(x, method){
  distances <- switch(method,
                      expectation = abs(colMeans(x)),
                      distance = sqrt(apply(x,2,crossprod)),
                      none = rep(1, ncol(x)),
                      wt_expect = abs(colMeans(x)/colSD(x))
  )
  penalties <- 1/distances
  penalties[1] <- 0
  return(penalties)
}

calc.lambdas <- function(aug, lambda.min.ratio, penfact, n.lambda) {
  xty <- aug$XtY
  xtx <- aug$XtX

  temp.xty <- abs(xty)/penfact
  temp.xty[penfact==0] <- 0
  max.lamb <- max(temp.xty)
  lambdas <- c(exp(seq(log(max.lamb), log(max.lamb) + log(lambda.min.ratio), length.out=n.lambda)),0)
  return(lambdas)
}

job.fun <- function(names, jobid) {
  return(grep(jobid,pattern))
}

date.fun <- function(names, date.start, date.end=NULL) {
  fn <- strsplit(names, c(".rds|_"))
  times <- sapply(fn, function(f) f[length(f)])
  dates <- sapply(fn, function(f) f[length(f)-1])
  times <- gsub("=",":", times)
  fulldates <- as.POSIXct(paste(dates,times), format="%Y-%m-%d %H:%M:%S")
  if (is.null(date.end)) {
    date.end <- Sys.time()
  }

  idx <- which(fulldates > date.start & fulldates < date.end)

  return(idx)
}

plot_time <- function(times, ylabs=NULL, ylim=NULL, alpha=NULL, base_size = 11, scale.time = 1) {
  if(is.null(alpha)) alpha <- 0.3
  if(is.null(base_size)) base_size <- 11



  times$method <- factor(times$method)
  # methods <- levels(times$method)
  n <- unique(times$n)
  p <- unique(times$p)

  plots <- vector("list", length(n) + length(p))
  names(plots) <- paste(c(n,p))

  if(length(n) ==1 & length(p)==1) {
    plots[[2]] <- NULL
    names(plots) <- paste(n,p)
  }
  # times <- times[times$method != "Simulated Annealing",]
  chartimes <- as.character(times$time)
  censored <- chartimes[grep(">", chartimes)]
  censored <- gsub(">","",censored)
  min.cens <- min(as.numeric(censored))
  times$time <- as.numeric(chartimes)
  times$time[is.na(times$time) & times$method == "Simulated Annealing"] <- min.cens
  times$time <- times$time/scale.time

  E <- tapply(times$time, INDEX = times$method, mean)
  hi <- tapply(times$time, INDEX = times$method, quantile, 0.975, na.rm=TRUE)
  low <- tapply(times$time, INDEX = times$method, quantile, 0.025, na.rm=TRUE)

  if(!is.matrix(E)) E<- t(as.matrix(E))
  if(!is.matrix(hi)) hi<- t(as.matrix(hi))
  if(!is.matrix(low)) low<- t(as.matrix(low))

  df <- data.frame(time = c(E), low = c(low), hi = c(hi),
                   groups = rep(colnames(E), each = nrow(E)),
                   row.names = NULL)
  df <- df[complete.cases(df),]
  # ylim <- c(0, max(hi[,colnames(hi) != "Simulated Annealing"], na.rm=TRUE) * 1.1)
  ylim <- c(0, max(hi, na.rm=TRUE) * 1.1)
  plots[[1]] <- ggplot2::ggplot( df,
                            ggplot2::aes(x=groups, y=time,
                                         color = groups, fill = groups,
                                         group=groups )) +
    ggplot2::geom_col(alpha = alpha) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = low, ymax = hi), width = 0.7, size = 1) +
    ggsci::scale_color_jama() +
    ggsci::scale_fill_jama() +
    ggplot2::labs(fill ="Method", color="Method") +
    ggplot2::ylab(ylabs[1]) + ggplot2::theme_bw(base_size) +
    # ggplot2::theme(legend.position = "none") +
    # ggplot2::scale_x_discrete(expand = c(0, 0)) /+
    # ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.9,
    #                                  hjust = 0.9)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank()) +
    ggplot2::xlab("") +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = ylim )
  return(plots)
}
