#### poisson data constructor ####
#' Constructs Poisson data for use in Survival models
#'
#' @param time observation times
#' @param event.times times at which events coccur
#' @param event indicator for censored (0) or event (1)
#' @param x covariates
#'
#' @return
#' @export
pois_dat <- function(time, event.times, event, x) {
  unique.times <- sort(unique(time))
  if(!(0 %in% unique.times)) unique.times <- c(0,unique.times)

  X <- NULL
  id <- NULL
  y <- y.temp <- time.temp <- offset <- time.save <- NULL
  for(i in seq_along(unique.times[-1])) {
    tt <- unique.times[i]
    which.idx <- (event.times > tt)
    id <- c(id, which(which.idx))
    time.temp <- event.times[which.idx]
    event.temp <- event[which.idx]
    X <- rbind(X, x[which.idx,])
    y.temp <- rep(0, sum(which.idx))
    y.temp[time.temp <= unique.times[i + 1] & event.temp == 1] <- 1
    y <- c(y, y.temp)
    offset.temp <- rep(unique.times[i + 1] - unique.times[i],
                       sum(which.idx))
    offset.temp[y.temp == 1] <- (time.temp - unique.times[i])[y.temp==1]
    offset <- c(offset, offset.temp)
    time.save <- c(time.save, time.temp)
  }
  cum.time.list <- tapply(offset, id, cumsum)
  cum.time<- rep(NA, length(y))
  for(i in id) {
    cum.time[id == i] <- cum.time.list[[i]]
  }
  return(list(x = X, y = y, offset = offset, time = cum.time, event.time = time.save, id = id))
}


gam_dat <- function(times, event.times, event, x, n.time.groups = 15) {
  poisdf <- pois_dat(time = seq(0,max(times), length.out=n.time.groups),
                     event.times = event.times,
                     event = event,
                     x = x)
  gam.df <- data.frame(y = poisdf$y, time = poisdf$time,
                       offset = (poisdf$offset), poisdf$x,
                       event.time = factor(poisdf$time),
                       id = poisdf$id)
  return(gam.df)
}
