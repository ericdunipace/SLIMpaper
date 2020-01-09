brierScore <- function(y, event, times,
                       surivalProb,
                       censoringProb,
                       returnRaw = TRUE) {
  #setup dimensions
  n <- length(y)
  nt <- length(times)
  dimS <- dim(survivalProb)
  dimC <- dim(censoringProb)
  S <- dimS[2]
  Sc <- dimC[2]

  #checks
  stopifnot(length(event) == n)
  stopifnot(dimS[1] == dimC[1])

  #sort times
  times <- sort(times)

  #generate things for C++ code
  survList <- lapply(1:nt, function(i) matrix(NA, nrow=n, ncol=S))
  censList <- lapply(1:nt, function(i) matrix(NA, nrow=Sc, ncol=n))
  censIndiv<- matrix(NA, nrow = Sc, ncol = n)

  for(i in 1:nt) {
    survList[[i]] <- t(survivalProb[i,,])
    censList[[i]] <- censoringProb[i,,]
  }
  time.idx <- as.numeric(cut(y, times, include.lowest = TRUE))
  for(i in 1:n) censIndiv[,i] <- censoringProb[time.idx[i],,i]

  bs <- brierScore_(as.double(y), as.integer(event),
                    as.double(times), survList, censList,
                    censIndiv,
                    as.integer(S),
                    as.integer(Sc))


  intbs <- (diff(c(0,times)) %*% bs)/max(times)


  output <- list(
    brier.score = list(
      mean = rowMeans(bs),
      low = apply(bs, 1, quantile, 0.025),
      high = apply(bs,1, quantile, 0.975),
      median = apply(bs,1, median),
      bscore = NULL
    ),
    int.BS = list(
      mean = mean(intbs),
      low = quantile(intbs, 0.025),
      high = quantile(intbs, 0.975),
      median = median(intbs),
      intBS = NULL)
  )
  if(returnRaw) {
    output$brier.score$bscore <- bs
    output$int.BS$intBS <-  c(intbs)
  }

  return(output)

}
