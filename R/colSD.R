colSD <- function(x){
  sds <- rep(NA, ncol(x))
  for(i in 1:ncol(x)){
    sds[i] <- sd(x[,i], na.rm=TRUE)
  }
  return(sds)
}
