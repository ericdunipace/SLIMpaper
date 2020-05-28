gam_iterate <- function(formula, y, x, extract, time, nT) {

  covar_vals_raw <- list()
  covar_vals <- array(NA, dim = c(nT, 7, ncol(y)))
  # foreach(s = 1:ncol(y)) %doPar% {
  # }
  for(s in 1:ncol(y)) {
    gamm_df   <- data.frame(y = y[,s], x, time = time)

    gamFit   <- gam(formula, data = gamm_df)
    covar_vals_raw[[1]] <- predict(gamFit, newdata = extract, type = "terms")
    covar_vals[,,s] <- covar_vals_raw[[1]][1:nT,2:8]
    covar_vals[,5,s] <- covar_vals[,5,s] + covar_vals_raw[[1]][1:nT,1]
    covar_vals[,6,s] <- rowSums(covar_vals_raw[[1]][extract$Resection=="Sub",])
    covar_vals[,7,s] <- rowSums(covar_vals_raw[[1]][extract$Resection=="Gross",])
  }
  return(covar_vals)
}
