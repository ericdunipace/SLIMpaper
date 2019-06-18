ridgePlot <- function(fit, index = 1, maxCoef = 10, scale = 1, alpha = 0.5, full = NULL, transform = function(x){x}) {
  require(ggplot2)
  require(ggridges)
  require(ggsci)

  idx <- index
  if(length(idx) > 1) {
    ridgeplot <- lapply(idx, function(i)
                        ridgePlot(fit,i, maxCoef, scale, alpha, full, transform)
                        )
  } else {
    if ( inherits(fit, "sparse-posterior") ) {
      n <- nrow(fit$eta[[1]])
      s <- ncol(fit$eta[[1]])
      whichModel <- which(fit$nzero <= maxCoef & fit$nzero>0)
      ncoef <- rep(fit$nzero[whichModel], each = ncol(eta[[1]]))
      eta <- lapply(fit$eta[whichModel], function(ee) transform(ee[idx,,drop=FALSE]))


      df_ridge <- data.frame(value = unlist(eta),
                             ncoef = factor(ncoef))
    }
    else {
      if (is.null(names(fit)) & length(fit) > 1) names(fit) <- paste0("Model ", 1:length(fit))
      if (!(all(sapply(fit, function(f) any(class(f) %in% "sparse-posterior"))))) {
        stop("Must be a fitted model from the 'SparsePosterior' package")
      }
      # if ( length(idx) > 1) stop("Can only do ridgeline plots for one observation at a time")

      n <- nrow(fit[[1]]$eta[[1]])
      s <- ncol(fit[[1]]$eta[[1]])

      whichModel <- lapply(fit, function(nn) which(nn$nzero <= maxCoef & nn$nzero>0))
      ncoef <- mapply(function(f,wm) {rep(f$nzero[wm], each=s)},
                      f = fit, wm = whichModel)

      eta <- mapply(function(f,wm) {lapply(f$eta[wm],
                                           function(e) transform(e[idx,,drop=FALSE]))},
                    f = fit, wm = whichModel)

      df_ridge <- data.frame(value = unlist(eta), ncoef = factor(unlist(ncoef)))
      df_ridge$Method <- factor(unlist(mapply(function(n,l) {rep(n,each=l*s)}, n = names(fit), l = sapply(eta, length))))
    }

    if(!is.null(full)) {
      if(is.matrix(full)) {
        if(nrow(full)>1) {
          if(nrow(full) != n) stop("'full' must be a vector of predictions for correct observation or matrix of predictions for every observation")
          if(ncol(full) != s) stop("'full' must have the same number of predictions per observation as are found in the coarse posterior version")
          full <- transform(full[idx,])
        }
      } else {
        if(length(full) != s) stop("'full' must have the same number of predictions as found in the coarse posterior version")
      }

      # ncoef_vec <- c(as.character(df_ridge$ncoef ), rep("Full", s))
      # ncoef_lev <- c(levels(df_ridge$ncoef), "Full")
      #
      # method_vec <- c(as.character(df_ridge$Method ), rep("Full", s))
      # method_lev <- c(levels(df_ridge$Method), "Full")

      ncoef_vec <- rep("Full", s)
      ncoef_lev <- "Full"

      method_vec <- rep("Full", s)
      method_lev <- "Full"
      full_df <- data.frame(value = full,
                            ncoef = factor(ncoef_vec, levels=ncoef_lev),
                            Method = factor(method_vec, levels=method_lev))
      df_ridge <- rbind(df_ridge, full_df)
    }

    # initialize ridgeplot
    ridgeplot <- ggplot(df_ridge, aes(y = factor(ncoef)))

    #is list of fits
    if(length(fit) > 1) {
      ridgeplot <- ridgeplot + geom_density_ridges(aes(x = value, fill=Method), scale=scale, alpha=alpha)
    } else { # if single fit
      ridgeplot <- ridgeplot + geom_density_ridges(aes(x = value), scale=scale, alpha=alpha)
    }

    #set themes and labels
    ridgeplot <- ridgeplot +
      theme_ridges() +
      ylab("Number of Active Coefficients") +
      xlab("")

    # remove full from legend
    if(!is.null(full)) {
      levs <- levels(df_ridge$Method)
      levs <- levs[levs != "Full"]
      cols <- c(pal_jama("default")(length(levs)), "#e41a1c")
      ridgeplot <- ridgeplot +
        scale_fill_manual(breaks=levs, values=cols)
    } else {
      ridgeplot <- ridgeplot + scale_fill_jama()
    }
  }

  return(ridgeplot)
}
