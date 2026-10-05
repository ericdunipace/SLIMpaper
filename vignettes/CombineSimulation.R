require(SLIMpaper)
library(dplyr)
library(forcats)
library(ggplot2)
library(gridExtra)

family <- "binomial"
penalty <- "mcp.net"
pf <- "none"
method <- "exact"
corrs <- c("Corr_0", "Corr_0.5", "Corr_0.9") #"Corr_0.9" #
date <- "2025-02-03"
n <- 131072 #1024, 131072
p <- 21
neighb <- c("Single_point")

neighb.plot.fun <- function(family, n, p, method, pf, penalty, corrs, date, neighborhood = c("Neighborhood", "Single_point"), basesize = 12){
  neighborhood <- match.arg(neighborhood)
  neighb <- switch(neighborhood,
                   "Neighborhood" = "neighb",
                   "Single_point" = "single"
  )
  w2_list <- mse_list <- w2r2_list <-
    w1r2_list <- deriv_list <- vector("list", length(corrs))
  names(w2_list) <- names(mse_list) <-
    names(w2r2_list) <- names(w1r2_list) <-
    names(deriv_list) <- corrs
  for(corr in corrs){
    folder <- file.path("Output", family, penalty, pf, method, corr,n,p)
    files <- list.files(folder, full.names = TRUE)
    files <- files[date.fun(files, date.start = date)]
    niter <- length(files)
    cat("Niter: ", niter, ", Correlation:", corr, ", Family: ", family,"\n", sep = "")
    outputs <- lapply(files, readRDS)

    label.df_sel <- data.frame(iter = rep(1:niter, each = 5),
                               transport.method = method,
                               method = factor(rep(c("Binary Programming","W2","Hahn-Carvalho","Stepwise","Simulated Annealing"),
                                                   niter))
    )

    label.df_proj <- data.frame(iter = rep(1:niter, each = 6),
                                transport.method = method,
                                method = factor(rep(c("W2","Hahn-Carvalho","Stepwise","Simulated Annealing", "W1", "WInfty"),
                                                    niter))
    )


    dist.list <- list(parameters = list(inSamp = NULL,
                                       single = NULL),
                      predictions = list(inSamp = NULL,
                                  single = NULL))
    mses <- dist.list
    distance <- dist.list

    for(i in c("inSamp", "newX","single")) {
      for(j in c("parameters","predictions")){
        mses[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$mse[[i]][[j]]))
        distance[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$W2_dist[[i]][[j]]))
      }
    }

    #### Single Point ####
    w2single <- msesingle <- list()
    # if (neighb == "single") {
    w2sellist <- lapply(outputs, function(o) o$W2_dist$single$selection  )
    msesellist<- lapply(outputs, function(o) o$mse$single$selection )

    if(!all(sapply(w2sellist, is.null))) {
      w2single$selection <- WpProj:::combine.distcompare(w2sellist)
    }
      w2single$projection <- WpProj:::combine.distcompare(lapply(outputs, function(o) o$W2_dist$single$projection ))
    if (!all(sapply(msesellist, is.null))) {
      msesingle$selection <- WpProj:::combine.distcompare(msesellist)
    }
      msesingle$projection <- WpProj:::combine.distcompare(lapply(outputs, function(o) o$mse$single$projection ))
    # }

    #
    w2_opt <- mse_opt <- list()
    # cols <- c("dist","nactive", "groups", "method")
    levs <- c("Binary Programming", "Relaxed B.P.","L1","L2","LInf")

    #Wass
    w2_opt$selection <- lapply(w2single$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
    w2_opt$selection$p <- w2single$selection$p
    w2_opt$projection <- lapply(w2single$projection[1:2], function(res) res[res$groups %in% c("L1", "Lasso", "LInf"),])
    w2_opt$projection$p <- w2single$projection$p

    for(pp in names(w2_opt)) {
      for(nn in names(w2_opt$selection)[1:2]) {
        if(is.null(w2_opt[[pp]][[nn]])) next
        w2_opt[[pp]][[nn]] <-  w2_opt[[pp]][[nn]] %>% filter(nactive > 0)
        w2_opt[[pp]][[nn]]$groups <- droplevels(w2_opt[[pp]][[nn]]$groups)

        w2_opt[[pp]][[nn]]$groups <- droplevels(w2_opt[[pp]][[nn]]$groups)
        if (pp == "projection") {
          levels(w2_opt[[pp]][[nn]]$groups) <- c("L1","L2","LInf")
          w2_opt[[pp]][[nn]]$groups <- factor(as.character(w2_opt[[pp]][[nn]]$groups), levels = levs)
        } else {
          levels(w2_opt[[pp]][[nn]]$groups) <- levs
        }

      }
    }
    class(w2_opt$selection) <-
      class(w2_opt$projection) <-
      unique(c(class(w2single$projection), class(w2single$selection)))
    w2names <- unique(c(
      names(w2_opt$selection)[1:2],
      names(w2_opt$projection)[1:2]
      ))
    w2_opt$combined <- lapply(w2names, function(nn) do.call("rbind", list(w2_opt$selection[[nn]], w2_opt$projection[[nn]])))
    w2_opt$combined[[3]] <- unique(c(w2_opt$selection$p, w2_opt$projection$p))
    names(w2_opt$combined) <- c("parameters", "predictions","p")
    class(w2_opt$combined) <- unique(c(class(w2single$selection),
                                       class(w2single$projection)))

    # mse
    msenames <- unique(c(
      names(msesingle$selection)[1:2],
      names(msesingle$projection)[1:2]
    ))
    mse_opt$selection <- lapply(msesingle$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
      mse_opt$selection$p <- msesingle$selection$p

    mse_opt$projection <- lapply(msesingle$projection[1:2], function(res) res[res$groups %in% c("L1","Lasso","LInf"),])
    mse_opt$projection$p <- msesingle$projection$p

    for(pp in names(mse_opt)) {
      for(nn in msenames) {
        if(is.null(mse_opt[[pp]][[nn]])) next
        mse_opt[[pp]][[nn]] <-  mse_opt[[pp]][[nn]] %>% filter(nactive > 0)
        mse_opt[[pp]][[nn]]$groups <- droplevels(mse_opt[[pp]][[nn]]$groups)

        mse_opt[[pp]][[nn]]$groups <- droplevels(mse_opt[[pp]][[nn]]$groups)
        if (pp == "projection") {
          levels(mse_opt[[pp]][[nn]]$groups) <- c("L1","L2","LInf")
          mse_opt[[pp]][[nn]]$groups <- factor(as.character(mse_opt[[pp]][[nn]]$groups), levels = levs)
        } else {
          levels(mse_opt[[pp]][[nn]]$groups) <- levs
        }
      }
    }


      class(mse_opt$projection) <-
      unique(c(
        class(msesingle$selection),
        class(msesingle$projection)
      ))
    if(!is.null(mse_opt$selection)) class(mse_opt$selection) <- class(mse_opt$projection)

    mse_opt$combined <- lapply(msenames[1:2], function(nn) do.call("rbind", list(mse_opt$selection[[nn]], mse_opt$projection[[nn]])))
    mse_opt$combined[[3]] <- unique(c(mse_opt$selection$p, mse_opt$projection$p))

    names(mse_opt$combined) <- c("parameters", "predictions","p")
    class(mse_opt$combined) <- unique(c(
      class(msesingle$selection),
      class(msesingle$projection)
    ))

    # debugonce(plot.combine.distcompare)
    # plot(w2single$selection, alpha = 0.5)
    # plot(msesingle$selection, alpha = 0.5)
    # plot(w2single$projection, alpha = 0.5)
    # plot(msesingle$projection, alpha = 0.5)
    w2_list[[corr]] <- w2_opt$combined
    mse_list[[corr]] <- mse_opt$combined

    for (i in 1:2) {
      if(!is.null(w2_list[[corr]][[i]])) w2_list[[corr]][[i]]$corr <- corr
      if(!is.null(mse_list[[corr]][[i]])) mse_list[[corr]][[i]]$corr <- corr
    }

    #### W2 r2 stuff ####
    w2r2single <- w1r2single <- list()
    fix.fun <- function(oo) {
      class(oo) <- c("WPR2", class(oo))
      return(oo)
    }

    # w2r2single$selection  <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W2_r2$null$single$selection)))
    # w2r2single$projection <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W2_r2$null$single$projection)))
    # w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W1_r2$null$single$selection)))
    # w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) fix.fun(o$W1_r2$null$single$projection)))


    # w2r2single$selection  <- combine.WPR2(lapply(outputs, function(o) (o$W2_r2$expectation$single$selection)))
    # w2r2single$projection <- combine.WPR2(lapply(outputs, function(o) (o$W2_r2$expectation$single$projection)))
    w1r2list <- lapply(outputs, function(o) (o$W1_r2$expectation$single$selection))
    w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$expectation$single$projection)))
    if(!all(sapply(w1r2list, is.null))) {
      w1r2single$selection <- WpProj:::combine.WPR2(w1r2list)
      w1r2_list[[corr]] <- combine.WPR2(w1r2single$selection %>% filter(groups %in% c("Binary Programming",  "Lasso")) %>%
                                          mutate(groups = fct_recode(groups, "Relaxed B.P." = "Lasso")) %>%
                                          mutate(groups = fct_drop(groups)),
                                        w1r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
                                          mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
                                          mutate(groups = fct_drop(groups)))

    } else {
      w1r2_list[[corr]] <-  w1r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
        mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
        mutate(groups = fct_drop(groups))
    }

    w2r2list <-lapply(outputs, function(o) (o$W2_r2$null$single$selection ))
    w2r2single$projection <- WpProj:::combine.WPR2(lapply(outputs, function(o) {
      temp <- o$W2_r2$null$single$projection %>%
        filter(nactive > 0)
      dists <- o$W2_dist$single$projection$predictions %>%
        filter(nactive > 0)
      old_base <- dists$dist
      # new_base <- o$W2_dist$single$projection$predictions %>%
      #   # group_by(groups) %>%
      #   mutate(base = max(dist)) %>%
      #   filter(nactive > 0)
      new_base <- dists %>% group_by(groups) %>% mutate(base = max(dist))
      temp$r2 <- 1- old_base^2 / new_base$base^2
      temp$base <- "dist.from.null"
      temp <- temp[,c("r2","nactive","groups","method","p","base")]
      return(temp)
    }))
    if(!all(sapply(w2r2list, is.null))) {
      w2r2single$selection <- WpProj:::combine.WPR2(lapply(outputs, function(o) (o$W2_r2$null$single$selection %>% filter(nactive > 0))))
      w2r2single$selection  <- combine.WPR2(lapply(outputs, function(o) {
        temp <- o$W2_r2$null$single$selection %>%
          filter(nactive > 0)
        dists <- o$W2_dist$single$selection$predictions %>%
          filter(nactive > 0)
        old_base <- dists$dist
        # new_base <-
          # max(o$W2_dist$single$projection$predictions %>%
          # summarize(md = max(dist)),
          # o$W2_dist$single$selection$predictions %>%
          #   summarize(md = max(dist)))
        # new_base <-
        # o$W2_dist$single$selection$predictions %>%
          # filter(nactive > 0)
          # group_by(groups) %>%
          # mutate(base = max(dist)) %>%

        new_base <- dists %>%
          group_by(groups) %>% mutate(base = max(dist))
        temp$r2 <- 1 - old_base^2 / new_base$base^2
        temp$base <- "dist.from.null"
        temp <- temp[,c("r2","nactive","groups","method","p","base")]
        return(temp)
      }))
      w2r2_list[[corr]] <- combine.WPR2(w2r2single$selection %>% filter(groups %in% c("Binary Programming",  "Lasso")) %>%
                                          mutate(groups = fct_recode(groups, "Relaxed B.P." = "Lasso")) %>%
                                          mutate(groups = fct_drop(groups)),
                                        w2r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
                                          mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
                                          mutate(groups = fct_drop(groups)))
    } else {
      w2r2_list[[corr]] <- w2r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
        mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
        mutate(groups = fct_drop(groups))
    }
    # w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$null$single$selection %>% filter(nactive >0))))
    # w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$null$single$projection %>% filter(nactive >0))))

    # plot(w2r2single$selection)
    # plot(w2r2single$projection)

    # w2r2single$projection <- combine.WPR2(lapply(outputs, function(o) {
    #   temp <- o$W2_r2$null$single$projection %>% filter(nactive > 0)
    #   dists <- o$W2_dist$single$projection$predictions %>%
    #     filter(nactive > 0)
    #   old_base <-  dists$dist
    #   new_base <- dists %>%
    #     # group_by(groups) %>% mutate(base = dist[which.min(nactive)])
    #     group_by(groups) %>% mutate(base = max(dist))
    #   temp$r2 <- 1- old_base^2 / new_base$base^2
    #   temp$base <- "dist.from.null"
    #   return(temp)
    # }))
    # w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) {
    #   temp <- o$W1_r2$null$single$selection
    #   old_base <- o$W2_dist$single$selection$predictions$dist
    #   new_base <- o$W2_dist$single$selection$predictions %>%
    #     group_by(groups) %>% mutate(base = dist[which.min(nactive)])
    #   temp$r2 <- 1- old_base / new_base$base
    #   return(temp)
    # }))
    # w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) {
    #   temp <- o$W1_r2$null$single$projection
    #   old_base <- o$W2_dist$single$projection$predictions$dist
    #   new_base <- o$W2_dist$single$selection$predictions %>%
    #     group_by(groups) %>% summarise(base = dist[which.min(nactive)])
    #   temp$r2 <- 1- old_base^2 / new_base$base^2
    #   return(temp)
    # }))

    # w


    w1r2_list[[corr]]$corr <- corr
    w2r2_list[[corr]]$corr <- corr

    #### derivatives ####
    if(family == 'binomial') {
      derivative.fun <- function(x) {
        if(!is.matrix(x)) {
          x <- as.matrix(x)
          if(dim(x)[2] == 1) x <- t(x)
        }
        if(all(x[,1] == 1)) x <- x[,2:ncol(x), drop = FALSE]
        p <- ncol(x)
        derivs <- matrix(0,nrow = nrow(x), ncol = p)
        names(derivs) <- names(x)

        derivs[,1] <- cos(pi/8 *  x[,1]* x[,2]) * pi/8 * x[,2]
        derivs[,2] <- cos(pi/8 *  x[,1]* x[,2]) * pi/8 * x[,1]
        derivs[,6] <- ifelse(x[,6]^2 > 3/8* pi, 0, -sin(2*x[,6]^2 + pi/2)* 4 * x[,6])
        derivs[,7] <- exp(-x[,7]^2/5) * (6/5 * x[,7]^2 - 1/4 * x[,7] - 2) -
          (2/5 * (x[,7]^3 - 1/8 * x[,7]^2 - 2 * x[,7]) * exp(- x[,7]^2/5) * 2/5 * x[,7])
        derivs[,11]<- -1.0/cosh(pi/8 * x[,11] * x[,15]^2) * pi/8 * x[,15]^2
        derivs[,15]<- -1.0/cosh(pi/8 * x[,11] * x[,15]^2) * pi/4 * x[,15] * x[,11]
        # derivs     <- cbind(0,derivs)
        return( derivs )
      }
      test.points <- lapply(outputs, function(o) o$data$test$X[,-1,drop=FALSE])
      proj.output <- lapply(outputs, function(o) o$models$projection)


      true.derivatives <- lapply(test.points, derivative.fun)
      nn.derivatives <- lapply(outputs, function(o) o$models$derivatives)
      nn.orig.derivatives <- lapply(outputs, function(o) o$models$estimation$original$derivatives[1,])
      d.dist.list   <- mapply(function(m, d, nn,nn.o) {
        rd <- sapply(m, function(ll) sapply(ll, function(tt) mean((tt[-1,,drop= FALSE] - c(d))^2)))
        # rnn <- sapply(m, function(ll) sapply(ll, function(tt)  mean((tt[-1,,drop= FALSE] - nn)^2)))
        dnn <- mean((c(d)- nn)^2)
        dnn.o <- mean((c(d) - nn.o)^2)
        rel <- unlist(rd)/dnn
        rel.o <- unlist(rd)/dnn.o

        w2.dist <- sapply(m, function(ll) sapply(ll, function(tt) WpProj::wasserstein(X = tt[-1,,drop= FALSE], Y = nn, p = 2,
                                                                                     observation.orientation = "colwise",
                                                                                     method = "exact")))

        boot <- data.frame(
          dist = c(rel),
          nactive = unlist(sapply(m, function(ll) sapply(ll, function(tt) c(unique(colSums(tt[-1,,drop= FALSE] !=0)))))),
          groups = unlist(sapply(names(m), function(nm) rep(nm, length(m[[nm]])))),
          method = "mse",
          ranks = NA,
          iter = NA)
        orig <- data.frame(
          dist = c(rel.o),
          nactive = unlist(sapply(m, function(ll) sapply(ll, function(tt) c(unique(colSums(tt[-1,,drop= FALSE] !=0)))))),
          groups = unlist(sapply(names(m), function(nm) rep(nm, length(m[[nm]])))),
          method = "mse",
          ranks = NA,
          iter = NA)
        w2 <- data.frame(
          dist = c(unlist(w2.dist)),
          nactive = unlist(sapply(m, function(ll) sapply(ll, function(tt) c(unique(colSums(tt[-1,,drop= FALSE] !=0)))))),
          groups = unlist(sapply(names(m), function(nm) rep(nm, length(m[[nm]])))),
          method = "w2",
          ranks = NA,
          iter = NA)
        return(list(boot = boot, orig = orig,
                    w2 = w2))
      },
      m = proj.output, d = true.derivatives, nn = nn.derivatives, nn.o = nn.orig.derivatives, SIMPLIFY = FALSE)
      for(ndd in 1:length(d.dist.list)) {
        d.dist.list[[ndd]]$boot$iter <- d.dist.list[[ndd]]$orig$iter <-
          d.dist.list[[ndd]]$w2$iter  <- ndd
        d.dist.list[[ndd]]$boot$corr <- d.dist.list[[ndd]]$orig$corr <-
          d.dist.list[[ndd]]$w2$corr <- corr
      }
      d.dist <- list(parameters = do.call("rbind",lapply(d.dist.list, function(d) d$boot)) %>%  filter(nactive >0),
                     predictions = do.call("rbind",lapply(d.dist.list, function(d) d$w2))  %>%  filter(nactive >0), p = w2_list[[1]]$p-1)
      class(d.dist) <- class(w2_list[[1]])
      deriv_list[[corr]] <- d.dist
      # plot(d.dist, alpha = 0.5, base_size = basesize, ylab = "MSE between coef. and deriv.", facet.group = "corr", CI = "none")
      # theta  <- lapply(proj, function(p) p)

    }
  }

  w2_stack <- mse_stack <- deriv_stack <- w2_deriv_stack <- vector("list", 2)
  names(w2_stack) <- names(mse_stack) <-
    names(deriv_stack) <- c("parameters", "predictions")
  for (i in 1:2) {
    checklist <- do.call( "rbind", lapply(w2_list, function(w) w[[i]]) )
    if(!is.null(checklist)) w2_stack[[i]]  <- checklist
    checklist <-  do.call( "rbind", lapply(mse_list, function(w) w[[i]]))
    if(!is.null(checklist)) mse_stack[[i]] <- checklist
    if(family == "binomial") {
      checklist <-  do.call( "rbind", lapply(deriv_list, function(w) w[[i]]))
      if(!is.null(checklist)) deriv_stack[[i]] <- checklist
    }
  }


  if (family == "binomial") {
    mse_stack$predictions <- mse_stack$predictions %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    w2_stack$predictions <-  w2_stack$predictions  %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    deriv_stack$predictions <- deriv_stack$predictions %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.",
                                                                    "Stepwise","Simulated Annealing","Hahn-Carvalho")))
    deriv_stack$parameters <- deriv_stack$parameters %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.",
                                                                              "Stepwise","Simulated Annealing","Hahn-Carvalho")))
    class(deriv_stack) <-  class(w2_list[[1]])

    deriv_stack$predictions <- deriv_stack$predictions %>% group_by(corr, method, groups, nactive) %>%
      summarise(lwr = quantile(dist, .25),
                upr = quantile(dist, .75),
                dist = mean(dist)) %>%
      relocate(dist, nactive, groups, method, corr, lwr, upr)
    deriv_stack$parameters <- deriv_stack$parameters %>% group_by(corr, method, groups, nactive) %>%
      summarise(lwr = quantile(dist, .25),
                upr = quantile(dist, .75),
                dist = mean(dist)) %>%
      relocate(dist, nactive, groups, method, corr, lwr, upr)
  } else {
    mse_stack$parameters <- mse_stack$parameters %>%
      group_by(corr, method, groups, nactive) %>%
      summarise(lwr = quantile(dist, .25),
                upr = quantile(dist, .75),
                dist = mean(dist)) %>%
      relocate(dist, nactive, groups, method, method, corr, lwr, upr) #%>%
    # group_by(corr, method, iter) %>%
    # summarise(dist = dist/min(dist), nactive = nactive,
    #           ranks = NA, groups = groups) %>%
    # relocate(dist, nactive, groups, method, ranks, iter, corr)
    w2_stack$parameters <- w2_stack$parameters %>%
      group_by(corr, method, groups, nactive) %>%
      summarise(lwr = quantile(dist, .25),
                upr = quantile(dist, .75),
                dist = mean(dist)) %>%
      relocate(dist, nactive, groups, method, method, corr,
               lwr, upr)
  }
  mse_stack$predictions <- mse_stack$predictions %>%
    group_by(corr, method, groups, nactive) %>%
    summarise(lwr = quantile(dist, .25),
              upr = quantile(dist, .75),
              dist = mean(dist) ) %>%
    relocate(dist, nactive, groups, method, corr, lwr, upr)
  w2_stack$predictions <- w2_stack$predictions %>%
    group_by(corr, method, groups, nactive) %>%
    summarise(lwr = quantile(dist, .25),
              upr = quantile(dist, .75),
              dist = mean(dist) ) %>%
    relocate(dist, nactive, groups, method, corr, lwr, upr)

  class(mse_stack) <- class(w2_stack) <- class(w2_list[[1]])

  if ("NULL" %in% class(mse_stack)) {
    cmse <- class(mse_stack)
    class(mse_stack) <- cmse[cmse != "NULL"]
  }
  if ("NULL" %in% class(w2_stack)) {
    cw2 <- class(w2_stack)
    class(w2_stack) <- cw2[cw2 != "NULL"]
  }
  msep <- plot(x = mse_stack,
               alpha = 0.5,
               base_size = basesize, ylab = "Relative MSE",
               facet.group = "corr", CI = "ribbon")
  w2p <- plot(w2_stack,  alpha = 0.5, base_size = basesize, ylab = "2-Wasserstein Distance", facet.group = "corr", CI = "none")

  if(family == "binomial") {
    msep$predictions <- msep$predictions + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    w2p$predictions  <- w2p$predictions  + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])

    if ("NULL" %in% class(deriv_stack)) {
      cderiv <- class(deriv_stack)
      class(deriv_stack) <- cderiv[cderiv != "NULL"]
    }
    deriv_stack_p <- plot(deriv_stack, alpha = 0.5, base_size = basesize, ylab = "Relative MSE", facet.group = "corr", CI = "none")

    deriv_stack_p$predictions <-  deriv_stack_p$predictions + ggplot2::ylab("2-Wasserstein Distance") + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])

    deriv_stack_p$parameters <-  deriv_stack_p$parameters + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
  }

  for (i in 1:2) {
    if(!is.null(msep[[i]])) {
      msep[[i]] <- msep[[i]] + theme(strip.background = element_blank(),
                                     strip.text.x = element_blank()) +
        scale_x_continuous(breaks = c(1, seq(5,25,5)))
      msep[[i]] <- remove_geoms(msep[[i]], "GeomRibbon", FALSE)
    }
    if(!is.null(w2p[[i]])) {
      w2p[[i]]  <- w2p[[i]] + theme(strip.background = element_blank(),
                                    strip.text.x = element_blank()) +
        scale_x_continuous(breaks = c(1, seq(5,25,5)))
      w2p[[i]]  <- remove_geoms(w2p[[i]], "GeomRibbon", FALSE)
    }
  }
  # pdf(file = "inst/figure/simulation/mse_normal_neighb.pdf", width = 7, height = 5)
  # msep$predictions
  # dev.off()
  #
  # pdf(file = "inst/figure/simulation/w2_normal_neighb.pdf", width = 7, height = 5)
  # w2p$predictions
  # dev.off()
  if(family == "binomial") {
    for (i in 1:2) {
      if(!is.null(deriv_stack_p[[i]])) {
        deriv_stack_p[[i]] <- deriv_stack_p[[i]] +
          scale_x_continuous(breaks = c(1, seq(5,25,5))) +
          theme(strip.background = element_blank(),
                strip.text.x = element_blank())
        deriv_stack_p[[i]] <- remove_geoms(deriv_stack_p[[i]], "GeomRibbon", FALSE)
      }
    }
    directory <- file.path("inst", "figure","rawplots")
    if(!dir.exists(directory)) dir.create(directory, recursive = TRUE)
    w2fn <- file.path(directory, paste0(paste0(c(neighb,family,gsub("[.]","_","mcp.net") , pf ,method, "w2plot"), collapse="_"), ".rds"))
    saveRDS(list(mse = msep, w2 = w2p, deriv = deriv_stack_p), file = w2fn)
  } else {
    directory <- file.path("inst", "figure","rawplots")
    if(!dir.exists(directory)) dir.create(directory, recursive = TRUE)
    w2fn <- file.path(directory, paste0(paste0(c(neighb,family,gsub("[.]","_","mcp.net") , pf ,method, "w2plot"), collapse="_"), ".rds"))
    saveRDS(list(mse = msep, w2 = w2p, deriv = NULL), file = w2fn)
  }
  # msefile <- file.path(directory, paste0(paste0(c(neighb,family,penalty , pf ,method, "mseplot"), collapse="_"), ".rds"))
  # pdf(file, width = width, height = height)
  # msep$predictions
  # dev.off()
  #
  # w2file <- file.path(directory, paste0(paste0(c(neighb,family,penalty , pf ,method, "w2plot"), collapse="_"), ".rds"))
  # pdf(file, width = width, height = height)
  # w2p$predictions
  # dev.off()

  # pdf(file = "inst/figure/simulation/w2r2_normal_neighb.pdf", width = 7, height = 5)
  # plot(w2r2, alpha = 0.5, base_size = 20)
  # dev.off()
  #
  # pdf(file = "inst/figure/simulation/w1r2_normal_neighb.pdf", width = 7, height = 5)
  # plot(w1r2, alpha = 0.5, base_size = 20)
  # dev.off()
  #
  #
  # pdf("for_lorenzo.pdf",  width = 12, height = 7)
  # msep$predictions
  # w2p$predictions
  # plot(w2r2, alpha = 0.5, base_size = 20)
  # plot(w1r2, alpha = 0.5, base_size = 20)
  # dev.off()
  w1r2 <- do.call("rbind", w1r2_list)
  w2r2 <- do.call("rbind", w2r2_list)
  if(family == "binomial") {
    if (is.data.frame(w1r2 ) ) {
      w1r2 <- w1r2 %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    }
    if ( is.data.frame(w2r2) ) {
      w2r2 <- w2r2 %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))

    }
    w2r2deriv <- deriv_stack$predictions %>% group_by(corr,groups, method) %>%
      reframe(r2 = 1-dist^2/max(dist)^2, nactive = nactive, p = 2, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      reframe(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    w1r2deriv <- deriv_stack$predictions %>% group_by(corr,groups, method) %>%
      reframe(r2 = 1-dist/max(dist), nactive = nactive, p = 1, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      reframe(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    class(w2r2deriv) <- class(w1r2deriv) <- c("WPR2", class(w1r2deriv))
    class(w2r2) <- class(w1r2) <- class(w1r2_list[[1]])
  }
  # w2r2.test <- w2_stack$predictions %>% group_by(groups, corr, nactive) %>%
  #   summarise(dist = mean(dist)) %>%
  #   group_by(groups, corr) %>%
  #   mutate(r2 = 1 - dist^2/max(dist^2))
  # class(w2r2.test) <-class(w1r2)
  w1r2p <-NULL
  if (inherits(w1r2, "WPR2")) {
    w1r2p <- plot(w1r2, alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE)
  }
  w2r2p <- NULL
  if (inherits(w2r2, "WPR2")) {
    w2r2p <- plot(w2r2 #%>% filter(groups != "Binary Programming")
                  , alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE,
                  ylim = c(0,1)) #+ ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[2:5])
  }

  # w1r2post <- w1r2
  # temp <- o$W2_r2$null$single$projection %>% filter(nactive > 0)
  # dists <- o$W2_dist$single$projection$predictions %>%
  #   filter(nactive > 0)
  # old_base <-  dists$dist
  # new_base <- dists %>%
  #   # group_by(groups) %>% mutate(base = dist[which.min(nactive)])
  #   group_by(groups) %>% mutate(base = max(dist))
  # temp$r2 <- 1- old_base^2 / new_base$base^2
  # temp$base <- "dist.from.null"
  w2r2derivp <- NULL
  if(family == "binomial") {
    if(!is.null(w1r2p)) {
      w1r2p <- w1r2p + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
        ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    }

    if (!is.null(w2r2p)) {
      w2r2p <-  w2r2p + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
        ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    }

    w1r2derivp <- plot(w1r2deriv, alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE)+
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
    w2r2derivp <- plot(w2r2deriv #%>% filter(groups != "Binary Programming")
                            , alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE,
                            ylim = c(0,1)) + #ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[2:5]) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
    w2r2derivp <- list(parameters = NULL, predictions = w2r2derivp)
  }
  if (!is.null(w2r2p)) {
  w2r2p <- w2r2p + theme(strip.background = element_blank(),
                         strip.text.x = element_blank()) +
    scale_x_continuous(breaks = c(1, seq(5,25,5)))
  }
  if (!is.null(w1r2p)) {
  w1r2p <- w1r2p + theme(strip.background = element_blank(),
                         strip.text.x = element_blank()) +
    scale_x_continuous(breaks = c(1, seq(5,25,5)))
  }

  w1r2postp <- w2r2postp <- NULL
  if(!is.null(w2_stack$parameters)) {
    w1r2post <- w2_stack$parameters %>% group_by(corr,groups, method) %>%
      reframe(r2 = 1-dist/max(dist), nactive = nactive, p = 1, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      reframe(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    w2r2post <- w2_stack$parameters %>% group_by(corr,groups, method) %>%
      reframe(r2 = 1-dist^2/max(dist)^2, nactive = nactive, p = 2, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      reframe(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    class(w1r2post) <- class(w2r2post) <- class(w2r2)
    w1r2postp <- plot(w1r2post, alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
    w2r2postp <- plot(w2r2post #%>% filter(groups != "Binary Programming")
                           , alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE,
                           ylim = c(0,1)) + #ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[2:5]) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
  }


  w2r2p <- list(parameters = w2r2postp, predictions = w2r2p)
  w1r2p <- list(parameters = w1r2postp, predictions = w1r2p)

  wpr2fn <- file.path("inst", "figure","rawplots", paste0(paste0(c(neighb, family, gsub("[.]","_",penalty) , pf ,method, "wpr2plot"), collapse="_"), ".rds"))
  saveRDS(list(w1r2 = w1r2p, w2r2 = w2r2p, w2r2deriv = w2r2derivp), file = wpr2fn)
}

#binomial
# debugonce(neighb.plot.fun)
neighb.plot.fun(family, n, p, method, pf, penalty, corrs, date, neighb[1], basesize = 11)
# neighb.plot.fun(family, n, p, method, pf, penalty, corrs, date, neighb[2], basesize = 11)

#gaussian
family <- "gaussian"
n <- 1024
p <- 21
# neighb <- c("Neighborhood", "Single_point")
#debugonce(neighb.plot.fun)
neighb.plot.fun(family, n, p, method, pf, penalty, corrs, date, neighb[1], basesize = 11)
# neighb.plot.fun(family, n, p, method, pf, penalty, corrs, date, neighb[2], basesize = 11)


#### generate plots ####
cube.root <- function(x) x^(1/3)
cube.root_trans <- function(){
  scales::trans_new("cube.root", function(x) x^(1/3), function(x)x^3,
                    domain = c(0, Inf))
}

combine_plot_temp <- function(family, which.plt = c("mse","w2", "deriv","w1r2","w2r2","w2r2deriv"),
                              height = 5, width = 7, breaks = NULL, trans = NULL, rm.layer = NULL,
                              limits = list(NULL),
                              # rel.heights = NULL,
                              # predictions = c("mse","w2", "deriv","w1r2","w2r2"),
                              parameters = NULL,
                              sep.parameters = FALSE,
                              label.map = NULL,
                              labels = NULL, pal = NULL) {
  remove_geom <- function(ggplot2_object, geom_type) {
    # Delete layers that match the requested type.
    layers <- lapply(ggplot2_object$layers, function(x) {
      if (class(x$geom)[1] == geom_type) {
        NULL
      } else {
        x
      }
    })
    # Delete the unwanted layers.
    layers <- layers[!sapply(layers, is.null)]
    ggplot2_object$layers <- layers
    ggplot2_object
  }
  lblr <- function(x) {
    ifelse(x == 0, "0",
           ifelse(x < 1, sprintf("%.2f", x),
                  sprintf("%.0f", x))
    )
  }
  penalty <- "mcp.net"
  pf <- "none"
  method <- "exact"
  # corr <- c("Corr_0", "Corr_0.5", "Corr_0.9")
  neighb <- c("single")
  which.plt <- match.arg(which.plt, several.ok = TRUE)
  if(is.null(breaks)) breaks <- list(NULL)
  if(is.null(trans)) trans <- list(NULL)

  # plts <- lapply(neighb, function(nn) vector("list",2))
  # names(plts) <- neighb
  plts <- list()
  for(cur_neighb in neighb) {
    # tempplts <-vector("list", length(corr))
    # names(tempplts) <- corr
    # for(cur_corr in corr){
    plt1 <- plt2 <- NULL
    if(any(which.plt %in% c("mse","w2", "deriv"))) {
      fn1 <- file.path("inst", "figure","rawplots", paste0(paste0(c(cur_neighb,family,gsub("[.]","_","mcp.net") , pf ,method, "w2plot"), collapse="_"), ".rds"))
      plt1 <- readRDS(fn1)
    }
    if (any(which.plt %in% c("w1r2","w2r2", "w2r2deriv"))) {
      fn2 <- file.path("inst", "figure","rawplots", paste0(paste0(c(cur_neighb, family,gsub("[.]","_","mcp.net") , pf ,method, "wpr2plot"), collapse="_"), ".rds"))
      plt2 <- readRDS(fn2)
    }
    plts <- c(plt1,
              plt2 )

    plts <- plts[which.plt]
    if(!is.null(parameters)) {
      stopifnot(parameters %in% which.plt)
      runparameters <- TRUE
    } else {
      runparameters <- FALSE
    }
    # }
    # plts[[cur_neighb]] <- do.call("rbind", )
    outfile <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, gsub("[.]","_","mcp.net"), pf, method, which.plt), collapse="_"), ".pdf"))
    # if(!is.null(parameters)) outfile <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, gsub("[.]","_","mcp.net"), pf, method, which.plt, "parameters"), collapse="_"), ".pdf"))
    if (sep.parameters) {
      outfile.predictions <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, gsub("[.]","_","mcp.net"), pf, method, which.plt), collapse="_"), "_predictions.pdf"))
      outfile.parameters <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, gsub("[.]","_","mcp.net"), pf, method, which.plt), collapse="_"), "_parameters.pdf"))

    }

    print.p.list <- vector("list", length(plts))
    names(print.p.list) <- names(plts)
    for (i in names(plts)) {
      pp <- list()
      if(i %in% parameters) {
        idx <- c("predictions","parameters")
      } else {
        idx <- "predictions"
      }
      if (i == "deriv") idx <- idx[2:1]
      for (j in idx) {
        express <- if( j == "parameters" &
                       runparameters & !is.null(plts[[i]][[j]]) &
                       i != "deriv") {
          ylab(expression(W[2](beta, theta)))
        } else if (i == "deriv" & j != "parameters" & !is.null(plts[[i]][[j]])) {
          ylab(expression(W[2](nabla[x]~mu, nabla[x]~nu)))
        } else {
          ylab(expression(W[2](mu, nu)))
        }
        pp[[j]] <- if( (i == "deriv" & j == "predictions") | i == "w2") {
          plts[[i]][[j]] + express
          # } else if (i != "w2r2" & i != "w1r2") {
          #   plts[[i]][[j]]
          # } else if (i %in% parameters) {

        } else {
          plts[[i]][[j]]
        }

        if(!is.null(trans[[i]][[j]]) & !isFALSE(trans[[i]][[j]])) {
          if(isTRUE(trans[[i]][[j]])) trans[[i]][[j]] <- "sqrt"
          if(is.null(breaks[[i]][[j]])) breaks[[i]][[j]] <- waiver()
          pp[[j]] <- pp[[j]] + expand_limits(y = 0) +
            scale_y_continuous(trans = trans[[i]][[j]],
                               expand = c(0.01,0.05),
                               breaks = breaks[[i]][[j]],
                               limits = limits[[i]][[j]],
                               labels = lblr)
        }

        if(!is.null(rm.layer)) pp[[j]] <- remove_geoms(pp[[j]], rm.layer, FALSE)

        if(!is.null(labels) & !is.null(pal)) {
          grps <- sort(unique(as.character(pp[[j]]$data$groups)))
          if(!is.null(label.map)) grps <- levels(forcats::fct_relevel(grps, label.map))
          pp[[j]] <- pp[[j]] +
            scale_color_manual(name = "Method:",
                               breaks = grps, #c("L1", "L2", "LInf"),
                               labels = labels,
                               values = pal) +
            scale_fill_manual(name = "Method:",
                              breaks = grps, #c("L1", "L2", "LInf"),
                              labels = labels,
                              values = pal) +
            scale_size_manual(name = "Method:",
                              breaks = grps, #c("L1", "L2", "LInf"),
                              labels = labels,
                              values = pal)
        }
      }


      print.p.list[[i]] <- pp

    }
    print.p <- list()
    print.p$predictions <- unlist(print.p.list, recursive = FALSE)
    if (sep.parameters) {
      temp <- unlist(print.p.list[names(plts) %in% parameters], recursive = FALSE)
      print.p$parameters <- temp[grepl("parameters", names(temp))]
      print.p$predictions <- print.p$predictions[!grepl("parameters", names(print.p$predictions))]
    }
    plot.legend <- get_legend(print.p$predictions[[1]] + ggplot2::theme(legend.position="bottom"))
    for(j in names(print.p)) {
      for(i in seq_along(print.p[[j]])) {
        if(i == 1 & j == "predictions") {
          rho0 <- data.frame(groups = "L1",
                             corr="Corr_0",
                             nactive = 0.5,
                             hi = 0.0,
                             low = 0.5,
                             dist = 0.05,
                             # .group = 1,
                             lab = "'Correlation:'~rho~'='~0.0"
          )
          rho05 <- data.frame(groups = "L1",
                              corr="Corr_0.5",
                              nactive = 0.5,
                              hi = 0.0,
                              low = 0.5,
                              dist = 0.05,
                              # .group = 2,
                              lab = "rho~'='~0.5"
          )
          rho09 <- data.frame(groups = "L1",
                              corr="Corr_0.9",
                              nactive = 0.5,
                              hi = 0.0,
                              low = 0.5,
                              dist = 0.05,
                              # .group = 3,
                              lab = "rho~'='~0.9"
          )
          print.p[[j]][[i]] <- print.p[[j]][[i]] + geom_text(data = rho0,
                                                             aes(label = lab), parse = TRUE,
                                                             color = "black",
                                                             hjust=0) +
            geom_text(data = rho05,
                      aes(label = lab), parse = TRUE,
                      color = "black",
                      hjust=0) +
            geom_text(data = rho09,
                      aes(label = lab), parse = TRUE,
                      color = "black",
                      hjust=0)
          # print.p[[i]] <-
          # print.p[[i]] + geom_text(aes(label="rho~'='~0.0", x=5, y=0.05), parse = TRUE)
        }
        print.p[[j]][[i]] <- print.p[[j]][[i]] + xlab("") + theme(panel.spacing.x = unit(4, "mm")) +
          ggplot2::theme(legend.position="none") +
          ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5))
        fg <- sapply(print.p[[j]][[i]]$facet$params$rows, rlang::as_label)
        if ( is.character(fg) ) {
          print.p[[j]][[i]] <- print.p[[j]][[i]] + ggplot2::facet_grid(cols =  vars(!!sym(fg)))
        }

        if(i != length(print.p[[j]])) {
          print.p[[j]][[i]] <- print.p[[j]][[i]] +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank())
        }
      }
    }
    # plot.length <- length(print.p)
    # if(is.null(rel.heights)) rel.heights <- rep(3, plot.length)
    if(sep.parameters) {
      pred <- arrangeGrob(do.call("rbind", lapply(print.p$predictions, ggplotGrob)),
                          bottom = grid::textGrob("Number of active coefficients",
                                                  vjust = -1.8, hjust=0.4))
      param <- arrangeGrob(do.call("rbind", lapply(print.p$parameters, ggplotGrob)),
                           bottom = grid::textGrob("Number of active coefficients",
                                                   vjust = -1.8, hjust=0.4))
      pdf(outfile.predictions, width = width, height = height)
      print(grid.arrange(pred,
                         plot.legend, nrow = 2, heights = c(10,.5)))
      dev.off()
      pdf(outfile.parameters, width = width, height = height)
      print(grid.arrange(param,
                         plot.legend, nrow = 2, heights = c(10,.5)))
      dev.off()

    } else {
      sim.plots <- arrangeGrob(do.call("rbind", lapply(print.p$predictions, ggplotGrob)),
                               bottom = grid::textGrob("Number of active coefficients",
                                                       vjust = -1.8, hjust=0.4))
      pdf(outfile, width = width, height = height)
      print(grid.arrange(sim.plots,
                         plot.legend, nrow = 2, heights = c(10,.5)))
      dev.off()
    }

    # sim.plots <- arrangeGrob(grobs = print.p,
    #                               nrow = plot.length, ncol = 1,
    #                          heights = rel.heights,
    #                          bottom = grid::textGrob("Number of active coefficients",
    #                                                  vjust = -2, hjust=0.4))

  }

}

combine_plot_temp("binomial", which.plt = c("deriv","w2r2deriv"),   height = 6, width = 7.5,
                  trans = list(deriv = list(predictions = "sqrt", parameters = "sqrt")),
                  breaks = list(deriv = list(predictions = c(0, 0.25, 1, 4),
                                             parameters = c(0,0.25, 1,4))),
                  # rel.heights = c(3,5),
                  limits = list(deriv = list(predictions = c(0,4),
                                             parameters = c(0,4))),
                  labels = expression(W[1], W[2], W[infinity]),
                  label.map = c("L1", "Lasso", "LInf"),
                  pal = ggsci::pal_jama()(5)[3:5], parameters = c("deriv"))

combine_plot_temp("gaussian", c("mse","w2","w2r2"),   height = 6, width = 7.5,
                  breaks = list(mse = list(predictions = c(0, 1, 8, 27, 64),
                                           parameters = c(0, 1, 16, 64,144)),
                                w2 = list(predictions = c(0, 1,4,9,16),
                                          parameters = c(0, 0.25, 1, 4))),
                  limits = list(mse = list(predictions = c(0,70),
                                           parameters = c(0,144)),
                                w2 = list(predictions = c(0,16),
                                          parameters = c(0,4))),
                  rm.layer = "GeomRibbon",
                  trans = list(mse = list(predictions = "cube.root",
                                          parameters = "sqrt"),
                               w2 = list(predictions= "sqrt",
                                         parameters = "sqrt")),
                  labels = expression("B.P.", "Relaxed B.P.", W[1], W[2],
                                      W[infinity]), pal = ggsci::pal_jama()(5),
                  label.map = c("Binary Programming","Relaxed B.P.","L1", "L2", "LInf"),
                  parameters = c("mse","w2","w2r2"),
                  sep.parameters = TRUE)

