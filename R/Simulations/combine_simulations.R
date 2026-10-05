require(SLIMpaper)
library(dplyr)
library(forcats)
library(ggplot2)

family <- "binomial"
penalty <- "mcp.net"
pf <- "none"
method <- "exact"
corrs <- c("Corr_0", "Corr_0.5", "Corr_0.9") #"Corr_0.9" #
date <- "2020-10-23 12:00:00"
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
    folder <- file.path("Output", neighborhood, family, penalty, pf, method, corr,n,p)
    files <- list.files(folder, full.names = TRUE)
    files <- files[date.fun(files, date.start = date)]
    niter <- length(files)
    cat("Niter: ", niter, " Correlation:", corr, ", Family: ", family,"\n")
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


    dist.list <- list(posterior = list(inSamp = NULL,
                                       single = NULL),
                      mean = list(inSamp = NULL,
                                  single = NULL))
    mses <- dist.list
    distance <- dist.list

    for(i in c("inSamp", "newX","single")) {
      for(j in c("posterior","mean")){
        mses[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$mse[[i]][[j]]))
        distance[[j]][[i]] <- do.call("rbind", lapply(outputs, function(o) o$W2_dist[[i]][[j]]))
      }
    }

    #### Single Point ####
    # browser()
    w2single <- msesingle <- list()
    w2single$selection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$single$selection  ))
    w2single$projection <- combine.dist.compare(lapply(outputs, function(o) o$W2_dist$single$projection ))
    msesingle$selection <- combine.dist.compare(lapply(outputs, function(o) o$mse$single$selection ))
    msesingle$projection <- combine.dist.compare(lapply(outputs, function(o) o$mse$single$projection ))

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
    class(w2_opt$selection) <- class(w2_opt$projection) <- class(w2single$selection)

    w2_opt$combined <- lapply(names(w2_opt$selection)[1:2], function(nn) do.call("rbind", list(w2_opt$selection[[nn]], w2_opt$projection[[nn]])))
    w2_opt$combined[[3]] <- w2_opt$selection$p
    names(w2_opt$combined) <- c("posterior", "mean","p")
    class(w2_opt$combined) <- class(w2single$selection)

    # mse
    mse_opt$selection <- lapply(msesingle$selection[1:2], function(res) res[res$groups %in% c("Binary Programming", "Lasso"),])
    mse_opt$selection$p <- msesingle$selection$p
    mse_opt$projection <- lapply(msesingle$projection[1:2], function(res) res[res$groups %in% c("L1","Lasso","LInf"),])
    mse_opt$projection$p <- msesingle$projection$p

    for(pp in names(mse_opt)) {
      for(nn in names(msesingle$selection)[1:2]) {
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

    class(mse_opt$selection) <- class(mse_opt$projection) <- class(msesingle$selection)

    mse_opt$combined <- lapply(names(mse_opt$selection)[1:2], function(nn) do.call("rbind", list(mse_opt$selection[[nn]], mse_opt$projection[[nn]])))
    mse_opt$combined[[3]] <- mse_opt$selection$p

    names(mse_opt$combined) <- c("posterior", "mean","p")
    class(mse_opt$combined) <- class(msesingle$selection)

    # debugonce(plot.combine.dist.compare)
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
    w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$expectation$single$selection)))
    w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$expectation$single$projection)))

    w2r2single$selection  <- combine.WPR2(lapply(outputs, function(o) (o$W2_r2$null$single$selection %>% filter(nactive >0))))
    w2r2single$projection <- combine.WPR2(lapply(outputs, function(o) (o$W2_r2$null$single$projection %>% filter(nactive >0))))
    # w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$null$single$selection %>% filter(nactive >0))))
    # w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) (o$W1_r2$null$single$projection %>% filter(nactive >0))))

    # plot(w2r2single$selection)
    # plot(w2r2single$projection)
    # browser()
    w2r2single$selection  <- combine.WPR2(lapply(outputs, function(o) {
      temp <- o$W2_r2$null$single$selection %>%
        filter(nactive > 0)
      dists <- o$W2_dist$single$selection$mean %>%
        filter(nactive > 0)
      old_base <- dists$dist
      new_base <- dists %>%
        # group_by(groups) %>% mutate(base = dist[which.min(nactive)])
      group_by(groups) %>% mutate(base = max(dist))
      # if(nrow(new_base) != length(old_base)) browser()
      temp$r2 <- 1- old_base^2 / new_base$base^2
      temp$base <- "dist.from.null"
      temp <- temp[,c("r2","nactive","groups","method","p","base")]
      return(temp)
      }))
    w2r2single$projection <- combine.WPR2(lapply(outputs, function(o) {
      temp <- o$W2_r2$null$single$projection %>% filter(nactive > 0)
      dists <- o$W2_dist$single$projection$mean %>%
        filter(nactive > 0)
      old_base <-  dists$dist
      new_base <- dists %>%
        # group_by(groups) %>% mutate(base = dist[which.min(nactive)])
        group_by(groups) %>% mutate(base = max(dist))
      # if(nrow(new_base) != length(old_base)) browser()
      temp$r2 <- 1- old_base^2 / new_base$base^2
      temp$base <- "dist.from.null"
      return(temp)
    }))
    # w1r2single$selection  <- combine.WPR2(lapply(outputs, function(o) {
    #   temp <- o$W1_r2$null$single$selection
    #   old_base <- o$W2_dist$single$selection$mean$dist
    #   new_base <- o$W2_dist$single$selection$mean %>%
    #     group_by(groups) %>% mutate(base = dist[which.min(nactive)])
    #   temp$r2 <- 1- old_base / new_base$base
    #   return(temp)
    # }))
    # w1r2single$projection <- combine.WPR2(lapply(outputs, function(o) {
    #   temp <- o$W1_r2$null$single$projection
    #   old_base <- o$W2_dist$single$projection$mean$dist
    #   new_base <- o$W2_dist$single$selection$mean %>%
    #     group_by(groups) %>% summarise(base = dist[which.min(nactive)])
    #   temp$r2 <- 1- old_base^2 / new_base$base^2
    #   return(temp)
    # }))

    # w

    w2r2_list[[corr]] <- combine.WPR2(w2r2single$selection %>% filter(groups %in% c("Binary Programming",  "Lasso")) %>%
                                        mutate(groups = fct_recode(groups, "Relaxed B.P." = "Lasso")) %>%
                                        mutate(groups = fct_drop(groups)),
                                      w2r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
                                        mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
                                        mutate(groups = fct_drop(groups)))

    w1r2_list[[corr]] <- combine.WPR2(w1r2single$selection %>% filter(groups %in% c("Binary Programming",  "Lasso")) %>%
                                        mutate(groups = fct_recode(groups, "Relaxed B.P." = "Lasso")) %>%
                                        mutate(groups = fct_drop(groups)),
                                      w1r2single$projection %>% filter(groups %in% c("L1",  "Lasso", "LInf")) %>%
                                        mutate(groups = fct_recode(groups, L2 = "Lasso")) %>%
                                        mutate(groups = fct_drop(groups)))
    w1r2_list[[corr]]$corr <- corr
    w2r2_list[[corr]]$corr <- corr

    #### derivatives ####
    if(family == 'binomial') {
      # browser()
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

                                w2.dist <- sapply(m, function(ll) sapply(ll, function(tt) limbs::wasserstein(X = tt[-1,,drop= FALSE], Y = nn, p = 2,
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
      d.dist <- list(posterior = do.call("rbind",lapply(d.dist.list, function(d) d$boot)) %>%  filter(nactive >0),
                     mean = do.call("rbind",lapply(d.dist.list, function(d) d$w2))  %>%  filter(nactive >0), p = w2_list[[1]]$p-1)
      class(d.dist) <- class(w2_list[[1]])
      deriv_list[[corr]] <- d.dist
      # plot(d.dist, alpha = 0.5, base_size = basesize, ylab = "MSE between coef. and deriv.", facet.group = "corr", CI = "none")
      # theta  <- lapply(proj, function(p) p)

    }
  }

  w2_stack <- mse_stack <- deriv_stack <- w2_deriv_stack <- vector("list", 2)
  names(w2_stack) <- names(mse_stack) <-
    names(deriv_stack) <- c("posterior", "mean")
  # browser()
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
    mse_stack$mean <- mse_stack$mean %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    w2_stack$mean <-  w2_stack$mean  %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    deriv_stack$mean <- deriv_stack$mean %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.",
                                                                    "Stepwise","Simulated Annealing","Hahn-Carvalho")))
    deriv_stack$posterior <- deriv_stack$posterior %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.",
                                                                          "Stepwise","Simulated Annealing","Hahn-Carvalho")))
    class(deriv_stack) <-  class(w2_list[[1]])
  } else {
    mse_stack$posterior <- mse_stack$posterior %>%
      group_by(corr, method, groups, iter, nactive) %>%
      mutate(dist = min(dist)) %>%
      select(-c(ranks)) %>%
      distinct() %>%
      mutate(ranks = NA) %>%
      relocate(dist, nactive, groups, method, ranks, iter, corr) #%>%
      # group_by(corr, method, iter) %>%
      # summarise(dist = dist/min(dist), nactive = nactive,
      #           ranks = NA, groups = groups) %>%
      # relocate(dist, nactive, groups, method, ranks, iter, corr)
    mse_stack$mean <- mse_stack$mean %>%
      group_by(corr, method, groups, iter, nactive) %>%
      mutate(dist = min(dist)) %>%
      select(-c(ranks)) %>%
      distinct() %>%
      mutate(ranks = NA) %>%
      relocate(dist, nactive, groups, method, ranks, iter, corr)
  }
  class(mse_stack) <- class(w2_stack) <- class(w2_list[[1]])

  msep <- plot.combine.dist.compare(mse_stack, alpha = 0.5, base_size = basesize, ylab = "Relative MSE", facet.group = "corr", CI = "none")
  w2p <- plot.combine.dist.compare(w2_stack,  alpha = 0.5, base_size = basesize, ylab = "2-Wasserstein Distance", facet.group = "corr", CI = "none")

  if(family == "binomial") {
    msep$mean <- msep$mean + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    w2p$mean  <- w2p$mean  + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    deriv_stack_p <- plot.combine.dist.compare(deriv_stack, alpha = 0.5, base_size = basesize, ylab = "Relative MSE", facet.group = "corr", CI = "none")

    deriv_stack_p$mean <-  deriv_stack_p$mean + ggplot2::ylab("2-Wasserstein Distance") + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])

    deriv_stack_p$posterior <-  deriv_stack_p$posterior + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
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
  # msep$mean
  # dev.off()
  #
  # pdf(file = "inst/figure/simulation/w2_normal_neighb.pdf", width = 7, height = 5)
  # w2p$mean
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
    w2fn <- file.path(directory, paste0(paste0(c(neighb,family,penalty , pf ,method, "w2plot"), collapse="_"), ".rds"))
    saveRDS(list(mse = msep, w2 = w2p, deriv = deriv_stack_p), file = w2fn)
  } else {
    directory <- file.path("inst", "figure","rawplots")
    w2fn <- file.path(directory, paste0(paste0(c(neighb,family,penalty , pf ,method, "w2plot"), collapse="_"), ".rds"))
    saveRDS(list(mse = msep, w2 = w2p, deriv = NULL), file = w2fn)
  }
  # msefile <- file.path(directory, paste0(paste0(c(neighb,family,penalty , pf ,method, "mseplot"), collapse="_"), ".rds"))
  # pdf(file, width = width, height = height)
  # msep$mean
  # dev.off()
  #
  # w2file <- file.path(directory, paste0(paste0(c(neighb,family,penalty , pf ,method, "w2plot"), collapse="_"), ".rds"))
  # pdf(file, width = width, height = height)
  # w2p$mean
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
  # msep$mean
  # w2p$mean
  # plot(w2r2, alpha = 0.5, base_size = 20)
  # plot(w1r2, alpha = 0.5, base_size = 20)
  # dev.off()
  w1r2 <- do.call("rbind", w1r2_list)
  w2r2 <- do.call("rbind", w2r2_list)

  if(family == "binomial") {
    w1r2 <- w1r2 %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    w2r2 <- w2r2 %>% filter(!(groups %in% c("Binary Programming", "Relaxed B.P.")))
    w2r2deriv <- deriv_stack$mean %>% group_by(corr,groups, method, iter) %>%
      summarise(r2 = 1-dist^2/max(dist)^2, nactive = nactive, p = 2, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      summarise(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    w1r2deriv <- deriv_stack$mean %>% group_by(corr,groups, method, iter) %>%
      summarise(r2 = 1-dist/max(dist), nactive = nactive, p = 1, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      summarise(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
   class(w2r2deriv) <- class(w1r2deriv) <- class(w2r2)
  }
  # browser()
  # w2r2.test <- w2_stack$mean %>% group_by(groups, corr, nactive) %>%
  #   summarise(dist = mean(dist)) %>%
  #   group_by(groups, corr) %>%
  #   mutate(r2 = 1 - dist^2/max(dist^2))
  # class(w2r2.test) <-class(w1r2)
  w1r2p <- plot.WPR2(w1r2, alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE)
  w2r2p <- plot.WPR2(w2r2 #%>% filter(groups != "Binary Programming")
                     , alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE,
                     ylim = c(0,1)) + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[2:5])

  # w1r2post <- w1r2

  # temp <- o$W2_r2$null$single$projection %>% filter(nactive > 0)
  # dists <- o$W2_dist$single$projection$mean %>%
  #   filter(nactive > 0)
  # old_base <-  dists$dist
  # new_base <- dists %>%
  #   # group_by(groups) %>% mutate(base = dist[which.min(nactive)])
  #   group_by(groups) %>% mutate(base = max(dist))
  # # if(nrow(new_base) != length(old_base)) browser()
  # temp$r2 <- 1- old_base^2 / new_base$base^2
  # temp$base <- "dist.from.null"
  w2r2derivp <- NULL
  if(family == "binomial") {
    w1r2p <- w1r2p + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    w2r2p <-  w2r2p + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[3:5]) +
      ggplot2::scale_fill_manual(values = ggsci::pal_jama()(5)[3:5])
    w1r2derivp <- plot.WPR2(w1r2deriv, alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE)+
      theme(strip.background = element_blank(),
      strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
    w2r2derivp <- plot.WPR2(w2r2deriv #%>% filter(groups != "Binary Programming")
                            , alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE,
                            ylim = c(0,1)) + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[2:5]) +
      theme(strip.background = element_blank(),
              strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
    w2r2derivp <- list(posterior = NULL, mean = w2r2derivp)
  }
  w2r2p <- w2r2p + theme(strip.background = element_blank(),
                         strip.text.x = element_blank()) +
    scale_x_continuous(breaks = c(1, seq(5,25,5)))
  w1r2p <- w1r2p + theme(strip.background = element_blank(),
                         strip.text.x = element_blank()) +
    scale_x_continuous(breaks = c(1, seq(5,25,5)))

  w1r2postp <- w2r2postp <- NULL
  if(!is.null(w2_stack$post)) {
    w1r2post <- w2_stack$post %>% group_by(corr,groups, method, iter) %>%
      summarise(r2 = 1-dist/max(dist), nactive = nactive, p = 1, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      summarise(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    w2r2post <- w2_stack$post %>% group_by(corr,groups, method, iter) %>%
      summarise(r2 = 1-dist^2/max(dist)^2, nactive = nactive, p = 2, base = "dist.from.null") %>%
      group_by(corr,groups, method, nactive,  p, base) %>%
      summarise(r2 = mean(r2)) %>%
      relocate(r2, nactive, groups, method, p, base, corr)
    class(w1r2post) <- class(w2r2post) <- class(w2r2)
    w1r2postp <- plot.WPR2(w1r2post, alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
    w2r2postp <- plot.WPR2(w2r2post #%>% filter(groups != "Binary Programming")
                           , alpha = 0.5, base_size = basesize, facet.group = "corr", ribbon = TRUE,
                           ylim = c(0,1)) + ggplot2::scale_color_manual(values = ggsci::pal_jama()(5)[2:5]) +
      theme(strip.background = element_blank(),
            strip.text.x = element_blank()) +
      scale_x_continuous(breaks = c(1, seq(5,25,5)))
  }


  w2r2p <- list(posterior = w2r2postp, mean = w2r2p)
  w1r2p <- list(posterior = w1r2postp, mean = w1r2p)

  wpr2fn <- file.path("inst", "figure","rawplots", paste0(paste0(c(neighb, family,penalty , pf ,method, "wpr2plot"), collapse="_"), ".rds"))
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
