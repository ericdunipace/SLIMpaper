require(SLIMpaper)
library(dplyr)
library(forcats)
library(ggplot2)
library(gridExtra)

cube.root <- function(x) x^(1/3)
cube.root_trans <- function(){
  scales::trans_new("cube.root", function(x) x^(1/3), function(x)x^3,
                    domain = c(0, Inf))
}

combine_plot_temp <- function(family, which.plt = c("mse","w2", "deriv","w1r2","w2r2","w2r2deriv"),
                              height = 5, width = 7, breaks = NULL, trans = NULL, rm.layer = NULL,
                              limits = list(NULL),
                              # rel.heights = NULL,
                              # mean = c("mse","w2", "deriv","w1r2","w2r2"),
                              posterior = NULL,
                              sep.posterior = FALSE,
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
      fn1 <- file.path("inst", "figure","rawplots", paste0(paste0(c(cur_neighb,family,penalty , pf ,method, "w2plot"), collapse="_"), ".rds"))
      plt1 <- readRDS(fn1)
    }
    if (any(which.plt %in% c("w1r2","w2r2", "w2r2deriv"))) {
      fn2 <- file.path("inst", "figure","rawplots", paste0(paste0(c(cur_neighb, family,penalty , pf ,method, "wpr2plot"), collapse="_"), ".rds"))
      plt2 <- readRDS(fn2)
    }
    plts <- c(plt1,
              plt2 )

    plts <- plts[which.plt]
    if(!is.null(posterior)) {
      stopifnot(posterior %in% which.plt)
      runPosterior <- TRUE
    } else {
      runPosterior <- FALSE
    }
    # }
    # plts[[cur_neighb]] <- do.call("rbind", )
    outfile <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt), collapse="_"), ".pdf"))
    # if(!is.null(posterior)) outfile <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt, "posterior"), collapse="_"), ".pdf"))
    if(sep.posterior) {
      outfile.mean <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt), collapse="_"), "_mean.pdf"))
      outfile.post <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt), collapse="_"), "_post.pdf"))

    }

    print.p.list <- vector("list", length(plts))
    names(print.p.list) <- names(plts)
    # browser()
    for(i in names(plts)) {
      pp <- list()
      if(i %in% posterior) {
        idx <- c("mean","posterior")
      } else {
        idx <- "mean"
      }
      if(i == "deriv") idx <- idx[2:1]
      for (j in idx) {
        express <- if( j == "posterior" &
                      runPosterior & !is.null(plts[[i]][[j]]) &
                      i != "deriv") {
          ylab(expression(W[2](beta, theta)))
        } else if (i == "deriv" & j != "posterior" & !is.null(plts[[i]][[j]])) {
          ylab(expression(W[2](nabla[x]~mu, nabla[x]~nu)))
        } else {
          ylab(expression(W[2](mu, nu)))
        }
        pp[[j]] <- if( (i == "deriv" & j == "mean") | i == "w2") {
          plts[[i]][[j]] + express
        # } else if (i != "w2r2" & i != "w1r2") {
        #   plts[[i]][[j]]
        # } else if (i %in% posterior) {

        } else {
          plts[[i]][[j]]
        }

        if(!is.null(trans[[i]][[j]]) & !isFALSE(trans[[i]][[j]])) {
          if(isTRUE(trans[[i]][[j]])) trans[[i]][[j]] <- "sqrt"
          if(is.null(breaks[[i]][[j]])) breaks[[i]][[j]] <- waiver()
          pp[[j]] <- pp[[j]] + scale_y_continuous(trans = trans[[i]][[j]],
                                        expand = c(0.01,0.05),
                                        breaks = breaks[[i]][[j]],
                                        limits = limits[[i]][[j]]) +
            expand_limits(y = 0)
        }
        if(!is.null(rm.layer)) pp[[j]] <- remove_geoms(pp[[j]], rm.layer, FALSE)

        if(!is.null(labels) & !is.null(pal)) {
          pp[[j]] <- pp[[j]] +
            scale_color_manual(name = "Method:",
                    breaks = levels(pp[[j]]$data$groups), #c("L1", "L2", "LInf"),
                    labels = labels,
                    values = pal) +
            scale_fill_manual(name = "Method:",
                              breaks = levels(pp[[j]]$data$groups), #c("L1", "L2", "LInf"),
                              labels = labels,
                              values = pal) +
            scale_size_manual(name = "Method:",
                              breaks = levels(pp[[i]][[j]]$data$groups), #c("L1", "L2", "LInf"),
                              labels = labels,
                              values = pal)
        }
      }


      print.p.list[[i]] <- pp

    }
    # browser()
    print.p <- list()
    print.p$mean <- unlist(print.p.list, recursive = FALSE)
    if(sep.posterior) {
      temp <- unlist(print.p.list[names(plts) %in% posterior], recursive = FALSE)
      print.p$post <- temp[grepl("posterior", names(temp))]
      print.p$mean <- print.p$mean[!grepl("posterior", names(print.p$mean))]
    }

    plot.legend <- get_legend(print.p$mean[[1]] + ggplot2::theme(legend.position="bottom"))
    for(j in names(print.p)) {
        for(i in seq_along(print.p[[j]])) {
          if(i == 1 & j == "mean") {
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
          if(i != length(print.p[[j]])) {
            print.p[[j]][[i]] <- print.p[[j]][[i]] +
              theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank())
          }
      }
      }
    # plot.length <- length(print.p)
    # if(is.null(rel.heights)) rel.heights <- rep(3, plot.length)
    if(sep.posterior) {
      pred <- arrangeGrob(do.call("rbind", lapply(print.p$mean, ggplotGrob)),
                  bottom = grid::textGrob("Number of active coefficients",
                                          vjust = -1.8, hjust=0.4))
      post <- arrangeGrob(do.call("rbind", lapply(print.p$post, ggplotGrob)),
                          bottom = grid::textGrob("Number of active coefficients",
                                                  vjust = -1.8, hjust=0.4))
      pdf(outfile.mean, width = width, height = height)
      print(grid.arrange(pred,
                         plot.legend, nrow = 2, heights = c(10,.5)))
      dev.off()
      pdf(outfile.post, width = width, height = height)
      print(grid.arrange(post,
                         plot.legend, nrow = 2, heights = c(10,.5)))
      dev.off()

    } else {
      sim.plots <- arrangeGrob(do.call("rbind", lapply(print.p$mean, ggplotGrob)),
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

# file.path("inst", "figure","simulations", paste0(paste0(c(cur_neighb,family,penalty , pf ,method, cur_corr, "w2plot"), collapse="_"), ".rds"))
# debugonce(combine_plot_temp)
combine_plot_temp("binomial", which.plt = c("deriv","w2r2deriv"),   height = 6, width = 7.5,
                  trans = list(deriv = list(mean = "sqrt", posterior = "sqrt")),
                  breaks = list(deriv = list(mean = c(0, 0.25, 1, 4),
                                             posterior = c(0,0.25, 1,4))),
                  # rel.heights = c(3,5),
                  limits = list(deriv = list(mean = NULL,
                                             posterior = c(0,4))),
                   labels = expression(W[1], W[2], W[infinity]),
                   pal = ggsci::pal_jama()(5)[3:5], posterior = c("deriv"))
combine_plot_temp("gaussian", c("mse","w2","w2r2"),   height = 6, width = 7.5,
                  breaks = list(mse = list(mean = c(0, 1, 8, 27, 64),
                                           posterior = c(0, 1, 4, 16, 36, 64, 100)),
                                w2 = list(mean = c(0, 0.25, 1),
                                          posterior = c(0, 0.25, 1, 2))),
                  rm.layer = "GeomRibbon",
                  trans = list(mse = list(mean = "cube.root",
                                          posterior = "sqrt"),
                               w2 = list(mean= "sqrt",
                                         posterior = "sqrt")),
                  labels = expression("B.P.", "Relaxed B.P.", W[1], W[2],
                                      W[infinity]), pal = ggsci::pal_jama()(5),
                  posterior = c("mse","w2","w2r2"),
                  sep.posterior = TRUE)
# combine_plot_temp("gaussian", "w2",   height = 2, width = 7.5, posterior = TRUE, rm.layer = "GeomRibbon",
#                   breaks = list(mse = c(0,0.1,0.2,0.3),
#                                 w2 = c(0:3)),
#                   labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity]), pal = ggsci::pal_jama()(5))
# combine_plot_temp("binomial", "wpr2", height = 2, width = 7.5, rm.layer = "GeomRibbon",
#                   labels = expression(W[1], W[2], W[infinity]), pal = ggsci::pal_jama()(5)[3:5])
# combine_plot_temp("gaussian", "wpr2", height = 2, width = 7.5, rm.layer = "GeomRibbon",
#                   labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity]), pal = ggsci::pal_jama()(5))
