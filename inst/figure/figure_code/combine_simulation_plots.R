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
    if(!is.null(parameters)) {
      stopifnot(parameters %in% which.plt)
      runparameters <- TRUE
    } else {
      runparameters <- FALSE
    }
    # }
    # plts[[cur_neighb]] <- do.call("rbind", )
    outfile <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt), collapse="_"), ".pdf"))
    # if(!is.null(parameters)) outfile <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt, "parameters"), collapse="_"), ".pdf"))
    if (sep.parameters) {
      outfile.predictions <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt), collapse="_"), "_predictions.pdf"))
      outfile.parameters <- file.path("inst", "figure","simulation", paste0(paste0(c(cur_neighb, family, penalty, pf, method, which.plt), collapse="_"), "_parameters.pdf"))

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

# file.path("inst", "figure","simulations", paste0(paste0(c(cur_neighb,family,penalty , pf ,method, cur_corr, "w2plot"), collapse="_"), ".rds"))
# debugonce(combine_plot_temp)
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

# debugonce(combine_plot_temp)
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
# combine_plot_temp("gaussian", "w2",   height = 2, width = 7.5, parameters = TRUE, rm.layer = "GeomRibbon",
#                   breaks = list(mse = c(0,0.1,0.2,0.3),
#                                 w2 = c(0:3)),
#                   labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity]), pal = ggsci::pal_jama()(5))
# combine_plot_temp("binomial", "wpr2", height = 2, width = 7.5, rm.layer = "GeomRibbon",
#                   labels = expression(W[1], W[2], W[infinity]), pal = ggsci::pal_jama()(5)[3:5])
# combine_plot_temp("gaussian", "wpr2", height = 2, width = 7.5, rm.layer = "GeomRibbon",
#                   labels = expression("B.P.", "Relaxed B.P.", W[1], W[2], W[infinity]), pal = ggsci::pal_jama()(5))
