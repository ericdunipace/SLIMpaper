# combine plot legends for gridExtra
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
get_legend<-function(plot, location = c("bottom","left", "right", "top")){
  loc <- match.arg(location)
  # tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot))
  # leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  # legend <- tmp$grobs[[leg]]
  # return(legend)
  g <- ggplot2::ggplotGrob(plot + ggplot2::theme(legend.position=loc))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  return(legend)
}

# grid_arrange_shared_legend <- function(...) {
#   plots <- list(...)
#   g <- ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position="bottom"))$grobs
#   legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
#   lheight <- sum(legend$height)
#   gridExtra::grid.arrange(
#     do.call(arrangeGrob, lapply(plots, function(x)
#       x + ggplot2::theme(legend.position="none"))),
#     legend,
#     ncol = 1,
#     heights = grid::unit.c(ggplot2::unit(1, "npc") - lheight, lheight))
# }

remove_geoms <- function(x, geom_type, last_only = TRUE) {
  # Find layers that match the requested type.
  selector <- sapply(x$layers,
                     function(y) {
                       class(y$geom)[1] == geom_type
                     })
  if(last_only)
    selector <- max(which(selector))
  # Delete the layers.
  x$layers[selector] <- NULL
  x
}
