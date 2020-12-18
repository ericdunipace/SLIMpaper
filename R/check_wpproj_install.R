.checkWpProj <- function() {
  if (!("WpProj" %in% utils::installed.packages())) {
    devtools::install_github("ericdunipace/WpProj", ref = "master")
  }
}
