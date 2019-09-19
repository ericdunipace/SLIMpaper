.checkSP <- function() {
  if(!("SparsePosterior" %in% utils::installed.packages())) {
    devtools::install_github("ericdunipace/SparsePosterior", ref="master")
  }
}
