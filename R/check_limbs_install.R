.checkLIMBS <- function() {
  if(!("limbs" %in% utils::installed.packages())) {
    devtools::install_github("ericdunipace/limbs", ref="master")
  }
}
