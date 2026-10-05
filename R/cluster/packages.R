#packages
if(!("RcppCGAL" %in% utils::installed.packages())) {
  devtools::install_github("ericdunipace/RcppCGAL", ref="master",force=TRUE)
}
if(!("limbs" %in% utils::installed.packages())) {
  devtools::install_github("ericdunipace/limbs", ref="master",force=TRUE)
} else {
    require(limbs)
}
if(!("SLIMpaper" %in% utils::installed.packages())) {
  devtools::install_github("ericdunipace/SLIMpaper", ref="master",force=TRUE)
} else {
  require(SLIMpaper)
}
# if(!("oem" %in% utils::installed.packages())){
#   install.packages("oem")
# } else {
#   require(oem)
# }
# if(!("mvtnorm" %in% utils::installed.packages())){
#   install.packages("mvtnorm")
# } else {
#   require(mvtnorm)
# }
# if(!("rstan" %in% utils::installed.packages())) {
#   install.packages("rstan")
#   rstan_options(auto_write = TRUE)
# } else {
#   require(rstan)
#   rstan_options(auto_write = TRUE)
# }
# if (!("rstanarm" %in% utils::installed.packages())) {
#   install.packages("rstanarm")
# } else {
#   require(rstanarm)
# }
# if (!("transport" %in% utils::installed.packages())) {
#   install.packages("transport")
# } else {
#   require(transport)
# }
# if (!("glmnet" %in% utils::installed.packages())) {
#   install.packages("glmnet")
# } else {
#   require(glmnet)
# }
