nn_fun_path <- torch <- NULL

.onLoad <- function(libname, pkgname) {
  # np <<- reticulate::import("numpy", delay_load = TRUE)
  # scipy <<- reticulate::import("scipy", delay_load = TRUE)
  torch <<- reticulate::import("torch", delay_load = TRUE)
  nn_fun_path <<- file.path(system.file(package="CoarsePosteriorSummary"), "python")
}


.onUnload <- function (libpath) {
  library.dynam.unload("CoarsePosteriorSummary", libpath)
}
