check_python_modules <- function(method = "auto", conda = "auto") {

  if(!reticulate::py_module_available("scipy")) {
    reticulate::py_install("scipy", method = method, conda = conda)
  }

  if(!reticulate::py_module_available("numpy")) {
    reticulate::py_install("numpy", method = method, conda = conda)
  }

  if(!reticulate::py_module_available("torch")) {
    reticulate::py_install("torch", method = method, conda = conda)
  }

}
