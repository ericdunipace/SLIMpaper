get_param <- function() {
  # gaussian_param <- function(corr, sigma2, ratio, ...) {
  #   # var_y <- sigma2/ratio
  #
  #   # ecorr *sigma2 +
  #
  #
  # }
  gaussian_param <- function() {

    theta_star <- c(1.7, #intercept
                    1.4, 1.2, 1, 1.4, 1.8,
                    -1.7, -1.2, -1.4, -1.8, -1.5,
                    0.2, 0.5, 0.09, .18, .1,
                    .29, .4, .05, .4, .06)
    sigma_star <- 1

    return(list(theta= theta_star, sigma2=sigma_star))

  }

  binary_param <- function() {
    theta_star <- c(-0.7,
                    .14, .12, .1, .14, .18,
                    -.17, -.12, -.14, -.18, -.15,
                    0.02, 0.05, 0.009, .018, .01,
                    .029, .04, .005, .04, .006)
    return(theta_star)
  }

  expo_param <- function() {
    theta_star <- c(-1.6,
                    .14, .12, .1, .14, .18,
                    -0.17, -0.12, -0.14, -0.18, -0.15,
                    0.02, 0.05, 0.009, 0.018, 0.01,
                    0.029, 0.04, 0.005, 0.04, 0.006)
    return(theta_star)
  }

  return(list(gaussian = gaussian_param,
              binary = binary_param,
              exponential = expo_param))
}
