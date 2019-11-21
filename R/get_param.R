get_param <- function() {
  # gaussian_param <- function(corr, sigma2, ratio, ...) {
  #   # var_y <- sigma2/ratio
  #
  #   # ecorr *sigma2 +
  #
  #
  # }
  theta_gen <- function() {
    theta_star <- c(1.7, #intercept
                    1.517764, 1.864788, 1.747564, 1.984903, 1.313780,
                    -c(1.386298, 1.471564, 1.277228, 1.786211, 1.169712),
                    0.28585294, 0.02419794, 0.14508495, 0.14131086, 0.12519745,
                    -0.03484582, -0.32138290, -0.14185643, -0.10779317, -0.38835672)
    theta_not_intercept <- theta_star[-1]
    theta_norm <- sum(theta_not_intercept^2 )
    theta_star[-1] <-  theta_not_intercept / sqrt( theta_norm )
    return(theta_star)
  }
  gaussian_param <- function() {

    theta_star <- theta_gen()
    sigma2_star <- 1
    # gives sigma^2/var(y) = 0.25
    return( list( theta = theta_star,
                  sigma2 = sigma2_star) )

  }

  binary_param <- function() {
    # theta_star <- c(-0.7,
    #                 .14, .12, .1, .14, .18,
    #                 -.17, -.12, -.14, -.18, -.15,
    #                 0.02, 0.05, 0.009, .018, .01,
    #                 .029, .04, .005, .04, .006)
    theta_star <- theta_gen()
    theta_star[1] <- -0.7
    return(theta_star)
  }

  expo_param <- function() {
    # theta_star <- c(-1.6,
    #                 .14, .12, .1, .14, .18,
    #                 -0.17, -0.12, -0.14, -0.18, -0.15,
    #                 0.02, 0.05, 0.009, 0.018, 0.01,
    #                 0.029, 0.04, 0.005, 0.04, 0.006)
    theta_star <- theta_gen()
    return(theta_star)
  }

  return(list(gaussian = gaussian_param,
              binary = binary_param,
              exponential = expo_param))
}
