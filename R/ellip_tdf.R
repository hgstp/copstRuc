#' Bivariate elliptical tail dependence function
#'
#' ellip_tdf computes the taildependence function corresponding to a
#' bivariate elliptical copula with parameter \eqn{\nu} and \eqn{\rho}.
#' It uses the definition with respect to the t-distribution
#'
#' @param x first	argument of the bivariate tail dependence function
#' @param y second	argument of the bivariate tail dependence function
#' @param nu parameter of regular variation of the generating variable
#' @param rho copula correlation
#' @return \code{tde}	value of the tail dependence function
#' with parameters \code{nu} and \code{rho} at a point (\code{x},\code{y})
#' @export
#' @example
#' ellip_tdf(1,1)



ellip_tdf <- function(x, y, nu = 2, rho = 0.5){
  x <- pmax(x, 1e-50)
  y <- pmax(y, 1e-50)
  tde <- 	x * (1 - pt(sqrt((nu + 1) / (1 - rho^2)) * ((x / y)^(1 / nu) - rho), nu + 1)) +
    y * (1 - pt(sqrt((nu + 1) / (1 - rho^2)) * ((y / x)^(1 / nu) - rho), nu + 1))
}


