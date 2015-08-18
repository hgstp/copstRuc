#' Copula structure analysis
#'
#' \code{copstruc} performs a regular, robust or extreme copula structure
#' analysis. Currently there is only an assumed factor analysis structure
#' implemented
#'
#' @param x data matrix containing n observations of the d-dimensional random variables.
#' @param method name of the estimaton method, which can be \code{tail}, \code{kendall}
#' or \code{standard}.
#' @param k	number of the largest order statistics, which is used for estimating the
#' tail dependence function. Default = 1.
#' @param k_select logical to decide if copstruc should be used to approximate \code{k}
#' @param asymp	if \code{asymp = TRUE}, the asymptotic covariance matrix of lower-diagonal
#' elements of the copula correlation estimator is calculated. Default = FALSE
#'
#' @return object of class copstruc, which is specified by:
#'
#'   \itemize{
#'      \item \code{R},	the estimated correlation matrix
#'      \item \code{Rvec}, the lower diagonal elements of R
#'      \item \code{Ga}, the asymptotic covariance of Rvec (if \code{asymp == TRUE})
#'      \item \code{nu}, the estimated coefficient of regular variation (if \code{method == 'tail'})
#'      \item \code{method}, the estimation method being used
#'      }
#' @export
#'



copstruc <- function(x, method = "tail", k = 100, asymp = FALSE, k_select = FALSE, ... ){
  x <- as.matrix(x)

  ret <- copstruc_est(x, method = method, k = k, asymp = asymp, k_select = k_select, ... )
  ret$call <- match.call()

  class(ret) <- "copstruc"
  return(ret)
}



print.copstruc <- function(x, ... ){
  cat("Call:\n")
  print(x$call)
  cat("\nNames of copstruc results:\n")
  print(names(x))
  cat("\nMethod for copula correlation estimator:\n")
  print(x$method)
  cat("\nAsymptotic covariance matrix available:\n")
  print(x$asymp)
}








