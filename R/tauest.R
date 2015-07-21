#' Estimation of Kendall's tau and asymptotic covariance matrix
#'
#' \code{tauest} computes an estimate of Kendall's tau using \code{pcaPP::cor.fk}.
#' Optionally the asymptotic covariance matrix of this estimator is computed as well
#'
#' @param x data matrix containing n observations of the d-dimensional random variables.
#' @param asymp	if \code{asymp = TRUE}, the asymptotic covariance matrix of lower-diagonal
#' elements of the copula correlation estimator is calculated. Default = FALSE
#' @param status if \code{status = TRUE}, the progress of calculation is printed.
#'
#' @export
#' @return list containing the estimated Kendall's tau matrix and in case \code{asymp = TRUE}
#' also the asymptotic covariance matrix


 tauest <- function(x, asymp = FALSE, status = FALSE){
   N <- dim(x)[1]
   D <- dim(x)[2]

   ret <- list()
   ret$tau <- pcaPP::cor.fk(x)
   ret$tauVec <- ret$tau[lower.tri(ret$tau)]

   if(asymp){
     DL <- D * (D - 1) / 2
     TAU4E <- rep(0, DL^2)

     out <- .C("tau4", data = as.double(x), n = as.integer(N),
              d = as.integer(D), tauest4 = as.double(TAU4E),
              status = as.integer(status))

     TAU4E <- matrix(out$tauest4,DL,DL)
     ret$tau4 <- TAU4E + t(TAU4E)
     diag(ret$tau4) <- diag(ret$tau4) / 2
     }
  return(ret)
 }

