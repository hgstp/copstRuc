#' Estimation of Kendall's tau and asymptotic covariance matrix
#'
#' \code{rho_tau_est} computes an estimate of Kendall's tau using \code{pcaPP::cor.fk}.
#' Optionally the asymptotic covariance matrix of this estimator is computed as well
#'
#' @param x data matrix containing n observations of the d-dimensional random variables.
#' @param asymp	if \code{asymp = TRUE}, the asymptotic covariance matrix of lower-diagonal
#' elements of the copula correlation estimator is calculated. Default = FALSE
#'
#' @export
#' @return list containing the estimated Kendall's tau matrix and in case \code{asymp = TRUE}
#' also the asymptotic covariance matrix


 rho_tau_est <- function(x, asymp = FALSE){
   N <- dim(x)[1]
   D <- dim(x)[2]

   ret <- list()
   ret$tau <- pcaPP::cor.fk(x)
   ret$tauVec <- ret$tau[lower.tri(ret$tau)]

   if(asymp){
     DL <- D * (D - 1) / 2
     var_tau <- rep(0, DL^2)
     x_vectorised <- as.double(x)

     var_tau <- matrix(asymp_var_tau(x_vectorised, D, N), DL, DL)
     ret$tau4 <- var_tau + t(var_tau)
     diag(ret$tau4) <- diag(ret$tau4) / 2
     }
  return(ret)
 }

