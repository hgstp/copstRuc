#' plot: k against estimate of nu
#'
#' plots the number of upper order statistics against \eqn{\widehat \nu}. k is varied
#' between 5 an 20% of the number of observations.
#'
#' @param x data matrix containing n observations of the d-dimensional random variables
#' @param prec	number of samples, which is used for smoothing the results of the estimation. Default = 31
#'
#' @export


 k_selection = function(x, prec = 31){
   n <- dim(x)[1]
   lower_k <- round(n * .05)
   upper_k <- round(n * .2)

   k <- seq(from = lower_k, to = upper_k, by = 10)
   nu_k <- seq_along(k)

   for(i in seq_along(k)){
     nu_k[i] <- rhoTailEst(x, k = k[i], prec = prec)$nu
   }

   plot(k, nu_k, type = "l", ylab = expression(paste("estimated ", nu)))
 }
