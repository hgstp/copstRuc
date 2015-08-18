#' Computes the asymptotic distribution free test statistic
#'
#' @param x object of class copstruc
#' @param m number of expected factors
#' @param weight_matrix_adf method to compute weight matrix
#' @param maxit_optim the maximum number of iterations in optim
#' @export
#'
#' @return object of class asymp_distr_free
#'

 asymp_distr_free <- function(x,  m = 1, weight_matrix_adf = "regular", ...){

  Rvec <- x$Rvec
  Ga <- x$Ga
  n <- x$obs
  d <- dim(x$R)[1]
  m <- as.integer(m)


  adf <- ADF(Rvec, Ga, n, d, m, weight_matrix_adf = weight_matrix_adf, ... )

  adf$call <- match.call()
  adf$nFac <- m
  adf$dim <- d

  class(adf) <- "asymp_distr_free"
  adf

 }


 #------------------------

 print.asymp_distr_free <- function(x, ... ){

   cat("Call:\n")
   print(x$call)
   cat("\nADF test statistic:\n")
   print(x$nF)
 }

 #------------------------

 summary.asymp_distr_free = function(object, ...){
    tab <- cbind(test_statistic = object$nF[1],
                 dimension = object$dim,
                 number_of_factors = object$nFac,
                 df = object$dof,
                 p_value = object$pValue[1])

   ret <- list(call = object$call, results = tab,
               factor_loadings = object$load,
               specific_loadings = object$uniq,
               weight_matrix = object$wm )

   return(ret)
 }

#------------------------
 print.summary.asymp_distr_free = function(object, ... ){
   cat("Call:\n")
   print(object$call)
   cat("\n")
   printCoefmat(object$results, P.values = TRUE, has.Pvalue = TRUE)
   cat("\n")
   printCoefmat(cbind(object$factor_loadings, object$specific_loadings))
   cat("\n")
   print(object$weight_matrix)
 }




