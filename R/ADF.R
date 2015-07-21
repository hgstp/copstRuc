# ADF computes the asymptotic distribution free test statistic

# arguments:
#    	Rvec:	  	lower diagonal elements of the estimated correlation matrix.
#    	Ga:  			asmptotic covariance matrix of Rvec.
#    	n:  			number of observations, used for the estimation in cop.struc
#    	d:  			dimension of the observation variables.
#    	m:  			number of expected factors.
#	weight_matrix_adf:	method to compute weight matrix
#	maxit_optim:		the maximum number of iterations in optim
#
# value:
#  dof:  	degree of freedom
#  nFml:  	calculated test statistics by using the estimated loadings from factanal.
#  nFP:  	a modified test statistic, which is estimated by minimizing over a discrepancy function with a certain penalty.
#  nF:  	test statistic using optimized factor loadings, estimated with optim.
#  mlStat:  	test statistic delivered by factanal.
#  con:  	convergence of optim, which is an integer code:
#
#			0: indicates successful completion.
#			1: indicates that the iteration limit maxit had been reached.
#			10: indicates degeneracy of the Nelder-Mead simplex.
#			51: indicates a warning from the L-BFGS-B method.
#			52: indicates an error from the L-BFGS-B method.


ADF  <-  function(Rvec, Ga, n, d, m, weight_matrix_adf = "regular", maxit_optim = 500){

  ret <- list()
  ret$dof  <-  d *(d - 1) / 2 - (d * m - m * (m - 1) / 2)

  if(ret$dof < 1){
    print('Error, too many factors - negative dof!')
  }

  if(ret$dof >= 1){
    dimIndex <- subsets(d)
    chol_Gamma <- chol(Ga)

    chol_Gamma_inv <- chol2inv(chol_Gamma)

    CR <- matrix(0, d, d)
    CR[dimIndex] <- Rvec
    CR <- CR + t(CR) + diag(1, d)

    cl <- list()
    cl$center <- rep(0, d)
    cl$n.obs <- n
    cl$cov <- CR
    # classical factor analysis
    fa <- factanal(factors = m, covmat = cl, method = 'mle', nstart = 5,
                   rotation = 'none', control = list(opt = list(maxit = 100000)))
    # loadings are starting values
    Lstart <- as.vector(fa$loadings[1 : (d * m)])

    if(m > 1){
      # use varimax rotated loading if more than one factor is used
      Lstart <- as.vector(varimax(matrix(fa$loadings, d, m))$loadings[1 : (d * m)])
      }

    tmp <- optim(Lstart, FQD, Rvec = Rvec, GaInv = chol_Gamma_inv,
                 d_FQD = d, m_FQD = m, weight_matrix = weight_matrix_adf,
                 control = list(maxit = maxit_optim))

    rot <- matrix(tmp$par, d, m)
    # varimax rotation of the final loadings
    if(m > 1){
      rot <- matrix(varimax(matrix(tmp$par, d, m))$loadings[1 : (d * m)], d, m)
      }
    # test staistik with penalty
    ret$nF <- n * FQDnonPenal(rot, R = Rvec, GaInv = chol_Gamma_inv,
                              weight_matrix = weight_matrix_adf,
                              d_FQD = d, m_FQD = m)
    # test staistik without penalty
    ret$nFP <- n * FQD(rot, R = Rvec, GaInv = chol_Gamma_inv,
                       weight_matrix = weight_matrix_adf,
                       d_FQD = d, m_FQD = m)

    ret$conv <- tmp$convergence
    ret$load <- rot

    ret$uniq  <-  1 - diag(matrix(tmp$par, d, m) %*% t(matrix(tmp$par, d, m)))

    # p-value based on asymptotic chi-2-distribution
    ret$pValue  <-  pchisq(ret$nF, ret$dof, lower.tail = FALSE)

    ret$wm <- weight_matrix_adf

    return(ret)
  }
}
