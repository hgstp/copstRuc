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
  subsets  <-  function(n){
    # computes a choose(n,2) times 2 matrix with all unique pair combinations
    L <- choose(n, 2)
    ret <- matrix(0, L, 2)

    pos <- 1
    for(i in 1 : (n - 1)){
      posnew  <-  pos + n - i - 1
      ret[pos:posnew,1]  <-  i
      ret[pos:posnew,2]  <-  (i + 1) : n
      pos  <-  posnew +1
    }
    return(ret)
  }

  #-----------------------------------------------------------------------


  FQD  <-  function(L, Rvec, GaInv, weight_matrix = "adjusted", d_FQD, m_FQD){

    # computes the quadratic discrepancy function (FQD) with penalty
    #
    #

    # arguments:
    #
    # L:			loading matrix
    # Rvec:			estimated copula correlations
    # GaInv:		weight matrix (information matrix)
    # weight_matrix:	method to compute weight matrix
    # d_FQD:		dimension
    # m_FQD:		number of factors



    Lmat  <-  matrix(L, d_FQD, m_FQD)

    # the factor model assumes:	R = L %*% t(L) + diag(V)

    Rtheta  <-  Lmat %*% t(Lmat)
    RL  <-  Rtheta[lower.tri(Rtheta)]

    tLVL  <-  t(Lmat) %*% diag(1/ diag( 1 - Rtheta ), d_FQD) %*% Lmat


    # penalty term due to the constraint: L^T * V^2 * L = diag
    # diagcheck is large if constraint is not fullfilled

    diagcheck  <-  sqrt(sum((tLVL[lower.tri(tLVL)])^2)) / 10

    # penalty term due to the constraint: Rtheta has to be correlation matrix
    # diagonecheck is large if constraint is not fullfilled

    diagonecheck  <-  sum(pmax(Rtheta, 1) - 1)^2 * 1000

    dL <- d_FQD * (d_FQD - 1)/2
    Delta <- matrix(0, dL, d_FQD * m_FQD)

    # computes the Jacobian matrix partial vecp( R(theta) )/partial theta
    # for the factor model R=L*L^T + diag(V^2)

    for(j in 0:(m_FQD - 1)){
      cc <- 1
      for(i in 1:(d_FQD - 2)){
        diag(Delta[cc : (d_FQD - i + cc - 1), j * d_FQD + ((i + 1) : d_FQD)] )  <-  L[j * d_FQD + i]
        Delta[cc : (d_FQD - i + cc - 1), j * d_FQD + i] <- L[j * d_FQD + ((i + 1) : d_FQD)]
        cc <- cc + d_FQD - i
      }
      Delta[cc, j * d_FQD + d_FQD] <- L[j * d_FQD + d_FQD - 1]
      Delta[cc, j * d_FQD + d_FQD-1] <- L[j * d_FQD + d_FQD]
    }



    Delta_GInv_Delta =  as.matrix( nearPD( t(Delta) %*% GaInv %*% Delta)$mat )


    chol_Delta_GInv_Delta = chol(Delta_GInv_Delta)

    Inv.Delta_GInv_Delta = chol2inv(chol_Delta_GInv_Delta)

    # weight matrix

    correction.term = GaInv%*%Delta%*% Inv.Delta_GInv_Delta %*%t(Delta)%*%GaInv
    correction.term = (correction.term + t(correction.term) )/2

    U = GaInv - correction.term
    U = (U + t(U))/2
    U = as.matrix(nearPD( U )$mat)



    if(weight_matrix=='adjusted'){ret = t(Rvec - RL)%*% U %*%(Rvec-RL)+ d_FQD^3*(diagcheck + diagonecheck)}

    if(weight_matrix=='regular'){ret = t(Rvec - RL)%*% GaInv %*%(Rvec-RL)+ d_FQD^3*(diagcheck + diagonecheck)}

    return( ret )
  }



  FQDnonPenal =  function(L, Rvec, GaInv, weight_matrix = "adjusted", d_FQD, m_FQD){

    # computes the quadratic discrepancy function (FQD) without penalty

    Lmat  =  matrix(L,d_FQD,m_FQD)

    Rtheta  =  Lmat %*% t(Lmat)
    RL  =  Rtheta[lower.tri(Rtheta)]
    dL = d_FQD * ( d_FQD - 1 ) /2

    Delta = matrix(0,dL,d_FQD* m_FQD)

    # computes the Jacobian matrix partial vecp( R(theta) )/partial theta  for the factor model R=L*L^T + diag(V^2)

    for(j in 0:(m_FQD-1)){
      cc = 1
      for(i in 1:(d_FQD -2)){
        diag( Delta[cc:(d_FQD -i+cc-1),j*d_FQD+ ((i+1):d_FQD )] )  =  L[j*d_FQD + i]
        Delta[cc:(d_FQD -i+cc-1),j*d_FQD + i] = L[j*d_FQD+ ((i+1):d_FQD ) ]
        cc = cc+d_FQD-i
      }
      Delta[cc,j*d_FQD + d_FQD]   =  L[j*d_FQD+ d_FQD - 1]
      Delta[cc,j*d_FQD+ d_FQD - 1] =  L[j*d_FQD + d_FQD ]
    }



    Delta_GInv_Delta <- as.matrix(nearPD(t(Delta) %*% GaInv %*% Delta)$mat)

    chol_Delta_GInv_Delta = chol(Delta_GInv_Delta)

    Inv.Delta_GInv_Delta = chol2inv(chol_Delta_GInv_Delta)


    # weight matrix

    correction.term = GaInv%*%Delta%*% Inv.Delta_GInv_Delta %*%t(Delta)%*%GaInv
    correction.term = (correction.term + t(correction.term) )/2

    U = GaInv - correction.term
    U = (U + t(U))/2
    U = as.matrix(nearPD( U )$mat)

    if(weight_matrix=='adjusted'){ret = t(Rvec - RL) %*% U %*% (Rvec-RL)}

    if(weight_matrix=='regular'){ret = t(Rvec - RL)%*% GaInv %*%(Rvec-RL)}


    return( ret )

  }
  #-----------------------------------------------------------------------

  ret = list()
  ret$dof  =  d*( d-1 )/2 - ( d*m - m*( m-1 )/2 )

  if(ret$dof<1){
    print('Error, too many factors - negative dof!')
  }

  if(ret$dof>=1){
    dimIndex = subsets(d)
    chol.Gamma = chol( Ga )

    GI = chol2inv( chol.Gamma )

    CR = matrix(0,d,d)
    CR[dimIndex] = Rvec
    CR = CR+t(CR) + diag(1,d)

    cl = list()
    cl$center = rep(0,d)
    cl$n.obs = n
    cl$cov = CR
    fa = factanal(factors=m,covmat=cl,method='mle',nstart = 5, rotation='none', control=list(opt=list(maxit=100000)))

    Lstart = as.vector(fa$loadings[1:(d*m)])

    if(m>1){ Lstart = as.vector( varimax( matrix(fa$loadings,d,m) )$loadings[1:(d*m)]) }

    tmp = optim(Lstart,FQD, Rvec=Rvec,GaInv=GI,d_FQD = d,m_FQD = m, weight_matrix = weight_matrix_adf, control = list( maxit = maxit_optim ) )




    rot  =  matrix(tmp$par,d,m)

    if(m>1){ rot  =  matrix(varimax(matrix(tmp$par,d,m))$loadings[1:(d*m)],d,m) }


    ret$nF  =  n*FQDnonPenal( rot,R=Rvec,GaInv=GI, weight_matrix = weight_matrix_adf, d_FQD=d,m_FQD=m)

    ret$nFP  =  n*FQD(rot,R=Rvec,GaInv=GI, weight_matrix = weight_matrix_adf, d_FQD=d,m_FQD=m)

    ret$conv  =  tmp$convergence

    ret$load  =  rot

    ret$uniq  =  1 - diag(matrix(tmp$par,d,m)%*%t(matrix(tmp$par,d,m)))

    ret$pValue  =  pchisq(ret$nF,ret$dof,lower.tail=FALSE)

    ret$wm = weight_matrix_adf

    return( ret )
  }
}
#------------------------------------------------------------------------------------
