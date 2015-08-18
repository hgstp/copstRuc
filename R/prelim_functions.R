# we first define some preliminary functions

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
 #Projection of a not positive (semi)-definit Correlation matrix to set of pos.def. Corr matrices
  PosDefCorr  <-  function(Rgiven){
  d <- dim(Rgiven)[1]
  X <- Rgiven
  ev <- eigen(X)
  deltaS <- X * 0
  Xold <- X * 0
  if(min(ev$values) < 1e-7 ){
    check <- TRUE
    count <- 0

    while(check){
      count <- count + 1
      R <- X - deltaS
      X <- ev$vectors %*% diag(pmax((d : 1) * 2e-7, ev$values), d) %*% t(ev$vectors)
      X <- (X + t(X)) / 2

      deltaS <- X - R
      diag(X) <- 1
      ev <- eigen(X)
      if(((min(ev$values) >= 0) & (sum((Xold - X)^2) < 1e-7)) | (count > 100)){
        check <- FALSE
      }

      Xold <- X
    }
    print(paste('Info: not positive definite matrix projected; steps=',count,', rel. distance=',signif(sum((Rgiven-X)^2)/sum(Rgiven^2),2),sep=''))

    if(count > 100){
      print('Warning: no convergence of projection')
      }
  }
  return(X)
}

#------------------------------------------
  FQD  <-  function(L, R, GaInv, weight_matrix = "adjusted", d_FQD, m_FQD){

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



    Delta_GInv_Delta <- as.matrix(Matrix::nearPD(t(Delta) %*% GaInv %*% Delta)$mat)


    chol_Delta_GInv_Delta <- chol(Delta_GInv_Delta)

    Inv_Delta_GInv_Delta <- chol2inv(chol_Delta_GInv_Delta)

    # weight matrix

    correction_term = GaInv%*%Delta%*% Inv_Delta_GInv_Delta %*%t(Delta)%*%GaInv
    correction_term = (correction_term + t(correction_term) )/2

    U = GaInv - correction_term
    U = (U + t(U))/2
    U = as.matrix(Matrix::nearPD( U )$mat)



    if(weight_matrix == 'adjusted'){
      ret <- t(Rvec - RL) %*% U %*% (Rvec - RL) + d_FQD^3 * (diagcheck + diagonecheck)
    }

    if(weight_matrix == 'regular'){
      ret <- t(Rvec - RL) %*% GaInv %*% (Rvec - RL) + d_FQD^3 * (diagcheck + diagonecheck)
    }

    return(ret)
  }

#---------------------------


  FQDnonPenal =  function(L, R, GaInv, weight_matrix = "adjusted", d_FQD, m_FQD){

    # computes the quadratic discrepancy function (FQD) without penalty

    Lmat <- matrix(L, d_FQD, m_FQD)

    Rtheta  <-  Lmat %*% t(Lmat)
    RL  <-  Rtheta[lower.tri(Rtheta)]
    dL <- d_FQD * (d_FQD - 1) / 2

    Delta <- matrix(0, dL, d_FQD * m_FQD)

    # computes the Jacobian matrix partial vecp( R(theta) )/partial theta
    # for the factor model R=L*L^T + diag(V^2)

    for(j in 0 : (m_FQD - 1)){
      cc <- 1
      for(i in 1 : (d_FQD - 2)){
        diag(Delta[cc : (d_FQD - i + cc - 1), j * d_FQD + ((i + 1) : d_FQD )]) <- L[j * d_FQD + i]
        Delta[cc : (d_FQD - i + cc - 1), j * d_FQD + i] <- L[j * d_FQD + ((i + 1) : d_FQD)]
        cc <- cc + d_FQD - i
      }
      Delta[cc, j * d_FQD + d_FQD] <- L[j * d_FQD + d_FQD - 1]
      Delta[cc, j * d_FQD + d_FQD - 1] <- L[j * d_FQD + d_FQD ]
    }


    Delta_GInv_Delta <- as.matrix(Matrix::nearPD(t(Delta) %*% GaInv %*% Delta)$mat)

    chol_Delta_GInv_Delta <- chol(Delta_GInv_Delta)

    Inv_Delta_GInv_Delta <- chol2inv(chol_Delta_GInv_Delta)


    # weight matrix

    correction_term <- GaInv %*% Delta %*% Inv_Delta_GInv_Delta %*% t(Delta) %*% GaInv
    correction_term <- (correction_term + t(correction_term)) / 2

    U <- GaInv - correction_term
    U <- (U + t(U)) / 2
    U <- as.matrix(Matrix::nearPD(U)$mat)

    if(weight_matrix == 'adjusted'){
      ret <- t(Rvec - RL) %*% U %*% (Rvec - RL)
    }

    if(weight_matrix == 'regular'){
      ret <- t(Rvec - RL) %*% GaInv %*% (Rvec - RL)
    }
    return(ret)
  }
#----------------------------------------
 # some function to compute partial derivatives of the tail dependence function
 derivative_laInverse_nu  =  function(LAE,x,y,nu,eps=1e-7){
   ret <- (rhosolve(LAE, x, y, nu * (1 + eps)) - rhosolve(LAE, x, y, nu * (1 - eps))) / eps / 2 / nu
 }


 derivative_lambda_x <- function(x, y, nu, rho){
   x <- pmax(x, 1e-50)
   y <- pmax(y, 1e-50)
   (1 - pt(sqrt((nu + 1) / (1 - rho^2)) * ((x / y)^(1 / nu) - rho), nu + 1)) -
     dt(sqrt((nu + 1) / (1 - rho^2)) * ((x / y)^(1 / nu) - rho), nu + 1) * sqrt((nu + 1) / (1 - rho^2)) * (x / y)^(1 / nu) / nu +
     dt(sqrt((nu + 1) / (1 - rho^2)) * ((y / x)^(1 / nu) - rho), nu + 1) * sqrt((nu + 1) / (1 - rho^2)) * (y / x)^(1 / nu + 1) / nu
 }

 derivative_lambda_y <- function(x, y, nu, rho){
   x <- pmax(x, 1e-50)
   y <- pmax(y, 1e-50)
   (1 - pt(sqrt((nu + 1) / (1 - rho^2)) * ((y / x)^(1/nu) - rho), nu + 1)) -
     dt(sqrt((nu + 1) / (1 - rho^2)) * ((y / x)^(1 / nu) - rho), nu + 1) * sqrt((nu + 1) / (1 - rho^2)) * (y / x)^(1 / nu) / nu +
     dt(sqrt((nu + 1) / (1 - rho^2)) * ((x / y)^(1 / nu) - rho), nu + 1) * sqrt((nu + 1) / (1 - rho^2)) * (x / y)^(1 / nu + 1) / nu
 }

 derivative_laRho <- function(x, y, nu, rho, eps = 1e-7){
   (ellip_tdf(x, y, nu, rho * (1 + eps)) - ellip_tdf(x, y, nu, rho * (1 - eps))) / eps / 2 / rho
 }

 derivative_laNu <- function(x, y, nu, rho, eps = 1e-07){
   (ellip_tdf(x, y, nu * (1 + eps), rho) - ellip_tdf(x, y, nu * (1 - eps), rho)) / eps / 2 / nu
 }

 #-----------------------------------------------------------------------
 # now we define functions to compute the inverse of the tail dependence function with respect to rho and nu

 lambda_rho = function(rho, X, Y, Nu, La){
   # computes difference of parametric and empirical tail copula estimate
   # rho has to be the first argument due to uniroot()

   X * (1 - pt(sqrt((Nu + 1) / (1 - rho^2)) * ((X / Y)^(1 / Nu ) - rho), Nu + 1)) +
     Y*(1 - pt(sqrt((Nu + 1) / (1 - rho^2)) * ((Y / X)^(1 / Nu) - rho), Nu + 1)) - La
 }

#--------------------------------------------

 rhosolve <- function(LAE, x, y, nu, rM = 1-1e-7){

   check <- lambda_rho(rM, x, y, nu, LAE)
   if(check < 0){
     ret <- NA
     }
   if(check >= 0){
     ret <- uniroot(lambda_rho, c(-1 + 1e-5, rM),
                    X = x, Y = y, Nu = nu, La = LAE)$root
     }
   return(ret)
 }

 #---------------------------------------------------
 # computes difference of parametric and empirical tail copula estimate
 # nu has to be the first argument due to uniroot()

 lambda_nu <- function(nu, X, Y, Rho, La){
   X * (1 - pt(sqrt((nu + 1) / (1 - Rho^2)) * ((X / Y)^(1 / nu) - Rho), nu + 1)) +
     Y*(1 - pt(sqrt((nu + 1) / (1 - Rho^2)) * ((Y / X)^(1 / nu) - Rho), nu + 1)) - La
 }

#-----------------------------------------------------

 nu_solve <- function(LAE, x, y, rho, am = 0){

   true <- 1
   # upper endpoint
   aM <- 1

   while(true){
     if(lambda_nu(aM, x, y, rho, LAE) > 0){
       aM <- aM * 2
       }
     if(lambda_nu(aM, x, y, rho, LAE) <= 0){
       true <- 0
       }
   }
   true <- 1
   count <- 0

   if(am == 0){
     am <- 1
     while(true){
       if(lambda_nu(am, x, y, rho, LAE) < 0){
         am <- am / 2
         count <- count + 1
       }
       if((lambda_nu(am, x, y, rho, LAE) >= 0) | (count > 12)){
         true <- 0
         }
     }
   }

   if(lambda_nu(am, x, y, rho, LAE) <= 0){
     ret <- 0
     }

   if(lambda_nu(am, x, y, rho, LAE) > 0){
     ret <- uniroot(lambda_nu, c(am, aM), X = x, Y = y, Rho = rho, La = LAE)$root
     }
   return(ret)
 }

 #-------------------------------------------------

 theta_function <- function(x, m){
   d <- (1 - 1 / m) / pi^2 * 16
   d * x^3 - 3 / 4 * pi * d * x^2 + (3 / 16 * pi^2 * d + 1 / m) * x
 }



