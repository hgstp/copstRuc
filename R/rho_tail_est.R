# param:
#
#  Xmult:          	dataset containing n observations of the d-dimensional random variables.
#  k:          		number of the largest order statistics, which is used for estimating the tail dependence function. Default = 1.
#  prec:          	number of samples, which is used for smoothing the results of the estimation. Default = 31
#  asymp:          	if asymp=TRUE, the asymptotic covariance matrix of lower-diagonal elements of the copula correlation estimator is calculated. Default = FALSE
#  status:          	if status=TRUE, the progress of calculation is printed.
#

 rho_tail_est <- function(Xmult, k, prec = 31, asymp = FALSE){

   n <- dim(Xmult)[1]
   d <- dim(Xmult)[2]

   # number of pairs
   dL <- d * (d - 1) / 2
   dInd <- subsets(d)

   RANKmult <- apply(Xmult, 2, rank)

   # correlation estimate based on kendall's tau
   Rtau <- sin(tauest(Xmult)$tau * pi / 2)

   thetavec <- theta_function(seq(0.05, pi / 2 - .05, length = prec), 2.5)
   thetavec[round(prec/2, 0)] <- pi / 4

   x <- sqrt(2) * cos(thetavec)
   y <- sqrt(2) * sin(thetavec)

   laE <- matrix(0, dL, prec)

   out <- .C('laEst', RANKM = as.integer(RANKmult), n = as.integer(n),
             k = as.integer(k), d = as.integer(d), prec = as.integer(prec),
             x = as.double(x), y = as.double(y), laE = as.double(laE))

   laE <- matrix(out$laE, dL, prec)

   # matrix of la-estimator: la[ij,phi]

   weight <- 1 - (thetavec * 4 / pi - 1)^2


   # estimated nu for x=1, y=1
   nu_at_tdc = 1:dL

   # index in thetavec corresponding to x=1, y=1
   piO4 <- as.integer(round(prec / 2, 0))


   for(i in 1:dL){
     nu_at_tdc[i] <- nu_solve(laE[i, piO4], 1, 1, Rtau[dInd[i, 1], dInd[i, 2]])
   }

   Qij <- matrix(0, dL, 2)
   Rvec <- Rtau[lower.tri(Rtau)]
   expo <- (1 - k^(-1 / 4)) * nu_at_tdc * abs(log(pmax(Rvec, 1e-10)))


   # range of admissible points for hat.Q.ij.star

   Qij[, 1] <- pmax(.1, atan(exp(-expo)))
   Qij[, 2] <- pmin(pi / 2 - .1, atan(exp(expo)))

   nu_tilde_ij <- matrix(0, dL, prec)
   nu_hat_ij  =  (1 : dL) * NA

   for(i in 1 : dL){
     if(nu_at_tdc[i]){
       a <- min(which(thetavec >= Qij[i,1]))
       b <- max(which(thetavec <= Qij[i,2]))
       ind <- (1 : prec) * 0

       # computes nu_tilde_ij if thetavec is element of the set hat.Q.ij.star
       for(p in a : b){
         if((laE[i,p] < ellip_tdf(x[p], y[p], abs(log(x[p] / y[p]) / log(max(Rvec[i], 0))), Rvec[i]))
            & (abs(log(x[p] / y[p]) / log(max(Rvec[i], 0))) < nu_at_tdc[i] * (1 - k^(-1 / 4)))){
           nu_tilde_ij[i,p] <- nu_solve(laE[i, p], x[p], y[p], Rvec[i],
                                        am = abs(log(x[p] / y[p]) / log(max(Rvec[i], 0))))
           ind[p] <- 1
         }
         }
       nu_hat_ij[i] <- sum(weight[a:b] * nu_tilde_ij[i, a : b]) / sum(weight * ind)
     }
     }
   # result of function rhoTailEst
   ret <- list()
   ret$nu <- mean(nu_hat_ij, na.rm = TRUE)
   nu <- ret$nu

   #now, nu is estimated; next we estimate rho for x=1, y=1

   rho_at_tdc <- 1 : dL

   for(i in 1 : dL){
     rho_at_tdc[i] <- rhosolve(laE[i, piO4], 1, 1, nu)
     }
   # range of admissible points for hat.U.ij.star

   Uij <- matrix(0, dL, 2)
   Uij[, 1] <- atan((pmax(rho_at_tdc, 0) * (1 - k^(-1 / 4)))^nu)
   Uij[, 2] <- atan((pmax(rho_at_tdc, 0) * (1 - k^(-1 / 4)))^(-nu))

   rho_tilde_ij <- matrix(0, dL, prec)
   rho_hat_ij <- (1 : dL) * NA

   for(i in 1 : dL){
     if(rho_at_tdc[i]){
       a <- min(which(thetavec >= Uij[i,1]))
       b <- max(which(thetavec <= Uij[i,2]))
       ind <- (1 : prec) * 0
       for(p in a : b){
         if((laE[i,p] < ellip_tdf(x[p], y[p], nu,
                                  min(exp(-abs(log(tan(thetavec[p]))) / nu), 1 - 1e-7)))
            & (rho_at_tdc[i] * (1 - k^(-1 / 4)) < exp(-abs(log(tan(thetavec[p]))) / nu))){
           rho_tilde_ij[i, p] <- rhosolve(laE[i, p], x[p], y[p], nu,
                                          rM = min(exp(-abs(log(tan(thetavec[p]))) / nu), 1 - 1e-7))
           ind[p] <- 1
         }
         }
       rho_hat_ij[i] <- sum(weight[a : b] * rho_tilde_ij[i, a : b]) / sum(weight * ind)
     }
   }

   # estimated copula correlation matrix
   Rho_lambda <- matrix(0, d, d)

   Rho_lambda[lower.tri(Rho_lambda)] <- rho_hat_ij
   Rho_lambda <- Rho_lambda + t(Rho_lambda)
   diag(Rho_lambda) <- 1

   # ret$R = PosDefCorr(Rho_lambda)

   ret$R <- nearPD(Rho_lambda, corr = TRUE )$mat
   Rvec <- ret$R[lower.tri(ret$R)]
   ret$Rvec <- Rvec
   #now, Rho_lambda is estimated

   if(asymp){
     a <- 1 : dL
     b <- 1 : dL
     for(i in 1 : dL){
       a[i] <- max(min(which(thetavec >= Uij[i, 1])), min(which(thetavec >= Qij[i, 1])))
       b[i] <- min(max(which(thetavec <= Uij[i, 2])), max(which(thetavec <= Qij[i, 2])))
       }
    #x <- sqrt(2)*cos(thetavec)
    #y = sqrt(2)*sin(thetavec)

    dlaX <- matrix(0, dL, prec)
    dlaY <- matrix(0, dL, prec)
    dlaR <- matrix(0, dL, prec)
    dlaNu <- matrix(0, dL, prec)
    dlaInNu <- (1 : dL) * 0

    # partial derivatives of tail dependence function for estimated nu and R

    for(ij in 1:dL){
      dlaX[ij, ] <- derivative_lambda_x(x, y, nu, Rvec[ij])
      dlaY[ij, ] <- derivative_lambda_y(x, y, nu, Rvec[ij])
      dlaR[ij, ] <- derivative_laRho(x, y, nu, Rvec[ij])
      dlaNu[ij, ] <- derivative_laNu(x, y, nu, Rvec[ij])
    }

    # part of Sigma_1 equation (3.12)

    for(ij in 1 : dL){
      for(p in a[ij] : b[ij]){
        tmp <- derivative_laInverse_nu(ellip_tdf(x[p], y[p], nu, Rvec[ij]),
                                       x[p], y[p], nu) * weight[p]
        dlaInNu[ij] <- dlaInNu[ij] + tmp
      }
      dlaInNu[ij] <- dlaInNu[ij] / sum(weight[a[ij] : b[ij]])
    }

    dimInd <- subsets(d)

    Sigma_1 <- matrix(0, dL, dL)
    Sigma_2_3 <- matrix(0, dL, dL)
    Sigma_4 <- matrix(0, dL, dL)

    out <- .C("La4EstN", RANKM = as.integer(RANKmult), n = as.integer(n), k = as.integer(k),
              d = as.integer(d), dimInd = as.integer(dimInd - 1),
              x = as.double(x), y = as.double(y), dlaX = as.double(dlaX), dlaY = as.double(dlaY),
              dlaR = as.double(dlaR), dlaNu = as.double(dlaNu), Weight = as.double(weight),
              a = as.double(a), b = as.double(b),
              Sigma4 = as.double(Sigma_4), Sigma23 = as.double(Sigma_2_3), Sigma1 = as.double(Sigma_1))

    Sigma_4 <- matrix(out$Sigma4, dL, dL)
    Sigma_4 <- Sigma_4 + t(Sigma_4)
    diag(Sigma_4) <- diag(Sigma_4) / 2

    Sigma_2_3 <- matrix(out$Sigma23, dL, dL)
    Sigma_2_3 <- Sigma_2_3 + t(Sigma_2_3)
    diag(Sigma_2_3) <- diag(Sigma_2_3) / 2
    Sigma_2_3 <- apply(Sigma_2_3, 1, sum)

    Sigma_1 <- matrix(out$Sigma1, dL, dL)
    Sigma_1 <- 2 * sum(Sigma_1) - sum(diag(Sigma_1))

    ret$Ga <- matrix(0, dL, dL)

    # equation (3.11)
    ret$Ga <- (dlaInNu %*% t(dlaInNu)) * Sigma_1 * 4 / d^2 / (d - 1)^2 +
      (dlaInNu %*% t(Sigma_2_3) + Sigma_2_3 %*% t(dlaInNu)) * 2 / d / (d - 1) + Sigma_4

    ev <- eigen(ret$Ga)
    evMin <-  min(ev$values)
    old <- ret$Ga
    if(min(ev$values) < 2e-7){
      ret$Ga <- as.matrix(nearPD(ret$Ga)$mat)
      print(paste('not positiv definit matrix projected. rel.distance=',signif(sum((ret$Ga-old)^2)/sum(old^2),2),' smallest EV=',evMin,sep=''))
      }
    }
   return(ret)
 }

