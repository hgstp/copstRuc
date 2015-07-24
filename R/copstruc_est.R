# arguments:
#
#  x:          	dataset containing n observations of the d-dimensional random variables.
#  method:         	name of the estimator, which can be 'tail', 'kendall' or 'standard'.
#  k:          		number of the largest order statistics, which is used for estimating the tail dependence function. Default = 1.
#  prec:          	number of samples, which is used for smoothing the results of the estimation. Default = 31
#  asymp:          	if asymp=TRUE, the asymptotic covariance matrix of lower-diagonal elements of the copula correlation estimator is calculated. Default = FALSE
#  status:          	if status=TRUE, the progress of calculation is printed.
#
# value:
#
#    R:          	the estimated correlation matrix
#    Rvec:          	lower diagonal elements of R
#    Ga:          	asymptotic Covariance of Rvec (if asymp==TRUE)
#    nu:          	the estimated tail index (if method=='tail')
#    method:          	the method being used





copstruc_est <- function(x, method = "kendall", k = 10, prec = 31,
                         asymp = FALSE, k_select = FALSE){

   # now the definition of copstruc_est starts


  if(!(method == 'kendall') & !(method == 'tail') & !(method == 'standard')){
    print('Admissible methods are kendall, tail or standard!')
  }

  if(method == 'kendall'){

    te <- rho_tau_est(x, asymp)
    ret <- list()

    # ret$R   =  PosDefCorr(sin(te$tau * pi/2))
    ret$R <- nearPD(sin(te$tau * pi/2), corr = TRUE)$mat
    ret$Rvec <- ret$R[lower.tri(ret$R)]

    if(asymp){
      ret$Ga <- pi^2 * (cos(te$tauVec * pi/2) %*% cos(t(te$tauVec) * pi/2)) * (te$tau4 - (te$tauVec%*%t(te$tauVec)))
      ev <- eigen(ret$Ga)
      old <- ret$Ga

      if(min(ev$values) < 2e-7){
        ret$Ga <- as.matrix(nearPD(ret$Ga)$mat)
        print(paste('not positiv definit matrix projected. rel.distance=',
                    signif(sum((ret$Ga - old)^2) / sum(old^2), 2), sep = ''))
      }
    }

    ret$method <- method
    ret$asymp <- asymp
    ret$obs <- dim(x)[1]
    return(ret)
  }

  if(method == 'tail' & k_select == TRUE){
    k_selection( x )
  }


  if(method == 'tail' & !(k_select == TRUE)){

    ret <- rho_tail_est(x, k, prec, asymp)
    ret$method <- method
    ret$asymp <- asymp
    ret$obs <- k
    return(ret)
  }

  if(method == 'standard'){

    ret <- list()

    #ret$R = PosDefCorr(cor(x))
    ret$R <- nearPD(cor(x), corr = TRUE)$mat
    ret$Rvec <- ret$R[lower.tri(ret$R)]

    if(asymp){
      n <- dim(x)[1]
      d <- dim(x)[2]

      dL <- d * (d - 1) / 2
      index <- subsets(n)
      dI <- subsets(d)
      Vest <- ret$Rvec
      mu <- apply(x,2,mean)

      S <- var(x)
      R4 <- matrix(0,dL,dL)
      ret$Ga <- matrix(0,dL,dL)

      for(ij in 1 : (dL - 1)){
        for(kl in (ij + 1) : dL){
          R4[ij,kl] <- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]]) * (x[, dI[ij, 2]] - mu[dI[ij, 2]]) *
                             (x[, dI[kl, 1]] - mu[dI[kl, 1]]) *
                             (x[, dI[kl, 2]] - mu[dI[kl,2]])) / sqrt(S[dI[ij, 1], dI[ij,1]] *
                                                                       S[dI[ij, 2], dI[ij, 2]] *
                                                                       S[dI[kl, 1], dI[kl, 1]] *
                                                                       S[dI[kl, 2], dI[kl, 2]])
        }
      }

      R4 <- (R4 + t(R4)) / n

      for(ij in 1 : dL){
        R4[ij, ij] <- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]])^2 *
                            (x[, dI[ij, 2]] - mu[dI[ij, 2]])^2) / n / (S[dI[ij, 1], dI[ij, 1]] *
                                                                         S[dI[ij, 2], dI[ij, 2]])
      }

      R <- cor(x)

      for(ij in 1 : dL){
        for(kl in ij : dL){
          Riikl <- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]])^2 * (x[, dI[kl, 1]] - mu[dI[kl, 1]])*
                         (x[, dI[kl, 2]] - mu[dI[kl, 2]])) / sqrt(S[dI[ij, 1], dI[ij, 1]]^2 *
                                                                    S[dI[kl, 1], dI[kl, 1]] *
                                                                    S[dI[kl, 2], dI[kl, 2]] ) / n

          Rjjkl <- sum((x[, dI[ij, 2]] - mu[dI[ij, 2]])^2 * (x[, dI[kl, 1]] - mu[dI[kl, 1]]) *
                         (x[, dI[kl, 2]] - mu[dI[kl, 2]])) / sqrt(S[dI[ij, 2], dI[ij, 2]]^2 *
                                                                    S[dI[kl, 1], dI[kl, 1]] *
                                                                    S[dI[kl, 2], dI[kl, 2]]) / n

          Rkkij <- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]]) * (x[, dI[ij, 2]] - mu[dI[ij, 2]]) *
                         (x[, dI[kl, 1]] - mu[dI[kl, 1]])^2) / sqrt(S[dI[ij, 1], dI[ij, 1]] *
                                                                      S[dI[ij, 2], dI[ij, 2]] *
                                                                      S[dI[kl, 1], dI[kl, 1]]^2 ) / n

          Rllij <- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]]) * (x[, dI[ij, 2]] - mu[dI[ij, 2]]) *
                         (x[, dI[kl, 2]] - mu[dI[kl, 2]])^2) / sqrt(S[dI[ij, 1], dI[ij, 1]] *
                                                                      S[dI[ij, 2], dI[ij, 2]] *
                                                                      S[dI[kl, 2], dI[kl, 2]]^2) / n

          Riikk<- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]])^2 *
                        (x[, dI[kl, 1]] - mu[dI[kl, 1]])^2) / S[dI[ij, 1], dI[ij, 1]] / S[dI[kl, 1], dI[kl, 1]] / n

          Riill <- sum((x[, dI[ij, 1]] - mu[dI[ij, 1]])^2 *
                         (x[, dI[kl, 2]] - mu[dI[kl, 2]])^2) / S[dI[ij, 1], dI[ij, 1]] / S[dI[kl, 2], dI[kl, 2]] / n

          Rjjkk <- sum((x[, dI[ij, 2]] - mu[dI[ij, 2]])^2 *
                         (x[, dI[kl, 1]] - mu[dI[kl, 1]])^2) / S[dI[ij, 2], dI[ij, 2]] / S[dI[kl, 1], dI[kl, 1]] / n

          Rjjll <- sum((x[, dI[ij, 2]] - mu[dI[ij, 2]])^2 *
                         (x[, dI[kl, 2]] - mu[dI[kl, 2]])^2) / S[dI[ij, 2], dI[ij, 2]] / S[dI[kl, 2], dI[kl, 2]] / n


          ret$Ga[ij,kl] <- R4[ij, kl] - R[dI[ij, 1], dI[ij, 2]] * (Riikl + Rjjkl) / 2 -
            R[dI[kl, 1], dI[kl, 2]] * (Rkkij + Rllij) / 2 + R[dI[ij, 1], dI[ij, 2]] *
            R[dI[kl, 1], dI[kl, 2]] * (Riikk + Riill + Rjjkk + Rjjll) / 4
        }
      }

      ret$Ga <- ret$Ga + t(ret$Ga)
      diag(ret$Ga) <- diag(ret$Ga / 2)
      old <- ret$Ga
      ev <- eigen(ret$Ga)

      if(min(ev$values) < 2e-7){
        ret$Ga <- as.matrix(nearPD( ret$Ga )$mat)
        print(paste('not positiv definit matrix projected. rel.distance=',
                    signif(sum((ret$Ga - old)^2) / sum(old^2), 2), sep = ''))

      }

    }

    ret$asymp <- asymp
    ret$obs <- dim(x)[1]
    return(ret)
  }
}



