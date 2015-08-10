#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix asymp_var_tail(NumericVector RANKM, int n, int k, int d, NumericVector dimInd,  NumericVector x, NumericVector y,
             NumericVector dlaX, NumericVector dlaY, NumericVector dlaR, NumericVector dlaNu, NumericVector Weight,
             NumericVector a, NumericVector b){
  int dL = d * (d - 1) / 2;
  double sumweight = 0;
  double EBtilde = 0;
  double sumSigmaijkl = 0;
  double sumSigma23 = 0;
  double sumSigma1 = 0;

  double weight = 0;

  double ki = 1 / k;

  int aij = 0;
  int akl = 0;
  int bij = 0;
  int bkl = 0;
  int L1 = 0;
  int L2 = 0;
  int L3 = 0;
  int L4 = 0;

  NumericMatrix Sigma(dL * dL, 3);

  for(int ij = 0; ij < dL; ++ij){
    aij = a[ij] - 1;
    bij = b[ij];
    for(int kl = ij; kl < dL; ++kl){
      akl = a[kl] - 1;
      bkl = b[kl];

      for(int p1 = aij; p1 < bij; ++p1){
        for(int p2 = akl; p2 < bkl; ++p2){
          for(int i = 0; i < n; ++i){
              if(RANKM[i + n * dimInd[ij]] >  (n - k * x[p1])){
                L1 = 1;
              } else{
                L1 = 0;
              }
              if(RANKM[i + n * dimInd[ij+dL]] > (n -  k * y[p1])){
                L2 = 1;
              } else{
                L2 = 0;
              }
              if(RANKM[i + n * dimInd[kl]]> (n - k * x[p2])){
                L3 = 1;
              } else{
                L3 = 0;
              }
              if(RANKM[i +  n * dimInd[kl+dL]] > (n - k * y[p2])){
                L4 = 1;
              } else{
                L4 = 0;
              }
              EBtilde += L1 * L2 * L3 * L4 * ki -
                L1 * L2 * L3 * dlaX[kl + p2 * dL] * ki -
                L1 * L2 * L4 * dlaY[kl + p2 * dL] * ki -
                L1 * L3 * L4 * dlaX[ij + p1 * dL] * ki -
                L2 * L3 * L4 * dlaY[ij + p1 * dL] * ki +
                L1 * L3 * dlaX[ij + p1 * dL] * dlaX[kl + p2 * dL] * ki +
                L1 * L4 * dlaX[ij + p1 * dL] * dlaY[kl + p2 * dL] * ki +
                L2 * L3 * dlaY[ij + p1 * dL] * dlaX[kl + p2 * dL] * ki +
                L2 * L4 * dlaY[ij + p1 * dL] * dlaY[kl + p2 * dL] * ki;
          }

          weight  = Weight[p1] * Weight[p2];
          sumweight += weight;
          EBtilde *= weight;

          sumSigmaijkl += (EBtilde / dlaR[ij + p1 * dL] / dlaR[kl + p2 * dL]);
          sumSigma23 += (EBtilde / dlaR[ij + p1 * dL] / dlaNu[kl + p2 * dL]);
          sumSigma1 += (EBtilde / dlaNu[ij + p1 * dL]/dlaNu[kl + p2 * dL]);
        }
      }
      Sigma(ij + kl*dL, 3) = sumSigmaijkl / sumweight;
      Sigma(ij + kl*dL, 2) =  sumSigma23 / sumweight;
      Sigma(ij + kl*dL,1) = sumSigma1 / sumweight;
    }
  }
  return Sigma;
}
