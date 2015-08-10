#include <Rcpp.h>
using namespace Rcpp;

// The function computes an estimate of the asymptotic covariance matrix
// of the copula correlation estimator based on Kendall's tau - see (4.1) in
// Klueppelber, Kuhn (2009). The entries of the asymptotic covariance matrix
// Gamma are given in (4.2)
//
//

// [[Rcpp::export]]
  NumericVector asymp_var_tau(NumericVector data, int d, int n){
    int dL = d * (d - 1) / 2;
    double N = 1 / (n - 1);

    NumericVector var_tau(dL*dL);

  for(int dim1b = 0; dim1b < d - 1; ++dim1b){

    for(int dim2b = dim1b + 1; dim2b < d; ++dim2b){

      for(int dim1 = 0; dim1 < d - 1; ++dim1){

        for(int dim2 = dim1 + 1; dim2 < d; ++dim2){

          int ij = dim1 * d -dim1 * (dim1 + 1) / 2 + dim2 - dim1 - 1;
          int kl = dim1b * d - dim1b * (dim1b + 1) / 2 + dim2b - dim1b - 1;
          var_tau[ij + kl*dL] = 0;

          if(ij >= kl){

            for(int i = 0; i < n; ++i){
              double tSij = N;
              double tSkl = N;

              for(int j = 0; j < n; ++j){
                if ((data[i + dim1 * n] - data[j + dim1 * n]) * (data[i + dim2 * n] - data[j + dim2 * n]) > 0) {
                  tSij += N;
                } else {
                  tSij += -N;
                }
                if ((data[i + dim1b * n] - data[j + dim1b * n]) * (data[i + dim2b * n] - data[j + dim2b * n]) > 0) {
                  tSkl += N;
                } else {
                  tSkl += -N;
                }
              }
              var_tau[ij + kl * dL] += (tSij * tSkl) / n;
            }
          }
        }
      }
    }
  }
  return var_tau;
}


