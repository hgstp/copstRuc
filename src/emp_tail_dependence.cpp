#include <Rcpp.h>
using namespace Rcpp;

// computes the empirical tail dependence function at prec values
// of the vector (x(theta),y(theta)) for all d *(d - 1) / 2 pairs
//
//

// [[Rcpp::export]]
NumericVector emp_tail_dependence(NumericMatrix RANKM, double n, double k, double d, int prec, NumericVector x, NumericVector y){
  int ij = 0;
  double dL = d *(d - 1) / 2;
  double tmp = 0;

  NumericVector laE(dL * prec);



  for(int p = 0;  p< prec; ++p){
    for(int dim1 = 0; dim1 < d - 1; ++dim1){
      for(int dim2 = dim1 + 1; dim2 < d; ++dim2){
        double count = 0;

        for(int nc = 0; nc < n; ++nc){
          if((RANKM(nc, dim1) > n - k * x[p]) & (RANKM(nc, dim2) > n - k * y[p])){
            count += 1;
          } else {
              count += 0;
            }
        }

        ij = dim1 * d - dim1 * (dim1 + 1) / 2 + dim2 - dim1 - 1;

        if(count == 0){
         tmp = 1 / (2 * k);
        }

        laE[p * dL + ij] =  count / k + tmp;
      }
    }
  }
  return laE;

}

