#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dist_squared(NumericMatrix X, NumericMatrix C) {
  int n = X.nrow();
  int p = X.ncol();
  int k = C.nrow();
  NumericMatrix dist(n,k);
  for (int i = 0; i < n; i++){
    for (int j = 0; j < k; j++){
      double sum = 0;
      for (int l = 0; l < p; l++) {
        sum += pow(X(i, l) - C(j, l), 2);
      }
      dist(i, j) = sum;
    }
  }
  return dist;
}
