#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp.
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @param a the parameter for gibbsc
//' @param b the parameter for gibbsc
//' @param n the parameter for gibbsc
//' @return a random sample from gibbs \code{mat}
//' @examples
//' \dontrun{
//'     gc <- gibbsC(100, 10, 2, 4, 16)
//'     print(gc)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N,  int thin, double a, double b, int n) {
  NumericMatrix mat(N, 2);
  double x = floor(n/2), y = 0.5;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y )[0];
      y = rbeta(1, x+a, n-x+b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
