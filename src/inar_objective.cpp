#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double inar_objective(double gamma, NumericVector y, NumericVector mu) {
  int n = y.size();
  double s = 0;
  for (int i = 1; i < n; ++i) {
    double a = y[i]*y[i-1] - mu[i]*mu[i-1]*(1.0 + 1.0/(mu[i-1] + gamma));
    s += a*a;
  }
  return s;
}
