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
double akkadianRootF(double x, double a)
{
  return (x+a/x)/2;
}


//' The inferred method of Akkadian-Babylonian mathematicians to approximate
//' square roots
//'
//' @param a the number to compute square root of
//' @param n the number of iterations to take in the method
//'
//' @description {For finding the square-root of whole numbers.}
//' @details {The algorithm is on wikipedia and quite straightforward.}
//' @export akkadianRoot
// [[Rcpp::export]]
double akkadianRoot(double a, int n = 100)
{
  NumericVector r(n);
  r[0] = a;
  for(int i = 1; i < n; i++)
  {
    r[i] = akkadianRootF(r[i-1], a);
  }
  return r[n-1];
}
