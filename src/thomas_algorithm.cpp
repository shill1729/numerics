#include <Rcpp.h>
using namespace Rcpp;


//' Thomas algorithm for solving tridiagonal linear systems. Stability only for diagonally dominant or symmetric positive definite.
//'
//' @param a lower diagonal of matrix
//' @param b diagonal of matrix
//' @param c upper diagonal of matrix
//' @param d vector on RHS of equation
//'
//' @description {Thomas algorithm for tridiagonal linear systems. This can be used to solve
//' PDEs arising from Feynman-Kac connections.}
//' @details {The algorithm is on wikipedia and quite straightforward.}
//' @export thomas_algorithm
// [[Rcpp::export]]
NumericVector thomas_algorithm(NumericVector a, NumericVector b, NumericVector c, NumericVector d)
{
  // On paper, b and d are from 1 to n, the vector a, the lower diag is from 2 to n, and c, the upper diag, is from 1 to n-1.
  int n = d.size();
  // To preserve the coefficients, we make new ones for the forward sweep and x for the solution to be returned.
  NumericVector cc(n - 1);
  NumericVector dd(n);
  NumericVector x(n);
  cc[0] = c[0] / b[0];
  // Forward sweep for c vector
  for (int i = 1; i < n - 1; ++i)
  {
    cc[i] = c[i] / (b[i] - a[i - 1] * cc[i - 1]);
  }
  // Repeating nearly the same for d
  dd[0] = d[0] / b[0];
  for (int i = 1; i < n; ++i)
  {
    // On paper a_i is written indexed from 2 to n, but since the values in the computer are indexed starting from 0, to n-1, and we are looping from 1 to n-1, we must decrement a_i to a_{i-1}
    // It actually doesn't matter for a_i=a constants
    dd[i] = (d[i] - a[i - 1] * dd[i - 1]) / (b[i] - a[i - 1] * cc[i - 1]);
  }
  // Back substitution for the solution
  x[n - 1] = dd[n - 1];
  for (int i = n - 2; i >= 0; --i)
  {
    x[i] = dd[i] - cc[i] * x[i + 1];
  }
  return x;
}
