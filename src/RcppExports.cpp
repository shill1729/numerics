// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// akkadianRootF
double akkadianRootF(double x, double a);
RcppExport SEXP _numerics_akkadianRootF(SEXP xSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(akkadianRootF(x, a));
    return rcpp_result_gen;
END_RCPP
}
// akkadianRoot
double akkadianRoot(double a, int n);
RcppExport SEXP _numerics_akkadianRoot(SEXP aSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(akkadianRoot(a, n));
    return rcpp_result_gen;
END_RCPP
}
// thomas_algorithm
NumericVector thomas_algorithm(NumericVector a, NumericVector b, NumericVector c, NumericVector d);
RcppExport SEXP _numerics_thomas_algorithm(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(thomas_algorithm(a, b, c, d));
    return rcpp_result_gen;
END_RCPP
}
// uniform_log_util
double uniform_log_util(std::vector<double> interval, double rate, unsigned int n);
RcppExport SEXP _numerics_uniform_log_util(SEXP intervalSEXP, SEXP rateSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(uniform_log_util(interval, rate, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_numerics_akkadianRootF", (DL_FUNC) &_numerics_akkadianRootF, 2},
    {"_numerics_akkadianRoot", (DL_FUNC) &_numerics_akkadianRoot, 2},
    {"_numerics_thomas_algorithm", (DL_FUNC) &_numerics_thomas_algorithm, 4},
    {"_numerics_uniform_log_util", (DL_FUNC) &_numerics_uniform_log_util, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_numerics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
