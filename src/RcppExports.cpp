// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// run_simulation
NumericMatrix run_simulation(double step, double final_t, int NH, int NM, NumericVector microProbs, double burnin, int seed, NumericVector pars);
RcppExport SEXP zikaABM_run_simulation(SEXP stepSEXP, SEXP final_tSEXP, SEXP NHSEXP, SEXP NMSEXP, SEXP microProbsSEXP, SEXP burninSEXP, SEXP seedSEXP, SEXP parsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    Rcpp::traits::input_parameter< double >::type final_t(final_tSEXP);
    Rcpp::traits::input_parameter< int >::type NH(NHSEXP);
    Rcpp::traits::input_parameter< int >::type NM(NMSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type microProbs(microProbsSEXP);
    Rcpp::traits::input_parameter< double >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    __result = Rcpp::wrap(run_simulation(step, final_t, NH, NM, microProbs, burnin, seed, pars));
    return __result;
END_RCPP
}
