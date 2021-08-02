// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// dependentLCM_fit_cpp
Rcpp::List dependentLCM_fit_cpp(Rcpp::IntegerMatrix& x_in, Rcpp::List hparams_list, Rcpp::List params_list, int nitr);
RcppExport SEXP _dependentLCM_dependentLCM_fit_cpp(SEXP x_inSEXP, SEXP hparams_listSEXP, SEXP params_listSEXP, SEXP nitrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix& >::type x_in(x_inSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type hparams_list(hparams_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type params_list(params_listSEXP);
    Rcpp::traits::input_parameter< int >::type nitr(nitrSEXP);
    rcpp_result_gen = Rcpp::wrap(dependentLCM_fit_cpp(x_in, hparams_list, params_list, nitr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dependentLCM_dependentLCM_fit_cpp", (DL_FUNC) &_dependentLCM_dependentLCM_fit_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_dependentLCM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
