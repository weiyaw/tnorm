// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// rmvtnorm
Eigen::MatrixXd rmvtnorm(const int n, const Eigen::VectorXd mean, const Eigen::MatrixXd cov, const Eigen::VectorXd initial, const Eigen::MatrixXd F, const Eigen::VectorXd g, const int burn);
RcppExport SEXP _tnorm_rmvtnorm(SEXP nSEXP, SEXP meanSEXP, SEXP covSEXP, SEXP initialSEXP, SEXP FSEXP, SEXP gSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type initial(initialSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type F(FSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type g(gSEXP);
    Rcpp::traits::input_parameter< const int >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvtnorm(n, mean, cov, initial, F, g, burn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tnorm_rmvtnorm", (DL_FUNC) &_tnorm_rmvtnorm, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_tnorm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
