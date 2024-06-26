// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// dist_SPD
double dist_SPD(const arma::mat& S1, const arma::mat& S2);
RcppExport SEXP _graph_cpt_dist_SPD(SEXP S1SEXP, SEXP S2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type S1(S1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type S2(S2SEXP);
    rcpp_result_gen = Rcpp::wrap(dist_SPD(S1, S2));
    return rcpp_result_gen;
END_RCPP
}
// ts_to_SPDts
List ts_to_SPDts(const arma::mat& data, int window_length);
RcppExport SEXP _graph_cpt_ts_to_SPDts(SEXP dataSEXP, SEXP window_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type window_length(window_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(ts_to_SPDts(data, window_length));
    return rcpp_result_gen;
END_RCPP
}
// SPDts_to_dists
arma::mat SPDts_to_dists(const List& SPDts, const List& states);
RcppExport SEXP _graph_cpt_SPDts_to_dists(SEXP SPDtsSEXP, SEXP statesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type SPDts(SPDtsSEXP);
    Rcpp::traits::input_parameter< const List& >::type states(statesSEXP);
    rcpp_result_gen = Rcpp::wrap(SPDts_to_dists(SPDts, states));
    return rcpp_result_gen;
END_RCPP
}
// ts_to_dists
arma::mat ts_to_dists(const arma::mat& data, const List& states, int window_length);
RcppExport SEXP _graph_cpt_ts_to_dists(SEXP dataSEXP, SEXP statesSEXP, SEXP window_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const List& >::type states(statesSEXP);
    Rcpp::traits::input_parameter< int >::type window_length(window_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(ts_to_dists(data, states, window_length));
    return rcpp_result_gen;
END_RCPP
}
// graph_cpt_manifold
List graph_cpt_manifold(const arma::mat& dists, const arma::mat& A, double beta);
RcppExport SEXP _graph_cpt_graph_cpt_manifold(SEXP distsSEXP, SEXP ASEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type dists(distsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_cpt_manifold(dists, A, beta));
    return rcpp_result_gen;
END_RCPP
}
// graph_cpt_mean
List graph_cpt_mean(const arma::vec& y, const arma::mat& A, const arma::vec& states, double beta);
RcppExport SEXP _graph_cpt_graph_cpt_mean(SEXP ySEXP, SEXP ASEXP, SEXP statesSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type states(statesSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_cpt_mean(y, A, states, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_graph_cpt_dist_SPD", (DL_FUNC) &_graph_cpt_dist_SPD, 2},
    {"_graph_cpt_ts_to_SPDts", (DL_FUNC) &_graph_cpt_ts_to_SPDts, 2},
    {"_graph_cpt_SPDts_to_dists", (DL_FUNC) &_graph_cpt_SPDts_to_dists, 2},
    {"_graph_cpt_ts_to_dists", (DL_FUNC) &_graph_cpt_ts_to_dists, 3},
    {"_graph_cpt_graph_cpt_manifold", (DL_FUNC) &_graph_cpt_graph_cpt_manifold, 3},
    {"_graph_cpt_graph_cpt_mean", (DL_FUNC) &_graph_cpt_graph_cpt_mean, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_graph_cpt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
