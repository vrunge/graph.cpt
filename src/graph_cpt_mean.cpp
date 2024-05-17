#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title graph.cpt.mean
//'
//' @description change point inference for univariate vector constraint by a graph structure and fixed means
//'
//' @name graph.cpt.mean
//'
//' @param y univariate data vector y
//' @param A transition matrix (filled with O and 1)
//' @param states vector of state means, a fixed mean value for each node
//'
//' @return A list
//'
//' @examples
//' n <- 5
//' p <- 10
//'
//' @export
// [[Rcpp::export]]
List graph_cpt_mean(const arma::vec& y,
                   const arma::mat& A,
                   const arma::vec& states,
                   double beta = 0)
{
   double res = 0;
   List z = List::create(y, res);
   return(z);
}

