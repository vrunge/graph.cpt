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
//' @param beta penalty value for transition between two different states/nodes
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
  int n = y.size();
  int d = states.size();
  arma::mat costQ(d, n + 1, arma::fill::zeros);
  arma::mat path(d, n + 1, arma::fill::zeros);

  costQ(0,0) = -beta;

  //costQ = costQ + 0.5;
  // Convert the resulting Armadillo vector to an R vector


  //NumericVector cost(costQ.size());
  //for (int i = 0; i < cost.size(); ++i){cost[i] = i + costQ(i);}


  double res = 0;
  List z = List::create(costQ, res);
  return(z);
}

