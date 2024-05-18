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
//' vec <- c(rnorm(5, sd = 1), rnorm(5, mean = 1, sd = 0), rnorm(5, mean = 1.9, sd = 1))
//' states <- c(0,1,2)
//' res <- graph_cpt_mean(y = vec, A = transition_matrix(3),states = c(0,1,2), beta = 0)
//' res
//' @export
// [[Rcpp::export]]
List graph_cpt_mean(const arma::vec& y,
                    const arma::mat& A,
                    const arma::vec& states,
                    double beta = 0)
{
  ///
  /// INITIALIZATION
  ///

  int n = y.size();
  int d = states.size();
  arma::mat costQ(d, n + 1, arma::fill::zeros);
  arma::mat path(d, n + 1, arma::fill::zeros);

  for (int i = 0; i < d; i++){ costQ(i,0) = 0;} // INITIAL COST

  ///
  /// DYNAMIC PROGRAMMING
  ///

  arma::mat costInter(d, d, arma::fill::zeros);

  for (int t = 0; t < n; t++)
  {
    for (int i = 0; i < d; i++)
    {
      costInter(i, i) = costQ(i,t);
      for (int j = 0; j < d; j++)
      {
        if(i != j)
        {
          if(A(j,i) == 1){costInter(i, j) = costQ(j,t) + beta;}
          else{costInter(i, j) = std::numeric_limits<double>::infinity();}
        }
        //std::cout << i << " - " << j << " : " << costInter(i, j) << std::endl;
      }
    }

    for (int j = 0; j < d; j++)
    {
      costQ(j, t + 1) = costInter.row(j).min() + pow(y[t] - states[j], 2);
    }
    for (int j = 0; j < d; j++)
    {
      path(j, t + 1) = index_min(costInter.row(j));
    }
  }


  ///
  /// BACKTRACKING
  ///

  NumericVector best_path(n);
  best_path[n-1] = index_min(costQ.col(n));

  for (int t = n-2; t >= 0; t--)
  {
    best_path[t] = path(best_path[t + 1], t + 2);
  }
  List res = List::create(_["costQ"] = costQ.cols(1, costQ.n_cols - 1),
                          _["path"] = path.cols(1, path.n_cols - 1) + 1,
                          _["best_path"] = best_path + 1);
  return(res);
}

