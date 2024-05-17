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
        std::cout << i << " - " << j << " : " << costInter(i, j) << std::endl;
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

  //arma::vec best_path(n, arma::fill::zeros);
  //best_path[n] = index_min(costQ.col(n));
  //std::cout << tau << std::endl;

  for (int t = n-1; t > 0; t--)
  {
   // best_path[t] = path(best_path[t+1], t);
  }

    //best_path <- s

    //for(i in (n+1):2)
    //{
    //  u <- path[s,i]
    //  s <- s - (u-1)
    //  if(s == 0){s <- nb}
    //  best_path <- c(s, best_path)
    //}

  //NumericVector vec_states(n);
  //for (int t = 0; t < n; t++){vec_states[t] = best_path[t];}

  // Optionally, initialize the elements of the vector (e.g., with zeros)
  //for(int i = 0; i < n; ++i)
  //{
  //  vec_states[i] = 0;
  //}

  double res = 0;
  List z = List::create(costQ, path, res);
  return(z);
}

