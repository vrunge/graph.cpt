#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @title graph_cpt_manifold
//'
//' @description change point inference for univariate vector constraint by a graph structure and distances computed on a manifold structure
//'
//' @name graph_cpt_manifold
//'
//' @param dists matrix of distances between data and states, size d x n for d states and n data points
//' @param A transition matrix (filled with O and 1)
//' @param beta penalty value for transition between two different states/nodes
//'
//' @return A list
//'
//' @export
// [[Rcpp::export]]
List graph_cpt_manifold(const arma::mat& dists,
                       const arma::mat& A,
                       double beta = 0)
{
  ///
  /// INITIALIZATION
  ///

   int n = dists.n_cols;
   int d = dists.n_rows;
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
       costQ(j, t + 1) = costInter.row(j).min() + pow(dists(j, t), 2);
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
