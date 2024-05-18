
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

#include "distances.h"

//' @title ts_to_SPDts
//'
//' @description Converting a time series into a sequence of SPD matrices (covariance matrices)
//'
//' @name ts_to_SPDts
//'
//' @param data matrix of data, size p x n for p dimensions and n data points
//' @param window_length the number of data points to take for computing each covariance matrix
//'
//' @return A list of SPD matrices
//'
//' @export
// [[Rcpp::export]]
List ts_to_SPDts(const arma::mat& data,
                 int window_length)
{
  int n = data.n_cols;
  int nb_cov =  n / window_length;
  List res(nb_cov);

  for (int t = 0; t < nb_cov; t++)
  {
    res[t] = cov(data.cols(t*window_length, (t+1)*window_length - 1));
  }
  return(res);
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


//' @title SPDts_to_dists
//'
//' @description Computing all the distance between a list of SPD matrices of size n (obtained by the raw data) and a list of reference matrices of size d.
//'
//' @name SPDts_to_dists
//'
//' @param SPDts list of SPD matrices of size n
//' @param states list of reference matrices of size d
//'
//' @return A matrix of distances between SPD matrices of size d x n
//'
//' @export
// [[Rcpp::export]]
arma::mat SPDts_to_dists(const List& SPDts,
                         const List& states)
{
  int n = SPDts.size();
  int d = states.size();
  arma::mat res(d, n, arma::fill::zeros);

  for (int t = 0; t < n; t++)
  {
    for (int i = 0; i < d; i++)
    {
      res(i,t) = dist_SPD(SPDts(t), states(i));
    }
  }
  return(res);
}


/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////



//' @title ts_to_dists
//'
//' @description function combining the functions ts_to_SPDts and SPDts_to_dists
//'
//' @name ts_to_dists
//'
//' @param data matrix of data, size p x n for p dimensions and n data points
//' @param states list of reference matrices of size d
//' @param window_length the number of data points to take for computing each covariance matrix
//'
//' @return A matrix of distances between SPD matrices of size d x n
//'
//' @export
// [[Rcpp::export]]
arma::mat ts_to_dists(const arma::mat& data,
                      const List& states,
                      int window_length)
{
  List  temp = ts_to_SPDts(data, window_length);
  arma::mat res = SPDts_to_dists(temp, states);
  return(res);
}





