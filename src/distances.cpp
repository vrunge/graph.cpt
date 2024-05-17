#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double dist_SPD(const arma::mat& S1, const arma::mat& S2)
{
  double res = 0;
  arma::mat S = arma::inv(arma::sqrtmat_sympd(S2)); // inverse of square root of S1 in S
  arma::mat P = S * S1 * S; // product
  arma::cx_vec eigval = arma::eig_gen(P); // eigenvalues od S * S1 * S

  for (std::size_t i = 0; i < eigval.size(); i++)
  {
    res = res + std::pow(log(eigval[i].real()), 2);
  }
  return(sqrt(res));
}


