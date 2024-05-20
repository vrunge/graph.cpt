#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


arma::mat raiz(const arma::mat& A)
{
  arma::vec D;
  arma::mat B;
  arma::eig_sym(D,B,A);

  unsigned int n = D.n_elem;
  arma::mat G(n,n, arma::fill::zeros);

  for(unsigned int i=0; i<n; i++)
  {
    if(D(i)<0 ) D(i)=0;
    G(i,i)=sqrt(D(i));
  }
  return B*G*B.t();
}



//' @title dist_SPD
//'
//' @description Computing distances between symmetric positive definite matrices (with invariant properties)
//'
//' @name dist_SPD
//'
//' @param S1 a SPD matrix
//' @param S2 a SPD matrix
//'
//' @return the square root of the sum of log(lambda_i)^2 with lambda_i the eigenvalues of "S1^(-1/2) * S2 * S1^(-1/2)"
//'
//' @examples
//' n <- 5
//' p <- 10
//' M1 <- matrix(rnorm(p*n), n, p)
//' M2 <- matrix(rnorm(p*n), n, p)
//' S1 <- M1%*%t(M1)
//' S2 <- M2%*%t(M2)
//' dist_SPD(S1,S2)
//'
//' @export
// [[Rcpp::export]]
double dist_SPD(const arma::mat& S1, const arma::mat& S2)
{
  double res = 0;
  arma::mat S = arma::inv(arma::sqrtmat_sympd(S2)); // inverse of square root of S2 in S
  arma::mat P = S * S1 * S; // product
  arma::cx_vec eigval = arma::eig_gen(P); // eigenvalues of S * S1 * S

  for (std::size_t i = 0; i < eigval.size(); i++)
  {
    //std::cout << "yyy " << i << " : "<< eigval[i] << " zzz " << std::endl;
    res = res + std::pow(log(eigval[i].real()), 2);
  }
  return(sqrt(res));
}


