

#' transition_matrix
#'
#' @description Generating a transition matrix for the change point problem constraint by this transition matrix
#'
#' @param d dimension of the transition matrix d x d
#' @param degree number of connections from one state if \code{type} is \code{"degree"}
#' @param type type of matrix to generate: \code{"cyclic"}, \code{"all"}, \code{"degree"}
#' @return A square transition matrix of size d by d. A one value in position (i,j) encodes a transition from state i to state j.
#' @examples
#' transition_matrix()
transition_matrix <- function(d = 5, degree = round(d/2), type = "cyclic")
{
  ##########  ##########  ##########  ##########  ##########
  if(type == "cyclic")
  {
    A <- matrix(0, d, d)
    for(k in 1:(d-1))
    {
      A[k,k+1] <- 1
      A[k,k] <- 1
    }
    A[d,1] <- 1
    A[d,d] <- 1
  }
  ##########  ##########  ##########  ##########  ##########
  if(type == "all")
  {
    A <- matrix(1, d, d)
  }
  ##########  ##########  ##########  ##########  ##########
  if(type == "degree")
  {
    A <- matrix(0, d, d)
    for(k in 1:d)
    {
      A[k,k] <- 1
      A[k, sample(setdiff(1:d, k), degree)] <- 1
    }
  }
  return(A)
}
