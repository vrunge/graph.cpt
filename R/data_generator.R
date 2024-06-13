

#' SPD_matrix
#'
#' @description Generating a symmetric positive definite (SPD) matrix.
#' A diagonal matrix with \code{eigenvalues} on its diagonal is generated.
#' The rotation is performed using an orthogonal matrix obtained by the QR decomposition of a random matrix
#' @param eigenvalues vector of eigenvalues
#' @return One SPD matrix
#' @examples
#' SPD_matrix(sample(5))
SPD_matrix <- function(eigenvalues)
{
  D <- diag(eigenvalues)
  n <- length(eigenvalues)
  random_matrix <- matrix(rnorm(n^2), n, n)
  Q <- qr.Q(qr(random_matrix)) # Q%*%t(Q) = Id
  W <- Q %*% D %*% t(Q)
  return(W)
}


################################################################################
################################################################################
################################################################################


#' dataGenerator_SPD
#'
#' @description Generating time series for change detection in SPD matrices.
#' Each data point is a matrix obtained by the function \code{SPD_matrix}.
#'
#' @param chpts a vector of increasing change-point indices (the last value is data length)
#' @param d dimension of the SPD matrix (d x d)
#' @param eigenvalues_mat matrix of eigenvalues for successive segments, dimension d x nb_seg with nb_seg the number of segments
#' @return a list of SPD matrices of length 'last value in chpts'
#' @examples
#' dataGenerator_SPD(chpts = c(20,40), d = 5, eigenvalues_mat = matrix(c(sample(5),2*sample(5)),5,2))
dataGenerator_SPD <- function(chpts,
                              d = 5,
                              eigenvalues_mat = sample(5))
{
  ############
  ### STOP ###
  ############

  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(!all(chpts > 0)){stop('chpts values are not all positives')}
  if(is.unsorted(chpts, strictly = TRUE)){stop('chpts should be a strictly increasing vector of change-point positions (indices)')}

  if(!is.numeric(eigenvalues_mat)){stop('parameters values are not all numeric')}
  if(length(chpts) != ncol(eigenvalues_mat)){stop('chpts has a length not equal to the number of columns in eigenvalues_mat')}

  ############  data generation   ############

  n <- chpts[length(chpts)]
  nb_segments <- length(chpts)
  repetition <- c(chpts[1], diff(chpts))
  segment <- rep(1:nb_segments, repetition)

  res <-  list()
  for(i in 1:n)
  {
    res[[i]] <- SPD_matrix(eigenvalues_mat[,segment[i]])
  }
  return(res)
}



################################################################################
################################################################################
################################################################################




#' dataGenerator_SPDw
#'
#' @description Generating time series for change detection in SPD matrices with Wishart distribution.
#' Each data point is a matrix obtained by the the realisation of a Wishart distribution.
#'
#' @param chpts a vector of increasing change-point indices (the last value is data length)
#' @param d dimension of the SPD matrix (d x d)
#' @param eigenvalues_mat matrix of eigenvalues for successive segments, dimension d x nb_seg with nb_seg the number of segments
#' @return a list of SPD matrices of length 'last value in chpts'
#' @examples
#' dataGenerator_SPDw(chpts = c(20,40), d = 5, eigenvalues_mat = matrix(c(sample(5),2*sample(5)),5,2))
dataGenerator_SPDw <- function(chpts,
                               d = 5,
                               eigenvalues_mat = sample(5))
{
  ############
  ### STOP ###
  ############

  if(!is.numeric(chpts)){stop('chpts values are not all numeric')}
  if(!all(chpts > 0)){stop('chpts values are not all positives')}
  if(is.unsorted(chpts, strictly = TRUE)){stop('chpts should be a strictly increasing vector of change-point positions (indices)')}

  if(!is.numeric(eigenvalues_mat)){stop('parameters values are not all numeric')}
  if(length(chpts) != ncol(eigenvalues_mat)){stop('chpts has a length not equal to the number of columns in eigenvalues_mat')}

  ############  data generation   ############

  n <- chpts[length(chpts)]
  nb_segments <- length(chpts)
  repetition <- c(chpts[1], diff(chpts))
  segment <- rep(1:nb_segments, repetition)

  res <-  list()
  for(i in 1:n)
  {
    res[[i]] <- SPD_matrix(eigenvalues_mat[,segment[i]])
  }
  return(res)
}
