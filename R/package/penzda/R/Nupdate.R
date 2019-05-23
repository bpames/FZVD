#' NUPDATE - Updates null basis using single Householder reflection.
#'
#' @param N matrix whose columns are orthonormal basis for null space of within-class covariance matrix W.
#' @param v new vector
#' @return N1 matrix whose columns are orthonormal basis for intersection of col(N) and orthogonal complement of v.
#' @export

Nupdate <- function(N,v){
  # Extract column of partial QR factorization (used to find N) to reflect
  # onto x-axis.
  x <- t(N) %*% v

  # dimension of current null-space.
  d <- length(x)

  # Update x to (implicitly) calculate Householder reflector.
  x[1] <- sign(x[1])*norm(x=x, type = 'F') + x[1];
  x <- x/norm(x=x, type = 'F');

  # Calculate N1 using Householder reflection.
  N1 <- N[, 2:d] - 2* (N %*% x) %*% t(x[2:d])

  # Output.
  return(N1)

}
