#' NORMALIZE.
#' Centers data set so that column means are equal to 0 and column variance is equal to 1.
#' @param x data set/
#' @return x normalized data set.
#' @return mu vector of column means of original data set.
#' @return sig vector of column standard deviations of original data set.
#' @export


normalize <- function(x){
  # Convert x to a matrix (if necessary).
  x <- as.matrix(x)

  # Get number of rows.
  n <- nrow(x)

  # Get column means.
  mu <- colMeans(x, na.rm = TRUE)

  # Center x.
  x <- x - rep(1,n) %*% t(as.matrix(mu))

  # Calculate column standard deviations.
  sig <- apply(x, 2, sd)

  # Scale ith columns of x by ith entry of sigma
  x <- x %*%  diag(1/sig)

  return(list(x=x, mu=mu, sig=sig))
}
