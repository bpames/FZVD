normalizetest <- function(x, mu, sig){
  # # Convert x and mu to matrices (necessary
  # x <- as.matrix(x)
  # mu <- as.matrix(mu)
  
  # Get number of rows.
  n <- nrow(x)
  
  # Center x according to mu.
  print(n)
  
  if (is.vector(x)){ # Have single testing observation, stored as vector.
    x <- x - mu
  }
  else{ # x is n x p data matrix.
    x <- x - rep(x=1, times=n) %*% t(mu)
  }
  
  # Scale ith columns of x by ith entry of sigma
  x <- x %*%  diag(1/sig) 
  
  return(x)
}