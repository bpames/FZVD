#' calcClassMeans - calculates class-means and factorization of between/within-class covariance matrices for use in penzda.
#' 
#' @param Xt training data set (as matrix)
#' @param Yt factor containing labels of training data.
#' @param D dictionary matrix.
#' @param tol stopping tolerance for ADMM scheme.
#' @param maxits maximum number of iterations performed by ADMM.
#' @param bta augmented Lagrangian penalty term.
#' @param quiet density outside planted submatrix
#' @param type constraint type: "ball" or "sphere".
#' @param gamscale scaling parameter controlling sparsity.
#' @return DVs discriminant vectors corresponding to the optimal solution, classMeans for purpose of prediction.
#' @export
#' 
calcClassMeans <- function(Xt, Yt){
  
  # Get data set dimensions.
  n <- nrow(Xt)
  p <- ncol(Xt)
  
  # Get number of classes.
  k <- nlevels(Yt)
  
  # Get class labels.
  labs <- levels(Yt)
  # Initialize classMeans.
  classMeans <- matrix(0, p, k)
  
  # Initialize R.
  R <- matrix(0, k,p)
  
  
  # Initialize gam.
  gam <- rep(0, k-1)
  
  # Initialize W factor (M).
  M <- matrix(0, n, p)
  
  #++++++++++++++++++++++++++++++++++++++++++
  # Make classMeans and M.
  for (i in 1:k){ # For each class.
    
    # Extract training observations in class i.
    classobs <- Xt[Yt == labs[i], ]
    
    # Get class-size
    ni <- nrow(classobs)
    
    # Compute within-class mean.
    classMeans[, i] <- colMeans(classobs)
    
    # Update W/M.
    M[Yt == labs[i], ] <- classobs - rep(1, ni)%*%t(classMeans[,i])
    
    # Update R.
    R[i,] <- sqrt(ni)*t(classMeans[,i])
    
  } # END for i in 1:k
  
  return(list(classMeans = classMeans, k = k, M = M, R = R))
}