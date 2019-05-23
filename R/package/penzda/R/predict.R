#' predict 
#'
#' Predicts class labels 
#' @param obj output from penzda function.
#' @param Xtest validation data set.
#' @param Ytest validation data labels.
#' @return preds 
#' @export

predict <- function(obj, Xtest, Ytest){
  
  # Get number of test observations.
  n <- nrow(Xtest)
  p <- ncol(Xtest)
  
  # Number of classes in test data.
  #k <- nlevels(Ytest)
  
  # Project the test data onto the space spanned by the discriminant vectors contained in obj.
  proj <- t(as.matrix(obj$DVs)) %*% t(Xtest)
  
  # Compute the centroids of the projected training data.
  cent <- t(as.matrix(obj$DVs)) %*% as.matrix(obj$classMeans)
  
  # Initialize matrix of distances to class-means
  dist <- matrix(0, n, obj$k)
  
  # Calculate distances to centroids of projected test observations.
  for (i in 1:n){
    for (j in 1:obj$k){
      dist[i,j] <- norm(x = as.matrix(proj[,i] - cent[,j]), type="F")
    }
  }
  
  # Label test_obs accoring to the closest centroid to its projection
  preds <- max.col(-dist)
  preds <- factor(x = preds, levels = levels(Ytest))
  
  # Calculate misclassification rate.
  mcrate <- sum(Ytest!= preds)/n
  
  # Count number of nonzero entries.
  l0 <- sum(colSums(obj$DVs != matrix(0, nrow = nrow(as.matrix(obj$DVs)), 
                                      ncol = ncol(as.matrix(obj$DVs)))))
  # l0 <- sum(apply(obj$DVs, 2, function(c) sum(c!=0) ) )
  
  # Output.
  return(list(mc = mcrate, l0 = l0, preds = preds, cent = cent, dist = dist))
}