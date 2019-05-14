#' PENZDA - PENalized Zero-variance Discriminant Analysis
#'
#' @param train training data set.
#' @param D dictionary matrix.
#' @param tol stopping tolerance for ADMM scheme.
#' @param maxits maximum number of iterations performed by ADMM.
#' @param beta augmented Lagrangian penalty term.
#' @param quiet density outside planted submatrix
#' @param consttype constraint type: "ball" or "sphere".
#' @param gamscale scaling parameter controlling sparsity.
#' @return DVs discriminant vectors corresponding to the optimal solution, classMeans for purpose of prediction.
#' @export

penzda <- function(train, D, tol=1e-4, maxits=1000, beta=5, quiet=false, consttype="ball", gamscale=0.5)
  
  # Calculate class-means of training data.
  classes <- train[,1] 

}