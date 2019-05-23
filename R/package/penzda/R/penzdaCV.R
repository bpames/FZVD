#' PENZDACV - PENalized Zero-variance Discriminant Analysis with Cross-Validation.
#'
#' @param Xt training data set (as matrix)
#' @param Yt factor containing labels of training data.
#' @param nfolds number of cross-validation folds.
#' @param D dictionary matrix.
#' @param gmults vector of multipliers defining potential regularization parameter values.
#' @param sparsity_level desired minimum sparsity level for validation scoring.
#' @param tol stopping tolerance for ADMM scheme.
#' @param maxits maximum number of iterations performed by ADMM.
#' @param bta augmented Lagrangian penalty term.
#' @param quiet true suppresses display of intermediate outputs.
#' @param type constraint type: "ball" or "sphere".
#' @return val_w discriminant vectors corresponding to best parameter.
#' @return DVs set of all discriminant vectors calculated.
#' @return gamma optimal regularization parameter choice.
#' @return bestind position of multiplier corresponding to best gamma.
#' @return valscores array of validation scores.
#' @return classMeans set of training data class-means.
#' @export

penzdaCV <- function(Xt, Yt, nfolds = 5, D = diag(p), gmults, 
                     sparsity_level = 0.25,
                     tol=1e-3, maxits=1000, bta=3, 
                     quiet=FALSE, type="ball")
{
  # Check for valid type.
  stopifnot(exprs = {
    (type=="ball" | type=="sphere")
  } ) 
  
  # Number of observations and features.
  n <- nrow(Xt)
  p <- ncol(Xt)
  
  # Number of gamma.
  ngam <- length(gmults)
  
  # Matrix of cross-validation scores.
  cvscores <- matrix(p+1, ngam, nfolds)
  
  # PERFORM NFOLDS-CROSS VALIDATION.
  for (f in 1:nfolds){ # nfolds CV.
    
    # Split training data into training/validation sets.
    if (nfolds == n){ # Leave-one-out CV.
      print('LOO-CV')
      Xvt <- Xt[-f,] # Exclude fth observation from training.
      Yvt <- Yt[-f] 
      Xval <- Xt[f,] # Use fth observation for validation.
      Yval <- Yt[f]
    }
    else if (nfolds > n){ # more folds than training observations.
      error('# of folds cannot exceed number of training observations.')
    }
    else{ # Use 1/nfolds fraction of training for validation, reset for training.
      # Sample validation indices.
      nval <- ceiling(n*1/nfolds)
      valinds <-  sample.int(n, size = nval, replace = FALSE)
      
      # Training/testing split.
      Xval <- Xtrain[valinds, ]
      Yval <- Ytrain[valinds]
      Xvt <- Xtrain[-valinds, ]
      Yvt <- Ytrain[-valinds]
    } # END IF: training/val split.
    
    # Normalize training data.
    trainlist <- normalize(x = Xvt)
    Xvt <- trainlist$x
    mu <- trainlist$mu
    sig <- trainlist$sig
    
    # Center/scale validation data.
    Xval <- normalizetest(x=Xval, mu = mu, sig = sig)
    
    # Call penzdaVal to calculate validation scores and discriminant vectors for each gamma.
    fres <- penzdaVAL(Xt = Xvt, Yt = Yvt, Xval = Xval, Yval = Yval, gmults = gmults,
                      sparsity_level = sparsity_level, tol = tol, maxits = maxits,
                      bta = bta, quiet = quiet, type = type)
    
    # Update cvscores.
    cvscores[,f] <- fres$val_score
    
    print(sprintf('Fold = %g', f))
    print(cvscores[,f])
    
  } # END FOR F.
  
  
  
  # DETERMINE BEST PARAMETERS AND SOLVE FOR OPTIMAL DISCRIMINANT VECTORS.
  
  # Calculate average validation scores across folds.
  best_ind <- which.min( rowMeans(cvscores) )
  
  # Calculate discriminant vectors with optimized choice of penalty.
  bestres <- penzda(Xt = Xt, Yt = Yt, D=D, tol = tol, maxits = maxits,
                    bta = bta, quiet = quiet, type=type, 
                    gamscale = gmults[best_ind])
  
  # OUTPUT.
  return(list(DVs = bestres$DVs,
              bestind = best_ind, cvscores = cvscores,
              classMeans = bestres$classMeans, k=bestres$k
  ))
  
}