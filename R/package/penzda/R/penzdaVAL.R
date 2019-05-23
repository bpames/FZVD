#' PENZDAVAL - PENalized Zero-variance Discriminant Analysis with Validation.
#'
#' @param Xt training data set (as matrix)
#' @param Yt factor containing labels of training data.
#' @param Xval validation data set (as matrix)
#' @param Yval factor containing labels of validation data.
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

penzdaVAL <- function(Xt, Yt, Xval, Yval, 
                      D = diag(p), gmults, sparsity_level = 0.25,
                      tol=1e-3, maxits=1000, bta=3, 
                      quiet=FALSE, type="ball")
{
  # Check for valid type.
  stopifnot(exprs = {
    (type=="ball" | type=="sphere")
  } ) 
  
  # Call calcClassMeans.
  cmnsres <- calcClassMeans(Xt = Xt, Yt = Yt)
  
  classMeans <- cmnsres$classMeans
  M <- cmnsres$M
  R <- cmnsres$R
  k <- cmnsres$k
  
  # Find null basis.
  N <- Null(t(M))
  
  # Form factor of objective matrix.
  RN <- R%*%N
  
  if(min(dim(RN)) <= 2){ # k or p <= 2
    # Take SVD (keep first right singular vector)
    SVDres <- svd(x=RN, nu=0, nv = 1)
  }
  else{ 
    # Use Lanczos method to find dominant singular value/right vector.
    SVDres <- rARPACK::svds(A = RN, k =1, nu = 0, nv =1)
  }
  
  # Extract singular value/vector.
  sig <- SVDres$d[1]
  w <- SVDres$v[,1]
  
  Nw <- N %*% w
  
  # Normalize R.
  R <- R/sig
  RN <- RN/sig
  
  # Define d operators.
  if (norm(x= (diag(diag(D)) - D), type='F') < 1e-12){ # D is diagonal.
    
    d <- diag(D) # Extract diagonal of D.
    
    if (norm(x= as.matrix(d - rep(1,p)), type='F') < 1e-12){ # D = I.
      Dx <- function(x){return(x)}
      Dtx <- function(x){return(x)}
    }
    else{ # D is diagonal, but not identity.
      Dx <- function(x){return(d*x)}
      Dtx <- function(x){return(d*x)}
    }
  }
  else{ # treat D as arbitrary matrix.
    Dx <- function(x){return(D%*%x)}
    Dtx <- function(x){return(t(D)%*%x)}
  }
  
  # Initialize scores.
  N0 <- N
  RN0 <- RN
  numgams <- length(gmults)
  val_score <- (p+1)*rep(k-1, numgams)
  best_ind <- 1
  triv <- 0
  quietADMM <- TRUE
  
  # Initialize discriminant vectors.
  DVs <- array(0,c(p,k-1, numgams))
  its <- matrix(0, k-1, numgams)
  
  # Initialize gammas.
  gammas <- matrix(0, numgams, k-1);
  gmax <- norm(RN %*% w,type = 'F')^2/(2*norm(Dx(Nw), type='1')); 
  
  # CALCULATE SETS OF DISCRIMINANT VECTORS AND VALIDATE.
  for (i in 1:numgams){
    
    # Reset N matrix and products.
    N <- N0
    RN <- RN0
    
    # Calculate each discriminant vector corresponding to gmult[i]
    for (j in 1:(k-1)){
      
      # Calculate initial solutions.
      if (i == 1){ # First discriminant vector.
        sols0 = list(x = w, y = Dx(N %*%w), z = rep(0,p))
      } # end if i =1.
      else{ # warm-start with last discriminant vector found for this multiplier.
        sols0 = list(x = t(N) %*% Dtx(as.matrix(DVs[,j, i-1])),
                     y = as.matrix(DVs[,j, i-1]), 
                     z = rep(0,p))
      } # end else
      
      # Choose gamma.
      gammas[i, j] <- gmults[i]*norm(x=RN%*% sols0$x, type = "F")^2/(2*norm(x=sols0$y, type="1"))
      
      # Call ADMM solver.
      ADMMres <- penzdaADMM(R=R, N=N, RN=RN, D=D, 
                            sols0=sols0, gam=gammas[i,j], bta = bta, 
                            tol = tol, maxits = maxits, 
                            type = type, quiet = quiet)
      
      # Save discriminant vector.
      if (norm(ADMMres$y)>1e-12){ # Normalize if not 0 vector.
        DVs[,j,i] <- ADMMres$y/norm(x=ADMMres$y, type="F")
      }
      else{print('Converged to zero')}
      
      
      # Record number of iterations.
      its[j,i] <- ADMMres$its
      
      # Update null-basis.
      if (j < (k-1)){
        # Extract vector to add to rows of W.
        v <- DVs[,j,i]
        
        # Call Nupdate.
        N <- Nupdate(N=N, v = v)
        
        # Update RN and initial solution.
        RN <- R %*% N
        
        if(min(dim(RN)) <= 2){ # k or p <= 2
          # Take SVD (keep first right singular vector)
          SVDres <- svd(x=RN, nu=0, nv = 1)
        }
        else{ 
          # Use Lanczos method to find dominant singular value/right vector.
          SVDres <- rARPACK::svds(A = RN, k =1, nu = 0, nv =1)
        }
        
        # Extract singular value/vector.
        sig <- SVDres$d[1]
        w <- SVDres$v[,1]
        
        # Normalize R.
        R <- R/sig
        RN <- RN/sig
        
        
      } # end if Nupdate
      
    } # End for j.
    
    # UPDATE VALIDATION SCORES FOR ITH SET OF DISCRIMINANT VECTORS.
    # Form penzda object for validation
    obji <- list(DVs = DVs[,,i], k=k, classMeans = classMeans)
    
    # Call predict fxn
    stts <- predict(obji, Xtest = Xval, Ytest = Yval)
    
    if (stts$l0 == 0){ # Found trivial solution 
      val_score[i] <- (p+1)*(k-1) 
      triv <- 1 
      break # Found trivial solution.
    } # End triv solution case.
    else if (stts$l0 > sparsity_level*(k-1)*p){ # Have nontrivial solution, but too dense. 
      val_score[i] <- stts$l0  
    } # End dense case.
    else{ # Sufficiently sparse. Use misclassification rate as score.
      val_score[i] <- stts$mc
    } # End sparse case and validation score update.
    
    # Check if validation score is best so far.
    if (val_score[i] < val_score[best_ind]){ # Min score so far.
      best_ind <- i # Update best_ind.
    }
    
  } # End for i.
  
  # OUTPUT.
  return(list(DVs = DVs[,,best_ind], allDVs = DVs, gam = gammas[i,],
              bestind = best_ind, val_score = val_score,
              its = its, classMeans = classMeans, k=k
  ))
  
} # END PENZDAVAL.