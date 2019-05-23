#' penzdaADMM - Alternating direction method of multipliers for PENalized Zero-variance Discriminant Analysis
#'
#' @param R factored between-class covariance matrix B= R'R.
#' @param N matrix whose columns are orthonormal basis for null space of within-class covariance matrix W.
#' @param RN saved product RN = R*N.
#' @param sols0 list of initial values of x, y, z.
#' @param gam regularization parameter controlling emphasis of 1-norm penalty.
#' @param maxits maximum number of iterations performed by ADMM.
#' @param bta augmented Lagrangian penalty term.
#' @param type constraint type: "ball" or "sphere".
#' @param quiet toggles whether to display intermediate iteration statistics.
#' @return x, y, z optimal solutions
#' @return its number of iterations required before convergence.
#' @export

penzdaADMM <- function(R, N, RN, D = diag(p), sols0, gam, bta = 3, 
                       tol = 1e-4, maxits = 1000, 
                       type = "ball", quiet = FALSE)
{  
  # Check for valid type.
  stopifnot(exprs = {
    (type=="ball" | type=="sphere")
  } ) 
  
  # Define D multiplication operators.
  if (norm(x= (diag(diag(D)) - D), type='F') < 1e-12){
    
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
  
  # Get number of classes.
  k <- nrow(R)
  
  # Initialize x,y,z.
  x <- sols0$x
  Nx <- N %*% x
  
  # Initialize y and z.
  y <- sols0$y
  z <- sols0$z
  
  # Calculate Cholesky factorization for x-update.
  V <- chol(diag(k) - 1/bta*(RN%*%t(RN)))
  
  # Iteratively update x,y,z.
  for (iter in 1:maxits){
    # Save previous iterate.
    yold <- y
    
    # Update iteration count.
    its <- iter
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update y using soft-thresholding.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (type == "ball"){ # shrinkage then optional normalization.
      # Apply vecshrink.
      y <- vecshrink(v=bta*Dx(Nx) + z, a=gam)
      
      # Normalize.
      tmp <- max(0, norm(x=y, type="F") - bta)
      y <- y/(bta + tmp)
    } 
    else if (type == "sphere"){ # Apply spherical constraint.
      
      # Calculate RHS.
      b <- Dx(Nx) + z
      
      # Calculate largest magnitude entry in b, to decide if will
      # threshold to 0.
      mx <- max(abs(b))
      
      # Update y.
      if (mx <= gam){ # shrinks to 0. Set y to have cardinality 1.
        y <- rep(0,p)
        ix <- which.max(abs(b))
        y[ix] <- sign(b[ix])
      }
      else{ # Otherwise shrink then normalize.
        y <- vecshrink(v = b, a= gam)
        y <- y/norm(x=y, type = "F")
      }
    }
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update x by solution of linear system of optimality conditions.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Save previous iterate
    xold <- x
    
    # Update right-hand side.
    b <- t(N) %*% Dtx(bta*y - z)
    
    # Update using the SMW identity.
    xtmp <- forwardsolve(l = t(V), x = RN %*% b)
    xtmp <- backsolve(r = V, x = xtmp)
    x <- 1/bta * b + 1/bta^2 * t(RN) %*% xtmp
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update z by dual ascent.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Calculate primal feasibility residual.
    Nx <- N %*% x
    r <- Dx(Nx) - y
    
    # Ascent step.
    z <- z + bta * r
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Check stopping criteria.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # Dual residual.
    s <- bta * (y - yold)
    
    # Residual norms.
    ds = norm(x = s, type = "F")
    dr = norm(x = r, type = "F")
    
    # Compute absolute and relative tolerances.
    ep = tol * (sqrt(p) + max(norm(x, type="F"), norm(as.matrix(y),type="F")) )
    es = tol * (sqrt(p) + norm(as.matrix(y), type="F"))
    
    # Check for convergence.
    if(dr < ep && ds < es) {
      if (quiet == FALSE) {
        print('DONEZO! Converged.')
      }
      
      # Output optimal solution.
      return(list(x=x, y=y, z=z, its=its))
      break  # STOP. Converged. 
    }
    
  }
  
  # For loop terminates after maximum number of iterations.
  # Output solution after maximum number of iterations.
  # print(c(its, dr, ep, ds, es))
  if (quiet == FALSE) {print('DONEZO! Didnt converge')}
  return(list(x=x, y=y, z=z, b=b, its=its))
  
  
}