#' PENZDA - PENalized Zero-variance Discriminant Analysis
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

penzda <- function(Xt, Yt, D = diag(p), tol=1e-3, maxits=1000, bta=3, quiet=FALSE, type="ball", gamscale=0.5)
{  
  # Check for valid type.
  stopifnot(exprs = {
    (type=="ball" | type=="sphere")
  } ) 
  # Get data set dimensions.
  n <- nrow(Xt)
  p <- ncol(Xt)
  
  # Get number of classes.
  k <- nlevels(Yt)
  print(x= k)
  
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
    print(x=i)
    
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
  
  # Find null basis.
  N <- Null(t(M))
  #N <- Null(M)
  print('Check null space')
  print(norm(M%*%N))
  
  # print(dim(M))
  # print(dim(N))
  # print(norm(t(M)%*%N) )
  # 
  # Compute leading eigenvector of N'*R'*R*N.
  RN <- R%*%N
  print(dim(RN))
  
  if(min(dim(RN)) <= 2){ # k or p <= 2
    # Take SVD (keep first right singular vector)
    SVDres <- svd(x=RN, nu=0, nv = 1)
  }
  else{ 
    # Use Lanczos method to find dominant singular value/right vector.
    SVDres <- svds(A = RN, k =1, nu = 0, nv =1)
  }
  
  # Extract singular value/vector.
  sig <- SVDres$d[1]
  w <- SVDres$v[,1]
  
  Nw <- N %*% w
  
  # Normalize R.
  R <- R/sig
  RN <- RN/sig
  
  # Define d operators.
  if (norm(x= (diag(diag(D)) - D), type='F') < 1e-12){
    print('D is diagonal')
    d <- diag(D) # Extract diagonal of D.
    if (norm(x= as.matrix(d - rep(1,p)), type='F') < 1e-12){ # D = I.
      print('D is identity')
      Dx <- function(x){return(x)}
    }
    else{ # D is diagonal, but not identity.
        Dx <- function(x){return(d*x)}
    }
  }
  else{ # treat D as arbitrary matrix.
    print('D is not diagonal')
    Dx <- function(x){return(D%*%x)}
  }
  
  
  # Initialize variables.
  sols0 = list(x=w, y=Dx(Nw), z=rep(0, p))
  
  # Initialize gam.
  gam[1] <- gamscale*norm(x=RN%*%w, type = "F")^2/norm(x=sols0$y, type="1")
  
  # Initialize discriminant vectors.
  DVs <- matrix(0,p,k-1)
  its <- rep(0, k-1)
  
  # Calculate discriminant vectors in sequence.
  for (i in 1:(k-1)){ # Calculate ith DV.
    
    # # Calculate remaining initial solutions.
    # sols0$x <- w
    # sols0$z <- rep(0,p)
    
    
    # Call ADMM solver.
    ADMMres <- penzdaADMM(R=R, N=N, RN=RN, D=D, 
                          sols0=sols0, gam=gam[i], bta = bta, 
                          tol = tol, maxits = maxits, 
                          type = type, quiet = quiet)
    
    print(norm(ADMMres$y, type= "F"))
    # Save discriminant vector.
    if (norm(ADMMres$y)>1e-12){ # Normalize if not 0 vector.
      print('Nonzero solution')
      DVs[,i] <- ADMMres$y/norm(x=ADMMres$y, type="F")
    }
    else{print('Converged to zero')}
    
    
    # Record number of iterations.
    its[i] <- ADMMres$its

  }
  
  
  
  return(list(DVs = DVs, its = its, classmns = classMeans, k=k, labels=labs, gam = gam))

}


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
  Nx <- N%*%x
  
  print('Gamma')
  print(gam)
  
  
  
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
    print('Iteration')
    print(its)
          
    
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
        y[ix] <- sign(b(ix))
      }
      else{ # Otherwise shrink then normalize.
        y <- vecshrink(v = b, a= gam)
        y <- y/norm(x=y, type = "F")
      }
      
    }
    
    
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # Update x by solution of linear system of optimality conditions.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    xold <- x
    # Update right-hand side.
    print(dim(N))
    b <- t(N) %*% Dtx(bta*y - z)
    
    # Update using the SMW identity.
    # xtmp <- forwardsolve(l = t(V), x = RN %*% b)
    # xtmp <- backsolve(r = V, x = xtmp)
    # x <- 1/bta * b + 1/bta^2 * t(RN) %*% xtmp
    # print(length(x))
    # print(dim(RN))
    x <- solve(a = (beta*diag(length(x)) - (t(RN) %*% RN)), b= b)
    
    print(norm(x - xold, type= "F"))
    
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
    ep = tol * (sqrt(p) + max(norm(x), norm(y)) )
    es = tol * (sqrt(p) + norm(y))
    
    # Check for convergence.
    if(dr < ep && ds < es) {
      if (quiet == FALSE) {
        print('DONEZO! Converged.')
        print(c(its, dr, ep, ds, es))
      }
      
      # Output optimal solution.
      return(list(x=x, y=y, z=z, its=its))
      break  # STOP. Converged. 
    }

  }
  
  # For loop terminates after maximum number of iterations.
  # Output solution after maximum number of iterations.
  # print(c(its, dr, ep, ds, es))
  return(list(x=x, y=y, z=z, b=b, its=its))
  if (quiet == FALSE) {print('DONEZO! Didnt converge')}

}



#' Soft threshholding operator.
#'
#' Applies the shrinkage operator to vector v with threshold a.
#' @param v vector
#' @param a threshold
#' @return s soft threshholded vector.
#' @export

vecshrink <- function(v,a){
  s <- sign(v)*pmax(abs(v)- a, 0 )
  return(s)
}

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
  N <- nrow(Xtest)
  
  # Number of classes in test data.
  #k <- nlevels(Ytest)
  
  # Project the test data onto the space spanned by the discriminant vectors contained in obj.
  proj <- t(obj$DVs) %*% t(Xtest)
  
  # Compute the centroids of the projected training data.
  cent <- t(obj$DVs) %*% obj$classmns
  
  # Initialize matrix of distances to class-means
  dist <- matrix(0, N, obj$k)
  
  # Calculate distances to centroids of projected test observations.
  for (i in 1:N){
    for (j in 1:obj$k){
      dist[i,j] <- norm(x = as.matrix(proj[,i] - cent[,j]), type="F")
    }
  }
  
  # Label test_obs accoring to the closest centroid to its projection
  preds <- max.col(-dist)
  
  # Calculate misclassification rate.
  mcrate <- sum(Ytest != preds)/N
  
  # Count number of nonzero entries.
  l0 <- sum(colSums(obj$DVs != 0))
  # l0 <- sum(apply(obj$DVs, 2, function(c) sum(c!=0) ) )
  
  # Output.
  return(list(mc = mcrate, l0 = l0, preds = preds, cent = cent, dist = dist))
}