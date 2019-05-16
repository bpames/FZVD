#' PENZDA - PENalized Zero-variance Discriminant Analysis
#'
#' @param Xt training data set (as matrix)
#' @param Yt factor containing labels of training data.
#' @param D dictionary matrix.
#' @param tol stopping tolerance for ADMM scheme.
#' @param maxits maximum number of iterations performed by ADMM.
#' @param beta augmented Lagrangian penalty term.
#' @param quiet density outside planted submatrix
#' @param type constraint type: "ball" or "sphere".
#' @param gamscale scaling parameter controlling sparsity.
#' @return DVs discriminant vectors corresponding to the optimal solution, classMeans for purpose of prediction.
#' @export

penzda <- function(Xt, Yt, D = diag(p), tol=1e-3, maxits=1000, beta=3, quiet=FALSE, type="ball", gamscale=0.5)
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
  
  # Initialize gamma.
  gamma <- rep(0, k-1)
  
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
  sols0 = list(x=w, y=Dx(Nw), z=rep(0,k-1))
  
  # Initialize gamma.
  gamma[1] <- gamscale*norm(x=RN%*%w, type = "F")^2/norm(x=sols0$y, type="1")
  
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
                          sols0=sols0, gamma=gamma(i), beta = beta, 
                          tol = tol, maxits = maxits, 
                          type = type, quiet = quiet)
    
    # Save discriminant vector.
    DVs[,i] <- ADMMres$y/norm(x=ADMMres$y, type="F")
    
    # Record number of iterations.
    its[i] <- ADMMres$its

  }
  
  
  
  return(list(cmns = classMeans, k=k, labs=labs, gamma = gamma))

}


#' PENZDA - PENalized Zero-variance Discriminant Analysis
#'
#' @param R factored between-class covariance matrix B= R'R.
#' @param N matrix whose columns are orthonormal basis for null space of within-class covariance matrix W.
#' @param RN saved product RN = R*N.
#' @param sols0 list of initial values of x, y, z.
#' @param gamma regularization parameter controlling emphasis of 1-norm penalty.
#' @param maxits maximum number of iterations performed by ADMM.
#' @param beta augmented Lagrangian penalty term.
#' @param type constraint type: "ball" or "sphere".
#' @param quiet toggles whether to display intermediate iteration statistics.
#' @return x, y, z optimal solutions
#' @return its number of iterations required before convergence.
#' @export

penzdaADMM <- function(R, N, RN, D = diag(p), sols0, gamma, beta = 3, 
                          tol = 1e-4, maxits = 1000, 
                          type = "ball", quiet = FALSE)
{  
}