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

  # Call calcClassMeans.
  cmnsres <- calcClassMeans(Xt = Xt, Yt = Yt)

  classMeans <- cmnsres$classMeans
  M <- cmnsres$M
  R <- cmnsres$R
  k <- cmnsres$k

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
  if (norm(x= (diag(diag(D)) - D), type='F') < 1e-12){ # D is diagonal

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

  # Initialize discriminant vectors.
  DVs <- matrix(0,p,k-1)
  print(k-1)
  its <- rep(0, k-1)

  # Initialize gam.
  gam <- rep(0, k-1)

  # Calculate discriminant vectors in sequence.
  for (i in 1:(k-1)){ # Calculate ith DV.

    # Initialize variables.
    sols0 = list(x=w, y=Dx(Nw), z=rep(0, p))

    # Initialize gam.
    gam[i] <- gamscale*norm(x=RN%*%w, type = "F")^2/(2*norm(x=sols0$y, type="1"))


    # Call ADMM solver.
    ADMMres <- penzdaADMM(R=R, N=N, RN=RN, D=D,
                          sols0=sols0, gam=gam[i], bta = bta,
                          tol = tol, maxits = maxits,
                          type = type, quiet = quiet)

    # Save discriminant vector.
    if (norm(ADMMres$y)>1e-12){ # Normalize if not 0 vector.
      DVs[,i] <- ADMMres$y/norm(x=ADMMres$y, type="F")
    }
    else{print('Converged to zero')}


    # Record number of iterations.
    its[i] <- ADMMres$its

    # Update null-basis.
    if (i < (k-1)){

      # Extract vector to add to rows of W.
      v <- DVs[,i]

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
        SVDres <- svds(A = RN, k =1, nu = 0, nv =1)
      }

      # Extract singular value/vector.
      sig <- SVDres$d[1]
      w <- SVDres$v[,1]

      Nw <- N %*% w

      # Normalize R.
      R <- R/sig
      RN <- RN/sig

    }

  }

  # Output.
  return(list(DVs = DVs, its = its, classMeans = classMeans, k=k, gam = gam))

}
