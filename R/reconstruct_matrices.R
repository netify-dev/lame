#' Reconstruct EZ and UVPM matrices from AME model output
#' 
#' @description
#' Helper functions to reconstruct EZ (linear predictor) and 
#' UVPM (posterior mean of UV) matrices that are no longer stored in the model
#' output to save memory.
#' 
#' Note: EZ returns the LINEAR PREDICTOR (η), not the response:
#' - For Gaussian: EZ = η = μ (identity link)
#' - For Poisson: EZ = η = log(λ) (can be negative)  
#' - For Binary: EZ = η = probit^-1(p) (can be any real value)
#' Use YPM for predictions on the response scale.
#' 
#' @param fit Fitted AME model object
#' @param X Covariate array (optional, will use fit$X if available)
#' 
#' @return Reconstructed matrix
#' 
#' @details
#' These matrices are not stored by default to save memory, but can be
#' reconstructed when needed for diagnostics or analysis.
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
reconstruct_EZ <- function(fit, X = NULL) {
  # Get dimensions
  if(!is.null(fit$YPM)) {
    n <- nrow(fit$YPM)
    m <- ncol(fit$YPM)
  } else {
    stop("Cannot determine dimensions from fit object")
  }
  
  # Start with intercept
  if(!is.null(fit$BETA)) {
    beta_mean <- colMeans(fit$BETA)
    
    # Check if X is stored in fit object (for bipartite models)
    if(!is.null(fit$X)) {
      X <- fit$X
    }
    
    if(!is.null(X)) {
      # X includes intercept as first layer
      EZ <- matrix(0, n, m)
      for(k in 1:length(beta_mean)) {
        if(k <= dim(X)[3]) {
          EZ <- EZ + X[,,k] * beta_mean[k]
        }
      }
    } else {
      # No covariates, just intercept
      EZ <- matrix(beta_mean[1], n, m)
    }
  } else {
    EZ <- matrix(0, n, m)
  }
  
  # Add additive effects
  if(!is.null(fit$APM)) {
    if(fit$mode == "bipartite") {
      # Bipartite: row and column effects are different
      for(i in 1:n) EZ[i,] <- EZ[i,] + fit$APM[i]
      if(!is.null(fit$BPM)) {
        for(j in 1:m) EZ[,j] <- EZ[,j] + fit$BPM[j]
      }
    } else if(fit$symmetric) {
      # Symmetric unipartite
      EZ <- EZ + outer(fit$APM, fit$APM, "+")
    } else {
      # Asymmetric unipartite
      if(!is.null(fit$BPM)) {
        EZ <- EZ + outer(fit$APM, fit$BPM, "+")
      }
    }
  }
  
  # Add multiplicative effects
  if(!is.null(fit$U)) {
    if(fit$mode == "bipartite" && !is.null(fit$G) && !is.null(fit$V)) {
      # Bipartite case with G matrix
      EZ <- EZ + fit$U %*% fit$G %*% t(fit$V)
    } else if(fit$symmetric && !is.null(fit$L)) {
      # Symmetric case: U %*% L %*% t(U)
      if(is.matrix(fit$L)) {
        EZ <- EZ + fit$U %*% fit$L %*% t(fit$U)
      } else if(is.numeric(fit$L)) {
        # L is a vector of eigenvalues
        EZ <- EZ + fit$U %*% diag(fit$L, nrow=length(fit$L)) %*% t(fit$U)
      }
    } else if(!is.null(fit$V)) {
      # Asymmetric case: U %*% t(V)
      EZ <- EZ + fit$U %*% t(fit$V)
    }
  }
  
  return(EZ)
}

#' @rdname reconstruct_EZ
#' @export
reconstruct_UVPM <- function(fit) {
  if(is.null(fit$U)) {
    return(NULL)
  }
  
  if(fit$symmetric) {
    # Symmetric case - always reconstruct from U and L
    if(!is.null(fit$L) && !is.null(fit$U)) {
      if(is.matrix(fit$L)) {
        return(fit$U %*% fit$L %*% t(fit$U))
      } else if(is.numeric(fit$L)) {
        # L is a vector of eigenvalues
        return(fit$U %*% diag(fit$L, nrow=length(fit$L)) %*% t(fit$U))
      }
    }
  } else {
    # Asymmetric case
    if(!is.null(fit$V)) {
      if(!is.null(fit$G) && fit$mode == "bipartite") {
        # Bipartite with G - check dimensions
        if(ncol(fit$U) == nrow(fit$G) && ncol(fit$G) == ncol(fit$V)) {
          return(fit$U %*% fit$G %*% t(fit$V))
        } else {
          # G dimensions don't match, use direct product
          return(fit$U %*% t(fit$V))
        }
      } else {
        # Regular asymmetric
        return(fit$U %*% t(fit$V))
      }
    }
  }
  
  return(NULL)
}