#' Use optimized C++ functions when available
#' 
#' This file contains wrapper functions that automatically use
#' optimized C++ implementations when available, falling back
#' to R implementations if needed.

#' @keywords internal
use_cpp_optimization <- function() {
  getOption("lame.use_cpp", TRUE)
}

#' Optimized precomputeX with C++ backend
#' @keywords internal
precomputeX_opt <- function(X) {
  if(use_cpp_optimization()) {
    tryCatch({
      result <- precomputeX_cpp(X)
      # Add attributes to X
      attributes(X)$Xr <- result$Xr
      attributes(X)$Xc <- result$Xc
      attributes(X)$mX <- result$mX
      attributes(X)$mXt <- result$mXt
      attributes(X)$XX <- result$XX
      attributes(X)$XXt <- result$XXt
      return(X)
    }, error = function(e) {
      # Fall back to R version
      return(precomputeX(X))
    })
  } else {
    return(precomputeX(X))
  }
}

#' Optimized llsrmRho with C++ backend
#' @keywords internal
llsrmRho_opt <- function(Y, Sab, rhos, s2 = 1) {
  if(use_cpp_optimization()) {
    tryCatch({
      return(llsrmRho_cpp(Y, Sab, rhos, s2))
    }, error = function(e) {
      # Fall back to R version
      return(llsrmRho(Y, Sab, rhos, s2))
    })
  } else {
    return(llsrmRho(Y, Sab, rhos, s2))
  }
}

#' Optimized rbeta_ab_fc with C++ backend
#' @keywords internal
rbeta_ab_fc_opt <- function(Z, Sab, rho, X = NULL, s2 = 1, offset = 0,
                           iV0 = NULL, m0 = NULL, g = length(Z)) {
  if(use_cpp_optimization()) {
    tryCatch({
      # Prepare inputs
      if(is.null(X)) {
        X <- design_array(intercept = FALSE, n = nrow(Z))
      }
      if(!is.array(offset)) {
        offset <- matrix(offset, nrow(Z), ncol(Z))
      }
      
      # Set priors
      p <- dim(X)[3]
      if(p > 0 & is.null(iV0)) {
        XX <- attributes(X)$XX
        if(is.null(XX)) {
          X <- precomputeX_opt(X)
          XX <- attributes(X)$XX
        }
        iV0 <- XX/g + diag(diag(XX), nrow = nrow(XX))/g^2
        
        if(all(attributes(X)$mX[,1] == 1)) {
          V0 <- solve(iV0)
          V0[1,1] <- V0[1,1] + sqrt(g) - g/nrow(Z)^2
          iV0 <- solve(V0)
        }
      }
      
      if(is.null(m0)) {
        m0 <- rep(0, dim(X)[3])
      }
      
      return(rbeta_ab_fc_cpp(Z, Sab, rho, X, s2, offset, iV0, m0, g))
    }, error = function(e) {
      # Fall back to R version
      return(rbeta_ab_fc(Z, Sab, rho, X, s2, offset, iV0, m0, g))
    })
  } else {
    return(rbeta_ab_fc(Z, Sab, rho, X, s2, offset, iV0, m0, g))
  }
}

#' Optimized ldZgbme with C++ backend
#' @keywords internal
ldZgbme_opt <- function(Z, Y, llYZ, EZ, rho, s2 = 1) {
  # For normal likelihood, use optimized C++ version
  if(use_cpp_optimization() && is.null(llYZ)) {
    tryCatch({
      return(ldZgbme_opt_cpp(Z, Y, EZ, rho, s2))
    }, error = function(e) {
      # Fall back to R version
      return(ldZgbme(Z, Y, function(y,z) { -0.5*(y-z)^2/s2 }, EZ, rho, s2))
    })
  } else {
    # For custom likelihood functions, use R version
    return(ldZgbme(Z, Y, llYZ, EZ, rho, s2))
  }
}

#' Optimized rrho_fc with C++ backend
#' @keywords internal
rrho_fc_opt <- function(Z, Sab, s2 = 1, offset = 0, ngp = 100, asp = TRUE) {
  if(use_cpp_optimization()) {
    tryCatch({
      if(!is.matrix(offset)) {
        offset <- matrix(offset, nrow(Z), ncol(Z))
      }
      return(rrho_fc_cpp(Z, Sab, s2, offset, ngp, asp))
    }, error = function(e) {
      # Fall back to R version
      return(rrho_fc(Z, Sab, s2, offset, ngp, asp))
    })
  } else {
    return(rrho_fc(Z, Sab, s2, offset, ngp, asp))
  }
}

#' Optimized array_to_list with C++ backend
#' @keywords internal
array_to_list_opt <- function(arr, actorByYr, pdLabs) {
  if(use_cpp_optimization()) {
    tryCatch({
      return(array_to_list_cpp(arr, actorByYr, pdLabs))
    }, error = function(e) {
      # Fall back to R version
      return(array_to_list(arr, actorByYr, pdLabs))
    })
  } else {
    return(array_to_list(arr, actorByYr, pdLabs))
  }
}