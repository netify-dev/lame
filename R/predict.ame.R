#' Predict method for AME models
#' 
#' @description
#' Generate predictions from fitted AME models, including point estimates
#' and predictive distributions.
#' 
#' @param object Fitted AME model object
#' @param newdata Optional new data (covariates) for prediction
#' @param type Character; type of prediction:
#'   - "response": predicted values on response scale (default)
#'   - "link": predicted values on link scale
#'   - "distribution": full posterior predictive distribution
#' @param n_samples For type="distribution", number of posterior samples
#' @param include_uncertainty Logical; include parameter uncertainty (default TRUE)
#' @param ... Additional arguments (not used)
#' 
#' @return Depending on type:
#'   - "response"/"link": Matrix of predictions
#'   - "distribution": Array of posterior predictive samples
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @examples
#' \dontrun{
#' # Fit model
#' fit <- ame(Y, X, R = 2)
#' 
#' # Point predictions
#' Y_pred <- predict(fit)
#' 
#' # Predictive distribution
#' Y_dist <- predict(fit, type = "distribution", n_samples = 100)
#' 
#' # Predictions with new covariates
#' Y_new <- predict(fit, newdata = X_new)
#' }
#' 
#' @export
predict.ame <- function(
  object, 
  newdata = NULL,
  type = c("response", "link", "distribution"),
  n_samples = 100,
  include_uncertainty = TRUE, ...
  ){
  
  type <- match.arg(type)
  
  # Extract model components
  n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
  m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))
  
  # Handle newdata if provided
  if(!is.null(newdata)) {
    if(!is.array(newdata) || length(dim(newdata)) != 3) {
      cli::cli_abort("newdata must be a 3-dimensional array")
    }
    if(dim(newdata)[1] != n || dim(newdata)[2] != m) {
      cli::cli_abort("newdata dimensions don't match model ({.val {n}} x {.val {m}})")
    }
    X <- newdata
  } else {
    # Use mean effects from fitted model
    X <- NULL
  }
  
  if(type == "distribution") {
    # Generate posterior predictive distribution
    return(predict_distribution(object, X, n_samples, include_uncertainty))
  } else {
    # Point predictions
    if(!is.null(X)) {
      # Compute linear predictor with new covariates
      beta_mean <- colMeans(object$BETA)
      XB <- array(0, dim = c(n, m))
      
      for(k in 1:dim(X)[3]) {
        XB <- XB + X[,,k] * beta_mean[k+1]  # +1 to skip intercept
      }
      XB <- XB + beta_mean[1]  # Add intercept
      
      # Add random effects if available
      if(!is.null(object$APM)) {
        for(i in 1:n) XB[i,] <- XB[i,] + object$APM[i]
      }
      if(!is.null(object$BPM)) {
        for(j in 1:m) XB[,j] <- XB[,j] + object$BPM[j]
      }
      
      # Add multiplicative effects
      UVPM <- reconstruct_UVPM(object)
      if(!is.null(UVPM)) {
        XB <- XB + UVPM
      }
      
      EZ <- XB
    } else {
      # Use fitted values
      if(type == "link") {
        # Reconstruct EZ from components
        EZ <- reconstruct_EZ(object)
        return(EZ)
      } else {
        # For response scale, transform EZ appropriately
        EZ <- reconstruct_EZ(object)
        # Transform to response scale based on family
        if(object$family == "binary") {
          return(pnorm(EZ))  # Probit link
        } else if(object$family == "normal") {
          return(EZ)  # Identity link
        } else if(object$family == "poisson") {
          return(exp(EZ))  # Log link
        } else if(object$family == "tobit") {
          return(pmax(0, EZ))  # Censored at 0
        } else {
          # Default: return YPM (posterior predictive mean)
          return(object$YPM)
        }
      }
    }
  }
}

#' Generate posterior predictive distribution
#' @noRd
predict_distribution <- function(object, X, n_samples, include_uncertainty) {
  
  n <- ifelse(object$mode == "bipartite", object$nA, nrow(object$YPM))
  m <- ifelse(object$mode == "bipartite", object$nB, ncol(object$YPM))
  
  # Initialize array for predictions
  Y_pred <- array(NA, dim = c(n, m, n_samples))
  
  # Get MCMC samples or simulate from approximate posterior
  if(include_uncertainty && !is.null(object$BETA)) {
    n_mcmc <- nrow(object$BETA)
    if(n_mcmc < n_samples) {
      cli::cli_warn("Only {.val {n_mcmc}} MCMC samples available, using all")
      n_samples <- n_mcmc
    }
    
    # Sample from MCMC iterations
    sample_idx <- sample(n_mcmc, n_samples, replace = (n_samples > n_mcmc))
    
    for(i in 1:n_samples) {
      # Get parameters for this sample
      beta <- object$BETA[sample_idx[i], ]
      
      # Compute linear predictor
      if(!is.null(X)) {
        XB <- array(0, dim = c(n, m))
        for(k in 1:dim(X)[3]) {
          XB <- XB + X[,,k] * beta[k+1]
        }
        XB <- XB + beta[1]
      } else {
        XB <- matrix(beta[1], n, m)  # Just intercept
      }
      
      # Add random effects (using posterior means - could be improved)
      if(!is.null(object$APM)) {
        XB <- XB + object$APM
      }
      if(!is.null(object$BPM)) {
        XB <- sweep(XB, 2, object$BPM, "+")
      }
      
      # Add multiplicative effects
      if(!is.null(object$U) && !is.null(object$V)) {
        # Could sample from U/V posterior if available
        if(!is.null(object$U_samples) && !is.null(object$V_samples)) {
          idx <- min(i, dim(object$U_samples)[3])
          UV <- object$U_samples[,,idx] %*% t(object$V_samples[,,idx])
        } else {
          # Add noise to posterior means
          s2 <- object$VC[sample_idx[i], ncol(object$VC)]
          R <- ncol(object$U)
          U_sim <- object$U + matrix(rnorm(n * R, 0, sqrt(s2/R/10)), n, R)
          V_sim <- object$V + matrix(rnorm(m * R, 0, sqrt(s2/R/10)), m, R)
          UV <- U_sim %*% t(V_sim)
        }
        XB <- XB + UV
      }
      
      # Add observation noise and transform
      s2 <- object$VC[sample_idx[i], ncol(object$VC)]
      Z <- XB + matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
      
      # Transform to response scale
      Y_pred[,,i] <- transform_to_response(Z, object$family)
    }
  } else {
    # No parameter uncertainty - just add observation noise
    EZ <- if(!is.null(X)) {
      # Recompute with new X
      beta_mean <- colMeans(object$BETA)
      XB <- array(0, dim = c(n, m))
      for(k in 1:dim(X)[3]) {
        XB <- XB + X[,,k] * beta_mean[k+1]
      }
      XB + beta_mean[1]
    } else {
      reconstruct_EZ(object)
    }
    
    s2 <- mean(object$VC[, ncol(object$VC)])
    
    for(i in 1:n_samples) {
      Z <- EZ + matrix(rnorm(n * m, 0, sqrt(s2)), n, m)
      Y_pred[,,i] <- transform_to_response(Z, object$family)
    }
  }
  
  return(Y_pred)
}

#' Transform linear predictor to response scale
#' @noRd
transform_to_response <- function(Z, family) {
  switch(family,
    normal = Z,
    binary = pnorm(Z),
    poisson = exp(Z),
    tobit = pmax(0, Z),
    ordinal = Z,  # Would need cutpoints for proper transformation
    cbin = pnorm(Z),
    frn = Z,
    rrl = Z,
    Z  # Default
  )
}

# simulate.ame is defined in simulate.ame.R
# This file contains predict, fitted, and residuals methods only

#' Extract fitted values from AME model
#' 
#' @param object Fitted AME model
#' @param ... Additional arguments
#' 
#' @return Matrix of fitted values
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
fitted.ame <- function(object, ...) {
  return(object$YPM)
}

#' Extract residuals from AME model
#' 
#' @param object Fitted AME model
#' @param type Type of residuals ("response" or "pearson")
#' @param ... Additional arguments
#' 
#' @return Matrix of residuals
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
residuals.ame <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  
  # Get observed values (need to be stored or passed)
  if(!is.null(object$Y_obs)) {
    Y <- object$Y_obs
  } else {
    cli::cli_warn("Original Y not stored in model object, residuals may be incorrect")
    Y <- object$YPM  # Fallback
  }
  
  if(type == "response") {
    return(Y - object$YPM)
  } else {
    # Pearson residuals
    s2 <- mean(object$VC[, ncol(object$VC)])
    return((Y - object$YPM) / sqrt(s2))
  }
}