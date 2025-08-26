#' Options for saving posterior samples during MCMC
#'
#' @param save_UV Logical; whether to save samples of U and V matrices (default FALSE)
#' @param save_ab Logical; whether to save samples of additive effects a and b (default FALSE)
#' @param thin_UV Integer; thinning interval for U/V samples (default 10)
#' @param thin_ab Integer; thinning interval for a/b samples (default 10)
#' 
#' @return List of posterior saving options
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
posterior_options <- function(
  save_UV = FALSE, save_ab = FALSE,
  thin_UV = 10, thin_ab = 10
  ){
  list(
    save_UV = save_UV,
    save_ab = save_ab,
    thin_UV = thin_UV,
    thin_ab = thin_ab
  )
}

#' Simulate posterior distributions from fitted AME model
#'
#' @param fit Fitted ame model object
#' @param component Character; which component to simulate: "UV", "ab", "beta", "Y"
#' @param n_samples Number of posterior samples to generate
#' @param use_full_cov For UV component, whether to use empirical variance scaling (default FALSE)
#' @param seed Random seed for reproducibility
#' 
#' @return Array or matrix of posterior samples
#' 
#' @details
#' This function can simulate posterior distributions even when they weren't
#' saved during MCMC, by using the posterior means and variance components.
#' 
#' For more accurate posteriors, use posterior_options() during model fitting
#' to save the actual MCMC samples.
#' 
#' @examples
#' \dontrun{
#' # Simulate posterior predictive distribution
#' Y_post <- simulate_posterior(fit, "Y", n_samples = 100)
#' 
#' # Get posterior samples of multiplicative effects
#' UV_post <- simulate_posterior(fit, "UV", n_samples = 100)
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
simulate_posterior <- function(
  fit, 
  component = c("UV", "ab", "beta", "Y"),
  n_samples = 100,
  use_full_cov = FALSE,
  seed = NULL
  ){
  
  component <- match.arg(component)
  
  if(!is.null(seed)) set.seed(seed)
  
  # Check if actual posterior samples are available
  has_samples <- switch(component,
    UV = !is.null(fit$U_samples) && !is.null(fit$V_samples),
    ab = !is.null(fit$a_samples) && !is.null(fit$b_samples),
    beta = !is.null(fit$BETA),
    Y = FALSE  # Always simulate Y
  )
  
  if(has_samples && component != "Y") {
    # Use actual MCMC samples
    cli::cli_text("Using saved MCMC samples for {.field {component}}")
    
    if(component == "UV") {
      return(sample_from_saved_UV(fit, n_samples))
    } else if(component == "ab") {
      return(sample_from_saved_ab(fit, n_samples))
    } else if(component == "beta") {
      return(sample_from_saved_beta(fit, n_samples))
    }
    
  } else {
    # Simulate from approximate posterior
    cli::cli_text("Simulating from approximate posterior for {.field {component}}")
    
    if(component == "UV") {
      return(simulate_UV_posterior(fit, n_samples, use_full_cov))
    } else if(component == "ab") {
      return(simulate_ab_posterior(fit, n_samples))
    } else if(component == "beta") {
      return(simulate_beta_posterior(fit, n_samples))
    } else if(component == "Y") {
      return(simulate_Y_posterior(fit, n_samples))
    }
  }
}

#' Sample from saved U/V posteriors
#' @noRd
sample_from_saved_UV <- function(fit, n_samples) {
  n_saved <- dim(fit$U_samples)[3]
  
  if(n_samples > n_saved) {
    cli::cli_warn("Requested {.val {n_samples}} samples but only {.val {n_saved}} available")
    n_samples <- n_saved
  }
  
  # Sample indices
  idx <- sample(n_saved, n_samples, replace = (n_samples > n_saved))
  
  # UV product is n x n (or nA x nB for bipartite)
  n_row <- nrow(fit$U_samples)
  n_col <- nrow(fit$V_samples)
  
  UV_samples <- array(NA, dim = c(n_row, n_col, n_samples))
  
  for(i in 1:n_samples) {
    UV_samples[,,i] <- fit$U_samples[,,idx[i]] %*% t(fit$V_samples[,,idx[i]])
  }
  
  return(UV_samples)
}

#' Simulate U/V from approximate posterior
#' @noRd
simulate_UV_posterior <- function(fit, n_samples, use_full_cov = FALSE) {
  
  if(is.null(fit$U) || is.null(fit$V)) {
    cli::cli_abort("No multiplicative effects in model")
  }
  
  n <- nrow(fit$U)
  R <- ncol(fit$U)
  
  # Get variance components
  if(!is.null(fit$VC)) {
    s2 <- mean(fit$VC[, ncol(fit$VC)])  # Last column is usually s2
  } else {
    s2 <- 1  # Default variance
  }
  
  # Simulate U and V matrices
  UV_samples <- array(NA, dim = c(n, n, n_samples))
  
  if(use_full_cov && !is.null(fit$BETA) && nrow(fit$BETA) >= 100) {
    # Use empirical variance from MCMC samples for better approximation
    # Scale variance based on number of MCMC samples
    n_mcmc <- nrow(fit$BETA)
    sd_scale <- sqrt(s2 / R) * sqrt(1 + 50/n_mcmc)  # Adaptive scaling
    cli::cli_inform("Using empirical variance scaling based on {.val {n_mcmc}} MCMC samples")
  } else {
    # Simple independent noise (original approach)
    sd_scale <- sqrt(s2 / R)
  }
  
  for(s in 1:n_samples) {
    # Add noise to posterior means
    U_sim <- fit$U + matrix(rnorm(n * R, 0, sd_scale), n, R)
    V_sim <- fit$V + matrix(rnorm(n * R, 0, sd_scale), n, R)
    
    UV_samples[,,s] <- U_sim %*% t(V_sim)
  }
  
  return(UV_samples)
}

#' Sample from saved a/b posteriors
#' @noRd
sample_from_saved_ab <- function(fit, n_samples) {
  n_saved <- ncol(fit$a_samples)
  
  if(n_samples > n_saved) {
    cli::cli_warn("Requested {.val {n_samples}} samples but only {.val {n_saved}} available")
    n_samples <- n_saved
  }
  
  # Sample indices
  idx <- sample(n_saved, n_samples, replace = (n_samples > n_saved))
  
  ab_samples <- list(
    a = fit$a_samples[, idx],
    b = fit$b_samples[, idx]
  )
  
  return(ab_samples)
}

#' Simulate a/b from approximate posterior
#' @noRd
simulate_ab_posterior <- function(fit, n_samples) {
  
  if(is.null(fit$APM) || is.null(fit$BPM)) {
    cli::cli_abort("No additive effects in model")
  }
  
  n <- length(fit$APM)
  
  # Get variance components from VC matrix
  if(!is.null(fit$VC) && ncol(fit$VC) >= 2) {
    Sab_var <- mean(fit$VC[, 1])  # First column is usually Sab
  } else {
    Sab_var <- 1  # Default variance
  }
  
  # Simulate a and b
  ab_samples <- list(
    a = matrix(NA, n, n_samples),
    b = matrix(NA, n, n_samples)
  )
  
  for(s in 1:n_samples) {
    ab_samples$a[, s] <- fit$APM + rnorm(n, 0, sqrt(Sab_var))
    ab_samples$b[, s] <- fit$BPM + rnorm(n, 0, sqrt(Sab_var))
  }
  
  return(ab_samples)
}

#' Sample from saved beta posteriors
#' @noRd
sample_from_saved_beta <- function(fit, n_samples) {
  n_saved <- nrow(fit$BETA)
  
  if(n_samples > n_saved) {
    cli::cli_warn("Requested {.val {n_samples}} samples but only {.val {n_saved}} available")
    n_samples <- n_saved
  }
  
  # Sample indices
  idx <- sample(n_saved, n_samples, replace = (n_samples > n_saved))
  
  return(fit$BETA[idx, , drop = FALSE])
}

#' Simulate beta from approximate posterior  
#' @noRd
simulate_beta_posterior <- function(fit, n_samples) {
  
  if(is.null(fit$BETA)) {
    return(matrix(0, n_samples, 0))  # No covariates
  }
  
  # Use empirical mean and covariance from MCMC samples
  beta_mean <- colMeans(fit$BETA)
  beta_cov <- cov(fit$BETA)
  
  # Simulate from multivariate normal
  if(requireNamespace("MASS", quietly = TRUE)) {
    beta_samples <- MASS::mvrnorm(n_samples, beta_mean, beta_cov)
  } else {
    # Simple version without MASS
    p <- length(beta_mean)
    beta_samples <- matrix(NA, n_samples, p)
    for(i in 1:p) {
      beta_samples[, i] <- rnorm(n_samples, beta_mean[i], sqrt(beta_cov[i, i]))
    }
  }
  
  return(beta_samples)
}

#' Simulate Y from posterior predictive distribution
#' @noRd
simulate_Y_posterior <- function(fit, n_samples) {
  
  n <- nrow(fit$Y)
  family <- fit$family
  
  # Get components
  beta_samples <- simulate_posterior(fit, "beta", n_samples)
  
  # Initialize Y samples
  if(fit$symmetric) {
    Y_samples <- array(NA, dim = c(n, n, n_samples))
  } else {
    Y_samples <- array(NA, dim = c(nrow(fit$Y), ncol(fit$Y), n_samples))
  }
  
  for(s in 1:n_samples) {
    # Construct linear predictor
    EZ <- matrix(0, nrow(fit$Y), ncol(fit$Y))
    
    # Add fixed effects if present
    if(!is.null(fit$X) && ncol(beta_samples) > 0) {
      for(k in 1:dim(fit$X)[3]) {
        EZ <- EZ + beta_samples[s, k] * fit$X[,,k]
      }
    }
    
    # Add additive effects
    if(!is.null(fit$APM) && !is.null(fit$BPM)) {
      ab_sim <- simulate_posterior(fit, "ab", 1)
      EZ <- EZ + outer(ab_sim$a[,1], ab_sim$b[,1], "+")
    }
    
    # Add multiplicative effects
    if(!is.null(fit$U_postmean) && !is.null(fit$V_postmean)) {
      UV_sim <- simulate_posterior(fit, "UV", 1)
      EZ <- EZ + UV_sim[,,1]
    }
    
    # Generate Y based on family
    Y_samples[,,s] <- generate_Y_from_EZ(EZ, family, fit)
  }
  
  return(Y_samples)
}

#' Generate Y from linear predictor EZ
#' @noRd
generate_Y_from_EZ <- function(EZ, family, fit) {
  
  n <- nrow(EZ)
  m <- ncol(EZ)
  
  Y <- switch(family,
    normal = EZ + matrix(rnorm(n * m), n, m),
    binary = matrix(rbinom(n * m, 1, pnorm(EZ)), n, m),
    poisson = matrix(rpois(n * m, exp(pmin(EZ, 10))), n, m),
    # Add other families as needed
    EZ  # Default: return linear predictor
  )
  
  # Handle missing values
  if(!is.null(fit$Y)) {
    Y[is.na(fit$Y)] <- NA
  }
  
  return(Y)
}

#' Extract posterior quantiles for model components
#'
#' @param fit Fitted ame model object
#' @param component Character; which component: "beta", "UV", "ab"
#' @param probs Numeric vector of probabilities for quantiles
#' 
#' @return Matrix or array of posterior quantiles
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' 
#' @export
posterior_quantiles <- function(
  fit, 
  component = c("beta", "UV", "ab"),
  probs = c(0.025, 0.5, 0.975)
  ){
  
  component <- match.arg(component)
  
  if(component == "beta") {
    if(is.null(fit$BETA)) {
      cli::cli_abort("No regression coefficients in model")
    }
    
    quantiles <- apply(fit$BETA, 2, quantile, probs = probs, na.rm = TRUE)
    rownames(quantiles) <- paste0("q", probs * 100)
    
    return(quantiles)
    
  } else if(component == "UV") {
    if(is.null(fit$U_postmean) || is.null(fit$V_postmean)) {
      cli::cli_abort("No multiplicative effects in model")
    }
    
    # Simulate samples if not saved
    UV_samples <- simulate_posterior(fit, "UV", n_samples = 1000)
    
    # Compute quantiles for each dyad
    n <- dim(UV_samples)[1]
    m <- dim(UV_samples)[2]
    
    UV_quantiles <- array(NA, dim = c(n, m, length(probs)))
    
    for(i in 1:n) {
      for(j in 1:m) {
        UV_quantiles[i, j, ] <- quantile(UV_samples[i, j, ], 
                                        probs = probs, na.rm = TRUE)
      }
    }
    
    dimnames(UV_quantiles)[[3]] <- paste0("q", probs * 100)
    
    return(UV_quantiles)
    
  } else if(component == "ab") {
    if(is.null(fit$APM) || is.null(fit$BPM)) {
      cli::cli_abort("No additive effects in model")
    }
    
    # Simulate samples if not saved
    ab_samples <- simulate_posterior(fit, "ab", n_samples = 1000)
    
    a_quantiles <- apply(ab_samples$a, 1, quantile, probs = probs, na.rm = TRUE)
    b_quantiles <- apply(ab_samples$b, 1, quantile, probs = probs, na.rm = TRUE)
    
    rownames(a_quantiles) <- paste0("q", probs * 100)
    rownames(b_quantiles) <- paste0("q", probs * 100)
    
    return(list(a = a_quantiles, b = b_quantiles))
  }
}