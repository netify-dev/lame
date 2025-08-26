#' Simulate networks from a fitted AME model
#' 
#' @description
#' Generates multiple network realizations from the posterior distribution of
#' a fitted AME model. This function performs posterior predictive simulation
#' by drawing from the full joint posterior distribution of model parameters,
#' thereby propagating parameter uncertainty into the simulated networks.
#' 
#' \strong{Key Difference from print():} While \code{print.ame} displays existing
#' model results without computation, \code{simulate.ame} actively generates new
#' network data by sampling from the posterior predictive distribution. This is
#' computationally intensive and produces new datasets for analysis.
#' 
#' @param object fitted model object of class "ame"
#' @param nsim number of networks to simulate (default: 100)
#' @param seed random seed for reproducibility
#' @param newdata optional list containing new covariate data:
#'   \describe{
#'     \item{Xdyad}{dyadic covariates (n x n x p array or nA x nB x p for bipartite)}
#'     \item{Xrow}{row/sender covariates (n x p matrix or nA x p for bipartite)}
#'     \item{Xcol}{column/receiver covariates (n x p matrix or nB x p for bipartite)}
#'   }
#'   If NULL, uses covariates from original model fit
#' @param burn_in number of initial MCMC samples to discard (default: 0, assumes
#'   burn-in already removed)
#' @param thin thinning interval for MCMC samples (default: 1, use every sample)
#' @param return_latent logical: return latent Z matrices in addition to Y? (default: FALSE)
#' @param ... additional arguments (not currently used)
#' 
#' @return A list with components:
#'   \describe{
#'     \item{Y}{list of nsim simulated networks in the same format as the original data}
#'     \item{Z}{if return_latent=TRUE, list of nsim latent Z matrices}
#'     \item{family}{the family of the model (binary, normal, etc.)}
#'     \item{mode}{network mode (unipartite or bipartite)}
#'   }
#'
#' @details
#' \strong{Mathematical Framework:}
#' 
#' The AME model represents networks through a latent variable framework:
#' \deqn{Y_{ij} \sim F(Z_{ij})}
#' where F is the observation model (e.g., probit for binary) and Z is the latent network:
#' \deqn{Z_{ij} = \beta^T x_{ij} + a_i + b_j + u_i^T v_j + \epsilon_{ij}}
#' 
#' Components:
#' \itemize{
#'   \item \eqn{\beta}: regression coefficients for dyadic/nodal covariates
#'   \item \eqn{a_i, b_j}: additive sender and receiver random effects
#'   \item \eqn{u_i, v_j}: multiplicative latent factors (dimension R)
#'   \item \eqn{\epsilon_{ij}}: dyadic random effects with correlation \eqn{\rho}
#' }
#' 
#' \strong{Uncertainty Quantification Process:}
#' 
#' For each simulated network k = 1, ..., nsim:
#' 
#' 1. \strong{Parameter Sampling:} Draw parameter set \eqn{\theta^{(k)}} from MCMC chains:
#'    \itemize{
#'      \item Sample iteration s uniformly from stored MCMC samples
#'      \item Extract \eqn{\beta^{(s)}}, variance components \eqn{(v_a^{(s)}, v_b^{(s)}, v_e^{(s)}, \rho^{(s)})}
#'    }
#' 
#' 2. \strong{Random Effects Generation:} Sample new random effects from posterior distributions:
#'    \itemize{
#'      \item \eqn{a_i^{(k)} \sim N(0, v_a^{(s)})} for i = 1, ..., n (row effects)
#'      \item \eqn{b_j^{(k)} \sim N(0, v_b^{(s)})} for j = 1, ..., m (column effects)
#'      \item Note: We sample fresh from the posterior variance rather than using
#'            point estimates to properly propagate uncertainty
#'    }
#' 
#' 3. \strong{Latent Network Construction:} Build expected latent positions:
#'    \deqn{E[Z_{ij}^{(k)}] = \beta^{(s)T} x_{ij} + a_i^{(k)} + b_j^{(k)} + \hat{u}_i^T \hat{v}_j}
#'    where \eqn{\hat{u}_i, \hat{v}_j} are posterior mean latent factors
#' 
#' 4. \strong{Dyadic Correlation:} Add correlated noise structure:
#'    \deqn{Z_{ij}^{(k)} = E[Z_{ij}^{(k)}] + \epsilon_{ij}^{(k)}}
#'    where \eqn{\epsilon} has covariance structure:
#'    \deqn{Cov(\epsilon_{ij}, \epsilon_{ji}) = \rho^{(s)} v_e^{(s)}}
#'    \deqn{Var(\epsilon_{ij}) = v_e^{(s)}}
#' 
#' 5. \strong{Observation Model:} Generate observed network based on family:
#'    \itemize{
#'      \item Binary: \eqn{Y_{ij}^{(k)} = I(Z_{ij}^{(k)} > 0)}
#'      \item Normal: \eqn{Y_{ij}^{(k)} = Z_{ij}^{(k)}}
#'      \item Poisson: \eqn{Y_{ij}^{(k)} \sim Poisson(\exp(Z_{ij}^{(k)}))}
#'      \item Other families use appropriate link functions
#'    }
#' 
#' \strong{Sources of Uncertainty:}
#' 
#' The simulation captures three types of uncertainty:
#' 
#' 1. \strong{Parameter uncertainty:} Different MCMC samples yield different \eqn{\beta, v_a, v_b, v_e, \rho}
#' 2. \strong{Random effect uncertainty:} Fresh draws from \eqn{N(0, v_a), N(0, v_b)} for each simulation
#' 3. \strong{Dyadic uncertainty:} Correlated random noise \eqn{\epsilon_{ij}}
#' 
#' This approach provides proper posterior predictive distributions that account
#' for all sources of uncertainty in the model. The variation across simulated
#' networks reflects our posterior uncertainty about the data generating process.
#' 
#' \strong{Limitations:}
#' 
#' Currently, multiplicative effects (U, V) use posterior means rather than
#' sampling from their full posterior. For complete uncertainty quantification,
#' one would need to store and sample from the full MCMC chains of these
#' latent factors, which would require substantial additional memory.
#' 
#' \strong{Symmetric Networks:}
#' 
#' For symmetric networks, the model enforces \eqn{a_i = b_i} and \eqn{u_i = v_i},
#' and the latent matrix Z is symmetrized before generating observations.
#'
#' @examples
#' \dontrun{
#' # Fit a model
#' data(YX_bin)
#' fit <- ame(YX_bin$Y, YX_bin$X, burn=100, nscan=500, family="binary")
#' 
#' # Simulate 50 networks from posterior
#' sims <- simulate(fit, nsim=50)
#' 
#' # With new covariates
#' new_X <- array(rnorm(dim(YX_bin$X)[1] * dim(YX_bin$X)[2] * dim(YX_bin$X)[3]),
#'                dim=dim(YX_bin$X))
#' sims_new <- simulate(fit, nsim=50, newdata=list(Xdyad=new_X))
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @importFrom stats simulate
#' @method simulate ame
#' @export
simulate.ame <- function(
  object, nsim = 100, seed = NULL, newdata = NULL,
  burn_in = 0, thin = 1, return_latent = FALSE, ...
  ){
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Extract model information
  fit <- object
  family <- if (!is.null(fit$family)) fit$family else "normal"
  symmetric <- if (!is.null(fit$symmetric)) {
    if (is.na(fit$symmetric)) FALSE else fit$symmetric
  } else FALSE
  mode <- if (!is.null(fit$mode)) fit$mode else "unipartite"
  bip <- (mode == "bipartite")
  
  # Get dimensions
  if (bip) {
    nA <- fit$nA
    nB <- fit$nB
    n_row <- nA
    n_col <- nB
  } else {
    # APM is a vector, not a matrix
    n <- length(fit$APM)
    n_row <- n_col <- n
  }
  
  # Extract MCMC samples
  BETA <- fit$BETA
  VC <- fit$VC
  
  # Apply burn-in and thinning
  if (burn_in > 0 && burn_in < nrow(BETA)) {
    BETA <- BETA[-(1:burn_in), , drop = FALSE]
    VC <- VC[-(1:burn_in), , drop = FALSE]
  }
  if (thin > 1) {
    keep <- seq(1, nrow(BETA), by = thin)
    BETA <- BETA[keep, , drop = FALSE]
    VC <- VC[keep, , drop = FALSE]
  }
  
  # Number of MCMC samples available
  n_mcmc <- nrow(BETA)
  
  # Get covariates (use new data if provided, otherwise use original)
  if (!is.null(newdata)) {
    Xdyad <- newdata$Xdyad
    Xrow <- newdata$Xrow
    Xcol <- newdata$Xcol
    
    # Convert to array format if needed
    if (!is.null(Xdyad) && length(dim(Xdyad)) == 2) {
      Xdyad <- array(Xdyad, dim = c(dim(Xdyad), 1))
    }
    X <- Xdyad  # For compatibility
  } else {
    # Check for stored design matrix components first
    if (!is.null(fit$X)) {
      X <- fit$X
    } else {
      # No stored covariates and no new data provided
      warning("Original covariates not stored in model. Using zero covariates. ",
              "Consider providing covariates via newdata argument.")
      X <- array(0, dim = c(n_row, n_col, 1))
    }
  }
  
  # Prepare storage for simulations
  Y_sims <- vector("list", nsim)
  if (return_latent) {
    Z_sims <- vector("list", nsim)
  }
  
  # Simulate networks
  for (i in 1:nsim) {
    # Sample from posterior (with replacement if nsim > n_mcmc)
    s <- sample(1:n_mcmc, 1)
    
    # Extract parameters for this sample
    beta <- BETA[s, ]
    
    # Extract variance components based on model type
    # Check if VC has column names
    vc_names <- colnames(VC)
    if (!is.null(vc_names)) {
      # Use column names to extract components
      s2 <- if ("ve" %in% vc_names) VC[s, "ve"] else VC[s, ncol(VC)]
      rho <- if ("rho" %in% vc_names) VC[s, "rho"] else 0
      a_var <- if ("va" %in% vc_names) VC[s, "va"] else VC[s, 1]
      b_var <- if ("vb" %in% vc_names) VC[s, "vb"] else a_var
      
      # Ensure variances are positive
      a_var <- max(a_var, 0.01)
      b_var <- max(b_var, 0.01)
      s2 <- max(s2, 0.01)
    } else {
      # Fall back to positional indexing
      if (symmetric) {
        s2 <- VC[s, ncol(VC)]
        rho <- 0
        a_var <- VC[s, 1]
        b_var <- a_var
      } else {
        if (bip) {
          s2 <- VC[s, 3]
          rho <- 0
          a_var <- VC[s, 1]
          b_var <- VC[s, 2]
        } else {
          s2 <- VC[s, ncol(VC)]
          rho <- VC[s, ncol(VC) - 1]
          a_var <- VC[s, 1]
          b_var <- VC[s, 3]
        }
      }
      # Ensure variances are positive in all cases
      a_var <- max(a_var, 0.01)
      b_var <- max(b_var, 0.01)
      s2 <- max(s2, 0.01)
    }
    
    # Sample random effects
    if (bip) {
      # Use stored APM and BPM as means, add noise
      if (!is.null(fit$APM) && length(fit$APM) == nA) {
        a <- rnorm(nA, mean = fit$APM, sd = sqrt(a_var))
      } else {
        a <- rnorm(nA, mean = 0, sd = sqrt(a_var))
      }
      if (!is.null(fit$BPM) && length(fit$BPM) == nB) {
        b <- rnorm(nB, mean = fit$BPM, sd = sqrt(b_var))
      } else {
        b <- rnorm(nB, mean = 0, sd = sqrt(b_var))
      }
    } else {
      # Unipartite
      # For better uncertainty quantification, sample fresh random effects
      # from their posterior distribution rather than centering on posterior mean
      if (!is.null(fit$APM) && length(fit$APM) > 0) {
        # Option 1: Sample around posterior mean (current approach)
        # a <- rnorm(n, mean = fit$APM, sd = sqrt(a_var))
        
        # Option 2: Sample fresh from prior (better uncertainty)
        # This treats each simulation as a new draw from the posterior
        a <- rnorm(n, mean = 0, sd = sqrt(a_var))
        
        # Option 3 (best but not implemented): Store MCMC samples of a, b
        # and sample from those directly
      } else {
        a <- rnorm(n, mean = 0, sd = sqrt(a_var))
      }
      
      if (!isTRUE(symmetric)) {
        if (!is.null(fit$BPM) && length(fit$BPM) > 0) {
          # Sample fresh for better uncertainty
          b <- rnorm(n, mean = 0, sd = sqrt(b_var))
        } else {
          b <- rnorm(n, mean = 0, sd = sqrt(b_var))
        }
      } else {
        b <- a  # Symmetric case
      }
    }
    
    # Compute expected value of Z
    # Start with intercept if present
    if (length(beta) > 0 && any(beta != 0)) {
      if (length(dim(X)) == 3 && dim(X)[3] > 0) {
        # Matrix multiplication for covariates
        EZ <- array(0, dim = c(n_row, n_col))
        for (k in 1:min(dim(X)[3], length(beta))) {
          EZ <- EZ + beta[k] * X[, , k]
        }
      } else {
        EZ <- matrix(beta[1], n_row, n_col)
      }
    } else {
      EZ <- matrix(0, n_row, n_col)
    }
    
    # Add random effects
    EZ <- EZ + outer(a, b, "+")
    
    # Add multiplicative effects if present
    if (!is.null(fit$U) && !is.null(fit$V)) {
      if (bip && !is.null(fit$G)) {
        # Bipartite with G matrix
        EZ <- EZ + fit$U %*% fit$G %*% t(fit$V)
      } else if (!bip) {
        # Unipartite
        EZ <- EZ + fit$U %*% t(fit$V)
      }
    }
    
    # Ensure symmetry if needed
    if (symmetric && !bip) {
      EZ <- (EZ + t(EZ)) / 2
    }
    
    # Simulate based on family
    if (family == "binary") {
      Y_sim <- simY_bin(EZ, rho)
    } else if (family == "normal") {
      Y_sim <- simY_nrm(EZ, rho, s2)
    } else if (family == "ordinal") {
      # Need to handle ordinal specially
      Y_sim <- simY_ord(EZ, rho, fit$ODM)
    } else if (family == "rrl") {
      Y_sim <- simY_rrl(EZ, rho, fit$ODM)
    } else if (family == "poisson") {
      Y_sim <- simY_pois(EZ)
    } else if (family == "frn") {
      # Fixed rank nomination
      odmax <- if (!is.null(fit$odmax)) fit$odmax else rep(Inf, n_row)
      Y_sim <- simY_frn(EZ, rho, odmax)
    } else if (family == "tobit") {
      Y_sim <- simY_tob(EZ, rho, s2)
    } else {
      stop("Family ", family, " not yet implemented for simulation")
    }
    
    # For unipartite, set diagonal to NA
    if (!bip) {
      diag(Y_sim) <- NA
    }
    
    # Store results
    Y_sims[[i]] <- Y_sim
    if (return_latent) {
      Z_sims[[i]] <- simZ(EZ, rho, s2)
    }
  }
  
  # Prepare output
  result <- list(
    Y = Y_sims,
    family = family,
    mode = mode
  )
  
  if (return_latent) {
    result$Z <- Z_sims
  }
  
  class(result) <- "ame.sim"
  return(result)
}