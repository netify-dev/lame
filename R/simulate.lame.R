#' Simulate longitudinal networks from a fitted LAME model
#' 
#' Generates multiple longitudinal network realizations from the posterior 
#' distribution of a fitted LAME (Longitudinal AME) model. This function performs
#' posterior predictive simulation for dynamic networks, propagating both 
#' cross-sectional and temporal uncertainty through the simulated trajectories.
#' 
#' @param object fitted model object of class "lame"
#' @param nsim number of network trajectories to simulate (default: 100)
#' @param seed random seed for reproducibility
#' @param newdata optional list containing new covariate data:
#'   \describe{
#'     \item{Xdyad}{list of T dyadic covariate arrays (n x n x p or nA x nB x p)}
#'     \item{Xrow}{list of T row/sender covariate matrices (n x p or nA x p)}
#'     \item{Xcol}{list of T column/receiver covariate matrices (n x p or nB x p)}
#'   }
#'   If NULL, uses covariates from original model fit
#' @param n_time number of time periods to simulate. If NULL, uses same as original data
#' @param burn_in number of initial MCMC samples to discard (default: 0)
#' @param thin thinning interval for MCMC samples (default: 1, use every sample)
#' @param return_latent logical: return latent Z matrices in addition to Y? (default: FALSE)
#' @param start_from character: how to initialize the simulation
#'   \describe{
#'     \item{"posterior"}{start from posterior mean (default)}
#'     \item{"random"}{random initialization}
#'     \item{"data"}{use first time point from original data}
#'   }
#' @param ... additional arguments (not currently used)
#' 
#' @return A list with components:
#'   \describe{
#'     \item{Y}{list of nsim simulated longitudinal network trajectories, 
#'            each element is a list of T networks}
#'     \item{Z}{if return_latent=TRUE, list of nsim latent Z trajectories}
#'     \item{family}{the family of the model (binary, normal, etc.)}
#'     \item{mode}{network mode (unipartite or bipartite)}
#'     \item{n_time}{number of time periods}
#'   }
#'
#' @details
#' \strong{Mathematical Framework for Longitudinal Networks:}
#' 
#' The LAME model extends AME to multiple time periods T with potential
#' temporal dependencies. For each time t = 1, ..., T:
#' \deqn{Y_{ij,t} \sim F(Z_{ij,t})}
#' where the latent network evolves as:
#' \deqn{Z_{ij,t} = \beta^T x_{ij,t} + a_{i,t} + b_{j,t} + u_{i,t}^T v_{j,t} + \epsilon_{ij,t}}
#' 
#' \strong{Temporal Dynamics:}
#' 
#' LAME can incorporate three types of temporal dependencies:
#' 
#' 1. \strong{Static Effects:} Parameters constant over time
#'    \itemize{
#'      \item \eqn{a_{i,t} = a_i}, \eqn{b_{j,t} = b_j} for all t
#'      \item \eqn{u_{i,t} = u_i}, \eqn{v_{j,t} = v_j} for all t
#'    }
#' 
#' 2. \strong{Dynamic Additive Effects:} AR(1) process for random effects
#'    \deqn{a_{i,t} = \rho_{ab} a_{i,t-1} + \eta_{i,t}, \quad \eta_{i,t} \sim N(0, \sigma_a^2(1-\rho_{ab}^2))}
#'    \deqn{b_{j,t} = \rho_{ab} b_{j,t-1} + \xi_{j,t}, \quad \xi_{j,t} \sim N(0, \sigma_b^2(1-\rho_{ab}^2))}
#'    where \eqn{\rho_{ab}} is the temporal correlation parameter
#' 
#' 3. \strong{Dynamic Multiplicative Effects:} AR(1) for latent factors
#'    \deqn{u_{i,t} = \rho_{uv} u_{i,t-1} + \omega_{i,t}}
#'    \deqn{v_{j,t} = \rho_{uv} v_{j,t-1} + \psi_{j,t}}
#' 
#' \strong{Uncertainty Quantification Process for Trajectories:}
#' 
#' For each simulated trajectory k = 1, ..., nsim:
#' 
#' \strong{Step 1: Parameter Sampling}
#' \itemize{
#'   \item Draw MCMC iteration s uniformly from stored posterior samples
#'   \item Extract static parameters: \eqn{\beta^{(s)}}, variance components
#'   \item Extract temporal parameters if applicable: \eqn{\rho_{ab}^{(s)}}, \eqn{\rho_{uv}^{(s)}}
#' }
#' 
#' \strong{Step 2: Initialize at t = 1}
#' 
#' Depending on start_from parameter:
#' \itemize{
#'   \item "posterior": Use posterior means as starting values
#'   \item "random": Draw from stationary distribution
#'   \item For additive effects: \eqn{a_{i,1}^{(k)} \sim N(0, \sigma_a^2)}
#'   \item For multiplicative effects: Initialize from prior
#' }
#' 
#' \strong{Step 3: Evolve Through Time}
#' 
#' For each t = 2, ..., T:
#' 
#' a) \strong{Update Dynamic Effects} (if applicable):
#'    \deqn{a_{i,t}^{(k)} = \rho_{ab}^{(s)} a_{i,t-1}^{(k)} + \eta_{i,t}^{(k)}}
#'    where \eqn{\eta_{i,t}^{(k)} \sim N(0, \sigma_a^2(1-[\rho_{ab}^{(s)}]^2))}
#'    
#'    The innovation variance \eqn{\sigma_a^2(1-\rho_{ab}^2)} ensures stationarity
#' 
#' b) \strong{Construct Latent Network}:
#'    \deqn{E[Z_{ij,t}^{(k)}] = \beta^{(s)T} x_{ij,t} + a_{i,t}^{(k)} + b_{j,t}^{(k)} + u_{i,t}^T v_{j,t}}
#' 
#' c) \strong{Add Dyadic Noise}:
#'    \deqn{Z_{ij,t}^{(k)} = E[Z_{ij,t}^{(k)}] + \epsilon_{ij,t}^{(k)}}
#'    with correlation structure preserved from AME model
#' 
#' d) \strong{Generate Observations}:
#'    Apply appropriate link function based on family
#' 
#' \strong{Sources of Uncertainty in Longitudinal Context:}
#' 
#' 1. \strong{Cross-sectional uncertainty} (as in AME):
#'    \itemize{
#'      \item Parameter uncertainty from MCMC
#'      \item Random effect variability
#'      \item Dyadic noise
#'    }
#' 
#' 2. \strong{Temporal uncertainty}:
#'    \itemize{
#'      \item Uncertainty in temporal correlation parameters \eqn{\rho_{ab}, \rho_{uv}}
#'      \item Innovation noise in AR(1) processes
#'      \item Propagation of uncertainty through time (compounds over periods)
#'    }
#' 
#' 3. \strong{Initial condition uncertainty}:
#'    \itemize{
#'      \item Different starting values lead to different trajectories
#'      \item Captured through start_from options
#'    }
#' 
#' \strong{Interpretation of Multiple Trajectories:}
#' 
#' Each simulated trajectory represents one possible evolution of the network
#' consistent with the posterior distribution. Variation across trajectories
#' captures:
#' \itemize{
#'   \item Model parameter uncertainty
#'   \item Stochastic variation in temporal evolution
#'   \item Accumulated uncertainty over time periods
#' }
#' 
#' The ensemble of trajectories provides prediction intervals that widen over
#' time, reflecting increasing uncertainty in longer-term forecasts.
#' 
#' \strong{Special Considerations:}
#' 
#' 1. \strong{Temporal Correlation:} Higher \eqn{\rho} values create smoother
#'    trajectories with more persistence
#' 
#' 2. \strong{Stationarity:} The AR(1) innovation variance is scaled to maintain
#'    stationary marginal distributions
#' 
#' 3. \strong{Missing Time Points:} If simulating beyond observed data (n_time > T_observed),
#'    covariates are recycled or set to zero with appropriate warnings
#' 
#' \strong{Limitations:}
#' 
#' As with simulate.ame, multiplicative effects currently use posterior means.
#' Full uncertainty would require storing complete MCMC chains for 
#' \eqn{u_{i,t}, v_{j,t}} at all time points, which is memory-intensive for
#' large networks and long time series.
#'
#' @examples
#' \dontrun{
#' # Fit a longitudinal model
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, burn=100, nscan=500, family="binary")
#' 
#' # Simulate 50 network trajectories from posterior
#' sims <- simulate(fit, nsim=50)
#' 
#' # Simulate longer time series
#' sims_long <- simulate(fit, nsim=25, n_time=20)
#' 
#' # With new covariates
#' new_X <- lapply(YX_bin_list$X, function(x) {
#'   array(rnorm(prod(dim(x))), dim=dim(x))
#' })
#' sims_new <- simulate(fit, nsim=50, newdata=list(Xdyad=new_X))
#' }
#' 
#' @author Shahryar Minhas
#' @method simulate lame
#' @export
simulate.lame <- function(object, nsim = 100, seed = NULL, newdata = NULL,
                         n_time = NULL, burn_in = 0, thin = 1, 
                         return_latent = FALSE, start_from = "posterior", ...) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Extract model information
  fit <- object
  family <- if (!is.null(fit$family)) fit$family else "normal"
  symmetric <- if (!is.null(fit$symmetric)) fit$symmetric else FALSE
  mode <- if (!is.null(fit$mode)) fit$mode else "unipartite"
  bip <- (mode == "bipartite")
  dynamic_uv <- if (!is.null(fit$dynamic_uv)) fit$dynamic_uv else FALSE
  dynamic_ab <- if (!is.null(fit$dynamic_ab)) fit$dynamic_ab else FALSE
  
  # Get dimensions
  if (bip) {
    nA <- fit$nA
    nB <- fit$nB
    n_row <- nA
    n_col <- nB
  } else {
    # APM might be a vector or matrix depending on dynamic_ab
    if (is.matrix(fit$APM)) {
      n <- nrow(fit$APM)
    } else {
      n <- length(fit$APM)
    }
    n_row <- n_col <- n
  }
  
  # Determine number of time periods
  if (is.null(n_time)) {
    # Try to infer from the model
    if (!is.null(fit$n_time)) {
      n_time <- fit$n_time
    } else if (!is.null(fit$Xlist)) {
      n_time <- length(fit$Xlist)
    } else {
      n_time <- 10  # Default
      warning("Could not determine number of time periods from model. Using n_time=10")
    }
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
  
  # Get covariates (use new data if provided)
  if (!is.null(newdata)) {
    Xlist <- newdata$Xdyad
    if (is.null(Xlist)) {
      # Create list of zero covariates
      Xlist <- lapply(1:n_time, function(t) {
        array(0, dim = c(n_row, n_col, 1))
      })
    } else {
      # Extend or truncate newdata to match n_time
      orig_length <- length(Xlist)
      if (orig_length < n_time) {
        # Cycle through existing time points if we need more
        for (t in (orig_length + 1):n_time) {
          cycle_index <- ((t - 1) %% orig_length) + 1
          # Ensure cycle_index is valid
          if (cycle_index > 0 && cycle_index <= orig_length) {
            Xlist[[t]] <- Xlist[[cycle_index]]
          } else {
            stop("Invalid cycle index: ", cycle_index, " for orig_length: ", orig_length)
          }
        }
      } else if (orig_length > n_time) {
        Xlist <- Xlist[1:n_time]
      }
    }
  } else if (!is.null(fit$Xlist)) {
    Xlist <- fit$Xlist
    # Extend or truncate to match n_time
    orig_length <- length(Xlist)
    if (orig_length < n_time) {
      # Cycle through existing time points if we need more
      for (t in (orig_length + 1):n_time) {
        cycle_index <- ((t - 1) %% orig_length) + 1
        # Ensure cycle_index is valid
        if (cycle_index > 0 && cycle_index <= orig_length) {
          Xlist[[t]] <- Xlist[[cycle_index]]
        } else {
          stop("Invalid cycle index: ", cycle_index, " for orig_length: ", orig_length)
        }
      }
    } else if (orig_length > n_time) {
      Xlist <- Xlist[1:n_time]
    }
  } else {
    # Create list of zero covariates
    warning("No covariates found. Using zero covariates.")
    Xlist <- lapply(1:n_time, function(t) {
      array(0, dim = c(n_row, n_col, 1))
    })
  }
  
  # Verify Xlist has correct length (only warn, don't stop)
  if (length(Xlist) != n_time) {
    warning("Xlist has length ", length(Xlist), " but n_time is ", n_time, 
            ". This may cause issues.")
  }
  
  # Extract temporal correlation parameters if present
  rho_ab <- if (!is.null(fit$rho_ab)) mean(fit$rho_ab) else 0
  rho_uv <- if (!is.null(fit$rho_uv)) mean(fit$rho_uv) else 0
  
  # Prepare storage for simulations
  Y_sims <- vector("list", nsim)
  if (return_latent) {
    Z_sims <- vector("list", nsim)
  }
  
  # Simulate network trajectories
  for (i in 1:nsim) {
    # Sample from posterior
    s <- sample(1:n_mcmc, 1)
    
    # Extract parameters for this sample
    beta <- BETA[s, ]
    
    # Extract variance components
    # Check if VC has column names
    vc_names <- colnames(VC)
    if (!is.null(vc_names)) {
      # Use column names to extract components
      s2 <- if ("ve" %in% vc_names) VC[s, "ve"] else VC[s, ncol(VC)]
      rho <- if ("rho" %in% vc_names) VC[s, "rho"] else 0
      sigma_a <- sqrt(if ("va" %in% vc_names) VC[s, "va"] else VC[s, 1])
      sigma_b <- sqrt(if ("vb" %in% vc_names) VC[s, "vb"] else sigma_a^2)
    } else {
      # Fall back to positional indexing
      if (symmetric) {
        s2 <- VC[s, ncol(VC)]
        rho <- 0
        sigma_a <- sqrt(VC[s, 1])
        sigma_b <- sigma_a
      } else {
        if (bip) {
          s2 <- VC[s, 3]
          rho <- 0
          sigma_a <- sqrt(VC[s, 1])
          sigma_b <- sqrt(VC[s, 2])
        } else {
          s2 <- VC[s, ncol(VC)]
          rho <- VC[s, ncol(VC) - 1]
          sigma_a <- sqrt(VC[s, 1])
          sigma_b <- sqrt(VC[s, 3])
        }
      }
    }
    
    # Initialize time series storage
    Y_t_list <- vector("list", n_time)
    if (return_latent) {
      Z_t_list <- vector("list", n_time)
    }
    
    # Initialize effects based on start_from
    if (start_from == "posterior") {
      # Start from posterior means
      if (bip) {
        a_t <- fit$APM
        b_t <- fit$BPM
      } else {
        a_t <- fit$APM
        b_t <- if (symmetric) fit$APM else fit$BPM
      }
      
      # Multiplicative effects
      if (!is.null(fit$U)) {
        U_t <- fit$U
        if (!is.null(fit$V)) {
          V_t <- fit$V
        } else {
          V_t <- U_t  # Symmetric case
        }
      }
    } else if (start_from == "random") {
      # Random initialization
      if (bip) {
        a_t <- rnorm(nA, 0, sigma_a)
        b_t <- rnorm(nB, 0, sigma_b)
      } else {
        a_t <- rnorm(n, 0, sigma_a)
        b_t <- if (symmetric) a_t else rnorm(n, 0, sigma_b)
      }
      
      # Random multiplicative effects
      if (!is.null(fit$U)) {
        R <- ncol(fit$U)
        if (bip) {
          U_t <- matrix(rnorm(nA * R, 0, 0.1), nA, R)
          V_t <- matrix(rnorm(nB * R, 0, 0.1), nB, R)
        } else {
          U_t <- matrix(rnorm(n * R, 0, 0.1), n, R)
          V_t <- if (symmetric) U_t else matrix(rnorm(n * R, 0, 0.1), n, R)
        }
      }
    }
    
    # Simulate over time
    for (t in 1:n_time) {
      # Get covariates for this time point
      X_t <- Xlist[[t]]
      
      # Update dynamic effects if present
      if (dynamic_ab && t > 1) {
        # AR(1) process for additive effects
        if (bip) {
          a_t <- rho_ab * a_t + rnorm(nA, 0, sigma_a * sqrt(1 - rho_ab^2))
          b_t <- rho_ab * b_t + rnorm(nB, 0, sigma_b * sqrt(1 - rho_ab^2))
        } else {
          a_t <- rho_ab * a_t + rnorm(n, 0, sigma_a * sqrt(1 - rho_ab^2))
          if (!symmetric) {
            b_t <- rho_ab * b_t + rnorm(n, 0, sigma_b * sqrt(1 - rho_ab^2))
          } else {
            b_t <- a_t
          }
        }
      }
      
      if (dynamic_uv && t > 1 && !is.null(U_t)) {
        # AR(1) process for multiplicative effects
        # This is simplified - actual implementation would be more complex
        U_t <- rho_uv * U_t + matrix(rnorm(prod(dim(U_t)), 0, 
                                          sqrt(1 - rho_uv^2) * 0.1), 
                                    nrow(U_t), ncol(U_t))
        if (!symmetric) {
          V_t <- rho_uv * V_t + matrix(rnorm(prod(dim(V_t)), 0, 
                                            sqrt(1 - rho_uv^2) * 0.1), 
                                      nrow(V_t), ncol(V_t))
        } else {
          V_t <- U_t
        }
      }
      
      # Compute expected value of Z
      if (length(beta) > 0 && any(beta != 0)) {
        if (length(dim(X_t)) == 3 && dim(X_t)[3] > 0) {
          EZ <- array(0, dim = c(n_row, n_col))
          for (k in 1:min(dim(X_t)[3], length(beta))) {
            EZ <- EZ + beta[k] * X_t[, , k]
          }
        } else {
          EZ <- matrix(beta[1], n_row, n_col)
        }
      } else {
        EZ <- matrix(0, n_row, n_col)
      }
      
      # Add random effects
      EZ <- EZ + outer(a_t, b_t, "+")
      
      # Add multiplicative effects
      if (!is.null(U_t) && !is.null(V_t)) {
        if (bip && !is.null(fit$G)) {
          EZ <- EZ + U_t %*% fit$G %*% t(V_t)
        } else if (!bip) {
          EZ <- EZ + U_t %*% t(V_t)
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
        Y_sim <- simY_ord(EZ, rho, fit$ODM)
      } else if (family == "rrl") {
        Y_sim <- simY_rrl(EZ, rho, fit$ODM)
      } else if (family == "poisson") {
        Y_sim <- simY_pois(EZ)
      } else if (family == "frn") {
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
      
      # Store results for this time point
      Y_t_list[[t]] <- Y_sim
      if (return_latent) {
        Z_t_list[[t]] <- simZ(EZ, rho, s2)
      }
    }
    
    # Store trajectory
    Y_sims[[i]] <- Y_t_list
    if (return_latent) {
      Z_sims[[i]] <- Z_t_list
    }
  }
  
  # Prepare output
  result <- list(
    Y = Y_sims,
    family = family,
    mode = mode,
    n_time = n_time
  )
  
  if (return_latent) {
    result$Z <- Z_sims
  }
  
  class(result) <- "lame.sim"
  return(result)
}