#' Compute GOF statistics from saved posterior samples
#' 
#' @description
#' Computes goodness-of-fit statistics after model estimation using saved 
#' posterior samples. This avoids the computational overhead of GOF calculation
#' during MCMC sampling. Requires that the model was run with appropriate
#' posterior sampling options.
#' 
#' @param fit An ame model object that was run with posterior sampling enabled
#' @param Y Original data matrix (if not stored in fit)
#' @param custom_gof Optional custom GOF function(s) - same format as for ame()
#' @param nsim Number of posterior predictive simulations to generate (default 100).
#'   If NULL, uses all available posterior samples.
#' @param verbose Logical; print progress information
#' 
#' @return A matrix of GOF statistics with the same format as if gof=TRUE was
#'   used during model estimation. First row contains observed statistics,
#'   subsequent rows contain posterior predictive statistics.
#'   
#' @details
#' This function requires that the model was estimated with posterior sampling
#' of the parameters needed to generate posterior predictive datasets. 
#' Specifically, it needs:
#' - BETA: regression coefficients
#' - VC: variance components  
#' - For models with random effects: samples of a, b
#' - For models with latent factors: U_samples, V_samples
#' 
#' To enable posterior sampling during model estimation, use:
#' \code{posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE)}
#' 
#' Computing GOF post-hoc has several advantages:
#' - Faster MCMC sampling (no GOF overhead)
#' - Can experiment with different GOF statistics without re-running model
#' - Can control number of posterior predictive simulations independently
#' 
#' @examples
#' \dontrun{
#' # Run model with posterior sampling but no GOF
#' opts <- posterior_options(save_UV = TRUE, save_ab = TRUE)
#' fit <- ame(Y, R = 2, gof = FALSE, posterior_opts = opts)
#' 
#' # Compute GOF post-hoc
#' gof_result <- gof(fit)
#' 
#' # Add custom statistics
#' custom <- function(Y) c(density = mean(Y > 0, na.rm = TRUE))
#' gof_custom <- gof(fit, custom_gof = custom)
#' }
#' 
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
gof <- function(
  fit, Y = NULL, custom_gof = NULL, 
  nsim = 100, verbose = TRUE) {
  
  # Get Y from fit if not provided
  if (is.null(Y)) {
    if (!is.null(fit$Y)) {
      Y <- fit$Y
    } else {
      stop("Y must be provided either as argument or stored in fit object")
    }
  }
  
  # Check that we have necessary posterior samples
  if (is.null(fit$BETA) || is.null(fit$VC)) {
    stop("Model must have BETA and VC posterior samples. Re-run with appropriate posterior_opts.")
  }
  
  # Determine network type
  is_bipartite <- fit$mode == "bipartite"
  
  # Determine number of simulations
  n_mcmc <- nrow(fit$BETA)
  if (n_mcmc < 2) {
    stop("Insufficient MCMC samples. Model must have at least 2 saved iterations. Check burn, nscan, and odens parameters.")
  }
  if (is.null(nsim)) {
    nsim <- min(100, n_mcmc - 1)  # Default to 100 or available samples
  } else {
    nsim <- min(nsim, n_mcmc - 1)
  }
  
  if (verbose) {
    cli::cli_alert_info("Computing GOF statistics post-hoc")
    cli::cli_alert_info("Using {nsim} posterior predictive simulations")
  }
  
  # Compute observed GOF using model metadata
  network_mode <- if (is_bipartite) "bipartite" else "unipartite"
  obs_gof <- gof_stats(Y, mode = network_mode)
  
  n_base_stats <- length(obs_gof)
  
  # Handle custom GOF
  custom_gof_names <- NULL
  n_custom_stats <- 0
  
  if (!is.null(custom_gof)) {
    if (is.function(custom_gof)) {
      # Test custom function
      tryCatch({
        test_output <- custom_gof(Y)
        n_custom_stats <- length(test_output)
        custom_gof_names <- names(test_output)
        if (is.null(custom_gof_names)) {
          custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
        }
        # Add custom stats to observed
        obs_gof <- c(obs_gof, test_output)
      }, error = function(e) {
        warning("Custom GOF function failed. Skipping.")
        custom_gof <- NULL
      })
    } else if (is.list(custom_gof)) {
      n_custom_stats <- length(custom_gof)
      custom_gof_names <- names(custom_gof)
      if (is.null(custom_gof_names)) {
        custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
      }
      # Compute custom stats for observed
      custom_obs <- numeric(n_custom_stats)
      for (i in seq_along(custom_gof)) {
        if (is.function(custom_gof[[i]])) {
          tryCatch({
            custom_obs[i] <- custom_gof[[i]](Y)
          }, error = function(e) {
            custom_obs[i] <- NA
          })
        }
      }
      obs_gof <- c(obs_gof, custom_obs)
    }
  }
  
  # Initialize GOF matrix
  n_total_stats <- length(obs_gof)  # Use actual length of obs_gof
  GOF <- matrix(NA, nsim + 1, n_total_stats)
  GOF[1, ] <- obs_gof
  
  # Set names
  base_names <- names(gof_stats(Y, mode = network_mode))
  all_names <- c(base_names, custom_gof_names)
  colnames(GOF) <- all_names
  rownames(GOF) <- c("obs", rep("", nsim))
  
  # Progress bar
  if (verbose) {
    cli::cli_progress_bar("Computing posterior predictive GOF", total = nsim)
  }
  
  # Generate posterior predictive datasets using simulate.ame
  # Generate all simulations at once for efficiency
  sims <- simulate(fit, nsim = nsim, seed = 6886)
  
  # Compute GOF for each simulated network
  for (i in seq_len(nsim)) {
    Y_sim <- sims$Y[[i]]
    
    # Compute GOF for simulated data
    GOF[i + 1, 1:n_base_stats] <- gof_stats(Y_sim, mode = network_mode)
    
    # Compute custom GOF
    if (!is.null(custom_gof)) {
      if (is.function(custom_gof)) {
        tryCatch({
          custom_vals <- custom_gof(Y_sim)
          GOF[i + 1, (n_base_stats + 1):(n_base_stats + length(custom_vals))] <- custom_vals
        }, error = function(e) {})
      } else if (is.list(custom_gof)) {
        for (j in seq_along(custom_gof)) {
          if (is.function(custom_gof[[j]])) {
            tryCatch({
              GOF[i + 1, n_base_stats + j] <- custom_gof[[j]](Y_sim)
            }, error = function(e) {})
          }
        }
      }
    }
    
    if (verbose) {
      cli::cli_progress_update()
    }
  }
  
  if (verbose) {
    cli::cli_progress_done()
    cli::cli_alert_success("GOF computation complete")
  }
  
  return(GOF)
}

