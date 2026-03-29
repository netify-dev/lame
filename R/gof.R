#' Compute GOF statistics from saved posterior samples
#' 
#' @description
#' Computes goodness-of-fit statistics after model estimation by generating
#' posterior predictive networks from the saved MCMC samples. This is useful
#' when the model was fitted with \code{gof = FALSE} to speed up MCMC sampling,
#' or when you want to evaluate custom GOF statistics without re-running the model.
#' 
#' @param fit An ame model object that was run with posterior sampling enabled
#' @param Y Original data. For ame objects, an n x n matrix (or nA x nB for
#'   bipartite). For lame objects, a list of matrices (one per time period).
#'   If NULL, extracted from the fit object.
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
#' \donttest{
#' # Run model without GOF during fitting
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2, gof = FALSE,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE)
#'
#' # Compute GOF post-hoc
#' gof_result <- gof(fit)
#' }
#' 
#' @seealso \code{\link{ame}}, \code{\link{lame}}, \code{\link{gof_plot}}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
gof <- function(
	fit, Y = NULL, custom_gof = NULL, 
	nsim = 100, verbose = TRUE) {
	
	if (is.null(Y)) {
		if (!is.null(fit$Y)) {
			Y <- fit$Y
		} else {
			stop("Y must be provided either as argument or stored in fit object")
		}
	}

	# lame stores Y as a 3D array; average across time for GOF baseline
	is_lame <- inherits(fit, "lame")
	if (is.array(Y) && length(dim(Y)) == 3) {
		Y_list <- lapply(seq_len(dim(Y)[3]), function(t) Y[, , t])
	} else if (is.list(Y)) {
		Y_list <- Y
	} else {
		Y_list <- list(Y)
	}

	if (is.null(fit$BETA) || is.null(fit$VC)) {
		stop("Model must have BETA and VC posterior samples. Re-run with appropriate posterior_opts.")
	}

	is_bipartite <- fit$mode == "bipartite"
	n_mcmc <- nrow(fit$BETA)
	if (n_mcmc < 2) {
		stop("Insufficient MCMC samples. Model must have at least 2 saved iterations. Check burn, nscan, and odens parameters.")
	}
	if (is.null(nsim)) {
		nsim <- min(100, n_mcmc - 1)
	} else {
		nsim <- min(nsim, n_mcmc - 1)
	}

	if (verbose) {
		cli::cli_alert_info("Computing GOF statistics post-hoc")
		cli::cli_alert_info("Using {nsim} posterior predictive simulations")
	}

	network_mode <- if (is_bipartite) "bipartite" else "unipartite"

	# compute observed GOF per time point, then average
	obs_gof_list <- lapply(Y_list, function(y) gof_stats(y, mode = network_mode))
	obs_gof <- Reduce("+", obs_gof_list) / length(obs_gof_list)
	
	n_base_stats <- length(obs_gof)
	
	custom_gof_names <- NULL
	n_custom_stats <- 0
	
	if (!is.null(custom_gof)) {
		if (is.function(custom_gof)) {
			tryCatch({
				test_output <- custom_gof(Y)
				n_custom_stats <- length(test_output)
				custom_gof_names <- names(test_output)
				if (is.null(custom_gof_names)) {
					custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
				}
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
	
	n_total_stats <- length(obs_gof)
	GOF <- matrix(NA, nsim + 1, n_total_stats)
	GOF[1, ] <- obs_gof
	
	base_names <- names(gof_stats(Y_list[[1]], mode = network_mode))
	all_names <- c(base_names, custom_gof_names)
	colnames(GOF) <- all_names
	rownames(GOF) <- c("obs", rep("", nsim))

	if (verbose) {
		cli::cli_progress_bar("Computing posterior predictive GOF", total = nsim)
	}

	sims <- simulate(fit, nsim = nsim, seed = 6886)

	# helper: compute average GOF across time periods for a simulated dataset
	avg_gof <- function(y_sim) {
		if (is.list(y_sim)) {
			gof_list <- lapply(y_sim, function(y) gof_stats(y, mode = network_mode))
			Reduce("+", gof_list) / length(gof_list)
		} else {
			gof_stats(y_sim, mode = network_mode)
		}
	}

	for (i in seq_len(nsim)) {
		Y_sim <- sims$Y[[i]]

		GOF[i + 1, 1:n_base_stats] <- avg_gof(Y_sim)

		# custom GOF
		if (!is.null(custom_gof)) {
			# for lame sims, pass first time period to custom functions
			Y_for_custom <- if (is.list(Y_sim)) Y_sim[[1]] else Y_sim
			if (is.function(custom_gof)) {
				tryCatch({
					custom_vals <- custom_gof(Y_for_custom)
					GOF[i + 1, (n_base_stats + 1):(n_base_stats + length(custom_vals))] <- custom_vals
				}, error = function(e) {})
			} else if (is.list(custom_gof)) {
				for (j in seq_along(custom_gof)) {
					if (is.function(custom_gof[[j]])) {
						tryCatch({
							GOF[i + 1, n_base_stats + j] <- custom_gof[[j]](Y_for_custom)
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

