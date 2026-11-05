# split a panel into period matrices, keeping singleton dimensions
.lame_gof_periods <- function(Y) {
	if (length(dim(Y)) == 3L) {
		z <- lapply(
			seq_len(dim(Y)[3L]),
			function(t) matrix(
				Y[, , t],
				nrow = dim(Y)[1L],
				ncol = dim(Y)[2L],
				dimnames = dimnames(Y)[1:2]
			)
		)
		names(z) <- dimnames(Y)[[3L]]
	} else z <- if (is.list(Y) && !is.data.frame(Y)) Y else list(Y)
	if (!length(z) || !all(vapply(z, function(y) is.matrix(y) && (is.numeric(y) || is.logical(y)), logical(1))))
		cli::cli_abort(
			"GOF data must be a numeric matrix, a three-dimensional array, or a list of numeric matrices.", class = "lame_gof_alignment"
		)
	z
}

# align a predictive replicate to the observed actors and cells
.lame_gof_align_replicate <- function(y, observed) {
	for (axis in 1:2) {
		want <- dimnames(observed)[[axis]]; have <- dimnames(y)[[axis]]
		if (!is.null(want) && !is.null(have)) {
			if (anyDuplicated(want) || anyDuplicated(have) || anyNA(want) || anyNA(have))
				cli::cli_abort("GOF actor labels must be unique and nonmissing.", class = "lame_gof_alignment")
			idx <- match(want, have)
			if (anyNA(idx)) cli::cli_abort(
				"Observed GOF actors are absent from the fitted predictive grid.", class = "lame_gof_alignment"
			)
			y <- if (axis == 1L) y[idx, , drop = FALSE] else y[, idx, drop = FALSE]
		} else if (dim(y)[axis] != dim(observed)[axis]) {
			cli::cli_abort(
				"GOF dimensions differ and actor labels are unavailable for alignment.", class = "lame_gof_alignment"
			)
		}
	}
	y[!is.finite(observed)] <- NA_real_
	y
}

# average a custom statistic over periods, checking it is consistent
.lame_gof_custom_average <- function(fun, periods, scalar = FALSE) {
	z <- lapply(periods, fun)
	valid <- vapply(
		z,
		function(x) is.numeric(x) && is.null(dim(x)) &&
			length(x) == length(z[[1L]]) && identical(names(x), names(z[[1L]])),
		logical(1)
	)
	if (!length(z[[1L]]) || !all(valid) || (scalar && length(z[[1L]]) != 1L))
		cli::cli_abort(
			"Custom GOF must return a numeric vector with consistent length and names for every period.", call = NULL
		)
	Reduce(`+`, z) / length(z)
}

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
#' @param custom_gof Optional custom GOF function(s), each accepting one
#'   period's matrix. A function can return a named numeric vector; a list
#'   contains scalar-valued functions. Observed and replicated statistics
#'   are both averaged across the same fitted periods.
#' @param nsim Number of posterior predictive simulations to generate (default 100).
#'   Posterior rows are sampled with replacement, so this may exceed the
#'   number of stored rows. If NULL, uses the number of post-initial rows.
#' @param verbose Logical; print progress information
#'
#' @return A \code{(nsim + 1) x n_stat} matrix of GOF statistics: the first
#'   row holds the observed statistics, the remaining rows the posterior
#'   predictive ones. For an \code{ame} fit this is the format of
#'   \code{fit$GOF}; a \code{lame} fit stores its own goodness of fit as a
#'   per-statistic list of period-by-draw matrices instead, so the two are
#'   not interchangeable. Attributes \code{replication_method} and
#'   \code{draw_index} record joint replication and the sampled MCMC rows.
#'
#' @details
#' Gaussian bipartite checks use weighted \code{four.cycle.dep}; a nonzero
#' count would be fixed by the missing-cell pattern. Each replicate is
#' evaluated on the same observed cells as its corresponding data slice,
#' aligning actor labels before applying the mask. Named data periods are
#' matched to the fit; unnamed periods must be in fitted order. All custom
#' functions receive individual matrices for both observed and replicated
#' data, and their values are averaged over periods with equal weight.
#'
#' For the ordinal family the posterior predictive draws are quantile-matched
#' to the observed category frequencies (the cutpoints are induced from the
#' data), so the simulated marginal distribution always equals the observed
#' one: the GOF statistics and \code{YPM} cannot show a marginal misfit for
#' ordinal outcomes and only assess the dependence structure (degree spread,
#' dyadic and triadic dependence).
#'
#' This function requires joint posterior replication and refuses fits that
#' lack aligned saved draws; it never substitutes conditional mean trajectories.
#' This function requires that the model was estimated with posterior sampling
#' of the parameters needed to generate posterior predictive datasets.
#' Specifically, it needs:
#' - BETA: regression coefficients
#' - VC: variance components
#' - For models with random effects: samples of a, b
#' - For models with latent factors: U_samples, V_samples
#'
#' To enable posterior sampling during model estimation, use:
#' \code{posterior_opts = posterior_options(save_UV = TRUE, save_ab = TRUE,
#' thin_UV = 1, thin_ab = 1)}. At rank zero, omit \code{save_UV}.
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
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2, gof = FALSE,
#'            nscan = 100, burn = 10, odens = 1, verbose = FALSE,
#'            posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE,
#'                                               thin_ab = 1, thin_UV = 1), chains = 1)
#'
#' # Compute GOF post-hoc
#' gof_result <- gof(fit)
#' }
#'
#' @seealso \code{\link{ame}}, \code{\link{lame}}, \code{\link{gof_plot}}
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
gof <- function(
	fit,
	Y = NULL,
	custom_gof = NULL,
	nsim = 100,
	verbose = TRUE
) {

	if (is.null(Y)) {
		if (!is.null(fit$Y)) {
			Y <- fit$Y
		} else {
			cli::cli_abort("Y must be provided either as argument or stored in fit object")
		}
	}

	if (!is.null(custom_gof) && !is.function(custom_gof) &&
			!(is.list(custom_gof) && all(vapply(custom_gof, is.function, logical(1)))))
		cli::cli_abort("custom_gof must be a function or a list of functions.", class = "lame_gof_custom")
	Y_list <- .lame_gof_periods(Y)
	if (!is.null(fit$Y)) {
		reference <- .lame_gof_periods(fit$Y)
		if (length(reference) != length(Y_list))
			cli::cli_abort("GOF data must contain the fitted number of periods.", class = "lame_gof_alignment")
		if (!is.null(names(reference)) && !is.null(names(Y_list))) {
			idx <- match(names(reference), names(Y_list))
			if (anyNA(idx) || anyDuplicated(names(reference)) || anyDuplicated(names(Y_list)))
				cli::cli_abort(
					"GOF period labels must match the fitted periods uniquely.", class = "lame_gof_alignment"
				)
			Y_list <- Y_list[idx]
		}
	}

	if (is.null(fit$BETA) || is.null(fit$VC)) {
		cli::cli_abort("Model must have BETA and VC posterior samples. Re-run with appropriate posterior_opts.")
	}

	is_bipartite <- identical(fit$mode, "bipartite")
	n_mcmc <- nrow(fit$BETA)
	if (n_mcmc < 2) {
		cli::cli_abort("Insufficient MCMC samples. Model must have at least 2 saved iterations. Check burn, nscan, and odens parameters.")
	}
	if (is.null(nsim)) nsim <- n_mcmc - 1L
	if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1L || nsim > .Machine$integer.max || nsim != as.integer(nsim))
		cli::cli_abort("nsim must be a positive integer.")

	if (verbose) {
		cli::cli_alert_info("Computing GOF statistics post-hoc")
		cli::cli_alert_info("Using {nsim} posterior predictive simulations")
	}

	network_mode <- if (is_bipartite) "bipartite" else "unipartite"

	# observed GOF averaged over periods
	obs_gof_list <- lapply(Y_list, function(y) gof_stats(y, mode = network_mode, family = fit$family))
	obs_gof <- Reduce("+", obs_gof_list) / length(obs_gof_list)

	n_base_stats <- length(obs_gof)

	custom_gof_names <- NULL
	n_custom_stats <- 0
	# collects custom_gof failures for a single warning at the end
	.cg_tracker <- .new_custom_gof_err_tracker()

	if (!is.null(custom_gof)) {
		if (is.function(custom_gof)) {
			tryCatch({
				test_output <- .lame_gof_custom_average(custom_gof, Y_list)
				n_custom_stats <- length(test_output)
				custom_gof_names <- names(test_output)
				if (is.null(custom_gof_names)) {
					custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
				}
				obs_gof <- c(obs_gof, test_output)
			}, error = function(e) {
				cli::cli_warn("Custom GOF function failed. Skipping.")
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
						custom_obs[i] <- .lame_gof_custom_average(custom_gof[[i]], Y_list, scalar = TRUE)
					}, error = function(e) {
						.record_custom_gof_err(.cg_tracker, e)
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

	base_names <- names(gof_stats(Y_list[[1]], mode = network_mode, family = fit$family))
	all_names <- c(base_names, custom_gof_names)
	colnames(GOF) <- all_names
	rownames(GOF) <- c("obs", rep("", nsim))

	if (verbose) {
		cli::cli_progress_bar("Computing posterior predictive GOF", total = nsim)
	}

	sims <- simulate(fit, nsim = nsim, seed = 6886, type = "posterior")

	# GOF of a simulated dataset averaged over periods
	avg_gof <- function(y_sim) {
		if (is.list(y_sim)) {
			gof_list <- lapply(y_sim, function(y) gof_stats(y, mode = network_mode, family = fit$family))
			Reduce("+", gof_list) / length(gof_list)
		} else {
			gof_stats(y_sim, mode = network_mode, family = fit$family)
		}
	}

	for (i in seq_len(nsim)) {
		Y_sim <- .lame_gof_periods(sims$Y[[i]])
		if (length(Y_sim) != length(Y_list))
			cli::cli_abort("Predictive and observed GOF periods differ.", class = "lame_gof_alignment")
		Y_sim <- lapply(
			seq_along(Y_sim),
			function(t)
				.lame_gof_align_replicate(Y_sim[[t]], Y_list[[t]])
		)

		GOF[i + 1, 1:n_base_stats] <- avg_gof(Y_sim)

		if (!is.null(custom_gof)) {
			if (is.function(custom_gof)) {
				tryCatch({
					custom_vals <- .lame_gof_custom_average(custom_gof, Y_sim)
					if (length(custom_vals) != n_custom_stats ||
							(!is.null(names(custom_vals)) && !identical(names(custom_vals), custom_gof_names)))
						cli::cli_abort(
							"Custom GOF output names and length must match the observed statistic.", call = NULL
						)
					GOF[i + 1, (n_base_stats + 1):(n_base_stats + length(custom_vals))] <- custom_vals
				}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
			} else if (is.list(custom_gof)) {
				for (j in seq_along(custom_gof)) {
					if (is.function(custom_gof[[j]])) {
						tryCatch({
							GOF[i + 1, n_base_stats + j] <- .lame_gof_custom_average(
								custom_gof[[j]],
								Y_sim,
								scalar = TRUE
							)
						}, error = function(e) .record_custom_gof_err(.cg_tracker, e))
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

	.maybe_warn_custom_gof(.cg_tracker, where = "GOF simulation")

	attr(GOF, "replication_method") <- sims$method
	attr(GOF, "draw_index") <- sims$draw_index
	return(GOF)
}
