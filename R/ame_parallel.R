#' Run AME model with multiple parallel chains
#'
#' \code{ame_parallel(Y, n_chains = k, ...)} is \code{ame(Y, chains = k, ...)}:
#' the same chains (seeds \code{seed + 1000 * (k - 1)}), the same drift check,
#' between-chain gate and extensions, and the same pooled fit. It is kept for the
#' older argument names; \code{combine_method = "list"} returns the fits of the
#' chains instead of the pooled fit.
#'
#' @param Y Network data matrix
#' @param n_chains Number of chains to run (default = 4)
#' @param cores Number of CPU cores to use: the chains run forked on Unix (one
#'   core per chain at most) and one after the other on Windows. Default is
#'   n_chains. Use 1 for sequential processing.
#' @param combine_method Method for combining chains: "pool" (default) or "list"
#' @param fitter Which fitter to use: \code{"auto"} (default; pick
#'   \code{\link{lame}} when \code{Y} is a list / 3-D array with >1 time slice,
#'   else \code{\link{ame}}), \code{"ame"} (force cross-sectional), or
#'   \code{"lame"} (force longitudinal -- required for \code{dynamic_*} flags).
#' @param ... Additional arguments passed to \code{ame()} or \code{lame()}.
#'
#' @return If combine_method = "pool": A single ame object with pooled chains
#'         If combine_method = "list": A list of ame objects, one per chain.
#'         A call whose chains are not settled after the extensions stops with
#'         an error of class \code{lame_chain_unsettled} (see \code{\link{ame}}).
#'
#' @seealso \code{\link{lame_multi}} for the unrelated multi-\emph{panel}
#'   wrapper that fits K \emph{distinct} networks with shared regression
#'   coefficients. \code{ame_parallel} / \code{lame_parallel} run K MCMC
#'   chains of the \emph{same} model (for R-hat / ESS diagnostics);
#'   \code{lame_multi} runs one chain across K panels with pooled beta.
#'
#' @examples
#' \donttest{
#' # Run 2 chains sequentially
#' data(YX_nrm)
#' fit_parallel <- ame_parallel(YX_nrm$Y, Xdyad = YX_nrm$X,
#'                              n_chains = 2, cores = 1,
#'                              nscan = 100, burn = 10, odens = 1,
#'                              verbose = FALSE)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
ame_parallel <- function(
	Y,
	n_chains = 4,
	cores = n_chains,
	combine_method = c("pool", "list"),
	fitter = c("auto", "ame", "lame"),
	...
) {

	combine_method <- match.arg(combine_method)
	fitter <- match.arg(fitter)

	# dispatch to lame() for longitudinal Y (a list or a 3-d array with more than one slice)
	if (fitter == "auto") {
		is_longitudinal <-
			(is.list(Y) && length(Y) > 1L) ||
			(is.array(Y) && length(dim(Y)) == 3L && dim(Y)[3] > 1L)
		fitter_name <- if (is_longitudinal) "lame" else "ame"
	} else fitter_name <- fitter

	# chains as in ame() / lame(): seeds seed + 1000 * (k - 1), the between-chain gate and the
	# single-chain extensions
	args <- c(list(Y = Y), list(...))
	args$chains <- as.integer(n_chains)
	args$cores <- cores
	if (identical(combine_method, "list")) {
		prior <- args$prior %||% list()
		attr(prior, "lame_chain_list") <- TRUE
		args$prior <- prior
	}
	cl <- as.call(c(list(as.name(fitter_name)), stats::setNames(lapply(names(args), as.name), names(args))))
	out <- eval(cl, list2env(args, parent = asNamespace("lame")))
	if (identical(combine_method, "list") && inherits(out, "ame") && !isTRUE(out$is_combined)) {
		attr(out$prior, "lame_chain_list") <- NULL
		out$chain_id <- 1L
		out <- list(out)
	}
	out
}

#' Run LAME (longitudinal AME) with multiple parallel chains
#'
#' Thin wrapper around \code{\link{ame_parallel}} that forces
#' \code{fitter = "lame"}. Use this for longitudinal data; in particular,
#' multi-chain \code{dynamic_beta} / \code{dynamic_uv} / \code{dynamic_ab}
#' fits are reached through here. The combined fit's \code{$BETA} is a 3-D
#' array when \code{dynamic_beta} is on, and
#' \code{posterior::summarise_draws(as_draws(fit), "rhat")}
#' gives R-hat across chains for every per-period coefficient.
#'
#' @param Y Longitudinal network: list of T relational matrices, or a 3-D
#'   array \code{[n, n, T]}.
#' @param n_chains Number of parallel chains (default 4).
#' @param cores CPU cores to use (default \code{n_chains}; 1 = sequential).
#' @param combine_method \code{"pool"} (default) or \code{"list"}.
#' @param ... Additional arguments forwarded to \code{\link{lame}}, including
#'   \code{dynamic_beta}, \code{dynamic_uv}, \code{dynamic_ab}, \code{family},
#'   \code{nscan}, \code{burn}, \code{odens}, etc.
#'
#' @return Same as \code{ame_parallel}: a combined \code{lame} fit (when
#'   \code{combine_method = "pool"}) or a list of fits.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit_pll <- lame_parallel(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), Xdyad = YX_bin_list$X,
#'                          family = "binary", n_chains = 2, cores = 1,
#'                          nscan = 50, burn = 10, odens = 5,
#'                          dynamic_beta = "dyad", verbose = FALSE)
#' dim(fit_pll$BETA)        # [iter * n_chains, p, T]
#' fit_pll$chain_indicator  # length iter * n_chains
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   posterior::summarise_draws(posterior::as_draws(fit_pll))
#' }
#' }
#'
#' @export
lame_parallel <- function(
	Y,
	n_chains = 4,
	cores = n_chains,
	combine_method = c("pool", "list"),
	...
) {
	combine_method <- match.arg(combine_method)
	ame_parallel(
		Y = Y,
		n_chains = n_chains,
		cores = cores,
		combine_method = combine_method,
		fitter = "lame",
		...
	)
}

#' Combine multiple AME chains
#'
#' @param chain_list List of ame fit objects from multiple chains
#' @param diagnostics Logical; whether to compute convergence diagnostics (default TRUE)
#'
#' @return A single ame object with combined chains. Conditional sampler
#'   health counts are summed when every chain supplies them. The original
#'   conditional counts and joint-move diagnostics are retained in
#'   \code{sampler_diagnostics_by_chain}; move diagnostics from a single
#'   chain are not presented as pooled values. If a legacy chain lacks the
#'   numerical/bound decomposition, the combined \code{sampler_health}
#'   is \code{NULL}; its available per-chain records remain accessible.
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # each fit is one chain (chains = 1); combine_ame_chains() pools them
#' fit1 <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, seed = 1, verbose = FALSE, chains = 1)
#' fit2 <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, seed = 2, verbose = FALSE, chains = 1)
#' pooled <- combine_ame_chains(list(fit1, fit2))
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
combine_ame_chains <- function(chain_list, diagnostics = TRUE) {

	n_chains <- length(chain_list)

	if (n_chains == 0) {
		cli::cli_abort("No chains provided")
	}

	if (n_chains == 1) {
		cli::cli_warn("Only one chain provided, returning as-is")
		return(chain_list[[1]])
	}

	if (any(vapply(chain_list, function(f) !is.null(f$log_lik_chunks), logical(1)))) {
		cli::cli_abort(c(
			"Combining chains with chunked likelihood storage is not implemented.",
			"i" = "Keep the chain objects separate, or refit with in-memory {.code save_log_lik = TRUE}."
		))
	}
	combined <- chain_list[[1]]
	combined <- .lame_pool_sampler_diagnostics(combined, chain_list)
	chain_sizes <- vapply(chain_list, function(f) dim(f$BETA)[1L], integer(1))
	chain_weights <- chain_sizes / sum(chain_sizes)
	# a chain that stored no draw (max_seconds during burn-in) has weight 0 and is left out of the averages
	weighted_mean <- function(xs) {
		use <- chain_weights > 0
		Reduce(`+`, Map(function(x, w) x * w, xs[use], chain_weights[use]))
	}

	# stack BETA and VC across chains; a dynamic_beta BETA [n_iter, p, T] is bound along the draws
	if (!is.null(combined$BETA)) {
		BETA_list <- lapply(chain_list, function(x) x$BETA)
		beta_is_dyn <- length(dim(BETA_list[[1]])) == 3L
		if (beta_is_dyn) {
			combined$BETA <- abind::abind(BETA_list, along = 1L)
		} else {
			combined$BETA <- do.call(rbind, BETA_list)
		}
		n_samples <- sapply(BETA_list, function(B) dim(B)[1])
		combined$chain_indicator <- rep(seq_len(n_chains), n_samples)
	}

	if (!is.null(combined$VC)) {
		VC_list <- lapply(chain_list, function(x) x$VC)
		combined$VC <- do.call(rbind, VC_list)
	}

	if (!is.null(combined$GOF)) {
		GOF_list <- lapply(chain_list, function(x) x$GOF)
		# one observed column or row, then the simulated draws of every chain
		if (is.list(GOF_list[[1]]) && !is.matrix(GOF_list[[1]])) {
			combined$GOF <- lapply(names(GOF_list[[1]]), function(nm) {
				do.call(
					cbind,
					c(
						list(GOF_list[[1]][[nm]][, 1L, drop = FALSE]),
						lapply(GOF_list, function(g) g[[nm]][, -1L, drop = FALSE])
					)
				)
			})
			names(combined$GOF) <- names(GOF_list[[1]])
		} else {
			combined$GOF <- do.call(
				rbind,
				c(
					list(GOF_list[[1]][1L, , drop = FALSE]),
					lapply(GOF_list, function(g) g[-1L, , drop = FALSE])
				)
			)
		}
	}

	# keep matching draw rows aligned across every saved posterior component
	for (nm in c(
		"U_samples", "V_samples", "G_samples", "U_draws", "V_draws",
		"a_samples", "b_samples", "a_draws", "b_draws",
		"a_level_samples", "b_level_samples"
	)) {
		xs <- lapply(chain_list, `[[`, nm)
		combined[[nm]] <- if (any(vapply(xs, is.null, logical(1)))) NULL else
			abind::abind(xs, along = length(dim(xs[[1L]])))
	}
	for (nm in c("rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "pi_snap")) {
		xs <- lapply(chain_list, `[[`, nm)
		combined[[nm]] <- if (any(vapply(xs, is.null, logical(1)))) NULL else
			unlist(xs, use.names = FALSE)
	}
	for (nm in c(
		"RHO", "RHO_BETA", "SIGMA_BETA", "RHO_ACTOR", "SIGMA_ACTOR",
		"snap_draws", "snap_draws_v", "THETA_ACTOR", "ALPHA", "log_lik"
	)) {
		xs <- lapply(chain_list, `[[`, nm)
		combined[[nm]] <- if (any(vapply(xs, is.null, logical(1)))) NULL else
		if (is.null(dim(xs[[1L]]))) unlist(xs, use.names = FALSE) else
			abind::abind(xs, along = 1L)
	}
	offsets <- c(0L, utils::head(cumsum(chain_sizes), -1L))
	combined$ab_draw_index <- if (all(vapply(chain_list, function(f) !is.null(f$ab_draw_index), logical(1))))
		unlist(Map(function(f, o) f$ab_draw_index + o, chain_list, offsets), use.names = FALSE) else NULL
	combined$uv_draw_index <- if (all(vapply(chain_list, function(f) !is.null(f$uv_draw_index), logical(1))))
		unlist(Map(function(f, o) f$uv_draw_index + o, chain_list, offsets), use.names = FALSE) else NULL
	combined$forecast_state <- .lame_pool_forecast_states(chain_list, offsets)
	for (nm in c("a_dynamic", "b_dynamic", "a_level", "b_level")) {
		xs <- lapply(chain_list, `[[`, nm)
		if (!any(vapply(xs, is.null, logical(1)))) combined[[nm]] <- weighted_mean(xs)
	}


	# U and V are identified only up to rotation, so the chains pool the rotation-invariant product UV'
	# rather than the factors themselves
	if (!is.null(combined$U) && !is.null(combined[["V"]]) &&
			ncol(combined$U) > 0L && ncol(combined[["V"]]) > 0L) {
		U_list <- lapply(chain_list, function(x) x$U)
		V_list <- lapply(chain_list, function(x) x[["V"]])
		# under dynamic_uv, U and V are 3-d and averaged element-wise to keep the dynamics
		u_is_dyn <- length(dim(U_list[[1]])) == 3L
		# bipartite fits pool U G V' and leave U, V and G per chain; exact [["G"]] avoids matching $GOF
		bip_G <- !is.null(chain_list[[1]][["G"]])
		G_list <- lapply(chain_list, function(x) x[["G"]])
		.slice2d <- function(A, t) matrix(A[, , t], nrow = dim(A)[1L], ncol = dim(A)[2L])
		if (u_is_dyn) {
			# positions are refactored from the pooled product below
			UVPM_list <- lapply(chain_list, `[[`, "UVPM")
			if (!any(vapply(UVPM_list, is.null, logical(1)))) {
				combined$UVPM <- weighted_mean(UVPM_list)
			} else {
				products <- Map(function(U, V, G) {
					out <- array(0, c(nrow(U), nrow(V), dim(U)[3L]))
					for (t in seq_len(dim(U)[3L])) {
						Ut <- .slice2d(U, t); Vt <- .slice2d(V, t)
						out[, , t] <- if (bip_G) Ut %*% G %*% t(Vt) else Ut %*% t(Vt)
					}
					out
				}, U_list, V_list, G_list)
				combined$UVPM <- weighted_mean(products)
			}

		} else if (bip_G) {
			# pool the stored per-chain UVPM; the last-draw product is a fallback for fits without one
			UVPM_list <- lapply(chain_list, function(x) x$UVPM)
			if (!any(vapply(UVPM_list, is.null, logical(1)))) {
				combined$UVPM <- weighted_mean(UVPM_list)
			} else {
				UVprod_list <- Map(function(U, V, G) U %*% G %*% t(V), U_list, V_list, G_list)
				combined$UVPM <- weighted_mean(UVprod_list)
			}
		} else {
			# pooled rotation-invariant product, from the stored full-rank UVPM when available, else the
			# rank-R
			# factor product; U and V are recovered from its SVD
			UVPM_list <- lapply(chain_list, function(x) x[["UVPM"]])
			UVprod <- if (!any(vapply(UVPM_list, is.null, logical(1)))) {
				weighted_mean(UVPM_list)
			} else {
				weighted_mean(Map(function(U, V) U %*% t(V), U_list, V_list))
			}
			combined$UVPM <- UVprod
			R_lat <- ncol(U_list[[1]])
			sv <- svd(UVprod, nu = R_lat, nv = R_lat)
			d_sqrt <- sqrt(pmax(sv$d[seq_len(R_lat)], 0))
			combined$U <- sweep(sv$u, 2, d_sqrt, "*")
			combined[["V"]] <- sweep(sv$v, 2, d_sqrt, "*")
			rownames(combined$U) <- rownames(U_list[[1]])
			rownames(combined[["V"]]) <- rownames(V_list[[1]])
		}
	} else if (!is.null(combined$UVPM)) {
		UVPM_list <- lapply(chain_list, function(x) x$UVPM)
		combined$UVPM <- weighted_mean(UVPM_list)
	}

	# symmetric fits store the multiplicative mean as ULUPM with V = NULL; pool it when every chain has
	# one and refresh U and L from its eigendecomposition so U L U' stays its rank-R truncation
	if (!is.null(combined$ULUPM)) {
		ULUPM_list <- lapply(chain_list, function(x) x[["ULUPM"]])
		if (!any(vapply(ULUPM_list, is.null, logical(1)))) {
			combined$ULUPM <- weighted_mean(ULUPM_list)
			# under dynamic_uv the positions are refactored from the pooled product below
			ulu_is_dyn <- length(dim(ULUPM_list[[1]])) == 3L
			if (!ulu_is_dyn && !is.null(combined$U)) {
				R_lat <- ncol(combined$U)
				eULU <- eigen(combined$ULUPM)
				eR <- which(rank(-abs(eULU$val), ties.method = "first") <= R_lat)
				combined$U <- eULU$vec[, eR, drop = FALSE]
				combined$L <- diag(eULU$val[eR], nrow = R_lat)
				rownames(combined$U) <- rownames(combined$ULUPM)
			}
		}
	}

	if (!is.null(combined$APM)) {
		APM_list <- lapply(chain_list, function(x) x$APM)
		combined$APM <- weighted_mean(APM_list)
	}

	if (!is.null(combined$BPM)) {
		BPM_list <- lapply(chain_list, function(x) x$BPM)
		combined$BPM <- weighted_mean(BPM_list)
	}

	# lame fits store YPM and EZ as per-period lists, ame fits as one matrix
	.avg_maybe_list <- function(items) {
		if (is.list(items[[1]]) && !is.matrix(items[[1]])) {
			T_ <- length(items[[1]])
			lapply(seq_len(T_), function(t) {
				weighted_mean(lapply(items, function(x) x[[t]]))
			})
		} else {
			weighted_mean(items)
		}
	}
	if (!is.null(combined$YPM)) {
		YPM_list <- lapply(chain_list, function(x) x$YPM)
		combined$YPM <- .avg_maybe_list(YPM_list)
		if (is.list(combined$YPM)) names(combined$YPM) <- names(YPM_list[[1]])
	}

	if (!is.null(combined$EZ)) {
		EZ_list <- lapply(chain_list, function(x) x$EZ)
		combined$EZ <- .avg_maybe_list(EZ_list)
		if (is.list(combined$EZ)) names(combined$EZ) <- names(EZ_list[[1]])
	}
	# Poisson fitted counts at the pooled posterior means of the predictor and residual variance
	if (identical(combined[["family"]], "poisson") && !is.null(combined$YPM)) {
		combined$YPM <- .poisson_fitted_from_fit(combined)
	}

	# convergence diagnostics
	if (diagnostics && !is.null(combined$BETA)) {
		combined$diagnostics <- compute_mcmc_diagnostics(chain_list)
	}

	# actor uncertainty combines within- and between-chain variation
	for (pair in list(
		c("APM", "APM_sd"),
		c("BPM", "BPM_sd"),
		c("a_dynamic", "a_dynamic_sd"),
		c("b_dynamic", "b_dynamic_sd")
	)) {
		means <- lapply(chain_list, `[[`, pair[1L]); sds <- lapply(chain_list, `[[`, pair[2L])
		combined[[pair[2L]]] <- if (any(vapply(c(means, sds), is.null, logical(1)))) NULL else
			sqrt(pmax(
				weighted_mean(Map(function(mu, sd) mu^2 + sd^2, means, sds)) -
					combined[[pair[1L]]]^2,
				0
			))
	}
	# coordinate summaries cannot inherit one chain's uncertainty after pooling
	combined$U_sd <- combined$V_sd <- NULL
	combined <- .lame_dynamic_positions(combined)
	if (!is.null(combined$ALPHA)) combined$alpha_post_mean <- c(0, colMeans(combined$ALPHA))
	if (!is.null(combined$THETA_ACTOR)) {
		combined$theta_actor_mean <- apply(combined$THETA_ACTOR, c(2, 3), mean)
		combined$theta_actor_sd <- apply(combined$THETA_ACTOR, c(2, 3), stats::sd)
	} else {
		means <- lapply(chain_list, `[[`, "theta_actor_mean")
		sds <- lapply(chain_list, `[[`, "theta_actor_sd")
		combined$theta_actor_mean <- if (any(vapply(means, is.null, logical(1)))) NULL else weighted_mean(means)
		combined$theta_actor_sd <- if (is.null(combined$theta_actor_mean) ||
				any(vapply(sds, is.null, logical(1))) || sum(chain_sizes) < 2L) NULL else {
					sums <- Map(
						function(mu, sd, n) (n - 1) * sd^2 +
							n * (mu - combined$theta_actor_mean)^2,
						means,
						sds,
						chain_sizes
					)
					sqrt(Reduce(`+`, sums) / (sum(chain_sizes) - 1))
				}
	}
	for (nm in c("theta_actor_mean", "theta_actor_sd")) {
		if (is.matrix(combined[[nm]])) dimnames(combined[[nm]]) <- dimnames(chain_list[[1L]]$theta_actor_mean)
	}
	for (pair in list(c("snap_draws", "snap_prob"), c("snap_draws_v", "snap_prob_v"))) {
		if (!is.null(combined[[pair[1L]]])) combined[[pair[2L]]] <- apply(combined[[pair[1L]]], c(2, 3), mean)
	}
	if (!is.null(combined$snap_draws_meta) && !is.null(combined$snap_draws))
		combined$snap_draws_meta$n_draws <- dim(combined$snap_draws)[1L]
	combined$n_chains <- n_chains
	combined$is_combined <- TRUE

	return(combined)
}

#' Compute MCMC convergence diagnostics for multiple chains
#'
#' @description
#' Computes R-hat and effective sample size (ESS) across two or
#' more independently-seeded chains, for both the regression coefficients
#' (\code{BETA}) and the variance components (\code{VC}). ESS is computed on
#' the pooled draws, so genuinely non-converged chains report a low ESS rather
#' than an inflated one.
#'
#' @param chain_list A list of independently seeded fits from the same model:
#'   either \code{ame}/\code{lame} objects carrying \code{BETA}, or
#'   \code{lame_robust} objects. Do not mix fit classes.
#'
#' @return A list with \code{rhat}, \code{rhat_classic}, \code{ess},
#'   \code{param_names} (over \code{BETA} and \code{VC}), \code{n_chains} and
#'   \code{n_samples}. For robust fits, additionally returns \code{ess_bulk},
#'   \code{ess_tail}, a full \code{diagnostics} data frame, and \code{pass}.
#'
#' @details
#' \code{rhat} is the split rank-normalised R-hat of Vehtari et al. (2021), the
#' larger of the bulk and tail values, as \code{posterior::rhat()} reports it.
#' It is the statistic of the between-chain gate of \code{\link{ame}} and
#' \code{\link{lame}}, computed on the first \code{min(n_samples)} draws of
#' every chain, so for a two-chain fit it equals \code{fit$rhat} for every
#' monitored coefficient and variance component, and it is the Rhat column of
#' \code{summary()}. It is \code{NA} for a quantity that is constant in every
#' chain and for chains with fewer than 8 draws. \code{rhat_classic} is the
#' Gelman-Rubin statistic of the whole chains, without splitting or rank
#' normalisation. The printed verdict compares \code{rhat} with the gate's
#' threshold, 1.05.
#'
#' For \code{lame_robust} fits, the \pkg{posterior} package is required.
#' Diagnostics include static coefficients, time intercepts, actor paths and
#' levels, estimated hyperparameters, every off-diagonal latent product, and
#' persistent pair precisions. Unused rank-zero factor hyperparameters and
#' fixed Gaussian precisions are excluded. Chains must have identical inputs
#' and priors, different seeds, and equal lengths of at least four retained
#' draws. Robust \code{pass} requires both \code{rhat} and \code{rhat_classic}
#' below 1.05 and finite positive ESS, bulk ESS and tail ESS for every
#' included parameter. The individual R-hat and effective sample size values
#' are reported for every included parameter.
#' The existing \code{ame}/\code{lame} diagnostic scope remains coefficients
#' and variance components.
#'
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # each fit is one chain (chains = 1); compute_mcmc_diagnostics() compares them
#' fit1 <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, seed = 1, verbose = FALSE, chains = 1)
#' fit2 <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1, seed = 2, verbose = FALSE, chains = 1)
#' diag_out <- compute_mcmc_diagnostics(list(fit1, fit2))
#' diag_out$rhat
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
compute_mcmc_diagnostics <- function(chain_list) {

	if (inherits(chain_list, "lame_robust"))
		cli::cli_abort("Supply a list of at least two independently seeded lame_robust fits.", call = NULL)
	if (is.list(chain_list) && any(vapply(chain_list, inherits, logical(1), "lame_robust")))
		return(.lame_robust_diagnostics(chain_list))

	# input contract: a list of per-chain fits, each carrying $BETA
	if (!is.list(chain_list) || inherits(chain_list, c("ame", "lame")) ||
			length(chain_list) == 0L ||
			!all(vapply(
				chain_list,
				function(x) is.list(x) && !is.null(x$BETA),
				logical(1)
			))) {
		cli::cli_abort(c(
			"{.arg chain_list} must be a list of fitted {.cls ame}/{.cls lame} objects, one per chain.",
			"i" = "Pass the per-chain fits (e.g. run {.fn ame} with {.arg n_chains} > 1, or supply a list of fits)."
		))
	}

	n_chains <- length(chain_list)
	if (n_chains < 2) {
		cli::cli_warn("Need at least 2 chains for convergence diagnostics")
		return(NULL)
	}
	# chains run with the same seed have zero between-chain variance, so R-hat would falsely signal
	# convergence
	b1 <- chain_list[[1]]$BETA
	chains_identical <- all(vapply(
		chain_list[-1],
		function(x) identical(x$BETA, b1),
		logical(1)
	))
	if (chains_identical) {
		cli::cli_warn(c(
			"All chains are identical; between-chain diagnostics are meaningless.",
			"i" = "Vary the sampler across chains by passing a different {.arg seed} to each fit.",
			"i" = "Note an external {.fn set.seed} does not affect the sampler; use the {.arg seed} argument."
		))
	}

	# R-hat and ESS for a block of per-chain matrices: rhat is the gate's split rank-normalised R-hat
	# on the first min(ns) draws (.lame_rhat_chains), rhat_classic the Gelman-Rubin statistic of the
	# whole chains
	diag_block <- function(mats) {
		mats <- lapply(mats, as.matrix)
		p <- ncol(mats[[1]])
		ns <- vapply(mats, nrow, integer(1))
		nm <- colnames(mats[[1]])
		if (is.null(nm)) nm <- paste0("p", seq_len(p))
		rh <- rc <- es <- stats::setNames(numeric(p), nm)
		for (j in seq_len(p)) {
			ch <- lapply(mats, function(x) x[, j])
			rh[j] <- .lame_rhat_chains(vapply(ch, function(x) x[seq_len(min(ns))], numeric(min(ns))))
			W <- mean(vapply(ch, stats::var, numeric(1)))
			B <- stats::var(vapply(ch, mean, numeric(1))) * min(ns)
			vp <- ((min(ns) - 1) * W + B) / min(ns)
			rc[j] <- if (is.finite(W) && W > 0) sqrt(vp / W) else NA_real_
			if (requireNamespace("coda", quietly = TRUE)) {
				# ESS on the pooled draws so between-chain disagreement counts, capped at the pooled draw count
				pooled <- unlist(ch, use.names = FALSE)
				es[j] <- min(
					as.numeric(coda::effectiveSize(coda::mcmc(pooled))),
					length(pooled)
				)
			}
		}
		list(rhat = rh, rhat_classic = rc, ess = es)
	}

	# a dynamic_beta BETA [n_stored, p, T] becomes [n_stored, p*T] so R-hat and ESS are per coefficient-period
	.flatten_beta <- function(B) {
		if (length(dim(B)) == 3L) {
			dn <- dimnames(B)
			M <- matrix(B, nrow = dim(B)[1])
			cn <- dn[[2]] %||% paste0("b", seq_len(dim(B)[2]))
			tn <- dn[[3]] %||% paste0("t", seq_len(dim(B)[3]))
			colnames(M) <- as.vector(outer(
				cn,
				tn,
				function(a, b)
					paste0(a, "[", b, "]")
			))
			M
		} else B
	}
	bd <- diag_block(lapply(chain_list, function(x) .flatten_beta(x$BETA)))
	rhat <- bd$rhat; rhat_classic <- bd$rhat_classic; ess <- bd$ess
	if (!is.null(chain_list[[1]]$VC)) {
		vd <- diag_block(lapply(chain_list, function(x) x$VC))
		rhat <- c(rhat, vd$rhat); rhat_classic <- c(rhat_classic, vd$rhat_classic); ess <- c(ess, vd$ess)
	}
	param_names <- names(rhat)
	# stored-draw count for both 2-D and 3-D BETA
	n_samples <- vapply(
		chain_list,
		function(x) dim(as.array(x$BETA))[1],
		integer(1)
	)

	diagnostics <- list(
		rhat = rhat,
		rhat_classic = rhat_classic,
		ess = ess,
		param_names = param_names,
		n_chains = n_chains,
		n_samples = n_samples
	)

	cli::cli_h3("MCMC convergence diagnostics")
	cli::cli_text("Number of chains: {.val {n_chains}}")
	cli::cli_text("Samples per chain: {.val {paste(n_samples, collapse = ', ')}}")
	cli::cli_text("Diagnostics cover regression coefficients and variance components.")

	rhat_max <- .lame_chain_rhat_max
	converged <- all(rhat < rhat_max, na.rm = TRUE)
	if (chains_identical) {
		# R-hat carries no information for identical chains (zero between-chain variance)
		cli::cli_text("Skipping the R-hat convergence verdict for identical chains.")
	} else if (converged) {
		cli::cli_alert_success("All parameters converged (R-hat < {rhat_max})")
	} else {
		problem_params <- which(rhat >= rhat_max)
		cli::cli_alert_warning("{.val {length(problem_params)}} parameter{?s} {?has/have} R-hat >= {rhat_max}")
		for (i in problem_params[1:min(5, length(problem_params))]) {
			cli::cli_text("  {param_names[i]}: R-hat = {.val {round(rhat[i], 3)}}")
		}
		if (length(problem_params) > 5) {
			cli::cli_text("  ... and {.val {length(problem_params) - 5}} more")
		}
	}

	return(diagnostics)
}

#' Helper function to handle NULL with default (base R has it from 4.4; the
#' package supports R 4.3)
#' @noRd
`%||%` <- function(x, y) {
	if (is.null(x)) y else x
}
