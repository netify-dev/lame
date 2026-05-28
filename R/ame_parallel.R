#' Run AME model with multiple parallel chains
#'
#' @param Y Network data matrix
#' @param n_chains Number of parallel chains to run (default = 4)
#' @param cores Number of CPU cores to use for parallel processing.
#'   Default is n_chains. Use 1 for sequential processing.
#' @param combine_method Method for combining chains: "pool" (default) or "list"
#' @param fitter Which fitter to use: \code{"auto"} (default; pick
#'   \code{\link{lame}} when \code{Y} is a list / 3-D array with >1 time slice,
#'   else \code{\link{ame}}), \code{"ame"} (force cross-sectional), or
#'   \code{"lame"} (force longitudinal -- required for \code{dynamic_*} flags).
#' @param ... Additional arguments passed to \code{ame()} or \code{lame()}.
#'
#' @return If combine_method = "pool": A single ame object with pooled chains
#'         If combine_method = "list": A list of ame objects, one per chain
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
	Y, n_chains = 4, cores = n_chains,
	combine_method = c("pool", "list"),
	fitter = c("auto", "ame", "lame"), ...
	){

	combine_method <- match.arg(combine_method)
	fitter <- match.arg(fitter)

	if(cores > 1 && !requireNamespace("parallel", quietly = TRUE)) {
		cli::cli_warn("Package 'parallel' not available. Running chains sequentially.")
		cores <- 1
	}
	# clamp cores to the machine's logical-core count; over-subscribing
	# is the most common silent foot-gun on shared workers
	if (cores > 1) {
		n_logical <- tryCatch(parallel::detectCores(logical = TRUE),
		                      error = function(e) NA_integer_)
		if (is.finite(n_logical) && cores > n_logical) {
			cli::cli_warn(c(
				"!" = "Requested {.val {cores}} cores but only {.val {n_logical}} are available; clamping.",
				"i" = "Set {.code cores = n_chains} explicitly to silence this warning."))
			cores <- n_logical
		}
	}

	ame_args <- list(...)

	####
	# dispatch to lame() for longitudinal data so that lame()-only args
	# (dynamic_beta etc.) reach the right fitter. auto-detect from Y:
	# a list or a 3-D array with > 1 slice is longitudinal.
	if (fitter == "auto") {
		is_longitudinal <-
			(is.list(Y) && length(Y) > 1L) ||
			(is.array(Y) && length(dim(Y)) == 3L && dim(Y)[3] > 1L)
		fitter_fn <- if (is_longitudinal) lame else ame
		fitter_name <- if (is_longitudinal) "lame" else "ame"
	} else if (fitter == "lame") {
		fitter_fn <- lame; fitter_name <- "lame"
	} else {
		fitter_fn <- ame; fitter_name <- "ame"
	}
	####

	####
	# generate seeds and define chain runner
	base_seed <- ame_args$seed %||% 6886
	chain_seeds <- base_seed + seq_len(n_chains) * 1000

	run_chain <- function(chain_id) {
		cli::cli_text("Starting chain {.val {chain_id}} ({.fn {fitter_name}})")

		ame_args$seed <- chain_seeds[chain_id]

		if(!is.null(ame_args$out_file)) {
			ame_args$out_file <- gsub("\\.RData$",
																paste0("_chain", chain_id, ".RData"),
																ame_args$out_file)
		}

		fit <- do.call(fitter_fn, c(list(Y = Y), ame_args))
		fit$chain_id <- chain_id

		cli::cli_text("Completed chain {.val {chain_id}}")

		return(fit)
	}
	####

	####
	# run chains
	if(cores > 1) {
		cli::cli_h2("Running {.val {n_chains}} chains in parallel on {.val {cores}} cores")

		cl <- parallel::makeCluster(cores)
		on.exit(parallel::stopCluster(cl))

		parallel::clusterEvalQ(cl, {
			library(lame)
		})

		parallel::clusterExport(cl,
			c("Y", "ame_args", "run_chain", "fitter_fn", "fitter_name"),
			envir = environment())

		chain_results <- parallel::parLapply(cl, seq_len(n_chains), run_chain)

	} else {
		cli::cli_h2("Running {.val {n_chains}} chains sequentially")

		chain_results <- lapply(seq_len(n_chains), run_chain)
	}
	####

	####
	# combine results
	if(combine_method == "pool") {
		cli::cli_text("Combining chains...")
		combined_fit <- combine_ame_chains(chain_results)
		return(combined_fit)
	} else {
		return(chain_results)
	}
	####
}

#' Run LAME (longitudinal AME) with multiple parallel chains
#'
#' Thin wrapper around \code{\link{ame_parallel}} that forces
#' \code{fitter = "lame"}. Use this for longitudinal data; in particular,
#' multi-chain \code{dynamic_beta} / \code{dynamic_uv} / \code{dynamic_ab}
#' fits are reached through here. The combined fit's \code{$BETA} is a 3-D
#' array when \code{dynamic_beta} is on, and \code{posterior::rhat(as_draws(fit))}
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
#' fit_pll <- lame_parallel(YX_bin_list$Y, Xdyad = YX_bin_list$X,
#'                          family = "binary", n_chains = 2, cores = 1,
#'                          nscan = 100, burn = 20, odens = 5,
#'                          dynamic_beta = "dyad", verbose = FALSE)
#' dim(fit_pll$BETA)        # [iter * n_chains, p, T]
#' fit_pll$chain_indicator  # length iter * n_chains
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   posterior::summarise_draws(posterior::as_draws(fit_pll))
#' }
#' }
#'
#' @export
lame_parallel <- function(Y, n_chains = 4, cores = n_chains,
                          combine_method = c("pool", "list"), ...) {
	combine_method <- match.arg(combine_method)
	ame_parallel(Y = Y, n_chains = n_chains, cores = cores,
	             combine_method = combine_method,
	             fitter = "lame", ...)
}

#' Combine multiple AME chains
#'
#' @param chain_list List of ame fit objects from multiple chains
#' @param diagnostics Logical; whether to compute convergence diagnostics (default TRUE)
#'
#' @return A single ame object with combined chains
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
combine_ame_chains <- function(chain_list, diagnostics = TRUE) {

	n_chains <- length(chain_list)

	if(n_chains == 0) {
		cli::cli_abort("No chains provided")
	}

	if(n_chains == 1) {
		cli::cli_warn("Only one chain provided, returning as-is")
		return(chain_list[[1]])
	}

	combined <- chain_list[[1]]

	####
	# combine BETA and VC matrices. BETA is 3-D when dynamic_beta is on
	# ([n_iter, p, T]); use abind along the first dim for that case so
	# the result has shape [sum(n_iter), p, T] with coefficient and period
	# dimnames preserved from the first chain. plain rbind handles 2-D.
	if(!is.null(combined$BETA)) {
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

	if(!is.null(combined$VC)) {
		VC_list <- lapply(chain_list, function(x) x$VC)
		combined$VC <- do.call(rbind, VC_list)
	}

	if(!is.null(combined$GOF)) {
		GOF_list <- lapply(chain_list, function(x) x$GOF)
		# GOF on lame fits is a NAMED LIST of [stat, T] matrices, not a 3-D
		# array. rbind() on the list-of-lists path fails; concatenate per-stat
		# across chains instead. ame fits store GOF as a 2-D matrix [iter,
		# 5/3 stats] and rbind works there.
		if (is.list(GOF_list[[1]]) && !is.matrix(GOF_list[[1]])) {
			combined$GOF <- lapply(names(GOF_list[[1]]), function(nm) {
				do.call(rbind, lapply(GOF_list, function(g) g[[nm]]))
			})
			names(combined$GOF) <- names(GOF_list[[1]])
		} else {
			combined$GOF <- do.call(rbind, GOF_list)
		}
	}
	####

	####
	# combine posterior samples if saved
	if(!is.null(combined$U_samples)) {
		U_samples_list <- lapply(chain_list, function(x) x$U_samples)
		combined$U_samples <- abind::abind(U_samples_list, along = 3)
	}

	if(!is.null(combined$V_samples)) {
		V_samples_list <- lapply(chain_list, function(x) x$V_samples)
		combined$V_samples <- abind::abind(V_samples_list, along = 3)
	}

	if(!is.null(combined$a_samples)) {
		a_samples_list <- lapply(chain_list, function(x) x$a_samples)
		combined$a_samples <- do.call(cbind, a_samples_list)
	}

	if(!is.null(combined$b_samples)) {
		b_samples_list <- lapply(chain_list, function(x) x$b_samples)
		combined$b_samples <- do.call(cbind, b_samples_list)
	}
	####

	####
	# average posterior means across chains. The multiplicative factors U, V
	# are identified only up to a rotation (UR, VR^{-T}) -- element-wise
	# averaging of rotation-equivalent matrices is meaningless. Instead, pool
	# the rotation-INVARIANT product U%*%V' across chains, then recover U/V
	# from its SVD. If the chains' U/V truly correspond up to small rotations
	# this is exact; if they sit in different basins, the singular spectrum
	# of the pooled product still has meaning while the individual U, V do
	# not -- the user is better off with the rotation-stable summary.
	if(!is.null(combined$U) && !is.null(combined$V) &&
	   ncol(combined$U) > 0L && ncol(combined$V) > 0L) {
		U_list <- lapply(chain_list, function(x) x$U)
		V_list <- lapply(chain_list, function(x) x$V)
		# under dynamic_uv U/V are 3-D; preserve the time dimension by
		# averaging element-wise instead of an SVD of the time-collapsed
		# product (which would lose the dynamics).
		u_is_dyn <- length(dim(U_list[[1]])) == 3L
		if (u_is_dyn) {
			combined$U <- Reduce("+", U_list) / n_chains
			combined$V <- Reduce("+", V_list) / n_chains
			# UVPM analogue under dynamic_uv -- per-time UV products averaged
			combined$UVPM <- Reduce("+",
				Map(function(U, V) {
					out <- array(0, dim = c(nrow(U), nrow(V), dim(U)[3]))
					for (t in seq_len(dim(U)[3])) out[,,t] <- U[,,t] %*% t(V[,,t])
					out
				}, U_list, V_list)) / n_chains
		} else {
			# pooled posterior-mean product (rotation-invariant)
			UVprod_list <- Map(function(U, V) U %*% t(V), U_list, V_list)
			UVprod <- Reduce("+", UVprod_list) / n_chains
			combined$UVPM <- UVprod
			R_lat <- ncol(U_list[[1]])
			sv <- svd(UVprod, nu = R_lat, nv = R_lat)
			d_sqrt <- sqrt(pmax(sv$d[seq_len(R_lat)], 0))
			combined$U <- sweep(sv$u, 2, d_sqrt, "*")
			combined$V <- sweep(sv$v, 2, d_sqrt, "*")
			rownames(combined$U) <- rownames(U_list[[1]])
			rownames(combined$V) <- rownames(V_list[[1]])
		}
	} else if (!is.null(combined$UVPM)) {
		UVPM_list <- lapply(chain_list, function(x) x$UVPM)
		combined$UVPM <- Reduce("+", UVPM_list) / n_chains
	}

	if(!is.null(combined$APM)) {
		APM_list <- lapply(chain_list, function(x) x$APM)
		combined$APM <- Reduce("+", APM_list) / n_chains
	}

	if(!is.null(combined$BPM)) {
		BPM_list <- lapply(chain_list, function(x) x$BPM)
		combined$BPM <- Reduce("+", BPM_list) / n_chains
	}

	# lame fits store YPM/EZ as a per-period list of matrices; ame fits
	# store them as a single matrix. average element-wise either way.
	.avg_maybe_list <- function(items, n) {
		if (is.list(items[[1]]) && !is.matrix(items[[1]])) {
			T_ <- length(items[[1]])
			lapply(seq_len(T_), function(t) {
				Reduce("+", lapply(items, function(x) x[[t]])) / n
			})
		} else {
			Reduce("+", items) / n
		}
	}
	if(!is.null(combined$YPM)) {
		YPM_list <- lapply(chain_list, function(x) x$YPM)
		combined$YPM <- .avg_maybe_list(YPM_list, n_chains)
		# preserve names on the per-period list for lame fits
		if (is.list(combined$YPM)) names(combined$YPM) <- names(YPM_list[[1]])
	}

	if(!is.null(combined$EZ)) {
		EZ_list <- lapply(chain_list, function(x) x$EZ)
		combined$EZ <- .avg_maybe_list(EZ_list, n_chains)
		if (is.list(combined$EZ)) names(combined$EZ) <- names(EZ_list[[1]])
	}
	####

	####
	# convergence diagnostics
	if(diagnostics && !is.null(combined$BETA)) {
		combined$diagnostics <- compute_mcmc_diagnostics(chain_list)
	}

	combined$n_chains <- n_chains
	combined$is_combined <- TRUE
	####

	return(combined)
}

#' Compute MCMC convergence diagnostics for multiple chains
#'
#' @description
#' Computes Gelman-Rubin R-hat and effective sample size (ESS) across two or
#' more independently-seeded chains, for both the regression coefficients
#' (\code{BETA}) and the variance components (\code{VC}). ESS is computed on
#' the pooled draws, so genuinely non-converged chains report a low ESS rather
#' than an inflated one.
#'
#' @param chain_list a list of fitted \code{ame}/\code{lame} objects, one per
#'   chain (each must carry a \code{BETA} matrix).
#'
#' @return A list with \code{rhat}, \code{ess}, \code{param_names} (over
#'   \code{BETA} and \code{VC}), \code{n_chains} and \code{n_samples}.
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#'
#' @export
compute_mcmc_diagnostics <- function(chain_list) {

	# input contract: a list of per-chain fits, each carrying $BETA
	if (!is.list(chain_list) || inherits(chain_list, c("ame", "lame")) ||
	    length(chain_list) == 0L ||
	    !all(vapply(chain_list,
	                function(x) is.list(x) && !is.null(x$BETA), logical(1)))) {
		cli::cli_abort(c(
			"{.arg chain_list} must be a list of fitted {.cls ame}/{.cls lame} objects, one per chain.",
			"i" = "Pass the per-chain fits (e.g. run {.fn ame} with {.arg n_chains} > 1, or supply a list of fits)."))
	}

	n_chains <- length(chain_list)
	if(n_chains < 2) {
		cli::cli_warn("Need at least 2 chains for convergence diagnostics")
		return(NULL)
	}

	# R-hat + ESS for one block of parameters (a list of per-chain matrices)
	diag_block <- function(mats) {
		mats <- lapply(mats, as.matrix)
		p  <- ncol(mats[[1]])
		ns <- vapply(mats, nrow, integer(1))
		nm <- colnames(mats[[1]])
		if (is.null(nm)) nm <- paste0("p", seq_len(p))
		rh <- es <- stats::setNames(numeric(p), nm)
		for (j in seq_len(p)) {
			ch <- lapply(mats, function(x) x[, j])
			W  <- mean(vapply(ch, stats::var, numeric(1)))
			B  <- stats::var(vapply(ch, mean, numeric(1))) * min(ns)
			vp <- ((min(ns) - 1) * W + B) / min(ns)
			rh[j] <- if (is.finite(W) && W > 0) sqrt(vp / W) else NA_real_
			if (requireNamespace("coda", quietly = TRUE)) {
				# ESS on the POOLED draws: summing per-chain ESS (the old
				# behaviour) ignores between-chain disagreement and reports a
				# hugely inflated ESS for non-converged chains
				es[j] <- as.numeric(coda::effectiveSize(
					coda::mcmc(unlist(ch, use.names = FALSE))))
			}
		}
		list(rhat = rh, ess = es)
	}

	bd <- diag_block(lapply(chain_list, function(x) x$BETA))
	rhat <- bd$rhat; ess <- bd$ess
	# variance components, when present
	if (!is.null(chain_list[[1]]$VC)) {
		vd <- diag_block(lapply(chain_list, function(x) x$VC))
		rhat <- c(rhat, vd$rhat); ess <- c(ess, vd$ess)
	}
	param_names <- names(rhat)
	n_samples <- vapply(chain_list, function(x) nrow(as.matrix(x$BETA)),
	                    integer(1))

	diagnostics <- list(
		rhat = rhat, ess = ess, param_names = param_names,
		n_chains = n_chains, n_samples = n_samples)

	cli::cli_h3("MCMC Convergence Diagnostics")
	cli::cli_text("Number of chains: {.val {n_chains}}")
	cli::cli_text("Samples per chain: {.val {paste(n_samples, collapse = ', ')}}")
	cli::cli_text("Diagnostics cover regression coefficients and variance components.")

	converged <- all(rhat < 1.1, na.rm = TRUE)
	if(converged) {
		cli::cli_alert_success("All parameters converged (R-hat < 1.1)")
	} else {
		problem_params <- which(rhat >= 1.1)
		cli::cli_alert_warning("{.val {length(problem_params)}} parameter{?s} have R-hat >= 1.1")
		for(i in problem_params[1:min(5, length(problem_params))]) {
			cli::cli_text("  {param_names[i]}: R-hat = {.val {round(rhat[i], 3)}}")
		}
		if(length(problem_params) > 5) {
			cli::cli_text("  ... and {.val {length(problem_params) - 5}} more")
		}
	}

	return(diagnostics)
}

#' Helper function to handle NULL with default
#' @noRd
`%||%` <- function(x, y) {
	if(is.null(x)) y else x
}
