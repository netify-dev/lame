# as_draws methods for ame / lame fits so they integrate cleanly with the
# stan ecosystem (posterior, bayesplot, tidybayes, loo).

#' Convert an AME / LAME fit to a posterior draws object
#'
#' Reshapes the BETA + VC posterior draws stored on an \code{ame} or
#' \code{lame} fit into a \code{posterior::draws_array} (or
#' \code{draws_df}-equivalent) suitable for downstream packages that
#' consume that representation. Multi-chain fits
#' (\code{ame(..., n_chains = N)}) are reshaped with chain as a separate
#' dimension; single-chain fits are returned as a single-chain
#' \code{draws_array}.
#'
#' The reshape uses \code{fit$chain_indicator} (set when \code{n_chains > 1})
#' to keep chain identity separate so within-chain Rhat / ESS estimates and
#' chain-coloured traceplots resolve correctly.
#'
#' @param x an \code{ame} or \code{lame} fit.
#' @param include character; which parameter blocks to include. Default
#'   \code{c("beta", "vc")}; can also include \code{"rho_uv"}, \code{"rho_ab"}
#'   for dynamic LAME fits.
#' @param ... ignored.
#' @return A \code{posterior::draws_array} (when \pkg{posterior} is
#'   available) or otherwise a plain 3-D array with named dimnames.
#' @export
as_draws.ame <- function(x, include = c("beta", "vc"), ...) {
	if (is.null(x$BETA) || dim(x$BETA)[1] == 0L) {
		cli::cli_abort("Fit has no stored posterior draws (BETA is empty).")
	}
	include <- intersect(include, c("beta", "vc", "rho_uv", "rho_ab"))
	# 3-d beta (dynamic_beta) -> flatten to [iter x (p*t)] with combined names
	beta_is_dyn <- length(dim(x$BETA)) == 3L
	beta_mat <- if (beta_is_dyn) {
		B <- x$BETA
		dn <- dimnames(B)
		nms_p <- if (!is.null(dn)) dn[[2]] else paste0("beta", seq_len(dim(B)[2]))
		nms_t <- if (!is.null(dn)) dn[[3]] else paste0("t", seq_len(dim(B)[3]))
		wide <- matrix(B, nrow = dim(B)[1], ncol = dim(B)[2] * dim(B)[3])
		colnames(wide) <- as.vector(outer(nms_p, nms_t,
			function(p, t) paste0(p, "[", t, "]")))
		wide
	} else {
		x$BETA
	}
	# infer n_chains / n_iter from chain_indicator if present
	ci <- x$chain_indicator
	n_iter_total <- dim(x$BETA)[1]
	if (is.null(ci) || length(ci) != n_iter_total) {
		n_chains <- 1L
		n_iter   <- n_iter_total
		idx_list <- list(seq_len(n_iter))
	} else {
		uc <- sort(unique(ci))
		n_chains <- length(uc)
		idx_list <- lapply(uc, function(c) which(ci == c))
		lens <- vapply(idx_list, length, integer(1))
		if (length(unique(lens)) != 1L) {
			cli::cli_warn("Unequal-length chains; truncating to the shortest.")
			n_iter <- min(lens)
			idx_list <- lapply(idx_list, head, n_iter)
		} else {
			n_iter <- lens[[1]]
		}
	}
	# assemble named draw matrices for each requested block
	blocks <- list()
	if ("beta" %in% include) blocks$beta <- beta_mat
	if ("vc"   %in% include && !is.null(x$VC)) blocks$vc <- x$VC
	if ("rho_uv" %in% include && !is.null(x$rho_uv))
		blocks$rho_uv <- matrix(x$rho_uv, ncol = 1, dimnames = list(NULL, "rho_uv"))
	if ("rho_ab" %in% include && !is.null(x$rho_ab))
		blocks$rho_ab <- matrix(x$rho_ab, ncol = 1, dimnames = list(NULL, "rho_ab"))
	# unify into one matrix (column-bind) -- assumes all blocks have nrow(beta) rows
	wide <- do.call(cbind, blocks)
	var_names <- colnames(wide)
	if (is.null(var_names)) var_names <- paste0("v", seq_len(ncol(wide)))
	# reshape into [iter, chain, var]
	arr <- array(NA_real_, dim = c(n_iter, n_chains, ncol(wide)),
	             dimnames = list(iteration = NULL,
	                             chain = paste0("c", seq_len(n_chains)),
	                             variable = var_names))
	for (c in seq_len(n_chains)) {
		arr[, c, ] <- wide[idx_list[[c]], , drop = FALSE]
	}
	if (requireNamespace("posterior", quietly = TRUE)) {
		return(posterior::as_draws_array(arr))
	}
	# graceful fallback: return the plain array; posterior is suggested, not required
	arr
}

#' @rdname as_draws.ame
#' @method as_draws lame
#' @export
as_draws.lame <- as_draws.ame

# as_draws for ame_als: route to bootstrap replicates when present,
# otherwise abort with a refit hint.
#' @rdname as_draws.ame
#' @method as_draws ame_als
#' @export
as_draws.ame_als <- function(x, ...) {
	if (is.null(x$bootstrap)) {
		cli::cli_abort(c(
			"{.fn posterior::as_draws} requires a draws-like object; an {.cls ame_als} point estimate has none.",
			"i" = "Refit with {.code ame_als(..., bootstrap = N)} or attach replicates via {.fn ame_als_bootstrap} to get a draws array."))
	}
	boot <- x$bootstrap
	# boot$coefs is the [r x p] matrix of per-replicate coefficient draws
	# (renamed from boot_coefs in the bootstrap assembly); treat as one
	# chain of r draws.
	mat <- boot$coefs
	if (is.null(mat) || !is.matrix(mat)) {
		cli::cli_abort("The bootstrap object on {.code fit$bootstrap} does not expose a per-draw coefficient matrix.")
	}
	# add iteration / chain attributes per posterior::as_draws_matrix contract
	if (requireNamespace("posterior", quietly = TRUE)) {
		return(posterior::as_draws_matrix(mat))
	}
	mat
}

#' Generic dispatcher for posterior::as_draws on lame fits
#'
#' Lightweight S3 generic so calls of the form \code{as_draws(fit)} dispatch
#' through R's S3 system even when the \pkg{posterior} package is not loaded.
#' When \pkg{posterior} \emph{is} loaded, its generic of the same name is
#' resolved first by R's namespace search; this fallback only fires for
#' bare-namespace use.
#'
#' @param x a fitted \code{ame} or \code{lame} object.
#' @param ... passed to the relevant method.
#' @return An object dispatched by the relevant method (typically a
#'   \code{posterior::draws_array} for an \code{ame} / \code{lame} fit).
#' @export
as_draws <- function(x, ...) UseMethod("as_draws")

# register the methods against posterior::as_draws when posterior is
# available so qualified posterior::as_draws(fit) dispatches to the
# methods defined here. called from the .onload hook in zzz.r.
.register_posterior_methods <- function() {
	if (requireNamespace("posterior", quietly = TRUE)) {
		registerS3method("as_draws", "ame",     as_draws.ame,
		                 envir = asNamespace("posterior"))
		registerS3method("as_draws", "lame",    as_draws.lame,
		                 envir = asNamespace("posterior"))
		registerS3method("as_draws", "ame_als", as_draws.ame_als,
		                 envir = asNamespace("posterior"))
	}
}
