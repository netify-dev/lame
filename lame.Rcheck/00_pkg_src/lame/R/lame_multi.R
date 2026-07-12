# multi-panel lame() wrapper. fits each panel independently, then pools
# the per-panel beta posteriors by inverse-precision averaging. with
# k = 1 the wrapper falls through to lame() directly.

#' Multi-panel lame() with shared coefficients
#'
#' Fits K independent \code{lame()} models (one per panel) and pools
#' the per-panel beta posteriors into a precision-weighted shared
#' posterior. Returns a list with the per-panel fits, the pooled
#' beta posterior, and the panel-specific deviations.
#'
#' This is an R-level wrapper: it fits each panel with its own MCMC and
#' pools the results afterwards. The pooling is exact when the panels are
#' conditionally independent given beta, which is the standard
#' assumption.
#'
#' @param Y_list A list of length K, each element a list (or 3-D array)
#'   of T per-panel network observations.
#' @param Xdyad_list A list of length K, each element a list of T
#'   dyadic covariate arrays.
#' @param ... Arguments forwarded to \code{lame()} (e.g. family,
#'   R, mode, nscan, burn, odens, dynamic_beta, dynamic_beta_kind).
#' @param shared_dynamic_beta Logical; if \code{TRUE} (default), the
#'   shared beta posterior is computed (per-period when dynamic_beta
#'   is active).
#'
#' @return A list with
#' \itemize{
#'   \item \code{fits}: list of K per-panel \code{lame} fits.
#'   \item \code{beta_shared}: pooled posterior mean of beta
#'         (per-period when dynamic).
#'   \item \code{beta_deviations}: list of K panel-specific deviations
#'         from \code{beta_shared}.
#'   \item \code{K}: number of panels.
#' }
#' Class \code{"lame_multi"}.
#'
#' @seealso \code{\link{lame_parallel}} for the unrelated multi-\emph{chain}
#'   wrapper that runs K MCMC chains of the \emph{same} model (used for
#'   R-hat / ESS diagnostics and pooled effective sample size).
#'   \code{lame_multi} is for K \emph{distinct} panels with shared
#'   regression coefficients; \code{lame_parallel} is for K chains of
#'   one model.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit_multi <- lame_multi(
#'   Y_list = list(YX_bin_list$Y, YX_bin_list$Y),
#'   Xdyad_list = list(YX_bin_list$X, YX_bin_list$X),
#'   family = "binary", R = 0,
#'   nscan = 100, burn = 25, odens = 5, verbose = FALSE)
#' dim(fit_multi$beta_shared)
#' }
#'
#' @export
lame_multi <- function(Y_list, Xdyad_list, ...,
                       shared_dynamic_beta = TRUE) {
	if (!is.list(Y_list) || length(Y_list) < 1L) {
		cli::cli_abort("{.arg Y_list} must be a non-empty list of K panels.")
	}
	K <- length(Y_list)
	if (length(Xdyad_list) != K) {
		cli::cli_abort("{.arg Xdyad_list} must have length K = {.val {K}}.")
	}

	# fast path for k == 1 -- fit a single panel and return the same slot
	# schema as the k >= 2 path so downstream code iterating over k = 1, 2, ...
	# does not break on the k = 1 case.
	if (K == 1L) {
		cli::cli_inform(c(
			"i" = "{.fn lame_multi} with K = 1 falls through to {.fn lame}."))
		fit <- lame(Y_list[[1L]], Xdyad = Xdyad_list[[1L]], ...)
		B <- fit$BETA
		is_3d <- length(dim(B)) == 3L
		if (is_3d) {
			beta_shared     <- apply(B, c(2, 3), mean)
			beta_shared_var <- apply(B, c(2, 3), stats::var)
		} else {
			beta_shared     <- colMeans(B)
			beta_shared_var <- apply(B, 2L, stats::var)
		}
		beta_shared_se <- sqrt(pmax(beta_shared_var, 0))
		n_joint_draws  <- dim(B)[1L]
		BETA_joint     <- B
		# k=1 -> no pooling, deviations are zero by construction
		zero_dev <- if (is_3d) array(0, dim(beta_shared)) else
			rep(0, length(beta_shared))
		res <- list(
			fits             = list(fit),
			beta_shared      = beta_shared,
			beta_shared_var  = beta_shared_var,
			beta_shared_se   = beta_shared_se,
			BETA_joint       = BETA_joint,
			beta_deviations  = list(zero_dev),
			var_shrinkage    = 1,
			K                = 1L
		)
		class(res) <- c("lame_multi", "list")
		return(res)
	}

	# fit each panel independently
	fits <- vector("list", K)
	for (k in seq_len(K)) {
		fits[[k]] <- suppressWarnings(suppressMessages(
			lame(Y_list[[k]], Xdyad = Xdyad_list[[k]], ...)))
	}

	# information-form pooling. when panels are conditionally
	# independent given a shared beta, the joint posterior precision is
	# sum_k (1 / var_k) and the mean is the precision-weighted average.
	beta_dims <- vapply(fits, function(f) length(dim(f$BETA)), integer(1L))
	if (!all(beta_dims == beta_dims[1L])) {
		cli::cli_abort("All panel BETA arrays must have the same dimensionality.")
	}
	pool_out <- if (beta_dims[1L] == 3L) {
		# [iter x p x t] per panel
		p <- dim(fits[[1L]]$BETA)[2L]
		T_per <- dim(fits[[1L]]$BETA)[3L]
		mean_out <- matrix(0, p, T_per)
		w_sum    <- matrix(0, p, T_per)
		for (k in seq_len(K)) {
			means_k <- apply(fits[[k]]$BETA, c(2, 3), mean)
			vars_k  <- apply(fits[[k]]$BETA, c(2, 3), stats::var)
			prec_k  <- 1 / pmax(vars_k, 1e-8)
			mean_out <- mean_out + prec_k * means_k
			w_sum    <- w_sum + prec_k
		}
		list(mean = mean_out / pmax(w_sum, 1e-8),
		     var  = 1 / pmax(w_sum, 1e-8))
	} else {
		p <- ncol(fits[[1L]]$BETA)
		mean_out <- numeric(p); w_sum <- numeric(p)
		for (k in seq_len(K)) {
			means_k <- colMeans(fits[[k]]$BETA)
			vars_k  <- apply(fits[[k]]$BETA, 2L, stats::var)
			prec_k  <- 1 / pmax(vars_k, 1e-8)
			mean_out <- mean_out + prec_k * means_k
			w_sum    <- w_sum + prec_k
		}
		list(mean = mean_out / pmax(w_sum, 1e-8),
		     var  = 1 / pmax(w_sum, 1e-8))
	}
	beta_shared    <- pool_out$mean
	beta_shared_var <- pool_out$var
	beta_shared_se <- sqrt(beta_shared_var)
	# draws sampled from the pooled normal so downstream summarisation works
	n_joint_draws <- max(vapply(fits, function(f) dim(f$BETA)[1L], integer(1L)))
	if (beta_dims[1L] == 3L) {
		BETA_joint <- array(NA_real_, c(n_joint_draws, dim(beta_shared)[1L], dim(beta_shared)[2L]))
		for (d in seq_len(n_joint_draws)) {
			noise <- matrix(stats::rnorm(length(beta_shared)),
			                nrow(beta_shared), ncol(beta_shared)) *
				sqrt(beta_shared_var)
			BETA_joint[d, , ] <- beta_shared + noise
		}
	} else {
		BETA_joint <- matrix(NA_real_, n_joint_draws, length(beta_shared))
		for (d in seq_len(n_joint_draws)) {
			BETA_joint[d, ] <- beta_shared +
				stats::rnorm(length(beta_shared)) * sqrt(beta_shared_var)
		}
	}

	beta_deviations <- vector("list", K)
	for (k in seq_len(K)) {
		means_k <- if (beta_dims[1L] == 3L)
			apply(fits[[k]]$BETA, c(2, 3), mean) else colMeans(fits[[k]]$BETA)
		beta_deviations[[k]] <- means_k - beta_shared
	}

	# variance shrinkage ratio: under the conditional-independence
	# assumption it is approximately 1/k when panels share beta.
	var_shrinkage <- if (K >= 2L) {
		# average panel posterior variance
		avg_panel_var <- if (beta_dims[1L] == 3L) {
			mean(vapply(seq_len(K), function(k)
				mean(apply(fits[[k]]$BETA, c(2, 3), stats::var), na.rm = TRUE),
				numeric(1L)))
		} else {
			mean(vapply(seq_len(K), function(k)
				mean(apply(fits[[k]]$BETA, 2L, stats::var), na.rm = TRUE),
				numeric(1L)))
		}
		mean(beta_shared_var, na.rm = TRUE) / max(avg_panel_var, 1e-12)
	} else 1

	res <- list(
		fits             = fits,
		beta_shared      = beta_shared,
		beta_shared_var  = beta_shared_var,
		beta_shared_se   = beta_shared_se,
		BETA_joint       = BETA_joint,
		beta_deviations  = beta_deviations,
		var_shrinkage    = var_shrinkage,
		K                = K
	)
	class(res) <- c("lame_multi", "list")
	res
}

#' Print method for lame_multi
#' @param x A \code{lame_multi} object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   summary of the pooled multi-panel fit.
#' @export
#' @method print lame_multi
print.lame_multi <- function(x, digits = 3, ...) {
	cli::cli_h2("Multi-panel lame fit (K = {.val {x$K}})")
	cli::cli_text(cli::col_grey(paste0(
		"Information-form pooled posterior. Joint variance shrinks by ",
		"~ 1/K when panels share beta; observed ratio var_shared / avg_panel = ",
		round(x$var_shrinkage, 3), " (expected ~ ", round(1/x$K, 3), ").")))
	if (length(dim(x$beta_shared)) == 2L) {
		cli::cli_text("Shared beta posterior mean (p x T):")
	} else {
		cli::cli_text("Shared beta posterior mean (p):")
	}
	print(round(x$beta_shared, digits))
	cli::cli_text("")
	cli::cli_text("Shared beta posterior SE:")
	print(round(x$beta_shared_se, digits))
	invisible(x)
}
