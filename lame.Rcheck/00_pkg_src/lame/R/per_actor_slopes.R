# post-mcmc per-actor time-varying slopes via ridge-penalised least
# squares on the residual after the main mcmc components. returns a
# point estimate, not a posterior summary.

#' Post-MCMC per-actor time-varying slopes
#'
#' Computes a smoothed per-actor time-varying slope coefficient on a
#' nodal covariate slice. For \code{kind = "row"}, each row-actor
#' \eqn{i} gets a length-\eqn{T} slope path \eqn{\beta_{i,t}} on
#' \code{Xrow[i, covariate_idx, t]}, fit by ridge-penalised least
#' squares on the residual after the main MCMC linear predictor.
#'
#' @param fit A fitted \code{lame} object.
#' @param kind \code{"row"} (default; per-row-actor slopes) or
#'   \code{"col"} (per-column-actor slopes).
#' @param covariate_idx Integer column index into the chosen nodal
#'   covariate (\code{Xrow} or \code{Xcol}). Default \code{1L}.
#' @param lambda Non-negative smoothing parameter (first-difference
#'   penalty across periods). Default \code{1}.
#'
#' @return A list with \code{slopes} (an \eqn{n_{actors} \times T}
#'   matrix), \code{kind}, \code{covariate_idx}, \code{lambda},
#'   and \code{label}. Class \code{"per_actor_slopes"}.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             nscan = 100, burn = 25, odens = 5, verbose = FALSE)
#' # post-MCMC per-row-actor slopes on the first dyadic covariate
#' pas <- per_actor_slopes(fit, kind = "row", lambda = 1)
#' dim(pas$slopes)
#' }
#'
#' @export
per_actor_slopes <- function(fit, kind = c("row", "col"),
                              covariate_idx = 1L, lambda = 1) {
	kind <- match.arg(kind)
	if (lambda < 0) cli::cli_abort("{.arg lambda} must be non-negative.")
	if (!inherits(fit, "lame")) {
		cli::cli_abort("{.arg fit} must be a fitted {.cls lame} object.")
	}

	# extract per-period residuals from the in-sample ez
	Y_list <- if (is.list(fit$Y)) fit$Y else if (length(dim(fit$Y)) == 3L)
		lapply(seq_len(dim(fit$Y)[3L]), function(t) fit$Y[, , t]) else
		list(fit$Y)
	EZ_list <- fit$EZ
	if (is.null(EZ_list)) {
		cli::cli_abort("Fit has no {.code $EZ}; cannot decompose residuals.")
	}
	T_per <- length(EZ_list)
	# aggregate the dyadic covariate slice to a per-actor value per period.
	# xlist is [n_a x n_b x p] per period; we pick covariate slice
	# `covariate_idx` and take the per-actor mean across the *other* index.
	X_dyad <- fit$Xlist
	if (is.null(X_dyad) || length(X_dyad) < T_per) {
		cli::cli_abort(c(
			"No {.code Xlist} found on the fit.",
			"i" = "{.fn per_actor_slopes} requires a dyadic covariate to slope on."))
	}

	n_actors <- if (kind == "row") {
		if (!is.null(fit$nA)) fit$nA else nrow(EZ_list[[1]])
	} else {
		if (!is.null(fit$nB)) fit$nB else ncol(EZ_list[[1]])
	}

	# residual per period (z - ez) summed across the other index
	r_per_actor <- matrix(0, n_actors, T_per)
	x_per_actor <- matrix(0, n_actors, T_per)
	for (t in seq_len(T_per)) {
		# observed cells
		Yt <- Y_list[[min(length(Y_list), t)]]
		EZt <- EZ_list[[t]]
		resid_t <- Yt - EZt
		# pick the dyadic covariate slice
		Xt <- X_dyad[[t]]
		if (length(dim(Xt)) < 3L) {
			Xt <- array(Xt, c(dim(Xt), 1L))
		}
		idx_pick <- min(covariate_idx, dim(Xt)[3L])
		X_slice <- Xt[, , idx_pick]
		if (kind == "row") {
			# row-actor mean of residual and covariate slice across columns
			r_per_actor[, t] <- rowMeans(resid_t, na.rm = TRUE)
			x_per_actor[, t] <- rowMeans(X_slice,  na.rm = TRUE)
		} else {
			r_per_actor[, t] <- colMeans(resid_t, na.rm = TRUE)
			x_per_actor[, t] <- colMeans(X_slice,  na.rm = TRUE)
		}
	}

	# ridge / first-difference penalised ls per actor
	slopes <- matrix(0, n_actors, T_per)
	for (i in seq_len(n_actors)) {
		x_i <- x_per_actor[i, ]
		r_i <- r_per_actor[i, ]
		ok <- is.finite(x_i) & is.finite(r_i)
		if (sum(ok) < 1L) next
		# build a [t x t] precision: data prec + lambda * d'd (1st-difference)
		if (T_per >= 2L) {
			D <- diff(diag(T_per))
			K <- crossprod(D)
		} else {
			K <- matrix(0, T_per, T_per)
		}
		# data precision: diag(x_i^2)
		# observed time positions only
		prec <- diag(ifelse(ok, x_i^2, 0), T_per) +
			lambda * K + diag(1e-8, T_per)
		rhs  <- ifelse(ok, x_i * r_i, 0)
		slopes[i, ] <- tryCatch(
			as.numeric(solve(prec, rhs)),
			error = function(e) rep(NA_real_, T_per))
	}

	res <- list(
		slopes        = slopes,
		kind          = kind,
		covariate_idx = covariate_idx,
		lambda        = lambda,
		label         = "post-MCMC per-actor time-varying slopes (penalised LS estimate, NOT a posterior summary)"
	)
	class(res) <- c("per_actor_slopes", "list")
	res
}

#' Print method for per_actor_slopes
#' @param x A \code{per_actor_slopes} result.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return \code{x}, invisibly. Called for its side effect of printing a
#'   summary of the per-actor time-varying slopes.
#' @export
#' @method print per_actor_slopes
print.per_actor_slopes <- function(x, digits = 3, ...) {
	cli::cli_h2("Per-actor time-varying slopes ({.val {x$kind}})")
	cli::cli_ul(c(
		"actors x periods: {.val {nrow(x$slopes)}} x {.val {ncol(x$slopes)}}",
		"covariate_idx: {.val {x$covariate_idx}}",
		"lambda: {.val {x$lambda}}",
		"label: {.emph {x$label}}"
	))
	# show a small head
	head_show <- min(6L, nrow(x$slopes))
	print(round(x$slopes[seq_len(head_show), , drop = FALSE], digits))
	if (nrow(x$slopes) > head_show) cli::cli_text("... ({nrow(x$slopes) - head_show} more rows)")
	invisible(x)
}
