#' Extract model coefficients from AME model
#'
#' Returns posterior means of regression coefficients from a fitted AME or
#' LAME model.
#'
#' For a static fit (the default, and any model with \code{dynamic_beta = FALSE}),
#' coefficients are returned as a named numeric vector computed as
#' \code{colMeans(fit$BETA)}.
#'
#' For a dynamic fit (\code{lame(..., dynamic_beta = ...)} where some
#' coefficient is time-varying), \code{fit$BETA} is a 3-dimensional array
#' \code{[n_stored, p, T]} and \code{coef.lame} returns a \code{[p, T]}
#' matrix of per-period posterior means. Rownames are the coefficient names;
#' colnames are the period labels (from \code{names(Y)} or \code{t1, t2, ...}).
#' Static coefficients in a dynamic fit are constant across the columns.
#'
#' For binary models, these are on the probit (latent) scale. Use
#' \code{\link{predict.ame}} with \code{type = "response"} to get predicted
#' probabilities.
#'
#' \strong{What \code{coef()} does not return.} The multiplicative latent
#' positions \eqn{U}, \eqn{V} are not part of the coefficient vector;
#' they live on \code{fit$U} and \code{fit$V} (or as 3-D arrays
#' \code{[n, R, T]} when \code{dynamic_uv} is on). The additive
#' sender / receiver effects \eqn{a, b} are on \code{fit$APM} and
#' \code{fit$BPM}. For a tidy frame of latent positions use
#' \code{\link{latent_positions}}; for sender / receiver lollipops use
#' \code{\link{ab_plot}}.
#'
#' @param object Fitted AME model (class \code{"ame"} or \code{"lame"}).
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector (static fit) or \code{p x T} matrix
#'   (\code{dynamic_beta} fit) of posterior mean coefficients.
#'
#' @seealso \code{\link{vcov.ame}} for the posterior covariance matrix,
#'   \code{\link{confint.ame}} for credible intervals,
#'   \code{\link{summary.ame}} for a full summary table
#'
#' @method coef ame
#' @export
coef.ame <- function(object, ...) {
	if(is.null(object$BETA)) {
		return(numeric(0))
	}
	# dynamic_beta storage is 3-d [iter x p x t]; return p x t matrix of
	# per-period posterior-mean coefficients with rownames = coef, colnames = period
	if (length(dim(object$BETA)) == 3L) {
		if (dim(object$BETA)[1] == 0L) return(numeric(0))
		out <- apply(object$BETA, c(2, 3), mean)
		# preserve dimnames from beta
		dn <- dimnames(object$BETA)
		if (!is.null(dn)) {
			rownames(out) <- dn[[2]]
			colnames(out) <- dn[[3]]
		}
		return(out)
	}
	if (nrow(object$BETA) == 0L) return(numeric(0))
	out <- colMeans(object$BETA)
	nms <- colnames(object$BETA)
	if(!is.null(nms)) names(out) <- nms
	out
}

#' @rdname coef.ame
#' @method coef lame
#' @export
coef.lame <- function(object, ...) {
	coef.ame(object, ...)
}

#' Posterior covariance of AME model coefficients
#'
#' Returns the posterior covariance matrix of regression coefficients.
#'
#' @param object fitted AME model (class "ame")
#' @param ... additional arguments (ignored)
#'
#' @return p x p covariance matrix of posterior BETA draws
#'
#' @method vcov ame
#' @export
vcov.ame <- function(object, ...) {
	if(is.null(object$BETA)) {
		return(matrix(NA, 0, 0))
	}
	# dynamic_beta storage: 3-d [iter x p x t] -- collapse to [iter x (p*t)]
	# so cov() returns a (p*t) x (p*t) covariance with combined coef[t] names
	if (length(dim(object$BETA)) == 3L) {
		B <- object$BETA
		if (dim(B)[1] < 2L) {
			pT <- dim(B)[2] * dim(B)[3]
			return(matrix(NA, pT, pT))
		}
		dn <- dimnames(B)
		nms_p <- if (!is.null(dn)) dn[[2]] else paste0("v", seq_len(dim(B)[2]))
		nms_t <- if (!is.null(dn)) dn[[3]] else paste0("t", seq_len(dim(B)[3]))
		wide_names <- as.vector(outer(nms_p, nms_t, function(p, t) paste0(p, "[", t, "]")))
		wide <- matrix(B, nrow = dim(B)[1], ncol = dim(B)[2] * dim(B)[3])
		out <- cov(wide)
		dimnames(out) <- list(wide_names, wide_names)
		return(out)
	}
	if (nrow(object$BETA) < 2L) {
		p <- ncol(object$BETA)
		return(matrix(NA, p, p))
	}
	out <- cov(object$BETA)
	nms <- colnames(object$BETA)
	if(!is.null(nms)) {
		rownames(out) <- colnames(out) <- nms
	}
	out
}

#' @rdname vcov.ame
#' @method vcov lame
#' @export
vcov.lame <- function(object, ...) {
	vcov.ame(object, ...)
}

#' Bayesian credible intervals for AME model parameters
#'
#' Returns posterior \strong{equal-tailed quantile-based} credible intervals
#' (not highest-posterior-density, HPD). Built directly from
#' \code{quantile(object$BETA, c(alpha/2, 1-alpha/2))} and
#' \code{quantile(object$VC, ...)}. These are Bayesian credible intervals, not
#' frequentist confidence intervals.
#'
#' @param object fitted AME / LAME model.
#' @param parm character vector of parameter names, or numeric indices. When
#'   \code{parm} is character, both BETA names (e.g. \code{"intercept"},
#'   \code{"x1_dyad"}) and variance-component names (\code{"va"}, \code{"vb"},
#'   \code{"cab"}, \code{"rho"}, \code{"ve"}) are accepted. \code{NULL}
#'   (default) returns intervals for all available parameters.
#' @param level credible level (default 0.95).
#' @param ... additional arguments (ignored).
#'
#' @return Matrix with one row per parameter and two columns
#'   (e.g. \code{"2.5\%"}, \code{"97.5\%"}).
#'
#' @section Note on interval type:
#' These are equal-tailed quantile intervals, not HPD. For an HPD interval use
#' e.g. \code{coda::HPDinterval} on the columns of \code{object$BETA} and
#' \code{object$VC} directly.
#'
#' @method confint ame
#' @export
confint.ame <- function(object, parm = NULL, level = 0.95, ...) {
	if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single number in (0, 1).")
	}
	if(is.null(object$BETA) || (length(dim(object$BETA)) >= 2L && dim(object$BETA)[1] == 0L)) {
		return(matrix(nrow = 0, ncol = 2,
									dimnames = list(NULL, c("lower", "upper"))))
	}

	alpha <- (1 - level) / 2
	probs <- c(alpha, 1 - alpha)
	pct <- paste0(round(probs * 100, 1), "%")

	# 3-d beta: each row of the ci table is one coef[t] combination.
	# loop time-outer, coef-inner so the row order matches what
	# vcov.ame() and as_draws.ame() produce (both built from
	# outer(nms_p, nms_t) in column-major / time-outer order), keeping
	# confint() rows aligned with vcov() diagonals.
	if (length(dim(object$BETA)) == 3L) {
		B <- object$BETA
		dn <- dimnames(B)
		nms_p <- if (!is.null(dn)) dn[[2]] else paste0("v", seq_len(dim(B)[2]))
		nms_t <- if (!is.null(dn)) dn[[3]] else paste0("t", seq_len(dim(B)[3]))
		ci_list <- list()
		for (t_ in seq_len(dim(B)[3])) {
			for (k in seq_len(dim(B)[2])) {
				q <- quantile(B[, k, t_], probs = probs, na.rm = TRUE)
				ci_list[[paste0(nms_p[k], "[", nms_t[t_], "]")]] <- q
			}
		}
		beta_ci <- do.call(rbind, ci_list)
		colnames(beta_ci) <- pct
	} else {
		beta_ci <- t(apply(object$BETA, 2, quantile, probs = probs, na.rm = TRUE))
		colnames(beta_ci) <- pct
	}

	if(!is.null(parm)) {
		vc_ci <- if (!is.null(object$VC) && ncol(object$VC) > 0) {
			m <- t(apply(object$VC, 2, quantile, probs = probs, na.rm = TRUE))
			colnames(m) <- pct
			m
		} else {
			matrix(numeric(0), 0, 2, dimnames = list(NULL, pct))
		}
		all_ci <- rbind(beta_ci, vc_ci)
		if(is.character(parm)) {
			bad <- setdiff(parm, rownames(all_ci))
			if (length(bad) > 0L) {
				cli::cli_abort(c(
					"Unknown parameter name{?s} in {.arg parm}: {.val {bad}}.",
					"i" = "Available: {.val {rownames(all_ci)}}."))
			}
			return(all_ci[parm, , drop = FALSE])
		} else if(is.numeric(parm)) {
			# positional indexing applies to beta rows only
			return(beta_ci[parm, , drop = FALSE])
		}
	}

	beta_ci
}

#' @rdname confint.ame
#' @method confint lame
#' @export
confint.lame <- function(object, parm = NULL, level = 0.95, ...) {
	confint.ame(object, parm = parm, level = level, ...)
}
