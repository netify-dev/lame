# ame_als_sandwich.r
#
# analytic (sandwich) covariance for the dyadic regression coefficients of a
# fast ame fit. this is a fast, closed-form alternative to the bootstrap for
# (intercept, dyadic-covariate) standard errors.
#
# it is a *conditional* sandwich: the additive effects a, b and the
# multiplicative term o are treated as fixed offsets, so it captures sampling
# variability of (mu, beta) given the estimated network structure but not the
# uncertainty in a, b, u, v themselves. it is therefore anti-conservative --
# the bootstrap (ame_als_bootstrap) remains the recommended uncertainty
# tool, and the only one for a, b, u, v and the variance components.

#' Sandwich covariance for the regression coefficients of a fast AME fit
#'
#' @description
#' Returns a heteroskedasticity-robust (optionally dyad-clustered) sandwich
#' covariance matrix for the intercept and \emph{dyadic}-covariate coefficients
#' of an \code{\link{ame_als}} fit. This is a fast analytic alternative to
#' the bootstrap for those coefficients.
#'
#' @details
#' The estimate is the conditional sandwich \eqn{B^{-} M B^{-}} with bread
#' \eqn{B = D'WD} (\eqn{D} the observed intercept + dyadic-covariate design,
#' \eqn{W} the fit's observation weights) and meat \eqn{M} the
#' heteroskedasticity-robust (\code{cluster = "none"}, an HC0 meat) or
#' dyad-clustered (\code{cluster = "dyad"}, the default) outer product of the
#' weighted score contributions \eqn{w_\ell e_\ell d_\ell}. For a normal or
#' transform fit the weights are unit, so this reduces to the ordinary
#' \eqn{D'D} sandwich; for an IRLS fit it uses the final IRLS weights, matching
#' the estimating equation the fit actually solved. Dyad clustering pools the
#' score across \eqn{(i,j)}, \eqn{(j,i)} and time, so it reflects dyadic
#' dependence (reciprocity, repeated observation).
#'
#' It is \strong{conditional}: the additive effects \code{a}, \code{b} and the
#' multiplicative term are held fixed, so it omits their estimation uncertainty
#' and is anti-conservative. Node-covariate, additive and multiplicative
#' standard errors are not returned -- use \code{\link{ame_als_bootstrap}}
#' for those and for fully-propagated inference.
#'
#' @param object an \code{ame_als} fit.
#' @param cluster \code{"dyad"} (default) for a dyad-clustered robust meat, or
#'   \code{"none"} for an HC0 (heteroskedasticity-only) meat. Ignored when the
#'   fit carries a bootstrap, since the bootstrap covariance is returned.
#' @param ... ignored.
#' @return A covariance matrix with matching row/column names. When the fit
#'   carries a \code{$bootstrap}, this is the bootstrap covariance over all
#'   estimated coefficients, matching \code{coef()}. Otherwise it is the
#'   conditional sandwich, covering \code{c(intercept, dyadic coefficients)}
#'   only.
#' @method vcov ame_als
#' @seealso \code{\link{ame_als_bootstrap}} for bootstrap uncertainty
#'   covering all parameters.
#' @export
vcov.ame_als <- function(object, cluster = c("dyad", "none"), ...) {
	cluster_given <- !missing(cluster)
	cluster <- match.arg(cluster)
	# the sandwich is tied to the surrogate gaussian working likelihood.
	# for the non-gaussian families (binary-transform, poisson-transform) it
	# is anti-conservative because it conditions on the transformed working
	# response and ignores the estimation uncertainty in the additive
	# effects. the bootstrap is the recommended uncertainty for those fits.
	# when a bootstrap is attached, return its covariance so that vcov(),
	# summary() and confint() all report the same uncertainty on one fit. the
	# replicate matrix reproduces $bootstrap$se exactly under cov().
	boot_coefs <- object$bootstrap[["coefs", exact = TRUE]]
	if (is.matrix(boot_coefs)) {
		# failed replicates are stored as NA rows; the se/ci machinery
		# already drops them, so the covariance must too
		boot_coefs <- boot_coefs[stats::complete.cases(boot_coefs), , drop = FALSE]
	}
	if (is.matrix(boot_coefs) && nrow(boot_coefs) > 1L) {
		if (cluster_given) {
			cli::cli_inform(c(
				"i" = "{.arg cluster} applies to the sandwich only; this fit carries a bootstrap, so its covariance is returned."))
		}
		V <- stats::cov(boot_coefs)
		dimnames(V) <- list(colnames(boot_coefs), colnames(boot_coefs))
		return(V)
	}
	if (object$family %in% c("binary", "poisson")) {
		cli::cli_warn(c(
			"!" = "{.fn vcov.ame_als} for {.val {object$family}} returns the conditional sandwich on the surrogate working likelihood -- anti-conservative.",
			"i" = "Refit with {.code ame_als(..., bootstrap = 200)} for fully propagated uncertainty.",
			"i" = "Or call {.code confint(fit)} after attaching a bootstrap with {.fn ame_als_bootstrap}."))
	}
	X  <- object$X
	p  <- dim(X)[3]
	pa <- p + 1L
	Tt <- object$n_time
	nr <- object$dims$n_row
	bip <- identical(object$mode, "bipartite")
	resid <- object$residuals
	if (is.null(resid)) {
		cli::cli_abort("This fit carries no residuals; refit with the current {.pkg lame}.")
	}

	# observation weights: unit for a normal/transform fit, the final irls
	# weights otherwise (so the sandwich matches the estimating equation solved)
	W_arr <- object$obs_weights

	# stack observed-cell design rows, residuals, weights and dyad cluster ids
	Dl <- vector("list", Tt); el <- vector("list", Tt)
	wl <- vector("list", Tt); cl <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		r <- resid[[t]]
		obs <- which(is.finite(r))
		if (length(obs) == 0L) next
		Dt <- matrix(1, length(obs), pa)
		if (p > 0) for (k in seq_len(p)) Dt[, k + 1L] <- X[, , k, t][obs]
		ij <- arrayInd(obs, dim(r))
		id <- if (bip) {
			paste(ij[, 1], ij[, 2], sep = "-")
		} else {                             # unordered dyad: pool (i,j),(j,i)
			paste(pmin(ij[, 1], ij[, 2]), pmax(ij[, 1], ij[, 2]), sep = "-")
		}
		wt <- if (is.null(W_arr)) rep(1, length(obs)) else W_arr[, , t][obs]
		wt[!is.finite(wt)] <- 0
		Dl[[t]] <- Dt; el[[t]] <- r[obs]; wl[[t]] <- wt; cl[[t]] <- id
	}
	D <- do.call(rbind, Dl)
	e <- unlist(el)
	w <- unlist(wl)
	g <- unlist(cl)
	if (is.null(D) || nrow(D) <= pa) {
		cli::cli_abort("Too few observed cells for a sandwich covariance.")
	}

	B  <- crossprod(D, w * D)                # bread d'wd
	score <- D * (w * e)                     # weighted score contributions
	M  <- if (cluster == "dyad") {
		Sg <- rowsum(score, g)               # per-dyad score sums
		crossprod(Sg)
	} else {
		crossprod(score)                     # hc0 meat
	}
	# b^- via symmetric eigendecomposition (robust to a rank-deficient design)
	eg  <- eigen((B + t(B)) / 2, symmetric = TRUE)
	thr <- max(eg$values, 0) * sqrt(.Machine$double.eps)
	keep <- eg$values > thr
	Qk  <- eg$vectors[, keep, drop = FALSE]
	Binv <- Qk %*% (t(Qk) / eg$values[keep])
	V <- Binv %*% M %*% Binv
	V <- (V + t(V)) / 2

	nm <- c("intercept",
	        if (p > 0) .ae_default(object$x_names[seq_len(p)],
	                               paste0("dyad", seq_len(p))))
	dimnames(V) <- list(nm, nm)

	# the conditional sandwich is defined for the dyadic-cell regression only;
	# node-covariate coefficients are identified from the additive effects and
	# are not covered. name them explicitly rather than dropping them silently.
	omitted <- setdiff(names(object$coefficients), nm)
	if (length(omitted) > 0) {
		cli::cli_inform(c(
			"i" = "{.fn vcov} covers the intercept and dyadic coefficients only.",
			"i" = "Node-covariate coefficients not covered: {.val {omitted}}.",
			"i" = "Use {.fn ame_als_bootstrap} for those."))
	}
	V
}

#' Confidence intervals for a fast AME fit
#'
#' @description
#' Wald confidence intervals for the intercept and dyadic-covariate
#' coefficients of an \code{\link{ame_als}} fit, built from the conditional
#' sandwich covariance (\code{\link{vcov.ame_als}}).
#'
#' @details
#' These intervals are a fast convenience. They are \strong{conditional} (the
#' additive and multiplicative effects are held fixed) and therefore
#' anti-conservative, and they cover only the regression coefficients the
#' sandwich covariance is defined for -- not the node-covariate, additive or
#' multiplicative parameters. For fully-propagated intervals on all parameters,
#' use \code{\link{ame_als_bootstrap}} and \code{\link{confint.boot_ame}}.
#'
#' @param object an \code{ame_als} fit.
#' @param parm character vector of parameter names, or integer indices; if
#'   \code{NULL} (default) all available coefficients are returned.
#' @param level confidence level (default 0.95).
#' @param ... passed to \code{\link{vcov.ame_als}} (e.g. \code{cluster}).
#' @return A matrix with one row per coefficient and lower/upper bound columns.
#' @method confint ame_als
#' @seealso \code{\link{ame_als_bootstrap}} for bootstrap intervals on all
#'   parameters.
#' @importFrom stats vcov
#' @export
confint.ame_als <- function(object, parm = NULL, level = 0.95, ...) {
	if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
		cli::cli_abort("{.arg level} must be a single number in (0, 1).")
	}
	# unified bootstrap path: if the fit was produced via
	# ame_als(..., bootstrap = n) or otherwise has a `$bootstrap` slot,
	# return the bootstrap intervals (fully propagated) instead of the
	# anti-conservative sandwich wald intervals.
	if (!is.null(object$bootstrap)) {
		return(stats::confint(object$bootstrap, parm = parm, level = level, ...))
	}
	V   <- stats::vcov(object, ...)
	est <- object$coefficients[rownames(V)]
	se  <- sqrt(diag(V))
	a   <- (1 - level) / 2
	z   <- stats::qnorm(c(a, 1 - a))
	ci  <- cbind(est + z[1] * se, est + z[2] * se)
	rownames(ci) <- rownames(V)
	colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), "%")
	cli::cli_inform(c(
		"!" = "Conditional sandwich (Wald) intervals for the regression coefficients only -- anti-conservative.",
		"i" = "For fully-propagated intervals on all parameters use {.fn ame_als_bootstrap}."))
	if (!is.null(parm)) {
		if (is.character(parm)) {
			bad <- setdiff(parm, rownames(ci))
			if (length(bad)) {
				cli::cli_abort(c(
					"Unknown parameter name in {.arg parm}: {.val {bad}}.",
					"i" = "Available: {.val {rownames(ci)}}."))
			}
		} else if (is.numeric(parm) &&
		           (any(parm < 1) || any(parm > nrow(ci)))) {
			cli::cli_abort("{.arg parm} index is out of range (1:{nrow(ci)}).")
		}
		ci <- ci[parm, , drop = FALSE]
	}
	ci
}
