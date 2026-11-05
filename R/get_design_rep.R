#' Create design array for replicate data
#'
#' @usage get_design_rep(Y, Xdyad, Xrow, Xcol, actorSet, intercept, n, N, pr, pc, pd,
#' symmetric = FALSE)
#' @param Y dependent variable in array format
#' @param Xdyad dyadic covariates in array format
#' @param Xrow sender covariates in array format
#' @param Xcol receiver covariates in array format
#' @param actorSet vector of actors
#' @param intercept logical indicating whether to include intercept
#' @param n number of actors
#' @param N number of replicates
#' @param pr number of receiver covariates
#' @param pc number of sender covariates
#' @param pd number of dyadic covariates
#' @param symmetric logical; in a symmetric network a dyad has one covariate
#' value, so a dyad is unobserved in a period when the dyadic covariate of
#' either of its two cells is missing (\code{\link{lame}} fills a hole in one
#' triangle from the other before the design is built, so only a dyad missing
#' in both reaches this rule)
#' @return returns list of design array values necessary for ame_repL, plus
#' \code{intercept}: the intercept flag actually used (FALSE when the
#' covariates span a constant across the observed panel). A cell whose
#' covariate is missing is set to \code{NA} in the returned \code{Y} and its
#' design value to 0.
#' @author Shahryar Minhas
#'
#' @keywords internal
#' @export get_design_rep

get_design_rep <- function(
	Y = Y,
	Xdyad,
	Xrow,
	Xcol,
	actorSet,
	intercept,
	n,
	N,
	pr,
	pc,
	pd,
	symmetric = FALSE
) {

	# covariate slices, before the panel-wide intercept decision
	Xt_list <- vector("list", N)
	for (t in 1:N) {
		# keep the covariate-name dimension when there is a single covariate
		xr <- if (is.null(Xrow)) NULL else Xrow[, , t]
		if (!is.null(xr) && is.null(dim(xr))) {
			xr <- matrix(xr, n, 1, dimnames = list(NULL, dimnames(Xrow)[[2]]))
		}
		xc <- if (is.null(Xcol)) NULL else Xcol[, , t]
		if (!is.null(xc) && is.null(dim(xc))) {
			xc <- matrix(xc, n, 1, dimnames = list(NULL, dimnames(Xcol)[[2]]))
		}
		xd <- if (is.null(Xdyad)) NULL else Xdyad[, , , t]
		if (!is.null(xd) && length(dim(xd)) == 2L) {
			xd <- array(
				xd,
				c(n, n, 1),
				dimnames = list(NULL, NULL, dimnames(Xdyad)[[3]])
			)
		}
		Xt_list[[t]] <- design_array_listwisedel(xr, xc, xd, FALSE, n)
	}

	# the intercept is redundant only if the whole panel design spans it
	keep_int <- isTRUE(as.logical(intercept)) &&
		!.lame_panel_intercept_redundant(Xt_list, Y)
	if (isTRUE(as.logical(intercept)) && !keep_int) {
		cli::cli_warn(c(
			"The covariates are collinear with the intercept across the observed panel.",
			"i" = "The intercept is dropped; the covariates carry the level."
		))
	}

	X <- array(dim = c(n, n, pr + pc + pd + keep_int, N))
	for (t in 1:N) {
		Xt <- Xt_list[[t]]

		# apply the panel-wide intercept decision to every slice
		if (keep_int) {
			tmp <- array(
				dim = dim(Xt) + c(0, 0, 1),
				dimnames = list(NULL, NULL, c('intercept', dimnames(Xt)[[3]]))
			)
			tmp[, , 1] <- 1 ; tmp[, , -1] <- Xt
			Xt <- tmp
		}

		X[, , , t] <- Xt
	}
	if ((pr + pc + pd + keep_int) > 0) {
		dimnames(X) <- list(actorSet, actorSet, dimnames(Xt)[[3]], dimnames(Y)[[3]])
	}
	dimnames(X)[[4]] <- dimnames(Y)[[3]]

	# NAs in the design array become NAs in Y
	if ((pr + pc + pd) > 0) {
		for (t in 1:N) {
			for (p in 1:dim(X)[3]) {
				na_mat <- X[, , p, t] ; na_mat[!is.na(na_mat)] <- 1 ; na_mat[is.nan(na_mat)] <- NA
				Y[, , t] <- Y[, , t] * na_mat
			}
		} ; rm(na_mat)
		# a symmetric dyad has one covariate value, so a missing value in either cell leaves it unobserved
		if (isTRUE(symmetric)) {
			for (tt in 1:N) {
				miss <- rowSums(is.na(X[, , , tt, drop = FALSE]), dims = 2L) > 0
				both <- miss | t(miss)
				if (any(both)) { Yt <- Y[, , tt] ; Yt[both] <- NA ; Y[, , tt] <- Yt }
			}
		}
	}

	# NAs in the design array become zero
	X[is.na(X)] <- 0

	# design array transformations for the beta update
	Xlist <- lapply(1:N, function(t) { array(X[, , , t], dim = dim(X)[1:3], dimnames = dimnames(X)[1:3]) })
	XrLong <- apply(X, c(1, 3, 4), sum)
	XcLong <- apply(X, c(2, 3, 4), sum)
	mXLong <- apply(X, c(3, 4), c)
	mXtLong <- apply(aperm(X, c(2, 1, 3, 4)), c(3, 4), c)
	if ((pr + pc + pd + keep_int) > 0) {
		xxLong <- array(apply(mXLong, 3, function(x) { t(x) %*% x }), dim = c(dim(X)[3], dim(X)[3], N))
		xxTLong <- array(unlist(lapply(1:N, function(t) { t(mXLong[, , t]) %*% mXtLong[, , t] })), dim = dim(xxLong))
	} else { xxLong <- NULL ; xxTLong <- NULL }

	return(
		list(
			Y = Y,
			X = X,
			Xlist = Xlist,
			XrLong = XrLong,
			XcLong = XcLong,
			mXLong = mXLong,
			mXtLong = mXtLong,
			xxLong = xxLong,
			xxTLong = xxTLong,
			intercept = keep_int
		)
	)
}

# whether the covariates span a constant on the observed panel, reducing one
# period at a time by an orthogonal transformation to bound memory
.lame_panel_intercept_redundant <- function(xlist, y) {
	p <- dim(xlist[[1L]])[3L]
	if (p == 0L) return(FALSE)
	blocks <- lapply(seq_along(xlist), function(t) {
		x <- matrix(xlist[[t]], ncol = p)
		keep <- !is.na(c(y[, , t])) & rowSums(!is.finite(x)) == 0L
		if (!any(keep)) return(NULL)
		q <- qr(cbind(x[keep, , drop = FALSE], 1), LAPACK = TRUE)
		qr.R(q)[, order(q$pivot), drop = FALSE]
	})
	x <- do.call(rbind, blocks)
	if (is.null(x)) return(FALSE)
	scale <- sqrt(colSums(x^2))
	scale[scale == 0] <- 1
	x <- sweep(x, 2L, scale, "/")
	qr(x)$rank == qr(x[, seq_len(p), drop = FALSE])$rank
}
