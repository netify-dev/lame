# linear solver for the regression block (mu, beta) of the ALS estimator, chosen by `linear_solver`:
#   "eigen" eigendecomposition of X'WX (default; minimum-norm on rank deficiency)
#   "qr"    QR of the weighted observed design, avoiding the squared condition number; a near-singular
#           design falls back to "eigen" with a warning
#   "auto"  "qr" for a moderately ill-conditioned full-rank design, "eigen" otherwise
# returns `solve(resid_arr)` giving (mu, beta) from a residual array, `solve_xty` (eigen only) from X'W r,
# and `XtX`; the factorisation is built once and reused across BCD iterations

.ae_linsolver <- function(
	linear_solver,
	X,
	mask,
	Warr,
	nr,
	nc,
	Tt,
	p,
	XtX = NULL
) {
	pa <- p + 1L
	# observed-cell design, built on first use; only the QR solver and the cell-level `solve` need it
	cells <- NULL
	observed_cells <- function() {
		if (is.null(cells)) {
			obs_idx <- which(mask) # linear indices into nr x nc x Tt
			Xcols <- matrix(0, length(obs_idx), p)
			if (p > 0) for (k in seq_len(p)) Xcols[, k] <- X[, , k, ][obs_idx]
			cells <<- list(obs_idx = obs_idx, Wv = Warr[obs_idx], Xcols = Xcols)
		}
		cells
	}

	# eigenvalues of X'WX give the numerical rank and the condition number of the weighted design
	if (is.null(XtX)) {
		oc <- observed_cells()
		D <- cbind(1, oc$Xcols)
		XtX <- crossprod(D, oc$Wv * D)
	}
	eg <- eigen((XtX + t(XtX)) / 2, symmetric = TRUE)
	lam <- eg$values
	lam_max <- max(lam, 0)
	rank_tol <- lam_max * max(pa, 1) * .Machine$double.eps * 100
	num_rank <- sum(lam > rank_tol)
	rank_deficient <- num_rank < pa
	# condition number of sqrt(W) X
	cond_x <- if (rank_deficient || lam_max <= 0) {
		Inf
	} else {
		sqrt(lam_max / min(lam[lam > rank_tol]))
	}

	# a near-singular full-rank design has a noise-dominated least-squares solution, so it gets the
	# truncated minimum-norm eigen solve like a rank-deficient one
	near_singular <- is.finite(cond_x) && cond_x > 1e5

	if (identical(linear_solver, "auto")) {
		linear_solver <- if (rank_deficient || near_singular) {
			"eigen" # minimum-norm if (near-)singular
		} else if (cond_x > 1e3) {
			"qr" # ill-conditioned, full rank
		} else {
			"eigen"
		}
	} else if (identical(linear_solver, "qr") && near_singular) {
		cli::cli_warn(c(
			"{.arg linear_solver} = {.val qr}: the regression design is near-singular (condition number {signif(cond_x, 3)}).",
			"i" = "The QR least-squares solution would be noise-dominated; using the minimum-norm {.val eigen} solver instead."
		))
		linear_solver <- "eigen"
	}

	if (identical(linear_solver, "qr")) {
		oc <- observed_cells()
		obs_idx <- oc$obs_idx
		sw <- sqrt(oc$Wv)
		qrD <- qr(cbind(1, oc$Xcols) * sw) # weighted observed design
		solve_fun <- function(resid_arr) {
			cf <- tryCatch(
				qr.coef(qrD, sw * resid_arr[obs_idx]),
				error = function(e) rep(NA_real_, pa)
			)
			cf[!is.finite(cf)] <- 0 # aliased columns -> 0
			as.numeric(cf)
		}
		solve_xty <- NULL
	} else {
		solve_fun <- function(resid_arr) {
			oc <- observed_cells()
			wy <- oc$Wv * resid_arr[oc$obs_idx]
			Xty <- c(sum(wy), if (p > 0) as.numeric(crossprod(oc$Xcols, wy)))
			.ae_safe_solve(XtX, Xty)
		}
		# the same solve from X'W r, for callers that accumulate it from sums (the BCD sweeps)
		solve_xty <- .ae_safe_solver(XtX)
	}

	list(
		solve = solve_fun,
		solve_xty = solve_xty,
		XtX = XtX,
		rank_deficient = rank_deficient,
		method = linear_solver
	)
}
