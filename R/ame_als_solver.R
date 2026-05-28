# ame_als_solver.R
#
# linear-solver factory for the regression block (mu, beta) of the fast AME
# estimator. `linear_solver` selects how the joint intercept + dyadic-covariate
# least-squares step is solved:
#   "eigen" - eigendecomposition of the weighted normal-equations matrix X'WX
#             (default; X'WX is precomputed once, minimum-norm on rank
#             deficiency). Forming X'WX squares the condition number of X.
#   "qr"    - QR factorisation of the weighted observed design, which solves the
#             least-squares step without squaring the condition number. For a
#             near-singular design (collinear covariates) the QR solution is
#             noise-dominated, so the solver falls back to the minimum-norm
#             "eigen" solution there and warns.
#   "auto"  - "qr" for a moderately ill-conditioned full-rank design, "eigen"
#             otherwise (minimum-norm when rank-deficient or near-singular).
#
# the factory returns a closure `solve(resid_arr)` that maps a residual array
# (Z - a - b - O, with the regression part removed) to the (p+1)-vector
# (mu, beta). The factorisation is fixed across BCD iterations, so it is built
# once; only the right-hand side changes.

.ae_linsolver <- function(linear_solver, X, mask, Warr, nr, nc, Tt, p) {
	pa <- p + 1L
	obs_idx <- which(mask)                         # linear indices, nr*nc*Tt
	Wv <- Warr[obs_idx]                            # weights at observed cells
	# observed-cell covariate columns (constant across iterations)
	Xcols <- matrix(0, length(obs_idx), p)
	if (p > 0) for (k in seq_len(p)) Xcols[, k] <- X[, , k, ][obs_idx]

	# weighted normal-equations matrix X'WX. its eigenvalues are the squared
	# singular values of sqrt(W) X, so they give both the numerical rank and the
	# true condition number of the design -- which drive the "auto" choice and
	# the rank-deficiency warning.
	D <- cbind(1, Xcols)
	XtX <- crossprod(D, Wv * D)
	eg  <- eigen((XtX + t(XtX)) / 2, symmetric = TRUE)
	lam <- eg$values
	lam_max <- max(lam, 0)
	rank_tol <- lam_max * max(pa, 1) * .Machine$double.eps * 100
	num_rank <- sum(lam > rank_tol)
	rank_deficient <- num_rank < pa
	# condition number of sqrt(W) X (= sqrt of the X'WX condition number)
	cond_x <- if (rank_deficient || lam_max <= 0) {
		Inf
	} else {
		sqrt(lam_max / min(lam[lam > rank_tol]))
	}

	# a design that is numerically full rank but near-singular has a unique
	# least-squares solution that is noise-dominated and can be enormous; the
	# QR path returns exactly that solution, whereas the eigen path truncates
	# and returns the stable minimum-norm solution. Treat such a design like a
	# rank-deficient one: prefer (and, for an explicit qr request, fall back to)
	# the minimum-norm solver.
	near_singular <- is.finite(cond_x) && cond_x > 1e5

	if (identical(linear_solver, "auto")) {
		linear_solver <- if (rank_deficient || near_singular) {
			"eigen"                                # minimum-norm if (near-)singular
		} else if (cond_x > 1e3) {
			"qr"                                   # ill-conditioned, full rank
		} else {
			"eigen"
		}
	} else if (identical(linear_solver, "qr") && near_singular) {
		cli::cli_warn(c(
			"{.arg linear_solver} = {.val qr}: the regression design is near-singular (condition number {signif(cond_x, 3)}).",
			"i" = "The QR least-squares solution would be noise-dominated; using the minimum-norm {.val eigen} solver instead."))
		linear_solver <- "eigen"
	}

	if (identical(linear_solver, "qr")) {
		sw  <- sqrt(Wv)
		qrD <- qr(D * sw)                          # weighted observed design
		solve_fun <- function(resid_arr) {
			cf <- tryCatch(qr.coef(qrD, sw * resid_arr[obs_idx]),
			               error = function(e) rep(NA_real_, pa))
			cf[!is.finite(cf)] <- 0                # aliased columns -> 0
			as.numeric(cf)
		}
	} else {
		solve_fun <- function(resid_arr) {
			wy  <- Wv * resid_arr[obs_idx]
			Xty <- c(sum(wy), if (p > 0) as.numeric(crossprod(Xcols, wy)))
			.ae_safe_solve(XtX, Xty)
		}
	}

	list(solve = solve_fun, rank_deficient = rank_deficient,
	     method = linear_solver)
}
