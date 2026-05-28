# ame_als.R
#
# fast, MCMC-free point estimation for the AME model via iterative block
# coordinate descent.
#
# the estimation algorithm is adapted from the Social Influence Regression
# (SIR) estimator of Hoff & Minhas --- the "iterative block coordinate descent"
# method described in Minhas & Hoff (2025), "Decomposing Network Dynamics:
# Social Influence Regression", Political Analysis --- and implemented in
# sir::sir_alsfit() (nicknamed "ALS" in that package). The estimator here is
# a port/adaptation of that algorithm to the additive-and-multiplicative
# effects (AME) model; it is not original lame methodology.

# --------------------------------------------------------------------------
# internal helpers
# --------------------------------------------------------------------------

# null-default
.ae_default <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# pairwise, NA-aware symmetrization of a square matrix. Unlike (M + t(M))/2,
# this does not propagate NA: a dyad observed in only one direction keeps its
# observed value (mirrored), and only dyads missing in both directions stay NA.
.ae_symmetrize <- function(M) {
	Mt <- t(M)
	out <- M
	both  <- is.finite(M) & is.finite(Mt)
	onlyt <- !is.finite(M) & is.finite(Mt)
	out[both]  <- (M[both] + Mt[both]) / 2
	out[onlyt] <- Mt[onlyt]
	out
}

# solve a symmetric normal-equations system A x = b by eigendecomposition,
# dropping directions with negligible eigenvalues. A is always a cross-product
# here (X'X, D'WD), so symmetric positive-semidefinite. returns the minimum-norm
# least-squares solution: stable for an exactly singular A (where solve() errors)
# and for a near-singular A (where solve() returns a finite but inflated answer).
# MASS::ginv() is a fallback only if eigen() itself fails on a non-finite matrix.
.ae_safe_solve <- function(A, bvec) {
	A <- (A + t(A)) / 2
	eg <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
	if (is.null(eg)) {
		out <- tryCatch(as.numeric(MASS::ginv(A) %*% bvec),
		                error = function(e) rep(0, length(bvec)))
		return(as.numeric(out))
	}
	lam <- eg$values
	thr <- max(lam, 0) * sqrt(.Machine$double.eps)
	keep <- is.finite(lam) & lam > thr
	if (!any(keep)) return(rep(0, length(bvec)))
	Qk <- eg$vectors[, keep, drop = FALSE]
	as.numeric(Qk %*% (crossprod(Qk, bvec) / lam[keep]))
}

# families the BCD estimator supports
.ae_supported_families <- c("normal", "binary", "poisson")

# identify a node-covariate effect by an explicit orthogonality constraint.
# A node covariate broadcasts to a per-actor constant, collinear with the
# additive sender/receiver effect, so its coefficient is not identified by the
# objective alone. We regress the raw additive effect `e` on [1, W] and report
# the slope coefficients (the between-actor regression coefficients) and the
# residual `e - [1,W]g` -- the identified additive heterogeneity, orthogonal to
# 1 and to W. This is the standard estimand under the assumption that the
# additive effect is uncorrelated with the node covariate.
#
# `active` (0/1 per actor) excludes actors with no observed dyads: their
# additive effect is not identified (the BCD pins it to 0), so they must not
# enter the between-actor regression, where they would otherwise pull the
# intercept and slopes toward an arbitrary value.
.ae_decompose <- function(e, W, active = NULL) {
	if (is.null(W) || ncol(W) == 0L) {
		return(list(intercept = 0, beta = numeric(0), resid = e))
	}
	if (is.null(active)) active <- rep(1, length(e))
	D <- cbind(1, W)
	Dw <- D * active
	g <- .ae_safe_solve(crossprod(Dw, D), crossprod(Dw, e))
	bcoef <- g[-1]
	names(bcoef) <- colnames(W)
	list(intercept = g[1], beta = bcoef, resid = as.numeric(e - D %*% g))
}

# connected components of the observed-dyad graph (union-find). Returns the
# component label of each row and column node. Unipartite: rows and columns are
# the same n actors, so the two label vectors coincide. Bipartite: a graph over
# the n_row + n_col distinct nodes.
.ae_components <- function(obs_count, bip) {
	nr <- nrow(obs_count); ncc <- ncol(obs_count)
	if (bip) {
		n  <- nr + ncc
		ij <- which(obs_count > 0, arr.ind = TRUE)
		edges <- if (nrow(ij) == 0L) matrix(0L, 0, 2) else
			cbind(ij[, 1], nr + ij[, 2])
	} else {
		n  <- nr
		adj <- (obs_count > 0) | t(obs_count > 0)
		adj[lower.tri(adj, diag = TRUE)] <- FALSE
		ij <- which(adj, arr.ind = TRUE)
		edges <- if (nrow(ij) == 0L) matrix(0L, 0, 2) else ij
	}
	parent <- seq_len(n)
	find <- function(x) {
		root <- x
		while (parent[root] != root) root <- parent[root]
		while (parent[x] != root) { nx <- parent[x]; parent[x] <<- root; x <- nx }
		root
	}
	for (e in seq_len(nrow(edges))) {
		ri <- find(edges[e, 1]); rj <- find(edges[e, 2])
		if (ri != rj) parent[ri] <- rj
	}
	root <- vapply(seq_len(n), find, integer(1))
	lab  <- match(root, unique(root))
	if (bip) {
		list(row = lab[seq_len(nr)], col = lab[nr + seq_len(ncc)],
		     n = length(unique(lab)))
	} else {
		list(row = lab, col = lab, n = length(unique(lab)))
	}
}

# number of connected components (thin wrapper, kept for call-site clarity)
.ae_n_components <- function(obs_count, bip) .ae_components(obs_count, bip)$n

# fit-preserving per-component additive gauge (directed networks only). On a
# disconnected observed-dyad graph each component has a free shift
# a_i += d, b_j -= d that leaves every fitted value unchanged. We pin d per
# component to the precision-weighted minimum-norm value, so the reported a, b
# -- and hence va, vb and the node-covariate regression -- are reproducible
# rather than an artdefact of the optimisation path. (Symmetric models have no
# such per-component freedom: a_i + a_j is not invariant to a one-sided shift.)
.ae_component_gauge <- function(a, b, comp, w_row, w_col) {
	if (comp$n <= 1L) return(list(a = a, b = b))
	for (k in seq_len(comp$n)) {
		ri <- which(comp$row == k); ci <- which(comp$col == k)
		den <- sum(w_row[ri]) + sum(w_col[ci])
		if (den <= 0) next
		d <- (sum(w_col[ci] * b[ci]) - sum(w_row[ri] * a[ri])) / den
		a[ri] <- a[ri] + d
		b[ci] <- b[ci] - d
	}
	list(a = a, b = b)
}

# map an observed-outcome array to a Gaussian working response on the latent
# scale. For "normal" this is exact; for the discrete families it uses the
# same latent transforms ame()/lame() use to initialize their samplers
# (see get_start_vals.R) --- a fast, documented approximation.
.ae_working_response <- function(Yarr, family) {
	d <- dim(Yarr)
	if (family == "normal") {
		return(Yarr)
	}
	if (family == "poisson") {
		return(log(Yarr + 1))
	}
	# binary: rank-based normal scores, computed within each slice. ordinal
	# is not supported by ALS; ordinal inference goes through the MCMC path.
	Z <- array(NA_real_, dim = d)
	for (t in seq_len(d[3])) {
		Yt <- Yarr[, , t]
		zt <- matrix(zscores(c(Yt)), d[1], d[2])
		if (family == "binary") {
			has0 <- any(Yt == 0, na.rm = TRUE)
			has1 <- any(Yt == 1, na.rm = TRUE)
			if (has0 && has1) {
				z01 <- 0.5 * (max(zt[Yt == 0], na.rm = TRUE) +
				              min(zt[Yt == 1], na.rm = TRUE))
				zt <- zt - z01
			}
		}
		Z[, , t] <- zt
	}
	Z
}

# coerce one dyadic covariate slice to an n_row x n_col x pd array
.ae_dyad_to_3d <- function(x, nr, nc) {
	if (is.null(x)) return(NULL)
	if (length(dim(x)) == 3L) return(x)
	if (is.matrix(x)) return(array(x, dim = c(nrow(x), ncol(x), 1L)))
	cli::cli_abort("Each {.arg Xdyad} element must be a matrix or a 3D array.")
}

# build the canonical representation of the covariates.
#
# dyadic covariates become the design array X (n_row x n_col x pd x T) used by
# the block coordinate descent. Row/column (node) covariates are kept SEPARATE
# as per-actor matrices W_row (n_row x pr) and W_col (n_col x pc): a node
# covariate broadcasts to a per-actor constant, which is collinear with the
# additive sender/receiver effect, so its coefficient is not identified by the
# objective. It is identified instead by an explicit orthogonality constraint
# (the additive effect is taken orthogonal to the node covariates) and is
# estimated post-hoc in .ae_assemble() as a between-actor regression. A
# time-varying node covariate is summarised by its per-actor mean over time
# (the static ALS model carries one coefficient per covariate).
.ae_build_design <- function(Xdyad, Xrow, Xcol, nr, nc, Tt) {
	# --- dyadic covariates -> design array X ---
	dcols <- list(); dnms <- character(0)
	if (!is.null(Xdyad)) {
		pd <- dim(Xdyad[[1]])[3]
		dn0 <- dimnames(Xdyad[[1]])
		dnm <- if (length(dn0) >= 3L) dn0[[3]] else NULL
		for (k in seq_len(pd)) {
			arr <- array(0, dim = c(nr, nc, Tt))
			for (t in seq_len(Tt)) arr[, , t] <- Xdyad[[t]][, , k]
			dcols[[length(dcols) + 1L]] <- arr
			# match the MCMC suffix convention: dyadic covariates get the
			# "_dyad" suffix idempotently so downstream code filtering by
			# name resolves the same way across ALS and MCMC fits. The
			# package-internal `.lame_apply_suffix()` handles both pre-suffixed
			# user inputs and legacy `.dyad` upgrades.
			base_nm <- .ae_default(dnm[k], paste0("dyad", k))
			dnms <- c(dnms, .lame_apply_suffix(base_nm, "dyad"))
		}
	}
	pd <- length(dcols)
	X <- array(0, dim = c(nr, nc, pd, Tt))
	for (k in seq_len(pd)) X[, , k, ] <- dcols[[k]]
	# NAs in covariates are treated as 0 (matching design_array); warn when any
	# appear off the diagonal, where 0 is a genuine modelling assumption
	if (pd > 0) {
		n_na <- sum(is.na(X))
		n_diag_na <- 0L
		if (nr == nc && n_na > 0L) {
			for (k in seq_len(pd)) for (t in seq_len(Tt)) {
				n_diag_na <- n_diag_na + sum(is.na(diag(X[, , k, t])))
			}
		}
		if (n_na > n_diag_na) {
			cli::cli_warn(c(
				"Missing values in dyadic covariates were replaced with 0.",
				"i" = "A missing covariate value is modelled as 0; supply an explicit imputation if that is not intended."))
		}
	}
	X[is.na(X)] <- 0
	if (pd > 0) dimnames(X) <- list(NULL, NULL, dnms, NULL)

	# --- node covariates -> per-actor means ---
	node_means <- function(Xn, n, suffix, lbl) {
		if (is.null(Xn)) return(NULL)
		p <- ncol(as.matrix(Xn[[1]]))
		nm0 <- colnames(as.matrix(Xn[[1]]))
		W <- matrix(0, n, p)
		varies <- FALSE
		for (k in seq_len(p)) {
			vals <- matrix(
				vapply(Xn, function(m) as.matrix(m)[, k], numeric(n)),
				n, length(Xn))
			W[, k] <- rowMeans(vals, na.rm = TRUE)
			# detect within-actor (over-time) variation that the per-actor mean
			# discards: the ALS model is static and carries one coefficient
			if (length(Xn) > 1L && !varies) {
				within <- mean(apply(vals, 1, stats::sd, na.rm = TRUE),
				               na.rm = TRUE)
				if (is.finite(within) && within > 1e-8) varies <- TRUE
			}
		}
		if (varies) {
			cli::cli_warn(c(
				"A {lbl} covariate varies within actor across time slices.",
				"i" = "The ALS estimator is static: each node covariate is summarised by its per-actor mean over time, so within-actor variation is discarded.",
				"i" = "For time-varying node effects use {.fn lame} (MCMC)."))
		}
		if (any(!is.finite(W))) {
			cli::cli_warn("Missing values in node covariates were replaced with 0.")
		}
		W[!is.finite(W)] <- 0
		# idempotent suffix application: pre-suffixed inputs and legacy
		# `.row`/`.col` are normalised rather than double-suffixed
		base_nm <- .ae_default(nm0, paste0(lbl, seq_len(p)))
		colnames(W) <- .lame_apply_suffix(base_nm, sub("^_", "", suffix))
		W
	}
	W_row <- node_means(Xrow, nr, "_row", "row")
	W_col <- node_means(Xcol, nc, "_col", "col")

	list(X = X, p_dyad = pd, x_names_dyad = dnms,
	     W_row = W_row, W_col = W_col)
}

# normalize raw user inputs into the canonical representation:
#   Yarr  : n_row x n_col x T   (unipartite diagonal blanked to NA)
#   X     : n_row x n_col x p x T
.ae_prepare <- function(Y, Xdyad, Xrow, Xcol, mode, longitudinal) {

	# --- Y to a list of matrices ---
	if (longitudinal) {
		if (is.array(Y) && length(dim(Y)) == 3L) {
			Ylist <- lapply(seq_len(dim(Y)[3]), function(t) Y[, , t])
			names(Ylist) <- dimnames(Y)[[3]]
		} else if (is.list(Y)) {
			Ylist <- Y
		} else {
			cli::cli_abort(c(
				"{.arg Y} for {.fn lame_als} must be a list of matrices or a 3D array.",
				"i" = "Received an object of class {.cls {class(Y)[1]}}."))
		}
	} else {
		if (!is.matrix(Y)) {
			cli::cli_abort(c(
				"{.arg Y} for {.fn ame_als} must be a matrix.",
				"i" = "For longitudinal data use {.fn lame_als}."))
		}
		Ylist <- list(Y)
	}
	Tt <- length(Ylist)
	if (Tt < 1L) cli::cli_abort("{.arg Y} contains no network slices.")

	nr <- nrow(Ylist[[1]]); nc <- ncol(Ylist[[1]])
	rn1 <- rownames(Ylist[[1]]); cn1 <- colnames(Ylist[[1]])
	for (t in seq_len(Tt)) {
		if (nrow(Ylist[[t]]) != nr || ncol(Ylist[[t]]) != nc) {
			cli::cli_abort(c(
				"All network slices must share the same dimensions.",
				"i" = "Slice 1 is {nr}x{nc} but slice {t} is {nrow(Ylist[[t]])}x{ncol(Ylist[[t]])}.",
				"i" = "Changing actor compositions are not supported by the ALS estimator; use {.fn lame} (MCMC)."))
		}
		# guard against silent actor misalignment: same-sized slices whose
		# named actors differ would be treated positionally and give wrong results
		if (longitudinal && t > 1L) {
			rnt <- rownames(Ylist[[t]]); cnt <- colnames(Ylist[[t]])
			if ((!is.null(rn1) && !is.null(rnt) && !identical(rnt, rn1)) ||
			    (!is.null(cn1) && !is.null(cnt) && !identical(cnt, cn1))) {
				cli::cli_abort(c(
					"All network slices must have the same actors in the same order.",
					"i" = "Slice {t} has different row/column names than slice 1.",
					"i" = "Changing actor compositions are not supported by the ALS estimator; use {.fn lame} (MCMC)."))
			}
		}
	}

	row_names <- .ae_default(rownames(Ylist[[1]]), paste0("row", seq_len(nr)))
	col_names <- .ae_default(colnames(Ylist[[1]]), paste0("col", seq_len(nc)))
	time_names <- .ae_default(names(Ylist), paste0("t", seq_len(Tt)))

	bip <- identical(mode, "bipartite")
	# a unipartite network has a single node set: row and column actors coincide
	if (!bip) col_names <- row_names

	# --- Y to canonical array, blank unipartite diagonal ---
	Yarr <- array(NA_real_, dim = c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		yt <- as.matrix(Ylist[[t]])
		storage.mode(yt) <- "double"
		if (!bip) diag(yt) <- NA_real_
		Yarr[, , t] <- yt
	}

	# warn on a time slice with no observed (non-NA, off-diagonal) entries: it
	# contributes nothing to the fit and usually signals a data problem
	if (Tt > 1L) {
		empty <- which(vapply(seq_len(Tt),
			function(t) !any(is.finite(Yarr[, , t])), logical(1)))
		if (length(empty) > 0L) {
			cli::cli_warn(c(
				"Time slice{?s} {.val {empty}} {?has/have} no observed entries.",
				"i" = "Empty slices contribute nothing to the fit; check for an all-missing network."))
		}
	}

	# --- covariates to per-slice lists ---
	norm_dyad <- function(Xd) {
		if (is.null(Xd)) return(NULL)
		if (!longitudinal) {
			Xd <- list(Xd)
		} else if (is.array(Xd) && length(dim(Xd)) == 4L) {
			Xd <- lapply(seq_len(dim(Xd)[4]), function(t) Xd[, , , t])
		} else if (is.array(Xd) && length(dim(Xd)) == 3L && !is.list(Xd)) {
			# a 3D array is ambiguous for a longitudinal fit (is the 3rd
			# dimension covariates or time?) -- ask for an unambiguous form
			cli::cli_abort(c(
				"{.arg Xdyad} for {.fn lame_als} is a 3D array, which is ambiguous.",
				"i" = "Supply a list of {Tt} matrices/arrays (one per time slice),",
				"i" = "or a 4D array {.code [n_row, n_col, p, T]}."))
		} else if (!is.list(Xd)) {
			Xd <- list(Xd)
		}
		if (length(Xd) != Tt) {
			cli::cli_abort("{.arg Xdyad} must have one entry per time slice ({Tt}).")
		}
		lapply(Xd, .ae_dyad_to_3d, nr = nr, nc = nc)
	}
	norm_node <- function(Xn, arg) {
		if (is.null(Xn)) return(NULL)
		if (!longitudinal || !is.list(Xn)) Xn <- list(Xn)
		if (length(Xn) != Tt) {
			cli::cli_abort("{.arg {arg}} must have one entry per time slice ({Tt}).")
		}
		lapply(Xn, as.matrix)
	}

	Xdyad <- norm_dyad(Xdyad)
	Xrow  <- norm_node(Xrow, "Xrow")
	Xcol  <- norm_node(Xcol, "Xcol")

	if (!is.null(Xdyad)) {
		pd1 <- dim(Xdyad[[1]])[3]
		for (t in seq_len(Tt)) {
			dd <- dim(Xdyad[[t]])
			if (dd[1] != nr || dd[2] != nc) {
				cli::cli_abort(c(
					"{.arg Xdyad} dimensions must match {.arg Y}.",
					"i" = "Y is {nr}x{nc}; Xdyad slice {t} is {dd[1]}x{dd[2]}."))
			}
			if (dd[3] != pd1) {
				cli::cli_abort(c(
					"All {.arg Xdyad} slices must hold the same number of covariates.",
					"i" = "Slice 1 has {pd1}; slice {t} has {dd[3]}."))
			}
		}
	}
	if (!is.null(Xrow)) {
		pr1 <- ncol(Xrow[[1]])
		for (t in seq_len(Tt)) {
			if (nrow(Xrow[[t]]) != nr || ncol(Xrow[[t]]) != pr1) {
				cli::cli_abort(c(
					"All {.arg Xrow} slices must be {nr} x {pr1} (rows = row/sender nodes).",
					"i" = "Slice {t} is {nrow(Xrow[[t]])} x {ncol(Xrow[[t]])}."))
			}
		}
	}
	if (!is.null(Xcol)) {
		pc1 <- ncol(Xcol[[1]])
		for (t in seq_len(Tt)) {
			if (nrow(Xcol[[t]]) != nc || ncol(Xcol[[t]]) != pc1) {
				cli::cli_abort(c(
					"All {.arg Xcol} slices must be {nc} x {pc1} (rows = column/receiver nodes).",
					"i" = "Slice {t} is {nrow(Xcol[[t]])} x {ncol(Xcol[[t]])}."))
			}
		}
	}

	des <- .ae_build_design(Xdyad, Xrow, Xcol, nr, nc, Tt)

	list(Yarr = Yarr, X = des$X, p_dyad = des$p_dyad,
	     x_names_dyad = des$x_names_dyad,
	     W_row = des$W_row, W_col = des$W_col,
	     nr = nr, nc = nc, Tt = Tt, bip = bip,
	     row_names = row_names, col_names = col_names, time_names = time_names)
}

# --------------------------------------------------------------------------
# core estimator: iterative block coordinate descent
# --------------------------------------------------------------------------

# exact solve of the two-way weighted additive normal equations
#   [ diag(cnt_row)  wsum         ] [a]   [base_rs]
#   [ wsum'          diag(cnt_col)] [b] = [base_cs]
# (symmetric models collapse to (diag(cnt_row) + wsum) a = base_rs, b = a).
# used as an exact fallback when the Gauss-Seidel sweep crawls under strongly
# unbalanced observation weights (an IRLS fit). the non-symmetric system is
# rank-deficient by one (the a_i + b_j gauge); a symmetric-eigen pseudo-inverse
# returns the minimum-norm solution, which the caller re-centres anyway and
# which stays well-defined on disconnected observed-dyad graphs.
.ae_additive_solve <- function(base_rs, base_cs, wsum,
                               cnt_row, cnt_col, symmetric) {
	psolve <- function(K, rhs) {
		eg   <- eigen((K + t(K)) / 2, symmetric = TRUE)
		thr  <- max(abs(eg$values)) * sqrt(.Machine$double.eps)
		keep <- abs(eg$values) > thr
		Qk   <- eg$vectors[, keep, drop = FALSE]
		as.vector(Qk %*% (crossprod(Qk, rhs) / eg$values[keep]))
	}
	if (symmetric) {
		K <- diag(x = cnt_row, nrow = length(cnt_row)) + wsum
		a <- psolve(K, base_rs)
		a[cnt_row == 0] <- 0
		return(list(a = a, b = a))
	}
	np <- length(cnt_row); nq <- length(cnt_col)
	K  <- matrix(0, np + nq, np + nq)
	K[seq_len(np), seq_len(np)]           <- diag(x = cnt_row, nrow = np)
	K[np + seq_len(nq), np + seq_len(nq)] <- diag(x = cnt_col, nrow = nq)
	K[seq_len(np), np + seq_len(nq)]      <- wsum
	K[np + seq_len(nq), seq_len(np)]      <- t(wsum)
	sol <- psolve(K, c(base_rs, base_cs))
	a <- sol[seq_len(np)]; b <- sol[np + seq_len(nq)]
	a[cnt_row == 0] <- 0; b[cnt_col == 0] <- 0
	list(a = a, b = b)
}

# Z       : n_row x n_col x T Gaussian working response (NA = unobserved)
# X       : n_row x n_col x p x T design array
# returns : list(mu, beta, a, b, U, V, L, Omat, sse, iterations, converged,
#                dev_history)
.ae_bcd_fit <- function(Z, X, R, symmetric, max_iter, tol,
                        init = NULL, lowrank_method = "mm",
                        weights = NULL, linear_solver = "eigen") {

	# guard against inputs that would otherwise raise an opaque low-level error
	# (a vector tol short-circuits the convergence test; max_iter < 1 returns a
	# degenerate, never-iterated fit)
	if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1) {
		cli::cli_abort("{.arg max_iter} must be a single positive integer.")
	}
	if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
		cli::cli_abort("{.arg tol} must be a single positive number.")
	}

	d  <- dim(Z)
	nr <- d[1]; nc <- d[2]; Tt <- d[3]
	p  <- dim(X)[3]

	mask <- is.finite(Z)
	if (!any(mask)) {
		cli::cli_abort(c(
			"No observed cells: the working response is entirely missing.",
			"i" = "Check that {.arg Y} has observed (non-NA, off-diagonal) entries."))
	}
	Zf <- Z; Zf[!mask] <- 0

	# per-cell observation weights: unit by default (so every block reduces
	# exactly to the unweighted estimator), family working weights under IRLS
	Warr <- array(0, dim = d)
	if (is.null(weights)) {
		Warr[mask] <- 1
	} else {
		ww <- weights
		ww[!mask | !is.finite(ww)] <- 0
		ww[ww < 0] <- 0
		Warr <- ww
	}
	# weighted observed-cell totals per dyad (the c_ij of the low-rank objective
	# and the denominators of the additive update; integer counts when unit)
	wsum <- apply(Warr, c(1, 2), sum)
	cnt_row <- rowSums(wsum)
	cnt_col <- colSums(wsum)
	cnt_row_pos <- pmax(cnt_row, 1)
	cnt_col_pos <- pmax(cnt_col, 1)

	# --- initialization (warm start when supplied) ---
	mu   <- .ae_default(init$mu, 0)
	beta <- .ae_default(init$beta, rep(0, p))
	if (length(beta) != p) beta <- rep(0, p)
	a <- .ae_default(init$a, rep(0, nr))
	b <- .ae_default(init$b, rep(0, nc))
	if (length(a) != nr) a <- rep(0, nr)
	if (length(b) != nc) b <- rep(0, nc)
	U <- init$U; V <- init$V; L <- init$L
	Omat <- matrix(0, nr, nc)
	if (R > 0 && !is.null(U)) {
		if (symmetric && !is.null(L)) {
			Omat <- U %*% (L * t(U))
		} else if (!is.null(V)) {
			Omat <- tcrossprod(U, V)
		}
	}

	xbeta_arr <- function(bt) {
		out <- array(0, dim = c(nr, nc, Tt))
		if (p > 0) for (t in seq_len(Tt)) {
			m <- matrix(0, nr, nc)
			for (k in seq_len(p)) m <- m + bt[k] * X[, , k, t]
			out[, , t] <- m
		}
		out
	}
	Xb <- xbeta_arr(beta)

	# precompute the regression-block solver for the augmented design
	# [intercept | dyadic covariates]. The design and the observation mask are
	# fixed across iterations, so the factorisation (X'WX eigen, or QR of the
	# observed design) is built once; only the right-hand side is recomputed.
	pa  <- p + 1L
	lin <- .ae_linsolver(linear_solver, X, mask, Warr, nr, nc, Tt, p)
	# warn once if the design is rank-deficient (constant / collinear covariate)
	if (lin$rank_deficient) {
		cli::cli_warn(c(
			"The regression design (intercept + dyadic covariates) is rank-deficient.",
			"i" = "Coefficients are not uniquely identified; the minimum-norm (pseudoinverse) solution is used.",
			"i" = "Check for a constant covariate (the intercept already covers it) or collinear covariates."))
	}

	sse_old <- Inf
	coef_mb_old <- rep(Inf, pa)
	dev_history <- numeric(0)
	converged <- FALSE
	add_converged <- TRUE          # did the final additive block solve converge?
	iter <- 0L

	for (iter in seq_len(max_iter)) {

		ab <- outer(a, b, "+")

		# --- block: intercept + regression coefficients (joint OLS) ---
		# solved jointly so a constant covariate is detected as collinear with
		# the intercept rather than aliased across two separate blocks. The
		# precomputed `lin` solver maps the residual target to (mu, beta).
		resid_arr <- Zf - c(ab) - c(Omat)          # broadcast ab, Omat over t
		coef_mb <- lin$solve(resid_arr)
		# change in the gauge-invariant coefficients (mu, beta) since last sweep;
		# used as a secondary convergence signal alongside the SSE plateau
		coef_delta <- max(abs(coef_mb - coef_mb_old))
		coef_mb_old <- coef_mb
		mu   <- coef_mb[1]
		beta <- coef_mb[-1]
		Xb   <- xbeta_arr(beta)

		# --- block: additive sender/receiver effects (a, b) ---
		# weighted residual not involving a or b (Warr is 0 off the mask)
		base_a <- array(0, dim = c(nr, nc, Tt))
		for (t in seq_len(Tt)) {
			base_a[, , t] <- Warr[, , t] *
				(Zf[, , t] - mu - Xb[, , t] - Omat)
		}
		base_rs <- apply(base_a, 1, sum)
		base_cs <- apply(base_a, 2, sum)
		# inner Gauss-Seidel sweeps to a tolerance. Each update averages over
		# OBSERVED cells only (a_i = mean_{obs j,t}(base - b_j)), so the missing
		# diagonal is handled exactly and every sub-step is an exact conditional
		# minimizer. The fixed point solves the two-way additive normal
		# equations exactly; the residual check below certifies it was reached.
		# symmetric models use the native single-vector recursion, whose fixed
		# point (diag(cnt) + wsum) a = base_rs is the exact minimizer of the
		# symmetric two-way sub-problem.
		for (inner in seq_len(200L)) {
			a_prev <- a; b_prev <- b
			if (symmetric) {
				a <- (base_rs - as.vector(wsum %*% a)) / cnt_row_pos
				a[cnt_row == 0] <- 0
				b <- a
			} else {
				a <- (base_rs - as.vector(wsum %*% b)) / cnt_row_pos
				a[cnt_row == 0] <- 0
				b <- (base_cs - as.vector(crossprod(wsum, a))) / cnt_col_pos
				b[cnt_col == 0] <- 0
			}
			if (max(abs(a - a_prev)) < 1e-10 &&
			    max(abs(b - b_prev)) < 1e-10) break
		}
		# residual certificate: the Gauss-Seidel fixed point solves the two-way
		# normal equations, so a non-negligible residual flags the sweep cap
		# being hit short of the exact additive minimizer
		ra <- cnt_row * a + as.vector(wsum %*% b) - base_rs
		rb <- cnt_col * b + as.vector(crossprod(wsum, a)) - base_cs
		ra[cnt_row == 0] <- 0; rb[cnt_col == 0] <- 0
		add_scale <- max(1, max(abs(base_rs)), max(abs(base_cs)))
		# Gauss-Seidel can crawl when the observation weights are strongly
		# unbalanced (an IRLS fit); fall back to an exact direct solve of the
		# additive normal equations so a, b are the true minimiser regardless
		if (max(max(abs(ra)), max(abs(rb))) > 1e-6 * add_scale) {
			adir <- .ae_additive_solve(base_rs, base_cs, wsum,
			                           cnt_row, cnt_col, symmetric)
			a <- adir$a; b <- adir$b
			ra <- cnt_row * a + as.vector(wsum %*% b) - base_rs
			rb <- cnt_col * b + as.vector(crossprod(wsum, a)) - base_cs
			ra[cnt_row == 0] <- 0; rb[cnt_col == 0] <- 0
		}
		# the flag reflects the FINAL outer iteration's additive block
		add_converged <- max(max(abs(ra)), max(abs(rb))) <= 1e-6 * add_scale
		# sum-to-zero identification: fold the removed grand mean into mu so
		# centering does not change the fitted values
		mu <- mu + mean(a) + mean(b)
		a <- a - mean(a)
		b <- b - mean(b)
		ab <- outer(a, b, "+")

		# --- block: multiplicative latent factors (U, V) ---
		# minimise the time-averaged weighted low-rank objective
		#   sum_ij c_ij (Rbar_ij - O_ij)^2 ,  rank(O) <= R,
		# where c_ij = wsum[i,j] is the (weighted) number of observed slices for
		# dyad (i,j) and Rbar is the weighted per-dyad mean residual.
		# `lowrank_method` selects the inner solver -- weighted majorise-
		# minimise ("mm", the default) or alternating least squares ("als" /
		# "hybrid"); each is a monotone, objective-non-increasing update. See
		# .ae_lowrank_update() in ame_als_lowrank.R.
		if (R > 0) {
			res_sum <- matrix(0, nr, nc)
			for (t in seq_len(Tt)) {
				res_sum <- res_sum +
					Warr[, , t] * (Zf[, , t] - mu - Xb[, , t] - ab)
			}
			obs_dyad <- wsum > 0
			Rbar <- matrix(0, nr, nc)
			Rbar[obs_dyad] <- res_sum[obs_dyad] / wsum[obs_dyad]
			lr <- .ae_lowrank_update(Omat, Rbar, wsum, obs_dyad,
			                         R, symmetric, lowrank_method)
			U <- lr$U; V <- lr$V; L <- lr$L; Omat <- lr$Omat
		}

		# --- convergence: relative change in (weighted) residual sum of squares
		sse <- 0
		for (t in seq_len(Tt)) {
			r <- Zf[, , t] - mu - Xb[, , t] - ab - Omat
			sse <- sse + sum(Warr[, , t] * r^2)
		}
		dev_history <- c(dev_history, sse)
		# skip the relative-change test on the first sweep (sse_old is still Inf).
		# require BOTH a residual-SSE plateau AND stable (mu, beta): the second
		# guard rules out a flat-SSE saddle along which the coefficients still
		# drift. The coefficient tolerance is relative and loose, so it never
		# blocks a genuine convergence.
		coef_stable <- coef_delta <= tol * 100 * (max(abs(coef_mb)) + 1e-8)
		if (iter > 1L && is.finite(sse) && coef_stable &&
		    abs(sse_old - sse) / (abs(sse_old) + 0.1) < tol) {
			converged <- TRUE
			sse_old <- sse
			break
		}
		sse_old <- sse
	}

	list(mu = mu, beta = beta, a = a, b = b, U = U, V = V, L = L,
	     Omat = Omat, sse = sse_old, iterations = iter,
	     converged = converged, add_converged = add_converged,
	     dev_history = dev_history)
}

# assemble the fitted-model object from a BCD solution
.ae_assemble <- function(fit, prep, family, mode, symmetric, R,
                         longitudinal, call, Zwork = NULL, link = NULL) {

	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt
	# `Zwork` (the IRLS final working response) overrides the one-shot transform
	Z  <- if (is.null(Zwork)) .ae_working_response(prep$Yarr, family) else Zwork
	# the symmetric fit minimises against a symmetrised working response;
	# symmetrise here too so residuals / ve match the fitted objective
	if (symmetric) for (t in seq_len(Tt)) {
		Z[, , t] <- .ae_symmetrize(Z[, , t])
	}
	mask <- is.finite(Z)

	mu <- fit$mu
	beta_dyad <- fit$beta              # dyadic-covariate coefficients
	a_raw <- fit$a; b_raw <- fit$b     # raw additive effects (absorb node structure)
	if (length(beta_dyad) > 0) names(beta_dyad) <- prep$x_names_dyad

	U <- fit$U; V <- fit$V; L <- fit$L
	Omat <- fit$Omat

	# --- identification: gauge the multiplicative term O to be double-centered ---
	# the additive term a_i + b_j and the rank-R multiplicative term u_i' v_j are
	# not separately identified for R > 0: a broadcast (row- or column-constant)
	# component can sit in either (a 1' is rank 1, so a pure sender effect can be
	# absorbed by a column of U with the matching column of V constant). We impose
	# the standard AME gauge by double-centering O (zero row and column means) and
	# folding its row means into the sender effect, its column means into the
	# receiver effect and its grand mean into mu. This is fit-preserving (EZ is
	# unchanged) and makes the reported a, b, U, V unique given the fit ---
	# independent of which BCD fixed point was reached.
	if (R > 0 && !is.null(Omat)) {
		rmn <- rowMeans(Omat); cmn <- colMeans(Omat); gmn <- mean(Omat)
		Omat <- Omat - outer(rmn, rep(1, nc)) - outer(rep(1, nr), cmn) + gmn
		a_raw <- a_raw + rmn
		b_raw <- b_raw + cmn
		mu <- mu - gmn
		# re-impose sum-to-zero, folding the grand mean back into mu
		mu <- mu + mean(a_raw) + mean(b_raw)
		a_raw <- a_raw - mean(a_raw)
		b_raw <- b_raw - mean(b_raw)
		# recompute U, V (or U, L) from the now double-centered O. Double-centering
		# cannot raise rank, so O still has rank <= R and this factorisation is
		# exact (O = U V', or U L U' when symmetric).
		if (symmetric) {
			eg <- eigen((Omat + t(Omat)) / 2, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
			U <- eg$vectors[, idx, drop = FALSE]
			L <- eg$values[idx]
			V <- U
		} else {
			sv <- svd(Omat)
			dd <- sqrt(pmax(sv$d[seq_len(R)], 0))
			U <- sweep(sv$u[, seq_len(R), drop = FALSE], 2, dd, "*")
			V <- sweep(sv$v[, seq_len(R), drop = FALSE], 2, dd, "*")
			L <- NULL
		}
	}
	if (!is.null(U)) dimnames(U) <- list(prep$row_names, paste0("dim", seq_len(R)))
	if (!is.null(V)) dimnames(V) <- list(prep$col_names, paste0("dim", seq_len(R)))

	# --- active actors + per-component additive gauge ---
	# observed-dyad counts per actor (symmetrised for undirected models so a
	# one-triangle input does not under-count); an actor with none is an
	# isolate whose additive effect is not identified.
	obs_count <- apply(!is.na(prep$Yarr), c(1, 2), sum)
	if (symmetric) obs_count <- obs_count + t(obs_count)
	cnt_row <- rowSums(obs_count); cnt_col <- colSums(obs_count)
	active_row <- as.numeric(cnt_row > 0)
	active_col <- as.numeric(cnt_col > 0)
	# on a disconnected observed graph each component has a free a_i += d,
	# b_j -= d shift; pin it to the precision-weighted minimum-norm value so the
	# reported a, b are reproducible. Directed networks only -- a symmetric
	# a_i + a_j has no per-component freedom.
	if (!symmetric) {
		comp <- .ae_components(obs_count, prep$bip)
		if (comp$n > 1L) {
			cg <- .ae_component_gauge(a_raw, b_raw, comp, cnt_row, cnt_col)
			a_raw <- cg$a; b_raw <- cg$b
		}
	}

	# --- linear predictor (link scale), from the RAW additive effects ---
	EZ <- vector("list", Tt)
	ab_raw <- outer(a_raw, b_raw, "+")
	resid_list <- vector("list", Tt)
	ve_num <- 0; ve_den <- 0
	for (t in seq_len(Tt)) {
		xb <- matrix(0, nr, nc)
		if (prep$p_dyad > 0) for (k in seq_len(prep$p_dyad)) {
			xb <- xb + beta_dyad[k] * prep$X[, , k, t]
		}
		ez <- mu + xb + ab_raw + Omat
		EZ[[t]] <- ez
		r <- Z[, , t] - ez
		resid_list[[t]] <- r
		ve_num <- ve_num + sum(r[mask[, , t]]^2)
		ve_den <- ve_den + sum(mask[, , t])
	}
	ve <- if (ve_den > 0) ve_num / ve_den else NA_real_

	# --- identify node-covariate coefficients by the orthogonality constraint ---
	# A node covariate broadcasts to a per-actor constant, collinear with the
	# additive effect; beta_row/beta_col are the between-actor regression
	# coefficients and a/b the residual heterogeneity orthogonal to them. This
	# is a fit-preserving reparameterisation of the raw additive effects.
	if (symmetric) {
		W_node <- prep$W_row
		if (!is.null(prep$W_col)) {
			W_node <- if (is.null(W_node)) prep$W_col else cbind(W_node, prep$W_col)
		}
		dec <- .ae_decompose(a_raw, W_node, active_row)
		a <- dec$resid; b <- dec$resid
		mu <- mu + 2 * dec$intercept       # the intercept enters g_i and g_j
		beta <- c(beta_dyad, dec$beta)
	} else {
		dec_r <- .ae_decompose(a_raw, prep$W_row, active_row)
		dec_c <- .ae_decompose(b_raw, prep$W_col, active_col)
		a <- dec_r$resid; b <- dec_c$resid
		mu <- mu + dec_r$intercept + dec_c$intercept
		beta <- c(beta_dyad, dec_r$beta, dec_c$beta)
	}
	names(a) <- prep$row_names
	names(b) <- prep$col_names

	# --- variance components (VC columns mirror ame()/lame()) ---
	# va/vb are empirical variances of the (covariate-orthogonal) additive
	# effects -- not identical to a random-effect variance component. Isolated
	# actors (no observed dyad) have an unidentified additive effect and are
	# excluded from the variance summaries.
	va <- if (sum(active_row) > 1) stats::var(a[active_row == 1]) else NA_real_
	vb <- if (symmetric) {
		NA_real_
	} else if (sum(active_col) > 1) {
		stats::var(b[active_col == 1])
	} else {
		NA_real_
	}
	# no sender/receiver covariance for symmetric or bipartite networks
	cab <- if (symmetric || prep$bip) {
		NA_real_
	} else {
		both <- active_row == 1 & active_col == 1
		if (sum(both) > 2) stats::cov(a[both], b[both]) else NA_real_
	}
	rho <- 0
	if (!prep$bip && !symmetric) {
		# pool (i,j) vs (j,i) residual pairs across time slices
		e_ij <- numeric(0); e_ji <- numeric(0)
		for (t in seq_len(Tt)) {
			rt <- resid_list[[t]]
			ut <- upper.tri(rt)
			e_ij <- c(e_ij, rt[ut])
			e_ji <- c(e_ji, t(rt)[ut])
		}
		ok <- is.finite(e_ij) & is.finite(e_ji)
		if (sum(ok) > 2) {
			rho <- suppressWarnings(stats::cor(e_ij[ok], e_ji[ok]))
			if (!is.finite(rho)) rho <- 0
		}
	}
	VC <- c(va = va, cab = cab, vb = vb, rho = rho, ve = ve)

	# --- response-scale fitted values; blank the unipartite diagonal ---
	# IRLS fits carry a genuine GLM `link`; the transform path keeps its own
	# inverse maps (poisson log1p, binary probit-type latent)
	fitted_list <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		fitted_list[[t]] <- if (!is.null(link)) {
			switch(paste(family, link),
				"poisson log"   = { m <- exp(EZ[[t]]); m[m > 1e6] <- 1e6; m },
				"binary logit"  = stats::plogis(EZ[[t]]),
				"binary probit" = stats::pnorm(EZ[[t]]),
				EZ[[t]])
		} else switch(family,
			normal  = EZ[[t]],
			poisson = { m <- exp(EZ[[t]]) - 1; m[m < 0] <- 0; m[m > 1e6] <- 1e6; m },
			binary  = stats::pnorm(EZ[[t]]))
		if (!prep$bip) {                  # self-ties are undefined
			diag(EZ[[t]]) <- NA_real_
			diag(fitted_list[[t]]) <- NA_real_
		}
		dimnames(fitted_list[[t]]) <- list(prep$row_names, prep$col_names)
		dimnames(EZ[[t]]) <- list(prep$row_names, prep$col_names)
	}
	names(EZ) <- prep$time_names
	names(fitted_list) <- prep$time_names

	coefficients <- c(intercept = mu, beta)

	approx_note <- if (family == "normal") {
		if (R > 0) {
			paste0("Least-squares (independent-error Gaussian ML) AME point ",
			       "estimate; with R > 0 the low-rank term is non-convex, so this ",
			       "is a local optimum.")
		} else {
			paste0("Least-squares (independent-error Gaussian ML) AME point ",
			       "estimate; the global solution for R = 0.")
		}
	} else if (!is.null(link)) {
		paste0("Approximate point estimate: iteratively reweighted least squares ",
		       "(", family, ", ", link, " link) -- a fast approximate GLM AME fit. ",
		       "For calibrated family-specific inference use ame()/lame() (MCMC).")
	} else {
		paste0("Approximate point estimate: iterative block coordinate descent ",
		       "run on a Gaussian working response (latent ", family,
		       " scores). For calibrated family-specific inference use ame()/lame() (MCMC).")
	}

	out <- list(
		call = call,
		family = family, mode = mode, symmetric = symmetric, R = R,
		longitudinal = longitudinal,
		mu = mu, beta = beta, a = a, b = b,
		U = U, V = V, L = L, G = if (prep$bip && R > 0) diag(R) else NULL,
		coefficients = coefficients,
		VC = VC, s2 = ve,
		EZ = EZ, fitted = fitted_list, YPM = fitted_list,
		residuals = resid_list,
		iterations = fit$iterations, converged = fit$converged,
		deviance = fit$sse, dev_history = fit$dev_history,
		n_time = Tt,
		dims = list(n_row = nr, n_col = nc, p = length(beta)),
		row_names = prep$row_names, col_names = prep$col_names,
		x_names = names(beta), time_names = prep$time_names,
		Y = prep$Yarr, X = prep$X, W_row = prep$W_row, W_col = prep$W_col,
		meta = list(
			sampler = "bcd",
			algorithm = "iterative block coordinate descent (Hoff & Minhas SIR estimator)",
			uncertainty_available = FALSE,
			approximation_note = approx_note)
	)
	# longitudinal fits inherit from "lame_als" so the lame_als-specific S3
	# methods (tidy.lame_als / glance.lame_als / update.lame_als) dispatch
	# correctly; both subclasses inherit from "ame_als" so the existing
	# generic S3 methods keep working
	class(out) <- if (isTRUE(longitudinal)) c("lame_als", "ame_als") else "ame_als"
	out
}

# attach a bootstrap result to a freshly-fit ame_als / lame_als object so
# the unified `bootstrap = N` argument can return point + intervals in one call
.ae_attach_bootstrap <- function(fit, R_boot, type, block_length, seed,
                                 verbose) {
	# warn early for tiny R; the existing post-fit n_valid<10 warning catches
	# replicate failures, this one warns just on the user-requested count
	if (R_boot < 30L) {
		cli::cli_warn(c(
			"{.arg bootstrap} = {R_boot} is below the recommended minimum of 30 replicates.",
			"i" = "Standard errors will be unreliable. Consider {.code bootstrap = 200}."))
	}
	# suppress the bootstrap function's own "too few replicates" warn when we
	# already pre-warned above -- otherwise the user sees two warnings about
	# the same condition. Other warnings from ame_als_bootstrap still surface.
	bt <- withCallingHandlers(
		ame_als_bootstrap(fit, R = R_boot, type = type,
		                  block_length = block_length, seed = seed,
		                  verbose = verbose),
		warning = function(w) {
			if (R_boot < 30L &&
			    grepl("bootstrap replicates? (total|succeeded)",
			          conditionMessage(w))) {
				invokeRestart("muffleWarning")
			}
		})
	fit$bootstrap <- bt
	# meta-bookkeeping so downstream accessors (confint/summary/print) know
	# uncertainty is now available without re-running the bootstrap function
	if (is.null(fit$meta)) fit$meta <- list()
	fit$meta$uncertainty_available <- TRUE
	fit$meta$uncertainty_source    <- "bootstrap"
	fit$meta$bootstrap_R           <- R_boot
	fit$meta$bootstrap_type        <- type
	fit
}

# shared driver for ame_als() and lame_als()
.ae_run <- function(Y, Xdyad, Xrow, Xcol, R, family, mode, symmetric,
                    max_iter, tol, verbose, seed, longitudinal, call,
                    lowrank_method = "mm", non_normal_method = "irls",
                    link = "probit", linear_solver = "eigen",
                    multistart = "none") {

	# accept short amen-era family spellings as aliases (parity with ame() / lame())
	amen_aliases <- c(nrm = "normal", bin = "binary", ord = "ordinal",
	                  pois = "poisson", tob = "tobit")
	if (length(family) == 1L && is.character(family) &&
	    family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
			"i" = "Update your script to {.code family = \"{new_family}\"} when convenient."))
		family <- new_family
	}
	family <- match.arg(family, c("normal", "binary", "poisson", "ordinal",
	                              "tobit", "cbin", "frn", "rrl"))
	if (!family %in% .ae_supported_families) {
		supported <- .ae_supported_families
		cli::cli_abort(c(
			"family = {.val {family}} is not supported by the ALS estimator.",
			"i" = "Supported families: {.val {supported}}.",
			"i" = "{.val {family}} is a rank/censoring likelihood with no well-defined Gaussian working response; use {.fn ame} or {.fn lame} (MCMC) instead."))
	}
	mode <- match.arg(mode, c("unipartite", "bipartite"))
	# the block coordinate descent is deterministic; `seed` is accepted only
	# for signature parity with ame()/lame() and is intentionally not used (so
	# the routine does not mutate the global RNG state).

	# Y must be numeric and free of Inf/NaN -- otherwise a non-numeric Y is
	# silently coerced to NA and Inf is silently treated as missing
	y_vals <- if (is.list(Y)) unlist(Y, use.names = FALSE) else Y
	if (!is.numeric(y_vals) && !is.logical(y_vals)) {
		cli::cli_abort("{.arg Y} must be numeric.")
	}
	if (any(is.infinite(y_vals) | is.nan(y_vals))) {
		cli::cli_abort(c(
			"{.arg Y} contains infinite or NaN values.",
			"i" = "Only finite values and {.val NA} (missing) are allowed."))
	}
	# reject Inf / NaN in covariates: otherwise node_means' rowMeans propagates
	# Inf into the design and the resulting non-finite cells get zero-filled
	# silently (mislabelled as "missing values" in the warning)
	for (cov_nm in c("Xdyad", "Xrow", "Xcol")) {
		cov_val <- get(cov_nm)
		if (!is.null(cov_val)) {
			vals <- if (is.list(cov_val)) {
				unlist(lapply(cov_val, function(z) as.numeric(z)), use.names = FALSE)
			} else as.numeric(cov_val)
			if (any(is.infinite(vals) | is.nan(vals))) {
				cli::cli_abort(c(
					"{.arg {cov_nm}} contains infinite or NaN values.",
					"i" = "Only finite values and {.val NA} (missing) are allowed."))
			}
		}
	}

	prep <- .ae_prepare(Y, Xdyad, Xrow, Xcol, mode, longitudinal)

	if (mode == "unipartite" && prep$nr != prep$nc) {
		cli::cli_abort(c(
			"Network is rectangular ({prep$nr}x{prep$nc}) but mode = {.val unipartite}.",
			"i" = "Use mode = {.val bipartite} for rectangular networks."))
	}
	if (symmetric && prep$bip) {
		cli::cli_warn("{.arg symmetric} is ignored for bipartite networks.")
		symmetric <- FALSE
	}
	if (symmetric && prep$nr != prep$nc) {
		cli::cli_abort("Symmetric models require a square network.")
	}
	# verify Y is actually symmetric -- otherwise the estimator silently
	# averages the upper and lower triangles and discards the asymmetry
	if (symmetric && prep$nr == prep$nc) {
		for (t in seq_len(prep$Tt)) {
			Yt <- prep$Yarr[, , t]
			off <- Yt - t(Yt)
			if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
				cli::cli_abort(c(
					"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
					"i" = "Symmetrize the input or set {.arg symmetric} = FALSE."))
			}
		}
	}

	if (length(R) != 1L || anyNA(R) || R < 0 ||
	    !isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	R <- as.integer(R)
	rmax <- min(prep$nr, prep$nc)
	if (R >= rmax) {
		cli::cli_abort(c(
			"{.arg R} = {R} must be smaller than min(n_row, n_col) = {rmax}."))
	}

	# symmetric model: symmetrise the dyadic covariates so the fitted linear
	# predictor beta' x_ij is itself symmetric
	if (symmetric && prep$p_dyad > 0) {
		for (k in seq_len(prep$p_dyad)) for (t in seq_len(prep$Tt)) {
			prep$X[, , k, t] <- .ae_symmetrize(prep$X[, , k, t])
		}
	}

	# warn if the observed-dyad graph is disconnected: the additive effects
	# then have a per-component level ambiguity beyond the global sum-to-zero
	ncomp <- .ae_n_components(apply(!is.na(prep$Yarr), c(1, 2), sum), prep$bip)
	if (ncomp > 1L) {
		cli::cli_warn(c(
			"The observed-dyad graph has {ncomp} disconnected components.",
			"i" = "Additive effects are identified only up to a per-component level shift; compare fitted values / predictions rather than raw a, b across components."))
	}

	# family / data sanity checks
	obs_y <- prep$Yarr[is.finite(prep$Yarr)]
	if (family == "poisson" && length(obs_y) && any(obs_y < 0)) {
		cli::cli_abort(c(
			"family = {.val poisson} requires non-negative counts.",
			"i" = "{.arg Y} has negative observed values."))
	}
	if (family == "binary" && length(obs_y) && !all(obs_y %in% c(0, 1))) {
		cli::cli_warn(c(
			"family = {.val binary} expects 0/1 data, but {.arg Y} has other values.",
			"i" = "The fit will proceed on rank-based scores; did you mean {.val ordinal} or {.val normal}?"))
	}

	# warn on a degenerate binary slice (a single distinct observed value):
	# it carries no within-slice rank information, so its working response
	# is constant and only the intercept/additive terms see it
	if (family == "binary") {
		for (t in seq_len(prep$Tt)) {
			yt <- prep$Yarr[, , t]
			if (length(unique(yt[is.finite(yt)])) == 1L) {
				cli::cli_warn(c(
					"Time slice {t} has a single distinct observed {family} value.",
					"i" = "A degenerate slice carries no within-slice rank information; its working-response scores are constant."))
			}
		}
	}

	Z <- .ae_working_response(prep$Yarr, family)
	if (symmetric) for (t in seq_len(prep$Tt)) {
		Z[, , t] <- .ae_symmetrize(Z[, , t])
	}

	p_total <- prep$p_dyad +
		(if (is.null(prep$W_row)) 0L else ncol(prep$W_row)) +
		(if (is.null(prep$W_col)) 0L else ncol(prep$W_col))
	if (verbose) {
		cli::cli_h3("AME exploration (iterative block coordinate descent)")
		cli::cli_text("Network: {prep$nr} x {prep$nc} {ifelse(prep$bip, 'bipartite', 'unipartite')}, {prep$Tt} time slice{?s}, R = {R}")
		cli::cli_text("Family: {.field {family}}{ifelse(symmetric, ' (symmetric)', '')} | Covariates: {p_total}")
	}

	# ALS / hybrid low-rank solvers are directed/bipartite only
	if (symmetric && lowrank_method != "mm" && R > 0) {
		cli::cli_warn(c(
			"{.arg lowrank_method} = {.val {lowrank_method}} is not available for symmetric models.",
			"i" = "The symmetric eigenmodel uses the majorise-minimise solver."))
		lowrank_method <- "mm"
	}

	# IRLS is the default for binary / poisson. tobit and ordinal are not
	# supported by ALS; their inference goes through the MCMC ame() / lame()
	# path.
	use_irls <- identical(non_normal_method, "irls") &&
		family %in% c("binary", "poisson")
	eff_link <- if (family == "poisson") "log" else link
	# IRLS produces strongly unbalanced working weights; the majorise-minimise
	# low-rank solver converges too slowly under them, so a directed R > 0 IRLS
	# fit uses the hybrid (ALS) solver instead
	if (use_irls && R > 0L && !symmetric && lowrank_method == "mm") {
		lowrank_method <- "hybrid"
	}

	if (use_irls) {
		ir <- .ae_irls_run(prep, family, eff_link, R, symmetric,
		                   max_iter, tol, lowrank_method = lowrank_method,
		                   linear_solver = linear_solver)
		fit <- ir$fit
		if (verbose) {
			cli::cli_alert_info("IRLS ({family}, {eff_link} link): {ir$irls_iters} reweighting iteration{?s}.")
		}
	} else if (multistart != "none" && R > 0) {
		n_starts <- if (multistart == "full") 8L else 4L
		fit <- .ae_multistart_fit(Z, prep$X, R, symmetric, max_iter, tol,
		                          lowrank_method, linear_solver,
		                          n_starts, seed)
		if (verbose) {
			cli::cli_alert_info("Multi-start: {n_starts} starts; kept start {fit$multistart_selected}.")
		}
		ss <- fit$multistart_sse
		spread <- (max(ss) - min(ss)) / max(1, abs(min(ss)))
		if (spread > 1e-4) {
			cli::cli_warn(c(
				"Multi-start found materially different local optima (relative SSE spread {signif(spread, 2)}).",
				"i" = "The lowest-SSE fit was kept; with R > 0 the optimum is not unique."))
		}
	} else {
		# nudge multi-start when R > 0 (nonconvex objective)
		if (R > 0L && multistart == "none" && verbose) {
			cli::cli_inform(c(
				"i" = "Single-start ALS at {.code R = {R}}. The latent objective is nonconvex.",
				"i" = "Rerun with {.code multistart = \"cheap\"} (4 starts) or {.code \"full\"} (8 starts) to check sensitivity."))
		}
		fit <- .ae_bcd_fit(Z, prep$X, R, symmetric, max_iter, tol,
		                   lowrank_method = lowrank_method,
		                   linear_solver = linear_solver)
	}

	if (verbose) {
		if (fit$converged) {
			cli::cli_alert_success("Converged in {fit$iterations} iteration{?s} (deviance/SSE = {round(fit$sse, 4)}).")
		} else {
			cli::cli_alert_warning("Did not converge in {max_iter} iterations; returning the last iterate.")
		}
	}
	if (!fit$converged) {
		cli::cli_warn(c(
			"Block coordinate descent did not converge within {.arg max_iter} = {max_iter} iterations.",
			"i" = "Increase {.arg max_iter} or relax {.arg tol}."))
	}
	if (!isTRUE(fit$add_converged)) {
		cli::cli_warn(c(
			"The additive-effects (Gauss-Seidel) block did not reach its exact solution within the sweep cap.",
			"i" = "This happens on weakly connected observed-dyad graphs; the additive effects are approximate.",
			"i" = "Treat {.code a}/{.code b} as approximate, or fit with {.fn lame} (MCMC)."))
	}

	out <- if (use_irls) {
		.ae_assemble(fit, prep, family, mode, symmetric, R, longitudinal, call,
		             Zwork = ir$Z, link = eff_link)
	} else {
		.ae_assemble(fit, prep, family, mode, symmetric, R, longitudinal, call)
	}
	out$lowrank_method <- lowrank_method
	out$meta$lowrank_method <- lowrank_method
	out$non_normal_method <- if (use_irls) "irls" else "transform"
	out$link <- if (use_irls) eff_link else NULL
	out$linear_solver <- linear_solver
	out$multistart <- multistart
	# final IRLS weights (NULL = unit weights); used by the weighted sandwich
	out$obs_weights <- if (use_irls) ir$W else NULL
	if (!is.null(fit$multistart_sse)) {
		out$multistart_sse <- fit$multistart_sse
		out$multistart_selected <- fit$multistart_selected
	}
	out
}

# --------------------------------------------------------------------------
# public estimators
# --------------------------------------------------------------------------

#' Fast (MCMC-free) AME estimation for a cross-sectional network
#'
#' @description
#' Fits an additive and multiplicative effects (AME) model to a single
#' cross-sectional network by \strong{iterative block coordinate descent},
#' producing a fast point estimate with no MCMC and no credible intervals.
#'
#' The estimation algorithm adapts the \emph{iterative block coordinate
#' descent} estimator of the Social Influence Regression (SIR) model of
#' Hoff & Minhas (the algorithm implemented in \code{sir::sir_alsfit()} and
#' nicknamed "ALS" in that package) to the AME model. It is a port and
#' adaptation, not original \pkg{lame} methodology.
#'
#' @details
#' The AME decomposition
#' \deqn{z_{ij} = \mu + \beta' x_{ij} + a_i + b_j + u_i' v_j + \epsilon_{ij}}
#' is conditionally linear: it is a linear regression in \eqn{(\mu, \beta, a, b)}
#' for fixed multiplicative factors \eqn{(U, V)}, and the optimal rank-\code{R}
#' \eqn{(U, V)} for fixed \eqn{(\mu, \beta, a, b)} is the truncated SVD of the
#' residual matrix. The estimator therefore cycles, until the residual sum of
#' squares stabilises, through three blocks --- each a monotone
#' (objective-non-increasing) update: a joint least-squares solve for the
#' intercept \eqn{\mu} and regression coefficients \eqn{\beta} (with a
#' \code{\link[MASS]{ginv}} pseudoinverse fallback for rank-deficient designs);
#' Gauss-Seidel sweeps to convergence for the additive effects \eqn{(a, b)};
#' and a weighted low-rank update for the multiplicative factors \eqn{(U, V)}
#' via SVD (eigen-decomposition when \code{symmetric = TRUE}). The weighting
#' makes the multiplicative step correct for unbalanced longitudinal panels
#' (dyads observed at unequal numbers of time points); severely unbalanced
#' panels may need more iterations to converge, so raise \code{max_iter} if the
#' fit reports non-convergence.
#'
#' For \code{family = "normal"} this is a least-squares (Gaussian maximum-
#' likelihood) fit --- the exact global solution when \code{R = 0}, and a
#' local optimum of the non-convex low-rank objective when \code{R > 0}.
#' For \code{"binary"}, \code{"poisson"} and
#' \code{"ordinal"} the same algorithm is run on a Gaussian working response
#' (the latent transforms \code{ame()} uses to initialise its sampler): a fast
#' \emph{approximation} suitable for exploration and starting values. For
#' \code{"tobit"} an EM outer loop replaces censored cells with their
#' truncated-normal conditional mean and the BCD M-step is warm-started from
#' the previous iterate (see the \dQuote{Tobit attenuation} section below for
#' the documented \eqn{\sigma^2} bias and the moderate-to-heavy censoring
#' regimes where the point estimates of \eqn{\beta} can drift noticeably from
#' the MCMC estimator). For calibrated family-specific inference, use the
#' MCMC estimator \code{\link{ame}}. The remaining censoring/rank families
#' (\code{"cbin"}, \code{"frn"}, \code{"rrl"}) are not supported and raise
#' an informative error.
#'
#' @section Tobit attenuation:
#' The EM-ALS path for \code{family = "tobit"} is a deterministic generalised
#' EM that imputes censored \eqn{Z_{ij}} on the truncated-normal mean
#' \eqn{\eta_{ij} - \sigma\,\phi(-\eta_{ij}/\sigma)/\Phi(-\eta_{ij}/\sigma)} and
#' refits the BCD low-rank fit on the imputed working response. The reported
#' \eqn{\sigma^2} carries a downward (attenuation) bias that scales with the
#' censoring fraction: the additive sender/receiver effects \eqn{a_i, b_j}
#' absorb part of the censored signal during the M-step, leaving the
#' residual scale smaller than the truth. In the package's own internal
#' simulations (\eqn{n = 80}, dyadic covariate, true \eqn{\sigma^2 = 1}) the
#' reported \eqn{\sigma^2} runs roughly 0.2-0.7 at 50\% censoring and below
#' 0.1 at 80\% censoring, while the MCMC tobit recovers \eqn{\sigma^2}
#' essentially without bias. The censored-cell EM also drifts the intercept
#' and the dyadic slope: in the same simulation EM-ALS gives
#' \eqn{(\hat\mu, \hat\beta_\mathrm{dyad}) \approx (-1.4, 0.31)} at 50\%
#' censoring and \eqn{(-0.9, 0.03)} at 80\% censoring (truth
#' \eqn{(0, 0.5)} and \eqn{(-1.5, 0.5)} respectively), while MCMC tobit on
#' the same draws gives \eqn{(0.11, 0.47)} and \eqn{(-1.28, 0.42)}. Treat the
#' EM-ALS tobit point estimates as fast exploratory summaries; use
#' \code{\link{ame}(..., family = "tobit")} for calibrated point and
#' interval estimates. The EM iteration count and \eqn{\sigma^2} are
#' surfaced on the fit as \code{$em_iters} and \code{$sigma2}.
#'
#' \strong{Row/column (node) covariates.} A node covariate broadcasts to a
#' per-actor constant, which is collinear with the additive sender/receiver
#' effect, so its coefficient is not identified by the objective alone. It is
#' identified here by an explicit constraint: the additive effects are taken
#' orthogonal to the node covariates, and \code{beta_row}/\code{beta_col} are
#' the corresponding between-actor regression coefficients (the additive
#' effects then carry only the residual heterogeneity). This is the standard
#' estimand under the assumption that the additive effects are uncorrelated
#' with the node covariates; if that assumption is doubtful the coefficient
#' absorbs the covariate-correlated part of the additive heterogeneity.
#' Time-varying node covariates are summarised by their per-actor mean; if a
#' node covariate varies within actor over time the discarded within-actor
#' variation triggers a warning.
#'
#' \strong{Identifiability of the additive and multiplicative terms.} For
#' \code{R > 0} the additive term \eqn{a_i + b_j} and the multiplicative term
#' \eqn{u_i' v_j} are not separately identified by the objective alone: a
#' broadcast (row- or column-constant) component can sit in either, since a
#' pure sender effect \eqn{a 1'} is itself rank one. The estimator imposes the
#' standard AME gauge --- the multiplicative term is double-centered (zero row
#' and column means), so all broadcast structure is carried by \eqn{a, b} ---
#' which makes the reported \eqn{a, b, U, V} unique given the fitted values and
#' independent of the optimisation path. For a unipartite network the
#' double-centering means are taken over the full matrix, which includes the
#' structurally unobserved self-tie diagonal that the model fills in by its
#' low-rank completion; the additive/multiplicative split therefore carries an
#' \eqn{O(1/n)} dependence on that completion. On a disconnected observed-dyad
#' graph the additive effects additionally have a per-component level shift,
#' which is pinned to a precision-weighted minimum-norm gauge so that
#' \code{a}, \code{b} and the variance components remain reproducible. A dyadic
#' covariate that is itself (near) low-rank can still be partially aliased with
#' the multiplicative term, so keep \code{R} modest relative to the covariate
#' structure. This residual aliasing is intrinsic to the AME model --- the MCMC
#' estimator resolves it only through its priors.
#'
#' \strong{Choosing R.} There is no automatic order-selection criterion (the
#' working-response objective has no likelihood, so AIC/BIC do not apply). Fit a
#' few values --- e.g. \code{lapply(0:4, function(r) ame_als(Y, R = r,
#' ...))} --- and inspect \code{deviance} (the residual sum of squares): it
#' falls steeply while real multiplicative signal is being captured and then
#' flattens, so the \dQuote{elbow} of that curve is a reasonable choice. With
#' \code{R > 0} the objective is non-convex; use \code{multistart} to guard
#' against local optima.
#'
#' Uncertainty is obtained separately, by the bootstrap; see
#' \code{\link{ame_als_bootstrap}}. (The SIR paper's own primary standard
#' errors are Hessian-based, classical and sandwich/robust. For AME a
#' Hessian-based variance is awkward on two counts: without an explicit gauge
#' fix the rotational invariance of the multiplicative factors leaves the joint
#' Hessian rank-deficient, and even with a gauge fixed the Gaussian
#' working-response approximation used for the non-normal families leaves a
#' Hessian-based variance miscalibrated. The bootstrap side-steps both, so it is
#' preferred here.)
#'
#' @param Y a square (unipartite) or rectangular (bipartite) relational matrix.
#' @param Xdyad an \code{n_row x n_col} matrix or \code{n_row x n_col x pd}
#'   array of dyadic covariates, or \code{NULL}.
#' @param Xrow an \code{n_row x pr} matrix of row/sender covariates, or \code{NULL}.
#' @param Xcol an \code{n_col x pc} matrix of column/receiver covariates, or \code{NULL}.
#' @param R integer dimension of the multiplicative effects (default \code{0}).
#'   \strong{The covariate coefficients are conditional on this choice.} The
#'   multiplicative term \eqn{u_i'v_j} is a flexible high-variance regressor
#'   that can correlate with the dyadic covariates, so the estimated
#'   \code{beta} can shift -- and occasionally change sign -- as \code{R}
#'   increases. Validate against an \code{R = 0} fit and the MCMC
#'   \code{\link{ame}}/\code{\link{lame}} before reporting.
#' @param family one of \code{"normal"}, \code{"binary"}, \code{"poisson"},
#'   \code{"ordinal"}.
#' @param mode \code{"unipartite"} (square) or \code{"bipartite"} (rectangular).
#' @param symmetric logical; fit a symmetric (undirected) model. Unipartite only.
#' @param max_iter maximum number of block coordinate descent iterations (default 200).
#' @param tol convergence tolerance on the relative change in residual sum of
#'   squares (default \code{1e-6}).
#' @param lowrank_method inner solver for the multiplicative (low-rank) block:
#'   \code{"mm"} (default) weighted majorise-minimise; \code{"als"} alternating
#'   least squares; \code{"hybrid"} runs both and keeps the lower-objective
#'   result. \code{"als"}/\code{"hybrid"} converge faster than \code{"mm"} on
#'   strongly unbalanced longitudinal panels and are available for directed and
#'   bipartite models (symmetric fits always use \code{"mm"}). All three
#'   minimise the same objective, so the point estimate is unchanged for
#'   balanced data.
#' @param non_normal_method for the non-normal families, \code{"irls"}
#'   (default for \code{binary}/\code{poisson}) runs iteratively reweighted
#'   least squares -- a fast approximate GLM AME fit with coefficients on
#'   the calibrated link scale (Poisson log, binary logit/probit).
#'   \code{"transform"} fits a single fixed Gaussian working response
#'   (\code{log(y+1)} for Poisson, rank-normal scores for binary/ordinal);
#'   its coefficients are on an uncalibrated rank scale -- good for
#'   direction/ranking of effects, not for magnitudes. \strong{For ordinal
#'   data the transform path attenuates \eqn{\beta} toward zero (classical
#'   discretisation bias); prefer \code{"em"} when \eqn{\beta} is an effect-size
#'   estimand.} \code{"em"} (for \code{ordinal}, opt-in; default for
#'   \code{tobit}) runs a deterministic EM outer loop that maximises the
#'   observed ordinal log-likelihood, jointly estimating cutpoints
#'   (returned on \code{fit$alpha}, anchored at \code{alpha[1] = 0}) and the
#'   regression block; the surrogate is monotone non-decreasing across EM
#'   iterations. \code{"irls"} is supported for \code{binary}/\code{poisson};
#'   \code{ordinal} accepts \code{"transform"} (default, back-compat) or
#'   \code{"em"}. A directed \code{R > 0} IRLS fit uses the hybrid low-rank
#'   solver internally (the IRLS weights are unbalanced). Uncertainty for
#'   any of these is the bootstrap.
#' @param link link for \code{non_normal_method = "irls"} with a \code{binary}
#'   family: \code{"probit"} (default; matches \code{\link{ame}}/\code{\link{lame}})
#'   or \code{"logit"}. \code{poisson} always uses the log link;
#'   ignored otherwise.
#' @param linear_solver solver for the regression block: \code{"eigen"}
#'   (default) eigendecomposes the normal equations; \code{"qr"} uses a QR
#'   factorisation of the observed design, which is more stable for
#'   ill-conditioned covariates; \code{"auto"} picks \code{"qr"} when the design
#'   is ill-conditioned but full rank. All give the same answer for
#'   well-conditioned designs.
#' @param multistart for \code{R > 0} (a non-convex objective), \code{"none"}
#'   (default) fits from a single deterministic start; \code{"cheap"} (4 starts)
#'   and \code{"full"} (8 starts) also try random low-rank starts and keep the
#'   lowest-SSE fit, warning when the starts reach materially different optima.
#'   Reproducible given \code{seed}; the global RNG stream is left unchanged.
#' @param verbose logical; print progress (default \code{TRUE}).
#' @param seed random seed (default \code{6886}). The block coordinate descent
#'   is deterministic, so the point estimate is reproducible regardless; the
#'   argument is retained for API consistency with \code{\link{ame}}.
#'
#' @return An object of class \code{"ame_als"}: a list with the point
#'   estimates \code{mu}, \code{beta}, \code{a}, \code{b}, \code{U}, \code{V}
#'   (\code{L} for symmetric models), the linear predictor \code{EZ},
#'   response-scale \code{fitted} values, working-scale \code{residuals},
#'   convergence information, and the variance-component vector \code{VC} with
#'   five descriptive entries:
#'   \describe{
#'     \item{\code{va}, \code{vb}}{empirical variances of the sender and
#'       receiver additive effects.}
#'     \item{\code{cab}}{covariance of the sender and receiver effects
#'       (\code{NA} for symmetric or bipartite models).}
#'     \item{\code{rho}}{dyadic residual reciprocity --- the correlation between
#'       the residuals of \eqn{(i,j)} and \eqn{(j,i)} --- not the sender/receiver
#'       correlation \code{cor(a, b)}.}
#'     \item{\code{ve}}{residual (working-scale) variance.}
#'   }
#'   These are descriptive summaries of the point estimates, not random-effect
#'   variance components. See \code{\link{ame_als_bootstrap}} for
#'   uncertainty and \code{\link{vcov.ame_als}} for a fast analytic
#'   covariance of the regression coefficients.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The iterative block
#' coordinate descent estimator adapted here originates with that work
#' (implemented in \code{sir::sir_alsfit()}).
#'
#' @seealso \code{\link{lame_als}} for longitudinal networks,
#'   \code{\link{ame_als_bootstrap}} for bootstrap uncertainty,
#'   \code{\link{ame}} for the full MCMC estimator.
#'
#' @examples
#' Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' coef(fit)
#'
#' @param bootstrap integer: if > 0, additionally run \code{bootstrap} replicates
#'   of the parametric or block bootstrap (via \code{\link{ame_als_bootstrap}})
#'   after the point fit and attach the result as \code{fit$bootstrap}. The
#'   downstream accessors (\code{\link{confint.ame_als}}, \code{summary},
#'   \code{print}) then surface bootstrap intervals instead of the
#'   anti-conservative sandwich Wald intervals. Default \code{0} (no
#'   bootstrap) -- bootstrap is expensive, so it is not imposed on a user
#'   who just wants a quick fit.
#' @param bootstrap_type character: \code{"parametric"} (default) or
#'   \code{"block"} -- the bootstrap scheme to use when \code{bootstrap > 0}.
#' @param bootstrap_block_length integer: block length for the block
#'   bootstrap; only used when \code{bootstrap_type = "block"}.
#' @param bootstrap_seed optional integer seed for the bootstrap (the point
#'   fit uses \code{seed}).
#'
#' @section Limitations vs. \code{\link{ame}} / \code{\link{lame}}:
#' The ALS estimator is a fast, frequentist point estimator. It is not a
#' drop-in replacement for the MCMC estimators -- the following features
#' apply to MCMC only:
#' \itemize{
#'   \item \strong{Families.} ALS supports \code{normal}, \code{binary},
#'     \code{poisson}, \code{ordinal} (transform or EM), and \code{tobit}
#'     (EM only; see the \dQuote{Tobit attenuation} section for the
#'     documented \eqn{\sigma^2} downward bias and the censoring regimes
#'     where the point estimates drift from the MCMC estimator). It does
#'     \strong{not} support \code{cbin}, \code{frn}, \code{rrl}; for those,
#'     fall back to \code{\link{ame}} / \code{\link{lame}}. The
#'     \code{bootstrap = N} option is \strong{not} calibrated for
#'     \code{family = "tobit"} -- the parametric simulator falls through to
#'     the binary/probit branch, producing degenerate 0/1 resamples and
#'     coefficient bootstraps collapsed near zero. Use the MCMC tobit for
#'     uncertainty quantification.
#'   \item \strong{Dynamic effects.} \code{\link{lame_als}} fits a STATIC
#'     model pooled across time slices. There is no \code{dynamic_uv},
#'     \code{dynamic_ab}, \code{dynamic_G}, or \code{dynamic_beta}.
#'     \code{lame(..., method = "als")} warns and ignores those arguments.
#'   \item \strong{Priors.} ALS has no priors. \code{prior = list(...)} and
#'     \code{g = ...} are MCMC-only; the dispatcher \code{ame(..., method = "als")}
#'     warns and ignores them.
#'   \item \strong{Posterior quantities.} No \code{$BETA} / \code{$VC} posterior
#'     draws on the point fit. Use \code{bootstrap = N} for a sampling
#'     distribution analogue; the draws are bootstrap, not Bayesian.
#'   \item \strong{Multi-chain.} Not applicable; \code{n_chains} is dropped.
#'   \item \strong{Convergence diagnostics.} No Rhat / ESS / \code{\link{trace_plot}}.
#' }
#'
#' Features that work the same on ALS fits: \code{coef},
#' \code{vcov} (sandwich on regression block), \code{confint} (auto-routes
#' bootstrap intervals when present, sandwich Wald otherwise), \code{predict},
#' \code{fitted}, \code{residuals}, \code{summary}, \code{print},
#' \code{nobs}, \code{\link{simulate.ame_als}},
#' \code{\link{gof_plot.ame_als}}, \code{\link{ab_plot.ame_als}},
#' \code{\link{uv_plot}}, \code{\link{latent_positions}}.
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
ame_als <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
                    R = 0, family = "normal",
                    mode = c("unipartite", "bipartite"),
                    symmetric = FALSE, max_iter = 200, tol = 1e-6,
                    lowrank_method = c("mm", "als", "hybrid"),
                    non_normal_method = c("irls", "transform"),
                    link = c("probit", "logit"),
                    linear_solver = c("eigen", "qr", "auto"),
                    multistart = c("none", "cheap", "full"),
                    bootstrap = 0L,
                    bootstrap_type = c("parametric", "block"),
                    bootstrap_block_length = 1L,
                    bootstrap_seed = NULL,
                    verbose = TRUE, seed = 6886) {
	mode <- match.arg(mode)
	lowrank_method <- match.arg(lowrank_method)
	non_normal_method <- match.arg(non_normal_method)
	link <- match.arg(link)
	linear_solver <- match.arg(linear_solver)
	multistart <- match.arg(multistart)
	bootstrap_type <- match.arg(bootstrap_type)
	if (!is.numeric(bootstrap) || length(bootstrap) != 1L || bootstrap < 0L ||
	    !isTRUE(all.equal(bootstrap, round(bootstrap)))) {
		cli::cli_abort("{.arg bootstrap} must be a non-negative integer (0 = no bootstrap).")
	}
	bootstrap <- as.integer(bootstrap)
	fit <- .ae_run(Y, Xdyad, Xrow, Xcol, R, family, mode, symmetric,
	               max_iter, tol, verbose, seed, longitudinal = FALSE,
	               call = match.call(), lowrank_method = lowrank_method,
	               non_normal_method = non_normal_method, link = link,
	               linear_solver = linear_solver, multistart = multistart)
	if (bootstrap > 0L) {
		fit <- .ae_attach_bootstrap(fit, bootstrap, bootstrap_type,
		                            bootstrap_block_length, bootstrap_seed,
		                            verbose)
	}
	fit
}


#' Fast (MCMC-free) AME estimation for a longitudinal network
#'
#' @description
#' Fits an additive and multiplicative effects (AME) model to a longitudinal
#' (replicated) network by \strong{iterative block coordinate descent},
#' producing a fast point estimate with no MCMC and no credible intervals.
#' This is the longitudinal counterpart of \code{\link{ame_als}}; the
#' estimated effects (\code{mu}, \code{beta}, \code{a}, \code{b}, \code{U},
#' \code{V}) are \emph{static} (pooled across time), as in a non-dynamic
#' \code{\link{lame}} fit.
#'
#' The estimation algorithm adapts the iterative block coordinate descent
#' estimator of the Social Influence Regression model of Hoff & Minhas
#' (\code{sir::sir_alsfit()}) to the AME model; it is a port and adaptation,
#' not original \pkg{lame} methodology. See \code{\link{ame_als}} for the
#' algorithm details.
#'
#' @param Y a list of \code{T} relational matrices, or a 3D array
#'   \code{[n_row, n_col, T]}. All slices must share the same dimensions and,
#'   if the matrices carry row/column names, the same actors in the same
#'   order (the ALS estimator does not handle changing actor
#'   compositions; use \code{\link{lame}} for that).
#' @param Xdyad a list of \code{T} dyadic covariate matrices/arrays, or \code{NULL}.
#' @param Xrow a list of \code{T} row/sender covariate matrices, or \code{NULL}.
#' @param Xcol a list of \code{T} column/receiver covariate matrices, or \code{NULL}.
#' @inheritParams ame_als
#'
#' @return An object of class \code{"ame_als"}; see \code{\link{ame_als}}.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The iterative block
#' coordinate descent estimator adapted here originates with that work
#' (implemented in \code{sir::sir_alsfit()}).
#'
#' @seealso \code{\link{ame_als}}, \code{\link{ame_als_bootstrap}},
#'   \code{\link{lame}} for the full MCMC estimator,
#'   \code{\link{als_dynamic_beta}} for a regression-only smoother
#'   that estimates \emph{only} a time-varying \eqn{\beta_t} (no
#'   \eqn{a, b, U, V}; not a special case of this function).
#'
#' @examples
#' Y <- replicate(4, { m <- matrix(rnorm(400), 20, 20); diag(m) <- NA; m },
#'                simplify = FALSE)
#' fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' coef(fit)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
lame_als <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
                     R = 0, family = "normal",
                     mode = c("unipartite", "bipartite"),
                     symmetric = FALSE, max_iter = 200, tol = 1e-6,
                     lowrank_method = c("mm", "als", "hybrid"),
                     non_normal_method = c("irls", "transform"),
                     link = c("probit", "logit"),
                     linear_solver = c("eigen", "qr", "auto"),
                     multistart = c("none", "cheap", "full"),
                     bootstrap = 0L,
                     bootstrap_type = c("parametric", "block"),
                     bootstrap_block_length = 1L,
                     bootstrap_seed = NULL,
                     verbose = TRUE, seed = 6886) {
	mode <- match.arg(mode)
	lowrank_method <- match.arg(lowrank_method)
	non_normal_method <- match.arg(non_normal_method)
	link <- match.arg(link)
	linear_solver <- match.arg(linear_solver)
	multistart <- match.arg(multistart)
	bootstrap_type <- match.arg(bootstrap_type)
	if (!is.numeric(bootstrap) || length(bootstrap) != 1L || bootstrap < 0L ||
	    !isTRUE(all.equal(bootstrap, round(bootstrap)))) {
		cli::cli_abort("{.arg bootstrap} must be a non-negative integer (0 = no bootstrap).")
	}
	bootstrap <- as.integer(bootstrap)
	fit <- .ae_run(Y, Xdyad, Xrow, Xcol, R, family, mode, symmetric,
	               max_iter, tol, verbose, seed, longitudinal = TRUE,
	               call = match.call(), lowrank_method = lowrank_method,
	               non_normal_method = non_normal_method, link = link,
	               linear_solver = linear_solver, multistart = multistart)
	if (bootstrap > 0L) {
		fit <- .ae_attach_bootstrap(fit, bootstrap, bootstrap_type,
		                            bootstrap_block_length, bootstrap_seed,
		                            verbose)
	}
	fit
}

# --------------------------------------------------------------------------
# S3 methods for ame_als objects
# --------------------------------------------------------------------------

#' Extract coefficients from a fast AME fit
#'
#' Returns the point-estimated regression coefficients (intercept first) of an
#' \code{ame_als} fit. There is no posterior here; for uncertainty use
#' \code{\link{ame_als_bootstrap}}.
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return A named numeric vector of coefficients.
#' @method coef ame_als
#' @export
coef.ame_als <- function(object, ...) {
	object$coefficients
}

#' Extract fitted values from a fast AME fit
#'
#' Returns response-scale fitted values: a single matrix for a cross-sectional
#' (\code{ame_als}) fit, or a list of matrices for a longitudinal
#' (\code{lame_als}) fit.
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return A matrix or list of matrices of fitted values.
#' @method fitted ame_als
#' @export
fitted.ame_als <- function(object, ...) {
	if (isTRUE(object$longitudinal)) object$fitted else object$fitted[[1]]
}

# inverse link: map a linear predictor EZ to the response scale
.ae_link_inverse <- function(ez, family, link) {
	if (!is.null(link)) {
		switch(paste(family, link),
			"poisson log"   = { m <- exp(ez); m[m > 1e6] <- 1e6; m },
			"binary logit"  = stats::plogis(ez),
			"binary probit" = stats::pnorm(ez),
			ez)
	} else switch(family,
		normal  = ez,
		poisson = { m <- exp(ez) - 1; m[m < 0] <- 0; m[m > 1e6] <- 1e6; m },
		binary  = stats::pnorm(ez),
		ez)
}

#' Residuals from a fast AME fit
#'
#' @description
#' Returns residuals from an \code{\link{ame_als}} or
#' \code{\link{lame_als}} fit.
#'
#' @details
#' \code{type = "response"} (default) returns observed \code{Y} minus the
#' response-scale fitted values, so \code{residuals()} reconciles with
#' \code{fitted()}. \code{type = "working"} returns the Gaussian working-scale
#' residuals the estimator actually minimised (identical to \code{"response"}
#' for the \code{normal} family).
#'
#' @param object an \code{ame_als} object.
#' @param type \code{"response"} (default) or \code{"working"}.
#' @param ... ignored.
#' @return A matrix (cross-sectional) or list of matrices (longitudinal).
#' @method residuals ame_als
#' @export
residuals.ame_als <- function(object, type = c("response", "working"), ...) {
	type <- match.arg(type)
	res <- if (type == "working") {
		object$residuals
	} else {
		lapply(seq_len(object$n_time), function(t) {
			r <- object$Y[, , t] - object$fitted[[t]]
			dimnames(r) <- dimnames(object$fitted[[t]])
			r
		})
	}
	if (isTRUE(object$longitudinal)) res else res[[1]]
}

#' Predictions from a fast AME fit
#'
#' @description
#' Returns predictions from an \code{\link{ame_als}} or
#' \code{\link{lame_als}} fit: the fitted linear predictor on the link
#' scale (\code{type = "link"}) or the response scale (\code{type = "response"},
#' the default).
#'
#' @details
#' With \code{newdata = NULL} the training-data fitted values are returned.
#' Supplying \code{newdata} -- a dyadic covariate array with the \emph{same}
#' actors and dimensions as the fitted design -- substitutes the dyadic
#' covariate contribution while holding the intercept, additive effects and
#' multiplicative term fixed; it predicts counterfactual dyadic covariates for
#' the same network, not out-of-sample actors.
#'
#' @param object an \code{ame_als} object.
#' @param newdata optional dyadic covariate array matching \code{object$X}
#'   (\code{[n_row, n_col, p, T]}); a matrix or 3D array is coerced.
#' @param type \code{"response"} (default) or \code{"link"}.
#' @param ... ignored.
#' @return A matrix (cross-sectional) or list of matrices (longitudinal).
#' @method predict ame_als
#' @export
predict.ame_als <- function(object, newdata = NULL,
                                type = c("response", "link"), ...) {
	type <- match.arg(type)
	# catch MCMC-only predict args and abort with a clear message. `h`
	# absorbed via `...` would otherwise be ignored on an ALS fit (which
	# has no dynamic state to propagate) and silently return fitted values.
	dots <- list(...)
	dot_names <- names(dots)
	if ("h" %in% dot_names && isTRUE(any(unlist(dots[dot_names == "h"]) > 0))) {
		cli::cli_abort(c(
			"{.fn predict.ame_als} does not support h-step forecasting.",
			"i" = "ALS is a point estimator with no dynamic state to propagate.",
			"i" = "Refit with {.code method = \"mcmc\", dynamic_beta = TRUE} (cross-sectional) or {.fn lame} (longitudinal) for forecasts."))
	}
	if (any(c("interval", "by_draw", "newexposure") %in% dot_names)) {
		cli::cli_warn(c(
			"Argument{?s} {.arg {intersect(c('interval','by_draw','newexposure'), dot_names)}} only apply to MCMC fits and are ignored here.",
			"i" = "For uncertainty on an ALS fit, run {.fn ame_als_bootstrap} or pass {.arg bootstrap} to {.fn ame_als}."))
	}
	EZ <- object$EZ
	if (!is.null(newdata)) {
		X_old  <- object$X
		p_dyad <- dim(X_old)[3]
		if (p_dyad == 0L) {
			cli::cli_abort("This fit has no dyadic covariates; {.arg newdata} cannot be applied.")
		}
		Xn <- newdata
		if (is.matrix(Xn)) Xn <- array(Xn, c(dim(Xn), 1L, 1L))
		if (length(dim(Xn)) == 3L) Xn <- array(Xn, c(dim(Xn), 1L))
		if (!identical(dim(Xn), dim(X_old))) {
			cli::cli_abort(c(
				"{.arg newdata} must match the fitted dyadic design dimensions.",
				"i" = "Expected [{paste(dim(X_old), collapse = ', ')}], got [{paste(dim(Xn), collapse = ', ')}]."))
		}
		beta_dyad <- object$beta[seq_len(p_dyad)]
		EZ <- lapply(seq_len(object$n_time), function(t) {
			xb_old <- xb_new <- matrix(0, object$dims$n_row, object$dims$n_col)
			for (k in seq_len(p_dyad)) {
				xb_old <- xb_old + beta_dyad[k] * X_old[, , k, t]
				xb_new <- xb_new + beta_dyad[k] * Xn[, , k, t]
			}
			ez <- object$EZ[[t]] - xb_old + xb_new
			dimnames(ez) <- dimnames(object$EZ[[t]])
			ez
		})
	}
	pred <- if (type == "link") EZ else {
		lapply(EZ, function(ez) .ae_link_inverse(ez, object$family, object$link))
	}
	if (isTRUE(object$longitudinal)) pred else pred[[1]]
}

#' Plot the convergence of a fast AME fit
#'
#' @description
#' Plots the block coordinate descent deviance/SSE history -- a convergence
#' diagnostic for an \code{\link{ame_als}} or \code{\link{lame_als}}
#' fit. For substantive plots see \code{\link{uv_plot}}, \code{\link{ab_plot}}
#' and \code{\link{gof_plot}}.
#'
#' @param x an \code{ame_als} object.
#' @param ... further arguments passed to the underlying plot.
#' @return \code{x}, invisibly.
#' @method plot ame_als
#' @export
plot.ame_als <- function(x, ...) {
	dh <- x$dev_history
	if (is.null(dh) || length(dh) < 2L) {
		cli::cli_abort(c(
			"No convergence history is stored for this fit.",
			"i" = "Use {.fn uv_plot}, {.fn ab_plot} or {.fn gof_plot} to visualise the fit."))
	}
	graphics::plot.default(seq_along(dh), dh, type = "b", pch = 19,
		xlab = "BCD iteration", ylab = "deviance / SSE",
		main = "ame_als convergence", ...)
	invisible(x)
}

#' Log-likelihood is not defined for a fast AME fit
#'
#' The fast estimator minimises a (working-response) least-squares objective,
#' not a family likelihood, so it has no log-likelihood and \code{AIC} / \code{BIC}
#' are undefined. Calling \code{logLik()} raises an informative error rather
#' than returning a misleading number.
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return Never returns; raises an error.
#' @method logLik ame_als
#' @export
logLik.ame_als <- function(object, ...) {
	cli::cli_abort(c(
		"A log-likelihood is not defined for the fast {.fn ame_als} estimator.",
		"i" = "It minimises a working-response least-squares objective, so {.fn AIC} / {.fn BIC} do not apply.",
		"i" = "For likelihood-based model comparison use the MCMC estimator {.fn ame} / {.fn lame}."))
}

#' Print an ame_als object
#'
#' @param x an \code{ame_als} object.
#' @param digits number of digits to display.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @method print ame_als
#' @export
print.ame_als <- function(x, digits = 4, ...) {
	cli::cli_text("{.strong AME fit via iterative block coordinate descent} (fast, MCMC-free)")
	cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}}{ifelse(x$symmetric, ' (symmetric)', '')} | R = {x$R}")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slice{?s}")
	cli::cli_text("Converged: {.val {x$converged}} in {x$iterations} iteration{?s}")
	# loud notice: ALS is always static, even on a longitudinal panel
	if (!is.null(x$n_time) && x$n_time > 1L) {
		cli::cli_alert_warning(c(
			"This is a {.strong STATIC} fit pooled across {x$n_time} time periods: ",
			"{.code U}, {.code V}, {.code a}, {.code b} do not vary by time. ",
			"The ALS estimator has no dynamic option; use {.fn lame} ",
			"(MCMC) with {.code dynamic_uv = TRUE} / {.code dynamic_ab = TRUE} ",
			"for time-varying latent positions."))
	}
	cli::cli_text("")
	cli::cli_text("{.strong Coefficients:}")
	print(round(x$coefficients, digits))
	cli::cli_text("")
	cli::cli_text("{.strong Variance components:}")
	print(round(x$VC, digits))
	cli::cli_text("")
	cli::cli_text(cli::col_grey("va/vb: variances of the sender/receiver effects | cab: their covariance | rho: dyadic residual reciprocity | ve: residual variance."))
	cli::cli_text(cli::col_grey("(Descriptive summaries, not random-effect variance components.)"))
	if (isTRUE(x$meta$uncertainty_available)) {
		nm <- x$meta$bootstrap_type %||% "parametric"
		Rb <- x$meta$bootstrap_R    %||% NA_integer_
		cli::cli_text(cli::col_grey(
			"Bootstrap intervals available ({nm}, R = {Rb}); use {.fn confint} for CIs."))
	} else {
		cli::cli_text(cli::col_grey("Point estimate only; no uncertainty. Re-fit with {.code bootstrap = N} (e.g. {.code ame_als(Y, ..., bootstrap = 200)}) to attach bootstrap intervals in one call."))
	}
	invisible(x)
}

#' Summarize an ame_als object
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return A list of class \code{"summary.ame_als"}, printed as a table.
#' @method summary ame_als
#' @export
summary.ame_als <- function(object, ...) {
	# Build the regression-coefficient table. When a bootstrap is
	# attached (i.e. uncertainty is available), merge the per-coefficient
	# SE / CI columns into the same table users see at the top of the
	# summary printout — without this, users have to call confint()
	# separately and miss the obvious join.
	bs <- object$bootstrap
	have_bs <- isTRUE(object$meta$uncertainty_available) &&
		!is.null(bs) && !is.null(bs$param_names) &&
		!is.null(bs$point_est)
	if (have_bs) {
		# `bs$param_names` lines up with bs$point_est / bs$se / bs$ci_lo /
		# bs$ci_hi. Match these to object$coefficients by name.
		coef_nms <- names(object$coefficients)
		match_idx <- match(coef_nms, bs$param_names)
		tab <- cbind(
			Estimate  = object$coefficients,
			Std.Error = bs$se[match_idx],
			`CI 2.5%`  = bs$ci_lo[match_idx],
			`CI 97.5%` = bs$ci_hi[match_idx])
		rownames(tab) <- coef_nms
	} else {
		tab <- cbind(Estimate = object$coefficients)
	}
	out <- list(
		call = object$call,
		family = object$family, mode = object$mode,
		symmetric = object$symmetric, R = object$R,
		n_time = object$n_time, dims = object$dims,
		coefficients = tab, VC = object$VC,
		converged = object$converged, iterations = object$iterations,
		deviance = object$deviance,
		approximation_note = object$meta$approximation_note,
		uncertainty_available = isTRUE(object$meta$uncertainty_available),
		uncertainty_source    = object$meta$uncertainty_source,
		bootstrap_R           = object$meta$bootstrap_R,
		bootstrap_type        = object$meta$bootstrap_type)
	class(out) <- "summary.ame_als"
	out
}

#' Print a fast AME summary
#'
#' @param x a \code{summary.ame_als} object.
#' @param digits number of digits to display.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @method print summary.ame_als
#' @export
print.summary.ame_als <- function(x, digits = 4, ...) {
	cli::cli_h2("AME exploration summary")
	cli::cli_text("Algorithm: iterative block coordinate descent (Hoff & Minhas SIR estimator)")
	cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}}{ifelse(x$symmetric, ' (symmetric)', '')} | R = {x$R}")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slice{?s}, {x$dims$p} covariate{?s}")
	cli::cli_text("Converged: {.val {x$converged}} ({x$iterations} iteration{?s}); deviance/SSE = {round(x$deviance, digits)}")
	if (!is.null(x$n_time) && x$n_time > 1L) {
		cli::cli_alert_warning(c(
			"This is a {.strong STATIC} fit pooled across {x$n_time} time periods: ",
			"{.code U}, {.code V}, {.code a}, {.code b} do not vary by time. ",
			"ALS has no dynamic option; use {.fn lame} (MCMC) with ",
			"{.code dynamic_uv = TRUE} / {.code dynamic_ab = TRUE} for time-varying effects."))
	}
	cli::cli_rule()
	cli::cli_text("{.strong Regression coefficients}")
	print(round(x$coefficients, digits))
	cli::cli_text("")
	cli::cli_text("{.strong Variance components}")
	print(round(x$VC, digits))
	cli::cli_text(cli::col_grey("va/vb: variances of the sender/receiver effects | cab: their covariance | rho: dyadic residual reciprocity | ve: residual variance."))
	cli::cli_text(cli::col_grey("(Descriptive summaries, not random-effect variance components.)"))
	cli::cli_rule()
	cli::cli_text(cli::col_grey(x$approximation_note))
	if (isTRUE(x$uncertainty_available)) {
		nm <- x$uncertainty_source %||% "bootstrap"
		Rb <- x$bootstrap_R %||% NA_integer_
		tp <- x$bootstrap_type %||% "parametric"
		cli::cli_text(cli::col_grey(
			"Uncertainty available: {nm} ({tp}, R = {Rb}). Use {.fn confint} for intervals."))
	} else {
		cli::cli_text(cli::col_grey("No standard errors: this is a point estimate. Re-fit with {.code bootstrap = N} (e.g. {.code ame_als(Y, ..., bootstrap = 200)}) to attach bootstrap intervals."))
	}
	invisible(x)
}
