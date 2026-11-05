# dynamic ALS/MAP point estimator for longitudinal AME models
#
# covers normal, binary and Poisson panels with optional changing composition,
# dynamic additive effects, selected dynamic regression coefficients, static
# low-rank UV when R > 0, smooth dynamic UV, bipartite dynamic G, and
# node-covariate coefficients identified by the ALS orthogonality
# decomposition. dynamic node coefficients use period-specific node values;
# static ones use the actor-level summaries of static ALS. snap-shift UV has no
# ALS route (see .lame_snap_als_refuse in lame_snap_als.R).

.lda_default <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.lda_max_abs <- function(x) {
	if (is.null(x) || length(x) == 0L) return(0)
	out <- max(abs(x), na.rm = TRUE)
	if (is.finite(out)) out else 0
}

.lda_presence_row_means <- function(mat, presence = NULL) {
	if (is.null(presence) || !all(dim(presence) == dim(mat))) {
		return(rowMeans(mat))
	}
	out <- numeric(nrow(mat))
	for (i in seq_len(nrow(mat))) {
		ok <- presence[i, ] & is.finite(mat[i, ])
		out[i] <- if (any(ok)) mean(mat[i, ok]) else 0
	}
	out
}

.lda_presence_col_vars <- function(mat, presence = NULL) {
	if (is.null(presence) || !all(dim(presence) == dim(mat))) {
		return(apply(mat, 2, stats::var))
	}
	vapply(seq_len(ncol(mat)), function(t) {
		ok <- presence[, t] & is.finite(mat[, t])
		if (sum(ok) > 1L) stats::var(mat[ok, t]) else NA_real_
	}, numeric(1))
}

.lda_presence_col_cov <- function(
	a_mat,
	b_mat,
	row_presence = NULL,
	col_presence = NULL
) {
	if (is.null(row_presence) || is.null(col_presence) ||
			!all(dim(row_presence) == dim(a_mat)) ||
			!all(dim(col_presence) == dim(b_mat)) ||
			nrow(a_mat) != nrow(b_mat)) {
		return(vapply(seq_len(ncol(a_mat)), function(t) {
			stats::cov(a_mat[, t], b_mat[, t])
		}, numeric(1)))
	}
	vapply(seq_len(ncol(a_mat)), function(t) {
		ok <- row_presence[, t] & col_presence[, t] &
			is.finite(a_mat[, t]) & is.finite(b_mat[, t])
		if (sum(ok) > 1L) stats::cov(a_mat[ok, t], b_mat[ok, t]) else NA_real_
	}, numeric(1))
}

.lda_transition_from_presence <- function(presence) {
	if (is.null(presence) || ncol(presence) < 2L) return(NULL)
	out <- matrix(
		1,
		nrow(presence),
		ncol(presence),
		dimnames = dimnames(presence)
	)
	p <- presence
	p[is.na(p)] <- FALSE
	out[, -1L] <- as.numeric(
		p[, -1L, drop = FALSE] &
			p[, -ncol(p), drop = FALSE]
	)
	out
}

.lda_combine_transition_weights <- function(weights, presence_weights) {
	if (is.null(presence_weights)) return(weights)
	if (is.null(weights)) return(presence_weights)
	if (!all(dim(weights) == dim(presence_weights))) return(weights)
	weights * presence_weights
}

.lda_clip_rho <- function(x, default = 0.95) {
	if (is.null(x) || length(x) == 0L || !is.finite(x[1L])) return(default)
	pmin(pmax(as.numeric(x[1L]), 0), 0.999)
}

.lda_path_penalty <- function(Tt, lambda, kind = "rw1", rho = 1) {
	K <- matrix(0, Tt, Tt)
	if (lambda <= 0 || Tt < 2L) return(K)
	if (identical(kind, "rw2")) {
		if (Tt < 3L) {
			cli::cli_abort("{.arg dynamic_beta_kind = \"rw2\"} requires at least 3 time periods.")
		}
		D <- matrix(0, Tt - 2L, Tt)
		for (t in seq_len(Tt - 2L)) {
			D[t, t] <- 1
			D[t, t + 1L] <- -2
			D[t, t + 2L] <- 1
		}
		return(lambda * crossprod(D))
	}
	D <- matrix(0, Tt - 1L, Tt)
	for (t in seq_len(Tt - 1L)) {
		D[t, t] <- if (identical(kind, "ar1")) -rho else -1
		D[t, t + 1L] <- 1
	}
	lambda * crossprod(D)
}

.lda_path_penalty_weighted <- function(Tt, lambda, weights = NULL, rho = 1) {
	K <- matrix(0, Tt, Tt)
	if (lambda <= 0 || Tt < 2L) return(K)
	w <- if (is.null(weights)) rep(1, Tt) else as.numeric(weights)
	if (length(w) != Tt) w <- rep(1, Tt)
	w[!is.finite(w) | w < 0] <- 0
	D <- matrix(0, Tt - 1L, Tt)
	for (t in seq_len(Tt - 1L)) {
		D[t, t] <- -rho
		D[t, t + 1L] <- 1
	}
	lambda * crossprod(sqrt(w[-1L]) * D)
}

# effective df of a bipartite fit's position paths: trace(H^-1 K) of the joint
# smoother of the (u, v) block, K the data information of the paths and H the
# Hessian of the penalized objective. the joint trace keeps the residual
# curvature of the bilinear term (w_ij r_ij g in d2 f / du_i dv_j) that
# conditional traces miss. H is block tridiagonal in time, so the diagonal
# blocks of its inverse come from one forward and one backward sweep
.lda_bip_joint_trace <- function(
	U,
	V,
	G,
	resid,
	obs_mask,
	wgt,
	lambda_uv,
	rho_uv,
	transition_u = NULL,
	transition_v = NULL,
	level_prec_g = 1,
	details = FALSE,
	row_unit = 1
) {
	nr <- dim(U)[1L]
	nc <- dim(V)[1L]
	RA <- dim(U)[2L]
	RB <- dim(V)[2L]
	Tt <- dim(U)[3L]
	m <- nr * RA + nc * RB
	ju <- seq_len(nr * RA)
	jv <- nr * RA + seq_len(nc * RB)
	pen_u <- lapply(seq_len(nr), function(i) {
		.lda_path_penalty_weighted(
			Tt,
			lambda_uv,
			if (is.null(transition_u)) rep(1, Tt) else transition_u[i, ],
			rho_uv
		)
	})
	pen_v <- lapply(seq_len(nc), function(j) {
		.lda_path_penalty_weighted(
			Tt,
			lambda_uv,
			if (is.null(transition_v)) rep(1, Tt) else transition_v[j, ],
			rho_uv
		)
	})
	pdiag <- function(t, s) {
		c(
			rep(vapply(pen_u, function(P) P[t, s], numeric(1)), each = RA),
			rep(vapply(pen_v, function(P) P[t, s], numeric(1)), each = RB)
		)
	}
	ng <- RA * RB
	KG <- matrix(0, ng, ng)
	KXG <- BXG <- vector("list", Tt)
	K <- vector("list", Tt)
	A <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		Ut <- matrix(U[, , t], nr, RA)
		Vt <- matrix(V[, , t], nc, RB)
		GV <- Vt %*% t(G)
		GU <- Ut %*% G
		w <- wgt[, , t]
		w[!obs_mask[, , t] | !is.finite(w) | w < 0] <- 0
		r <- resid[[t]]
		r[!is.finite(r)] <- 0
		Kt <- matrix(0, m, m)
		Ct <- matrix(0, m, m)
		for (i in seq_len(nr)) {
			rows <- (i - 1L) * RA + seq_len(RA)
			Kt[rows, rows] <- crossprod(GV * sqrt(w[i, ]))
			Kt[rows, jv] <- kronecker(t(GV * w[i, ]), matrix(GU[i, ], 1L, RB))
			Ct[rows, jv] <- kronecker(matrix(w[i, ] * r[i, ], 1L, nc), G)
		}
		for (j in seq_len(nc)) {
			cols <- nr * RA + (j - 1L) * RB + seq_len(RB)
			Kt[cols, cols] <- crossprod(GU * sqrt(w[, j]))
		}
		Kt[jv, ju] <- t(Kt[ju, jv])
		Ct[jv, ju] <- t(Ct[ju, jv])
		kxg <- cxg <- matrix(0, m, ng)
		for (i in seq_len(nr)) {
			rows <- (i - 1L) * RA + seq_len(RA)
			xg <- kronecker(Vt, matrix(Ut[i, ], 1L, RA))
			kxg[rows, ] <- crossprod(GV * sqrt(w[i, ]), xg * sqrt(w[i, ]))
			KG <- KG + crossprod(xg * sqrt(w[i, ]))
			mixed <- colSums(Vt * (w[i, ] * r[i, ]))
			for (b in seq_len(RB)) cxg[rows, (b - 1L) * RA + seq_len(RA)] <- diag(mixed[b], RA)
		}
		for (j in seq_len(nc)) {
			rows <- nr * RA + (j - 1L) * RB + seq_len(RB)
			xg <- kronecker(matrix(Vt[j, ], 1L, RB), Ut)
			kxg[rows, ] <- crossprod(GU * sqrt(w[, j]), xg * sqrt(w[, j]))
			mixed <- colSums(Ut * (w[, j] * r[, j]))
			for (b in seq_len(RB)) cxg[rows[b], (b - 1L) * RA + seq_len(RA)] <- mixed
		}
		KXG[[t]] <- kxg
		BXG[[t]] <- kxg - cxg
		K[[t]] <- Kt
		A[[t]] <- Kt - Ct + diag(pdiag(t, t) + (.lda_level_prec + 1e-6) * row_unit, m)
	}
	inv <- function(M) tryCatch(solve(M), error = function(e) .lda_safe_pinv(M))
	Sinv <- vector("list", Tt)
	B <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		S <- A[[t]]
		if (t > 1L) {
			B[[t - 1L]] <- pdiag(t - 1L, t)
			S <- S - B[[t - 1L]] * t(B[[t - 1L]] * Sinv[[t - 1L]])
		}
		Sinv[[t]] <- inv(S)
	}
	D <- Sinv[[Tt]]
	total <- sum(D * K[[Tt]])
	if (Tt > 1L) for (t in (Tt - 1L):1L) {
		SB <- Sinv[[t]] * rep(B[[t]], each = m)
		D <- Sinv[[t]] + SB %*% D %*% t(SB)
		total <- total + sum(D * K[[t]])
	}
	# eliminate the temporal paths, then solve the small static interaction block
	z <- x <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		rhs <- BXG[[t]]
		if (t > 1L) rhs <- rhs - B[[t - 1L]] * z[[t - 1L]]
		z[[t]] <- Sinv[[t]] %*% rhs
	}
	x[[Tt]] <- z[[Tt]]
	if (Tt > 1L) for (t in (Tt - 1L):1L) {
		x[[t]] <- z[[t]] - Sinv[[t]] %*% (B[[t]] * x[[t + 1L]])
	}
	schur <- KG + diag(level_prec_g, ng)
	information <- KG
	for (t in seq_len(Tt)) {
		schur <- schur - crossprod(BXG[[t]], x[[t]])
		information <- information + crossprod(x[[t]], K[[t]] %*% x[[t]]) -
			crossprod(KXG[[t]], x[[t]]) - crossprod(x[[t]], KXG[[t]])
	}
	answer <- total + sum(inv(schur) * information)
	if (details) return(list(
		trace = answer,
		conditional = total,
		A = A,
		K = K,
		B = B,
		KXG = KXG,
		BXG = BXG,
		KG = KG,
		schur = schur
	))
	answer
}

.lda_path_hessian <- function(
	obs_prec,
	lambda,
	kind = "rw1",
	rho = 1,
	transition_weights = NULL,
	level_prec = 0,
	ridge = 1e-8
) {
	Tt <- length(obs_prec)
	P <- if (!is.null(transition_weights) && !identical(kind, "rw2")) {
		.lda_path_penalty_weighted(
			Tt,
			lambda,
			transition_weights,
			rho = if (identical(kind, "ar1")) rho else 1
		)
	} else {
		.lda_path_penalty(Tt, lambda, kind, rho)
	}
	diag(obs_prec, nrow = Tt) + P + diag(ridge + level_prec, Tt)
}

# eigen pseudo-inverse on the tolerance .ae_safe_solve() uses, so a path
# sub-problem with an unobserved period is handled the same alone or inside the
# constrained block
.lda_safe_pinv <- function(A) {
	A <- (A + t(A)) / 2
	eg <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
	if (is.null(eg)) return(matrix(0, nrow(A), ncol(A)))
	lam <- eg$values
	thr <- max(lam, 0) * sqrt(.Machine$double.eps)
	keep <- is.finite(lam) & lam > thr
	if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
	Qk <- eg$vectors[, keep, drop = FALSE]
	Qk %*% (t(Qk) / lam[keep])
}

# exact minimiser of separable penalized path sub-problems subject to every
# period's mean being zero: the level of a_t (or b_t) is not identified against
# the intercept, and shared Lagrange multipliers keep this an exact block
# minimiser (projecting the mean out afterwards would not be a descent step)
.lda_solve_paths_centred <- function(
	obs_prec,
	rhs,
	lambda,
	kind = "ar1",
	rho = 1,
	transition_weights = NULL,
	level_prec = 0,
	ridge = 1e-8
) {
	n <- nrow(rhs)
	Tt <- ncol(rhs)
	out <- matrix(0, n, Tt)
	if (n == 0L) return(out)
	# an actor never observed and without a level prior has a free path level: it
	# is held at zero and left out of the gauge, as in the static additive block
	active <- if (level_prec > 0) rep(TRUE, n) else rowSums(obs_prec) > 0
	if (!any(active)) return(out)
	Kinv <- vector("list", n)
	S <- matrix(0, Tt, Tt)
	m <- numeric(Tt)
	for (i in which(active)) {
		tw <- if (is.null(transition_weights)) NULL else transition_weights[i, ]
		Ki <- .lda_safe_pinv(.lda_path_hessian(
			obs_prec[i, ],
			lambda,
			kind,
			rho,
			tw,
			level_prec,
			ridge
		))
		Kinv[[i]] <- Ki
		S <- S + Ki
		m <- m + as.numeric(Ki %*% rhs[i, ])
	}
	mult <- -.ae_safe_solve(S, m)
	for (i in which(active)) {
		out[i, ] <- as.numeric(Kinv[[i]] %*% (rhs[i, ] + mult))
	}
	out
}

.lda_unit_weights <- function(Z) {
	W <- array(0, dim = dim(Z))
	W[is.finite(Z)] <- 1
	W
}

.lda_eta_array <- function(E) {
	array(unlist(E, use.names = FALSE), dim = c(nrow(E[[1L]]), ncol(E[[1L]]), length(E)))
}

.lda_irls_work <- function(Yarr, eta_arr, family, link) {
	work_fun <- switch(
		paste(family, link),
		"poisson log" = function(eta) .ae_irls_poisson(Yarr, eta),
		"binary probit" = function(eta) .ae_irls_binary_probit(Yarr, eta),
		"binary logit" = function(eta) .ae_irls_binary_logit(Yarr, eta),
		cli::cli_abort("IRLS is not available for family {.val {family}} with link {.val {link}}.")
	)
	wk <- work_fun(eta_arr)
	Z <- wk$z
	W <- wk$w
	Z[!is.finite(Yarr)] <- NA_real_
	W[!is.finite(Yarr)] <- NA_real_
	W[!is.finite(W) | W < 0] <- 0
	list(Z = Z, W = W)
}

.lda_o_array <- function(O, Tt) {
	if (is.null(O)) return(NULL)
	d <- dim(O)
	if (length(d) == 3L) return(O)
	array(rep(O, Tt), dim = c(nrow(O), ncol(O), Tt))
}

.lda_o_slice <- function(Oarr, t) {
	if (length(dim(Oarr)) == 3L) Oarr[, , t] else Oarr
}

.lda_symmetric_factor_slice <- function(
	O,
	R,
	row_names = NULL,
	dim_names = NULL
) {
	O <- (O + t(O)) / 2
	n <- nrow(O)
	if (R <= 0L) {
		U <- V <- matrix(0, n, 0L)
		return(list(U = U, V = V, O = matrix(0, n, n)))
	}
	eg <- eigen(O, symmetric = TRUE)
	rr <- min(R, length(eg$values))
	idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(rr)]
	vals <- eg$values[idx]
	scale <- sqrt(pmax(abs(vals), 0))
	U <- sweep(eg$vectors[, idx, drop = FALSE], 2, scale, "*")
	V <- sweep(U, 2, sign(vals), "*")
	if (rr < R) {
		U <- cbind(U, matrix(0, n, R - rr))
		V <- cbind(V, matrix(0, n, R - rr))
	}
	U <- U[, seq_len(R), drop = FALSE]
	V <- V[, seq_len(R), drop = FALSE]
	if (is.null(dim_names)) dim_names <- paste0("dim", seq_len(R))
	dimnames(U) <- dimnames(V) <- list(row_names, dim_names)
	O_rank <- tcrossprod(U, V)
	dimnames(O_rank) <- list(row_names, row_names)
	list(U = U, V = V, O = O_rank)
}

.lda_symmetric_o_to_uv_cube <- function(
	Oarr,
	R,
	row_names = NULL,
	dim_names = NULL,
	time_names = NULL
) {
	if (length(dim(Oarr)) != 3L) Oarr <- .lda_o_array(Oarr, 1L)
	d <- dim(Oarr)
	if (is.null(row_names)) row_names <- dimnames(Oarr)[[1L]]
	if (is.null(time_names)) time_names <- dimnames(Oarr)[[3L]]
	if (is.null(dim_names)) dim_names <- paste0("dim", seq_len(R))
	U <- array(
		0,
		dim = c(d[1L], R, d[3L]),
		dimnames = list(row_names, dim_names, time_names)
	)
	V <- array(
		0,
		dim = c(d[1L], R, d[3L]),
		dimnames = list(row_names, dim_names, time_names)
	)
	O_rank <- array(
		0,
		dim = d,
		dimnames = list(row_names, row_names, time_names)
	)
	for (t in seq_len(d[3L])) {
		f <- .lda_symmetric_factor_slice(
			Oarr[, , t],
			R,
			row_names,
			dim_names
		)
		U[, , t] <- f$U
		V[, , t] <- f[["V"]]
		O_rank[, , t] <- f$O
	}
	list(U = U, V = V, Oarr = O_rank)
}

.lda_rect_identity <- function(nr, nc) {
	out <- matrix(0, nr, nc)
	if (nr > 0L && nc > 0L) {
		for (k in seq_len(min(nr, nc))) out[k, k] <- 1
	}
	out
}

.lda_coerce_factor <- function(M, n, R, scale = 0.1, seed = NULL) {
	if (R <= 0L) return(matrix(0, n, 0L))
	if (!is.null(seed)) set.seed(seed)
	out <- matrix(0, n, R)
	if (!is.null(M) && length(dim(M)) == 2L && nrow(M) == n) {
		k <- min(ncol(M), R)
		if (k > 0L) out[, seq_len(k)] <- M[, seq_len(k), drop = FALSE]
		if (R > k) out[, (k + 1L):R] <- stats::rnorm(n * (R - k), 0, scale)
	} else {
		out[, ] <- stats::rnorm(n * R, 0, scale)
	}
	out
}

.lda_uv_to_o_array <- function(U, V) {
	if (is.null(U) || is.null(V)) return(NULL)
	du <- dim(U)
	dv <- dim(V)
	if (length(du) == 2L) {
		return(array(tcrossprod(U, V), dim = c(du[1L], dv[1L], 1L)))
	}
	Tt <- du[3L]
	out <- array(0, dim = c(du[1L], dv[1L], Tt))
	for (t in seq_len(Tt)) {
		Ut <- matrix(U[, , t], nrow = du[1L], ncol = du[2L])
		Vt <- matrix(V[, , t], nrow = dv[1L], ncol = dv[2L])
		out[, , t] <- tcrossprod(Ut, Vt)
	}
	dimnames(out) <- list(dimnames(U)[[1L]], dimnames(V)[[1L]], dimnames(U)[[3L]])
	out
}

.lda_uvg_to_o_array <- function(U, V, G, Tt = NULL, time_names = NULL) {
	if (is.null(U) || is.null(V) || is.null(G)) return(NULL)
	du <- dim(U)
	dv <- dim(V)
	dg <- dim(G)
	Tt <- max(
		if (length(du) == 3L) du[3L] else 1L,
		if (length(dv) == 3L) dv[3L] else 1L,
		if (length(dg) == 3L) dg[3L] else 1L,
		.lda_default(Tt, 1L)
	)
	out <- array(0, dim = c(du[1L], dv[1L], Tt))
	for (t in seq_len(Tt)) {
		Ut <- if (length(du) == 3L) {
			matrix(U[, , t], nrow = du[1L], ncol = du[2L])
		} else U
		Vt <- if (length(dv) == 3L) {
			matrix(V[, , t], nrow = dv[1L], ncol = dv[2L])
		} else V
		Gt <- if (length(dg) == 3L) {
			matrix(G[, , t], nrow = dg[1L], ncol = dg[2L])
		} else G
		out[, , t] <- Ut %*% Gt %*% t(Vt)
	}
	time_names <- .lda_default(time_names, if (length(du) == 3L) {
		dimnames(U)[[3L]]
	} else if (length(dv) == 3L) {
		dimnames(V)[[3L]]
	} else if (length(dg) == 3L) {
		dimnames(G)[[3L]]
	} else NULL)
	dimnames(out) <- list(dimnames(U)[[1L]], dimnames(V)[[1L]], time_names)
	out
}

.lda_procrustes_q <- function(A, ref) {
	if (ncol(A) == 1L) {
		sgn <- sign(sum(A * ref, na.rm = TRUE))
		if (!is.finite(sgn) || sgn == 0) sgn <- 1
		return(matrix(sgn, 1L, 1L))
	}
	M <- crossprod(A, ref)
	if (!all(is.finite(M))) return(diag(ncol(A)))
	sv <- tryCatch(svd(M), error = function(e) NULL)
	if (is.null(sv)) return(diag(ncol(A)))
	sv$u %*% t(sv$v)
}

# a common rotation of u_t and v_t leaves u_t v_t' and the fit unchanged and
# only moves the transition penalty, so aligning each period to the previous
# one is a descent step when the rotation uses the penalty's own metric; the
# per-actor transition weights therefore enter the cross-product
.lda_align_uv_cube <- function(U, V, weights_u = NULL, weights_v = NULL) {
	if (is.null(U) || is.null(V) || length(dim(U)) != 3L || dim(U)[3L] < 2L) {
		return(list(U = U, V = V))
	}
	Tt <- dim(U)[3L]
	use_w <- function(w, n) {
		if (is.null(w) || !is.matrix(w) || nrow(w) != n || ncol(w) != Tt) NULL else w
	}
	wu <- use_w(weights_u, dim(U)[1L])
	wv <- use_w(weights_v, dim(V)[1L])
	for (t in 2:Tt) {
		Ut <- matrix(U[, , t], nrow = dim(U)[1L], ncol = dim(U)[2L])
		Vt <- matrix(V[, , t], nrow = dim(V)[1L], ncol = dim(V)[2L])
		Ur <- matrix(U[, , t - 1L], nrow = dim(U)[1L], ncol = dim(U)[2L])
		Vr <- matrix(V[, , t - 1L], nrow = dim(V)[1L], ncol = dim(V)[2L])
		A <- rbind(Ut, Vt)
		ref <- rbind(Ur, Vr)
		w <- if (is.null(wu) && is.null(wv)) NULL else {
			c(
				if (is.null(wu)) rep(1, nrow(Ut)) else pmax(0, wu[, t]),
				if (is.null(wv)) rep(1, nrow(Vt)) else pmax(0, wv[, t])
			)
		}
		Q <- .lda_procrustes_q(if (is.null(w)) A else A * w, ref)
		U[, , t] <- Ut %*% Q
		V[, , t] <- Vt %*% Q
	}
	list(U = U, V = V)
}

.lda_coef_index <- function(mask, Tt) {
	p <- length(mask)
	static_idx <- which(!mask)
	dynamic_idx <- which(mask)
	n_static <- length(static_idx)
	n_dynamic <- length(dynamic_idx) * Tt
	static_map <- rep(NA_integer_, p)
	if (n_static > 0L) static_map[static_idx] <- seq_len(n_static)
	dynamic_map <- matrix(NA_integer_, p, Tt)
	pos <- n_static
	for (k in dynamic_idx) {
		dynamic_map[k, ] <- pos + seq_len(Tt)
		pos <- pos + Tt
	}
	list(
		n = n_static + n_dynamic,
		static_idx = static_idx,
		dynamic_idx = dynamic_idx,
		static_map = static_map,
		dynamic_map = dynamic_map
	)
}

.lda_coef_layout <- function(prep, include_nodes = TRUE) {
	nms <- c("intercept", prep$x_names_dyad)
	block <- c("intercept", rep("dyad", prep$p_dyad))
	if (isTRUE(include_nodes)) {
		if (!is.null(prep$W_row) && ncol(prep$W_row) > 0L) {
			nms <- c(nms, colnames(prep$W_row))
			block <- c(block, rep("row", ncol(prep$W_row)))
		}
		if (!is.null(prep$W_col) && ncol(prep$W_col) > 0L) {
			nms <- c(nms, colnames(prep$W_col))
			block <- c(block, rep("col", ncol(prep$W_col)))
		}
	}
	list(names = nms, block = block)
}

# the covariate-unit multiplier of the smoothness weight of coefficient k
.lda_beta_ms <- function(beta_dyn, k) if (is.null(beta_dyn$ms)) 1 else beta_dyn$ms[[k]]

.lda_beta_dyn_from_mask <- function(mask, coef_names, coef_block) {
	mask <- as.logical(mask)
	groups <- ifelse(mask, coef_block, "")
	names(mask) <- coef_names
	names(groups) <- coef_names
	group_names <- sort(unique(groups[groups != ""]))
	list(
		any = any(mask),
		mask = mask,
		static_idx = which(!mask),
		dynamic_idx = which(mask),
		groups = groups,
		group_names = group_names,
		n_groups = length(group_names)
	)
}

.lda_dynamic_node_names <- function(beta_dyn) {
	if (is.null(beta_dyn) || !isTRUE(beta_dyn$any)) return(character(0))
	mask <- beta_dyn$mask
	groups <- beta_dyn$groups
	if (is.null(names(mask)) || is.null(names(groups))) return(character(0))
	names(mask)[mask & groups %in% c("row", "col")]
}

.lda_node_matrix <- function(prep, margin, t = NULL, dynamic = FALSE) {
	W <- if (identical(margin, "row")) prep$W_row else prep$W_col
	W_arr <- if (identical(margin, "row")) prep$W_row_arr else prep$W_col_arr
	if (isTRUE(dynamic) && !is.null(W_arr) && !is.null(t)) {
		out <- W_arr[, , t, drop = FALSE]
		out <- matrix(out, nrow = dim(W_arr)[1L], ncol = dim(W_arr)[2L])
		dimnames(out) <- list(NULL, dimnames(W_arr)[[2L]])
		return(out)
	}
	W
}

.lda_bind_node_matrices <- function(W_row, W_col) {
	if (is.null(W_row)) return(W_col)
	if (is.null(W_col)) return(W_row)
	cbind(W_row, W_col)
}

.lda_bind_node_arrays <- function(A, B) {
	if (is.null(A)) return(B)
	if (is.null(B)) return(A)
	n <- dim(A)[1L]
	p <- dim(A)[2L] + dim(B)[2L]
	Tt <- dim(A)[3L]
	out <- array(0, dim = c(n, p, Tt))
	out[, seq_len(dim(A)[2L]), ] <- A
	out[, dim(A)[2L] + seq_len(dim(B)[2L]), ] <- B
	dimnames(out) <- list(
		NULL,
		c(dimnames(A)[[2L]], dimnames(B)[[2L]]),
		dimnames(A)[[3L]]
	)
	out
}

.lda_node_varies_names <- function(prep, margin) {
	W_arr <- if (identical(margin, "row")) prep$W_row_arr else prep$W_col_arr
	if (is.null(W_arr) || dim(W_arr)[3L] < 2L) return(character(0))
	nms <- dimnames(W_arr)[[2L]]
	varies <- logical(dim(W_arr)[2L])
	for (k in seq_along(varies)) {
		vals <- W_arr[, k, , drop = FALSE]
		vals <- matrix(vals, nrow = dim(W_arr)[1L], ncol = dim(W_arr)[3L])
		within <- mean(apply(vals, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
		varies[k] <- is.finite(within) && within > 1e-8
	}
	nms[varies]
}

.lda_design_obs <- function(prep, t, obs, coef_names, symmetric = FALSE) {
	p <- length(coef_names)
	Xt <- matrix(1, nrow = sum(obs), ncol = p)
	colnames(Xt) <- coef_names
	if (sum(obs) == 0L) return(Xt)
	ij <- which(obs, arr.ind = TRUE)
	ii <- ij[, 1L]
	jj <- ij[, 2L]
	W_node <- NULL
	if (isTRUE(symmetric)) {
		W_node <- .lda_bind_node_matrices(
			.lda_node_matrix(prep, "row", t, dynamic = TRUE),
			.lda_node_matrix(prep, "col", t, dynamic = TRUE)
		)
	}
	for (k in seq_along(coef_names)) {
		nm <- coef_names[k]
		if (identical(nm, "intercept")) next
		if (prep$p_dyad > 0L && nm %in% prep$x_names_dyad) {
			Xt[, k] <- prep$X[, , match(nm, prep$x_names_dyad), t][obs]
		} else if (!is.null(W_node) && isTRUE(symmetric) &&
				nm %in% colnames(W_node)) {
			Xt[, k] <- W_node[ii, nm] + W_node[jj, nm]
		} else {
			W_row_t <- .lda_node_matrix(prep, "row", t, dynamic = TRUE)
			W_col_t <- .lda_node_matrix(prep, "col", t, dynamic = TRUE)
			if (!is.null(W_row_t) && nm %in% colnames(W_row_t)) {
				Xt[, k] <- W_row_t[ii, nm]
			} else if (!is.null(W_col_t) && nm %in% colnames(W_col_t)) {
				Xt[, k] <- W_col_t[jj, nm]
			} else {
				Xt[, k] <- 0
			}
		}
	}
	Xt
}

.lda_update_beta_path <- function(
	Z,
	W,
	prep,
	coef_path_old,
	a_mat,
	b_mat,
	Oarr,
	beta_dyn,
	lambda_beta,
	kind,
	rho_beta,
	symmetric = FALSE
) {
	Tt <- prep$Tt
	p <- nrow(coef_path_old)
	coef_names <- rownames(coef_path_old)
	mask <- beta_dyn$mask
	idx <- .lda_coef_index(mask, Tt)
	H <- matrix(0, idx$n, idx$n)
	g <- numeric(idx$n)

	for (t in seq_len(Tt)) {
		obs <- is.finite(Z[, , t])
		if (!any(obs)) next
		target <- Z[, , t] - outer(a_mat[, t], b_mat[, t], "+") -
			.lda_o_slice(Oarr, t)
		y <- target[obs]
		wt <- W[, , t][obs]
		keep <- is.finite(wt) & wt > 0
		if (!any(keep)) next
		y <- y[keep]
		wt <- wt[keep]
		Xt <- .lda_design_obs(
			prep,
			t,
			obs,
			coef_names,
			symmetric = symmetric
		)
		Xt <- Xt[keep, , drop = FALSE]
		Xt_w <- Xt * sqrt(wt)
		XtX <- crossprod(Xt_w)
		Xty <- crossprod(Xt, wt * y)
		param_idx <- integer(p)
		for (k in seq_len(p)) {
			param_idx[k] <- if (mask[k]) idx$dynamic_map[k, t] else idx$static_map[k]
		}
		for (k in seq_len(p)) {
			g[param_idx[k]] <- g[param_idx[k]] + Xty[k]
			for (l in seq_len(p)) {
				H[param_idx[k], param_idx[l]] <-
					H[param_idx[k], param_idx[l]] + XtX[k, l]
			}
		}
	}

	if (length(idx$dynamic_idx) > 0L && lambda_beta > 0) {
		for (k in idx$dynamic_idx) {
			idk <- idx$dynamic_map[k, ]
			H[idk, idk] <- H[idk, idk] +
				.lda_path_penalty(Tt, lambda_beta * .lda_beta_ms(beta_dyn, k), kind, rho_beta)
		}
	}
	H <- H + diag(1e-8, idx$n)
	theta <- .ae_safe_solve(H, g)

	out <- matrix(0, p, Tt, dimnames = dimnames(coef_path_old))
	for (k in seq_len(p)) {
		if (mask[k]) {
			out[k, ] <- theta[idx$dynamic_map[k, ]]
		} else {
			out[k, ] <- theta[idx$static_map[k]]
		}
	}
	out
}

# refuse by name a coefficient that is a linear combination of the others over
# the observed cells, pooled over periods (see .lame_check_design_rank() in
# R/design_rank.R)
.lda_check_design_rank <- function(prep, coef_names, symmetric = FALSE) {
	nr <- prep$nr
	nc <- prep$nc
	p <- length(coef_names)
	if (p == 0L) return(invisible(NULL))
	every <- matrix(TRUE, nr, nc)
	xlist <- lapply(seq_len(prep$Tt), function(t) {
		array(
			.lda_design_obs(prep, t, every, coef_names, symmetric = symmetric),
			c(nr, nc, p),
			dimnames = list(NULL, NULL, coef_names)
		)
	})
	.lame_check_design_rank(xlist, prep$Yarr)
}

# mean square of every regression column over the observed cells (one for the
# intercept): the unit of a covariate, so that rescaling a covariate rescales
# only its coefficient path under the smoothness penalty
.lda_design_ms <- function(prep, coef_names, symmetric = FALSE) {
	tot <- numeric(length(coef_names))
	n <- 0
	for (t in seq_len(prep$Tt)) {
		ok <- is.finite(prep$Yarr[, , t])
		if (!any(ok)) next
		x <- .lda_design_obs(prep, t, ok, coef_names, symmetric = symmetric)
		tot <- tot + colSums(x^2)
		n <- n + nrow(x)
	}
	ms <- if (n > 0) tot / n else rep(1, length(coef_names))
	ms[!is.finite(ms) | ms <= 0] <- 1
	stats::setNames(ms, coef_names)
}

.lda_xbeta_array <- function(
	prep,
	coef_path,
	symmetric = FALSE,
	beta_dyn = NULL
) {
	nr <- prep$nr
	nc <- prep$nc
	Tt <- prep$Tt
	p <- nrow(coef_path)
	cn <- rownames(coef_path)
	dyn_node_names <- .lda_dynamic_node_names(beta_dyn)
	out <- array(0, dim = c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		xt <- matrix(coef_path[1L, t], nr, nc)
		if (prep$p_dyad > 0L && p > 1L) {
			kmax <- min(prep$p_dyad, p - 1L)
			for (k in seq_len(kmax)) {
				xt <- xt + coef_path[k + 1L, t] * prep$X[, , k, t]
			}
		}
		if (!is.null(cn) && isTRUE(symmetric)) {
			W_node_static <- .lda_bind_node_matrices(prep$W_row, prep$W_col)
			W_node_dyn <- .lda_bind_node_matrices(
				.lda_node_matrix(prep, "row", t, dynamic = TRUE),
				.lda_node_matrix(prep, "col", t, dynamic = TRUE)
			)
			W_node <- W_node_static
			if (!is.null(W_node)) {
				for (nm in intersect(colnames(W_node), cn)) {
					W_use <- if (nm %in% dyn_node_names &&
							!is.null(W_node_dyn) &&
							nm %in% colnames(W_node_dyn)) {
						W_node_dyn
					} else W_node_static
					g <- as.numeric(coef_path[nm, t] * W_use[, nm])
					xt <- xt + outer(g, g, "+")
				}
			}
		} else if (!is.null(cn)) {
			if (!is.null(prep$W_row)) {
				for (nm in intersect(colnames(prep$W_row), cn)) {
					W_row_t <- .lda_node_matrix(
						prep,
						"row",
						t,
						dynamic = nm %in% dyn_node_names
					)
					g <- as.numeric(coef_path[nm, t] * W_row_t[, nm])
					xt <- xt + matrix(g, nr, nc)
				}
			}
			if (!is.null(prep$W_col)) {
				for (nm in intersect(colnames(prep$W_col), cn)) {
					W_col_t <- .lda_node_matrix(
						prep,
						"col",
						t,
						dynamic = nm %in% dyn_node_names
					)
					g <- as.numeric(coef_path[nm, t] * W_col_t[, nm])
					xt <- xt + matrix(g, nr, nc, byrow = TRUE)
				}
			}
		}
		out[, , t] <- xt
	}
	out
}

.lda_update_dynamic_ab <- function(
	Z,
	W,
	Xb,
	Oarr,
	a_mat,
	b_mat,
	lambda_ab,
	rho_ab,
	kind = "ar1",
	transition_a = NULL,
	transition_b = NULL,
	level_prec = 0
) {
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)

	# sender paths are separable given the receiver paths, so the whole block is
	# solved at once under the per-period sum-to-zero gauge
	den_a <- matrix(0, nr, Tt)
	rhs_a <- matrix(0, nr, Tt)
	for (i in seq_len(nr)) {
		for (t in seq_len(Tt)) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			Ot <- .lda_o_slice(Oarr, t)
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			ok_idx <- which(ok)[keep]
			wt <- wt[keep]
			base <- Z[i, ok_idx, t] - Xb[i, ok_idx, t] - Ot[i, ok_idx] -
				b_mat[ok_idx, t]
			rhs_a[i, t] <- sum(wt * base)
			den_a[i, t] <- sum(wt)
		}
	}
	a_mat <- .lda_solve_paths_centred(
		den_a,
		rhs_a,
		lambda_ab,
		kind,
		rho_ab,
		transition_weights = transition_a,
		level_prec = level_prec
	)
	den_b <- matrix(0, nc, Tt)
	rhs_b <- matrix(0, nc, Tt)
	for (j in seq_len(nc)) {
		for (t in seq_len(Tt)) {
			ok <- mask[, j, t]
			if (!any(ok)) next
			Ot <- .lda_o_slice(Oarr, t)
			wt <- W[ok, j, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			ok_idx <- which(ok)[keep]
			wt <- wt[keep]
			base <- Z[ok_idx, j, t] - Xb[ok_idx, j, t] - Ot[ok_idx, j] -
				a_mat[ok_idx, t]
			rhs_b[j, t] <- sum(wt * base)
			den_b[j, t] <- sum(wt)
		}
	}
	b_mat <- .lda_solve_paths_centred(
		den_b,
		rhs_b,
		lambda_ab,
		kind,
		rho_ab,
		transition_weights = transition_b,
		level_prec = level_prec
	)
	list(a = a_mat, b = b_mat)
}

.lda_update_dynamic_ab_symmetric <- function(
	Z,
	W,
	Xb,
	Oarr,
	a_mat,
	lambda_ab,
	rho_ab,
	kind = "ar1",
	transition_a = NULL,
	level_prec = 0
) {
	n <- dim(Z)[1L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)

	# actor i enters both the row-i and column-i residuals, so the block is not
	# separable: every sub-problem is built from one snapshot of the incoming paths
	# and solved together under the sum-to-zero gauge. with M the block-diagonal
	# curvature and H the full one, 2M - H is the weighted graph Laplacian plus a
	# positive penalty diagonal, so this step always lowers the objective
	a_old <- a_mat
	den <- matrix(0, n, Tt)
	rhs <- matrix(0, n, Tt)
	for (i in seq_len(n)) {
		for (t in seq_len(Tt)) {
			Ot <- .lda_o_slice(Oarr, t)

			ok_r <- mask[i, , t]
			if (any(ok_r)) {
				wt <- W[i, ok_r, t]
				keep <- is.finite(wt) & wt > 0
				if (any(keep)) {
					idx <- which(ok_r)[keep]
					wt <- wt[keep]
					base <- Z[i, idx, t] - Xb[i, idx, t] - Ot[i, idx] -
						a_old[idx, t]
					rhs[i, t] <- rhs[i, t] + sum(wt * base)
					den[i, t] <- den[i, t] + sum(wt)
				}
			}

			ok_c <- mask[, i, t]
			if (any(ok_c)) {
				wt <- W[ok_c, i, t]
				keep <- is.finite(wt) & wt > 0
				if (any(keep)) {
					idx <- which(ok_c)[keep]
					wt <- wt[keep]
					base <- Z[idx, i, t] - Xb[idx, i, t] - Ot[idx, i] -
						a_old[idx, t]
					rhs[i, t] <- rhs[i, t] + sum(wt * base)
					den[i, t] <- den[i, t] + sum(wt)
				}
			}
		}
	}
	a_mat <- .lda_solve_paths_centred(
		den,
		rhs,
		lambda_ab,
		kind,
		rho_ab,
		transition_weights = transition_a,
		level_prec = level_prec
	)
	list(a = a_mat, b = a_mat)
}

.lda_update_static_ab <- function(Z, W, Xb, Oarr, symmetric = FALSE) {
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)
	base <- array(0, dim = dim(Z))
	for (t in seq_len(Tt)) {
		Ot <- .lda_o_slice(Oarr, t)
		base[, , t] <- ifelse(
			mask[, , t],
			W[, , t] * (Z[, , t] - Xb[, , t] - Ot),
			0
		)
	}
	wsum <- apply(ifelse(mask, W, 0), c(1, 2), sum)
	base_rs <- apply(base, 1, sum)
	base_cs <- apply(base, 2, sum)
	cnt_row <- rowSums(wsum)
	cnt_col <- colSums(wsum)
	ab <- .ae_additive_solve(
		base_rs,
		base_cs,
		wsum,
		cnt_row,
		cnt_col,
		symmetric = symmetric
	)
	# sum-to-zero identification; the removed grand level is returned so the
	# caller can fold it into the intercept, which leaves the fit and objective
	# unchanged
	level <- mean(ab$a) + mean(ab$b)
	ab$a <- ab$a - mean(ab$a)
	ab$b <- ab$b - mean(ab$b)
	list(a = matrix(ab$a, nr, Tt), b = matrix(ab$b, nc, Tt), level = level)
}

.lda_update_static_lowrank <- function(
	Z,
	W,
	Xb,
	a_mat,
	b_mat,
	Omat,
	R,
	symmetric,
	lowrank_method,
	ridge = 0,
	rel = FALSE
) {
	if (R <= 0L) {
		return(list(
			Omat = matrix(0, dim(Z)[1L], dim(Z)[2L]),
			U = NULL,
			V = NULL,
			L = NULL
		))
	}
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)
	wsum <- apply(ifelse(mask, W, 0), c(1, 2), sum)
	res_sum <- matrix(0, nr, nc)
	for (t in seq_len(Tt)) {
		Wt <- W[, , t]
		Wt[!mask[, , t] | !is.finite(Wt)] <- 0
		res_t <- Z[, , t] - Xb[, , t] - outer(a_mat[, t], b_mat[, t], "+")
		res_t[!mask[, , t]] <- 0
		res_sum <- res_sum + Wt * res_t
	}
	obs_dyad <- wsum > 0
	Rbar <- matrix(0, nr, nc)
	Rbar[obs_dyad] <- res_sum[obs_dyad] / wsum[obs_dyad]
	# count and binary families carry the ridge of the level prior on the factors (see dynamic_als_lowrank_ridge.R)
	if (ridge > 0) {
		return(.lda_lowrank_update_ridge(
			Omat,
			Rbar,
			wsum,
			obs_dyad,
			R,
			symmetric,
			lowrank_method,
			ridge
		))
	}
	# the unweighted least-squares fit of a normal response measures its tolerances in the response unit
	.ae_lowrank_update(
		Omat,
		Rbar,
		wsum,
		obs_dyad,
		R,
		symmetric,
		lowrank_method,
		rel = rel
	)
}

# marginal level prior precision for the multiplicative terms. the difference
# penalty only constrains u_t - rho u_{t-1}, so without a level prior a binary
# or Poisson IRLS fit can push the factor paths to arbitrarily large levels
# (quasi-separation). this is the ALS counterpart of the sampler's N(0, s2 I)
# factor prior on the working scale, held fixed so block minimization stays a
# descent method
.lda_level_prec <- 1

.lda_update_uv_path_block <- function(
	target,
	mask,
	W,
	other_cube,
	margin,
	lambda_uv,
	rho_uv,
	transition_weights = NULL,
	ridge = 1e-6,
	level_prec = 0
) {
	d <- dim(other_cube)
	R <- d[2L]
	Tt <- d[3L]
	n_path <- if (identical(margin, "row")) dim(target)[1L] else dim(target)[2L]
	out <- array(0, dim = c(n_path, R, Tt))

	for (idx_path in seq_len(n_path)) {
		tw <- if (is.null(transition_weights)) rep(1, Tt) else
			transition_weights[idx_path, ]
		P <- kronecker(
			.lda_path_penalty_weighted(Tt, lambda_uv, tw, rho_uv),
			diag(R)
		)
		H <- P
		g <- numeric(R * Tt)
		for (t in seq_len(Tt)) {
			block <- ((t - 1L) * R + 1L):(t * R)
			if (identical(margin, "row")) {
				ok <- mask[idx_path, , t]
				if (!any(ok)) next
				wt <- W[idx_path, ok, t]
				keep <- is.finite(wt) & wt > 0
				if (!any(keep)) next
				Xt <- matrix(other_cube[ok, , t], ncol = R)
				yt <- target[idx_path, ok, t]
			} else {
				ok <- mask[, idx_path, t]
				if (!any(ok)) next
				wt <- W[ok, idx_path, t]
				keep <- is.finite(wt) & wt > 0
				if (!any(keep)) next
				Xt <- matrix(other_cube[ok, , t], ncol = R)
				yt <- target[ok, idx_path, t]
			}
			Xt <- Xt[keep, , drop = FALSE]
			yt <- yt[keep]
			wt <- wt[keep]
			H[block, block] <- H[block, block] + crossprod(Xt * sqrt(wt))
			g[block] <- g[block] + as.numeric(crossprod(Xt, wt * yt))
		}
		# level prior on every time block plus the numerical floor
		H <- H + diag(level_prec + ridge, R * Tt)
		sol <- .ae_safe_solve(H, g)
		sol[!is.finite(sol)] <- 0
		for (t in seq_len(Tt)) {
			block <- ((t - 1L) * R + 1L):(t * R)
			out[idx_path, , t] <- sol[block]
		}
	}
	out
}

# sign pattern s of a symmetric factor pair, read off the incoming cubes:
# v = u s, so sum_it u_itr v_itr carries the sign of s_r
.lda_symmetric_signs <- function(U_cube, V_cube) {
	R <- dim(U_cube)[2L]
	if (R == 0L) return(numeric(0))
	vapply(seq_len(R), function(r) {
		cs <- sum(U_cube[, r, ] * V_cube[, r, ], na.rm = TRUE)
		if (!is.finite(cs) || cs >= 0) 1 else -1
	}, numeric(1))
}

# symmetric dynamic multiplicative block: o_t = u_t s u_t' with s = diag(+/-1),
# so actor i enters both its row and column residuals. actors are swept in
# sequence, each solve the exact minimiser in that actor's whole path given the
# others, which keeps the sweep monotone
.lda_update_dynamic_uv_symmetric <- function(
	Z,
	W,
	Xb,
	a_mat,
	b_mat,
	U_cube,
	V_cube,
	lambda_uv,
	rho_uv,
	transition_u = NULL,
	transition_v = NULL,
	ridge = 1e-6,
	level_prec = 0
) {
	udn <- dimnames(U_cube)
	vdn <- dimnames(V_cube)
	d <- dim(U_cube)
	n <- d[1L]
	R <- d[2L]
	Tt <- d[3L]
	sgn <- .lda_symmetric_signs(U_cube, V_cube)
	U <- U_cube
	if (R > 0L && n > 1L) {
		target <- array(0, dim = dim(Z))
		for (t in seq_len(Tt)) {
			target[, , t] <- Z[, , t] - Xb[, , t] -
				outer(a_mat[, t], b_mat[, t], "+")
		}
		obs <- is.finite(Z) & is.finite(W) & W > 0
		Wz <- W
		Wz[!obs] <- 0
		Tz <- target
		Tz[!obs] <- 0
		for (i in seq_len(n)) {
			others <- seq_len(n)[-i]
			tw_u <- if (is.null(transition_u)) rep(1, Tt) else transition_u[i, ]
			tw_v <- if (is.null(transition_v)) rep(1, Tt) else transition_v[i, ]
			# both factors carry the same path, so both smoothness penalties and both
			# level priors act on it
			P <- kronecker(
				.lda_path_penalty_weighted(Tt, lambda_uv, tw_u, rho_uv) +
					.lda_path_penalty_weighted(Tt, lambda_uv, tw_v, rho_uv),
				diag(R)
			)
			H <- P + diag(2 * level_prec + ridge, R * Tt)
			g <- numeric(R * Tt)
			for (t in seq_len(Tt)) {
				block <- ((t - 1L) * R + 1L):(t * R)
				Xt <- matrix(U[others, , t], ncol = R)
				Xt <- sweep(Xt, 2, sgn, "*")
				w_row <- Wz[i, others, t]
				w_col <- Wz[others, i, t]
				y_row <- Tz[i, others, t]
				y_col <- Tz[others, i, t]
				w_tot <- w_row + w_col
				if (!any(w_tot > 0)) next
				H[block, block] <- H[block, block] +
					crossprod(Xt * sqrt(w_tot))
				g[block] <- g[block] +
					as.numeric(crossprod(Xt, w_row * y_row + w_col * y_col))
			}
			sol <- .ae_safe_solve(H, g)
			sol[!is.finite(sol)] <- 0
			for (t in seq_len(Tt)) {
				U[i, , t] <- sol[((t - 1L) * R + 1L):(t * R)]
			}
		}
	}
	V <- U
	if (R > 0L) for (r in seq_len(R)) V[, r, ] <- U[, r, ] * sgn[r]
	dimnames(U) <- udn
	dimnames(V) <- vdn
	list(U = U, V = V, Oarr = .lda_uv_to_o_array(U, V))
}

.lda_update_dynamic_uv_directed <- function(
	Z,
	W,
	Xb,
	a_mat,
	b_mat,
	U_cube,
	V_cube,
	lambda_uv,
	rho_uv,
	transition_u = NULL,
	transition_v = NULL,
	align = TRUE,
	level_prec = 0,
	ridge = 1e-6
) {
	Tt <- dim(Z)[3L]
	udn <- dimnames(U_cube)
	vdn <- dimnames(V_cube)
	target <- array(0, dim = dim(Z))
	for (t in seq_len(Tt)) {
		target[, , t] <- Z[, , t] - Xb[, , t] -
			outer(a_mat[, t], b_mat[, t], "+")
	}
	mask <- is.finite(Z)
	U_cube <- .lda_update_uv_path_block(
		target,
		mask,
		W,
		V_cube,
		margin = "row",
		lambda_uv = lambda_uv,
		rho_uv = rho_uv,
		transition_weights = transition_u,
		level_prec = level_prec,
		ridge = ridge
	)
	dimnames(U_cube) <- udn
	V_cube <- .lda_update_uv_path_block(
		target,
		mask,
		W,
		U_cube,
		margin = "col",
		lambda_uv = lambda_uv,
		rho_uv = rho_uv,
		transition_weights = transition_v,
		level_prec = level_prec,
		ridge = ridge
	)
	dimnames(V_cube) <- vdn
	if (isTRUE(align)) {
		uv <- .lda_align_uv_cube(
			U_cube,
			V_cube,
			weights_u = transition_u,
			weights_v = transition_v
		)
		U_cube <- uv$U
		V_cube <- uv[["V"]]
	}
	Oarr <- .lda_uv_to_o_array(U_cube, V_cube)
	list(U = U_cube, V = V_cube, Oarr = Oarr)
}

.lda_update_dynamic_uv_bip_g <- function(
	Z,
	W,
	Xb,
	a_mat,
	b_mat,
	U_cube,
	V_cube,
	G,
	lambda_uv,
	rho_uv,
	transition_u = NULL,
	transition_v = NULL,
	level_prec = 0,
	ridge = 1e-6
) {
	Tt <- dim(Z)[3L]
	RA <- dim(U_cube)[2L]
	RB <- dim(V_cube)[2L]
	target <- array(0, dim = dim(Z))
	for (t in seq_len(Tt)) {
		target[, , t] <- Z[, , t] - Xb[, , t] -
			outer(a_mat[, t], b_mat[, t], "+")
	}
	mask <- is.finite(Z)
	V_eff <- array(
		0,
		dim = c(dim(V_cube)[1L], RA, Tt),
		dimnames = list(
			dimnames(V_cube)[[1L]],
			dimnames(U_cube)[[2L]],
			dimnames(V_cube)[[3L]]
		)
	)
	for (t in seq_len(Tt)) {
		Vt <- matrix(V_cube[, , t], nrow = dim(V_cube)[1L], ncol = RB)
		V_eff[, , t] <- Vt %*% t(G)
	}
	udn <- dimnames(U_cube)
	vdn <- dimnames(V_cube)
	gdn <- dimnames(G)
	U_cube <- .lda_update_uv_path_block(
		target,
		mask,
		W,
		V_eff,
		margin = "row",
		lambda_uv = lambda_uv,
		rho_uv = rho_uv,
		transition_weights = transition_u,
		level_prec = level_prec,
		ridge = ridge
	)
	dimnames(U_cube) <- udn
	U_eff <- array(
		0,
		dim = c(dim(U_cube)[1L], RB, Tt),
		dimnames = list(
			dimnames(U_cube)[[1L]],
			dimnames(V_cube)[[2L]],
			dimnames(U_cube)[[3L]]
		)
	)
	for (t in seq_len(Tt)) {
		Ut <- matrix(U_cube[, , t], nrow = dim(U_cube)[1L], ncol = RA)
		U_eff[, , t] <- Ut %*% G
	}
	V_cube <- .lda_update_uv_path_block(
		target,
		mask,
		W,
		U_eff,
		margin = "col",
		lambda_uv = lambda_uv,
		rho_uv = rho_uv,
		transition_weights = transition_v,
		level_prec = level_prec,
		ridge = ridge
	)
	dimnames(V_cube) <- vdn
	dimnames(G) <- gdn
	list(
		U = U_cube,
		V = V_cube,
		G = G,
		Oarr = .lda_uvg_to_o_array(U_cube, V_cube, G)
	)
}

.lda_update_static_uv_bip_g <- function(
	target,
	mask,
	W,
	U,
	V,
	G,
	ridge = 1e-6,
	level_prec = 0
) {
	dg <- dim(G)
	RA <- ncol(U)
	RB <- ncol(V)
	Tt <- dim(target)[3L]
	udn <- dimnames(U)
	vdn <- dimnames(V)

	for (i in seq_len(nrow(U))) {
		H <- diag(level_prec + ridge, RA)
		g <- numeric(RA)
		for (t in seq_len(Tt)) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			jj <- which(ok)[keep]
			wt <- wt[keep]
			Gt <- if (length(dg) == 3L) matrix(G[, , t], RA, RB) else G
			Xt <- matrix(0, length(jj), RA)
			for (m in seq_along(jj)) Xt[m, ] <- as.numeric(V[jj[m], ] %*% t(Gt))
			Xtw <- Xt * sqrt(wt)
			H <- H + crossprod(Xtw)
			g <- g + as.numeric(crossprod(Xt, wt * target[i, jj, t]))
		}
		U[i, ] <- .ae_safe_solve(H, g)
	}
	dimnames(U) <- udn

	for (j in seq_len(nrow(V))) {
		H <- diag(level_prec + ridge, RB)
		g <- numeric(RB)
		for (t in seq_len(Tt)) {
			ok <- mask[, j, t]
			if (!any(ok)) next
			wt <- W[ok, j, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			ii <- which(ok)[keep]
			wt <- wt[keep]
			Gt <- if (length(dg) == 3L) matrix(G[, , t], RA, RB) else G
			Xt <- matrix(0, length(ii), RB)
			for (m in seq_along(ii)) Xt[m, ] <- as.numeric(U[ii[m], ] %*% Gt)
			Xtw <- Xt * sqrt(wt)
			H <- H + crossprod(Xtw)
			g <- g + as.numeric(crossprod(Xt, wt * target[ii, j, t]))
		}
		V[j, ] <- .ae_safe_solve(H, g)
	}
	dimnames(V) <- vdn
	list(
		U = U,
		V = V,
		Oarr = .lda_uvg_to_o_array(U, V, G, Tt = dim(target)[3L])
	)
}

.lda_update_static_g_bip <- function(
	target,
	mask,
	W,
	U,
	V,
	G,
	ridge = 1e-6,
	level_prec = 0
) {
	du <- dim(U)
	dv <- dim(V)
	RA <- if (length(du) == 3L) du[2L] else ncol(U)
	RB <- if (length(dv) == 3L) dv[2L] else ncol(V)
	Tt <- dim(target)[3L]
	p <- RA * RB
	H <- diag(level_prec + ridge, p)
	g <- numeric(p)
	for (t in seq_len(Tt)) {
		Ut <- if (length(du) == 3L) matrix(U[, , t], nrow = du[1L], ncol = RA) else U
		Vt <- if (length(dv) == 3L) matrix(V[, , t], nrow = dv[1L], ncol = RB) else V
		for (i in seq_len(nrow(target))) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			jj <- which(ok)[keep]
			wt <- wt[keep]
			yt <- target[i, jj, t]
			for (m in seq_along(jj)) {
				x <- as.vector(outer(Ut[i, ], Vt[jj[m], ]))
				H <- H + wt[m] * tcrossprod(x)
				g <- g + wt[m] * x * yt[m]
			}
		}
	}
	out <- matrix(.ae_safe_solve(H, g), RA, RB)
	out[!is.finite(out)] <- 0
	dimnames(out) <- dimnames(G)
	out
}

# degrees of freedom of the heavy-tailed transition kind
.lda_t_nu <- 4

# the t kind is not fitted on bipartite panels: the interaction matrix carries
# the scale of the positions, so the innovations the weights read stay small and
# the fit is in effect an ar1 fit with a stiffer penalty. lame() and
# lame_dynamic_als() refuse it with this message
.lda_refuse_bipartite_t <- function(caller = c("lame", "lame_dynamic_als")) {
	caller <- match.arg(caller)
	head <- if (identical(caller, "lame")) {
		"{.code dynamic_uv_kind = \"t\"} with {.code method = \"als\"} is not available for bipartite networks."
	} else {
		"{.fn lame_dynamic_als} does not fit {.code dynamic_uv_kind = \"t\"} on bipartite networks."
	}
	cli::cli_abort(
		c(
			head,
			"i" = "On simulated bipartite t panels its transition weights never marked a transition as heavy (all stayed between 1.2 and 1.5), so it fitted the latent product no differently from the {.val ar1} kind (correlation with the generating product 0.68 against 0.68).",
			"i" = "Use {.code lame(..., dynamic_uv = TRUE, dynamic_uv_kind = \"ar1\")} with {.code method = \"als\"} or {.code method = \"mcmc\"}."
		),
		call = parent.frame()
	)
}

# squared innovation length per actor and period of a factor path cube
.lda_path_innov_ss <- function(cube, rho) {
	d <- dim(cube)
	if (is.null(d) || length(d) != 3L || d[3L] < 2L) return(NULL)
	innov <- cube[, , -1L, drop = FALSE] -
		rho * cube[, , -d[3L], drop = FALSE]
	apply(innov^2, c(1, 3), sum)
}

# the "t" transition kind replaces the Gaussian AR(1) penalty lambda * ss by the
# multivariate-t(nu) penalty
#   ((nu + R) / 2) * log(1 + 2 * lambda * ss / nu),
# which is concave in ss, so its tangent at the current iterate, the weighted
# quadratic
#   lambda * w * ss  with  w = (nu + R) / (nu + 2 * lambda * ss),
# majorises it. holding the weights at the incoming iterate makes the sweep an
# MM step on one fixed objective; the scale is a fixed hyperparameter
.lda_t_transition_weights <- function(cube, rho, lambda, nu = .lda_t_nu) {
	d <- dim(cube)
	dn <- dimnames(cube)
	out <- matrix(
		1,
		d[1L],
		d[3L],
		dimnames = if (!is.null(dn)) list(dn[[1L]], dn[[3L]]) else NULL
	)
	ss <- .lda_path_innov_ss(cube, rho)
	if (is.null(ss)) return(out)
	if (!is.finite(lambda) || lambda <= 0) return(out)
	out[, -1L] <- (nu + d[2L]) / (nu + 2 * lambda * ss)
	out[!is.finite(out) | out <= 0] <- 1
	out
}

# marginal t(nu) transition penalty of one factor path cube, the quantity the
# weighted quadratic majorises
.lda_t_path_penalty <- function(
	cube,
	rho,
	lambda,
	weights = NULL,
	nu = .lda_t_nu
) {
	ss <- .lda_path_innov_ss(cube, rho)
	if (is.null(ss) || !is.finite(lambda) || lambda <= 0) return(0)
	R <- dim(cube)[2L]
	term <- ((nu + R) / 2) * log1p(2 * lambda * ss / nu)
	if (!is.null(weights) && all(dim(weights) == c(nrow(ss), ncol(ss) + 1L))) {
		term <- term * pmax(0, weights[, -1L, drop = FALSE])
	}
	sum(term)
}

# a warm start replaces the internal static seed pass; it must carry the
# static-fit pieces at this problem's dimensions, otherwise the cold start is used
.lda_usable_start <- function(
	start,
	prep,
	R_eff,
	symmetric,
	R_row = R_eff,
	R_col = R_eff
) {
	if (is.null(start)) return(NULL)
	if (!is.list(start)) return(NULL)
	ok_vec <- function(x, n) is.numeric(x) && length(x) == n && all(is.finite(x))
	if (!ok_vec(start$a, prep$nr)) return(NULL)
	if (!ok_vec(start$b, prep$nc)) return(NULL)
	if (!is.numeric(start$coefficients) || is.null(names(start$coefficients))) {
		return(NULL)
	}
	if (R_eff > 0L) {
		ok_mat <- function(x, n, rank) is.matrix(x) && nrow(x) == n &&
			ncol(x) == rank && all(is.finite(x))
		if (!ok_mat(start$U, prep$nr, R_row)) return(NULL)
		if (!(symmetric && is.null(start[["V"]])) &&
				!ok_mat(start[["V"]], prep$nc, R_col)) {
			return(NULL)
		}
	}
	start
}

.lda_fit_surface_delta <- function(base, other) {
	if (is.null(base$EZ) || is.null(other$EZ)) return(NA_real_)
	max(vapply(seq_len(base$n_time), function(t) {
		max(abs(base$EZ[[t]] - other$EZ[[t]]), na.rm = TRUE)
	}, numeric(1)), na.rm = TRUE)
}

.lda_dynamic_stability <- function(
	fit,
	starts = c("quick", "validation"),
	seed = 6886,
	call_args = list()
) {
	starts <- match.arg(starts)
	n_start <- if (identical(starts, "quick")) 2L else 4L
	refits <- vector("list", n_start)
	for (s in seq_len(n_start)) {
		args <- call_args
		args$seed <- seed + 1000L + s
		args$verbose <- FALSE
		args$bootstrap <- 0L
		args$stability <- "none"
		args$start_jitter <- 0.05
		args$call <- NULL
		refits[[s]] <- tryCatch(
			do.call(lame_dynamic_als, args),
			error = function(e) e
		)
	}
	ok <- vapply(refits, function(x) inherits(x, "lame_dynamic_als"), logical(1))
	valid <- refits[ok]
	surface_delta <- if (length(valid) > 0L) {
		vapply(valid, function(x) .lda_fit_surface_delta(fit, x), numeric(1))
	} else numeric(0)
	coef_delta <- if (length(valid) > 0L) {
		vapply(
			valid,
			function(x) max(abs(x$coef_path - fit$coef_path), na.rm = TRUE),
			numeric(1)
		)
	} else numeric(0)
	list(
		preset = starts,
		n_start = n_start,
		n_success = sum(ok),
		converged = vapply(valid, function(x) isTRUE(x$converged), logical(1)),
		max_surface_delta = if (length(surface_delta)) max(surface_delta, na.rm = TRUE) else NA_real_,
		max_coef_delta = if (length(coef_delta)) max(coef_delta, na.rm = TRUE) else NA_real_,
		refits = valid,
		failures = refits[!ok]
	)
}

.lda_dynamic_bootstrap <- function(
	fit,
	R_boot,
	type,
	block_length,
	seed,
	call_args = list(),
	verbose = FALSE
) {
	if (R_boot <= 0L) return(NULL)
	.ae_check_bootstrap_n(R_boot)
	if (!isTRUE(fit$converged)) {
		cli::cli_abort(
			c(
				"Bootstrap intervals require a converged original fit.",
				"i" = "Refit with a larger iteration budget and inspect the convergence diagnostics."
			),
			class = "lame_als_bootstrap_convergence"
		)
	}
	if (!identical(type, "parametric")) {
		cli::cli_warn(c(
			"{.arg bootstrap_type = \"block\"} is not yet available for dynamic ALS.",
			"i" = "Using parametric dynamic refits instead."
		))
		type <- "parametric"
	}
	sims <- .lda_boot_sims(fit, R_boot, seed)
	obs_mask <- is.finite(fit$Y)
	mask_sim <- function(Ysim) {
		if (is.list(Ysim) && !is.data.frame(Ysim)) {
			out <- Ysim
			for (t in seq_along(out)) {
				out[[t]][!obs_mask[, , t]] <- NA_real_
				dimnames(out[[t]]) <- dimnames(fit$Y)[1:2]
				if (isTRUE(fit$changing_composition)) {
					row_keep <- (fit$row_presence %||% fit$actor_presence)[, t]
					col_keep <- (fit$col_presence %||% fit$actor_presence)[, t]
					out[[t]] <- out[[t]][row_keep, col_keep, drop = FALSE]
				}
			}
			names(out) <- dimnames(fit$Y)[[3L]]
			return(out)
		}
		if (is.array(Ysim) && length(dim(Ysim)) == 3L) {
			Ysim[!obs_mask] <- NA_real_
			dimnames(Ysim) <- dimnames(fit$Y)
			if (isTRUE(fit$changing_composition)) {
				out <- lapply(seq_len(dim(Ysim)[3L]), function(t) {
					row_keep <- (fit$row_presence %||% fit$actor_presence)[, t]
					col_keep <- (fit$col_presence %||% fit$actor_presence)[, t]
					Ysim[row_keep, col_keep, t, drop = FALSE][, , 1L]
				})
				names(out) <- dimnames(fit$Y)[[3L]]
				return(out)
			}
			return(Ysim)
		}
		Ysim
	}
	# the replicate loop and its accounting are in ame_als_bootstrap_dynamic.R
	.lda_boot_run(fit, sims, mask_sim, call_args, seed, R_boot, verbose)
}

.lda_eta <- function(
	prep,
	coef_path,
	a_mat,
	b_mat,
	Oarr,
	symmetric = FALSE,
	beta_dyn = NULL
) {
	Xb <- .lda_xbeta_array(
		prep,
		coef_path,
		symmetric = symmetric,
		beta_dyn = beta_dyn
	)
	E <- vector("list", prep$Tt)
	for (t in seq_len(prep$Tt)) {
		E[[t]] <- Xb[, , t] + outer(a_mat[, t], b_mat[, t], "+") +
			.lda_o_slice(Oarr, t)
	}
	E
}

.lda_objective <- function(
	Z,
	W,
	prep,
	coef_path,
	a_mat,
	b_mat,
	Oarr,
	beta_dyn,
	lambda_beta,
	beta_kind,
	rho_beta,
	dynamic_ab,
	lambda_ab,
	rho_ab,
	transition_a = NULL,
	transition_b = NULL,
	symmetric = FALSE,
	dynamic_uv = FALSE,
	U = NULL,
	V = NULL,
	lambda_uv = 0,
	rho_uv = 1,
	transition_u = NULL,
	transition_v = NULL,
	dynamic_uv_kind = "ar1",
	G_static = NULL,
	level_prec_uv = 0,
	level_prec_g = 0,
	level_prec_ab = 0,
	L = NULL
) {
	E <- .lda_eta(
		prep,
		coef_path,
		a_mat,
		b_mat,
		Oarr,
		symmetric = symmetric,
		beta_dyn = beta_dyn
	)
	sse <- 0
	for (t in seq_len(prep$Tt)) {
		ok <- is.finite(Z[, , t])
		r <- Z[, , t] - E[[t]]
		sse <- sse + sum(W[, , t][ok] * r[ok]^2)
	}
	pen_ab <- 0
	if (isTRUE(dynamic_ab) && lambda_ab > 0 && prep$Tt > 1L) {
		for (i in seq_len(nrow(a_mat))) {
			d <- a_mat[i, -1L] - rho_ab * a_mat[i, -prep$Tt]
			wa <- if (is.null(transition_a)) rep(1, prep$Tt - 1L) else
				pmax(0, transition_a[i, -1L])
			pen_ab <- pen_ab + lambda_ab * sum(wa * d^2)
		}
		if (!isTRUE(symmetric)) {
			for (j in seq_len(nrow(b_mat))) {
				d <- b_mat[j, -1L] - rho_ab * b_mat[j, -prep$Tt]
				wb <- if (is.null(transition_b)) rep(1, prep$Tt - 1L) else
					pmax(0, transition_b[j, -1L])
				pen_ab <- pen_ab + lambda_ab * sum(wb * d^2)
			}
		}
	}
	pen_beta <- 0
	if (isTRUE(beta_dyn$any) && lambda_beta > 0) {
		for (k in which(beta_dyn$mask)) {
			P <- .lda_path_penalty(prep$Tt, lambda_beta * .lda_beta_ms(beta_dyn, k), beta_kind, rho_beta)
			pen_beta <- pen_beta + as.numeric(t(coef_path[k, ]) %*% P %*% coef_path[k, ])
		}
	}
	pen_uv <- 0
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") &&
			lambda_uv > 0 && prep$Tt > 1L &&
			!is.null(U) && !is.null(V) && length(dim(U)) == 3L) {
		# the heavy-tailed kind carries its own marginal penalty; of the majoriser's
		# weights in transition_u / transition_v only the presence part belongs here
		pres_u <- .lda_transition_from_presence(prep$row_presence)
		pres_v <- .lda_transition_from_presence(prep$col_presence)
		pen_uv <- .lda_t_path_penalty(U, rho_uv, lambda_uv, pres_u) +
			.lda_t_path_penalty(V, rho_uv, lambda_uv, pres_v)
	} else if (isTRUE(dynamic_uv) && lambda_uv > 0 && prep$Tt > 1L &&
			!is.null(U) && !is.null(V) && length(dim(U)) == 3L) {
		for (i in seq_len(dim(U)[1L])) {
			d <- U[i, , -1L, drop = FALSE] -
				rho_uv * U[i, , -prep$Tt, drop = FALSE]
			wu <- if (is.null(transition_u)) rep(1, prep$Tt - 1L) else
				pmax(0, transition_u[i, -1L])
			pen_uv <- pen_uv + lambda_uv * sum(wu * apply(d^2, 3, sum))
		}
		for (j in seq_len(dim(V)[1L])) {
			d <- V[j, , -1L, drop = FALSE] -
				rho_uv * V[j, , -prep$Tt, drop = FALSE]
			wv <- if (is.null(transition_v)) rep(1, prep$Tt - 1L) else
				pmax(0, transition_v[j, -1L])
			pen_uv <- pen_uv + lambda_uv * sum(wv * apply(d^2, 3, sum))
		}
	}
	# fixed-precision marginal level prior on the multiplicative terms, matching
	# the diag(level_prec) of the factor block solves
	pen_level <- 0
	if (level_prec_ab > 0) {
		pen_level <- pen_level + level_prec_ab * (sum(a_mat^2, na.rm = TRUE) +
				if (is.null(b_mat) || isTRUE(symmetric)) 0 else sum(b_mat^2, na.rm = TRUE))
	}
	if (level_prec_uv > 0 && !is.null(U) && !is.null(V)) {
		if (isTRUE(symmetric) && !is.null(L) && length(dim(U)) == 2L) {
			# static eigenmodel U diag(L) U' has the factor U |L|^(1/2) twice (U orthonormal): 2 sum |L|
			pen_level <- pen_level + level_prec_uv * 2 * sum(abs(L))
		} else {
			pen_level <- pen_level + level_prec_uv * (sum(U^2) + sum(V^2))
		}
	}
	if (level_prec_g > 0 && !is.null(G_static)) {
		pen_level <- pen_level + level_prec_g * sum(G_static^2)
	}
	sse + pen_ab + pen_beta + pen_uv + pen_level
}

# deviance of a non-normal fit on the data scale plus the penalty terms the
# block updates minimize; the outer IRLS gate accepts a reweighted iterate only
# when this falls, which keeps sparse panels from running away under
# near-separation
.lda_penalized_deviance <- function(prep, family, link, eta_list, pen) {
	eta <- .lda_eta_array(eta_list)
	y <- prep$Yarr
	dev <- .ae_irls_deviance(y, eta, family, link)
	dev + pen
}

.lda_default_lambdas <- function(prep, lambda_ab, lambda_beta, lambda_uv) {
	mask <- is.finite(prep$Yarr)
	row_deg <- numeric(0)
	for (t in seq_len(prep$Tt)) {
		row_deg <- c(row_deg, rowSums(mask[, , t]), colSums(mask[, , t]))
	}
	deg_pos <- row_deg[row_deg > 0]
	deg_med <- if (length(deg_pos) > 0L) stats::median(deg_pos, na.rm = TRUE) else 1
	if (is.null(lambda_ab)) {
		lambda_ab <- max(1, deg_med / 4)
	}
	if (is.null(lambda_beta)) {
		nobs_t <- vapply(seq_len(prep$Tt), function(t) sum(mask[, , t]), numeric(1))
		nobs_pos <- nobs_t[nobs_t > 0]
		nobs_med <- if (length(nobs_pos) > 0L) stats::median(nobs_pos, na.rm = TRUE) else 1
		lambda_beta <- max(1, sqrt(nobs_med))
	}
	if (is.null(lambda_uv)) {
		lambda_uv <- max(1, deg_med / 4)
	}
	list(ab = lambda_ab, beta = lambda_beta, uv = lambda_uv)
}

.lda_decompose_margin_node_paths <- function(
	coef_path,
	e_mat,
	W,
	active,
	time_names,
	intercept_multiplier = 1,
	W_arr = NULL,
	presence = NULL
) {
	if (is.null(W) || ncol(W) == 0L) {
		return(list(
			coef_path = coef_path,
			e_mat = e_mat,
			beta_names = character(0)
		))
	}
	Tt <- ncol(e_mat)
	n_actor <- nrow(e_mat)
	node_names <- colnames(W)
	dynamic_names <- intersect(node_names, rownames(coef_path))
	static_names <- setdiff(node_names, dynamic_names)
	beta_names <- character(0)

	if (length(dynamic_names) > 0L) {
		for (t in seq_len(Tt)) {
			W_t <- if (!is.null(W_arr)) {
				tmp <- W_arr[, , t, drop = FALSE]
				tmp <- matrix(
					tmp,
					nrow = dim(W_arr)[1L],
					ncol = dim(W_arr)[2L]
				)
				dimnames(tmp) <- list(NULL, dimnames(W_arr)[[2L]])
				tmp
			} else W
			W_dyn <- W_t[, dynamic_names, drop = FALSE]
			active_t <- if (is.matrix(active) && ncol(active) >= t) {
				as.numeric(active[, t])
			} else active
			dec <- .ae_decompose(e_mat[, t], W_dyn, active_t)
			e_mat[, t] <- dec$resid
			coef_path[1L, t] <- coef_path[1L, t] +
				intercept_multiplier * dec$intercept
			coef_path[names(dec$beta), t] <-
				coef_path[names(dec$beta), t] + dec$beta
		}
		beta_names <- c(beta_names, dynamic_names)
	}

	if (length(static_names) > 0L) {
		W_static <- W[, static_names, drop = FALSE]
		active_static <- if (is.matrix(active)) {
			as.numeric(rowSums(active) > 0)
		} else active
		dec <- .ae_decompose(
			.lda_presence_row_means(e_mat, presence),
			W_static,
			active_static
		)
		if (length(dec$beta) > 0L) {
			node_level <- as.numeric(dec$intercept + W_static %*% dec$beta)
			e_mat <- e_mat - matrix(node_level, n_actor, Tt)
			coef_path[1L, ] <- coef_path[1L, ] +
				intercept_multiplier * dec$intercept
			node_path <- matrix(
				rep(dec$beta, Tt),
				nrow = length(dec$beta),
				dimnames = list(names(dec$beta), time_names)
			)
			coef_path <- rbind(coef_path, node_path)
			beta_names <- c(beta_names, names(dec$beta))
		}
	}

	list(
		coef_path = coef_path,
		e_mat = e_mat,
		beta_names = unique(beta_names)
	)
}

.lda_decompose_node_paths <- function(prep, symmetric, coef_path, a_mat, b_mat) {
	has_node <- !is.null(prep$W_row) || !is.null(prep$W_col)
	if (!isTRUE(has_node)) {
		return(list(
			coef_path = coef_path,
			a_mat = a_mat,
			b_mat = b_mat,
			beta_names = character(0)
		))
	}

	Tt <- prep$Tt
	obs_count <- apply(!is.na(prep$Yarr), c(1, 2), sum)
	obs_row_time <- matrix(FALSE, prep$nr, Tt)
	obs_col_time <- matrix(FALSE, prep$nc, Tt)
	for (t in seq_len(Tt)) {
		yt <- prep$Yarr[, , t]
		obs_row_time[, t] <- rowSums(!is.na(yt)) > 0
		obs_col_time[, t] <- colSums(!is.na(yt)) > 0
	}
	dimnames(obs_row_time) <- list(prep$row_names, prep$time_names)
	dimnames(obs_col_time) <- list(prep$col_names, prep$time_names)
	if (isTRUE(symmetric)) obs_count <- obs_count + t(obs_count)
	active_row_time <- obs_row_time
	active_col_time <- obs_col_time

	if (isTRUE(symmetric)) {
		W_node <- prep$W_row
		if (!is.null(prep$W_col)) {
			W_node <- if (is.null(W_node)) prep$W_col else cbind(W_node, prep$W_col)
		}
		W_node_arr <- .lda_bind_node_arrays(prep$W_row_arr, prep$W_col_arr)
		active_sym_time <- obs_row_time | obs_col_time
		dec <- .lda_decompose_margin_node_paths(
			coef_path,
			a_mat,
			W_node,
			active_sym_time,
			prep$time_names,
			intercept_multiplier = 2,
			W_arr = W_node_arr,
			presence = prep$row_presence
		)
		coef_path <- dec$coef_path
		a_mat <- dec$e_mat
		b_mat <- a_mat
		return(list(
			coef_path = coef_path,
			a_mat = a_mat,
			b_mat = b_mat,
			beta_names = dec$beta_names
		))
	}

	dec_r <- .lda_decompose_margin_node_paths(
		coef_path,
		a_mat,
		prep$W_row,
		active_row_time,
		prep$time_names,
		W_arr = prep$W_row_arr,
		presence = prep$row_presence
	)
	coef_path <- dec_r$coef_path
	a_mat <- dec_r$e_mat
	dec_c <- .lda_decompose_margin_node_paths(
		coef_path,
		b_mat,
		prep$W_col,
		active_col_time,
		prep$time_names,
		W_arr = prep$W_col_arr,
		presence = prep$col_presence
	)
	coef_path <- dec_c$coef_path
	b_mat <- dec_c$e_mat
	list(
		coef_path = coef_path,
		a_mat = a_mat,
		b_mat = b_mat,
		beta_names = unique(c(dec_r$beta_names, dec_c$beta_names))
	)
}

# innovation scale implied by fitted AR(1) paths: the sd of x_t - rho x_{t-1}
# pooled over every series. the penalty shrinks the paths, so this is a starting
# value for the MCMC rather than an estimate. only objects whose last dimension
# is `periods` are differenced
.lda_innov_sd <- function(paths, rho, periods) {
	if (!is.list(paths)) paths <- list(paths)
	innov <- numeric(0)
	for (path in paths) {
		if (is.null(path)) next
		d <- dim(path)
		if (is.null(d) || periods < 2L || d[length(d)] != periods) next
		M <- if (length(d) == 3L) matrix(path, ncol = d[3L]) else path
		innov <- c(
			innov,
			as.numeric(
				M[, -1L, drop = FALSE] -
					rho * M[, -ncol(M), drop = FALSE]
			)
		)
	}
	innov <- innov[is.finite(innov)]
	if (length(innov) < 2L) return(NULL)
	out <- stats::sd(innov)
	if (!is.finite(out) || out <= 0) return(NULL)
	out
}

.lda_as_fit <- function(
	prep,
	family,
	mode,
	symmetric,
	R,
	R_row,
	R_col,
	call,
	coef_path,
	beta_dyn,
	dynamic_beta_kind,
	a_mat,
	b_mat,
	Omat,
	Oarr,
	U,
	V,
	L,
	G_static,
	objective_trace,
	converged,
	iterations,
	lambda_ab,
	lambda_beta,
	lambda_uv,
	rho_ab,
	rho_beta,
	rho_uv,
	dynamic_ab,
	dynamic_uv,
	dynamic_uv_kind,
	lowrank_method,
	link = NULL,
	Zwork = NULL,
	obs_weights = NULL,
	irls_trace = NULL,
	stability = NULL,
	transition_u = NULL,
	transition_v = NULL,
	convergence_trace = NULL,
	tol = NULL,
	reweight_steps = NULL
) {
	# the level precision and the position ridge of a normal response carry the unit of its rows, as in the fit
	uv_unit_f <- if (!identical(family, "normal")) 1 else if (prep$bip) .lame_unit_var(prep$Yarr, family) else .lame_unit_sd(prep$Yarr, family)
	nr <- prep$nr
	nc <- prep$nc
	Tt <- prep$Tt
	if (is.null(Oarr)) Oarr <- .lda_o_array(Omat, Tt)
	objective_diff <- diff(objective_trace)
	objective_increase_tol <- 1e-10 * (abs(utils::head(objective_trace, -1L)) + 1)
	if (length(reweight_steps)) {
		keep_diff <- setdiff(seq_along(objective_diff), reweight_steps)
		objective_diff <- objective_diff[keep_diff]
		objective_increase_tol <- objective_increase_tol[keep_diff]
	}
	objective_increases <- if (length(objective_diff) > 0L) {
		sum(objective_diff > objective_increase_tol, na.rm = TRUE)
	} else {
		0L
	}
	max_objective_increase <- if (length(objective_diff) > 0L) {
		max(c(0, objective_diff), na.rm = TRUE)
	} else {
		0
	}
	objective_label <- if (identical(family, "normal")) {
		"penalized least-squares objective"
	} else {
		"penalized IRLS working objective"
	}
	# the non-normal working objective is comparable only within one IRLS working
	# problem, so the fit reports and stores the penalized deviance on the data
	# scale; accepted deviances are non-increasing and the running minimum ends on
	# the returned state
	pdev_trace <- if (!identical(family, "normal") && length(irls_trace) > 0L) {
		cummin(as.numeric(irls_trace))
	} else NULL
	reported_objective <- if (!is.null(pdev_trace)) {
		pdev_trace[length(pdev_trace)]
	} else {
		objective_trace[length(objective_trace)]
	}
	reported_label <- if (!is.null(pdev_trace)) {
		"penalized deviance"
	} else {
		objective_label
	}
	descent_history <- if (!is.null(pdev_trace)) pdev_trace else objective_trace
	Omat_store <- if (isTRUE(dynamic_uv)) {
		apply(Oarr, c(1, 2), mean)
	} else Omat
	dec <- .lda_decompose_node_paths(prep, symmetric, coef_path, a_mat, b_mat)
	coef_path <- dec$coef_path
	a_mat <- dec$a_mat
	b_mat <- dec$b_mat
	EZ <- .lda_eta(
		prep,
		coef_path,
		a_mat,
		b_mat,
		Oarr,
		symmetric = symmetric,
		beta_dyn = beta_dyn
	)
	fitted_list <- vector("list", Tt)
	resid_list <- vector("list", Tt)
	Zresid <- if (is.null(Zwork)) prep$Yarr else Zwork
	ve_num <- 0
	ve_den <- 0
	vew_num <- 0 # weighted (working-scale) numerator
	for (t in seq_len(Tt)) {
		fitted_list[[t]] <- .ae_link_inverse(EZ[[t]], family, link)
		if (!prep$bip && nr == nc) {
			diag(EZ[[t]]) <- NA_real_
			diag(fitted_list[[t]]) <- NA_real_
		}
		dimnames(EZ[[t]]) <- list(prep$row_names, prep$col_names)
		dimnames(fitted_list[[t]]) <- list(prep$row_names, prep$col_names)
		ok <- is.finite(Zresid[, , t])
		r <- Zresid[, , t] - EZ[[t]]
		resid_list[[t]] <- r
		ve_num <- ve_num + sum(r[ok]^2)
		ve_den <- ve_den + sum(ok)
		if (!is.null(obs_weights)) {
			wt <- obs_weights[, , t]
			okw <- ok & is.finite(wt)
			vew_num <- vew_num + sum(wt[okw] * r[okw]^2)
		}
	}
	names(EZ) <- prep$time_names
	names(fitted_list) <- prep$time_names
	names(resid_list) <- prep$time_names
	ve <- if (ve_den > 0) ve_num / ve_den else NA_real_
	pois_latent_var <- 0
	obs_mask <- is.finite(prep$Yarr)
	# dispersion of the weighted working problem of the additive block: the
	# weighted mean square (quasi-likelihood dispersion), since the IRLS weights
	# carry the heteroskedasticity. df_static counts the static mean-model blocks;
	# path df are added below
	df_static <- .ae_df_mean(
		prep$p_dyad,
		sum(rowSums(apply(obs_mask, c(1, 2), sum)) > 0),
		sum(colSums(apply(obs_mask, c(1, 2), sum)) > 0),
		max(R, R_row, R_col, 0L),
		symmetric
	)
	if (symmetric) df_static <- 2L * df_static
	# effective df of the penalized paths: each path's df is the trace of its
	# smoother, trace((K + P)^-1 K), K the per-period observation information and
	# P the transition penalty plus the level prior. the per-period identification
	# directions are subtracted on the same effective-period scale: one per margin
	# for the additive effects, and 2r + r^2 (centring and rotation of u_t v_t')
	# for the positions. without a penalty this is the per-period parameter count
	df_paths <- 0
	if (isTRUE(dynamic_ab) && Tt > 1L && !isTRUE(beta_dyn$any)) {
		level_prec_path <- if (!identical(family, "normal")) .lda_level_prec else 0
		transition_a_df <- .lda_transition_from_presence(prep$row_presence)
		transition_b_df <- if (symmetric) NULL else .lda_transition_from_presence(prep$col_presence)
		path_df <- function(den_mat, transition) {
			n <- nrow(den_mat)
			if (n == 0L) return(0)
			total <- 0
			for (i in seq_len(n)) {
				tw <- if (is.null(transition)) NULL else transition[i, ]
				Hi <- .lda_path_hessian(
					den_mat[i, ],
					lambda_ab,
					"ar1",
					rho_ab,
					tw,
					level_prec_path
				)
				Si <- tryCatch(solve(Hi), error = function(e) .lda_safe_pinv(Hi))
				total <- total + sum(diag(Si) * den_mat[i, ])
			}
			total
		}
		wgt <- if (is.null(obs_weights)) {
			array(as.numeric(obs_mask), dim(obs_mask))
		} else {
			obs_weights
		}
		den_a <- matrix(0, nr, Tt)
		den_b <- matrix(0, nc, Tt)
		for (t in seq_len(Tt)) {
			wt <- wgt[, , t]
			wt[!obs_mask[, , t]] <- 0
			wt[!is.finite(wt) | wt < 0] <- 0
			if (symmetric) {
				den_a[, t] <- rowSums(wt) + colSums(wt)
			} else {
				den_a[, t] <- rowSums(wt)
				den_b[, t] <- colSums(wt)
			}
		}
		tr_a <- path_df(den_a, transition_a_df)
		n_act_a <- max(sum(rowSums(den_a) > 0), 1)
		df_paths <- tr_a - tr_a / n_act_a
		if (!symmetric) {
			tr_b <- path_df(den_b, transition_b_df)
			n_act_b <- max(sum(rowSums(den_b) > 0), 1)
			df_paths <- df_paths + tr_b - tr_b / n_act_b
		}
		# the paths replace the static additive effects in df_static
		df_paths <- df_paths - if (symmetric) {
			max(sum(rowSums(den_a) > 0) - 1L, 0L)
		} else {
			max(sum(rowSums(den_a) > 0) - 1L, 0L) + max(sum(rowSums(den_b) > 0) - 1L, 0L)
		}
		# a symmetric fit pools the row-i and column-i cells into one path, so tr_a is
		# on the unique-pair scale while ve_den and df_static count each pair twice
		if (symmetric) df_paths <- 2 * df_paths
	}
	# regression paths share information with the additive paths: count their
	# joint constrained smoother in place of the static count
	if (isTRUE(beta_dyn$any)) {
		mean_weights <- if (is.null(obs_weights)) {
			array(as.numeric(obs_mask), dim(obs_mask))
		} else obs_weights
		mean_df <- .lda_joint_mean_df(
			prep,
			names(beta_dyn$mask),
			beta_dyn,
			mean_weights,
			dynamic_ab,
			lambda_ab,
			rho_ab,
			lambda_beta,
			dynamic_beta_kind,
			rho_beta,
			level_prec_ab = if (!identical(family, "normal") && dynamic_ab) .lda_level_prec else 0,
			symmetric = symmetric
		)
		base_mean_df <- .ae_df_mean(
			prep$p_dyad,
			sum(rowSums(apply(obs_mask, c(1, 2), sum)) > 0),
			sum(colSums(apply(obs_mask, c(1, 2), sum)) > 0),
			0L,
			symmetric
		)
		if (symmetric) base_mean_df <- 2 * base_mean_df
		df_paths <- mean_df - base_mean_df
	}
	# position paths (directed and bipartite, ar1 and t kinds): each actor's path
	# is the penalized solve of .lda_update_uv_path_block() given the other margin.
	# for the t kind the curvature at the fit adds 2 phi''(ss) (D'e)(D'e)' per
	# transition to the weighted quadratic, with phi''(ss) = -2 lambda^2 w^2 / (nu + r)
	R_uv <- if (is.array(U) && length(dim(U)) == 3L) dim(U)[2L] else 0L
	R_uv_col <- if (is.array(V) && length(dim(V)) == 3L) dim(V)[2L] else 0L
	G_uv <- if (prep$bip) G_static else NULL
	bip_uv_ok <- !prep$bip || (is.matrix(G_uv) && nrow(G_uv) == R_uv &&
			ncol(G_uv) == R_uv_col && all(is.finite(G_uv)))
	if (isTRUE(dynamic_uv) && dynamic_uv_kind %in% c("ar1", "t") && !symmetric &&
			bip_uv_ok && R_uv > 0L && R_uv_col > 0L && Tt > 1L) {
		wgt_uv <- if (is.null(obs_weights)) {
			array(as.numeric(obs_mask), dim(obs_mask))
		} else {
			obs_weights
		}
		t_kind <- identical(dynamic_uv_kind, "t")
		uv_trace <- function(other_cube, margin, transition, own_cube) {
			n_path <- if (identical(margin, "row")) nr else nc
			Rm <- dim(other_cube)[2L]
			total <- 0
			n_act <- 0L
			for (i in seq_len(n_path)) {
				tw <- if (is.null(transition)) rep(1, Tt) else transition[i, ]
				Kd <- matrix(0, Rm * Tt, Rm * Tt)
				for (t in seq_len(Tt)) {
					blk <- ((t - 1L) * Rm + 1L):(t * Rm)
					if (identical(margin, "row")) {
						ok <- obs_mask[i, , t]; wt <- wgt_uv[i, , t]
					} else {
						ok <- obs_mask[, i, t]; wt <- wgt_uv[, i, t]
					}
					keep <- ok & is.finite(wt) & wt > 0
					if (!any(keep)) next
					Xt <- matrix(other_cube[keep, , t], ncol = Rm)
					Kd[blk, blk] <- crossprod(Xt * sqrt(wt[keep]))
				}
				if (!any(Kd != 0)) next
				n_act <- n_act + 1L
				P <- kronecker(
					.lda_path_penalty_weighted(Tt, lambda_uv, tw, rho_uv),
					diag(Rm)
				)
				if (t_kind) {
					x <- matrix(own_cube[i, , ], Rm, Tt)
					for (t in 2:Tt) {
						w <- tw[t]
						if (!is.finite(w) || w <= 0) next
						e <- x[, t] - rho_uv * x[, t - 1L]
						g <- numeric(Rm * Tt)
						g[((t - 1L) * Rm + 1L):(t * Rm)] <- e
						g[((t - 2L) * Rm + 1L):((t - 1L) * Rm)] <- -rho_uv * e
						phi2 <- -2 * lambda_uv^2 * w^2 / (.lda_t_nu + Rm)
						P <- P + 2 * phi2 * tcrossprod(g)
					}
				}
				H <- Kd + P + diag((.lda_level_prec + 1e-6) * uv_unit_f, Rm * Tt)
				S <- tryCatch(solve(H, Kd), error = function(e) .lda_safe_pinv(H) %*% Kd)
				total <- total + sum(diag(S))
			}
			list(tr = total, n = n_act)
		}
		# row paths are solved against v_t g' and column paths against u_t g on a
		# bipartite fit, against v_t and u_t otherwise
		V_other <- V; U_other <- U
		if (prep$bip) {
			V_other <- array(0, dim = c(dim(V)[1L], R_uv, Tt))
			U_other <- array(0, dim = c(dim(U)[1L], R_uv_col, Tt))
			for (t in seq_len(Tt)) {
				V_other[, , t] <- matrix(V[, , t], ncol = R_uv_col) %*% t(G_uv)
				U_other[, , t] <- matrix(U[, , t], ncol = R_uv) %*% G_uv
			}
		}
		tu <- uv_trace(V_other, "row", transition_u, U)
		tv <- uv_trace(U_other, "col", transition_v, V)
		if (tu$n > 0L && tv$n > 0L) {
			t_eff <- (tu$tr + tv$tr) / (tu$n * R_uv + tv$n * R_uv_col)
			# per effective period: the centring of each margin and the directions the two
			# conditional traces share (rotation of u_t v_t', or the r_row x r_col span of
			# u_t k v_t' on a bipartite fit)
			df_uv <- if (prep$bip) {
				# the joint trace already leaves out the directions u_t a, v_t b with
				# a g b' = g, so only the centring of each margin is removed
				.lda_bip_joint_trace(
					U,
					V,
					G_uv,
					resid_list,
					obs_mask,
					wgt_uv,
					lambda_uv,
					rho_uv,
					transition_u,
					transition_v,
					row_unit = uv_unit_f
				) -
					(R_uv + R_uv_col) * t_eff
			} else {
				tu$tr + tv$tr - (R_uv + R_uv_col + R_uv * R_uv_col) * t_eff
			}
			# the paths replace the static rank-r count in df_static
			R_static <- max(R, R_row, R_col, 0L)
			df_paths <- df_paths + df_uv - R_static * (tu$n + tv$n - R_static)
		}
	}
	# position paths of a symmetric fit: one path per actor, solved by
	# .lda_update_dynamic_uv_symmetric() with a Hessian twice the unique-pair one,
	# so its trace is on the unique-pair scale. the shared per-period directions are
	# the centring (r) and the transformations that keep u_t s u_t' (r (r - 1) / 2),
	# r (r + 1) / 2 in all; ve_den and df_static count each pair twice, so the net
	# df enters twice
	if (isTRUE(dynamic_uv) && dynamic_uv_kind %in% c("ar1", "t") && symmetric &&
			!prep$bip && R_uv > 0L && Tt > 1L &&
			is.array(V) && length(dim(V)) == 3L) {
		wgt_uv <- if (is.null(obs_weights)) {
			array(as.numeric(obs_mask), dim(obs_mask))
		} else {
			obs_weights
		}
		t_kind <- identical(dynamic_uv_kind, "t")
		sgn <- .lda_symmetric_signs(U, V)
		tr_sym <- 0
		n_act <- 0L
		for (i in seq_len(nr)) {
			others <- seq_len(nr)[-i]
			tw_u <- if (is.null(transition_u)) rep(1, Tt) else transition_u[i, ]
			tw_v <- if (is.null(transition_v)) rep(1, Tt) else transition_v[i, ]
			Kd <- matrix(0, R_uv * Tt, R_uv * Tt)
			for (t in seq_len(Tt)) {
				blk <- ((t - 1L) * R_uv + 1L):(t * R_uv)
				w_row <- wgt_uv[i, others, t]
				w_col <- wgt_uv[others, i, t]
				w_row[!obs_mask[i, others, t] | !is.finite(w_row) | w_row < 0] <- 0
				w_col[!obs_mask[others, i, t] | !is.finite(w_col) | w_col < 0] <- 0
				w_tot <- w_row + w_col
				if (!any(w_tot > 0)) next
				Xt <- sweep(matrix(U[others, , t], ncol = R_uv), 2, sgn, "*")
				Kd[blk, blk] <- crossprod(Xt * sqrt(w_tot))
			}
			if (!any(Kd != 0)) next
			n_act <- n_act + 1L
			P <- kronecker(
				.lda_path_penalty_weighted(Tt, lambda_uv, tw_u, rho_uv) +
					.lda_path_penalty_weighted(Tt, lambda_uv, tw_v, rho_uv),
				diag(R_uv)
			)
			if (t_kind) {
				x <- matrix(U[i, , ], R_uv, Tt)
				for (t in 2:Tt) {
					e <- x[, t] - rho_uv * x[, t - 1L]
					g <- numeric(R_uv * Tt)
					g[((t - 1L) * R_uv + 1L):(t * R_uv)] <- e
					g[((t - 2L) * R_uv + 1L):((t - 1L) * R_uv)] <- -rho_uv * e
					for (w in c(tw_u[t], tw_v[t])) {
						if (!is.finite(w) || w <= 0) next
						phi2 <- -2 * lambda_uv^2 * w^2 / (.lda_t_nu + R_uv)
						P <- P + 2 * phi2 * tcrossprod(g)
					}
				}
			}
			H <- Kd + P + diag((2 * .lda_level_prec + 1e-6) * uv_unit_f, R_uv * Tt)
			S <- tryCatch(solve(H, Kd), error = function(e) .lda_safe_pinv(H) %*% Kd)
			tr_sym <- tr_sym + sum(diag(S))
		}
		if (n_act > 0L) {
			t_eff <- tr_sym / (n_act * R_uv)
			df_uv <- tr_sym - (R_uv * (R_uv + 1L) / 2) * t_eff
			R_static <- max(R, R_row, R_col, 0L)
			df_mult_static <- R_static * n_act - R_static * (R_static + 1L) %/% 2L
			df_paths <- df_paths + 2 * (df_uv - df_mult_static)
		}
	}
	df_mean <- df_static + df_paths
	ve_noise <- if (!is.null(obs_weights) && ve_den > 0) {
		vew_num / max(ve_den - df_static, 1)
	} else if (ve_den > 0) {
		ve_num / max(ve_den - df_static, 1)
	} else {
		ve
	}
	# binary IRLS fits report the latent-scale error variance fixed by the probit
	# identification, as the static ALS route and the MCMC do; the working
	# dispersion is kept as ve_working
	ve_working <- NULL
	if (identical(family, "binary") && !is.null(link)) {
		ve_working <- ve
		ve <- 1
	} else if (identical(family, "poisson") && identical(link, "log")) {
		# the IRLS fit solves the marginal log mean of the lognormal Poisson model
		# (y ~ Pois(exp(z)), z ~ N(eta, ve)); ve comes from the count overdispersion,
		# and the level moves to the latent scale by subtracting ve / 2
		ve_working <- ve
		ve_lat <- .ae_pois_latent_var(prep$Yarr, EZ, obs_mask, df_mean)
		icpt_row <- match("intercept", rownames(coef_path))
		if (is.finite(ve_lat) && !is.na(icpt_row)) {
			ve <- ve_lat
			pois_latent_var <- ve_lat
			coef_path[icpt_row, ] <- coef_path[icpt_row, ] - ve_lat / 2
			for (t in seq_len(Tt)) {
				EZ[[t]] <- EZ[[t]] - ve_lat / 2
				fitted_list[[t]] <- .ae_link_inverse(
					EZ[[t]],
					family,
					link,
					latent_var = ve_lat
				)
				dimnames(fitted_list[[t]]) <- list(prep$row_names, prep$col_names)
			}
		} else if (is.finite(ve_lat)) {
			ve <- ve_lat
		}
	} else {
		# normal family: residual mean square on the df left after the mean model,
		# counting penalized paths by their effective df
		ve <- if (ve_den > 0) ve_num / max(ve_den - df_mean, 1) else ve
	}

	dimnames(a_mat) <- list(prep$row_names, prep$time_names)
	dimnames(b_mat) <- list(prep$col_names, prep$time_names)
	a_bar <- .lda_presence_row_means(a_mat, prep$row_presence)
	b_bar <- .lda_presence_row_means(b_mat, prep$col_presence)
	# reported va / vb are variances of the actor levels, as in the MCMC VC: under
	# dynamic_ab they come from the cross-period covariance profile of the paths,
	# otherwise from the empirical variance minus the first-order sampling
	# variance. uncorrected values are kept as vc_raw and vc_period
	va_period <- if (nr > 1L) mean(.lda_presence_col_vars(a_mat, prep$row_presence), na.rm = TRUE) else NA_real_
	vb_period <- if (symmetric) NA_real_ else
	if (nc > 1L) mean(.lda_presence_col_vars(b_mat, prep$col_presence), na.rm = TRUE) else NA_real_
	va_raw <- if (nr > 1L) stats::var(a_bar) else NA_real_
	vb_raw <- if (symmetric) NA_real_ else if (nc > 1L) stats::var(b_bar) else NA_real_
	vc_noise <- list(va = NA_real_, vb = NA_real_)
	if (isTRUE(dynamic_ab)) {
		sides <- list(a = .ae_path_lag_profile(a_mat, prep$row_presence))
		if (!symmetric) sides$b <- .ae_path_lag_profile(b_mat, prep$col_presence)
		rho_within <- .ae_path_decay_rho(sides)
		if (!is.finite(rho_within)) rho_within <- rho_ab
		lev_a <- .ae_path_level_var(sides$a, rho_within)
		va <- if (is.finite(lev_a$va)) lev_a$va else va_raw
		vb <- if (symmetric) NA_real_ else {
			lev_b <- .ae_path_level_var(sides$b, rho_within)
			if (is.finite(lev_b$va)) lev_b$va else vb_raw
		}
	} else {
		vc_noise <- .ae_effect_noise_var(
			obs_mask,
			obs_weights,
			ve_noise,
			symmetric,
			as.numeric(rowSums(apply(obs_mask, c(1, 2), sum)) > 0),
			as.numeric(colSums(apply(obs_mask, c(1, 2), sum)) > 0)
		)
		va <- if (is.finite(va_raw) && is.finite(vc_noise$va)) max(va_raw - vc_noise$va, 0) else va_raw
		vb <- if (is.finite(vb_raw) && is.finite(vc_noise$vb)) max(vb_raw - vc_noise$vb, 0) else vb_raw
	}
	cab <- if (symmetric || prep$bip) NA_real_ else
		mean(
			.lda_presence_col_cov(a_mat, b_mat, prep$row_presence, prep$col_presence),
			na.rm = TRUE
		)
	rho <- 0
	rho_path <- NULL
	if (!prep$bip && !symmetric) {
		rho_path <- rep(NA_real_, Tt)
		names(rho_path) <- prep$time_names
		e_ij <- numeric(0)
		e_ji <- numeric(0)
		for (t in seq_len(Tt)) {
			rt <- resid_list[[t]]
			ut <- upper.tri(rt)
			e_t <- rt[ut]
			et_t <- t(rt)[ut]
			ok_t <- is.finite(e_t) & is.finite(et_t)
			if (sum(ok_t) > 2L) {
				rho_path[t] <- suppressWarnings(stats::cor(e_t[ok_t], et_t[ok_t]))
				if (!is.finite(rho_path[t])) rho_path[t] <- NA_real_
			}
			e_ij <- c(e_ij, e_t)
			e_ji <- c(e_ji, et_t)
		}
		ok <- is.finite(e_ij) & is.finite(e_ji)
		if (sum(ok) > 2L) {
			rho <- suppressWarnings(stats::cor(e_ij[ok], e_ji[ok]))
			if (!is.finite(rho)) rho <- 0
		}
	}
	VC <- c(va = va, cab = cab, vb = vb, rho = rho, ve = ve)

	coef_names <- rownames(coef_path)
	coef_overall <- rowMeans(coef_path)
	names(coef_overall) <- coef_names
	beta_mask <- beta_dyn$mask
	if (length(beta_mask) < length(coef_names)) {
		beta_mask <- c(beta_mask, rep(FALSE, length(coef_names) - length(beta_mask)))
	}
	beta_mask <- beta_mask[seq_along(coef_names)]
	names(beta_mask) <- coef_names
	beta_groups <- beta_dyn$groups
	if (length(beta_groups) < length(coef_names)) {
		beta_groups <- c(
			beta_groups,
			rep("", length(coef_names) - length(beta_groups))
		)
	}
	beta_groups <- beta_groups[seq_along(coef_names)]
	names(beta_groups) <- coef_names
	BETA <- if (isTRUE(beta_dyn$any)) {
		array(
			coef_path,
			dim = c(1L, nrow(coef_path), ncol(coef_path)),
			dimnames = list("point", coef_names, prep$time_names)
		)
	} else {
		matrix(
			coef_path[, 1L],
			nrow = 1L,
			dimnames = list("point", coef_names)
		)
	}

	out <- list(
		call = call,
		family = family,
		mode = mode,
		symmetric = symmetric,
		R = R,
		R_row = if (prep$bip) R_row else NULL,
		R_col = if (prep$bip) R_col else NULL,
		longitudinal = TRUE,
		mu = coef_overall["intercept"],
		beta = coef_overall[-1L],
		coef_path = coef_path,
		beta_dynamic_per_t = coef_path,
		beta_dynamic_mask = beta_mask,
		beta_dynamic_groups = beta_groups,
		BETA = BETA,
		a = a_bar,
		b = b_bar,
		a_dynamic = a_mat,
		b_dynamic = b_mat,
		APM = a_bar,
		BPM = b_bar,
		U = U,
		V = V,
		L = L,
		G = if (prep$bip && R_row > 0L && R_col > 0L) G_static else NULL,
		Omat = Omat_store,
		Oarr = Oarr,
		UVPM = if (isTRUE(dynamic_uv)) Oarr else Omat_store,
		coefficients = coef_overall,
		VC = VC,
		s2 = ve,
		ve_working = ve_working,
		# additive variances before the level and noise correction, and the mean
		# per-period variance of the fitted paths
		vc_raw = c(va = va_raw, vb = vb_raw),
		vc_period = c(va = va_period, vb = vb_period),
		vc_noise = c(va = vc_noise$va, vb = vc_noise$vb),
		df_mean = df_mean,
		pois_latent_var = if (pois_latent_var > 0) pois_latent_var else NULL,
		EZ = EZ,
		fitted = fitted_list,
		YPM = fitted_list,
		residuals = resid_list,
		iterations = iterations,
		converged = converged,
		deviance = reported_objective,
		objective = reported_objective,
		objective_label = reported_label,
		# one entry per accepted working problem for the IRLS families, one per block
		# sweep for the normal family; monotone either way
		dev_history = descent_history,
		objective_history = objective_trace,
		penalized_deviance = pdev_trace,
		n_time = Tt,
		dims = list(n_row = nr, n_col = nc, p = length(coef_overall)),
		row_names = prep$row_names,
		col_names = prep$col_names,
		x_names = coef_names[-1L],
		time_names = prep$time_names,
		Y = prep$Yarr,
		X = prep$X,
		W_row = prep$W_row,
		W_col = prep$W_col,
		W_row_arr = prep$W_row_arr,
		W_col_arr = prep$W_col_arr,
		changing_composition = isTRUE(prep$changing_composition),
		row_presence = prep$row_presence,
		col_presence = prep$col_presence,
		actor_presence = prep$actor_presence,
		dynamic_uv = isTRUE(dynamic_uv),
		dynamic_uv_kind = if (isTRUE(dynamic_uv)) dynamic_uv_kind else NULL,
		dynamic_ab = isTRUE(dynamic_ab),
		dynamic_beta = isTRUE(beta_dyn$any),
		# the estimator and simulate.ame_als use the pooled VC["rho"]; rho_path is a
		# labelled diagnostic
		dynamic_rho = FALSE,
		rho_path = rho_path,
		residual_rho_path = rho_path,
		rho_path_role = if (is.null(rho_path)) NULL else "residual_diagnostic",
		rho_uv = if (isTRUE(dynamic_uv)) rho_uv else NULL,
		sigma_uv = if (isTRUE(dynamic_uv)) .lda_innov_sd(list(U, V), rho_uv, Tt) else NULL,
		rho_ab = if (isTRUE(dynamic_ab)) rho_ab else NULL,
		sigma_ab = if (isTRUE(dynamic_ab)) {
			.lda_innov_sd(if (isTRUE(symmetric)) list(a_mat) else list(a_mat, b_mat), rho_ab, Tt)
		} else NULL,
		rho_beta = if (isTRUE(beta_dyn$any)) rho_beta else NULL,
		# only evolving coefficients carry innovations; a static row would add
		# (1 - rho) times its level to the scale
		sigma_beta = if (isTRUE(beta_dyn$any) && any(beta_mask) &&
				nrow(coef_path) == length(beta_mask)) {
			.lda_innov_sd(coef_path[beta_mask, , drop = FALSE], rho_beta, Tt)
		} else NULL,
		lambda_u = if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) transition_u else NULL,
		lambda_v = if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) transition_v else NULL,
		lambda_uv_als = if (isTRUE(dynamic_uv)) lambda_uv else NULL,
		lambda_ab_als = if (isTRUE(dynamic_ab)) lambda_ab else NULL,
		lambda_beta_als = if (isTRUE(beta_dyn$any)) lambda_beta else NULL,
		lowrank_method = lowrank_method,
		non_normal_method = if (identical(family, "normal")) "none" else "irls",
		link = link,
		obs_weights = obs_weights,
		irls_trace = irls_trace,
		stability = stability,
		dynamic = list(
			uv = isTRUE(dynamic_uv),
			uv_kind = if (isTRUE(dynamic_uv)) dynamic_uv_kind else NULL,
			ab = isTRUE(dynamic_ab),
			ab_kind = "ar1_penalty",
			beta = isTRUE(beta_dyn$any),
			beta_kind = dynamic_beta_kind,
			snap_model_estimated = FALSE,
			t_model_estimated = isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t"),
			mode = mode
		),
		diagnostics = list(
			converged = converged,
			objective = objective_trace,
			convergence_trace = convergence_trace,
			objective_label = objective_label,
			# for the IRLS families the trace is comparable only within one working
			# problem, so reweight boundaries are excluded from the increase count
			objective_scope = if (identical(family, "normal")) {
				"whole fit"
			} else {
				"within one IRLS working problem"
			},
			reweight_steps = reweight_steps,
			penalized_deviance = pdev_trace,
			reported_objective = reported_objective,
			reported_objective_label = reported_label,
			objective_increases = objective_increases,
			max_objective_increase = max_objective_increase,
			last_objective_delta = if (length(objective_trace) > 1L) {
				utils::tail(diff(objective_trace), 1L)
			} else NA_real_,
			last_relative_objective_delta = if (length(objective_trace) > 1L) {
				abs(utils::tail(diff(objective_trace), 1L)) /
					max(1, abs(objective_trace[length(objective_trace) - 1L]))
			} else NA_real_,
			max_objective_delta = if (length(objective_trace) > 1L)
				max(abs(diff(objective_trace))) else NA_real_,
			last_max_fitted_delta = if (!is.null(convergence_trace) &&
					nrow(convergence_trace) > 0L) {
				utils::tail(convergence_trace$max_fitted_delta, 1L)
			} else NA_real_,
			last_max_parameter_delta = if (!is.null(convergence_trace) &&
					nrow(convergence_trace) > 0L) {
				utils::tail(convergence_trace$max_parameter_delta, 1L)
			} else NA_real_,
			convergence_tolerance = tol,
			fitted_value_tolerance = if (!is.null(tol)) sqrt(tol) else NA_real_,
			fitted_value_delta_scope = "observed_cells",
			message = if (converged) "Converged" else "Stopped at max_iter before the convergence criteria were reached."
		),
		meta = list(
			sampler = "dynamic_als",
			algorithm = "dynamic penalized least-squares block coordinate descent",
			uncertainty_available = FALSE,
			approximation_note = paste0(
				if (identical(family, "normal")) "Gaussian" else "IRLS",
				" dynamic ALS/MAP point estimate with RW/AR smoothing penalties. ",
				"It provides fitted paths and diagnostics, not posterior draws."
			)
		)
	)
	class(out) <- c("lame_dynamic_als", "lame_als", "ame_als")
	.ae_update_response_rho(out, link)
}

lame_dynamic_als <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	R = 0,
	R_row = NULL,
	R_col = NULL,
	family = "normal",
	mode = c("unipartite", "bipartite"),
	symmetric = FALSE,
	dynamic_ab = FALSE,
	dynamic_uv = FALSE,
	dynamic_uv_kind = c("ar1", "snap", "t"),
	dynamic_beta = FALSE,
	dynamic_beta_kind = c("ar1", "rw1", "rw2", "matern32"),
	dynamic_G = FALSE,
	prior = list(),
	lambda_ab = NULL,
	lambda_beta = NULL,
	lambda_uv = NULL,
	max_iter = 200,
	tol = 5e-5,
	link = c("probit", "logit"),
	stability = c("none", "quick", "validation"),
	lowrank_method = c("hybrid", "mm", "als"),
	verbose = TRUE,
	seed = 6886,
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	call = NULL,
	start = NULL,
	start_jitter = 0,
	bias_correct = NULL,
	rvar = NULL,
	cvar = NULL,
	nvar = NULL,
	dcor = NULL,
	dynamic_beta_pool = NULL,
	...
) {
	bias_correct <- .ae_check_bias_correct(bias_correct)
	lowrank_user <- !missing(lowrank_method)
	mode <- match.arg(mode)
	dynamic_uv_kind <- match.arg(dynamic_uv_kind)
	dynamic_beta_kind <- match.arg(dynamic_beta_kind)
	link <- match.arg(link)
	stability <- match.arg(stability)
	lowrank_method <- match.arg(lowrank_method)
	bootstrap_type <- match.arg(bootstrap_type)
	# the sampler's switches (NULL when unset): the EB estimator honours rvar,
	# cvar, nvar, dcor and dynamic_beta_pool, and a route that cannot honour a
	# switch says so
	for (sw in c("rvar", "cvar", "nvar", "dcor")) {
		v <- get(sw)
		if (!is.null(v) && !(is.logical(v) && length(v) == 1L && !is.na(v))) {
			cli::cli_abort("{.arg {sw}} must be {.code TRUE} or {.code FALSE}.")
		}
	}
	if (!is.null(dynamic_beta_pool)) {
		dynamic_beta_pool <- match.arg(dynamic_beta_pool, c("none", "rho", "sigma", "both"))
	}
	.lda_warn_ignored_switches <- function(names_) {
		if (!length(names_)) return(invisible(NULL))
		ignored_text <- paste(names_, collapse = ", ")
		cli::cli_warn(c(
			"{.fn lame} {.arg method} = {.val als}: ignored MCMC-only arguments on the dynamic ALS route: {ignored_text}.",
			"i" = "Use {.code method = \"mcmc\"} if you need those sampler controls."
		))
	}
	# nothing is read from `...`, so an argument there was meant for a control
	# this route does not have; say so
	dots <- list(...)
	if (length(dots) > 0L) {
		dot_names <- names(dots)
		if (is.null(dot_names)) dot_names <- rep("", length(dots))
		dot_names[!nzchar(dot_names)] <- "<unnamed>"
		static_only <- intersect(
			dot_names,
			c(
				"linear_solver", "multistart",
				"non_normal_method"
			)
		)
		hint <- if (length(static_only) > 0L) {
			"{.arg {static_only}} {?is/are} a static ALS control; the dynamic route solves whole-path blocks with its own solver."
		} else {
			"Check the argument names against {.fn lame_dynamic_als}."
		}
		cli::cli_abort(c(
			"{.fn lame_dynamic_als} does not accept {.arg {dot_names}}.",
			"i" = hint,
			"i" = "Use {.fn lame} with {.code method = \"als\"} for the documented dynamic ALS entry point."
		))
	}
	# seed locally and restore the global RNG stream on exit
	.had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	.old_seed <- if (.had_seed) get(".Random.seed", envir = globalenv()) else NULL
	on.exit({
		if (.had_seed) {
			assign(".Random.seed", .old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)
	amen_aliases <- c(nrm = "normal", bin = "binary", pois = "poisson")
	if (length(family) == 1L && is.character(family) &&
			family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
			"i" = "Update your script to {.code family = \"{new_family}\"} when convenient."
		))
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
			!family %in% c("normal", "binary", "poisson")) {
		cli::cli_abort(c(
			"{.fn lame_dynamic_als} currently supports {.val normal}, {.val binary}, and {.val poisson}.",
			"i" = "Use {.code method = \"mcmc\"} for dynamic {.val {family}} models."
		))
	}
	# Poisson responses must be non-negative integer counts, as in ame()/lame()
	if (identical(family, "poisson")) .lame_check_poisson_counts(Y)
	.lame_als_check_covariate_scale(Xdyad, Xrow, Xcol)
	if (identical(bias_correct, "eb") && !family %in% c("binary", "poisson", "normal")) {
		cli::cli_abort(
			c(
				"The empirical-Bayes estimator ({.code bias_correct = \"eb\"}) is defined for the normal, binary and Poisson likelihoods; {.arg family} is {.val {family}}.",
				"i" = "Dynamic {.val {family}} ALS fits take {.code bias_correct = NULL} or {.val none}."
			),
			class = "lame_als_eb"
		)
	}
	if (identical(bias_correct, "bootstrap")) {
		cli::cli_abort(
			c(
				"The parametric-bootstrap bias correction is available for static binary ALS fits only.",
				"i" = "Dynamic binary ALS fits take {.code bias_correct = NULL} (the empirical-Bayes estimator) or {.val none}; dynamic Poisson fits take {.code bias_correct = NULL} (the empirical-Bayes estimator), {.val jackknife} (directed panels) or {.val none}."
			),
			class = "lame_jk_split"
		)
	}
	if (identical(bias_correct, "jackknife") &&
			(!identical(family, "poisson") || !identical(mode, "unipartite") ||
					isTRUE(symmetric) || (isTRUE(dynamic_uv) && dynamic_uv_kind != "ar1"))) {
		cli::cli_abort(
			c(
				"The split-sample jackknife is not available for this dynamic ALS specification.",
				"i" = "The dynamic correction estimates directed Poisson variance and reciprocal covariance."
			),
			class = "lame_jk_split"
		)
	}
	if (identical(family, "poisson")) link <- "log"
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "snap")) {
		cli::cli_abort(c(
			"{.fn lame_dynamic_als} supports smooth {.arg dynamic_uv_kind} values {.val ar1} and {.val t}.",
			"i" = "The snap kind has no ALS route; use {.code dynamic_uv_kind = \"snap\"} with {.code method = \"mcmc\"} for directed panels, or {.val ar1} here."
		))
	}
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") &&
			identical(mode, "bipartite")) {
		.lda_refuse_bipartite_t(caller = "lame_dynamic_als")
	}
	# as in the MCMC route: a time-varying interaction matrix is not identified
	# apart from the latent positions
	.lame_refuse_dynamic_G(dynamic_G)
	if (identical(dynamic_beta_kind, "matern32")) {
		cli::cli_abort(c(
			"{.arg dynamic_beta_kind = \"matern32\"} is not yet supported by {.fn lame_dynamic_als}.",
			"i" = "Use {.val ar1}, {.val rw1}, or {.val rw2}, or fit with {.code method = \"mcmc\"}."
		))
	}
	if (!is.numeric(bootstrap) || bootstrap < 0L ||
			!isTRUE(all.equal(bootstrap, round(bootstrap)))) {
		cli::cli_abort("{.arg bootstrap} must be a non-negative integer.")
	}
	bootstrap <- as.integer(bootstrap)
	# bootstrap replicates of a dynamic non-normal fit are not calibrated, so
	# the route refuses them
	if (bootstrap > 0L && family %in% c("binary", "poisson") &&
			(is.null(bias_correct) || identical(bias_correct, "eb"))) {
		cli::cli_abort(
			c(
				"A bootstrap is not offered for the dynamic empirical-Bayes estimator.",
				"i" = "It carries Laplace standard errors: {.code vcov(fit)}, {.code confint(fit)}, {.code summary(fit)} and {.code fit$eb$se}.",
				"i" = "For posterior summaries of the paths and variance components use {.fn lame} with {.code method = \"mcmc\"}, optionally warm-started from this fit with {.code start_vals = als_start_vals(fit)}."
			),
			class = "lame_dynamic_eb_unsupported"
		)
	}
	if (bootstrap > 0L && !identical(family, "normal")) {
		fam <- family
		cli::cli_abort(
			c(
				"A bootstrap is not available for a {.val {fam}} dynamic ALS fit.",
				"i" = "Each replicate refits through the IRLS reweighting loop and lands at a different local optimum: on a 60-actor binary {.arg dynamic_uv} panel the replicate standard deviations came out 8 to 11 times the MCMC posterior standard deviation, with 14% of replicates not converging.",
				"i" = "For uncertainty on this fit use MCMC warm-started from it: {.code lame(..., start_vals = als_start_vals(fit))}.",
				"i" = "The parametric bootstrap is available for normal-family dynamic fits and for every static ALS fit."
			),
			class = "lame_als_bootstrap_unsupported"
		)
	}
	# the penalised normal fit runs the bootstrap: refuse a count below the
	# minimum before fitting
	if (bootstrap > 0L && identical(bias_correct, "none")) .ae_check_bootstrap_n(bootstrap)
	if (!is.numeric(start_jitter) || length(start_jitter) != 1L ||
			!is.finite(start_jitter) || start_jitter < 0) {
		cli::cli_abort("{.arg start_jitter} must be a non-negative finite number.")
	}

	if (is.null(call)) call <- match.call()
	call_args <- list(
		Y = Y,
		Xdyad = Xdyad,
		Xrow = Xrow,
		Xcol = Xcol,
		R = R,
		R_row = R_row,
		R_col = R_col,
		family = family,
		mode = mode,
		symmetric = symmetric,
		dynamic_ab = dynamic_ab,
		dynamic_uv = dynamic_uv,
		dynamic_uv_kind = dynamic_uv_kind,
		dynamic_beta = dynamic_beta,
		dynamic_beta_kind = dynamic_beta_kind,
		prior = prior,
		lambda_ab = lambda_ab,
		lambda_beta = lambda_beta,
		lambda_uv = lambda_uv,
		max_iter = max_iter,
		tol = tol,
		link = link,
		stability = "none",
		bias_correct = bias_correct,
		rvar = rvar,
		cvar = cvar,
		nvar = nvar,
		dcor = dcor,
		dynamic_beta_pool = dynamic_beta_pool,
		lowrank_method = lowrank_method,
		bootstrap_type = bootstrap_type,
		bootstrap_block_length = bootstrap_block_length,
		call = call
	)
	prep <- .ae_prepare(
		Y,
		Xdyad,
		Xrow,
		Xcol,
		mode,
		longitudinal = TRUE,
		node_time_policy = "quiet",
		symmetric = isTRUE(symmetric)
	)
	if (mode == "unipartite" && prep$nr != prep$nc) {
		cli::cli_abort("Unipartite dynamic ALS requires square network slices.")
	}
	if (symmetric && prep$bip) {
		cli::cli_warn("{.arg symmetric} is ignored for bipartite networks.")
		symmetric <- FALSE
	}
	if (symmetric && prep$nr != prep$nc) {
		cli::cli_abort("Symmetric dynamic ALS requires a square network.")
	}
	if (symmetric && prep$nr == prep$nc) {
		for (t in seq_len(prep$Tt)) {
			Yt <- prep$Yarr[, , t]
			off <- Yt - t(Yt)
			if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
				cli::cli_abort(c(
					"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
					"i" = "Symmetrize the input or set {.arg symmetric} = FALSE."
				))
			}
		}
	}
	if (prep$Tt < 2L) {
		cli::cli_abort("{.fn lame_dynamic_als} requires at least 2 time periods.")
	}
	if (length(R) != 1L || anyNA(R) || R < 0 ||
			!isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	R <- as.integer(R)
	.validate_rank <- function(x, nm) {
		if (length(x) != 1L || anyNA(x) || x < 0 ||
				!isTRUE(all.equal(x, round(x)))) {
			cli::cli_abort("{.arg {nm}} must be a single non-negative integer.")
		}
		as.integer(x)
	}
	if (prep$bip) {
		R_row_eff <- .validate_rank(.lda_default(R_row, R), "R_row")
		R_col_eff <- .validate_rank(.lda_default(R_col, R), "R_col")
	} else {
		if (!is.null(R_row) || !is.null(R_col)) {
			cli::cli_warn("{.arg R_row} and {.arg R_col} are ignored for unipartite dynamic ALS.")
		}
		R_row_eff <- R
		R_col_eff <- R
	}
	R_eff <- if (prep$bip) max(R_row_eff, R_col_eff) else R
	if (isTRUE(dynamic_uv) &&
			(R_row_eff <= 0L || R_col_eff <= 0L)) {
		cli::cli_abort(c(
			"Dynamic multiplicative ALS requires a positive multiplicative latent rank.",
			"i" = if (prep$bip) {
				"Set positive {.arg R_row} and {.arg R_col} values, or set {.arg R} > 0 for both."
			} else {
				"Set {.arg R} to at least 1, or turn off {.arg dynamic_uv}."
			}
		))
	}
	if (prep$bip) {
		if (R_row_eff >= prep$nr) {
			cli::cli_abort("{.arg R_row} = {R_row_eff} must be smaller than n_row = {prep$nr}.")
		}
		if (R_col_eff >= prep$nc) {
			cli::cli_abort("{.arg R_col} = {R_col_eff} must be smaller than n_col = {prep$nc}.")
		}
	} else {
		rmax <- min(prep$nr, prep$nc)
		if (R >= rmax) {
			cli::cli_abort("{.arg R} = {R} must be smaller than min(n_row, n_col) = {rmax}.")
		}
	}
	if (symmetric && prep$p_dyad > 0L) {
		for (k in seq_len(prep$p_dyad)) for (t in seq_len(prep$Tt)) {
			prep$X[, , k, t] <- .ae_symmetrize(prep$X[, , k, t])
		}
	}
	if (symmetric && lowrank_method != "mm" && R > 0L) {
		# the symmetric eigenmodel takes the majorise-minimise solver silently; only
		# an explicit choice of an unavailable solver warns
		if (isTRUE(lowrank_user)) {
			cli::cli_warn(c(
				"{.arg lowrank_method} = {.val {lowrank_method}} is not available for symmetric models.",
				"i" = "The symmetric eigenmodel uses the majorise-minimise solver."
			))
		}
		lowrank_method <- "mm"
	}
	# stability restarts and bootstrap refits replay call_args, so record the
	# solver this fit used
	call_args$lowrank_method <- lowrank_method

	base_layout <- .lda_coef_layout(prep, include_nodes = FALSE)
	full_layout <- .lda_coef_layout(prep, include_nodes = TRUE)
	dynamic_beta_parse <- dynamic_beta
	if (is.logical(dynamic_beta) && length(dynamic_beta) > 1L &&
			length(dynamic_beta) == length(base_layout$names) &&
			length(full_layout$names) > length(base_layout$names)) {
		dynamic_beta_parse <- c(
			dynamic_beta,
			rep(
				FALSE,
				length(full_layout$names) -
					length(base_layout$names)
			)
		)
	}
	beta_dyn_full <- parse_dynamic_beta(
		dynamic_beta_parse,
		full_layout$names,
		full_layout$block,
		intercept = TRUE,
		Tn = prep$Tt,
		family = family,
		mode = mode
	)
	node_dynamic <- beta_dyn_full$mask &
		full_layout$block %in% c("row", "col")
	if (any(node_dynamic) && "intercept" %in% full_layout$names) {
		beta_dyn_full$mask[match("intercept", full_layout$names)] <- TRUE
		beta_dyn_full <- .lda_beta_dyn_from_mask(
			beta_dyn_full$mask,
			full_layout$names,
			full_layout$block
		)
	}
	base_idx <- which(!full_layout$block %in% c("row", "col"))
	node_idx <- which(full_layout$block %in% c("row", "col"))
	include_idx <- sort(unique(c(base_idx, node_idx[beta_dyn_full$mask[node_idx]])))
	coef_names <- full_layout$names[include_idx]
	coef_block <- full_layout$block[include_idx]
	beta_dyn <- .lda_beta_dyn_from_mask(
		beta_dyn_full$mask[include_idx],
		coef_names,
		coef_block
	)
	# refuse exactly collinear design columns (node covariates included) for
	# both the EB and the penalised estimators
	.lda_check_design_rank(prep, full_layout$names, symmetric = isTRUE(symmetric))
	# the smoothness weight of a moving coefficient is per unit of its covariate: lambda_beta times its mean square
	if (beta_dyn$any) beta_dyn$ms <- .lda_design_ms(prep, coef_names, isTRUE(symmetric))
	# binary, Poisson and normal dynamic fits are empirical-Bayes fits by default
	# (lame_dynamic_als_eb.R), treating the paths as random to remove the
	# incidental-parameter bias; bias_correct = "none" fits the penalised model
	# below, and unsupported specifications are refused by name
	pln_eb <- identical(family, "poisson")
	nrm_eb <- identical(family, "normal")
	if ((identical(family, "binary") || pln_eb || nrm_eb) &&
			(is.null(bias_correct) || identical(bias_correct, "eb")) &&
			(isTRUE(dynamic_ab) || isTRUE(beta_dyn$any) || isTRUE(dynamic_uv))) {
		if (nrm_eb) {
			# the Gaussian cell: identity link (the penalised fit's link argument is unused)
			link <- "identity"
			if (bootstrap > 0L) {
				cli::cli_abort(
					c(
						"A bootstrap is not offered for the dynamic empirical-Bayes estimator.",
						"i" = "It carries the exact Gaussian standard errors given the hyperparameters: {.code vcov(fit)}, {.code confint(fit)}, {.code summary(fit)} and {.code fit$eb$se}.",
						"i" = "The parametric bootstrap of the penalised normal fit is {.code bias_correct = \"none\"}; for posterior summaries of the paths and variance components use {.fn lame} with {.code method = \"mcmc\"}, optionally warm-started from this fit with {.code start_vals = als_start_vals(fit)}."
					),
					class = "lame_dynamic_eb_unsupported"
				)
			}
		}
		# the paths use period-specific node values, so the static-mean warnings
		# below do not apply
		.lame_inform_cov_na(prep$cov_na$n)
		.lde_check_supported(
			prep,
			symmetric,
			R_eff,
			R_row_eff,
			R_col_eff,
			dynamic_ab,
			dynamic_uv,
			dynamic_uv_kind,
			beta_dyn,
			dynamic_beta_kind,
			lambda_ab,
			lambda_beta,
			lambda_uv,
			prior,
			link
		)
		if (isTRUE(lowrank_user) && R_eff > 0L) {
			cli::cli_warn(c(
				"{.arg lowrank_method} = {.val {lowrank_method}} does not apply to the empirical-Bayes estimator.",
				"i" = "It solves the latent factors jointly with the other effects by Gauss-Newton; {.code bias_correct = \"none\"} uses the alternating low-rank solvers."
			))
		}
		# a bipartite interaction U G V' has rank min(R_row, R_col); the sampler
		# takes both ranks or neither
		R_eb <- R_eff
		if (isTRUE(prep$bip)) {
			if (xor(R_row_eff > 0L, R_col_eff > 0L)) {
				cli::cli_abort(c(
					"Only one of {.arg R_row} ({R_row_eff}) / {.arg R_col} ({R_col_eff}) is positive.",
					"i" = "Set both {.arg R_row} and {.arg R_col} (or pass {.arg R}) to include latent factors."
				))
			}
			R_eb <- min(R_row_eff, R_col_eff)
			if (R_row_eff != R_col_eff) {
				cli::cli_inform(c(
					"i" = "{.arg R_row} = {R_row_eff} and {.arg R_col} = {R_col_eff}: the interaction U G V' has rank {R_eb}, which is what the empirical-Bayes estimator fits; the factor of the larger rank carries zero columns beyond it."
				))
			}
		}
		# a directed network's dyadic correlation is the sampler's dcor: a dyad's two
		# cells form one bivariate probit observation with one correlation for all
		# periods (lame_dynamic_als_eb_pair.R); undirected and bipartite networks have
		# no cell pairs
		dcor_eb <- if (is.null(dcor)) !symmetric else isTRUE(dcor)
		dcor_eb <- dcor_eb && !symmetric && !isTRUE(prep$bip)
		if (dcor_eb && identical(link, "logit") && !pln_eb) {
			if (!is.null(dcor)) {
				cli::cli_abort(
					c(
						"{.code dcor = TRUE} needs the probit link.",
						"i" = "The dyadic correlation is the correlation of the latent normal errors of the two cells of a dyad, as in {.fn lame}; use {.code link = \"probit\"} or {.code dcor = FALSE}."
					),
					class = c("lame_dynamic_eb_unsupported", "lame_als_eb")
				)
			}
			dcor_eb <- FALSE
		}
		# an undirected network has one effect per actor, switched by nvar (the
		# sampler overrides rvar and cvar with it)
		rvar_eff <- if (symmetric) isTRUE(.lda_default(nvar, TRUE)) else isTRUE(.lda_default(rvar, TRUE))
		cvar_eff <- if (symmetric) rvar_eff else isTRUE(.lda_default(cvar, TRUE))
		# the moving coefficients over the whole design, node covariates included
		beta_mask <- if (isTRUE(beta_dyn$any)) {
			parse_dynamic_beta(
				dynamic_beta_parse,
				full_layout$names,
				full_layout$block,
				intercept = TRUE,
				Tn = prep$Tt,
				family = family,
				mode = mode
			)$mask
		} else NULL
		# a replay of this call (ame_als_refit, stability restarts) must not read as a
		# user choice of solver, which the EB estimator does not use
		if (!isTRUE(lowrank_user)) call_args$lowrank_method <- NULL
		# EB iteration budget of lame(method = "als"), 1000 passes, unless set; replays
		# keep the budget the fit used
		eb_budget <- if (missing(max_iter)) 1000L else max_iter
		call_args$max_iter <- eb_budget
		return(.lde_dynamic_fit(
			prep,
			family,
			mode,
			symmetric,
			link,
			prior,
			eb_budget,
			tol,
			verbose,
			seed,
			call,
			call_args,
			bias_requested = bias_correct,
			start = start,
			start_jitter = start_jitter,
			stability = stability,
			R = R_eb,
			beta_mask = beta_mask,
			beta_kind = dynamic_beta_kind,
			dynamic_ab = isTRUE(dynamic_ab),
			rvar = rvar_eff,
			cvar = cvar_eff,
			R_row = if (isTRUE(prep$bip)) R_row_eff,
			R_col = if (isTRUE(prep$bip)) R_col_eff,
			beta_pool = .lda_default(dynamic_beta_pool, "none"),
			dynamic_uv = isTRUE(dynamic_uv),
			dynamic_uv_kind = dynamic_uv_kind,
			dcor = dcor_eb
		))
	}
	# the penalised routes honour none of these switches
	.lda_warn_ignored_switches(names(Filter(
		Negate(is.null),
		list(
			rvar = rvar,
			cvar = cvar,
			nvar = nvar,
			dcor = dcor,
			dynamic_beta_pool = dynamic_beta_pool
		)
	)))

	dyn_node_names <- .lda_dynamic_node_names(beta_dyn)
	vary_row_static <- setdiff(
		.lda_node_varies_names(prep, "row"),
		dyn_node_names
	)
	vary_col_static <- setdiff(
		.lda_node_varies_names(prep, "col"),
		dyn_node_names
	)
	if (length(vary_row_static) > 0L) {
		cli::cli_warn(c(
			"Some row node covariates vary within actor across time slices.",
			"i" = "{.fn lame_dynamic_als} uses period-specific row values only for row coefficients selected by {.arg dynamic_beta}.",
			"i" = "Using per-actor means for static row coefficient{?s}: {.val {vary_row_static}}."
		))
	}
	if (length(vary_col_static) > 0L) {
		cli::cli_warn(c(
			"Some column node covariates vary within actor across time slices.",
			"i" = "{.fn lame_dynamic_als} uses period-specific column values only for column coefficients selected by {.arg dynamic_beta}.",
			"i" = "Using per-actor means for static column coefficient{?s}: {.val {vary_col_static}}."
		))
	}
	if (!isTRUE(dynamic_ab) && !isTRUE(beta_dyn$any) &&
			!isTRUE(dynamic_uv)) {
		return(lame_als(
			Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			R = R_eff,
			family = family,
			mode = mode,
			symmetric = symmetric,
			verbose = verbose,
			seed = seed
		))
	}
	# a static fit informs from its own preparation, a dynamic one here
	.lame_inform_cov_na(prep$cov_na$n)

	# unit of a normal response: directed and undirected position rows carry its
	# square root, so their penalties carry the unit; bipartite rows are unit free
	# (the interaction matrix carries the unit), so theirs carry its square. the
	# convergence tests are relative to the unit
	unit_sd <- .lame_unit_sd(prep$Yarr, family)
	uv_unit <- if (!identical(family, "normal")) 1 else if (prep$bip) unit_sd^2 else unit_sd
	lambda_uv_default <- is.null(.lda_default(lambda_uv, prior$lambda_uv_als))
	lams <- .lda_default_lambdas(
		prep,
		lambda_ab = .lda_default(lambda_ab, prior$lambda_ab_als),
		lambda_beta = .lda_default(lambda_beta, prior$lambda_beta_als),
		lambda_uv = .lda_default(lambda_uv, prior$lambda_uv_als)
	)
	lambda_ab <- as.numeric(lams$ab)
	lambda_beta <- as.numeric(lams$beta)
	lambda_uv <- as.numeric(lams$uv)
	if (lambda_uv_default) lambda_uv <- lambda_uv * uv_unit
	# AR(1) persistence of each smoothness penalty: not estimated, it sits at the
	# prior mean the MCMC uses, so an ALS warm start hands the chain its usual start
	rho_ab <- .lda_clip_rho(prior$rho_ab_mean, 0.8)
	rho_beta <- .lda_clip_rho(prior$rho_beta_mean, 0.8)
	rho_uv <- .lda_clip_rho(prior$rho_uv_mean, 0.9)

	if (verbose) {
		cli::cli_h3("Dynamic ALS / MAP point fit")
		rank_msg <- if (prep$bip) {
			paste0("R_row = ", R_row_eff, ", R_col = ", R_col_eff)
		} else paste0("R = ", R)
		cli::cli_text("Network: {prep$nr} x {prep$nc}, {prep$Tt} time slices, {rank_msg}")
		cli::cli_text("Family: {.field {family}}{ifelse(identical(family, 'normal'), '', paste0(' / ', link, ' IWLS'))}")
		cli::cli_text("Dynamic UV: {.val {isTRUE(dynamic_uv)}} | Dynamic additive effects: {.val {isTRUE(dynamic_ab)}} | Dynamic beta: {.val {isTRUE(beta_dyn$any)}}")
	}

	# seed the paths from a static fit: a warm start handed in by a refit, or a
	# fresh internal pass over this data
	fit0 <- .lda_usable_start(
		start,
		prep,
		R_eff,
		symmetric,
		R_row_eff,
		R_col_eff
	)
	# only a warm start handed in by a refit is a state the fit may prefer over
	# what its reweights find
	warm_start_used <- !is.null(fit0)
	if (!is.null(start) && is.null(fit0)) {
		cli::cli_warn(c(
			"{.arg start} does not match this model's actor set or latent rank.",
			"i" = "Seeding the paths from a fresh static pass instead."
		))
	}
	# the static seed pass is only an initializer: its iteration-limit condition
	# is not shown (it is recorded in fit$diagnostics$static_seed), and the dynamic
	# fit must pass its own score check
	seed_pass <- NULL
	if (is.null(fit0)) {
		seed_stopped <- FALSE
		fit0 <- withCallingHandlers(
			.ae_run(
				Y,
				Xdyad,
				Xrow,
				Xcol,
				R_eff,
				family,
				mode,
				symmetric,
				max_iter = min(100L, max_iter),
				tol = tol,
				verbose = FALSE,
				seed = seed,
				longitudinal = TRUE,
				call = call,
				lowrank_method = lowrank_method,
				non_normal_method = "transform",
				link = "probit",
				linear_solver = "eigen",
				multistart = "none",
				node_time_policy = "quiet",
				certify_gaussian = FALSE,
				report_cov_na = FALSE
			),
			lame_als_convergence = function(w) {
				if (inherits(w, "warning")) {
					seed_stopped <<- TRUE
					invokeRestart("muffleWarning")
				}
			}
		)
		seed_pass <- list(converged = !seed_stopped, iterations = fit0$iterations)
	}
	coef0 <- fit0$coefficients[coef_names]
	coef0[is.na(coef0)] <- 0
	coef0[!is.finite(coef0)] <- 0
	coef_path <- matrix(
		rep(coef0, prep$Tt),
		nrow = length(coef_names),
		dimnames = list(coef_names, prep$time_names)
	)
	a_mat <- matrix(fit0$a, prep$nr, prep$Tt)
	b_mat <- matrix(fit0$b, prep$nc, prep$Tt)
	# a start that carries fitted trajectories seeds the paths with them rather
	# than with their time mean
	fits_path <- function(x, nr) is.matrix(x) && nrow(x) == nr &&
		ncol(x) == prep$Tt && all(is.finite(x))
	if (fits_path(fit0$a_path, prep$nr)) a_mat[] <- fit0$a_path
	if (fits_path(fit0$b_path, prep$nc)) b_mat[] <- fit0$b_path
	if (is.matrix(fit0$coef_path) && ncol(fit0$coef_path) == prep$Tt &&
			all(coef_names %in% rownames(fit0$coef_path)) &&
			all(is.finite(fit0$coef_path[coef_names, ]))) {
		coef_path[] <- fit0$coef_path[coef_names, , drop = FALSE]
	}
	Omat <- matrix(0, prep$nr, prep$nc)
	U <- NULL
	V <- NULL
	L <- NULL
	G_static <- NULL

	if (R_eff > 0L && !is.null(fit0$U) &&
			(isTRUE(symmetric) || prep$bip || !is.null(fit0[["V"]]))) {
		if (prep$bip) {
			row_dim_names <- if (R_row_eff > 0L) {
				paste0("row_dim", seq_len(R_row_eff))
			} else character(0)
			col_dim_names <- if (R_col_eff > 0L) {
				paste0("col_dim", seq_len(R_col_eff))
			} else character(0)
			U <- .lda_coerce_factor(
				fit0$U,
				prep$nr,
				R_row_eff,
				seed = seed + 11L
			)
			V <- .lda_coerce_factor(
				fit0[["V"]],
				prep$nc,
				R_col_eff,
				seed = seed + 17L
			)
			dimnames(U) <- list(prep$row_names, row_dim_names)
			dimnames(V) <- list(prep$col_names, col_dim_names)
			G_static <- .lda_rect_identity(R_row_eff, R_col_eff)
			if (is.matrix(fit0$G) &&
					identical(dim(fit0$G), c(R_row_eff, R_col_eff)) &&
					all(is.finite(fit0$G))) G_static[] <- fit0$G
			if (identical(family, "normal")) {
				# the static seed splits the response unit evenly over the rows of both sides;
				# here the rows are unit free and the interaction matrix carries the unit
				U <- U / sqrt(unit_sd)
				V <- V / sqrt(unit_sd)
				G_static <- G_static * unit_sd
			}

			dimnames(G_static) <- list(row_dim_names, col_dim_names)
			if (isTRUE(dynamic_uv)) {
				U <- array(
					rep(U, prep$Tt),
					dim = c(prep$nr, R_row_eff, prep$Tt),
					dimnames = list(
						prep$row_names,
						row_dim_names,
						prep$time_names
					)
				)
				V <- array(
					rep(V, prep$Tt),
					dim = c(prep$nc, R_col_eff, prep$Tt),
					dimnames = list(
						prep$col_names,
						col_dim_names,
						prep$time_names
					)
				)
				Oarr <- .lda_uvg_to_o_array(
					U,
					V,
					G_static,
					Tt = prep$Tt,
					time_names = prep$time_names
				)
			} else {
				Oarr <- .lda_uvg_to_o_array(
					U,
					V,
					G_static,
					Tt = prep$Tt,
					time_names = prep$time_names
				)
			}
			Omat <- apply(Oarr, c(1, 2), mean)
		} else {
			uv_names <- paste0("dim", seq_len(R))
			if (isTRUE(dynamic_uv) && isTRUE(symmetric) && !is.null(fit0$L)) {
				O0 <- fit0$U %*% (fit0$L * t(fit0$U))
				sf <- .lda_symmetric_factor_slice(O0, R, prep$row_names, uv_names)
				U <- sf$U
				V <- sf[["V"]]
				L <- NULL
				Omat <- sf$O
			} else {
				U <- fit0$U
				V <- if (isTRUE(symmetric) && !is.null(fit0$U)) {
					if (!is.null(fit0$L)) {
						sweep(fit0$U, 2, sign(fit0$L), "*")
					} else fit0$U
				} else fit0[["V"]]
				L <- fit0$L
				Omat <- if (isTRUE(symmetric) && !is.null(fit0$L)) {
					fit0$U %*% (fit0$L * t(fit0$U))
				} else {
					tcrossprod(U, V)
				}
			}
			if (isTRUE(dynamic_uv)) {
				U <- array(
					rep(U, prep$Tt),
					dim = c(prep$nr, R, prep$Tt),
					dimnames = list(prep$row_names, uv_names, prep$time_names)
				)
				V <- array(
					rep(V, prep$Tt),
					dim = c(prep$nc, R, prep$Tt),
					dimnames = list(prep$col_names, uv_names, prep$time_names)
				)
				uv <- .lda_align_uv_cube(U, V)
				U <- uv$U
				V <- uv[["V"]]
				Oarr <- .lda_uv_to_o_array(U, V)
				if (isTRUE(symmetric)) {
					sf <- .lda_symmetric_o_to_uv_cube(
						Oarr,
						R,
						row_names = prep$row_names,
						dim_names = uv_names,
						time_names = prep$time_names
					)
					U <- sf$U
					V <- sf[["V"]]
					Oarr <- sf$Oarr
				}
				Omat <- apply(Oarr, c(1, 2), mean)
				L <- NULL
			} else {
				if (!is.null(U)) dimnames(U) <- list(prep$row_names, uv_names)
				if (!is.null(V)) dimnames(V) <- list(prep$col_names, uv_names)
				Oarr <- .lda_o_array(Omat, prep$Tt)
			}
		}
	} else {
		Oarr <- .lda_o_array(Omat, prep$Tt)
	}
	if (R_eff > 0L && isTRUE(dynamic_uv) && !is.null(U) && !is.null(V) &&
			length(dim(U)) == 3L) {
		up <- fit0$U_path
		vp <- fit0[["V_path"]]
		if (length(dim(up)) == 3L && length(dim(vp)) == 3L &&
				all(dim(up) == dim(U)) && all(dim(vp) == dim(V)) &&
				all(is.finite(up)) && all(is.finite(vp))) {
			U[] <- up
			V[] <- vp
			Oarr <- if (prep$bip) {
				.lda_uvg_to_o_array(
					U,
					V,
					G_static,
					Tt = prep$Tt,
					time_names = prep$time_names
				)
			} else {
				.lda_uv_to_o_array(U, V)
			}
			Omat <- apply(Oarr, c(1, 2), mean)
		}
	}
	if (start_jitter > 0) {
		set.seed(seed)
		coef_path <- coef_path +
			matrix(
				stats::rnorm(length(coef_path), 0, start_jitter / 5),
				nrow = nrow(coef_path),
				dimnames = dimnames(coef_path)
			)
		a_mat <- a_mat + matrix(
			stats::rnorm(length(a_mat), 0, start_jitter),
			nrow = nrow(a_mat),
			dimnames = dimnames(a_mat)
		)
		b_mat <- b_mat + matrix(
			stats::rnorm(length(b_mat), 0, start_jitter),
			nrow = nrow(b_mat),
			dimnames = dimnames(b_mat)
		)
		if (R_eff > 0L && !is.null(U) && !is.null(V)) {
			U <- U + array(
				stats::rnorm(length(U), 0, start_jitter),
				dim = dim(U),
				dimnames = dimnames(U)
			)
			V <- V + array(
				stats::rnorm(length(V), 0, start_jitter),
				dim = dim(V),
				dimnames = dimnames(V)
			)
			if (prep$bip && !is.null(G_static)) {
				Oarr <- .lda_uvg_to_o_array(
					U,
					V,
					G_static,
					Tt = prep$Tt,
					time_names = prep$time_names
				)
				Omat <- apply(Oarr, c(1, 2), mean)
			} else if (length(dim(U)) == 3L) {
				uv <- .lda_align_uv_cube(U, V)
				U <- uv$U
				V <- uv[["V"]]
				Oarr <- .lda_uv_to_o_array(U, V)
				if (isTRUE(symmetric)) {
					sf <- .lda_symmetric_o_to_uv_cube(
						Oarr,
						R,
						row_names = prep$row_names,
						dim_names = dimnames(U)[[2L]],
						time_names = prep$time_names
					)
					U <- sf$U
					V <- sf[["V"]]
					Oarr <- sf$Oarr
				}
				Omat <- apply(Oarr, c(1, 2), mean)
			} else {
				Omat <- tcrossprod(U, V)
				Oarr <- .lda_o_array(Omat, prep$Tt)
			}
		}
	}

	Z <- if (identical(family, "normal")) prep$Yarr else .ae_working_response(prep$Yarr, family)
	W <- .lda_unit_weights(prep$Yarr)
	if (isTRUE(symmetric)) {
		for (t in seq_len(prep$Tt)) {
			Z[, , t] <- .ae_symmetrize(Z[, , t])
			W[, , t] <- .ae_symmetrize(W[, , t])
		}
	}
	use_irls <- !identical(family, "normal")
	irls_trace <- numeric(0)
	transition_u <- NULL
	transition_v <- NULL
	transition_a <- .lda_transition_from_presence(prep$row_presence)
	transition_b <- .lda_transition_from_presence(prep$col_presence)
	fit_mask <- is.finite(Z) & is.finite(W) & W > 0
	objective_trace <- numeric(0)
	convergence_trace <- data.frame(
		iter = integer(0),
		objective = numeric(0),
		relative_objective_delta = numeric(0),
		max_fitted_delta = numeric(0),
		max_parameter_delta = numeric(0)
	)
	converged <- FALSE
	old_eta <- NULL
	# fixed marginal level prior on the multiplicative terms (see .lda_level_prec):
	# active for dynamic factor paths, static bipartite factor blocks and the
	# static unipartite low-rank block of the count and binary families
	# (dynamic_als_lowrank_ridge.R); the normal unipartite static block has a
	# minimiser without it
	uv_level_active <- R_eff > 0L && (isTRUE(dynamic_uv) || prep$bip || isTRUE(use_irls))
	g_level_active <- prep$bip && R_eff > 0L
	level_prec_uv <- if (uv_level_active) .lda_level_prec * uv_unit else 0
	level_prec_g <- if (g_level_active) .lda_level_prec else 0
	# under a saturating likelihood an actor-period with no ties has only the
	# smoothness penalty holding its effect finite, so non-normal families add the
	# level prior to the a/b paths too
	level_prec_ab <- if (isTRUE(use_irls) && isTRUE(dynamic_ab)) .lda_level_prec else 0
	# non-normal families run IRLS as an outer loop: the working problem is rebuilt
	# once the block sweeps converge, and a reweighted iterate is kept only if it
	# lowers the penalized deviance. directed count paths damp an overshoot, and an
	# unresolved overshoot restores the best state
	irls_maxit <- as.integer(max_iter)
	# a short inner run per working problem acts as damped IRLS, so the working
	# response is not over-fitted before the next reweight
	irls_inner_cap <- min(40L, as.integer(max_iter))
	# a normal fit's sweep limit is raised once by .ae_budget_factor
	# (ame_als_refusals.R) before the fit is refused; the IRLS limit is on
	# reweighting steps (up to 40 sweeps each) and is not raised
	max_iter_run <- if (isTRUE(use_irls)) as.integer(max_iter) else as.integer(
		min(.Machine$integer.max, .ae_budget_factor * max_iter)
	)
	sweep_budget <- if (isTRUE(use_irls)) irls_inner_cap * irls_maxit else max_iter_run
	stop_reason <- "iteration limit"
	inner_cap_now <- irls_inner_cap
	outer_k <- 1L
	iter_inner <- 0L
	reweight_steps <- integer(0)
	best_state <- NULL
	best_pdev <- Inf
	irls_reweight <- function() {
		eta_start <- .lda_eta(
			prep,
			coef_path,
			a_mat,
			b_mat,
			Oarr,
			symmetric = symmetric,
			beta_dyn = beta_dyn
		)
		work <- .lda_irls_work(
			prep$Yarr,
			.lda_eta_array(eta_start),
			family = family,
			link = link
		)
		Zn <- work$Z
		Wn <- work$W
		if (isTRUE(symmetric)) {
			for (t in seq_len(prep$Tt)) {
				Zn[, , t] <- .ae_symmetrize(Zn[, , t])
				Wn[, , t] <- .ae_symmetrize(Wn[, , t])
			}
		}
		list(Z = Zn, W = Wn)
	}
	snapshot_state <- function() list(
		coef_path = coef_path,
		a_mat = a_mat,
		b_mat = b_mat,
		U = U,
		V = V,
		L = L,
		G_static = G_static,
		Oarr = Oarr,
		Omat = Omat,
		Z = Z,
		W = W
	)
	stationarity_settings <- function() {
		list(
			family = family,
			link = link,
			symmetric = symmetric,
			beta_dyn = beta_dyn,
			lambda_beta = lambda_beta,
			dynamic_beta_kind = dynamic_beta_kind,
			rho_beta = rho_beta,
			dynamic_ab = dynamic_ab,
			lambda_ab = lambda_ab,
			rho_ab = rho_ab,
			level_prec_ab = level_prec_ab,
			dynamic_uv_kind = dynamic_uv_kind,
			lambda_uv = lambda_uv,
			rho_uv = rho_uv,
			level_prec_uv = level_prec_uv,
			level_prec_g = level_prec_g,
			dynamic_uv = dynamic_uv,
			unit = unit_sd
		)
	}
	stationarity_at <- function() {
		.lda_stationarity(prep, snapshot_state(), stationarity_settings(), tol)
	}
	# penalized objective of the current state and working problem, which every
	# block update minimizes in its own parameters
	objective_now <- function() {
		.lda_objective(
			Z,
			W,
			prep,
			coef_path,
			a_mat,
			b_mat,
			Oarr,
			beta_dyn,
			lambda_beta,
			dynamic_beta_kind,
			rho_beta,
			dynamic_ab,
			lambda_ab,
			rho_ab,
			transition_a = transition_a,
			transition_b = transition_b,
			symmetric = symmetric,
			dynamic_uv = dynamic_uv,
			U = U,
			V = V,
			lambda_uv = lambda_uv,
			rho_uv = rho_uv,
			transition_u = transition_u,
			transition_v = transition_v,
			dynamic_uv_kind = dynamic_uv_kind,
			G_static = G_static,
			level_prec_uv = level_prec_uv,
			level_prec_g = level_prec_g,
			level_prec_ab = level_prec_ab,
			L = L
		)
	}
	restore_state <- function(st) {
		coef_path <<- st$coef_path
		a_mat <<- st$a_mat
		b_mat <<- st$b_mat
		U <<- st$U
		V <<- st[["V"]]
		L <<- st$L
		G_static <<- st$G_static
		Oarr <<- st$Oarr
		Omat <<- st$Omat
		Z <<- st$Z
		W <<- st$W
		invisible(NULL)
	}
	# a block's parameters a fraction f of the way from `from` to `to`, with the
	# derived arrays rebuilt; FALSE when the block is not a vector space (a mix of
	# two rank-R matrices has higher rank)
	blend_block <- function(block, from, to, f) {
		mix <- function(x0, x1) if (is.null(x0) || is.null(x1)) x1 else x0 + f * (x1 - x0)
		if (identical(block, "beta")) {
			coef_path <<- mix(from$coef_path, to$coef_path)
			return(TRUE)
		}
		if (identical(block, "ab")) {
			coef_path <<- mix(from$coef_path, to$coef_path)
			a_mat <<- mix(from$a_mat, to$a_mat)
			b_mat <<- mix(from$b_mat, to$b_mat)
			return(TRUE)
		}
		if (R_eff <= 0L) return(TRUE)
		if (!isTRUE(dynamic_uv) && !prep$bip) return(FALSE)
		U <<- mix(from$U, to$U)
		V <<- mix(from[["V"]], to[["V"]])
		G_static <<- mix(from$G_static, to$G_static)
		Oarr <<- if (prep$bip) {
			.lda_uvg_to_o_array(U, V, G_static, Tt = prep$Tt, time_names = prep$time_names)
		} else {
			.lda_uv_to_o_array(U, V)
		}
		Omat <<- apply(Oarr, c(1L, 2L), mean)
		L <<- NULL
		TRUE
	}
	# descent guard: each block minimizes the objective exactly, so a sweep that
	# raises it (a ridge fallback, non-finite solve or a rotation that broke the
	# penalty) is undone and redone with every block checked, stepping a block back
	# along its displacement until it stops raising the objective. once needed,
	# the guard stays on for the rest of the fit. its tolerance is relative to
	# |objective| + obj_unit (the response variance for a normal fit, one otherwise)
	obj_unit <- if (identical(family, "normal")) {
		local({
			vy <- stats::var(as.numeric(prep$Yarr[is.finite(prep$Yarr)]))
			if (is.finite(vy) && vy > 0) vy else 1
		})
	} else 1
	guard_on <- FALSE
	guard_log <- list(
		redone = 0L,
		halved = 0L,
		rolled_back = 0L,
		blocks = character(0),
		stalled = FALSE,
		iteration = NA_integer_
	)
	sweep_rollbacks <- 0L
	guarded <- function(block, expr) {
		if (!guard_on) return(eval.parent(substitute(expr)))
		start <- snapshot_state()
		before <- objective_now()
		eval.parent(substitute(expr))
		after <- objective_now()
		slack <- 1e-10 * (abs(before) + obj_unit)
		if (is.finite(after) && after <= before + slack) return(invisible(NULL))
		end <- snapshot_state()
		accepted <- FALSE
		for (k in seq_len(12L)) {
			if (!blend_block(block, start, end, 2^(-k))) break
			trial <- objective_now()
			if (is.finite(trial) && trial <= before + slack) {
				accepted <- TRUE
				break
			}
		}
		if (accepted) {
			guard_log$halved <<- guard_log$halved + 1L
		} else {
			restore_state(start)
			guard_log$rolled_back <<- guard_log$rolled_back + 1L
			sweep_rollbacks <<- sweep_rollbacks + 1L
		}
		guard_log$blocks <<- union(guard_log$blocks, block)
		invisible(NULL)
	}
	# the seeded additive paths enter on the per-period sum-to-zero gauge, so every
	# later iterate is feasible and the constrained solves are descent steps
	if (isTRUE(dynamic_ab) && prep$Tt > 0L) {
		for (t in seq_len(prep$Tt)) {
			a_mat[, t] <- a_mat[, t] - mean(a_mat[, t])
			b_mat[, t] <- b_mat[, t] - mean(b_mat[, t])
		}
		if (isTRUE(symmetric)) b_mat <- a_mat
	}
	# penalized deviance of the starting state, on the reweight gate's scale;
	# seeding the gate with it means a restart can never return a worse state
	penalized_deviance_at <- function() {
		etal <- .lda_eta(
			prep,
			coef_path,
			a_mat,
			b_mat,
			Oarr,
			symmetric = symmetric,
			beta_dyn = beta_dyn
		)
		earr <- .lda_eta_array(etal)
		tu <- .lda_transition_from_presence(prep$row_presence)
		tv <- .lda_transition_from_presence(prep$col_presence)
		if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) {
			tu <- .lda_combine_transition_weights(
				.lda_t_transition_weights(U, rho_uv, lambda_uv),
				tu
			)
			tv <- .lda_combine_transition_weights(
				.lda_t_transition_weights(V, rho_uv, lambda_uv),
				tv
			)
		}
		pen0 <- .lda_objective(
			earr,
			array(1, dim(earr)),
			prep,
			coef_path,
			a_mat,
			b_mat,
			Oarr,
			beta_dyn,
			lambda_beta,
			dynamic_beta_kind,
			rho_beta,
			dynamic_ab,
			lambda_ab,
			rho_ab,
			transition_a = transition_a,
			transition_b = transition_b,
			symmetric = symmetric,
			dynamic_uv = dynamic_uv,
			U = U,
			V = V,
			lambda_uv = lambda_uv,
			rho_uv = rho_uv,
			transition_u = tu,
			transition_v = tv,
			dynamic_uv_kind = dynamic_uv_kind,
			G_static = G_static,
			level_prec_uv = level_prec_uv,
			level_prec_g = level_prec_g,
			level_prec_ab = level_prec_ab,
			L = L
		)
		.lda_penalized_deviance(prep, family, link, etal, pen0)
	}
	# warm restarts keep the supplied state as a candidate; directed count paths
	# also keep their cold seed and damp overshooting IRLS proposals
	count_damping <- identical(family, "poisson") && !prep$bip && !symmetric &&
		isTRUE(dynamic_uv) &&
		identical(dynamic_uv_kind, "ar1") && !isTRUE(beta_dyn$any)
	line_search_trace <- list()
	if (isTRUE(use_irls) && (isTRUE(warm_start_used) || count_damping)) {
		rw0 <- irls_reweight()
		Z <- rw0$Z
		W <- rw0$W
		start_pdev <- penalized_deviance_at()
		if (is.finite(start_pdev)) {
			best_pdev <- start_pdev
			best_state <- snapshot_state()
			# the start is a candidate the fit can return, so it belongs in the
			# penalized-deviance trace
			irls_trace <- start_pdev
		}
	}
	for (iter in seq_len(sweep_budget)) {
		iter_inner <- iter_inner + 1L
		if (isTRUE(use_irls) && iter == 1L) {
			rw <- irls_reweight()
			Z <- rw$Z
			W <- rw$W
		}
		transition_u <- .lda_transition_from_presence(prep$row_presence)
		transition_v <- .lda_transition_from_presence(prep$col_presence)
		if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) {
			transition_u <- .lda_combine_transition_weights(
				.lda_t_transition_weights(U, rho_uv, lambda_uv),
				transition_u
			)
			transition_v <- .lda_combine_transition_weights(
				.lda_t_transition_weights(V, rho_uv, lambda_uv),
				transition_v
			)
		}
		coef_old <- coef_path
		a_old <- a_mat
		b_old <- b_mat
		O_old <- Oarr
		G_old <- G_static
		sweep_start <- snapshot_state()
		sweep_rollbacks <- 0L
		# the objective the sweep must not raise: the last sweep's within a working
		# problem, the start's on its first sweep. the fit's first sweep has none, as
		# its start need not be on the gauge of the block solves
		sweep_ref <- if (guard_on || iter == 1L) NA_real_ else if (iter_inner > 1L && length(objective_trace) > 0L) {
			objective_trace[length(objective_trace)]
		} else objective_now()

		guarded("beta", {
			coef_path <- .lda_update_beta_path(
				Z,
				W,
				prep,
				coef_path,
				a_mat,
				b_mat,
				Oarr,
				beta_dyn = beta_dyn,
				lambda_beta = lambda_beta,
				kind = dynamic_beta_kind,
				rho_beta = rho_beta,
				symmetric = symmetric
			)
		})
		Xb <- .lda_xbeta_array(
			prep,
			coef_path,
			symmetric = symmetric,
			beta_dyn = beta_dyn
		)
		guarded("ab", {
			if (isTRUE(dynamic_ab)) {
				ab <- if (isTRUE(symmetric)) {
					.lda_update_dynamic_ab_symmetric(
						Z,
						W,
						Xb,
						Oarr,
						a_mat,
						lambda_ab = lambda_ab,
						rho_ab = rho_ab,
						kind = "ar1",
						transition_a = transition_a,
						level_prec = level_prec_ab
					)
				} else {
					.lda_update_dynamic_ab(
						Z,
						W,
						Xb,
						Oarr,
						a_mat,
						b_mat,
						lambda_ab = lambda_ab,
						rho_ab = rho_ab,
						kind = "ar1",
						transition_a = transition_a,
						transition_b = transition_b,
						level_prec = level_prec_ab
					)
				}
				a_mat <- ab$a
				b_mat <- ab$b
			} else {
				ab <- if (isTRUE(beta_dyn$mask["intercept"]) ||
						!"intercept" %in% rownames(coef_path)) {
					.lda_update_static_ab_centered(Z, W, Xb, Oarr, symmetric)
				} else {
					.lda_update_static_ab(Z, W, Xb, Oarr, symmetric)
				}
				a_mat <- ab$a
				b_mat <- ab$b
				# the additive block centres itself; carrying the removed level into the
				# intercept keeps the fitted values and the objective unchanged
				if (is.finite(ab$level) && ab$level != 0 &&
						"intercept" %in% rownames(coef_path)) {
					coef_path["intercept", ] <- coef_path["intercept", ] + ab$level
					Xb <- Xb + ab$level
				}
			}
		})
		if (guard_on) Xb <- .lda_xbeta_array(
			prep,
			coef_path,
			symmetric = symmetric,
			beta_dyn = beta_dyn
		)
		target_mult <- array(0, dim = dim(Z))
		for (t in seq_len(prep$Tt)) {
			target_mult[, , t] <- Z[, , t] - Xb[, , t] -
				outer(a_mat[, t], b_mat[, t], "+")
		}
		guarded("mult", {
			if (isTRUE(dynamic_uv)) {
				if (prep$bip) {
					uv <- .lda_update_dynamic_uv_bip_g(
						Z,
						W,
						Xb,
						a_mat,
						b_mat,
						U,
						V,
						G_static,
						lambda_uv = lambda_uv,
						rho_uv = rho_uv,
						transition_u = transition_u,
						transition_v = transition_v,
						level_prec = level_prec_uv,
						ridge = 1e-6 * uv_unit
					)
					U <- uv$U
					V <- uv[["V"]]
					G_static <- .lda_update_static_g_bip(
						target_mult,
						is.finite(Z),
						W,
						U,
						V,
						G_static,
						level_prec = level_prec_g
					)
					Oarr <- .lda_uvg_to_o_array(
						U,
						V,
						G_static,
						Tt = prep$Tt,
						time_names = prep$time_names
					)
				} else if (isTRUE(symmetric)) {
					uv <- .lda_update_dynamic_uv_symmetric(
						Z,
						W,
						Xb,
						a_mat,
						b_mat,
						U,
						V,
						lambda_uv = lambda_uv,
						rho_uv = rho_uv,
						transition_u = transition_u,
						transition_v = transition_v,
						level_prec = level_prec_uv,
						ridge = 1e-6 * uv_unit
					)
					U <- uv$U
					V <- uv[["V"]]
					Oarr <- uv$Oarr
				} else {
					uv <- .lda_update_dynamic_uv_directed(
						Z,
						W,
						Xb,
						a_mat,
						b_mat,
						U,
						V,
						lambda_uv = lambda_uv,
						rho_uv = rho_uv,
						transition_u = transition_u,
						transition_v = transition_v,
						level_prec = level_prec_uv,
						ridge = 1e-6 * uv_unit
					)
					U <- uv$U
					V <- uv[["V"]]
					Oarr <- uv$Oarr
				}
				Omat <- apply(Oarr, c(1, 2), mean)
				L <- NULL
			} else if (R_eff > 0L) {
				if (prep$bip) {
					lr <- .lda_update_static_uv_bip_g(
						target_mult,
						is.finite(Z),
						W,
						U,
						V,
						G_static,
						level_prec = level_prec_uv,
						ridge = 1e-6 * uv_unit
					)
					U <- lr$U
					V <- lr[["V"]]
					G_static <- .lda_update_static_g_bip(
						target_mult,
						is.finite(Z),
						W,
						U,
						V,
						G_static,
						level_prec = level_prec_g
					)
					Oarr <- .lda_uvg_to_o_array(
						U,
						V,
						G_static,
						Tt = prep$Tt,
						time_names = prep$time_names
					)
					Omat <- apply(Oarr, c(1, 2), mean)
					L <- NULL
				} else {
					lr <- .lda_update_static_lowrank(
						Z,
						W,
						Xb,
						a_mat,
						b_mat,
						Omat,
						R,
						symmetric = symmetric,
						lowrank_method = lowrank_method,
						ridge = level_prec_uv,
						rel = identical(family, "normal")
					)
					Omat <- lr$Omat
					U <- lr$U
					V <- lr[["V"]]
					L <- lr$L
					Oarr <- .lda_o_array(Omat, prep$Tt)
				}
			}

			if (prep$bip && R_eff > 0L) {
				penalty_side <- function(cube, transition) {
					value <- level_prec_uv * sum(cube^2)
					if (isTRUE(dynamic_uv) && prep$Tt > 1L) {
						for (i in seq_len(dim(cube)[1L])) {
							P <- .lda_path_penalty_weighted(
								prep$Tt,
								lambda_uv,
								if (is.null(transition)) NULL else transition[i, ],
								rho_uv
							)
							for (r in seq_len(dim(cube)[2L])) value <- value +
								drop(crossprod(cube[i, r, ], P %*% cube[i, r, ]))
						}
					}
					value
				}
				pu <- penalty_side(U, transition_u)
				pv <- penalty_side(V, transition_v)
				pg <- level_prec_g * sum(G_static^2)
				if (all(is.finite(c(pu, pv, pg))) && min(pu, pv, pg) > 0) {
					common <- exp(mean(log(c(pu, pv, pg))))
					su <- sqrt(common / pu); sv <- sqrt(common / pv)
					U <- U * su; V <- V * sv; G_static <- G_static / (su * sv)
					Oarr <- .lda_uvg_to_o_array(
						U,
						V,
						G_static,
						Tt = prep$Tt,
						time_names = prep$time_names
					)
					Omat <- apply(Oarr, c(1L, 2L), mean)
				}
			}
		})

		obj <- objective_now()
		if (!guard_on && (!is.finite(obj) ||
				(is.finite(sweep_ref) && obj > sweep_ref + 1e-10 * (abs(sweep_ref) + obj_unit)))) {
			# a block raised the objective: undo the sweep and redo it with every block checked
			restore_state(sweep_start)
			guard_on <- TRUE
			guard_log$redone <- guard_log$redone + 1L
			if (is.na(guard_log$iteration)) guard_log$iteration <- as.integer(iter)
			iter_inner <- iter_inner - 1L
			next
		}
		objective_trace <- c(objective_trace, obj)
		eta_now <- .lda_eta(
			prep,
			coef_path,
			a_mat,
			b_mat,
			Oarr,
			symmetric = symmetric,
			beta_dyn = beta_dyn
		)
		max_eta_delta <- if (is.null(old_eta)) Inf else {
			max(vapply(seq_len(prep$Tt), function(t) {
				ok <- fit_mask[, , t]
				if (!any(ok)) return(0)
				max(abs(eta_now[[t]][ok] - old_eta[[t]][ok]), na.rm = TRUE)
			}, numeric(1)))
		}
		old_eta <- eta_now
		rel_obj <- NA_real_
		param_delta <- NA_real_
		if (length(objective_trace) > 1L) {
			prev_obj <- objective_trace[length(objective_trace) - 1L]
			rel_obj <- abs(prev_obj - obj) / max(obj_unit, abs(prev_obj))
			param_delta <- max(
				.lda_max_abs(coef_path - coef_old),
				.lda_max_abs(a_mat - a_old),
				.lda_max_abs(b_mat - b_old),
				.lda_max_abs(Oarr - O_old),
				if (!is.null(G_static) && !is.null(G_old))
					.lda_max_abs(G_static - G_old) else 0
			)
		}
		convergence_trace <- rbind(
			convergence_trace,
			data.frame(
				iter = iter,
				objective = obj,
				relative_objective_delta = rel_obj,
				max_fitted_delta = max_eta_delta,
				max_parameter_delta = param_delta
			)
		)
		# the relative objective change is only meaningful within one working problem
		if (iter_inner <= 1L) rel_obj <- NA_real_
		inner_conv <- iter_inner > 1L && is.finite(rel_obj) && rel_obj < tol &&
			is.finite(max_eta_delta) && max_eta_delta < sqrt(tol) * sqrt(obj_unit)
		# a guarded sweep in which no block could lower the objective has stalled
		stalled <- guard_on && sweep_rollbacks > 0L && is.finite(rel_obj) && rel_obj < tol
		if (!isTRUE(use_irls)) {
			if (inner_conv && isTRUE(stationarity_at()$pass)) {
				converged <- TRUE
				break
			}
			if (stalled) {
				converged <- isTRUE(stationarity_at()$pass)
				guard_log$stalled <- !converged
				break
			}
			next
		}
		if (inner_conv || iter_inner >= inner_cap_now || stalled) {
			eta_now_list <- .lda_eta(
				prep,
				coef_path,
				a_mat,
				b_mat,
				Oarr,
				symmetric = symmetric,
				beta_dyn = beta_dyn
			)
			eta_arr <- .lda_eta_array(eta_now_list)
			pen <- .lda_objective(
				eta_arr,
				array(1, dim(eta_arr)),
				prep,
				coef_path,
				a_mat,
				b_mat,
				Oarr,
				beta_dyn,
				lambda_beta,
				dynamic_beta_kind,
				rho_beta,
				dynamic_ab,
				lambda_ab,
				rho_ab,
				transition_a = transition_a,
				transition_b = transition_b,
				symmetric = symmetric,
				dynamic_uv = dynamic_uv,
				U = U,
				V = V,
				lambda_uv = lambda_uv,
				rho_uv = rho_uv,
				transition_u = transition_u,
				transition_v = transition_v,
				dynamic_uv_kind = dynamic_uv_kind,
				G_static = G_static,
				level_prec_uv = level_prec_uv,
				level_prec_g = level_prec_g,
				level_prec_ab = level_prec_ab,
				L = L
			)
			pdev <- .lda_penalized_deviance(prep, family, link, eta_now_list, pen)
			if (count_damping && !is.null(best_state) &&
					(!is.finite(pdev) || pdev > best_pdev)) {
				proposal <- snapshot_state()
				proposal_value <- pdev
				accepted_fraction <- 0
				for (half in seq_len(30L)) {
					fraction <- 2^(-half)
					coef_path <- best_state$coef_path + fraction *
						(proposal$coef_path - best_state$coef_path)
					a_mat <- best_state$a_mat + fraction * (proposal$a_mat - best_state$a_mat)
					b_mat <- best_state$b_mat + fraction * (proposal$b_mat - best_state$b_mat)
					U <- best_state$U + fraction * (proposal$U - best_state$U)
					V <- best_state$V + fraction * (proposal$V - best_state$V)
					Oarr <- .lda_uv_to_o_array(U, V)
					Omat <- apply(Oarr, c(1L, 2L), mean)
					trial_value <- penalized_deviance_at()
					if (is.finite(trial_value) && trial_value <= best_pdev) {
						pdev <- trial_value
						accepted_fraction <- fraction
						break
					}
				}
				if (accepted_fraction == 0) pdev <- Inf
				line_search_trace[[length(line_search_trace) + 1L]] <- list(
					objective_before = best_pdev,
					proposal = proposal_value,
					objective_after = pdev,
					fraction = accepted_fraction,
					trials = half
				)
			}
			irls_trace <- c(irls_trace, pdev)
			if (!is.finite(pdev) || pdev > best_pdev) {
				# overshoot: restore the best accepted state and stop
				st <- best_state
				if (!is.null(st)) {
					coef_path <- st$coef_path; a_mat <- st$a_mat; b_mat <- st$b_mat
					U <- st$U; V <- st[["V"]]; L <- st$L
					G_static <- st$G_static
					Oarr <- st$Oarr; Omat <- st$Omat; Z <- st$Z; W <- st$W
				}
				# an overshoot does not certify the restored state, so check its score; a
				# non-stationary state is continued with a shorter run (half the sweeps) down
				# to a single sweep
				converged <- isTRUE(stationarity_at()$pass)
				if (!converged && inner_cap_now > 1L && outer_k < irls_maxit && !is.null(best_state)) {
					inner_cap_now <- max(1L, inner_cap_now %/% 2L)
					outer_k <- outer_k + 1L
					rw <- irls_reweight()
					Z <- rw$Z
					W <- rw$W
					reweight_steps <- c(reweight_steps, length(objective_trace))
					iter_inner <- 0L
					old_eta <- NULL
					next
				}
				stop_reason <- "overshoot"
				break
			}
			best_state <- snapshot_state()
			best_pdev <- pdev
			if (isTRUE(stationarity_at()$pass)) {
				converged <- TRUE
				break
			}
			# a damped step is lengthened again after each step that lowered the deviance
			if (inner_cap_now < irls_inner_cap) inner_cap_now <- min(irls_inner_cap, 2L * inner_cap_now)
			if (outer_k >= irls_maxit) {
				converged <- FALSE
				break
			}
			prev_pdev <- pdev
			outer_k <- outer_k + 1L
			rw <- irls_reweight()
			Z <- rw$Z
			W <- rw$W
			reweight_steps <- c(reweight_steps, length(objective_trace))
			iter_inner <- 0L
			old_eta <- NULL
		}
	}

	refinement <- NULL
	if (!converged && !use_irls && prep$bip && isTRUE(dynamic_uv) &&
			!isTRUE(beta_dyn$any) &&
			identical(dynamic_uv_kind, "ar1") && length(U) > 0L && length(V) > 0L) {
		refined <- .lda_refine_bip_normal(
			prep,
			snapshot_state(),
			stationarity_settings(),
			tol,
			as.integer(max_iter)
		)
		refinement <- refined$diagnostics
		if (isTRUE(refined$accepted)) {
			old_state <- snapshot_state()
			st <- refined$state
			coef_path <- st$coef_path
			a_mat <- st$a_mat
			b_mat <- st$b_mat
			U <- st$U
			V <- st$V
			G_static <- st$G_static
			Oarr <- st$Oarr
			Omat <- st$Omat
			converged <- TRUE
			objective_trace <- c(objective_trace, refinement$objective_after)
			old_eta <- .lda_eta_array(.lda_eta(
				prep,
				old_state$coef_path,
				old_state$a_mat,
				old_state$b_mat,
				old_state$Oarr,
				symmetric = symmetric,
				beta_dyn = beta_dyn
			))
			new_eta <- .lda_eta_array(.lda_eta(
				prep,
				coef_path,
				a_mat,
				b_mat,
				Oarr,
				symmetric = symmetric,
				beta_dyn = beta_dyn
			))
			convergence_trace <- rbind(
				convergence_trace,
				data.frame(
					iter = iter + 1L,
					objective = refinement$objective_after,
					relative_objective_delta = abs(
						refinement$objective_after -
							refinement$objective_before
					) / max(1, abs(refinement$objective_before)),
					max_fitted_delta = max(abs(new_eta[fit_mask] - old_eta[fit_mask])),
					max_parameter_delta = max(
						.lda_max_abs(coef_path - old_state$coef_path),
						.lda_max_abs(a_mat - old_state$a_mat),
						.lda_max_abs(b_mat - old_state$b_mat),
						.lda_max_abs(Oarr - old_state$Oarr),
						.lda_max_abs(G_static - old_state$G_static)
					)
				)
			)
		}
	}

	if (!converged && use_irls && identical(family, "poisson") &&
			!prep$bip && !symmetric && isTRUE(dynamic_uv) &&
			!isTRUE(beta_dyn$any) &&
			identical(dynamic_uv_kind, "ar1") && length(U) > 0L && length(V) > 0L) {
		refined <- .lda_refine_directed_count(
			prep,
			snapshot_state(),
			stationarity_settings(),
			tol,
			as.integer(max_iter)
		)
		refinement <- refined$diagnostics
		if (isTRUE(refined$accepted)) {
			st <- refined$state
			coef_path <- st$coef_path
			a_mat <- st$a_mat
			b_mat <- st$b_mat
			U <- st$U
			V <- st$V
			Oarr <- st$Oarr
			Omat <- st$Omat
			converged <- TRUE
			irls_trace <- c(irls_trace, refinement$objective_after)
		}
	}
	# a stall that a certified joint refinement resolved is a converged fit
	if (converged) guard_log$stalled <- FALSE

	if (verbose) {
		if (converged) {
			cli::cli_alert_success("Dynamic ALS converged in {iter} iteration{?s}.")
		} else if (isTRUE(guard_log$stalled)) {
			cli::cli_alert_danger("Dynamic ALS stopped in {iter} iteration{?s}: a block update could not lower the objective.")
		} else {
			cli::cli_alert_danger("Dynamic ALS stopped ({stop_reason}) without reaching its convergence criterion.")
		}
	}
	# a fit that did not converge is refused below, with the fit in the condition
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") &&
			!is.null(U) && length(dim(U)) == 3L) {
		transition_u <- .lda_combine_transition_weights(
			.lda_t_transition_weights(U, rho_uv, lambda_uv),
			.lda_transition_from_presence(prep$row_presence)
		)
		transition_v <- if (!is.null(V) && length(dim(V)) == 3L) {
			.lda_combine_transition_weights(
				.lda_t_transition_weights(V, rho_uv, lambda_uv),
				.lda_transition_from_presence(prep$col_presence)
			)
		} else NULL
	}

	stationarity <- stationarity_at()
	if (isTRUE(use_irls)) { rw <- irls_reweight(); Z <- rw$Z; W <- rw$W }

	fit_out <- .lda_as_fit(
		prep,
		family,
		mode,
		symmetric,
		R_eff,
		R_row_eff,
		R_col_eff,
		call,
		coef_path,
		beta_dyn,
		dynamic_beta_kind,
		a_mat,
		b_mat,
		Omat,
		Oarr,
		U,
		V,
		L,
		G_static,
		objective_trace,
		converged,
		iter,
		lambda_ab,
		lambda_beta,
		lambda_uv,
		rho_ab,
		rho_beta,
		rho_uv,
		dynamic_ab,
		dynamic_uv,
		dynamic_uv_kind,
		lowrank_method,
		link = if (identical(family, "normal")) NULL else link,
		Zwork = Z,
		obs_weights = W,
		irls_trace = irls_trace,
		transition_u = transition_u,
		transition_v = transition_v,
		convergence_trace = convergence_trace,
		tol = tol,
		reweight_steps = reweight_steps
	)
	fit_out$diagnostics$stationarity <- stationarity[c("pass", "max_score", "threshold", "by_block")]
	fit_out$diagnostics$joint_refinement <- refinement
	fit_out$diagnostics$irls_line_search <- line_search_trace
	fit_out$diagnostics$static_seed <- seed_pass
	# descent guard record (active = FALSE when no block solve raised the objective)
	fit_out$diagnostics$descent_guard <- list(
		active = guard_on,
		first_iteration = guard_log$iteration,
		sweeps_redone = guard_log$redone,
		blocks_halved = guard_log$halved,
		blocks_rolled_back = guard_log$rolled_back,
		blocks = guard_log$blocks,
		stalled = guard_log$stalled
	)
	if (isTRUE(guard_log$stalled)) {
		fit_out$diagnostics$message <- "Stopped: a block update could not lower the penalized objective and the others had nothing left to lower."
		cli::cli_abort(
			c(
				"{.fn lame_dynamic_als} stopped without converging: a block update could not lower the penalized objective.",
				"x" = "The {.field {guard_log$blocks}} block update raised the objective at every step the descent guard tried (it first intervened at iteration {guard_log$iteration}), and the other blocks had nothing left to lower.",
				"i" = "The last state that lowered the objective is in the {.field fit} field of the condition; its score is in {.code fit$diagnostics$stationarity}."
			),
			class = "lame_als_convergence",
			fit = fit_out
		)
	}
	if (!converged) {
		st_fin <- fit_out$diagnostics$stationarity
		blk <- if (length(st_fin$by_block)) names(st_fin$by_block)[which.max(st_fin$by_block)] else NA_character_
		limit <- identical(stop_reason, "iteration limit")
		fit_out$diagnostics$message <- if (limit) {
			"Stopped at the iteration limit before the convergence criteria were reached."
		} else {
			"Stopped: the reweighted solve could not lower the penalized deviance."
		}
		cli::cli_abort(
			c(
				"{.fn lame_dynamic_als} did not reach its convergence criterion within {.arg max_iter} = {max_iter}{if (limit && max_iter_run > max_iter) paste0(' (raised once to ', max_iter_run, ')') else ''}.",
				"x" = "It stopped on the {stop_reason} after {fit_out$iterations} sweep{?s}: the largest scaled first-order score is {signif(st_fin$max_score, 3)} against a threshold of {signif(st_fin$threshold, 3)}, in the block {.field {blk}}.",
				"i" = if (limit) {
					"Increase {.arg max_iter}; for {.code R > 0} the objective is non-convex and a flat one converges slowly, so a lower {.arg R} may also finish it."
				} else {
					"The reweighted solve could not lower the penalized deviance, so a larger {.arg max_iter} does not help: lower {.arg R}, or use the default empirical-Bayes estimator ({.code bias_correct = NULL})."
				},
				"i" = "The last iterate is in the {.field fit} field of the condition."
			),
			class = "lame_als_convergence",
			kind = "bcd",
			state = st_fin,
			max_iter = max_iter,
			retry = FALSE,
			fit = fit_out
		)
	}
	count_correction <- identical(family, "poisson") && !prep$bip &&
		!symmetric && !isTRUE(beta_dyn$any) && prep$nr >= 4L &&
		(R_eff == 0L || (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "ar1")))
	if (identical(bias_correct, "jackknife") && !count_correction) {
		cli::cli_abort(
			c(
				"The split-sample jackknife is not available for this dynamic ALS specification.",
				"i" = "It requires a directed Poisson panel with at least four actors, static coefficients, and rank zero or AR(1) positions."
			),
			class = "lame_jk_split"
		)
	}
	if (count_correction && !identical(bias_correct, "none") && converged) {
		fit_out <- .lda_count_jackknife(fit_out, call_args, max_iter, tol)
		fit_out$bias_correct$requested <- bias_correct
	} else {
		fit_out$bias_correct <- list(
			method = "none",
			applied = FALSE,
			requested = bias_correct,
			reason = if (!converged) "the full fit has not converged" else
			if (identical(bias_correct, "none")) "explicitly disabled" else
				"no automatic correction for this specification"
		)
	}
	# everything a refit needs except the data, which the fit already carries
	fit_out$als_refit <- list(
		fitter = "lame_dynamic_als",
		args = call_args[setdiff(
			names(call_args),
			c("Y", "Xdyad", "call")
		)],
		seed = seed
	)
	if (!identical(stability, "none")) {
		fit_out$stability <- .lda_dynamic_stability(
			fit_out,
			starts = stability,
			seed = seed,
			call_args = call_args
		)
	}
	if (bootstrap > 0L) {
		fit_out$bootstrap_dynamic <- .lda_dynamic_bootstrap(
			fit_out,
			R_boot = bootstrap,
			type = bootstrap_type,
			block_length = bootstrap_block_length,
			seed = if (is.null(bootstrap_seed)) seed else bootstrap_seed,
			call_args = call_args,
			verbose = verbose
		)
		fit_out$meta$uncertainty_available <- TRUE
		fit_out$meta$uncertainty_source <- "dynamic_bootstrap"
	}
	fit_out
}

#' @export
coef.lame_dynamic_als <- function(object, ...) {
	if (isTRUE(object$dynamic_beta)) object$coef_path else object$coefficients
}

#' @export
predict.lame_dynamic_als <- function(
	object,
	newdata = NULL,
	type = c("response", "link"),
	...
) {
	# a withheld product (.ae_withhold_product) is restored to its variational means, which this reads
	object <- .ae_product_restore(object)
	type <- match.arg(type)
	# a future horizon passed through ... must not return training-period fitted
	# values: this method changes covariates on the fitted grid only and does not
	# propagate dynamic paths beyond it
	dots <- list(...)
	if ("h" %in% names(dots)) {
		h <- dots[["h"]]
		if (!is.numeric(h) || length(h) != 1L || !is.finite(h) ||
				h < 0 || h != floor(h))
			cli::cli_abort("{.arg h} must be a single non-negative integer.")
		if (h > 0) cli::cli_abort(c(
			"{.fn predict.lame_dynamic_als} does not support h-step forecasting.",
			"i" = "It predicts on the fitted periods; future ALS path propagation is not implemented.",
			"i" = "Use {.fn lame} with {.code method = \"mcmc\"} and supported dynamic components for forecasts."
		))
	}
	mcmc_args <- intersect(c("interval", "by_draw", "newexposure"), names(dots))
	if (length(mcmc_args)) cli::cli_warn(c(
		"Argument{?s} {.arg {mcmc_args}} only apply to MCMC fits and are ignored here.",
		"i" = "For uncertainty on an ALS fit, use {.fn ame_als_bootstrap} or the {.arg bootstrap} fitting argument."
	))
	EZ <- object$EZ
	if (!is.null(newdata)) {
		p_dyad <- dim(object[["X"]])[3L]
		if (p_dyad == 0L) {
			cli::cli_abort("This fit has no dyadic covariates; {.arg newdata} cannot be applied.")
		}
		Xn <- newdata
		if (is.matrix(Xn)) Xn <- array(Xn, c(dim(Xn), 1L, 1L))
		if (length(dim(Xn)) == 3L) Xn <- array(Xn, c(dim(Xn), 1L))
		if (!identical(dim(Xn), dim(object[["X"]]))) {
			cli::cli_abort(c(
				"{.arg newdata} must match the fitted dyadic design dimensions.",
				"i" = "Expected [{paste(dim(object[['X']]), collapse = ', ')}], got [{paste(dim(Xn), collapse = ', ')}]."
			))
		}
		EZ <- vector("list", object$n_time)
		coef_names <- rownames(object$coef_path)
		for (t in seq_len(object$n_time)) {
			xb_new <- .lda_xbeta_cells(
				object,
				Xn,
				t,
				stats::setNames(object$coef_path[, t], coef_names)
			)
			Ot <- if (!is.null(object$Oarr)) object$Oarr[, , t] else object$Omat
			ez <- xb_new +
				outer(object$a_dynamic[, t], object$b_dynamic[, t], "+") +
				Ot
			if (!identical(object$mode, "bipartite") &&
					object$dims$n_row == object$dims$n_col) {
				diag(ez) <- NA_real_
			}
			dimnames(ez) <- list(object$row_names, object$col_names)
			EZ[[t]] <- ez
		}
		names(EZ) <- object$time_names
	}
	# the response estimand is that of fitted() (R/dynamic_als_predictive.R): the
	# predictive probability of an EB binary fit at the requested design, the
	# plug-in value otherwise
	pred <- if (type == "link") EZ else if (.lde_is_binary_eb(object)) {
		.lde_predictive_response(object, EZ, if (is.null(newdata)) object[["X"]] else Xn)
	} else {
		lapply(
			EZ,
			function(ez) .ae_link_inverse(
				ez,
				object$family,
				object$link,
				latent_var = .ae_default(object$pois_latent_var, 0)
			)
		)
	}
	pred
}

#' @method print lame_dynamic_als
#' @export
print.lame_dynamic_als <- function(x, digits = 4, ...) {
	.lame_cli_stdout({
		if (identical(x$bias_correct$method, "eb") && isTRUE(x$bias_correct$dynamic)) {
			cli::cli_text("{.strong Dynamic AME fit via empirical-Bayes ALS} (fast, MCMC-free)")
		} else {
			cli::cli_text("{.strong Dynamic AME fit via penalized ALS/MAP} (fast, MCMC-free)")
		}
		cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}} | R = {x$R}")
		cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slices")
		cli::cli_text("Dynamic UV: {.val {isTRUE(x$dynamic_uv)}} | Dynamic additive effects: {.val {isTRUE(x$dynamic_ab)}} | Dynamic beta: {.val {isTRUE(x$dynamic_beta)}}")
		cli::cli_text("Converged: {.val {x$converged}} in {x$iterations} iteration{?s}")
		objective_label <- x$objective_label %||% "objective"
		if (!is.null(x$objective) && is.finite(x$objective)) {
			cli::cli_text("Reported {objective_label}: {round(x$objective, digits)} (lower is better)")
		}
		cli::cli_text("")
		cli::cli_text("{.strong Coefficients:}")
		if (isTRUE(x$dynamic_beta)) {
			cli::cli_verbatim(utils::capture.output(print(round(x$coef_path, digits))))
		} else {
			cli::cli_verbatim(utils::capture.output(print(round(x$coefficients, digits))))
		}
		cli::cli_text("")
		rho_note <- .ae_rho_not_estimated_text(x)
		if (!is.null(rho_note)) cli::cli_text("{rho_note}")
		.ae_print_bias_correct(x$bias_correct, x$family)
		.lde_print_eb(x, digits)
		cli::cli_text(cli::col_grey("Point estimate only; no posterior draws."))
	})
	invisible(x)
}

#' @method summary lame_dynamic_als
#' @export
summary.lame_dynamic_als <- function(object, ...) {
	out <- list(
		call = object$call,
		family = object$family,
		mode = object$mode,
		R = object$R,
		n_time = object$n_time,
		dims = object$dims,
		coefficients = if (isTRUE(object$dynamic_beta)) object$coef_path else
			cbind(Estimate = object$coefficients),
		VC = object$VC,
		rho_note = .ae_rho_not_estimated_text(object),
		bias_correct = object$bias_correct,
		eb = object$eb,
		dynamic_uv = object$dynamic_uv,
		dynamic_ab = object$dynamic_ab,
		dynamic_beta = object$dynamic_beta,
		converged = object$converged,
		iterations = object$iterations,
		objective = .ae_default(object$objective, object$deviance),
		objective_label = object$objective_label,
		diagnostics = object$diagnostics
	)
	class(out) <- "summary.lame_dynamic_als"
	out
}

#' @method print summary.lame_dynamic_als
#' @export
print.summary.lame_dynamic_als <- function(x, digits = 4, ...) {
	.lame_cli_stdout({
		cli::cli_h2("Dynamic ALS summary")
		cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}} | R = {x$R}")
		cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slices")
		cli::cli_text("Dynamic UV: {.val {isTRUE(x$dynamic_uv)}} | Dynamic additive effects: {.val {isTRUE(x$dynamic_ab)}} | Dynamic beta: {.val {isTRUE(x$dynamic_beta)}}")
		objective_label <- x$objective_label %||%
			(x$diagnostics$reported_objective_label %||% "objective")
		cli::cli_text("Converged: {.val {x$converged}} ({x$iterations} iteration{?s}); {objective_label} = {round(x$objective, digits)}")
		monitor_label <- x$diagnostics$objective_label %||% "objective"
		if (!identical(monitor_label, objective_label)) {
			cli::cli_text(
				cli::col_grey("Block sweeps are monitored on the {monitor_label} ({x$diagnostics$objective_scope %||% 'whole fit'}); {.code fit$dev_history} holds the {objective_label}.")
			)
		}
		if (isTRUE((x$diagnostics$objective_increases %||% 0L) > 0L)) {
			cli::cli_text(
				cli::col_yellow("{monitor_label} increased in {x$diagnostics$objective_increases} step{?s}; max increase = {round(x$diagnostics$max_objective_increase, digits)}.")
			)
		}
		dg <- x$diagnostics$descent_guard
		if (isTRUE(dg$active)) {
			cli::cli_text(
				cli::col_grey("Descent guard from iteration {dg$first_iteration}: {dg$sweeps_redone} sweep{?s} redone, {dg$blocks_halved} block update{?s} stepped back, {dg$blocks_rolled_back} kept at its start.")
			)
		}
		cli::cli_rule()
		cli::cli_text("{.strong Coefficients}")
		cli::cli_verbatim(utils::capture.output(print(round(x$coefficients, digits))))
		cli::cli_text("")
		cli::cli_text("{.strong Variance summaries}")
		cli::cli_verbatim(utils::capture.output(print(round(x$VC, digits))))
		if (!is.null(x$rho_note)) cli::cli_text("{x$rho_note}")
		.ae_print_bias_correct(x$bias_correct, x$family)
		.lde_print_eb(x, digits)
		cli::cli_text(cli::col_grey("Point estimate only; no posterior draws or credible intervals."))
	})
	invisible(x)
}
