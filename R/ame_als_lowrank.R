# multiplicative (low-rank) block solvers for the fast ame estimator. the
# unpenalized solvers (.ae_lowrank_mm, .ae_lowrank_als) minimise the weighted
# objective sum_ij c_ij (rbar_ij - o_ij)^2 with rank(o) <= r and serve the
# normal-family and fixed-transform paths. the IRLS families (binary, poisson)
# use the penalized MAP solvers (.ae_lowrank_em_als, .ae_lowrank_em_mm), whose
# ridge comes from the MCMC prior [u_i, v_i] ~ N(0, Suv), Suv ~ IW(suv0 I, kappa0);
# without it a binary rank-r factor block diverges by quasi-separation

# weighted low-rank objective over observed dyads
.ae_lowrank_obj <- function(O, Rbar, C, obs_dyad) {
	sum(C[obs_dyad] * (Rbar[obs_dyad] - O[obs_dyad])^2)
}

# weighted majorize-minimize update. a cell's step is its weight over the
# largest weight, floored at `cfloor`: binary and count weights are at most one,
# while a normal cell weighs 1 / ve, so the normal fit passes cfloor = 0 and the
# largest weight takes a full step whatever the unit of the response
.ae_lowrank_mm <- function(
	Omat,
	Rbar,
	C,
	obs_dyad,
	R,
	symmetric,
	maxit = 25L,
	cfloor = 1,
	rel = FALSE
) {
	cmax <- max(C); if (cmax < cfloor) cmax <- cfloor
	cmax <- max(cmax, .Machine$double.xmin)
	w_obs <- C[obs_dyad] / cmax
	U <- V <- L <- NULL
	for (mm in seq_len(maxit)) {
		A <- Omat
		A[obs_dyad] <- Omat[obs_dyad] +
			w_obs * (Rbar[obs_dyad] - Omat[obs_dyad])
		if (symmetric) {
			A <- (A + t(A)) / 2
			eg <- eigen(A, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
			U <- eg$vectors[, idx, drop = FALSE]
			L <- eg$values[idx]
			V <- U
			Omat_new <- U %*% (L * t(U))
		} else {
			sv <- svd(A)
			dd <- sqrt(pmax(sv$d[seq_len(R)], 0))
			U <- sweep(sv$u[, seq_len(R), drop = FALSE], 2, dd, "*")
			V <- sweep(sv$v[, seq_len(R), drop = FALSE], 2, dd, "*")
			L <- NULL
			Omat_new <- tcrossprod(U, V)
		}
		delta <- max(abs(Omat_new - Omat))
		Omat <- Omat_new
		# a normal response (`rel`) has an arbitrary unit, so its change is measured
		# against the size of the product; other families use the absolute tolerance
		if (if (rel) delta <= 1e-8 * max(abs(Omat_new)) else delta < 1e-8) break
	}
	list(U = U, V = V, L = L, Omat = Omat)
}

# rebalance U and V to equal scale without changing U V'; ALS can otherwise
# drift along the U G, V G^-T gauge and worsen the row systems' conditioning
.ae_rebalance_uv <- function(U, V) {
	R <- ncol(U)
	if (R == 0L) return(list(U = U, V = V))
	qu <- qr(U); qv <- qr(V)
	sv <- svd(qr.R(qu) %*% t(qr.R(qv)))
	dd <- sqrt(pmax(sv$d, 0))
	list(
		U = qr.Q(qu) %*% sweep(sv$u, 2, dd, "*"),
		V = qr.Q(qv) %*% sweep(sv$v, 2, dd, "*")
	)
}

# weighted alternating least squares (directed / bipartite): with V held fixed each
# row is the exact weighted LS minimiser
#   u_i = (sum_j c_ij v_j v_j')^+ (sum_j c_ij rbar_ij v_j)
# and likewise for V, so the objective is non-increasing across sweeps
.ae_lowrank_als <- function(
	Rbar,
	C,
	obs_dyad,
	R,
	U_init,
	V_init,
	maxit = 25L,
	tol = 1e-6,
	rel = FALSE
) {
	nr <- nrow(Rbar); nc <- ncol(Rbar)
	U <- U_init; V <- V_init
	Rb0 <- Rbar; Rb0[!obs_dyad] <- 0
	ss0 <- sum(C[obs_dyad] * Rbar[obs_dyad]^2)
	obj_old <- .ae_lowrank_obj(tcrossprod(U, V), Rbar, C, obs_dyad)
	it <- 0L
	for (it in seq_len(maxit)) {
		for (i in seq_len(nr)) {
			ci <- C[i, ]
			if (all(ci == 0)) { U[i, ] <- 0; next }
			U[i, ] <- .ae_safe_solve(
				crossprod(V, ci * V),
				crossprod(V, ci * Rb0[i, ])
			)
		}
		for (j in seq_len(nc)) {
			cj <- C[, j]
			if (all(cj == 0)) { V[j, ] <- 0; next }
			V[j, ] <- .ae_safe_solve(
				crossprod(U, cj * U),
				crossprod(U, cj * Rb0[, j])
			)
		}
		rb <- .ae_rebalance_uv(U, V); U <- rb$U; V <- rb[["V"]]
		obj_new <- .ae_lowrank_obj(tcrossprod(U, V), Rbar, C, obs_dyad)
		# relative change of the objective; a normal response (`rel`) measures it
		# against the objective, floored at 1e-12 of the target's weighted sum of squares
		done <- (obj_old - obj_new) / (if (rel) max(
			obj_old,
			1e-12 * ss0,
			.Machine$double.xmin
		) else max(1, obj_old)) < tol
		obj_old <- obj_new
		if (done) break
	}
	list(
		U = U,
		V = V,
		L = NULL,
		Omat = tcrossprod(U, V),
		obj = obj_old,
		iters = it
	)
}

# dispatch the multiplicative-block update on `method`; "als" and "hybrid" are
# directed / bipartite only, symmetric fits use MM
.ae_lowrank_update <- function(Omat, Rbar, C, obs_dyad, R, symmetric, method, rel = FALSE) {
	if (method == "mm" || symmetric) {
		return(.ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, symmetric, rel = rel))
	}
	# zero factors are a fixed point of ALS, so seed it with one MM step
	seed <- .ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, symmetric, maxit = 1L, rel = rel)
	als <- .ae_lowrank_als(Rbar, C, obs_dyad, R, seed$U, seed[["V"]], rel = rel)
	if (method == "als") return(als)
	# hybrid: also run the full MM and keep the lower objective
	mm <- .ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, symmetric, rel = rel)
	if (als$obj <= .ae_lowrank_obj(mm$Omat, Rbar, C, obs_dyad)) als else mm
}

# penalized (MAP / empirical-Bayes) low-rank solvers for the IRLS families

# penalized ALS with an EM hyper-step (directed / bipartite). each row solve
# carries the prior precision of u_i ~ N(0, diag(sig2_u)):
#   u_i = (sum_j c_ij v_j v_j' + diag(1/sig2_u))^{-1} sum_j c_ij rbar_ij v_j
# the per-column scale is updated once per call from the incoming factors,
# before the ridged solves shrink them (ECM ordering),
#   sig2_r = (suv0 + sum_i u_ir^2 + sum_i Var(u_ir | data)) / (nr + kappa0);
# updating it after every shrink collapses to Suv -> 0 when suv0 is small.
# suv0 is the per-coordinate IW prior scale (1 on the probit scale, the
# moment-scaled vscale of prior$Suv0 for poisson)
.ae_lowrank_em_als <- function(
	Rbar,
	C,
	obs_dyad,
	R,
	U,
	V,
	sig2_u,
	sig2_v,
	suv0,
	kappa0,
	maxit = 15L,
	tol = 1e-6,
	off = which(!obs_dyad)
) {
	# cells off the observed set carry no weight, so zero their target and weight
	Rb0 <- Rbar; C0 <- C
	if (length(off)) {
		if (any(Rb0[off] != 0)) Rb0[off] <- 0
		if (any(C0[off] != 0)) C0[off] <- 0
	}
	# the compiled update also balances the factors when the variance order and the
	# factor order agree and its QR is unpivoted; otherwise the code below does it
	out <- als_factor_balanced_cpp(
		Rb0,
		C0,
		U,
		V,
		sig2_u,
		sig2_v,
		suv0,
		kappa0,
		maxit,
		tol
	)
	balanced_done <- isTRUE(out$balanced)
	out$balanced <- NULL
	if (balanced_done) return(out)
	variance_order <- order(sqrt(out$sig2_u * out$sig2_v), decreasing = TRUE)
	factor_order <- order(sqrt(colSums(out$U^2) * colSums(out[["V"]]^2)), decreasing = TRUE)
	balanced <- if (identical(variance_order, factor_order)) {
		.ae_factor_balance(out$U, out[["V"]], out$sig2_u, out$sig2_v)
	} else {
		# near equal variances, exchanging singular directions can cycle the
		# second-moment update; column rescaling lowers the same penalty without it
		uu <- colSums(out$U^2)
		vv <- colSums(out[["V"]]^2)
		ratio <- rep(1, length(uu))
		positive <- uu > 0 & vv > 0
		ratio[positive] <- (vv[positive] * out$sig2_u[positive] /
				(uu[positive] * out$sig2_v[positive]))^0.25
		# a log-scale ratio where tiny column norms and variances would underflow
		unstable <- positive & (!is.finite(ratio) | ratio <= 0)
		if (any(unstable)) {
			ratio[unstable] <- exp(
				(log(vv[unstable]) +
						log(out$sig2_u[unstable]) - log(uu[unstable]) -
						log(out$sig2_v[unstable])) / 4
			)
		}
		list(
			U = out$U * rep(ratio, each = nrow(out$U)),
			V = out[["V"]] / rep(ratio, each = nrow(out[["V"]]))
		)
	}
	out$U <- balanced$U
	out[["V"]] <- balanced$V
	out$Omat <- tcrossprod(out$U, out[["V"]])
	out
}

# sum over actors of the trace of the posterior covariance of an actor's
# factor row given the others; row i has the Gram
#   sum_j c_ij f_j f_j' + I / sig2
# `fm` is nr x R, `c_mat` the symmetric dyad weights, each unordered dyad
# entering the row Gram once
.ae_sym_trace <- function(fm, c_mat, sig2) {
	R <- ncol(fm)
	ridge <- diag(1 / max(sig2, 1e-8), R)
	total <- 0
	for (i in seq_len(nrow(fm))) {
		g <- crossprod(fm, c_mat[i, ] * fm) + ridge
		total <- total + sum(diag(.ae_safe_inverse(g)))
	}
	total
}

# inverse of a small symmetric positive definite matrix, with the generalised
# inverse when it is numerically singular
.ae_safe_inverse <- function(g) {
	ch <- tryCatch(chol(g), error = function(e) NULL)
	if (!is.null(ch)) return(chol2inv(ch))
	eg <- eigen((g + t(g)) / 2, symmetric = TRUE)
	keep <- eg$values > max(eg$values) * 1e-12
	eg$vectors[, keep, drop = FALSE] %*% (t(eg$vectors[, keep, drop = FALSE]) / eg$values[keep])
}

# penalized MM update for the symmetric eigenmodel: the surrogate's
# eigenvalues are soft-thresholded by lam / cmax (lam = 1 / sig2) before rank-r
# truncation, the exact minimiser under the nuclear-norm penalty the factor
# ridge implies. sig2 is updated once per call from the incoming unpenalized
# eigenvalues (ECM ordering),
#   sig2 = (kappa0 suv0 + 2 sum|l_r|) / (kappa0 + 2 nr R)
# or, with `second_moment`, by the EM update with the posterior second moment,
#   sig2 = (suv0 + ||F||^2 + tr) / (nr R + kappa0),
# where F = U |L|^{1/2} and `tr` comes from .ae_sym_trace at `sig2_in`; this
# keeps L = 0 from being absorbing and is used by the empirical-Bayes fit
# (ame_als_eb.R)
.ae_lowrank_em_mm <- function(
	Omat,
	Rbar,
	C,
	obs_dyad,
	R,
	suv0,
	kappa0,
	maxit = 25L,
	cold_start = TRUE,
	second_moment = FALSE,
	sig2_in = NULL,
	cfloor = 1
) {
	nr <- nrow(Rbar)
	cmax <- max(C); if (cmax < cfloor) cmax <- cfloor
	cmax <- max(cmax, .Machine$double.xmin)
	w_obs <- C[obs_dyad] / cmax

	# E-step: sig2 from the incoming factors before shrinking; the incoming Omat
	# carries the previous BCD iteration's shrunk estimate, and a cold start uses
	# the unpenalized surrogate's spectrum
	if (isTRUE(cold_start) && max(abs(Omat)) < 1e-10) {
		A0 <- Omat
		A0[obs_dyad] <- Omat[obs_dyad] +
			w_obs * (Rbar[obs_dyad] - Omat[obs_dyad])
		A0 <- (A0 + t(A0)) / 2
		eg0 <- eigen(A0, symmetric = TRUE)
	} else {
		eg0 <- eigen((Omat + t(Omat)) / 2, symmetric = TRUE)
	}
	idx0 <- order(abs(eg0$values), decreasing = TRUE)[seq_len(R)]
	lref <- eg0$values[idx0]
	if (second_moment) {
		fref <- sweep(eg0$vectors[, idx0, drop = FALSE], 2L, sqrt(abs(lref)), "*")
		tr_cond <- .ae_sym_trace(fref, C, if (is.null(sig2_in)) suv0 else sig2_in)
		sig2 <- (suv0 + sum(abs(lref)) + tr_cond) / (nr * R + kappa0)
	} else {
		sig2 <- (kappa0 * suv0 + 2 * sum(abs(lref))) / (kappa0 + 2 * nr * R)
	}
	lam <- 1 / max(sig2, 1e-8)

	# M-step: soft-thresholded MM sweeps at the fixed sig2
	U <- L <- NULL
	for (mm in seq_len(maxit)) {
		A <- Omat
		A[obs_dyad] <- Omat[obs_dyad] +
			w_obs * (Rbar[obs_dyad] - Omat[obs_dyad])
		A <- (A + t(A)) / 2
		eg <- eigen(A, symmetric = TRUE)
		idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
		lraw <- eg$values[idx]
		lsh <- sign(lraw) * pmax(abs(lraw) - lam / cmax, 0)
		U <- eg$vectors[, idx, drop = FALSE]
		L <- lsh
		Omat_new <- U %*% (lsh * t(U))
		delta <- max(abs(Omat_new - Omat))
		Omat <- Omat_new
		if (delta < 1e-8) break
	}
	list(
		U = U,
		V = U,
		L = L,
		Omat = Omat,
		sig2_u = rep(sig2, R),
		sig2_v = rep(sig2, R)
	)
}

# dispatch the penalized multiplicative-block update (IRLS families):
# directed / bipartite fits use the ridged ALS solver seeded by one
# unpenalized MM step, symmetric fits the soft-thresholded MM; the incoming
# sig2 only seeds the ALS ridge on the cold-start sweep
.ae_lowrank_update_pen <- function(
	Omat,
	Rbar,
	C,
	obs_dyad,
	R,
	symmetric,
	U,
	V,
	sig2_u,
	sig2_v,
	suv0,
	kappa0,
	off = which(!obs_dyad),
	second_moment = FALSE,
	cfloor = 1
) {
	if (symmetric) {
		return(.ae_lowrank_em_mm(
			Omat,
			Rbar,
			C,
			obs_dyad,
			R,
			suv0,
			kappa0,
			cold_start = is.null(sig2_u),
			second_moment = second_moment,
			sig2_in = if (is.null(sig2_u)) NULL else sig2_u[1L],
			cfloor = cfloor
		))
	}
	if (is.null(U) || is.null(V)) {
		seed <- .ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, FALSE, maxit = 1L, cfloor = cfloor)
		U <- seed$U; V <- seed[["V"]]
	}
	nr <- nrow(Rbar); nc <- ncol(Rbar)
	if (is.null(sig2_u)) sig2_u <- (suv0 + colSums(U^2)) / (nr + kappa0)
	if (is.null(sig2_v)) sig2_v <- (suv0 + colSums(V^2)) / (nc + kappa0)
	.ae_lowrank_em_als(
		Rbar,
		C,
		obs_dyad,
		R,
		U,
		V,
		sig2_u,
		sig2_v,
		suv0,
		kappa0,
		off = off
	)
}

# minimize the factor penalty at a fixed product: pair the largest singular
# values with the largest geometric-mean prior variances, then balance the two
# factors under their own variances, via a small QR core (pivoting undone)
.ae_factor_balance <- function(U, V, sig2_u, sig2_v) {
	rank <- ncol(U)
	qr_u <- qr(U)
	qr_v <- qr(V)
	r_u <- qr.R(qr_u)[, order(qr_u$pivot), drop = FALSE]
	r_v <- qr.R(qr_v)[, order(qr_v$pivot), drop = FALSE]
	spectrum <- svd(tcrossprod(r_u, r_v), nu = rank, nv = rank)
	index <- order(sqrt(sig2_u * sig2_v), decreasing = TRUE)
	ratio <- (sig2_u[index] / sig2_v[index])^0.25
	u <- matrix(0, nrow(U), rank)
	v <- matrix(0, nrow(V), rank)
	qu <- qr.Q(qr_u) %*% spectrum$u
	qv <- qr.Q(qr_v) %*% spectrum$v
	u[, index] <- qu * rep(sqrt(spectrum$d) * ratio, each = nrow(qu))
	v[, index] <- qv * rep(sqrt(spectrum$d) / ratio, each = nrow(qv))
	list(U = u, V = v)
}
