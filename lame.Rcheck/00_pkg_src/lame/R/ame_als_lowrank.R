# ame_als_lowrank.r
#
# multiplicative (low-rank) block solvers for the fast ame estimator.
#
# the default is the weighted majorize-minimize (mm) update. `lowrank_method`
# = "als" or "hybrid" select an alternating-least-squares solver that converges
# faster than mm when the per-dyad weights c_ij are very unbalanced (a sparse
# or strongly unbalanced longitudinal panel). als is available for directed and
# bipartite models; symmetric (eigenmodel) fits fall back to mm, where the
# same-actor-on-both-sides indefinite structure makes a plain als unstable.
#
# the unpenalized solvers (.ae_lowrank_mm / .ae_lowrank_als) minimise the same
# weighted low-rank objective
#   sum_ij c_ij (rbar_ij - o_ij)^2 ,  rank(o) <= r
# (the time-averaged form of the multiplicative sub-problem), so each is a
# monotone, objective-non-increasing low-rank update. they are used by the
# normal-family and fixed-transform paths.
#
# the irls families (binary / poisson) instead use the penalized (map)
# solvers further down (.ae_lowrank_em_als / .ae_lowrank_em_mm), which carry
# the ridge implied by the mcmc prior [u_i, v_i] ~ N(0, Suv),
# Suv ~ IW(suv0 I, kappa0), with the per-column prior scale sig2 estimated by
# an em step that includes the posterior-variance trace term. without the
# ridge the rank-r factor block is an unpenalized rank-r glm mle, which for a
# binary likelihood diverges by quasi-separation (the deviance falls
# monotonically as ||uv'|| grows), inflating every coefficient. note the em
# hyper-step optimises a lower bound, so the penalized solvers are monotone
# in the penalized objective only at fixed sig2; the irls outer loop keeps
# the best penalized-deviance iterate, which is safe either way.

# weighted low-rank objective over observed dyads
.ae_lowrank_obj <- function(O, Rbar, C, obs_dyad) {
	sum(C[obs_dyad] * (Rbar[obs_dyad] - O[obs_dyad])^2)
}

# weighted majorize-minimize update (the baseline low-rank algorithm)
.ae_lowrank_mm <- function(Omat, Rbar, C, obs_dyad, R, symmetric,
                           maxit = 25L) {
	cmax <- max(C); if (cmax < 1) cmax <- 1
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
		if (delta < 1e-8) break
	}
	list(U = U, V = V, L = L, Omat = Omat)
}

# rebalance u, v so the two factors carry equal scale, without changing u v'.
# als can otherwise let one factor grow and the other shrink (the u g, v g^-t
# gauge freedom), which makes the per-row systems progressively worse scaled.
.ae_rebalance_uv <- function(U, V) {
	R <- ncol(U)
	if (R == 0L) return(list(U = U, V = V))
	qu <- qr(U); qv <- qr(V)
	sv <- svd(qr.R(qu) %*% t(qr.R(qv)))
	dd <- sqrt(pmax(sv$d, 0))
	list(U = qr.Q(qu) %*% sweep(sv$u, 2, dd, "*"),
	     V = qr.Q(qv) %*% sweep(sv$v, 2, dd, "*"))
}

# weighted alternating least squares (directed / bipartite). holding v fixed,
# each row u_i is the exact weighted-ls minimiser
#   u_i = (sum_j c_ij v_j v_j')^+ (sum_j c_ij rbar_ij v_j) ;
# symmetrically for v. every row/column block is an exact conditional minimiser,
# so the weighted low-rank objective is non-increasing across sweeps.
.ae_lowrank_als <- function(Rbar, C, obs_dyad, R, U_init, V_init,
                            maxit = 25L, tol = 1e-6) {
	nr <- nrow(Rbar); nc <- ncol(Rbar)
	U <- U_init; V <- V_init
	Rb0 <- Rbar; Rb0[!obs_dyad] <- 0
	obj_old <- .ae_lowrank_obj(tcrossprod(U, V), Rbar, C, obs_dyad)
	it <- 0L
	for (it in seq_len(maxit)) {
		for (i in seq_len(nr)) {
			ci <- C[i, ]
			if (all(ci == 0)) { U[i, ] <- 0; next }
			U[i, ] <- .ae_safe_solve(crossprod(V, ci * V),
			                         crossprod(V, ci * Rb0[i, ]))
		}
		for (j in seq_len(nc)) {
			cj <- C[, j]
			if (all(cj == 0)) { V[j, ] <- 0; next }
			V[j, ] <- .ae_safe_solve(crossprod(U, cj * U),
			                         crossprod(U, cj * Rb0[, j]))
		}
		# rebalance the factor scales (leaves u v' unchanged)
		rb <- .ae_rebalance_uv(U, V); U <- rb$U; V <- rb$V
		obj_new <- .ae_lowrank_obj(tcrossprod(U, V), Rbar, C, obs_dyad)
		done <- (obj_old - obj_new) / max(1, obj_old) < tol
		obj_old <- obj_new
		if (done) break
	}
	list(U = U, V = V, L = NULL, Omat = tcrossprod(U, V),
	     obj = obj_old, iters = it)
}

# dispatch the multiplicative-block update on `method`. "als"/"hybrid" are
# directed/bipartite only -- symmetric fits always use mm.
.ae_lowrank_update <- function(Omat, Rbar, C, obs_dyad, R, symmetric, method) {
	if (method == "mm" || symmetric) {
		return(.ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, symmetric))
	}
	# als from o = 0 is a fixed point (zero factors stay zero), so seed it with
	# one mm step for a robust spectral start
	seed <- .ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, symmetric, maxit = 1L)
	als  <- .ae_lowrank_als(Rbar, C, obs_dyad, R, seed$U, seed$V)
	if (method == "als") return(als)
	# hybrid: also run the full mm and keep whichever objective is lower, so the
	# hybrid result is never worse than the baseline mm
	mm <- .ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, symmetric)
	if (als$obj <= .ae_lowrank_obj(mm$Omat, Rbar, C, obs_dyad)) als else mm
}

# --------------------------------------------------------------------------
# penalized (map / empirical-bayes) low-rank solvers for the irls families
# --------------------------------------------------------------------------

# penalized alternating least squares with an em hyper-step (directed /
# bipartite). each row solve carries the prior precision implied by the mcmc
# prior u_i ~ N(0, diag(sig2_u)):
#   u_i = (sum_j c_ij v_j v_j' + diag(1/sig2_u))^{-1} sum_j c_ij rbar_ij v_j
# (the irls working weights c_ij already put the working response on the
# unit-dispersion link scale, so the prior precision enters unscaled).
#
# ordering (ECM): the per-column prior scale is updated once per call, from
# the incoming factors, before the ridged row solves shrink them,
#   sig2_r = (suv0 + sum_i u_ir^2 + sum_i Var(u_ir | data)) / (nr + kappa0),
# where the posterior-variance trace term Var(u_ir | data) is the honest EM
# second moment (a free byproduct of the ridged Gram inverse). the outer BCD
# loop then supplies the EM iteration. this ordering matters: updating sig2
# repeatedly within the block (after every shrink) spirals into the degenerate
# Suv -> 0 fixed point whenever the prior scale suv0 is small -- which crushes
# genuine poisson structure to var(UV') ~ 0 (binary is spared only because its
# suv0 = 1 floors sig2 at a healthy value). seeding sig2 from the unpenalized
# spectral factors and updating it before shrinking keeps both the null case
# (shrinks to ~0, matching mcmc) and the signal case (recovers real structure)
# in the correct basin. suv0 is the per-coordinate iw prior scale (1 on the
# binary probit scale; the moment-scaled vscale used for prior$Suv0 in the
# mcmc poisson path).
.ae_lowrank_em_als <- function(Rbar, C, obs_dyad, R, U, V, sig2_u, sig2_v,
                               suv0, kappa0, maxit = 15L, tol = 1e-6) {
	nr <- nrow(Rbar); nc <- ncol(Rbar)
	Rb0 <- Rbar; Rb0[!obs_dyad] <- 0
	C0 <- C; C0[!obs_dyad] <- 0

	# --- E-step: update sig2 from the incoming (not-yet-shrunk) factors ---
	# the trace uses the ridged Gram at the incoming sig2, the current best
	# estimate of the posterior precision
	Lu <- diag(1 / pmax(sig2_u, 1e-8), R)
	Vu <- rep(0, R)
	for (i in seq_len(nr)) {
		Mi <- tryCatch(solve(crossprod(V, C0[i, ] * V) + Lu),
		               error = function(e) MASS::ginv(crossprod(V, C0[i, ] * V) + Lu))
		Vu <- Vu + diag(Mi)
	}
	sig2_u <- (suv0 + colSums(U^2) + Vu) / (nr + kappa0)
	Lv <- diag(1 / pmax(sig2_v, 1e-8), R)
	Vv <- rep(0, R)
	for (j in seq_len(nc)) {
		Mi <- tryCatch(solve(crossprod(U, C0[, j] * U) + Lv),
		               error = function(e) MASS::ginv(crossprod(U, C0[, j] * U) + Lv))
		Vv <- Vv + diag(Mi)
	}
	sig2_v <- (suv0 + colSums(V^2) + Vv) / (nc + kappa0)

	# --- M-step: penalized (ridged) ALS sweeps at the fixed sig2 above ---
	Lu <- diag(1 / pmax(sig2_u, 1e-8), R)
	Lv <- diag(1 / pmax(sig2_v, 1e-8), R)
	obj_old <- Inf
	for (it in seq_len(maxit)) {
		for (i in seq_len(nr)) {
			ci <- C0[i, ]
			U[i, ] <- .ae_safe_solve(crossprod(V, ci * V) + Lu,
			                         crossprod(V, ci * Rb0[i, ]))
		}
		for (j in seq_len(nc)) {
			cj <- C0[, j]
			V[j, ] <- .ae_safe_solve(crossprod(U, cj * U) + Lv,
			                         crossprod(U, cj * Rb0[, j]))
		}
		O <- tcrossprod(U, V)
		obj_new <- sum(C0[obs_dyad] * (Rbar[obs_dyad] - O[obs_dyad])^2) +
			sum(sweep(U^2, 2, 1 / pmax(sig2_u, 1e-8), "*")) +
			sum(sweep(V^2, 2, 1 / pmax(sig2_v, 1e-8), "*"))
		done <- is.finite(obj_old) &&
			abs(obj_old - obj_new) / max(1, abs(obj_old)) < tol
		obj_old <- obj_new
		if (done) break
	}
	list(U = U, V = V, L = NULL, Omat = tcrossprod(U, V),
	     sig2_u = sig2_u, sig2_v = sig2_v, obj = obj_old, iters = maxit)
}

# penalized majorize-minimize update for the symmetric eigenmodel. the
# penalized rank-r step under the mm surrogate is soft-thresholding: the
# eigenvalues of the surrogate are shrunk toward zero by lam/cmax with
# lam = 1/sig2 before rank-r truncation (the exact minimizer of the
# majorized objective plus the nuclear-norm penalty implied by the factor
# ridge, since min over factorizations of the ridge is 2*lam*||O||_*). the
# isotropic prior scale sig2 is updated by the iw posterior-mode analog
#   sig2 = (kappa0*suv0 + 2*sum|l_r|) / (kappa0 + 2*nr*R)
# using the ECM ordering (as for the directed als path): sig2 is updated
# once per call, from the incoming *unpenalized* eigenvalues, before the
# soft-threshold sweeps shrink them. this avoids the degenerate Suv -> 0
# spiral that repeated within-block updates hit when suv0 is small (poisson).
.ae_lowrank_em_mm <- function(Omat, Rbar, C, obs_dyad, R,
                              suv0, kappa0, maxit = 25L) {
	nr <- nrow(Rbar)
	cmax <- max(C); if (cmax < 1) cmax <- 1
	w_obs <- C[obs_dyad] / cmax

	# --- E-step: sig2 from the incoming factors, before shrinking ---
	# once past the cold start, the incoming Omat is the previous BCD
	# iteration's *shrunk* low-rank estimate, so its spectrum carries the
	# adapted prior scale across BCD iterations (as the directed ALS path
	# threads sig2 via the incoming factors); on a cold start (Omat == 0) we
	# fall back to the unpenalized surrogate's spectrum for a spectral start.
	if (max(abs(Omat)) < 1e-10) {
		A0 <- Omat
		A0[obs_dyad] <- Omat[obs_dyad] +
			w_obs * (Rbar[obs_dyad] - Omat[obs_dyad])
		A0 <- (A0 + t(A0)) / 2
		eg0 <- eigen(A0, symmetric = TRUE)
		idx0 <- order(abs(eg0$values), decreasing = TRUE)[seq_len(R)]
		lref <- eg0$values[idx0]
	} else {
		eg0 <- eigen((Omat + t(Omat)) / 2, symmetric = TRUE)
		idx0 <- order(abs(eg0$values), decreasing = TRUE)[seq_len(R)]
		lref <- eg0$values[idx0]
	}
	sig2 <- (kappa0 * suv0 + 2 * sum(abs(lref))) / (kappa0 + 2 * nr * R)
	lam <- 1 / max(sig2, 1e-8)

	# --- M-step: soft-thresholded MM sweeps at the fixed sig2 above ---
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
	list(U = U, V = U, L = L, Omat = Omat,
	     sig2_u = rep(sig2, R), sig2_v = rep(sig2, R))
}

# dispatch the *penalized* multiplicative-block update (irls families).
# directed / bipartite fits use the ridged als row solver seeded by one
# unpenalized mm step (a spectral start: als from zero factors is a fixed
# point, and the prior scale must be set from unshrunk factors); symmetric
# fits use the soft-thresholded mm. both use the ECM ordering internally
# (update the prior scale from the incoming, unshrunk factors before
# shrinking), so the incoming sig2 state only seeds the als ridge for the
# very first (cold-start) sweep.
.ae_lowrank_update_pen <- function(Omat, Rbar, C, obs_dyad, R, symmetric,
                                   U, V, sig2_u, sig2_v, suv0, kappa0) {
	if (symmetric) {
		return(.ae_lowrank_em_mm(Omat, Rbar, C, obs_dyad, R, suv0, kappa0))
	}
	if (is.null(U) || is.null(V)) {
		seed <- .ae_lowrank_mm(Omat, Rbar, C, obs_dyad, R, FALSE, maxit = 1L)
		U <- seed$U; V <- seed$V
	}
	nr <- nrow(Rbar); nc <- ncol(Rbar)
	if (is.null(sig2_u)) sig2_u <- (suv0 + colSums(U^2)) / (nr + kappa0)
	if (is.null(sig2_v)) sig2_v <- (suv0 + colSums(V^2)) / (nc + kappa0)
	.ae_lowrank_em_als(Rbar, C, obs_dyad, R, U, V, sig2_u, sig2_v,
	                   suv0, kappa0)
}
