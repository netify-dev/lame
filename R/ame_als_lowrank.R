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
# all three solvers minimise the same weighted low-rank objective
#   sum_ij c_ij (rbar_ij - o_ij)^2 ,  rank(o) <= r
# (the time-averaged form of the multiplicative sub-problem), so each is a
# monotone, objective-non-increasing low-rank update.

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
