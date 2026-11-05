# static low-rank block of the penalised dynamic fit under a ridge on its factors
# the ridge keeps the binary and sparse count objectives bounded below (quasi-separation)

# penalised weighted low-rank objective of the product U V'
.lda_lowrank_obj_ridge <- function(U, V, Rbar, C, obs_dyad, ridge) {
	.ae_lowrank_obj(tcrossprod(U, V), Rbar, C, obs_dyad) + ridge * (sum(U^2) + sum(V^2))
}

# soft-thresholded weighted MM update; at a fixed product the ridge is a nuclear-norm penalty
.lda_lowrank_mm_ridge <- function(
	Omat,
	Rbar,
	C,
	obs_dyad,
	R,
	symmetric,
	ridge,
	maxit = 25L
) {
	cmax <- max(C)
	if (cmax < 1) cmax <- 1
	w_obs <- C[obs_dyad] / cmax
	thr <- ridge / cmax
	U <- V <- L <- NULL
	for (mm in seq_len(maxit)) {
		A <- Omat
		A[obs_dyad] <- Omat[obs_dyad] + w_obs * (Rbar[obs_dyad] - Omat[obs_dyad])
		if (symmetric) {
			A <- (A + t(A)) / 2
			eg <- eigen(A, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
			U <- eg$vectors[, idx, drop = FALSE]
			lraw <- eg$values[idx]
			L <- sign(lraw) * pmax(abs(lraw) - thr, 0)
			V <- U
			Omat_new <- U %*% (L * t(U))
		} else {
			sv <- svd(A)
			dd <- sqrt(pmax(sv$d[seq_len(R)] - thr, 0))
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

# weighted alternating ridge regressions (directed), rebalancing the factors after each sweep
.lda_lowrank_als_ridge <- function(
	Rbar,
	C,
	obs_dyad,
	R,
	U_init,
	V_init,
	ridge,
	maxit = 25L,
	tol = 1e-6
) {
	nr <- nrow(Rbar)
	nc <- ncol(Rbar)
	U <- U_init
	V <- V_init
	Rb0 <- Rbar
	Rb0[!obs_dyad] <- 0
	rid <- diag(ridge, R)
	obj_old <- .lda_lowrank_obj_ridge(U, V, Rbar, C, obs_dyad, ridge)
	it <- 0L
	for (it in seq_len(maxit)) {
		for (i in seq_len(nr)) {
			ci <- C[i, ]
			if (all(ci == 0)) { U[i, ] <- 0; next }
			U[i, ] <- .ae_safe_solve(crossprod(V, ci * V) + rid, crossprod(V, ci * Rb0[i, ]))
		}
		for (j in seq_len(nc)) {
			cj <- C[, j]
			if (all(cj == 0)) { V[j, ] <- 0; next }
			V[j, ] <- .ae_safe_solve(crossprod(U, cj * U) + rid, crossprod(U, cj * Rb0[, j]))
		}
		rb <- .ae_rebalance_uv(U, V)
		U <- rb$U
		V <- rb[["V"]]
		obj_new <- .lda_lowrank_obj_ridge(U, V, Rbar, C, obs_dyad, ridge)
		done <- (obj_old - obj_new) / max(1, obj_old) < tol
		obj_old <- obj_new
		if (done) break
	}
	list(U = U, V = V, L = NULL, Omat = tcrossprod(U, V), obj = obj_old, iters = it)
}

# ridge version of .ae_lowrank_update(): "mm" and the eigenmodel soft-threshold,
# "als" runs the ridge regressions from one such step, "hybrid" keeps the better
.lda_lowrank_update_ridge <- function(Omat, Rbar, C, obs_dyad, R, symmetric, method, ridge) {
	if (method == "mm" || symmetric) {
		return(.lda_lowrank_mm_ridge(Omat, Rbar, C, obs_dyad, R, symmetric, ridge))
	}
	seed <- .lda_lowrank_mm_ridge(Omat, Rbar, C, obs_dyad, R, symmetric, ridge, maxit = 1L)
	als <- .lda_lowrank_als_ridge(Rbar, C, obs_dyad, R, seed$U, seed[["V"]], ridge)
	if (method == "als") return(als)
	mm <- .lda_lowrank_mm_ridge(Omat, Rbar, C, obs_dyad, R, symmetric, ridge)
	if (als$obj <= .lda_lowrank_obj_ridge(mm$U, mm[["V"]], Rbar, C, obs_dyad, ridge)) als else mm
}
