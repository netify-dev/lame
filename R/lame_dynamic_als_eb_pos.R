# moving latent positions (dynamic_uv = TRUE) of the dynamic binary empirical-Bayes
# estimator of lame_dynamic_als_eb.R. each coordinate follows the stationary AR(1) of lame(),
#   u_i1 ~ N(0, S), u_it = rho u_i,t-1 + e_it, e_it ~ N(0, s2), S = s2 / (1 - rho^2),
# with s2 ~ IG(2, 1) and rho ~ TN(0.9, 0.1); an undirected network has one position per actor
# and the product F S F' with the start's sign pattern. point estimates of a bilinear term
# overfit a sparse binary panel, so the positions get a mean-field Gaussian posterior
# q(U) q(V) with one full path covariance per actor. the two sides are updated in turn and
# balanced along the gauge after each sweep, and the product's variance enters the additive
# block's working problem by Gauss-Hermite quadrature (.lde_work_gh()). the stationary
# variance S minimises a profile objective by search over log S, since EM on the second
# moments of q collapses the scale; rho is the EM update at the held S

# hyperprior and the AR(1) precision

# hyperprior of the positions as lame() sets it; the user's `prior` list overrides the defaults
.lde_pos_prior <- function(prior = list(), gaussian_bip_ar = FALSE, scale_default = 1) {
	prior <- if (is.list(prior)) prior else list()
	# a normal bipartite panel with the Gaussian AR(1) transition takes the sampler's wider,
	# lower persistence prior (its interaction matrix carries the scale)
	list(
		rho_mean = .lda_default(prior$rho_uv_mean, if (gaussian_bip_ar) 0.5 else 0.9),
		rho_sd = .lda_default(prior$rho_uv_sd, if (gaussian_bip_ar) 0.2 else 0.1),
		shape = .lda_default(prior$sigma_uv_shape, 2),
		scale = .lda_default(prior$sigma_uv_scale, scale_default),
		# the support of the sampler's truncated normal
		rho_lo = -0.99,
		rho_hi = 0.99
	)
}

# the hyperparameter record of the positions at stationary variance S and persistence rho
.lde_pos_hp <- function(S, rho) list(S = S, rho = rho, s2 = S * (1 - rho^2))

# precision of a stationary AR(1) of length Tt with unit innovation variance
.lde_pos_q <- function(Tt, rho) {
	if (Tt == 1L) return(matrix(1 - rho^2, 1L, 1L))
	Q <- matrix(0, Tt, Tt)
	diag(Q) <- c(1, rep(1 + rho^2, max(Tt - 2L, 0L)), 1)
	k <- seq_len(Tt - 1L)
	Q[cbind(k, k + 1L)] <- -rho
	Q[cbind(k + 1L, k)] <- -rho
	Q
}

# the same for a scale-mixture transition: the step into period t (t >= 2) has the weight w[t] (the expected mixing
# weight, the innovation variance being s2 / w[t]); with all weights one it is .lde_pos_q(). tridiagonal:
#   d_1 = (1 - rho^2) + rho^2 w_2 ,  d_t = w_t + rho^2 w_t+1 (2 <= t < Tt) ,  d_Tt = w_Tt ,  q_t,t+1 = -rho w_t+1 ,
# and its determinant is (1 - rho^2) times the product of the weights
.lde_pos_qw <- function(Tt, rho, w) {
	if (Tt == 1L) return(matrix(1 - rho^2, 1L, 1L))
	wt <- w[-1L]
	d <- numeric(Tt)
	d[1L] <- (1 - rho^2) + rho^2 * wt[1L]
	if (Tt > 2L) d[2:(Tt - 1L)] <- wt[1:(Tt - 2L)] + rho^2 * wt[2:(Tt - 1L)]
	d[Tt] <- wt[Tt - 1L]
	Q <- diag(d, Tt)
	k <- seq_len(Tt - 1L)
	Q[cbind(k, k + 1L)] <- -rho * wt
	Q[cbind(k + 1L, k)] <- -rho * wt
	Q
}

# the diagonal (n x Tt) and first off-diagonal (n x Tt-1, between t and t+1) of the weighted precision of every
# actor's path, for the weights `lam` (n x Tt; NULL: all one)
.lde_pos_qcoef <- function(n, Tt, rho, lam = NULL) {
	if (is.null(lam)) lam <- matrix(1, n, Tt)
	d <- matrix(0, n, Tt)
	e <- matrix(0, n, max(Tt - 1L, 1L))
	if (Tt == 1L) {
		d[, 1L] <- 1 - rho^2
		return(list(d = d, e = e))
	}
	d[, 1L] <- (1 - rho^2) + rho^2 * lam[, 2L]
	if (Tt > 2L) d[, 2:(Tt - 1L)] <- lam[, 2:(Tt - 1L), drop = FALSE] + rho^2 * lam[, 3:Tt, drop = FALSE]
	d[, Tt] <- lam[, Tt]
	e[] <- -rho * lam[, -1L, drop = FALSE]
	list(d = d, e = e)
}

# position update of one side

# q of one side's positions given the other's. W and Rr are per-period lists of the (n x m)
# weights and working residuals (zero off the observed cells), Vh and Om the other side's
# (m x R) means and (R x R x m) covariances. returns the means (n x R x Tt), the diagonal
# covariance blocks of each actor's stacked path (R x R x n x Tt), the lag-one blocks
# (R x R x n x (Tt - 1)) and the KL pieces sum_i log|P_i| (ldp), sum_i tr(P0 Sigma_i) (trp)
# and n log|P0| (ldp0)
.lde_pos_side <- function(W, Rr, Vh, Om, s2, rho, n, R, Tt, lam = NULL) {
	P0 <- kronecker(.lde_pos_q(Tt, rho) / s2, diag(R))
	Lam <- array(0, dim = c(R, R, n, Tt))
	H <- array(0, dim = c(n, R, Tt))
	for (t in seq_len(Tt)) {
		for (k in seq_len(R)) for (l in k:R) {
			ev <- Vh[[t]][, k] * Vh[[t]][, l] + Om[[t]][k, l, ]
			v <- as.vector(W[[t]] %*% ev)
			Lam[k, l, , t] <- v
			Lam[l, k, , t] <- v
		}
		H[, , t] <- (W[[t]] * Rr[[t]]) %*% Vh[[t]]
	}
	mean_ <- array(0, dim = c(n, R, Tt))
	sig <- array(0, dim = c(R, R, n, Tt))
	cl <- array(0, dim = c(R, R, n, max(Tt - 1L, 1L)))
	ldp <- 0
	trp <- 0
	ldp0 <- n * as.numeric(determinant(P0)$modulus)
	for (i in seq_len(n)) {
		# the prior precision of this actor's path (its transition weights, if the transition is a scale mixture)
		P0i <- if (is.null(lam)) P0 else kronecker(.lde_pos_qw(Tt, rho, lam[i, ]) / s2, diag(R))
		Pi <- P0i
		for (t in seq_len(Tt)) {
			ix <- (t - 1L) * R + seq_len(R)
			Pi[ix, ix] <- Pi[ix, ix] + Lam[, , i, t]
		}
		ch <- .lde_chol(Pi)
		ldp <- ldp + 2 * sum(log(diag(ch)))
		Sg <- chol2inv(ch)
		trp <- trp + sum(P0i * Sg)
		mean_[i, , ] <- matrix(as.vector(Sg %*% as.vector(H[i, , ])), R, Tt)
		for (t in seq_len(Tt)) {
			ix <- (t - 1L) * R + seq_len(R)
			sig[, , i, t] <- Sg[ix, ix]
			if (t < Tt) cl[, , i, t] <- Sg[ix, ix + R]
		}
	}
	list(mean = mean_, sig = sig, cl = cl, ldp = ldp, trp = trp, ldp0 = ldp0, P0 = P0)
}

# one sweep over both sides: q(U) given q(V), then q(V) given the fresh q(U), by the compiled
# lde_pos_sweep_cpp (src/dynamic_eb_positions.cpp)
.lde_pos_sweep <- function(ctx, pos, hp, wk, eta_add) {
	if (isTRUE(ctx$sym)) return(.lde_pos_sweep_sym(ctx, pos, hp, wk, eta_add))
	out <- lde_pos_sweep_cpp(
		wk$w,
		wk$z,
		eta_add,
		ctx$obs,
		pos$qu$mean,
		pos$qu$sig,
		pos$qv$mean,
		pos$qv$sig,
		ctx$na,
		ctx$nb,
		ctx$R,
		ctx$Tt,
		hp$s2,
		hp$rho,
		if (is.null(pos$qu$lam)) numeric(0) else pos$qu$lam,
		if (is.null(pos$qv$lam)) numeric(0) else pos$qv$lam
	)
	# the transition weights are those the sweep was taken at
	out$qu$lam <- pos$qu$lam
	out$qv$lam <- pos$qv$lam
	out
}

# R reference of the sweep, used by the tests
.lde_pos_sweep_ref <- function(ctx, pos, hp, wk, eta_add) {
	if (isTRUE(ctx$sym)) return(.lde_pos_sweep_sym(ctx, pos, hp, wk, eta_add))
	Tt <- ctx$Tt; na <- ctx$na; nb <- ctx$nb; R <- ctx$R
	obs <- ctx$obs
	W <- lapply(seq_len(Tt), function(t) {
		w <- wk$w[, , t]
		w[!obs[, , t] | !is.finite(w)] <- 0
		w
	})
	Rr <- lapply(seq_len(Tt), function(t) {
		r <- wk$z[, , t] - eta_add[, , t]
		r[!obs[, , t] | !is.finite(r)] <- 0
		r
	})
	Vh <- lapply(seq_len(Tt), function(t) matrix(pos$qv$mean[, , t], nb, R))
	Om <- lapply(seq_len(Tt), function(t) array(pos$qv$sig[, , , t], c(R, R, nb)))
	qu <- .lde_pos_side(W, Rr, Vh, Om, hp$s2, hp$rho, na, R, Tt, lam = pos$qu$lam)
	Uh <- lapply(seq_len(Tt), function(t) matrix(qu$mean[, , t], na, R))
	Omu <- lapply(seq_len(Tt), function(t) array(qu$sig[, , , t], c(R, R, na)))
	qv <- .lde_pos_side(lapply(W, t), lapply(Rr, t), Uh, Omu, hp$s2, hp$rho, nb, R, Tt, lam = pos$qv$lam)
	qu$lam <- pos$qu$lam
	qv$lam <- pos$qv$lam
	list(qu = qu, qv = qv)
}

# undirected networks

# the receiver-side q of an undirected network: V = F S, the means and the covariances
# of the coordinates carry the signs of the product's fixed sign pattern
.lde_pos_mirror <- function(qs, sg) {
	R <- length(sg)
	mean_ <- qs$mean
	for (k in seq_len(R)) mean_[, k, ] <- sg[k] * qs$mean[, k, ]
	sgm <- outer(sg, sg)
	sig <- qs$sig
	cl <- qs$cl
	for (k in seq_len(R)) for (l in seq_len(R)) {
		sig[k, l, , ] <- sgm[k, l] * qs$sig[k, l, , ]
		cl[k, l, , ] <- sgm[k, l] * qs$cl[k, l, , ]
	}
	list(mean = mean_, sig = sig, cl = cl, lam = qs$lam)
}

# one sweep of an undirected network: the positions of an actor enter the cells of its
# row only (every dyad once), against the other actors' current q, so the actors are
# updated one after the other (a coordinate ascent over the actors)
.lde_pos_sweep_sym <- function(ctx, pos, hp, wk, eta_add) {
	Tt <- ctx$Tt; n <- ctx$na; R <- ctx$R
	sg <- ctx$sgn
	W <- wk$w
	W[!ctx$obs | !is.finite(W)] <- 0
	Rr <- wk$z - eta_add
	Rr[!ctx$obs | !is.finite(Rr)] <- 0
	qs <- pos$qu
	P0 <- kronecker(.lde_pos_q(Tt, hp$rho) / hp$s2, diag(R))
	# the second moments of the mapped coordinates s_k u_k of every actor, by pair and period
	pairs <- which(upper.tri(diag(R), diag = TRUE), arr.ind = TRUE)
	EV <- lapply(seq_len(nrow(pairs)), function(m) {
		k <- pairs[m, 1L]; l <- pairs[m, 2L]
		sg[k] * sg[l] * (matrix(qs$mean[, k, ], n, Tt) * matrix(qs$mean[, l, ], n, Tt) + matrix(qs$sig[k, l, , ], n, Tt))
	})
	MS <- lapply(seq_len(R), function(k) sg[k] * matrix(qs$mean[, k, ], n, Tt))
	ldp <- 0
	trp <- 0
	for (i in seq_len(n)) {
		Wi <- matrix(W[i, , ], n, Tt)
		Ri <- matrix(Rr[i, , ], n, Tt)
		P0i <- if (is.null(qs$lam)) P0 else kronecker(.lde_pos_qw(Tt, hp$rho, qs$lam[i, ]) / hp$s2, diag(R))
		Pi <- P0i
		h <- numeric(Tt * R)
		for (m in seq_len(nrow(pairs))) {
			k <- pairs[m, 1L]; l <- pairs[m, 2L]
			lam <- colSums(Wi * EV[[m]])
			for (t in seq_len(Tt)) {
				a <- (t - 1L) * R + k; b <- (t - 1L) * R + l
				Pi[a, b] <- Pi[a, b] + lam[t]
				if (k != l) Pi[b, a] <- Pi[b, a] + lam[t]
			}
		}
		for (k in seq_len(R)) {
			hk <- colSums(Wi * Ri * MS[[k]])
			h[(seq_len(Tt) - 1L) * R + k] <- hk
		}
		ch <- .lde_chol(Pi)
		ldp <- ldp + 2 * sum(log(diag(ch)))
		Sg <- chol2inv(ch)
		trp <- trp + sum(P0i * Sg)
		mu <- matrix(as.vector(Sg %*% h), R, Tt)
		qs$mean[i, , ] <- mu
		for (t in seq_len(Tt)) {
			ix <- (t - 1L) * R + seq_len(R)
			qs$sig[, , i, t] <- Sg[ix, ix]
			if (t < Tt) qs$cl[, , i, t] <- Sg[ix, ix + R]
		}
		for (m in seq_len(nrow(pairs))) {
			k <- pairs[m, 1L]; l <- pairs[m, 2L]
			EV[[m]][i, ] <- sg[k] * sg[l] * (mu[k, ] * mu[l, ] + qs$sig[k, l, i, ])
		}
		for (k in seq_len(R)) MS[[k]][i, ] <- sg[k] * mu[k, ]
	}
	qs$ldp <- ldp
	qs$trp <- trp
	qs$ldp0 <- n * as.numeric(determinant(P0)$modulus)
	qs$P0 <- P0
	list(qu = qs, qv = .lde_pos_mirror(qs, sg))
}

# gauge: the balance of the two sides' scales

# the likelihood sees the positions only through the product and its variance, which
# u_it -> B u_it, v_jt -> B^-T v_jt leaves unchanged; the prior does not, and the sweeps alone
# move along that valley very slowly. with G = B'B and M_u, M_v the prior-weighted second
# moments of the two sides, the KL along the valley is, up to a constant,
#     tr(G M_u) + tr(G^-1 M_v) - T (n_u - n_v) log det G ,
# minimised at G = M_u^-1/2 H M_u^-1/2, H = c/2 I + (M_u^1/2 M_v M_u^1/2 + c^2/4 I)^1/2,
# c = T (n_u - n_v). the state moves there at once (B = G^1/2 for the senders, G^-1/2 for the
# receivers), which keeps the fixed points of the sweeps; rotations and undirected networks
# are left alone
.lde_pos_moments <- function(qs, rho, s2) {
	n <- dim(qs$mean)[1L]; R <- dim(qs$mean)[2L]; Tt <- dim(qs$mean)[3L]
	# the entries of the (weighted) precision of every actor's path, up to the innovation variance
	qc <- .lde_pos_qcoef(n, Tt, rho, qs$lam)
	Cm <- matrix(0, R, R)
	Mn <- matrix(0, R, R)
	for (t in seq_len(Tt)) {
		mt <- matrix(qs$mean[, , t], n, R)
		Cm <- Cm + matrix(matrix(qs$sig[, , , t], R * R, n) %*% qc$d[, t], R, R)
		Mn <- Mn + crossprod(mt * qc$d[, t], mt)
		if (t < Tt) {
			m1 <- matrix(qs$mean[, , t + 1L], n, R)
			cl <- matrix(matrix(qs$cl[, , , t], R * R, n) %*% qc$e[, t], R, R)
			Cm <- Cm + cl + t(cl)
			Mn <- Mn + crossprod(mt * qc$e[, t], m1) + crossprod(m1 * qc$e[, t], mt)
		}
	}
	list(C = (Cm + t(Cm)) / (2 * s2), M = (Cm + Mn + t(Cm + Mn)) / (2 * s2))
}

# the matrix G of the valley's minimum for the two sides' moments
.lde_pos_gauge_g <- function(Mu, Mv, cc) {
	eu <- eigen(Mu, symmetric = TRUE)
	Mh <- eu$vectors %*% (sqrt(eu$values) * t(eu$vectors))
	Mih <- eu$vectors %*% (t(eu$vectors) / sqrt(eu$values))
	K <- Mh %*% Mv %*% Mh
	ek <- eigen((K + t(K)) / 2, symmetric = TRUE)
	H <- ek$vectors %*% ((cc / 2 + sqrt(pmax(ek$values + cc^2 / 4, 0))) * t(ek$vectors))
	G <- Mih %*% H %*% Mih
	(G + t(G)) / 2
}

# B q for the blocks of one side (means n x R x T, covariance blocks R x R x n x T)
.lde_pos_apply <- function(qs, B) {
	n <- dim(qs$mean)[1L]; R <- dim(qs$mean)[2L]
	for (t in seq_len(dim(qs$mean)[3L])) qs$mean[, , t] <- matrix(qs$mean[, , t], n, R) %*% t(B)
	both <- function(a) {
		d <- dim(a)
		a1 <- array(B %*% matrix(a, R, prod(d[-1L])), d)
		a2 <- aperm(a1, c(2L, 1L, 3L, 4L))
		aperm(array(B %*% matrix(a2, R, prod(d[-1L])), d), c(2L, 1L, 3L, 4L))
	}
	qs$sig <- both(qs$sig)
	qs$cl <- both(qs$cl)
	qs
}

.lde_pos_gauge <- function(pos, hp, R, Tt) {
	mu <- .lde_pos_moments(pos$qu, hp$rho, hp$s2)
	mv <- .lde_pos_moments(pos$qv, hp$rho, hp$s2)
	nu <- dim(pos$qu$mean)[1L]; nv <- dim(pos$qv$mean)[1L]
	G <- .lde_pos_gauge_g(mu$M, mv$M, Tt * (nu - nv))
	eg <- eigen(G, symmetric = TRUE)
	Bu <- eg$vectors %*% (sqrt(eg$values) * t(eg$vectors))
	Bv <- eg$vectors %*% (t(eg$vectors) / sqrt(eg$values))
	ldg <- sum(log(eg$values))
	qu <- .lde_pos_apply(pos$qu, Bu)
	qv <- .lde_pos_apply(pos$qv, Bv)
	# the KL pieces that depend on the gauge: tr(P0 Sigma) and the log determinant of the precision
	qu$trp <- sum(G * mu$C); qu$ldp <- pos$qu$ldp - nu * Tt * ldg
	qv$trp <- sum(solve(G) * mv$C); qv$ldp <- pos$qv$ldp + nv * Tt * ldg
	list(qu = qu, qv = qv)
}

# moments of the product

# E_q[u_it' v_jt] for every cell: na x nb x Tt
.lde_pos_product <- function(pos, Tt) {
	na <- dim(pos$qu$mean)[1L]; nb <- dim(pos$qv$mean)[1L]
	O <- array(0, dim = c(na, nb, Tt))
	for (t in seq_len(Tt)) {
		O[, , t] <- tcrossprod(matrix(pos$qu$mean[, , t], na), matrix(pos$qv$mean[, , t], nb))
	}
	O
}

# var_q(u_it' v_jt) = sum_kl (uhat_k uhat_l Sv_kl + Su_kl vhat_k vhat_l + Su_kl Sv_kl), one product of
# (cells by R^2) matrices per period: with the pairs (k, l) as columns it is
#   [Su_kl | uhat_k uhat_l] [Sv_kl + vhat_k vhat_l | Sv_kl]'
.lde_pos_prodvar <- function(pos, R, Tt) {
	Um <- pos$qu$mean; Vm <- pos$qv$mean
	na <- dim(Um)[1L]; nb <- dim(Vm)[1L]
	out <- array(0, dim = c(na, nb, Tt))
	kk <- rep(seq_len(R), times = R)
	ll <- rep(seq_len(R), each = R)
	for (t in seq_len(Tt)) {
		su <- t(matrix(pos$qu$sig[, , , t], R * R, na))
		sv <- t(matrix(pos$qv$sig[, , , t], R * R, nb))
		Ut <- matrix(Um[, , t], na, R)
		Vt <- matrix(Vm[, , t], nb, R)
		out[, , t] <- tcrossprod(
			cbind(su, Ut[, kk, drop = FALSE] * Ut[, ll, drop = FALSE]),
			cbind(sv + Vt[, kk, drop = FALSE] * Vt[, ll, drop = FALSE], sv)
		)
	}
	out
}

# R reference, one outer product per pair and period (used by the tests)
.lde_pos_prodvar_ref <- function(pos, R, Tt) {
	Um <- pos$qu$mean; Vm <- pos$qv$mean
	na <- dim(Um)[1L]; nb <- dim(Vm)[1L]
	out <- array(0, dim = c(na, nb, Tt))
	for (t in seq_len(Tt)) {
		for (k in seq_len(R)) for (l in seq_len(R)) {
			su <- pos$qu$sig[k, l, , t]; sv <- pos$qv$sig[k, l, , t]
			out[, , t] <- out[, , t] + outer(su, sv) + outer(Um[, k, t] * Um[, l, t], sv) +
				outer(su, Vm[, k, t] * Vm[, l, t])
		}
	}
	out
}

# hyperparameters

# expected sufficient statistics of the position paths of both sides (the sums of
# squares and lag-one products of the coordinates under q)
.lde_pos_stats <- function(pos, R, Tt, sym = FALSE) {
	stack <- function(qs) {
		n <- dim(qs$mean)[1L]
		d <- matrix(0, n * R, Tt)
		V <- d
		C <- matrix(0, n * R, Tt - 1L)
		for (k in seq_len(R)) {
			rows <- (k - 1L) * n + seq_len(n)
			d[rows, ] <- qs$mean[, k, ]
			for (t in seq_len(Tt)) V[rows, t] <- qs$sig[k, k, , t]
			for (t in seq_len(Tt - 1L)) C[rows, t] <- qs$cl[k, k, , t]
		}
		# the transition weights of a scale-mixture transition, once for every coordinate of an actor
		list(d = d, V = V, C = C, w = if (!is.null(qs$lam)) qs$lam[rep(seq_len(n), R), , drop = FALSE])
	}
	a <- stack(pos$qu)
	if (sym) return(.lde_ar_stats_core(a$d, a$V, a$C, a$w))
	b <- stack(pos$qv)
	.lde_ar_stats_core(rbind(a$d, b$d), rbind(a$V, b$V), rbind(a$C, b$C), if (!is.null(a$w)) rbind(a$w, b$w))
}

# EM update of the persistence at the held stationary variance S (s2 = S (1 - rho^2) follows
# rho), under the hyperprior of the positions
.lde_pos_rho <- function(St, S, pr, rho_ref = NULL) {
	cnt <- St$nser * St$Tt
	obj <- function(rho) {
		s2 <- S * (1 - rho^2)
		ll <- -(cnt / 2) * log(S) - ((cnt - St$nser) / 2) * log(1 - rho^2) -
			(St$S1 / (2 * S) + (St$Sn - 2 * rho * St$Sc + rho^2 * St$Sl) / (2 * S * (1 - rho^2)))
		ll - (pr$shape + 1) * log(s2) - pr$scale / s2 - (rho - pr$rho_mean)^2 / (2 * pr$rho_sd^2)
	}
	# a coarse grid (steps of 0.025) locates the maximum before refining, so a second local
	# maximum cannot capture the optimiser
	grid <- seq(pr$rho_lo, pr$rho_hi, length.out = 81L)
	k <- which.max(vapply(grid, obj, numeric(1)))
	lo <- grid[max(k - 1L, 1L)]
	hi <- grid[min(k + 1L, length(grid))]
	opt <- stats::optimize(obj, c(lo, hi), maximum = TRUE, tol = 1e-8)
	# the gain, on the deviance scale, of the update over the reference persistence
	gain <- if (is.null(rho_ref)) 0 else 2 * (opt$objective - obj(rho_ref))
	list(rho = opt$maximum, gain = gain)
}

# the expected squared innovation of every actor and period, E_q |u_t - rho u_t-1|^2 (n x Tt - 1), and the
# expected squared first position, E_q |u_1|^2 (n), from the means, the variances and the lag-one
# covariances of q (the diagonals of its blocks)
.lde_pos_eesq <- function(qs, rho) {
	n <- dim(qs$mean)[1L]; R <- dim(qs$mean)[2L]; Tt <- dim(qs$mean)[3L]
	e1 <- numeric(n)
	ee <- matrix(0, n, max(Tt - 1L, 1L))
	for (k in seq_len(R)) {
		m <- matrix(qs$mean[, k, ], n, Tt)
		v <- matrix(qs$sig[k, k, , ], n, Tt)
		e1 <- e1 + m[, 1L]^2 + v[, 1L]
		if (Tt > 1L) {
			cv <- matrix(qs$cl[k, k, , ], n, Tt - 1L)
			ee <- ee + (m[, -1L, drop = FALSE] - rho * m[, -Tt, drop = FALSE])^2 + v[, -1L, drop = FALSE] +
				rho^2 * v[, -Tt, drop = FALSE] - 2 * rho * cv
		}
	}
	list(e1 = e1, ee = ee)
}

# scale-mixture transition: the step into period t has innovation variance s2 / lambda_it with
# lambda_it ~ Gamma(nu / 2, nu / 2), a Student-t transition (the first period is the Gaussian
# stationary law). q(lambda_it) = Gamma((nu + R) / 2, b_it), b_it = (nu + E_q |step|^2 / s2) / 2,
# so the path-precision weight is E lambda_it = (nu + R) / (nu + E |step|^2 / s2) and
# E log lambda_it = digamma((nu + R) / 2) - log b_it enters the objective. returns the weights
# (n x Tt, one in the first period), E log lambda and b
.lde_pos_lambda <- function(qs, hp, nu, R) {
	n <- dim(qs$mean)[1L]; Tt <- dim(qs$mean)[3L]
	lam <- matrix(1, n, Tt)
	elog <- matrix(0, n, Tt)
	b <- matrix(NA_real_, n, Tt)
	if (Tt > 1L) {
		ee <- .lde_pos_eesq(qs, hp$rho)$ee
		a <- (nu + R) / 2
		bb <- (nu + ee / hp$s2) / 2
		lam[, -1L] <- a / bb
		elog[, -1L] <- digamma(a) - log(bb)
		b[, -1L] <- bb
	}
	list(lam = lam, elog = elog, b = b)
}

# twice the KL divergence of q from the prior over both sides: the energy of q under the prior,
# minus the dimension and log |P0|, plus the log determinant of q's precision. scale-mixture
# weights enter the energy at E lambda and log |P0| at E log lambda, and add the KL of their
# gamma posterior from the Gamma(nu / 2, nu / 2) prior:
#   energy = sum_i [(1 - rho^2) E|u_i1|^2 + sum_t E lambda_it E|u_it - rho u_i,t-1|^2] / s2
#   log |P0| = R [Tt log(1 / s2) + log(1 - rho^2) + sum_t E log lambda_it]
.lde_pos_kl2 <- function(pos, R, Tt, hp, sym = FALSE, t_nu = 4) {
	one <- function(qs) {
		n <- dim(qs$mean)[1L]
		mix <- !is.null(qs$lam)
		ee <- .lde_pos_eesq(qs, hp$rho)
		w <- if (mix) qs$lam[, -1L, drop = FALSE] else 1
		energy <- ((1 - hp$rho^2) * sum(ee$e1) + if (Tt > 1L) sum(w * ee$ee) else 0) / hp$s2
		ldp0 <- n * R * (-Tt * log(hp$s2) + log(1 - hp$rho^2)) + if (mix) R * sum(qs$elog[, -1L]) else 0
		kl <- energy - n * R * Tt - ldp0 + qs$ldp
		if (mix && Tt > 1L) {
			a <- (t_nu + R) / 2
			b <- qs$b[, -1L]
			kl <- kl + 2 * sum(
				(a - t_nu / 2) * digamma(a) - lgamma(a) + lgamma(t_nu / 2) +
					(t_nu / 2) * (log(b) - log(t_nu / 2)) + a * (t_nu / 2 - b) / b
			)
		}
		kl
	}
	if (sym) one(pos$qu) else one(pos$qu) + one(pos$qv)
}

# minus two times the log hyperprior density of (s2, rho) at the held values
.lde_pos_prior_dev <- function(hp, pr) {
	2 * (pr$shape + 1) * log(hp$s2) + 2 * pr$scale / hp$s2 + (hp$rho - pr$rho_mean)^2 / pr$rho_sd^2
}

# working problem under the position noise

# quadrature rules of the cell problem: the static estimator's 9-node Gauss-Hermite rule for a
# nearly known predictor, and for a larger predictor variance the smallest rule of a ladder with
# at least `per_sd` nodes per sd, so the probit's unit-width transition is resolved.
# per_sd = 12 gives 9 nodes up to a variance of 0.56, 15 up to 1.6, 31 up to 6.7 and 61 from 12 on
.lde_gh_per_sd <- 12
.lde_gh_cache <- new.env(parent = emptyenv())
.lde_gh_rules <- function(K = 9L) {
	key <- as.character(K)
	hit <- .lde_gh_cache[[key]]
	if (!is.null(hit)) return(hit)
	rungs <- c(9L, 15L, 21L, 31L, 41L, 61L)
	lad <- c(K, rungs[rungs > K])
	rules <- lapply(lad, .ae_mf_gh)
	out <- list(K = lad, z = lapply(rules, function(r) r$z), w = lapply(rules, function(r) r$w))
	.lde_gh_cache[[key]] <- out
	out
}

# expected working problem of a binary (probit or logit) cell whose predictor is N(eta, vnoise),
# vnoise the variance of the product under q, as in .ae_mf_cells(): E log p and its derivative in
# eta, E psi, and the expected Fisher weight with its first two derivatives (for the Laplace
# correction), all by Gauss-Hermite quadrature on the ladder above. the working response is
# eta + E psi / E w (the increment capped as in .lde_work()); vnoise = 0 gives the plain problem.
#   dev:      -2 E log p, the deviance of the expected log likelihood
#   dev_marg: -2 log E p = -2 log Phi(s c eta), c = (1 + vnoise)^-1/2, the marginal deviance used by
#             the profile over the position scale, since E log p is optimised by the collapse
#             (by quadrature for the logit)
.lde_work_gh_ref <- function(
	y,
	eta,
	vnoise,
	K = 9L,
	delta_cap = formals(.lde_work)$delta_cap,
	wfloor = formals(.lde_work)$wfloor,
	block = 20000L,
	per_sd = .lde_gh_per_sd,
	link = "probit"
) {
	logit <- identical(link, "logit")
	ix <- which(is.finite(y))
	nobs <- length(ix)
	rules <- .lde_gh_rules(K)
	es <- ew <- eh <- ehh <- el <- dmarg <- numeric(nobs)
	vv <- pmax(vnoise[ix], 0)
	# the rule of every cell: the smallest rung with at least per_sd * sd nodes, the last one beyond it
	lev <- pmin(findInterval(per_sd * sqrt(vv), rules$K, left.open = TRUE) + 1L, length(rules$K))
	lev[is.na(lev)] <- 1L
	for (r in seq_along(rules$K)) {
		cells <- which(lev == r)
		nd <- list(z = rules$z[[r]], w = rules$w[[r]])
		# cells in blocks so the cells x nodes matrices stay small; the node terms are those of
		# .lde_probit_raw(), with the outcome picked by index
		for (b0 in seq.int(1L, max(length(cells), 1L), by = block)) {
			if (!length(cells)) break
			bi <- cells[b0:min(b0 + block - 1L, length(cells))]
			yy <- y[ix[bi]]
			m <- eta[ix[bi]]
			v <- vv[bi]
			e <- outer(sqrt(v), nd$z) + m
			one <- yy == 1
			if (logit) {
				# the score of an outcome of 1 is p0 and minus that of an outcome of 0 is p1 (the roles of lam1 and
				# lam0), so the weight is p0 p1 = lam1 lam0; w' = w (p0 - p1), w'' = w' (p0 - p1) - 2 w^2
				lp1 <- stats::plogis(e, log.p = TRUE)
				lp0 <- stats::plogis(e, lower.tail = FALSE, log.p = TRUE)
				lam1 <- stats::plogis(e, lower.tail = FALSE)
				lam0 <- stats::plogis(e)
				wr <- lam1 * lam0
				mm <- lam1 - lam0
				hpr <- wr * mm
				hppr <- hpr * mm - 2 * wr * wr
			} else {
				lp1 <- stats::pnorm(e, log.p = TRUE)
				lp0 <- stats::pnorm(e, lower.tail = FALSE, log.p = TRUE)
				ld <- -0.5 * e * e - 0.91893853320467267
				lam1 <- exp(ld - lp1)
				lam0 <- exp(ld - lp0)
				wr <- lam1 * lam0
				mm <- lam0 - lam1 - 2 * e
				hpr <- wr * mm
				hppr <- hpr * mm + wr * (lam0 * lam0 + lam1 * lam1 + e * (lam1 - lam0) - 2)
			}
			sc <- lam0
			sc[one, ] <- lam1[one, ]
			lg <- lp0
			lg[one, ] <- lp1[one, ]
			es[bi] <- as.vector(sc %*% nd$w) * (2 * one - 1)
			ew[bi] <- as.vector(pmax(wr, wfloor) %*% nd$w)
			eh[bi] <- as.vector(hpr %*% nd$w)
			ehh[bi] <- as.vector(hppr %*% nd$w)
			el[bi] <- as.vector(lg %*% nd$w)
			# marginal deviance: the noise integrated out of the likelihood, Phi(c eta) (by quadrature for the logit)
			if (logit) {
				pm <- lam1
				pm[one, ] <- lam0[one, ]
				dmarg[bi] <- log(as.vector(pm %*% nd$w))
			} else {
				dmarg[bi] <- .lde_probit_raw(yy, m / sqrt(1 + v))$ell
			}
		}
	}
	delta <- es / ew
	big <- abs(delta) > delta_cap
	if (any(big)) {
		delta[big] <- sign(delta[big]) * delta_cap
		ew[big] <- abs(es[big]) / delta_cap
	}
	z <- w <- score <- hprime <- hpp <- array(0, dim(y))
	z[ix] <- eta[ix] + delta
	w[ix] <- ew
	score[ix] <- es
	hprime[ix] <- ifelse(is.finite(eh), eh, 0)
	hpp[ix] <- ifelse(is.finite(ehh), ehh, 0)
	list(z = z, w = w, score = score, hprime = hprime, hpp = hpp, dev = -2 * sum(el), dev_marg = -2 * sum(dmarg))
}

# the same problem by the compiled kernel lde_gh_cells_cpp() (src/dynamic_eb_kernels.cpp), equal
# to the reference to rounding. `block` is ignored; `what` selects the outputs (0 the working
# response and weight, 1 also the score, the weight derivatives and the marginal deviance, 2 everything)
.lde_work_gh <- function(
	y,
	eta,
	vnoise,
	K = 9L,
	delta_cap = formals(.lde_work)$delta_cap,
	wfloor = formals(.lde_work)$wfloor,
	block = NULL,
	what = 2L,
	per_sd = .lde_gh_per_sd,
	link = "probit"
) {
	rules <- .lde_gh_rules(K)
	if (length(vnoise) != length(y)) vnoise <- rep_len(as.numeric(vnoise), length(y))
	lde_gh_cells_cpp(y, eta, vnoise, rules$z, rules$w, wfloor, delta_cap, what, per_sd, identical(link, "logit"))
}

# start

# start of q from a fit without positions: time-constant means from the leading singular triplets
# of each dyad's weighted mean working residual, shrunk by `scale` (the posterior mean of a bilinear
# term is smaller than its least-squares fit), and the prior covariance at S0 and rho0 (default
# the prior mean). S0 and scale only size the first iterate; the 1e-6 and 1e-12 below keep a
# square root and a division finite
.lde_pos_start <- function(ctx, fit0, S0 = 0.5, rho0 = NULL, scale = 0.3, jitter = 0) {
	if (is.null(rho0)) rho0 <- ctx$uvp$rho_mean
	na <- ctx$na; nb <- ctx$nb; R <- ctx$R; Tt <- ctx$Tt
	sv <- .lde_pos_svd(ctx, fit0, scale, jitter)
	U0 <- sv$U0; V0 <- sv$V0; sgn <- sv$sgn
	mk <- function(M0, n) {
		list(
			mean = array(rep(M0, Tt), dim = c(n, R, Tt)),
			sig = array(rep(diag(S0, R), n * Tt), dim = c(R, R, n, Tt)),
			cl = array(
				rep(diag(S0 * rho0, R), n * max(Tt - 1L, 1L)),
				dim = c(R, R, n, max(Tt - 1L, 1L))
			),
			lam = if (identical(ctx$uv_kind, "t")) matrix(1, n, Tt)
		)
	}
	qu <- mk(U0, na)
	qv <- if (ctx$sym) .lde_pos_mirror(qu, sgn) else mk(V0, nb)
	list(pos = list(qu = qu, qv = qv), hp = .lde_pos_hp(S0, rho0), sgn = sgn)
}

# the time-constant means of the start (the shrunk leading singular triplets of the weighted mean working
# residual of every dyad, or the leading eigenpairs of an undirected network, whose signs are the fixed sign
# pattern of the product F S F'), shared by the moving positions and the static factors
.lde_pos_svd <- function(ctx, fit0, scale = 0.3, jitter = 0, mm = FALSE) {
	ev <- fit0$ev
	na <- ctx$na; nb <- ctx$nb; R <- ctx$R; Tt <- ctx$Tt
	cs <- matrix(0, na, nb)
	rb <- matrix(0, na, nb)
	for (t in seq_len(Tt)) {
		w <- ev$wk$w[, , t]
		r <- (ev$wk$z[, , t] - ev$eta[, , t]) * ctx$obs[, , t]
		cs <- cs + w
		rb <- rb + w * r
	}
	# normal-panel weights are 1 / ve, so the floors are relative to the weights and singular values
	nrm <- isTRUE(ctx$nrm)
	M <- ifelse(cs > 0, rb / pmax(cs, 1e-12 * if (nrm) max(cs) else 1), 0)
	dfl <- 1e-6 * sqrt(ctx$prior$vunit)
	sgn <- NULL
	if (mm) {
		# weighted least-squares rank-R fit by the static estimator's MM iteration (.ae_lowrank_mm());
		# the one-shot spectrum of an undirected network can rank a noise eigenvalue above a signal one
		# and fix the wrong sign. weights are scaled by their maximum so a normal panel still moves
		lr <- .ae_lowrank_mm(matrix(0, na, nb), M, if (nrm) cs / max(cs) else cs, cs > 0, R, isTRUE(ctx$sym))
		if (ctx$sym) {
			sgn <- ifelse(lr$L < 0, -1, 1)
			U0 <- sweep(lr$U, 2L, sqrt(pmax(abs(lr$L), dfl)), "*") * scale
			V0 <- sweep(U0, 2L, sgn, "*")
		} else {
			U0 <- lr$U * scale
			V0 <- lr$V * scale
		}
	} else if (ctx$sym) {
		eg <- eigen((M + t(M)) / 2, symmetric = TRUE)
		idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
		sgn <- ifelse(eg$values[idx] < 0, -1, 1)
		dd <- sqrt(pmax(abs(eg$values[idx]), dfl))
		U0 <- sweep(eg$vectors[, idx, drop = FALSE], 2L, dd, "*") * scale
		V0 <- sweep(U0, 2L, sgn, "*")
	} else {
		sv <- svd(M, nu = R, nv = R)
		dd <- sqrt(pmax(sv$d[seq_len(R)], dfl))
		U0 <- sweep(sv$u[, seq_len(R), drop = FALSE], 2L, dd, "*") * scale
		V0 <- sweep(sv$v[, seq_len(R), drop = FALSE], 2L, dd, "*") * scale
	}
	if (jitter > 0) {
		U0 <- U0 + matrix(stats::rnorm(length(U0), 0, jitter), nrow(U0))
		V0 <- if (ctx$sym) sweep(U0, 2L, sgn, "*") else V0 + matrix(stats::rnorm(length(V0), 0, jitter), nrow(V0))
	}
	list(U0 = U0, V0 = V0, sgn = sgn)
}

# one evaluation

# evaluate a state at the stationary variance S in its hyperparameter record: one sweep of q at
# the incoming working problem, then the additive block at the new offset E_q[u'v] and noise
# var_q(u'v), its certificate and the next state. the certificate holds the additive block's
# scaled scores, the product's move in posterior sd and the hyperparameters' move under one EM
# step; the objective is the profile objective of S (see .lde_pos_profile)
.lde_eval_dyn <- function(ctx, st, tol = 5e-5, freeze_ab = FALSE, move_ve = TRUE) {
	q <- ctx$q; Tt <- ctx$Tt; R <- ctx$R
	sta <- st
	sta$O <- NULL
	eta_add <- .lde_eta(ctx, sta)
	static <- identical(ctx$fmode, "fac")
	hp <- st$hyp$uv
	ve_gap <- NULL
	if (isTRUE(ctx$vehyp) && move_ve) {
		# count or normal panel: the noise variance first moves to its EM image at the incoming iterate
		# (Gauss-Seidel, as in .lde_eval()); the size of the move is the residual of its EM equation
		vl <- if (is.null(st$vcell)) 0 else st$vcell
		ve1 <- .lde_ve_image_at(ctx, eta_add + st$O, st$Vpos + vl, st$cvpair, st$dshift, st$hyp)
		ve_gap <- abs(log(ve1 / st$hyp$ve))
		st$hyp$ve <- ve1
		sta$hyp$ve <- ve1
	}
	wk0 <- .lde_work_noise(ctx, eta_add + st$O, st$Vpos, st$hyp, 0L)
	if (static) {
		# static factors: the static estimator's sweeps, gauge and prior update at the period-summed problem
		fs <- .lde_fac_step(ctx, st, wk0, eta_add)
		pos <- fs$pos; O <- fs$O; Vpos <- fs$Vpos
	} else {
		pos <- .lde_pos_sweep(ctx, st$pos, hp, wk0, eta_add)
		if (!isTRUE(ctx$sym)) pos <- .lde_pos_gauge(pos, hp, R, Tt)
		# scale-mixture transition: the mixing weights take their variational update from the new q,
		# read by the next sweep
		if (identical(ctx$uv_kind, "t")) {
			lu <- .lde_pos_lambda(pos$qu, hp, ctx$t_nu, R)
			pos$qu[c("lam", "elog", "b")] <- lu[c("lam", "elog", "b")]
			if (isTRUE(ctx$sym)) {
				pos$qv <- .lde_pos_mirror(pos$qu, ctx$sgn)
			} else {
				lv <- .lde_pos_lambda(pos$qv, hp, ctx$t_nu, R)
				pos$qv[c("lam", "elog", "b")] <- lv[c("lam", "elog", "b")]
			}
		}
		O <- .lde_pos_product(pos, Tt)
		Vpos <- .lde_pos_prodvar(pos, R, Tt)
	}
	wk <- .lde_work_noise(ctx, eta_add + O, Vpos, st$hyp, 1L)
	wk_add <- wk
	wk_add$z <- wk$z - O
	sysm <- .lde_system(ctx, sta, wk_add)
	fac <- .lde_factor(sysm$A, sysm$b, sysm$C, sysm$Hzz, q, q, sysm$cdiag, unit = .lde_cunit(ctx))
	# with the additive hyperparameters held (the scale search) their posterior moments are not formed
	post <- .lde_posterior(fac, cond = TRUE, moments = !freeze_ab)
	lap <- .lde_laplace(ctx, fac, post$cond, wk, sta)
	em_self <- if (freeze_ab) NULL else .lde_em(ctx, sta, post)
	# the dyadic correlation takes its scoring step at the state reached, in the scale search too
	rstep <- .lde_pair_rho_step(ctx, st$hyp, lap, eta_add + O, Vpos)
	rho_new <- rstep$rho
	if (!is.null(rho_new) && !freeze_ab) em_self$rho_d <- rho_new
	# a state evaluated in place (a fit's record): the residual of its EM equation at this
	# evaluation's predictor variances
	dshift <- if (isTRUE(ctx$pln)) .lde_pln_dshift(ctx, sta, eta_add, lap, fac)
	if (isTRUE(ctx$vehyp) && is.null(ve_gap)) ve_gap <- abs(
		log(.lde_ve_image_at(ctx, eta_add + O, Vpos + lap$vcell, lap$cvpair, dshift, st$hyp) / st$hyp$ve)
	)
	if (!static) {
		rr <- .lde_pos_rho(.lde_pos_stats(pos, R, Tt, ctx$sym), hp$S, ctx$uvp, rho_ref = hp$rho)
		rho_em <- rr$rho
		hp_new <- .lde_pos_hp(hp$S, rho_em)
	} else if (identical(ctx$fscale, "cavity")) {
		# cavity rule: the factor prior scale moves toward its plug-in value from the cells' cavity
		# weights at the state reached, unless it is held
		fs$prior <- st$hyp$fac
		fs$change <- 0
		if (!isTRUE(ctx$hold_scale)) {
			nu_c <- .lde_cavity_w(ctx, eta_add + O, Vpos + lap$vcell, lap, st$hyp)$w
			fsu <- .lde_fac_scale_update(ctx, pos, st$hyp$fac, rowSums(nu_c, dims = 2L))
			fs$prior <- fsu$prior
			# the scale is held to a tenth of the other tolerances, as in the static estimator
			fs$change <- fsu$change * sqrt(tol) / (10 * tol)
		}
	}
	cert <- .lde_certificate(ctx, sta, wk, sysm, lap, em_self, tol)
	# move of the product in units of its posterior sd (sqrt(2) times a Newton step, see
	# .lde_certificate()); the 1e-6 keeps a collapsed cell from dividing by zero
	d_o <- sqrt(2) * max(abs(O - st$O)[ctx$obs] / sqrt(Vpos[ctx$obs] + 1e-6 * ctx$prior$vunit))
	gap_pos <- if (static) fs$change else abs(rho_em - hp$rho)
	# the dyadic correlation is held to a tenth of the other tolerances; with free additive
	# hyperparameters the certificate already counts its move
	if (freeze_ab && !is.null(rho_new)) gap_pos <- max(gap_pos, abs(rho_new - st$hyp$rho_d) * sqrt(tol) / (10 * tol))
	# the residual of a count cell's noise-variance EM equation counts with the positions' own gap
	if (!is.null(ve_gap)) gap_pos <- max(gap_pos, ve_gap)
	gap <- if (freeze_ab) gap_pos else max(cert$scale_change, gap_pos)
	# the mixing weights of a scale-mixture transition: their relative move
	if (!static && !is.null(pos$qu$lam) && Tt > 1L) {
		gap <- max(
			gap,
			abs(log(pos$qu$lam[, -1L]) - log(st$pos$qu$lam[, -1L])),
			if (!isTRUE(ctx$sym)) abs(log(pos$qv$lam[, -1L]) - log(st$pos$qv$lam[, -1L]))
		)
	}
	by_block <- c(cert$by_block, positions = d_o)
	thr <- sqrt(tol)
	mx <- max(by_block)
	cert$by_block <- by_block
	cert$max_score <- mx
	cert$scale_change <- gap
	cert$pass <- all(is.finite(c(by_block, gap))) && mx <= thr && gap <= thr
	cert$merit <- max(mx, gap)
	ste <- sta
	ste$pos <- pos; ste$O <- O; ste$Vpos <- Vpos
	# the persistence enters the objective only through the prior of q, so the objective is taken
	# at the EM persistence for this q
	wk_obj <- wk
	wk_obj$dev <- wk$dev_marg
	obj <- .lde_objective(ctx, sta, wk_obj, fac) + if (static) {
		fs$obj / ctx$nu
	} else {
		(.lde_pos_kl2(pos, R, Tt, hp, ctx$sym, ctx$t_nu) + .lde_pos_prior_dev(hp, ctx$uvp) - rr$gain) / ctx$nu
	}
	ev <- list(
		st = ste,
		eta = eta_add + O,
		wk = wk,
		sysm = sysm,
		fac = fac,
		post = post,
		lap = lap,
		em_self = em_self,
		cert = cert,
		objective = obj,
		pos_in = st$pos,
		rho_step = rstep
	)
	if (isTRUE(ctx$vehyp)) ev$ve_gap <- ve_gap
	rz <- c(numeric(q), sysm$rth - lap$corr)
	facs <- .lde_step_factor(ctx, sta, fac, lap)
	rz <- rz + as.vector(facs$dH %*% c(sta$lev, sta$theta))
	sol <- .lde_solve(facs, sysm$rs, rz)
	cst <- .lde_candidate(ctx, sta, sol, NULL)
	cst$hyp <- if (freeze_ab) st$hyp else .lde_em(ctx, cst, post)
	if (!is.null(rho_new)) cst$hyp$rho_d <- rho_new
	if (static) cst$hyp$fac <- fs$prior else cst$hyp$uv <- hp_new
	cst$pos <- pos; cst$O <- O; cst$Vpos <- Vpos
	if (isTRUE(ctx$vehyp)) { cst$vcell <- lap$vcell; cst$cvpair <- lap$cvpair; cst$dshift <- dshift }
	ev$cand <- cst
	ev
}

# the positions' q mixed linearly, (1 - a) q0 + a q1, in the means and the diagonal and lag-one
# covariance blocks; the KL pieces stay those of q1 until the next sweep, and an undirected
# network's receiver side mirrors the sender side
.lde_pos_mix <- function(pos0, pos1, a, sym = FALSE, sgn = NULL) {
	mix <- function(q0, q1) {
		out <- q1
		for (f in c(
			"mean",
			"sig",
			"cl",
			if (!is.null(q0$lam) && !is.null(q1$lam)) "lam"
		)) out[[f]] <- (1 - a) * q0[[f]] + a * q1[[f]]
		out
	}
	qu <- mix(pos0$qu, pos1$qu)
	list(qu = qu, qv = if (sym) .lde_pos_mirror(qu, sgn) else mix(pos0$qv, pos1$qv))
}

# blend two states (a damped step): the additive part and the hyperparameters mix, and the
# positions mix with `pos0`, the incoming q of the evaluation that produced s1
.lde_blend_dyn <- function(ctx, s0, s1, a, pos0 = NULL) {
	out <- .lde_blend(ctx, s0, s1, a)
	out$pos <- s1$pos; out$O <- s1$O; out$Vpos <- s1$Vpos
	static <- identical(ctx$fmode, "fac")
	if (static) out$hyp$fac <- .lde_fac_mix(s0$hyp$fac, s1$hyp$fac, a)
	if (!is.null(pos0) && a < 1) {
		out$pos <- .lde_pos_mix(pos0, s1$pos, a, isTRUE(ctx$sym), ctx$sgn)
		if (static) {
			pp <- .lde_fac_product(ctx, .lde_fac_state(ctx, out$pos, out$hyp$fac))
			out$O <- pp$O
			out$Vpos <- pp$Vpos
		} else {
			out$O <- .lde_pos_product(out$pos, ctx$Tt)
			out$Vpos <- .lde_pos_prodvar(out$pos, ctx$R, ctx$Tt)
		}
	}
	if (!static) out$hyp$uv <- .lde_pos_hp(s1$hyp$uv$S, (1 - a) * s0$hyp$uv$rho + a * s1$hyp$uv$rho)
	out
}

# the iteration at a held position scale

# the numeric leaves of a (nested) list as one vector, and the list with such a vector put back
.lde_leaves <- function(x) {
	if (is.list(x)) unlist(lapply(x, .lde_leaves), use.names = FALSE)
	else if (is.numeric(x)) as.numeric(x)
	else numeric(0)
}
.lde_leaves_put <- function(x, v, pos) {
	if (is.list(x)) {
		for (k in seq_along(x)) x[k] <- list(.lde_leaves_put(x[[k]], v, pos))
		x
	} else if (is.numeric(x)) {
		n <- length(x)
		x[] <- v[pos$i + seq_len(n)]
		pos$i <- pos$i + n
		x
	} else x
}

# the continuous part of a state as one vector: coefficients, levels, paths, hyperparameters and
# both sides' position means and covariances (the extrapolation moves the covariances with the
# means, or the next sweep pulls the tail back)
.lde_dyn_vec <- function(st) {
	c(
		st$theta,
		st$lev,
		as.vector(st$x),
		.lde_leaves(st$hyp),
		unlist(lapply(c("qu", "qv"), function(side) {
			lapply(
				c("mean", "sig", "cl", if (!is.null(st$pos[[side]]$lam)) "lam"),
				function(f) as.vector(st$pos[[side]][[f]])
			)
		}), use.names = FALSE)
	)
}

# the state with such a vector put back; the position hyperparameters stay consistent
# (s2 = S (1 - rho^2)), the product and its variance follow the positions (the receiver side of an
# undirected network is the sign-mapped mirror of the sender side)
.lde_dyn_put <- function(ctx, st, v) {
	i <- 0L
	take <- function(n) {
		out <- v[i + seq_len(n)]
		i <<- i + n
		out
	}
	st$theta <- take(length(st$theta))
	st$lev <- take(length(st$lev))
	st$x[] <- take(length(st$x))
	pos <- new.env()
	pos$i <- i
	st$hyp <- .lde_leaves_put(st$hyp, v, pos)
	i <- pos$i
	if (!is.null(st$hyp$uv)) st$hyp$uv <- .lde_pos_hp(st$hyp$uv$S, st$hyp$uv$rho)
	if (!is.null(st$hyp$fac)) st$hyp$fac <- .lde_fac_fix(st$hyp$fac)
	# an extrapolated correlation beyond the open interval is not a state
	if (!is.null(st$hyp$rho_d) && !(abs(st$hyp$rho_d) < 0.98)) stop("invalid dyadic correlation")
	# nor is a variance of the noise that is not positive
	if (!is.null(st$hyp$ve) && !(st$hyp$ve > 0)) stop("invalid noise variance")
	for (side in c("qu", "qv")) {
		for (f in c(
			"mean",
			"sig",
			"cl",
			if (!is.null(st$pos[[side]]$lam)) "lam"
		)) st$pos[[side]][[f]][] <- take(length(st$pos[[side]][[f]]))
		# a mixing weight is positive (an extrapolation past zero is not a state)
		if (!is.null(st$pos[[side]]$lam) && any(!(st$pos[[side]]$lam > 0))) stop("invalid transition weights")
	}
	if (isTRUE(ctx$sym)) st$pos$qv <- .lde_pos_mirror(st$pos$qu, ctx$sgn)
	if (identical(ctx$fmode, "fac")) {
		pp <- .lde_fac_product(ctx, .lde_fac_state(ctx, st$pos, st$hyp$fac))
		st$O <- pp$O
		st$Vpos <- pp$Vpos
	} else {
		st$O <- .lde_pos_product(st$pos, ctx$Tt)
		st$Vpos <- .lde_pos_prodvar(st$pos, ctx$R, ctx$Tt)
	}
	st
}

# squared extrapolation of a fixed-point iteration (SQUAREM, Varadhan and Roland 2008) from x0,
# x1 = F(x0), x2 = F(x1): with r = x1 - x0 and v = x2 - 2 x1 + x0 the step length is
# alpha = -|r| / |v|, at most -1, and the new point is x0 - 2 alpha r + alpha^2 v
.lde_squarem <- function(x0, x1, x2, shrink = 0) {
	r <- x1 - x0
	v <- x2 - 2 * x1 + x0
	nv <- sqrt(sum(v^2))
	if (!is.finite(nv) || nv <= 0) return(NULL)
	alpha <- min(-sqrt(sum(r^2)) / nv, -1)
	if (shrink > 0) alpha <- -1 + (alpha + 1) * (1 - shrink)
	x0 - 2 * alpha * r + alpha^2 * v
}

# objective decrease still to come: the last change times r / (1 - r), r the larger of the last
# two ratios of successive changes (capped at 0.999); zero when the last change is zero
.lde_dyn_remaining <- function(obj) {
	if (length(obj) < 4L) return(Inf)
	d <- abs(diff(utils::tail(obj, 4L)))
	if (d[3L] == 0) return(0)
	ratio <- function(a, b) if (a > 0) b / a else if (b > 0) Inf else 0
	r <- min(max(ratio(d[1L], d[2L]), ratio(d[2L], d[3L])), 0.999)
	d[3L] * r / (1 - r)
}

# the reweighting loop at a held position scale (the loop of .lde_fit_core() with the evaluation above)
#   settle = TRUE: converged when the certificate passes, the objective has stayed within tol / 10
#   over the last `window` passes and an extrapolation moved it by less than 10 times that range;
#   the certificate alone misses the slow modes of a weakly identified bilinear term. the range is
#   relaxed for a search point far above the best fixed point so far (obj_ref, `excess`)
#   relaxation: omega damps the step of the whole state, halved when successive moves reverse and
#   raised when they align, which breaks the 2-cycles of strongly identified panels
#   acceleration: after `gate` plain passes that have not halved the merit, a squared extrapolation
#   (.lde_squarem()) is tried; it is kept when finite, within `kappa` of the plain iterate's merit
#   and without a jump in the objective, and it must lower the merit a cycle later, else the gate
#   doubles and three failures in a row end extrapolation. it changes the speed, not the fixed point
#   settle = FALSE: the certificate alone, the starting state included
.lde_fit_core_dyn <- function(
	ctx,
	st0,
	max_iter = 200L,
	tol = 5e-5,
	freeze_ab = FALSE,
	settle = TRUE,
	window = 10L,
	gate = 5L,
	kappa = 100,
	obj_ref = NULL,
	excess = 0.01
) {
	st <- st0
	ev <- .lde_eval_dyn(ctx, st, tol, freeze_ab = freeze_ab, move_ve = !(isTRUE(st$seeded) && !settle))
	trace <- data.frame(
		iter = 0L,
		merit = ev$cert$merit,
		score = ev$cert$max_score,
		gap = ev$cert$scale_change,
		objective = ev$objective,
		step = NA_real_,
		jump = FALSE
	)
	converged <- !settle && isTRUE(ev$cert$pass)
	pln_t <- isTRUE(ctx$vehyp)
	thr_t <- sqrt(tol)
	if (pln_t && converged) {
		# a count fit is also certified by its noise variance's distance to the limit of its iteration
		# (.lde_tail()), from the EM residual and the tail rate the record carries
		r0 <- if (isTRUE(is.finite(st$ve_rate))) min(st$ve_rate, 0.9999) else NA_real_
		converged <- isTRUE(is.finite(r0)) && isTRUE(ev$ve_gap * r0 / (1 - r0) <= thr_t / 2)
	}
	termination <- if (converged) "converged" else "iteration limit"
	it <- 0L
	omega <- 1
	# the scalar hyperparameters at the accepted plain passes since the last extrapolation (a count fit)
	hs <- if (pln_t) matrix(.lde_tail_scalars(ctx, st), nrow = 1L)
	v_old <- .lde_hvec(ctx, ev$cand)
	move_old <- NULL
	best <- NULL
	hist_obj <- ev$objective
	jumped <- FALSE
	rem_jump <- Inf
	# the first extrapolation also waits `gate` plain passes
	since_jump <- 0L
	gate_eff <- gate
	kappa_eff <- kappa
	jump_ref <- NA_real_
	n_fail <- 0L
	while (!converged && it < max_iter) {
		it <- it + 1L
		g <- ev$cand
		cand <- if (omega < 1) .lde_blend_dyn(ctx, ev$st, g, omega, pos0 = ev$pos_in) else g
		jumped <- FALSE
		rem_jump <- Inf
		evc <- tryCatch(.lde_eval_dyn(ctx, cand, tol, freeze_ab = freeze_ab), error = function(e) NULL)
		if (is.null(evc) || !is.finite(evc$cert$merit) || !is.finite(evc$objective)) {
			termination <- "no improving step"
			break
		}
		# squared extrapolation from x0, x1 = F(x0) (cand) and x2 = F(x1) when the plain pass has not
		# halved the merit. first check that the last extrapolation paid off; if not, later ones are
		# rarer and must improve the merit, and three failures in a row stop them
		if (!is.na(jump_ref) && since_jump >= gate_eff) {
			if (ev$cert$merit < jump_ref) {
				n_fail <- 0L
				gate_eff <- gate
				kappa_eff <- kappa
			} else {
				n_fail <- n_fail + 1L
				gate_eff <- if (n_fail >= 3L) Inf else min(2 * gate_eff, 64L)
				kappa_eff <- 1
			}
			jump_ref <- NA_real_
		}
		if (settle && it < max_iter && evc$cert$merit > 0.5 * ev$cert$merit) {
			x0 <- .lde_dyn_vec(st)
			x1 <- .lde_dyn_vec(cand)
			# the second iterate of the map being iterated: the plain one, or the damped one
			x2 <- .lde_dyn_vec(
				if (omega < 1) .lde_blend_dyn(ctx, evc$st, evc$cand, omega, pos0 = evc$pos_in) else evc$cand
			)
			# two moves that go back and forth (cosine below -1/2) are no tail to extrapolate; the
			# relaxation takes over
			d1 <- x1 - x0
			d2 <- x2 - x1
			n1 <- sqrt(sum(d1^2))
			n2 <- sqrt(sum(d2^2))
			nn <- n1 * n2
			osc <- is.finite(nn) && nn > 0 && sum(d1 * d2) / nn < -0.5
			if (osc) omega <- max(omega / 2, 1 / 64)
			# a tail is a contraction: when the second move is not the shorter, the iteration is still
			# growing into its state and an extrapolation would run off
			grow <- !(is.finite(n1) && is.finite(n2) && n2 < n1)
			# the gate holds back the extrapolation, not the detection of an oscillating map
			for (shrink in if (osc || grow || since_jump < gate_eff) numeric(0) else c(0, 0.5)) {
				xs <- .lde_squarem(x0, x1, x2, shrink)
				if (is.null(xs) || any(!is.finite(xs))) break
				# an extrapolated state with an invalid hyperparameter fails or warns when evaluated and is rejected
				cx <- tryCatch(.lde_dyn_put(ctx, cand, xs), error = function(e) NULL, warning = function(w) NULL)
				ex <- if (is.null(cx)) NULL else tryCatch(
					.lde_eval_dyn(ctx, cx, tol, freeze_ab = freeze_ab),
					error = function(e) NULL,
					warning = function(w) NULL
				)
				# a state whose objective is far above the plain iterate's is rejected: the room is 10 deviance
				# units plus the range the objective has covered in this call
				room <- 10 + diff(range(c(hist_obj, evc$objective)))
				if (!is.null(ex) && is.finite(ex$objective) && is.finite(ex$cert$merit) && ex$cert$merit < kappa_eff * evc$cert$merit &&
						ex$objective <= evc$objective + room) {
					# objective change from the plain iterate to its extrapolation: what the extrapolation claims is
					# still to come
					rem_jump <- if (shrink == 0) abs(ex$objective - evc$objective) else Inf
					jump_ref <- evc$cert$merit
					cand <- cx; evc <- ex; jumped <- TRUE
					break
				}
			}
		}
		since_jump <- if (jumped) 0L else since_jump + 1L
		st <- cand
		ev <- evc
		if (pln_t) hs <- if (jumped) matrix(
			.lde_tail_scalars(ctx, ev$st),
			nrow = 1L
		) else rbind(hs, .lde_tail_scalars(ctx, ev$st))
		if (is.null(best) || ev$cert$merit < best$ev$cert$merit) best <- list(st = st, ev = ev)
		hist_obj <- c(hist_obj, ev$objective)
		if (jumped) {
			v_old <- .lde_hvec(ctx, ev$cand)
			move_old <- NULL
		} else {
			v_new <- .lde_hvec(ctx, ev$cand)
			move <- v_new - v_old
			v_old <- v_new
			if (!is.null(move_old)) {
				nn <- sqrt(sum(move^2)) * sqrt(sum(move_old^2))
				if (is.finite(nn) && nn > 0) {
					cs <- sum(move * move_old) / nn
					if (cs < -0.5) omega <- max(omega / 2, 1 / 64)
					else if (cs > 0.5) omega <- min(1, omega * 1.5)
				}
			}
			move_old <- move
		}
		trace <- rbind(
			trace,
			data.frame(
				iter = it,
				merit = ev$cert$merit,
				score = ev$cert$max_score,
				gap = ev$cert$scale_change,
				objective = ev$objective,
				step = omega,
				jump = jumped
			)
		)
		n_h <- length(hist_obj)
		# how far the objective may still move: tol / 10, or a fraction `excess` of a search point's
		# distance above the best fixed point so far (obj_ref)
		slack <- tol / 10
		if (!is.null(obj_ref) && is.finite(obj_ref)) slack <- max(slack, excess * max(0, ev$objective - obj_ref))
		quiet <- n_h > window && diff(range(hist_obj[(n_h - window):n_h])) <= slack
		if (!pln_t) {
			if (isTRUE(ev$cert$pass) && (!settle || (quiet && (!jumped || rem_jump <= 10 * slack)))) {
				converged <- TRUE
				termination <- "converged"
			}
		} else if (isTRUE(ev$cert$pass) && (!settle || .lde_tail_ok(hs, thr_t))) {
			# a count fit: the distance of the noise variance (and static factor scales) to the limit of the
			# iteration, read from the tails of their moves, replaces the quiet-objective test for static
			# factors, which cannot see the end of a slow tail; the moving positions keep it
			if (!settle || identical(ctx$fmode, "fac") || (quiet && (!jumped || rem_jump <= 10 * slack))) {
				converged <- TRUE
				termination <- "converged"
			}
		}
	}
	if (pln_t) {
		ve_rate <- if (nrow(hs) >= 5L) .lde_tail(hs[, 1L])$rate else NA_real_
		st$ve_rate <- ve_rate
		ev$st$ve_rate <- ve_rate
	}
	if (!converged && !is.null(best) && best$ev$cert$merit < ev$cert$merit) { st <- best$st; ev <- best$ev }
	# objective at the fixed point: extrapolated from a clean tail (the last three changes of one sign
	# and shrinking) when the iteration stopped short, else the objective of the last state
	rem_end <- if (converged) 0 else .lde_dyn_remaining(hist_obj)
	d_last <- if (length(hist_obj) >= 2L) utils::tail(diff(hist_obj), 1L) else 0
	d3 <- if (length(hist_obj) >= 4L) diff(utils::tail(hist_obj, 4L)) else NULL
	clean <- !is.null(d3) && all(d3 != 0) && all(sign(d3) == sign(d3[3L])) &&
		abs(d3[2L]) < abs(d3[1L]) && abs(d3[3L]) < abs(d3[2L])
	list(
		st = ev$st,
		ev = ev,
		converged = converged,
		iterations = it,
		termination = termination,
		trace = trace,
		objective_limit = ev$objective + if (clean && is.finite(rem_end)) sign(d_last) * rem_end else 0
	)
}

# scale of the positions

# golden-section search of a function of S over [lower, upper] on the log scale, resolving S to the
# factor exp(width * phi^n_iter). a best evaluation within that resolution of an outer edge extends
# the range by its own width on that side, at most `max_ext` times. returns every evaluation and
# whether the best one is still at an edge
.lde_golden_log <- function(at, lower, upper, n_iter = 7L, max_ext = 2L) {
	gr <- (sqrt(5) - 1) / 2
	S <- numeric(0)
	v <- numeric(0)
	ev <- function(x) {
		val <- at(exp(x))
		S <<- c(S, exp(x))
		v <<- c(v, val)
		val
	}
	search <- function(a, b) {
		c1 <- b - gr * (b - a); d1 <- a + gr * (b - a)
		fc <- ev(c1); fd <- ev(d1)
		for (k in seq_len(n_iter)) {
			if (fc < fd) {
				b <- d1; d1 <- c1; fd <- fc; c1 <- b - gr * (b - a); fc <- ev(c1)
			} else {
				a <- c1; c1 <- d1; fc <- fd; d1 <- a + gr * (b - a); fd <- ev(d1)
			}
		}
	}
	a_all <- log(lower); b_all <- log(upper)
	width <- b_all - a_all
	tol_edge <- width * gr^n_iter
	search(a_all, b_all)
	for (ext in seq_len(max_ext)) {
		lb <- log(S[which.min(v)])
		if (lb > b_all - tol_edge) {
			search(b_all, b_all + width)
			b_all <- b_all + width
		} else if (lb < a_all + tol_edge) {
			search(a_all - width, a_all)
			a_all <- a_all - width
		} else break
	}
	lb <- log(S[which.min(v)])
	list(
		S = S,
		value = v,
		edge = lb > b_all - tol_edge || lb < a_all + tol_edge,
		range = exp(c(a_all, b_all))
	)
}

# quadratic search of a function of S on the log scale about a centre: the centre and the points w
# either side, then, while the parabola through the three evaluations about the best one has its
# vertex outside them (or is not convex), the vertex itself, at most max_extra times and within
# `span` of the centre. `edge` is TRUE when no vertex was found inside three evaluations
.lde_quad_search_log <- function(at, centre, w = 0.12, max_extra = 3L, span = 0.7) {
	x0 <- log(centre)
	S <- numeric(0)
	v <- numeric(0)
	ev <- function(x) {
		val <- at(exp(x))
		S <<- c(S, exp(x))
		v <<- c(v, val)
		val
	}
	for (x in x0 + c(0, -w, w)) ev(x)
	edge <- TRUE
	for (k in 0:max_extra) {
		ls <- log(S)
		fin <- is.finite(v)
		if (!any(fin)) break
		ib <- which(fin)[which.min(v[fin])]
		o <- order(ls)
		pos <- which(o == ib)
		lo <- max(1L, min(pos - 1L, length(o) - 2L))
		idx <- o[lo:(lo + 2L)]
		if (!all(fin[idx]) || length(idx) < 3L) break
		xs <- ls[idx]
		fs <- v[idx]
		# the parabola through three points by divided differences
		d1 <- (fs[2L] - fs[1L]) / (xs[2L] - xs[1L])
		d2 <- ((fs[3L] - fs[2L]) / (xs[3L] - xs[2L]) - d1) / (xs[3L] - xs[1L])
		if (is.finite(d2) && d2 > 0) {
			vtx <- (xs[1L] + xs[2L]) / 2 - d1 / (2 * d2)
			# a vertex inside the three evaluations, or within 1 percent of S of the outer one, counts as found
			if (vtx >= min(xs) - 0.01 && vtx <= max(xs) + 0.01) {
				edge <- FALSE
				break
			}
			xnew <- vtx
		} else {
			# not convex: a step of 2 w beyond the best evaluation, away from the worse neighbour
			xnew <- ls[ib] + (if (ib == idx[1L]) -2 * w else if (ib == idx[3L]) 2 * w else 0)
			if (xnew == ls[ib]) break
		}
		xnew <- min(max(xnew, x0 - span), x0 + span)
		if (k == max_extra || any(abs(xnew - ls) < 1e-6)) break
		ev(xnew)
	}
	list(S = S, value = v, edge = edge, range = exp(range(log(S))))
}

# unit of the stationary variance of the positions: one where the link fixes the latent scale, the
# sd of the response for a normal panel (the product of two coordinates is added to the response)
.lde_pos_scale <- function(ctx) if (isTRUE(ctx$nrm)) sqrt(ctx$prior$ve_prior[["s20"]]) else 1

# profile over the stationary variance S of the positions: at each S the mean-field iteration
# (persistence by EM, additive variance components held so the objective is smooth in S) runs from
# the solution at the nearest uncollapsed scale, and its value is the profile objective. the
# evaluations nearest the best that hit the iteration cap are continued to convergence, a quadratic
# in log S through them refines the minimiser, the refined scale is kept only if no worse than the
# best evaluation, and the fit there is polished with every hyperparameter free. returns the
# polished fit and the search record.
#   * the starting bracket [0.08, 4] runs from a product that barely moves a probit predictor to
#     one that dominates the unit noise; a minimiser at an edge extends it (.lde_golden_log)
#   * n_iter = 7 golden steps narrow it to 3.4% of its log-width; fewer can end on the flat second
#     branch of fixed points of a weakly identified panel
#   * inner fits settle to tol_inner = 1e-3 deviance units (a point above the best to a fraction
#     `excess` of its distance), with max_inner iterations at first and the whole budget for the
#     deciding evaluations; one that still does not settle is reported as `unsettled`
#   * the quadratic uses the five evaluations nearest the best
#   * `quick` (fits with the dyadic correlation, .lde_pos_profile_pair()) uses .lde_quad_search_log()
#     about a centre and the three nearest evaluations, and returns `failed` before any polish when
#     that search ends without a vertex or with an unsettled evaluation
.lde_pos_profile <- function(
	ctx,
	st0,
	tol = 5e-5,
	max_iter = 300L,
	verbose = FALSE,
	lower = 0.08 * .lde_pos_scale(ctx),
	upper = 4 * .lde_pos_scale(ctx),
	n_iter = 7L,
	tol_inner = 1e-3,
	max_inner = 100L,
	quick = NULL
) {
	cache <- list()
	key <- function(S) sprintf("%.8f", S)
	# with the dyadic correlation the inner fits of the search use the coarse tier of the pair
	# quadrature (.lde_pair_tier1); the polish at the scale found uses the full accuracy
	ctx_s <- ctx
	if (isTRUE(ctx$dcor)) ctx_s$pair_tier <- 1L
	# a state whose position means are numerically zero (a unit of a position is of order one)
	collapsed <- function(st) max(abs(st$pos$qu$mean)) < 1e-6
	at <- function(S) {
		# start from the solution at the nearest uncollapsed scale evaluated so far (in log S), so a short
		# move of the scale is a short transient
		cur <- st0
		ok <- vapply(cache, function(z) !collapsed(z$st), logical(1))
		if (any(ok)) {
			cc <- cache[ok]
			cur <- cc[[which.min(abs(log(vapply(cc, function(z) z$S, numeric(1))) - log(S)))]]$st
		}
		cur$hyp$uv <- .lde_pos_hp(S, cur$hyp$uv$rho)
		vals <- vapply(cache, function(z) z$value, numeric(1))
		fit <- .lde_fit_core_dyn(
			ctx_s,
			cur,
			max_iter = min(max_inner, max_iter),
			tol = tol_inner,
			freeze_ab = TRUE,
			obj_ref = if (any(is.finite(vals))) min(vals[is.finite(vals)])
		)
		# the objective at the fixed point of this S (extrapolated from a clean tail if the cap was reached)
		val <- fit$objective_limit
		if (!is.finite(val)) val <- Inf
		cache[[key(S)]] <<- list(
			S = S,
			st = fit$st,
			converged = fit$converged,
			iterations = fit$iterations,
			value = val
		)
		if (verbose) {
			cli::cli_text("Position scale S = {format(S, digits = 3)}: objective {format(val, digits = 7)}, persistence {format(fit$st$hyp$uv$rho, digits = 3)}, {fit$iterations} passes{if (!fit$converged) ' (cap)' else ''}")
		}
		val
	}
	gs <- if (is.null(quick)) {
		.lde_golden_log(at, lower, upper, n_iter = n_iter)
	} else {
		.lde_quad_search_log(at, quick$centre, w = quick$w, max_extra = quick$max_extra, span = quick$span)
	}
	ss <- gs$S
	values <- function() vapply(ss, function(S) cache[[key(S)]]$value, numeric(1))
	# the quadratic search leaves a parabola through three evaluations, which the refinement reads
	kq <- if (is.null(quick)) 5L else 3L
	nearest <- function(vals, k = kq) {
		ib <- which.min(vals)
		list(ib = ib, near = order(abs(log(ss) - log(ss[ib])))[seq_len(min(k, length(ss)))])
	}
	# the evaluations that decide the minimiser must be fixed points: those among the five nearest the best
	# that stopped at the cap are continued, up to the whole budget, until one fails to settle
	unsettled <- FALSE
	for (round in 1:3) {
		nr <- nearest(values())
		todo <- nr$near[!vapply(nr$near, function(k) isTRUE(cache[[key(ss[k])]]$converged), logical(1))]
		if (!length(todo) || unsettled) break
		for (k in todo) {
			ent <- cache[[key(ss[k])]]
			f2 <- .lde_fit_core_dyn(ctx_s, ent$st, max_iter = max_iter, tol = tol_inner, freeze_ab = TRUE)
			val <- f2$objective_limit
			if (!is.finite(val)) val <- Inf
			cache[[key(ss[k])]] <- list(
				S = ent$S,
				st = f2$st,
				converged = f2$converged,
				iterations = ent$iterations + f2$iterations,
				value = val
			)
			if (!f2$converged) { unsettled <- TRUE; break }
		}
	}
	vals <- values()
	nr <- nearest(vals)
	unsettled <- unsettled || any(!vapply(nr$near, function(k) isTRUE(cache[[key(ss[k])]]$converged), logical(1)))
	# the deciding evaluations and how many did not settle, for the refusal message (an evaluation seen
	# to fail counts even after it leaves the nearest set)
	n_near <- length(nr$near)
	n_unsettled <- max(
		sum(!vapply(nr$near, function(k) isTRUE(cache[[key(ss[k])]]$converged), logical(1))),
		as.integer(unsettled)
	)
	# a quadratic search with no vertex inside three evaluations, or with one unsettled, hands the
	# panel to the golden search before any polish
	if (!is.null(quick) && (isTRUE(gs$edge) || unsettled)) {
		return(
			list(
				fit = NULL,
				S = NA_real_,
				grid = NULL,
				range = gs$range,
				at_edge = gs$edge,
				unsettled = unsettled,
				failed = TRUE
			)
		)
	}
	ib <- nr$ib
	s_star <- ss[ib]
	# a quadratic in log S through the evaluations nearest the best one
	near <- nr$near
	if (length(near) >= (if (is.null(quick)) 4L else 3L) && all(is.finite(vals[near]))) {
		ls <- log(ss[near])
		cf <- tryCatch(stats::coef(stats::lm(vals[near] ~ ls + I(ls^2))), error = function(e) NULL)
		if (!is.null(cf) && all(is.finite(cf)) && cf[3L] > 0) {
			v <- unname(-cf[2L] / (2 * cf[3L]))
			if (v > min(ls) && v < max(ls)) s_star <- exp(v)
		}
	}
	# the refined scale is an evaluation of its own, kept only if it is no worse than the best one
	if (s_star != ss[ib]) {
		v_star <- at(s_star)
		if (!(v_star <= vals[ib])) s_star <- ss[ib]
	}
	cl <- cache[[which.min(abs(log(vapply(cache, function(z) z$S, numeric(1))) - log(s_star)))]]
	start <- cl$st
	start$hyp$uv <- .lde_pos_hp(s_star, start$hyp$uv$rho)
	fit <- .lde_fit_core_dyn(ctx, start, max_iter = max_iter, tol = tol)
	grid <- data.frame(
		S = vapply(cache, function(z) z$S, numeric(1)),
		objective = vapply(cache, function(z) z$value, numeric(1)),
		iterations = vapply(cache, function(z) z$iterations, numeric(1)),
		converged = vapply(cache, function(z) isTRUE(z$converged), logical(1))
	)
	rownames(grid) <- NULL
	list(
		fit = fit,
		S = s_star,
		grid = grid[order(grid$S), ],
		range = gs$range,
		at_edge = gs$edge,
		unsettled = unsettled,
		n_near = n_near,
		n_unsettled = n_unsettled
	)
}

# profile over the scale for a fit with the dyadic correlation: the cheap independent-cell fit
# first, then the pair likelihood by the quadratic search about the scale it found
# (.lde_quad_search_log(), 4 to 6 inner fits). a panel where that search finds no vertex inside
# three evaluations or leaves one unsettled runs the search over [0.08, 4] with the pair
# likelihood; options(lame.dcor_full_search = TRUE) runs that search throughout (see
# lame_dynamic_als_eb.R) and options(lame.dcor_golden_search = TRUE) a golden search over
# S_i exp(+-width)
.lde_pos_profile_pair <- function(ctx, st0, tol = 5e-5, max_iter = 300L, verbose = FALSE, width = 0.37, n_iter = 4L) {
	ctx_i <- ctx
	ctx_i$dcor <- FALSE
	ctx_i$pc <- NULL
	st_i <- st0
	st_i$hyp$rho_d <- NULL
	pi <- .lde_pos_profile(ctx_i, st_i, tol = tol, max_iter = max_iter, verbose = verbose)
	s_i <- pi$S
	# every inner fit starts from the independent solution (positions, paths, levels), with the
	# additive hyperparameters held at the start's values as in the full search and the correlation at zero
	st1 <- pi$fit$st
	keep <- setdiff(names(st1$hyp), c("uv", "rho_d"))
	st1$hyp[keep] <- st0$hyp[keep]
	st1$hyp$rho_d <- 0
	# quadratic search about the independent cells' scale; with no vertex, an unsettled evaluation or
	# an unsettled independent search, the full search over [0.08, 4] runs from the start of the fit
	# (or the golden search over S_i exp(+-width) under options(lame.dcor_golden_search = TRUE))
	golden <- isTRUE(getOption("lame.dcor_golden_search", FALSE))
	pr <- NULL
	if (!golden && !isTRUE(pi$at_edge) && !isTRUE(pi$unsettled)) {
		pr <- .lde_pos_profile(
			ctx,
			st1,
			tol = tol,
			max_iter = max_iter,
			verbose = verbose,
			quick = list(centre = s_i, w = 0.12, max_extra = 3L, span = 0.7)
		)
		if (isTRUE(pr$failed)) pr <- NULL
	}
	if (is.null(pr)) {
		# the pair quadrature's plans depend on the calls since the cache was emptied; reset it so the
		# search runs as it would alone
		.lde_pair_reset()
		pr <- if (golden) {
			.lde_pos_profile(
				ctx,
				st1,
				tol = tol,
				max_iter = max_iter,
				verbose = verbose,
				lower = s_i * exp(-width),
				upper = s_i * exp(width),
				n_iter = n_iter
			)
		} else {
			.lde_pos_profile(ctx, st0, tol = tol, max_iter = max_iter, verbose = verbose)
		}
	}
	pr$independent <- list(S = s_i, grid = pi$grid, iterations = pi$fit$iterations)
	pr
}

# reported positions

# position means and posterior sd in the gauge of the time-averaged product: a common rotation of
# every period leaves the product, the likelihood and the isotropic AR(1) prior unchanged, so the
# coordinates are rotated to the principal axes of the average product, as the sampler does
.lde_pos_report <- function(ctx, pos) {
	if (isTRUE(ctx$sym)) return(.lde_pos_report_sym(ctx, pos))
	Tt <- ctx$Tt; R <- ctx$R
	na <- dim(pos$qu$mean)[1L]; nb <- dim(pos$qv$mean)[1L]
	Ubar <- matrix(0, na, R)
	Obar <- matrix(0, na, nb)
	for (t in seq_len(Tt)) {
		Ut <- matrix(pos$qu$mean[, , t], na, R)
		Vt <- matrix(pos$qv$mean[, , t], nb, R)
		Ubar <- Ubar + Ut / Tt
		Obar <- Obar + tcrossprod(Ut, Vt) / Tt
	}
	Qr <- diag(R)
	sv <- svd(Obar, nu = R, nv = R)
	# a numerically zero product (collapsed positions) has no principal axes; the gauge stays the identity
	if (all(is.finite(sv$d[seq_len(R)])) && sv$d[1L] > 1e-10) {
		Uref <- sweep(sv$u[, seq_len(R), drop = FALSE], 2L, sqrt(pmax(sv$d[seq_len(R)], 0)), "*")
		sr <- svd(crossprod(Ubar, Uref))
		Qr <- sr$u %*% t(sr$v)
		# the sign of a column is a gauge as well: the largest entry of the reference column is positive
		sg <- sign(Uref[cbind(apply(abs(Uref), 2L, which.max), seq_len(R))])
		sg[sg == 0] <- 1
		Qr <- sweep(Qr, 2L, sg, "*")
	}
	rot <- function(qs, n) {
		m <- array(0, dim = c(n, R, Tt))
		sd <- array(0, dim = c(n, R, Tt))
		for (t in seq_len(Tt)) {
			m[, , t] <- matrix(qs$mean[, , t], n, R) %*% Qr
			for (i in seq_len(n)) {
				sd[i, , t] <- sqrt(pmax(diag(crossprod(Qr, qs$sig[, , i, t] %*% Qr)), 0))
			}
		}
		list(m = m, sd = sd)
	}
	ru <- rot(pos$qu, na)
	rv <- rot(pos$qv, nb)
	list(U = ru$m, V = rv$m, U_sd = ru$sd, V_sd = rv$sd, rotation = Qr)
}

# the reported positions of an undirected network: the rotations that leave the product
# F S F' and the isotropic prior unchanged are the orthogonal matrices that preserve S,
# one in each sign group, so every group of columns is rotated to the principal axes
# (the leading eigenvectors of that sign) of the time-averaged product; V = F S
.lde_pos_report_sym <- function(ctx, pos) {
	Tt <- ctx$Tt; R <- ctx$R; sg <- ctx$sgn
	n <- dim(pos$qu$mean)[1L]
	Fbar <- matrix(0, n, R)
	Obar <- matrix(0, n, n)
	for (t in seq_len(Tt)) {
		Ft <- matrix(pos$qu$mean[, , t], n, R)
		Fbar <- Fbar + Ft / Tt
		Obar <- Obar + tcrossprod(sweep(Ft, 2L, sg, "*"), Ft) / Tt
	}
	Qr <- diag(R)
	eg <- eigen((Obar + t(Obar)) / 2, symmetric = TRUE)
	for (g in unique(sg)) {
		cc <- which(sg == g)
		ord <- if (g > 0) order(eg$values, decreasing = TRUE) else order(eg$values)
		take <- ord[seq_along(cc)]
		Uref <- sweep(eg$vectors[, take, drop = FALSE], 2L, sqrt(abs(eg$values[take])), "*")
		# (a numerically zero reference, as above, leaves the gauge of the group at the identity)
		if (all(is.finite(Uref)) && max(abs(Uref)) > 1e-10) {
			sr <- svd(crossprod(Fbar[, cc, drop = FALSE], Uref))
			Q <- sr$u %*% t(sr$v)
			sgc <- sign(Uref[cbind(apply(abs(Uref), 2L, which.max), seq_along(cc))])
			sgc[sgc == 0] <- 1
			Qr[cc, cc] <- sweep(Q, 2L, sgc, "*")
		}
	}
	m <- array(0, dim = c(n, R, Tt))
	sd <- array(0, dim = c(n, R, Tt))
	for (t in seq_len(Tt)) {
		m[, , t] <- matrix(pos$qu$mean[, , t], n, R) %*% Qr
		for (i in seq_len(n)) {
			sd[i, , t] <- sqrt(pmax(diag(crossprod(Qr, pos$qu$sig[, , i, t] %*% Qr)), 0))
		}
	}
	V <- m
	for (k in seq_len(R)) V[, k, ] <- sg[k] * m[, k, ]
	list(U = m, V = V, U_sd = sd, V_sd = sd, rotation = Qr)
}
