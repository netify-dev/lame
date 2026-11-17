# standard errors of the reported effects, positions, product and variance components of an
# empirical-bayes ALS fit. the effects, positions and product take the delta method under the
# joint Laplace covariance of (mu, beta, effects, gamma, factors) that .ae_eb_vcov forms,
# carried through the identified report (.ae_report_fold); the variance components take the
# curvature of the marginal likelihood of each hyperparameter at the fitted state

# columns of the factor block of the stacked vector (.ae_eb_layout), one per rank and side
.ae_eb_factor_cols <- function(m0, nr, nc, R, sym) {
	iu <- lapply(seq_len(R), function(r) m0 + (r - 1L) * nr + seq_len(nr))
	iv <- if (sym) iu else lapply(seq_len(R), function(r) m0 + R * nr + (r - 1L) * nc + seq_len(nc))
	list(iu = iu, iv = iv)
}

# directions of the factor block along which the product does not change: U -> U A,
# V -> V A^-T for a directed or bipartite fit, F -> F K S with K skew when symmetric
.ae_eb_gauge_dirs <- function(Fu, Fv, sg, sym) {
	nr <- nrow(Fu); R <- ncol(Fu)
	nf <- if (sym) nr * R else (nr + nrow(Fv)) * R
	at_u <- function(r) (r - 1L) * nr + seq_len(nr)
	at_v <- function(r) R * nr + (r - 1L) * nrow(Fv) + seq_len(nrow(Fv))
	out <- list()
	for (k in seq_len(R)) for (l in seq_len(R)) {
		if (sym && k >= l) next
		t_kl <- numeric(nf)
		if (sym) {
			t_kl[at_u(l)] <- Fu[, k] * sg[l]
			t_kl[at_u(k)] <- -Fu[, l] * sg[k]
		} else {
			t_kl[at_u(l)] <- Fu[, k]
			t_kl[at_v(k)] <- -Fv[, l]
		}
		out[[length(out) + 1L]] <- t_kl
	}
	if (!length(out)) return(NULL)
	do.call(cbind, out)
}

# covariance of the factor block with the gauge directions projected out and, with `centre`,
# each rank block centred over its actors, as the reported positions are
.ae_eb_factor_cov <- function(Sff, Td, sizes, centre = TRUE) {
	if (!is.null(Td)) {
		Q <- qr.Q(qr(Td))
		SQ <- Sff %*% Q
		Sff <- Sff - SQ %*% t(Q) - Q %*% t(SQ) + Q %*% (crossprod(Q, SQ) %*% t(Q))
	}
	ends <- cumsum(sizes)
	starts <- ends - sizes + 1L
	centre_rows <- function(X) {
		for (k in seq_along(sizes)) {
			ix <- starts[k]:ends[k]
			X[ix, ] <- sweep(X[ix, , drop = FALSE], 2L, colMeans(X[ix, , drop = FALSE]))
		}
		X
	}
	Sc <- if (centre) centre_rows(t(centre_rows(Sff))) else Sff
	(Sc + t(Sc)) / 2
}

# variance of every cell of the centred product x_i' S y_j for Gaussian rows x_i, y_j with
# means Xm, Ym: the delta-method term plus the second-order term of the product
.ae_eb_cell_var <- function(Xm, Ym, cov_xx, cov_yy, cov_xy, sg) {
	R <- ncol(Xm)
	out <- matrix(0, nrow(Xm), nrow(Ym))
	for (r in seq_len(R)) for (s in seq_len(R)) {
		w <- sg[r] * sg[s]
		dxx <- diag(cov_xx[[r]][[s]]); dyy <- diag(cov_yy[[r]][[s]])
		out <- out + w * (
			outer(dxx, Ym[, r] * Ym[, s]) +
				outer(Xm[, r] * Xm[, s], dyy) +
				2 * outer(Xm[, s], Ym[, r]) * cov_xy[[r]][[s]] +
				outer(dxx, dyy) +
				cov_xy[[r]][[s]] * cov_xy[[s]][[r]]
		)
	}
	pmax(out, 0)
}

# sds of the eigenvectors and eigenvalues of a symmetric product M = Fc S Fc' (Fc = C F, the
# rows centred by C or C the identity) under the covariance Sfu of the factor block F (rank-major):
# first-order perturbation du_r = (l_r I - M)^+ dM u_r and dl_r = u_r' dM u_r with
# dM = C dF S Fc' + Fc S dF' C. with the gauge-projected covariance Sc (blocks cov_uu) it also
# returns the shift of each eigenvalue from the product of the posterior means to the
# posterior mean of the eigenvalue
.ae_eb_eigen_sd <- function(Fc, sg, Sfu, C, cov_uu = NULL, Sc = NULL) {
	nr <- nrow(Fc); R <- ncol(Fc)
	M <- Fc %*% (sg * t(Fc))
	eg <- eigen((M + t(M)) / 2, symmetric = TRUE)
	idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
	Ue <- eg$vectors[, idx, drop = FALSE]; lam <- eg$values[idx]
	FS <- sweep(Fc, 2L, sg, "*")
	U_sd <- matrix(NA_real_, nr, R)
	L_sd <- L_shift <- rep(NA_real_, R)
	for (r in seq_len(R)) {
		Pr <- (diag(nr) - tcrossprod(Ue)) / lam[r]
		for (s in setdiff(seq_len(R), r)) Pr <- Pr + tcrossprod(Ue[, s]) / (lam[r] - lam[s])
		cr <- as.vector(crossprod(FS, Ue[, r]))
		K <- matrix(0, nr, nr * R)
		g <- numeric(nr * R)
		for (k in seq_len(R)) {
			ix <- (k - 1L) * nr + seq_len(nr)
			K[, ix] <- cr[k] * C + outer(FS[, k], Ue[, r])
			g[ix] <- 2 * cr[k] * Ue[, r]
		}
		J <- Pr %*% K
		U_sd[, r] <- sqrt(pmax(rowSums((J %*% Sfu) * J), 0))
		L_sd[r] <- sqrt(max(sum(g * as.vector(Sfu %*% g)), 0))
		# the eigenvalue of the product of the posterior means falls short of the posterior
		# mean of the eigenvalue by the second-order terms u' E[dF S dF'] u + tr((l I - M)^+ dM u u' dM')
		if (!is.null(Sc)) {
			L_shift[r] <- sum(vapply(seq_len(R), function(k) sg[k] * sum(Ue[, r] * (cov_uu[[k]][[k]] %*% Ue[, r])), 0)) +
				sum(J * (K %*% Sc))
		}
	}
	list(U = U_sd, L = L_sd, L_shift = if (!is.null(Sc)) L_shift, vectors = Ue, values = lam)
}

# below the signal threshold of the product (`info`, .ae_product_statistic()) the position
# posterior is far from gaussian along the directions the product barely constrains, and the
# linearised covariance is too narrow; the variational posterior of each row keeps that spread,
# so the position sds of a directed or bipartite fit are at least the variational row sds
# `U_var`, `V_var` (same shape as U_sd, V_sd); a static fit passes them with the spread each row
# takes from the other side (.ae_eb_weak_row_var())
.ae_eb_widen_weak_positions <- function(out, info, U_var, V_var, align = FALSE) {
	stat <- .ae_product_statistic(out$family, info)
	if (!isTRUE(out$symmetric) && is.finite(stat) && stat < info[["threshold"]] &&
			identical(dim(U_var), dim(out$U_sd)) && identical(dim(V_var), dim(out$V_sd))) {
		out$U_sd <- pmax(out$U_sd, U_var)
		out$V_sd <- pmax(out$V_sd, V_var)
		if (isTRUE(align)) out <- .ae_eb_scale_align_sd(out)
	}
	out
}

# the reported positions of a weak product are compared in their own gauge, where the scale of each
# column is split between the sides by the prior: a least-squares alignment of a posterior draw
# (U + e_U, V + e_V) rescales its column by c = e^t, t = (E_V - E_U) / (N_U + N_V + 2 E_U + 2 E_V) to
# first order (E the summed variances of a side, N its squared norm), so the side that carries the noise
# is shrunk and the other grows: the error of an entry is e^t e_U + (e^t - 1) U and e^-t e_V +
# (e^-t - 1) V, whose mean squares replace the sds. a balanced split (directed fits) gives t near 0.
# taken by the dynamic report (`align`); the static floor, from the joint covariance of the fit, already
# matched the aligned spread of the exact posterior of a weak bipartite product, and the step widens it
.ae_eb_scale_align_sd <- function(out) {
	U <- out$U; V <- out[["V"]]
	if (!is.matrix(U) || !is.matrix(V) || !identical(dim(out$U_sd), dim(U)) || !identical(dim(out$V_sd), dim(V)) ||
			!all(is.finite(U)) || !all(is.finite(V))) return(out)
	# unequal ranks (bipartite): only the columns the two sides share form the product
	for (r in seq_len(min(ncol(U), ncol(V)))) {
		eu <- sum(out$U_sd[, r]^2); ev <- sum(out$V_sd[, r]^2)
		den <- sum(U[, r]^2) + sum(V[, r]^2) + 2 * eu + 2 * ev
		if (!is.finite(den) || !(den > 0)) next
		t <- (ev - eu) / den
		out$U_sd[, r] <- sqrt(exp(2 * t) * out$U_sd[, r]^2 + (exp(t) - 1)^2 * U[, r]^2)
		out$V_sd[, r] <- sqrt(exp(-2 * t) * out$V_sd[, r]^2 + (1 - exp(-t))^2 * V[, r]^2)
	}
	out
}

# variances of the rows of a weak directed or bipartite product beside the linearised ones: the
# variational covariance `Sig_i` of a row conditions on the other side at its mean, while the
# row's conditional mean Sig_i sum_j w_ij z_ij v_j moves with the other side; over the
# variational covariances `Om_j` of that side the move has variance Sig_i (sum_j w_ij Om_j) Sig_i
# (w the cell information `W`, E[w z^2] = 1), which is added. the rows are taken independent,
# and the directions that leave the product unchanged are projected out and each rank block is
# centred, as for the linearised covariance (.ae_eb_factor_cov()); only the diagonal is formed,
# from the block-diagonal covariance D and the gauge basis Q: diag(C P D P C), P = I - Q Q' (C the
# identity with `centre = FALSE`).
# `Fu`, `Fv` the factors of the fitted state; returns the variances of the factor block (rank-major)
.ae_eb_weak_row_var <- function(Fu, Fv, Sig, Om, W, centre = TRUE) {
	R <- ncol(Fu); nr <- nrow(Fu); nc <- nrow(Fv)
	if (length(dim(W)) == 3L) W <- rowSums(W, dims = 2L)
	W[!is.finite(W)] <- 0
	# a column held at zero has no conditional mean to move: it keeps its variational covariance
	held <- colSums(Fu^2) == 0 | colSums(Fv^2) == 0
	row_cov <- function(S, O, Wt) {
		M <- matrix(O, R * R) %*% Wt
		vapply(seq_len(dim(S)[3L]), function(i) {
			Si <- matrix(S[, , i], R, R)
			Mi <- matrix(M[, i], R, R)
			Mi[held, ] <- 0; Mi[, held] <- 0
			Si + Si %*% Mi %*% Si
		}, matrix(0, R, R))
	}
	Cu <- array(row_cov(Sig, Om, t(W)), c(R, R, nr))
	Cv <- array(row_cov(Om, Sig, W), c(R, R, nc))
	at_u <- function(r) (r - 1L) * nr + seq_len(nr)
	at_v <- function(r) R * nr + (r - 1L) * nc + seq_len(nc)
	dD <- c(vapply(seq_len(R), function(r) Cu[r, r, ], numeric(nr)), vapply(seq_len(R), function(r) Cv[r, r, ], numeric(nc)))
	# a held column gives no gauge direction
	gd <- .ae_eb_gauge_dirs(Fu, Fv, rep(1, R), FALSE)
	if (!is.null(gd)) gd <- gd[, colSums(gd^2) > 0, drop = FALSE]
	sizes <- c(rep(nr, R), rep(nc, R))
	blk <- rep(seq_along(sizes), sizes)
	ones <- outer(blk, seq_along(sizes), "==") * 1
	if (is.null(gd) || !ncol(gd)) {
		if (!centre) return(dD)
		# centred block-diagonal D: the row sums of D within a block are its diagonal
		return(dD - 2 * dD / sizes[blk] + as.vector(rowsum(dD, blk))[blk] / sizes[blk]^2)
	}
	Q <- qr.Q(qr(gd))
	# G = D Q, D block diagonal by row
	G <- matrix(0, nrow(Q), ncol(Q))
	for (r in seq_len(R)) for (s in seq_len(R)) {
		G[at_u(r), ] <- G[at_u(r), ] + Cu[r, s, ] * Q[at_u(s), , drop = FALSE]
		G[at_v(r), ] <- G[at_v(r), ] + Cv[r, s, ] * Q[at_v(s), , drop = FALSE]
	}
	H <- crossprod(Q, G)
	# A = D - Q G' - G Q' + Q H Q': its diagonal, its row sums within each block and the block totals
	dA <- dD - 2 * rowSums(Q * G) + rowSums((Q %*% H) * Q)
	# positions reported uncentred (a dynamic fit) keep the projected diagonal
	if (!centre) return(dA)
	Qb <- crossprod(Q, ones); Gb <- crossprod(G, ones)
	# A 1_b on the rows of block b (D 1_b there is the diagonal of D)
	rA <- dD - rowSums(Q * t(Gb)[blk, , drop = FALSE]) - rowSums(G * t(Qb)[blk, , drop = FALSE]) +
		rowSums((Q %*% H) * t(Qb)[blk, , drop = FALSE])
	tA <- as.vector(rowsum(rA, blk))
	dA - 2 * rA / sizes[blk] + tA[blk] / sizes[blk]^2
}

# sds and intervals of the eigenvalues of a symmetric product: the linearised sds `L_sd`, the
# interval centred on the posterior mean of each eigenvalue (L + L_shift, .ae_eb_eigen_sd()); for
# rank one, and for a column of a higher rank whose eigenvalue or linearised sd vanishes, the
# eigenvalue is +-|C f|^2 with the rows of f gaussian under the variational posterior (row
# covariances `Sig`, R x R x n), a quadratic form whose mean and variance give a gamma interval.
# `centred`: the rows of f are centred over the actors (C the centring matrix), else C = I;
# `every`: every column takes the gamma interval; `L_extra`: variance each eigenvalue takes from
# the uncertainty of the scale of the factor prior, added to either
.ae_eb_eigen_interval <- function(out, L_sd, L_shift = NULL, Sig = NULL, centred = TRUE, every = FALSE, L_extra = NULL) {
	ex <- if (length(L_extra) == length(L_sd) && all(is.finite(L_extra))) pmax(L_extra, 0) else rep(0, length(L_sd))
	L_sd <- sqrt(L_sd^2 + ex)
	out$L_sd <- L_sd
	zq <- stats::qnorm(0.975)
	if (!is.null(L_shift) && length(out$L) == length(L_sd)) {
		out$L_ci <- cbind(lower = out$L + L_shift - zq * L_sd, upper = out$L + L_shift + zq * L_sd)
	}
	R <- length(out$L)
	zero <- abs(out$L) <= sqrt(.Machine$double.eps) * max(1, abs(out$L))
	if (R >= 1L && !is.null(Sig) && identical(dim(Sig)[3L], nrow(out$U))) {
		n <- nrow(out$U)
		if (is.null(out$L_ci)) out$L_ci <- cbind(lower = out$L - zq * out$L_sd, upper = out$L + zq * out$L_sd)
		for (r in seq_len(R)) {
			if (!every && R > 1L && !zero[r] && isTRUE(out$L_sd[r] > 0)) next
			sv <- pmax(as.vector(Sig[r, r, ]), 0) * (if (centred) 1 - 1 / n else 1)
			fc <- out$U[, r] * sqrt(abs(out$L[r]))
			m <- abs(out$L[r]) + sum(sv)
			v <- 2 * sum(sv^2) + 4 * sum(fc^2 * sv) + ex[r]
			if (!is.finite(m) || !is.finite(v) || !(v > 0)) next
			q <- stats::qgamma(c(0.025, 0.975), shape = m^2 / v, scale = v / m) * (if (out$L[r] < 0) -1 else 1)
			out$L_sd[r] <- sqrt(v)
			out$L_ci[r, ] <- sort(q)
		}
	}
	out
}

# below the signal threshold of the product (`info`, .ae_product_statistic()) an undirected fit
# determines neither the sign of an eigenvalue nor the direction of its eigenvector: the sampler's
# posterior of the leading eigenvalue of a weak product puts close to half its mass on each sign,
# and the reported eigenvector is the leading direction of the noise as much as of the product.
# with the sign not determined the posterior of an eigenvalue is the even mixture of +-|C f|^2, the
# quadratic form of .ae_eb_eigen_interval() (variational row covariances `Sig`), whose central 95%
# interval is +-q with q the 95% quantile of |C f|^2, taken for every column. with the direction not
# determined the eigenvector is uniform on the unit sphere of centred vectors, signed as the
# estimate, and the error of the reported entry u_i has mean square u_i^2 (1 - 2 E|c|) + 1 / n (c the
# cosine of the two directions, E|c| = Gamma((n - 1) / 2) / (sqrt(pi) Gamma(n / 2))): the position
# sds are at least its root. the statistic is read per column: the fit's statistic split over the
# columns by their squared eigenvalues (the squared norm of the double-centred product is their sum),
# so a column below the threshold takes this rule beside one the data determine, whichever side of
# the threshold the fit as a whole is. `centred = FALSE`: positions reported uncentred (a dynamic
# fit), the direction uniform over all unit vectors (E|c| with n - 1 replaced by n) and the rows of f
# uncentred; `L_extra` the variance each eigenvalue takes from the scale of the factor prior;
# `margin`: the multiple of the threshold below which a column's sign counts as undetermined (1 for a
# static fit: of 2822 simulated columns above the threshold, normal, binary and count, one had the wrong sign)
.ae_eb_weak_symmetric <- function(out, info, Sig, centred = TRUE, L_extra = NULL, margin = 1, cells = TRUE) {
	stat <- .ae_product_statistic(out$family, info)
	if (!isTRUE(out$symmetric) || !is.finite(stat) || !is.matrix(out$U)) return(out)
	n <- nrow(out$U)
	R <- ncol(out$U)
	if (n < 3L || length(out$L) != R) return(out)
	l2 <- ifelse(is.finite(out$L), out$L^2, 0)
	stat_r <- if (sum(l2) > 0) stat * R * l2 / sum(l2) else rep(stat, R)
	weak <- stat_r < margin * info[["threshold"]]
	if (!any(weak)) return(out)
	if (identical(dim(out$U_sd), dim(out$U))) {
		d <- if (centred) n - 1 else n
		ec <- exp(lgamma(d / 2) - lgamma((d + 1) / 2)) / sqrt(pi)
		fl <- sqrt(out$U^2 * (1 - 2 * ec) + 1 / n)
		out$U_sd[, weak] <- pmax(out$U_sd[, weak], fl[, weak])
		if (identical(dim(out$V_sd), dim(out$U))) out$V_sd <- out$U_sd
	}
	# its product cells: with sign and direction independent of the truth's, the error of the cell
	# l u_i u_j has mean square (L U_i U_j)^2 + l^2 E[u_i^2 u_j^2], E[u_i^2 u_j^2] = 1 / (d (d + 2))
	# off the diagonal and 3 / (d (d + 2)) on it for a direction uniform on the unit sphere of d
	# dimensions, l^2 at the reported L^2 (the quadratic form's second moment carries the noise the
	# weak fit absorbs); the cell sds are at least its root
	if (isTRUE(cells) && identical(dim(out$UVPM_sd), c(n, n))) {
		d <- if (centred) n - 1 else n
		eu <- (1 + 2 * diag(n)) / (d * (d + 2))
		cf <- matrix(0, n, n)
		for (r in which(weak)) {
			if (!is.finite(out$L[r])) next
			cf <- cf + (out$L[r] * tcrossprod(out$U[, r]))^2 + out$L[r]^2 * eu
		}
		if (all(is.finite(cf))) {
			out$UVPM_sd <- pmax(out$UVPM_sd, sqrt(cf))
			if (!is.null(out$ULUPM_sd)) out$ULUPM_sd <- out$UVPM_sd
		}
	}
	if (!is.null(out$L_ci) && !is.null(Sig) && identical(dim(Sig)[3L], n)) {
		for (r in which(weak)) {
			sv <- pmax(as.vector(Sig[r, r, ]), 0) * (if (centred) 1 - 1 / n else 1)
			fc <- out$U[, r] * sqrt(abs(out$L[r]))
			m <- abs(out$L[r]) + sum(sv)
			v <- 2 * sum(sv^2) + 4 * sum(fc^2 * sv) + (if (length(L_extra) == length(out$L) && is.finite(L_extra[r])) max(L_extra[r], 0) else 0)
			if (!is.finite(m) || !is.finite(v) || !(v > 0)) next
			q <- stats::qgamma(0.95, shape = m^2 / v, scale = v / m)
			out$L_sd[r] <- sqrt(v)
			out$L_ci[r, ] <- c(-q, q)
		}
	}
	out
}

# sds of the reported a, b (cross-sectional and centred over the observed actors), positions,
# eigenvalues and product cells, from the joint covariance `S` (the inverse of the extended
# information) at the fitted state `fit`; `vb` is the variational posterior, if any
# second-order term of the expected log likelihood of a count fit, 0.5 sum_ij w_ij t_ij, returned
# as the cell weights 0.5 w (w the cell information) and the variances t: t_ij is the part of the
# posterior variance of the product cell that the Laplace (linear) term omits, the Isserlis terms
# sum_rs [C(u_ir, u_is) C(v_jr, v_js) + C(u_ir, v_js) C(u_is, v_jr)] of the factor covariance S
# with the gauge projected out. where the data inform the product weakly, t_ij is the whole of
# its variance, and without it the noise variance of the count drops by about its mean
.ae_pln_second_order <- function(S, ctx, ab, fit, R, wsum) {
	m0 <- ctx$m; nr <- ctx$nr; nc <- ctx$nc
	if (R < 1L || is.null(fit$U) || nrow(S) <= m0) return(NULL)
	sym <- identical(ab$kind, "symmetric")
	fc <- .ae_eb_factor_cols(m0, nr, nc, R, sym)
	if (sym) {
		Fu <- sweep(fit$U, 2L, sqrt(abs(fit$L)), "*")
		sg <- ifelse(fit$L < 0, -1, 1)
		Fv <- Fu
	} else {
		Fu <- fit$U; Fv <- fit[["V"]]
		sg <- rep(1, R)
	}
	jf <- (m0 + 1L):nrow(S)
	sizes <- if (sym) rep(nr, R) else c(rep(nr, R), rep(nc, R))
	Sc <- .ae_eb_factor_cov(S[jf, jf, drop = FALSE], .ae_eb_gauge_dirs(Fu, Fv, sg, sym), sizes)
	ju <- lapply(fc$iu, function(ix) ix - m0)
	jv <- if (sym) ju else lapply(fc$iv, function(ix) ix - m0)
	t2 <- matrix(0, nr, nc)
	for (r in seq_len(R)) for (s in seq_len(R)) {
		dxx <- diag(Sc[ju[[r]], ju[[s]], drop = FALSE])
		dyy <- diag(Sc[jv[[r]], jv[[s]], drop = FALSE])
		t2 <- t2 + sg[r] * sg[s] * (outer(dxx, dyy) + Sc[ju[[r]], jv[[s]], drop = FALSE] * Sc[ju[[s]], jv[[r]], drop = FALSE])
	}
	t2 <- pmax(t2, 0)
	# noise the weak columns take up beyond the regular count (.ae_pln_bulk_charge)
	if (!sym) t2 <- t2 + .ae_pln_bulk_charge(fit$score_sum, tcrossprod(Fu, Fv), wsum, R)
	list(w = 0.5 * ab$nu * wsum, t = t2)
}

# noise variance the factor columns of a directed or bipartite count fit take up beyond the regular
# nr + nc degrees of freedom each, spread over the observed cells in proportion to their working
# variance. the product is fitted to M = sqrt(w) (s / w + P), the mean working residual plus the
# product in units of the cell working sd, whose noise has the variance tau of its bulk (the mean
# square of the singular values past the first R and below the edge of that noise). a column of
# signal strength t (its sum of squares) has squared singular value x = (t + nr tau)(t + nc tau) / t
# and a truncation to it takes up x - t = (nr + nc) tau + nr nc tau^2 / t of the noise, so its excess
# over the regular count is x - (nr + nc) tau - t: large near the edge x = (sqrt(nr) + sqrt(nc))^2 tau,
# vanishing for a strong column. below the edge the strength is not identified and is held at its
# edge value sqrt(nr nc) tau, where the excess falls with x. a column the fit shrinks to f of the
# singular value sigma takes up the share 2 f sigma - f^2 of the sigma^2 a truncation takes, and the
# excess with it (none for a column held at zero); 0 without the cell scores
.ae_pln_bulk_charge <- function(score, P, wsum, R) {
	nr <- nrow(P); nc <- ncol(P)
	if (R < 1L || is.null(score) || !identical(dim(score), dim(P)) || !identical(dim(wsum), dim(P))) return(0)
	obs <- is.finite(wsum) & wsum > 0
	if (!any(obs) || !all(is.finite(score[obs]))) return(0)
	M <- PM <- matrix(0, nr, nc)
	M[obs] <- sqrt(wsum[obs]) * (score[obs] / wsum[obs] + P[obs])
	PM[obs] <- sqrt(wsum[obs]) * P[obs]
	if (!all(is.finite(M)) || !all(is.finite(PM))) return(0)
	ex <- .ae_bulk_excess(M, PM, R, hold = 1)
	out <- matrix(0, nr, nc)
	out[obs] <- ex / (sum(obs) * wsum[obs])
	out
}

# the excess of .ae_pln_bulk_charge() in the sum-of-squares units of M, the fitted matrix (the
# product plus its residual), from M and the fitted product PM; 0 when it cannot be formed. below
# the edge the strength is bounded by sqrt(nr nc) but not located, and is held at `hold` times that
# bound: a count fit holds it at the bound (1); a normal fit at the middle of the range (1/2), as its
# variational posterior charges a weak column less than the regular count
.ae_bulk_excess <- function(M, PM, R, hold = 1) {
	nr <- nrow(M); nc <- ncol(M)
	sv2 <- svd(M, nu = 0L, nv = 0L)$d^2
	f <- svd(PM, nu = 0L, nv = 0L)$d
	edge <- (sqrt(nr) + sqrt(nc))^2
	# the bulk: the components past the first R that lie below the edge of its own noise
	K <- R; tau <- NA_real_
	for (it in seq_len(10L)) {
		dof <- (nr - K - 1) * (nc - K - 1)
		if (K >= length(sv2) || !(dof > 0)) return(0)
		tau <- sum(sv2[-seq_len(K)]) / dof
		if (!is.finite(tau) || !(tau > 0)) return(0)
		K1 <- max(R, sum(sv2 > edge * tau))
		if (K1 == K) break
		K <- K1
	}
	x <- sv2[seq_len(R)] / tau
	d <- x - nr - nc
	e <- d - ifelse(x >= edge, (d + sqrt(pmax(d^2 - 4 * nr * nc, 0))) / 2, hold * sqrt(nr * nc))
	q <- pmin(f[seq_len(R)] / sqrt(sv2[seq_len(R)]), 1)
	ex <- tau * sum(e * (2 * q - q^2))
	if (!is.finite(ex)) return(0)
	ex
}

# the excess (.ae_bulk_excess()) of a normal fit at its evaluation `ev`, in the
# units of the noise sum of squares: M is the period mean of the residual plus the product, times the
# square root of the periods a cell is observed in, so its noise has the variance ve of one cell
.ae_gauss_bulk_excess <- function(cx, ev, R) {
	n <- rowSums(cx$obs, dims = 2L)
	obs <- n > 0
	r <- cx$y - ev$eta
	r[!cx$obs] <- 0
	rs <- rowSums(r, dims = 2L)
	M <- PM <- matrix(0, cx$nr, cx$nc)
	M[obs] <- sqrt(n[obs]) * (rs[obs] / n[obs] + ev$Pm[obs])
	PM[obs] <- sqrt(n[obs]) * ev$Pm[obs]
	if (!any(obs) || !all(is.finite(M)) || !all(is.finite(PM))) return(0)
	.ae_bulk_excess(M, PM, R, hold = 0.5)
}

.ae_eb_report_se <- function(S, ctx, ab, fit, vb = NULL, R = 0L, W = NULL, X = NULL) {
	nr <- ctx$nr; nc <- ctx$nc; m0 <- ctx$m
	sym <- identical(ab$kind, "symmetric")
	d <- nrow(S)
	use_f <- R > 0L && !is.null(fit$U) && d > m0
	ia <- ctx$i_s[seq_len(nr)]
	ib <- if (sym) ia else ctx$i_s[nr + seq_len(nc)]
	pr <- if (is.null(ab$Wr)) 0L else ncol(ab$Wr)
	pc <- if (is.null(ab$Wc)) 0L else ncol(ab$Wc)
	# effects net of their covariate part, a = s - M gamma
	Ga <- matrix(0, nr, d)
	Ga[cbind(seq_len(nr), ia)] <- 1
	Gb <- NULL
	if (sym) {
		if (ctx$pg > 0L) Ga[, ctx$i_g] <- -ctx$M
	} else {
		Gb <- matrix(0, nc, d)
		Gb[cbind(seq_len(nc), ib)] <- 1
		if (pr) Ga[, ctx$i_g[seq_len(pr)]] <- -ab$Wr
		if (pc) Gb[, ctx$i_g[pr + seq_len(pc)]] <- -ab$Wc
	}
	# a longitudinal fit reports the effects of the observed actors at mean zero: their gradient
	# loses the mean gradient of the observed actors, taken before the product's means are added
	mean_act <- function(G, act) if (is.null(G) || !any(act)) NULL else colMeans(G[act, , drop = FALSE])
	ma <- mean_act(Ga, ab$act_a)
	mb <- if (sym) NULL else mean_act(Gb, ab$act_b)
	out <- list()
	if (use_f) {
		fc <- .ae_eb_factor_cols(m0, nr, nc, R, sym)
		if (sym) {
			Fu <- if (is.null(vb)) sweep(fit$U, 2L, sqrt(abs(fit$L)), "*") else vb$S$F
			sg <- if (is.null(vb)) ifelse(fit$L < 0, -1, 1) else vb$S$sg
			Fv <- Fu
		} else {
			Fu <- fit$U; Fv <- fit[["V"]]
			sg <- rep(1, R)
		}
		ubar <- colMeans(Fu); vbar <- colMeans(Fv)
		Uc <- sweep(Fu, 2L, ubar); Vc <- sweep(Fv, 2L, vbar)
		# the row and column means of the product the report moves into a and b
		add_fold <- function(G, Xc, ybar, ix, iy, n_x, n_y) {
			for (r in seq_len(R)) {
				G[, ix[[r]]] <- G[, ix[[r]]] + sg[r] * ybar[r] * (diag(n_x) - 1 / n_x)
				G[, iy[[r]]] <- G[, iy[[r]]] + sg[r] * outer(Xc[, r], rep(1 / n_y, n_y))
			}
			G
		}
		if (sym) {
			# symmetric: the row and column parts agree, so the effect takes one of them
			Ga <- add_fold(Ga, Uc, ubar, fc$iu, fc$iu, nr, nr)
		} else {
			Ga <- add_fold(Ga, Uc, vbar, fc$iu, fc$iv, nr, nc)
			Gb <- add_fold(Gb, Vc, ubar, fc$iv, fc$iu, nc, nr)
		}
		# factor block with the gauge projected out, centred as reported
		jf <- (m0 + 1L):d
		sizes <- if (sym) rep(nr, R) else c(rep(nr, R), rep(nc, R))
		Sc <- .ae_eb_factor_cov(
			S[jf, jf, drop = FALSE],
			.ae_eb_gauge_dirs(Fu, Fv, sg, sym),
			sizes
		)
		ju <- lapply(fc$iu, function(ix) ix - m0)
		jv <- lapply(fc$iv, function(ix) ix - m0)
		blocks <- function(ix, iy) lapply(seq_len(R), function(r) lapply(seq_len(R), function(s) Sc[ix[[r]], iy[[s]], drop = FALSE]))
		cov_uu <- blocks(ju, ju)
		cov_vv <- if (sym) cov_uu else blocks(jv, jv)
		cov_uv <- if (sym) cov_uu else blocks(ju, jv)
		out$cells <- sqrt(.ae_eb_cell_var(Uc, Vc, cov_uu, cov_vv, cov_uv, sg))
		# a directed binary column held at zero: its scale reads the rows' full second moments
		# (.ae_mf_scale_target), under which a column below the bulk edge falls to the floor of its
		# prior, so the prior no longer carries its cell variance. the strength of such a column (its
		# cell variance lambda) is not identified below the edge, where lambda reaches the edge of the
		# working noise 1 / w, lambda_c = 1 / (w sqrt(nr nc)); at the middle of that range, as the bulk
		# rule holds a normal fit's weak column (.ae_bulk_excess), it adds lambda_c / 2 to every cell
		if (identical(ab$kind, "directed") && is.null(ab$family)) {
			held <- colSums(Fu^2) == 0 | colSums(Fv^2) == 0
			wbar <- mean(ctx$wsum[is.finite(ctx$wsum) & ctx$wsum > 0])
			if (any(held) && is.finite(wbar) && wbar > 0) out$held_cell_var <- sum(held) / (2 * wbar * sqrt(nr * nc))
		}
		if (sym) {
			# positions are the eigenvectors of the centred product (.ae_eb_eigen_sd)
			es <- .ae_eb_eigen_sd(Uc, sg, S[jf, jf, drop = FALSE], diag(nr) - 1 / nr, cov_uu, Sc)
			out$U <- es$U
			out$L <- es$L
			out$L_shift <- es$L_shift
		} else {
			out$U <- sqrt(pmax(vapply(seq_len(R), function(r) diag(cov_uu[[r]][[r]]), numeric(nr)), 0))
			out[["V"]] <- sqrt(pmax(vapply(seq_len(R), function(r) diag(cov_vv[[r]][[r]]), numeric(nc)), 0))
			dim(out$U) <- c(nr, R); dim(out[["V"]]) <- c(nc, R)
			# below the signal threshold, the rows with the spread they take from the other side
			# (.ae_eb_weak_row_var()), for the position sds of .ae_eb_attach_se()
			# (the variational row covariances the fit reports, fit$eb$factor_posterior)
			info <- fit$eb_mf$record$information
			fp <- fit$eb_mf$posterior
			stat <- .ae_product_statistic(ab$family, info)
			if (is.finite(stat) && stat < info[["threshold"]] && !is.null(W) &&
					identical(dim(fp$Sig), c(R, R, nr)) && identical(dim(fp$Om), c(R, R, nc))) {
				vw <- tryCatch(.ae_eb_weak_row_var(Fu, Fv, fp$Sig, fp$Om, W), error = function(e) NULL)
				if (length(vw) == (nr + nc) * R && all(is.finite(vw))) {
					dw <- sqrt(pmax(vw, 0))
					out$U_weak <- matrix(dw[seq_len(nr * R)], nr, R)
					out$V_weak <- matrix(dw[nr * R + seq_len(nc * R)], nc, R)
				}
			}
		}
	}
	out$sab_info <- .ae_eb_sab_info(S[seq_len(m0), seq_len(m0), drop = FALSE], ctx, ab, fit$sab)
	if (!is.null(W)) {
		Fu_s <- if (use_f) Fu else NULL
		Fv_s <- if (use_f) Fv else NULL
		out$df_eff <- .ae_eb_cell_df(S, ctx, ab, if (use_f) R else 0L, W, X, Fu_s, Fv_s, if (use_f) sg)
		# the coefficients and effects alone, beside which the product's own df is counted
		out$df_linear <- if (use_f) .ae_eb_cell_df(S, ctx, ab, 0L, W, X) else out$df_eff
	}
	# variances of G theta and, for the observed rows, of (G - 1 m') theta
	quad <- function(G, m, act) {
		if (is.null(G)) return(NULL)
		GS <- G %*% S
		v <- rowSums(GS * G)
		vc <- v
		if (!is.null(m)) {
			vc[act] <- v[act] - 2 * as.vector(GS[act, , drop = FALSE] %*% m) + sum(m * as.vector(S %*% m))
		}
		list(sqrt(pmax(v, 0)), sqrt(pmax(vc, 0)))
	}
	qa <- quad(Ga, ma, ab$act_a)
	qb <- if (sym) qa else quad(Gb, mb, ab$act_b)
	out$a <- qa[[1L]]; out$a_centred <- qa[[2L]]
	out$b <- qb[[1L]]; out$b_centred <- qb[[2L]]
	out
}

# fisher information of the effect covariance (va, cab, vb) in the marginal likelihood of
# the data with the effects integrated out: 0.5 tr(B D_j B D_k) with B = Lam - Lam P Lam,
# Lam the prior precision of the effects, P their posterior covariance marginal over the
# coefficients and factors and D_j the derivative of the prior covariance; the inverse-wishart
# hyperprior adds eta0 / 2 tr(Sab^-1 d_j Sab^-1 d_k)
.ae_eb_sab_info <- function(S0, ctx, ab, sab) {
	P <- .ae_eb_css(ctx, S0)
	Lam <- .ae_eb_precision(sab, ab)
	attr(Lam, "ae_prior") <- NULL
	B <- Lam - Lam %*% P %*% Lam
	nr <- ab$nr
	q <- nrow(B)
	ia <- seq_len(nr)
	ib <- if (identical(ab$kind, "symmetric")) ia else nr + seq_len(ab$nc)
	# B D_j: the columns of B that D_j selects, swapped for the covariance
	BD <- list(va = { M <- matrix(0, q, q); M[, ia] <- B[, ia]; M })
	dd <- list(va = matrix(c(1, 0, 0, 0), 2L))
	if (identical(ab$kind, "directed")) {
		BD$cab <- { M <- matrix(0, q, q); M[, ia] <- B[, ib]; M[, ib] <- B[, ia]; M }
		dd$cab <- matrix(c(0, 1, 1, 0), 2L)
	}
	if (!identical(ab$kind, "symmetric")) {
		BD$vb <- { M <- matrix(0, q, q); M[, ib] <- B[, ib]; M }
		dd$vb <- matrix(c(0, 0, 0, 1), 2L)
	}
	k <- names(BD)
	Si <- if (identical(ab$kind, "symmetric")) diag(1 / sab[1L, 1L], 2L) else solve(sab)
	info <- matrix(0, length(k), length(k), dimnames = list(k, k))
	for (j in k) for (l in k) {
		info[j, l] <- 0.5 * sum(BD[[j]] * t(BD[[l]])) +
			0.5 * ab$eta0 * sum(diag(Si %*% dd[[j]] %*% Si %*% dd[[l]]))
	}
	info
}

# effective number of mean parameters, sum over observed cells of the information weight
# times the posterior variance of the cell's predictor (the trace of the hat matrix)
.ae_eb_cell_df <- function(S, ctx, ab, R, W, X, Fu = NULL, Fv = NULL, sg = NULL) {
	nr <- ctx$nr; nc <- ctx$nc; m0 <- ctx$m
	sym <- identical(ab$kind, "symmetric")
	ia <- ctx$i_s[seq_len(nr)]
	ib <- if (sym) ia else ctx$i_s[nr + seq_len(nc)]
	fc <- if (R > 0L) .ae_eb_factor_cols(m0, nr, nc, R, sym)
	Tt <- dim(W)[3L]
	p <- ctx$pa - 1L
	total <- 0
	for (t in seq_len(Tt)) {
		wt <- W[, , t]
		ok <- which(is.finite(wt) & wt > 0)
		if (!length(ok)) next
		I <- (ok - 1L) %% nr + 1L; J <- (ok - 1L) %/% nr + 1L
		pos <- cbind(ctx$i_th[1L], if (p > 0L) matrix(ctx$i_th[-1L], length(ok), p, byrow = TRUE), ia[I], ib[J])
		one <- rep(1, length(ok))
		g <- cbind(one, if (p > 0L) matrix(X[, , , t], nr * nc, p)[ok, , drop = FALSE], one, one)
		for (r in seq_len(R)) {
			pos <- cbind(pos, fc$iu[[r]][I], fc$iv[[r]][J])
			g <- cbind(g, sg[r] * Fv[J, r], sg[r] * Fu[I, r])
		}
		v <- numeric(length(ok))
		for (x in seq_len(ncol(g))) for (y in seq_len(ncol(g))) v <- v + g[, x] * g[, y] * S[cbind(pos[, x], pos[, y])]
		total <- total + sum(wt[ok] * pmax(v, 0))
	}
	total
}

# sd and interval of each estimated variance component: the inverse fisher information of the
# effect covariance (.ae_eb_sab_info), the inverse-gamma posterior of a normal noise variance
# less the effective number of mean parameters, the curvature of the marginal likelihood
# profile of a count noise variance and the pair-likelihood sd of the dyadic correlation
# inflated by the same loss of degrees of freedom.
# intervals are on the log scale for variances, the fisher z scale for rho and score-type
# (.ae_eb_cab_score_ci), inside its bound, for cab
.ae_eb_vc_se <- function(out, ab, rs = NULL) {
	VC <- out$VC
	se <- setNames(rep(NA_real_, length(VC)), names(VC))
	lo <- hi <- se
	z <- stats::qnorm(0.975)
	# `mean`: the value is a posterior mean, and the interval is centred on the log scale at the
	# median v exp(-s^2 / 2) of the lognormal with that mean and log-scale sd s = sd / v
	log_ci <- function(k, sd, mean = FALSE) {
		v <- VC[[k]]
		if (!is.finite(v) || v <= 0 || !is.finite(sd) || sd <= 0) return(invisible())
		s <- sd / v
		m <- if (isTRUE(mean)) v * exp(-s^2 / 2) else v
		se[[k]] <<- sd
		lo[[k]] <<- m * exp(-z * s)
		hi[[k]] <<- m * exp(z * s)
	}
	z_ci <- function(k, r, sd) {
		if (!is.finite(r) || abs(r) >= 1 || !is.finite(sd) || sd <= 0) return(invisible())
		sd_z <- sd / (1 - r^2)
		se[[k]] <<- sd
		lo[[k]] <<- tanh(atanh(r) - z * sd_z)
		hi[[k]] <<- tanh(atanh(r) + z * sd_z)
	}
	sab <- out$eb$sab
	info <- rs$sab_info
	if (!is.null(sab) && is.matrix(info)) {
		V <- tryCatch(solve(info), error = function(e) NULL)
		if (!is.null(V)) {
			sd_of <- function(k) if (k %in% rownames(V) && V[k, k] > 0) sqrt(V[k, k]) else NA_real_
			# the effect variances are inverse-wishart posterior means (.ae_eb_update_sab)
			log_ci("va", sd_of("va"), mean = TRUE)
			if ("vb" %in% rownames(V)) log_ci("vb", sd_of("vb"), mean = TRUE)
			# the covariance takes a score-type interval kept inside |cab| < sqrt(va vb): its sampling
			# variance grows with the covariance itself, (va vb + cab^2) / n for normal effects, so
			# the interval is the set of values c with |cab - c| <= z sd(c), sd(c) that variance read
			# at c. a wald interval reads it at the estimate and is widest where the estimate is far
			# from the truth, which over-covers
			sd_c <- sd_of("cab")
			if ("cab" %in% rownames(V) && is.finite(sd_c)) {
				sc <- sqrt(sab[1L, 1L] * sab[2L, 2L])
				ci <- .ae_eb_cab_score_ci(sab[1L, 2L], sd_c, sc^2, z)
				se[["cab"]] <- sd_c
				lo[["cab"]] <- max(ci[1L], -sc)
				hi[["cab"]] <- min(ci[2L], sc)
			}
		}
	}
	# a normal or pair noise parameter loses the effective number of mean parameters
	nu_cell <- if (isTRUE(out$symmetric)) 0.5 else 1
	n_cell <- nu_cell * sum(is.finite(out$Y))
	df <- nu_cell * .ae_default(rs$df_eff, 0)
	ve <- if ("ve" %in% names(VC)) VC[["ve"]] else NA_real_
	if (identical(out$family, "normal") && is.finite(ve) && !is.null(out$eb$ve_prior)) {
		# a product below the signal threshold (.ae_withhold_product) takes up more noise than its
		# variational posterior shows; a determined product is counted by the laplace df
		info <- out$eb$variational$information
		weak <- .ae_default(out$R, 0L) > 0L && !is.null(info) && isTRUE(info[["per_loading"]] < info[["threshold"]])
		vi <- if (weak) .ae_gauss_ve_interval(out, rs, z)
		if (!is.null(vi)) {
			se[["ve"]] <- vi[["sd"]]; lo[["ve"]] <- vi[["lower"]]; hi[["ve"]] <- vi[["upper"]]
		} else {
			nu <- n_cell + out$eb$ve_prior[["nu0"]] - df
			if (nu > 4) log_ci("ve", ve * sqrt(2 / (nu - 4)))
		}
	}
	if (identical(out$family, "poisson") && is.finite(ve) && !is.null(out$eb$ve_search)) {
		log_ci("ve", ve * .ae_pln_ve_sd_log(out$eb$ve_search$curve))
	}
	rho_sd <- if (!is.null(out$rho_estimation$sd)) out$rho_estimation$sd else out$eb$pair$sd
	if ("rho" %in% names(VC) && !is.null(out$eb$pair) && is.finite(VC[["rho"]]) && length(rho_sd)) {
		z_ci("rho", VC[["rho"]], rho_sd * sqrt(n_cell / max(n_cell - df, 1)))
	}
	list(se = se, ci = cbind(lower = lo, upper = hi))
}

# score-type interval of a covariance estimate `c` with sd `sd` at the estimate, `k` = va vb: the
# values x with (c - x)^2 <= z^2 sd0^2 (1 + x^2 / k), sd0^2 = sd^2 / (1 + c^2 / k) the variance at
# zero; the wald interval where the quadratic does not close (fewer than about z^2 effective actors)
.ae_eb_cab_score_ci <- function(c, sd, k, z) {
	wald <- c(c - z * sd, c + z * sd)
	if (!is.finite(k) || k <= 0) return(wald)
	q <- z^2 * sd^2 / (1 + c^2 / k)
	A <- 1 - q / k
	disc <- c^2 - A * (c^2 - q)
	if (!(A > 0) || !(disc >= 0)) return(wald)
	c((c - sqrt(disc)) / A, (c + sqrt(disc)) / A)
}

# interval of a normal noise variance: the residual sum of squares over its degrees of freedom,
# n - df_effects - df_product, with the inverse-gamma posterior of that variance. the variational
# posterior of a product the data barely determine understates how much noise the product takes
# up, so its df is the divergence of a spectral estimator of the double-centred partial residual
# (candes, sing-long and trzasko 2013): sum_k (1 + |m - n|) psi_k / s_k + 2 sum_k sum_{j != k}
# s_k psi_k / (s_k^2 - s_j^2), with s the singular values of the residual and psi those of the
# fitted product. an undirected fit keeps the laplace df
.ae_gauss_ve_interval <- function(out, rs, z) {
	res <- unlist(out$residuals, use.names = FALSE)
	ok <- is.finite(res)
	if (!any(ok) || isTRUE(out$symmetric) || is.null(out$a_sd) || is.null(out$b_sd)) return(NULL)
	rss <- sum(res[ok]^2)
	n_obs <- sum(ok)
	# the coefficients and the reported effects, which carry the row and column means of the
	# product: each effect counts its posterior variance over the noise variance per cell
	cells <- rowSums(is.finite(out$Y), dims = 2L)
	ve <- out$VC[["ve"]]
	df_lin <- length(out$coefficients) + (sum(out$a_sd^2 * rowSums(cells)) + sum(out$b_sd^2 * colSums(cells))) / ve
	df_prod <- 0
	R <- .ae_default(out$R, 0L)
	if (R > 0L) {
		if (!is.matrix(out$UVPM) || !all(is.finite(out$UVPM))) return(NULL)
		P <- out$UVPM
		# the period mean of the partial residual: the product is shared by the periods
		pr <- Reduce(`+`, lapply(seq_along(out$EZ), function(t) {
			y <- out$Y[, , t]
			x <- y - (out$EZ[[t]] - P)
			x[!is.finite(x)] <- P[!is.finite(x)]
			x
		})) / length(out$EZ)
		pr <- pr - rowMeans(pr) - rep(colMeans(pr), each = nrow(pr)) + mean(pr)
		sv <- svd(pr, nu = 0L, nv = 0L)$d
		psi <- svd(P, nu = 0L, nv = 0L)$d[seq_len(R)]
		for (k in seq_len(R)) {
			if (!(sv[k] > 0)) return(NULL)
			gap <- sv[k]^2 - sv[-k]^2
			if (any(gap == 0)) return(NULL)
			df_prod <- df_prod + (1 + abs(nrow(P) - ncol(P))) * psi[k] / sv[k] + 2 * sum(sv[k] * psi[k] / gap)
		}
	}
	nu <- n_obs - df_lin - df_prod
	if (!is.finite(nu) || nu <= 4) return(NULL)
	v <- rss / nu
	sd <- v * sqrt(2 / nu)
	c(sd = sd, lower = v * exp(-z * sd / v), upper = v * exp(z * sd / v), df_product = df_prod)
}

# sd of log ve from the marginal likelihood profile of a count fit: the parabola through the
# maximum and its two neighbours, or through the three points nearest the maximum
.ae_pln_ve_sd_log <- function(curve) {
	ok <- is.finite(curve$objective) & curve$ve > 0
	x <- log(curve$ve[ok]); y <- curve$objective[ok]
	if (length(x) < 3L) return(NA_real_)
	o <- order(x); x <- x[o]; y <- y[o]
	b <- min(max(which.max(y), 2L), length(x) - 1L)
	x3 <- x[(b - 1L):(b + 1L)]; y3 <- y[(b - 1L):(b + 1L)]
	cc <- ((y3[3L] - y3[2L]) / (x3[3L] - x3[2L]) - (y3[2L] - y3[1L]) / (x3[2L] - x3[1L])) / (x3[3L] - x3[1L])
	if (!is.finite(cc) || cc >= 0) return(NA_real_)
	1 / sqrt(-2 * cc)
}

# attach the standard errors of the reported quantities to an assembled EB fit
.ae_eb_attach_se <- function(out, rs, ab, longitudinal) {
	if (!is.null(rs)) {
		nm_a <- names(out$a); nm_b <- names(out$b)
		out$a_sd <- setNames(if (isTRUE(longitudinal)) rs$a_centred else rs$a, nm_a)
		out$b_sd <- setNames(if (isTRUE(longitudinal)) rs$b_centred else rs$b, nm_b)
		if (!is.null(rs$cells) && is.matrix(out$UVPM)) {
			out$UVPM_sd <- rs$cells
			# a probit product with fewer than three events per loading, where the product is weakly
			# identified: at least the sd of each cell under the variational factor moments, which keep
			# the spread of the posterior the linearised covariance misses
			fp <- out$eb$factor_posterior
			epl <- out$eb$variational$information[["per_loading"]]
			if (identical(out$family, "binary") && isTRUE(epl < 3) && !is.null(fp$Sig) && is.matrix(out$U)) {
				mf <- tryCatch(
					if (isTRUE(out$symmetric)) {
						.ae_mf_prod_var_sym(sweep(out$U, 2L, sqrt(abs(out$L)), "*"), sign(out$L), fp$Sig)
					} else .ae_mf_prod_var(out$U, out[["V"]], fp$Sig, fp$Om),
					error = function(e) NULL
				)
				if (identical(dim(mf), dim(out$UVPM_sd))) out$UVPM_sd <- pmax(out$UVPM_sd, sqrt(pmax(mf, 0)))
			}
			if (isTRUE(rs$held_cell_var > 0)) out$UVPM_sd <- sqrt(out$UVPM_sd^2 + rs$held_cell_var)
			dimnames(out$UVPM_sd) <- dimnames(out$UVPM)
			if (!is.null(out$ULUPM)) out$ULUPM_sd <- out$UVPM_sd
		}
		if (!is.null(rs$U) && is.matrix(out$U)) {
			out$U_sd <- rs$U
			out$V_sd <- if (isTRUE(out$symmetric)) out$U_sd else rs[["V"]]
			# below the signal threshold the rows keep their variational spread (fit$eb$factor_posterior)
			fp <- out$eb$factor_posterior
			if (!isTRUE(out$symmetric) && !is.null(fp$Sig) && !is.null(fp$Om) &&
					identical(dim(fp$Sig)[3L], nrow(out$U)) && identical(dim(fp$Om)[3L], nrow(out[["V"]]))) {
				row_sd <- function(A) t(sqrt(pmax(apply(A, 3L, diag), 0)))
				# with the spread from the other side where the report formed it (rs$U_weak)
				weak <- identical(dim(rs$U_weak), dim(out$U)) && identical(dim(rs$V_weak), dim(out[["V"]]))
				out <- .ae_eb_widen_weak_positions(
					out,
					out$eb$variational$information,
					if (weak) rs$U_weak else matrix(row_sd(fp$Sig), nrow(out$U)),
					if (weak) rs$V_weak else matrix(row_sd(fp$Om), nrow(out[["V"]]))
				)
			}
			dimnames(out$U_sd) <- dimnames(out$U)
			if (!is.null(out$V_sd)) dimnames(out$V_sd) <- dimnames(out[["V"]])
		}
		if (!is.null(rs$L)) {
			# the eigenvalues respond to the scale of the factor prior as the cells do (rs$cells_scale,
			# the scaled derivatives of the centred product): first order u_r' J u_r
			L_extra <- NULL
			if (is.matrix(out$U) && length(rs$cells_scale)) {
				L_extra <- vapply(seq_len(ncol(out$U)), function(r) {
					u <- out$U[, r]
					sum(vapply(rs$cells_scale, function(J) {
						if (!identical(dim(J), c(length(u), length(u)))) return(0)
						J[!is.finite(J)] <- 0
						sum(u * (J %*% u))^2
					}, 0))
				}, 0)
			}
			out <- .ae_eb_eigen_interval(out, rs$L, rs$L_shift, out$eb$factor_posterior$Sig, L_extra = L_extra)
			out <- .ae_eb_weak_symmetric(out, out$eb$variational$information, out$eb$factor_posterior$Sig, L_extra = L_extra)
			R <- length(out$L)
			zero <- abs(out$L) <= sqrt(.Machine$double.eps) * max(1, abs(out$L))
			# a column held at zero: the fit estimates its positions at zero, and its eigenvector is not
			# determined, uniform on the unit sphere of centred vectors as below the threshold
			# (.ae_eb_weak_symmetric()), whose entries have mean square 1 / n: the sd of the reported
			# zero is 1 / sqrt(n), the floor u_i^2 (1 - 2 E|c|) + 1 / n there at u_i = 0
			if (any(zero) && is.matrix(out$U)) {
				out$U[, zero] <- 0; out[["V"]][, zero] <- 0
				if (identical(dim(out$U_sd), dim(out$U))) {
					out$U_sd[, zero] <- 1 / sqrt(nrow(out$U))
					out$V_sd <- out$U_sd
				}
			}
		}
	}
	vs <- .ae_eb_vc_se(out, ab, rs)
	out$eb$vc_information <- rs$sab_info
	out$eb$df_effective <- rs$df_eff
	out$eb$df_linear <- rs$df_linear
	out$VC_se <- vs$se
	out$VC_ci <- vs$ci
	out
}
