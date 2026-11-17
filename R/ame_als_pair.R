# dyadic correlation of the binary probit EB fit of a directed network (dcor = TRUE in ame() and lame()):
# the two latent errors of a dyad are bivariate normal with correlation rho, so a pair of outcomes has
#   P(y_ij, y_ji | eta_ij, eta_ji, rho) = Phi2(s_ij eta_ij, s_ji eta_ji; s_ij s_ji rho),  s = 2 y - 1,
# and a dyad with one missing cell has the probit probability of the other. the mean-field fit
# (R/ame_als_meanfield.R) uses the expected pair log probability, its score and 2 x 2 information by tensor
# Gauss-Hermite quadrature; rho takes the posterior mean of a quadratic model under the samplers' prior.
# compiled kernels are in src/ae_kernels.cpp (their R references are in the tests)

# log of the bivariate normal cdf Phi2(a, b; r) by the compiled kernel: pbivnorm's scheme, with
# probabilities below 1e-6 integrated in the log domain by Gauss-Legendre over (a - 14, a), a the smaller limit
.ae_pair_logphi2 <- function(a, b, r) {
	ae_pair_logphi2_cpp(a, b, r, .ae_threads())
}

# log P, P = Phi2(a, b; r), and its derivatives: gradient (La, Lb), Hessian (Laa, Lbb, Lab),
# score of r (Lr = phi2 / P) and its derivative (Lrr)
.ae_pair_logp <- function(a, b, r) {
	ae_pair_logp_cpp(a, b, r, .ae_threads())
}

# expected log probability of a pair, its score psi and information W in the two predictors, and the score
# and curvature of rho, for (eta1, eta2) bivariate normal (m1, m2, v1, v2, cv), by a tensor Gauss-Hermite
# rule of K nodes a side; s1, s2 the signs 2 y - 1. compiled kernel (mode 0 of ae_pair_quad_cpp), one row
# per pair; `chunk` is ignored
.ae_pair_quad <- function(m1, m2, v1, v2, cv, s1, s2, rho, K = 7L, chunk = 8192L) {
	ae_pair_quad_cpp(m1, m2, v1, v2, cv, s1, s2, rho, as.integer(K), 0L, FALSE, FALSE, .ae_threads())
}

# expected Fisher information problem of a pair under q, (eta_1, eta_2) bivariate normal, by a tensor
# Gauss-Hermite rule of K[i] nodes a side (K may differ by pair); returns ell, the expected score and information,
# with `deriv` their derivatives in the two means and with `rho_deriv` (needs `deriv`) those in the correlation.
# compiled kernel (mode 1 of ae_pair_quad_cpp): vectorised over the pairs, one thread unless
# options(lame.als_threads) says more
.ae_pair_quad_fisher <- function(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv = FALSE, rho_deriv = FALSE) {
	ae_pair_quad_cpp(m1, m2, v1, v2, cv, s1, s2, rho, as.integer(K), 1L, deriv, rho_deriv, .ae_threads())
}

# the pairs of a directed panel: the dyads (i < j) of every period with both cells observed
.ae_pair_context <- function(y) {
	n <- dim(y)[1L]; Tt <- dim(y)[3L]
	ut <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	i <- ut[, 1L]; j <- ut[, 2L]
	off <- (seq_len(Tt) - 1L) * n * n
	one <- as.vector(outer(i + (j - 1L) * n, off, "+"))
	two <- as.vector(outer(j + (i - 1L) * n, off, "+"))
	ii <- rep(i, Tt); jj <- rep(j, Tt)
	keep <- is.finite(y[one]) & is.finite(y[two])
	list(
		one = one[keep],
		two = two[keep],
		i = ii[keep],
		j = jj[keep],
		s1 = 2 * y[one[keep]] - 1,
		s2 = 2 * y[two[keep]] - 1,
		y1 = y[one[keep]],
		y2 = y[two[keep]],
		n = n,
		Tt = Tt,
		n_pairs = sum(keep)
	)
}

# the covariance of the two predictors of a dyad that the additive part carries (the factors of
# the two actors are independent under q): cov(a_i + b_j, a_j + b_i) from the joint covariance of
# the effects, cov(a_i, a_j) + cov(a_i, b_i) + cov(b_j, a_j) + cov(b_j, b_i)
.ae_pair_cov <- function(cx, S) {
	if (is.null(S$Sc)) return(rep(0, cx$pair$n_pairs))
	nr <- cx$nr
	i <- cx$pair$i; j <- cx$pair$j
	Sc <- S$Sc
	Sc[cbind(i, j)] + Sc[cbind(i, nr + i)] + Sc[cbind(nr + j, j)] + Sc[cbind(nr + j, nr + i)]
}

# expected working problem of every cell when dyads are pairs: the .ae_mf_cells arrays (a pair's log
# likelihood sits in its first cell), the cross weight cw, the vector G = W m + psi of the additive normal
# equations, and the score and curvature of rho
.ae_pair_cells <- function(cx, S, eta, V3, rho) {
	pc <- cx$pair
	yy <- cx$y
	yy[c(pc$one, pc$two)] <- NA_real_
	wk <- .ae_mf_cells(yy, eta, V3, cx$link, cx$quad)
	m1 <- eta[pc$one]; m2 <- eta[pc$two]
	cvv <- .ae_pair_cov(cx, S)
	q <- .ae_pair_quad(m1, m2, V3[pc$one], V3[pc$two], cvv, pc$s1, pc$s2, rho, cx$pair_quad)
	# cap the working-response increment at 30, as in the IRLS step
	cap <- function(psi, w) {
		delta <- psi / pmax(w, 1e-12)
		big <- abs(delta) > 30
		if (any(big)) w[big] <- abs(psi[big]) / 30
		w
	}
	w1 <- cap(q$psi1, q$w11); w2 <- cap(q$psi2, q$w22)
	z <- w <- ell <- cw <- G <- array(0, dim(cx$y))
	z[pc$one] <- m1 + q$psi1 / pmax(w1, 1e-12)
	z[pc$two] <- m2 + q$psi2 / pmax(w2, 1e-12)
	w[pc$one] <- w1; w[pc$two] <- w2
	ell[pc$one] <- q$ell
	cw[pc$one] <- q$w12; cw[pc$two] <- q$w12
	G[pc$one] <- w1 * m1 + q$w12 * m2 + q$psi1
	G[pc$two] <- w2 * m2 + q$w12 * m1 + q$psi2
	rest <- which(is.finite(yy))
	z[rest] <- wk$z[rest]; w[rest] <- wk$w[rest]; ell[rest] <- wk$ell[rest]
	G[rest] <- wk$w[rest] * wk$z[rest]
	list(
		wk = list(z = z, w = w, ell = ell),
		pair = list(cw = cw, G = G, rs = sum(q$rs), rc = sum(q$rc), rho = rho, cv = cvv)
	)
}

# cavity weights when dyads are pairs: the curvature in one predictor's mean of the log pair probability
# integrated over the Gaussian predictors (a closed-form bivariate normal probability); returns diagonal
# weights (the single-cell weight where the partner is missing) and cross weights, as arrays
.ae_pair_info_w <- function(cx, eta, V3, cv, rho) {
	pc <- cx$pair
	yy <- cx$y
	yy[c(pc$one, pc$two)] <- NA_real_
	w <- .ae_mf_info_w(yy, eta, V3, cx$link)
	cw <- array(0, dim(cx$y))
	v1 <- V3[pc$one]; v2 <- V3[pc$two]
	sg1 <- sqrt(1 + v1); sg2 <- sqrt(1 + v2)
	A <- pc$s1 * eta[pc$one] / sg1; B <- pc$s2 * eta[pc$two] / sg2
	r <- pc$s1 * pc$s2 * (rho + cv) / (sg1 * sg2)
	d <- .ae_pair_logp(A, B, r)
	w[pc$one] <- pmax(-d$Laa / sg1^2, 0)
	w[pc$two] <- pmax(-d$Lbb / sg2^2, 0)
	x12 <- -pc$s1 * pc$s2 * d$Lab / (sg1 * sg2)
	cw[pc$one] <- x12; cw[pc$two] <- x12
	list(w = w, cw = cw)
}

# the sums of the cross weights that the additive block needs: the dyad-level sums over time
# of cw, of cw [1, x] (gmc: the cross analogue of gm), the (theta, theta) block
# sum cw x_ij x_ji' and the (theta, actor) blocks, from the compiled kernel (ae_pair_cross_cpp).
# cw is symmetric in (i, j) within a period
.ae_pair_cross <- function(cx, cw) {
	nr <- cx$nr; Tt <- cx$Tt; p <- cx$p
	pa <- p + 1L
	cw[!is.finite(cw)] <- 0
	Csum <- rowSums(cw, dims = 2L)
	gmc <- .ae_eb_hprime_sums(cw, cx$X, nr, nr, Tt, p)
	cr <- ae_pair_cross_cpp(cw, cx$X, nr, Tt, p)
	Mka <- array(cr$Mk, c(nr, nr, pa))
	Mk <- lapply(seq_len(pa), function(k) Mka[, , k])
	list(
		Csum = Csum,
		gmc = gmc,
		xx = cr$xx,
		Mk = Mk,
		d = rowSums(Csum),
		th_a = t(vapply(Mk, colSums, numeric(nr))),
		th_b = t(vapply(Mk, rowSums, numeric(nr)))
	)
}

# add the cross blocks to the likelihood part of the normal equations of
# (mu, beta | actor effects | gamma): H += sum_{ordered (i, j), t} cw_ijt g_ij g_ji'
# with g_ij the gradient of the predictor of the cell in the parameters
.ae_pair_add_cross <- function(ctx, cross) {
	nr <- ctx$nr
	ia <- ctx$i_s[seq_len(nr)]; ib <- ctx$i_s[nr + seq_len(nr)]
	H <- ctx$H
	H[ctx$i_th, ctx$i_th] <- H[ctx$i_th, ctx$i_th] + cross$xx
	H[ia, ctx$i_th] <- H[ia, ctx$i_th] + t(cross$th_a)
	H[ctx$i_th, ia] <- H[ctx$i_th, ia] + cross$th_a
	H[ib, ctx$i_th] <- H[ib, ctx$i_th] + t(cross$th_b)
	H[ctx$i_th, ib] <- H[ctx$i_th, ib] + cross$th_b
	H[ia, ia] <- H[ia, ia] + cross$Csum
	H[ib, ib] <- H[ib, ib] + cross$Csum
	H[cbind(ia, ib)] <- H[cbind(ia, ib)] + cross$d
	H[cbind(ib, ia)] <- H[cbind(ib, ia)] + cross$d
	ctx$H <- H
	ctx$cross <- cross
	ctx
}

# cross blocks of the information involving the factor rows, added to `H` of .ae_eb_extend; the factors
# of two actors are independent under q, so the second moments are products of means
.ae_pair_extend <- function(H, ctx, ab, fit) {
	cross <- ctx$cross
	nr <- ctx$nr; nc <- ctx$nc; m0 <- ctx$m
	U <- fit$U; V <- fit[["V"]]
	R <- ncol(U)
	C <- cross$Csum
	iu <- function(r) m0 + (r - 1L) * nr + seq_len(nr)
	iv <- function(r) m0 + R * nr + (r - 1L) * nc + seq_len(nc)
	ia <- ctx$i_s[seq_len(nr)]; ib <- ctx$i_s[nr + seq_len(nc)]
	cr <- matrix(0, nrow(H), ncol(H))
	for (r in seq_len(R)) {
		for (k in seq_along(ctx$i_th)) {
			cr[ctx$i_th[k], iu(r)] <- as.vector(crossprod(cross$Mk[[k]], V[, r]))
			cr[ctx$i_th[k], iv(r)] <- as.vector(cross$Mk[[k]] %*% U[, r])
		}
		cr[ia, iu(r)] <- V[, r] * C
		cr[cbind(ia, iv(r))] <- as.vector(C %*% U[, r])
		cr[cbind(ib, iu(r))] <- as.vector(C %*% V[, r])
		cr[ib, iv(r)] <- U[, r] * C
		for (r2 in seq_len(R)) {
			cr[cbind(iu(r), iv(r2))] <- as.vector(C %*% (V[, r] * U[, r2]))
		}
	}
	cr <- cr + t(cr)
	for (r in seq_len(R)) for (r2 in seq_len(R)) {
		cr[iu(r), iu(r2)] <- t(t(C * V[, r2]) * V[, r])
		cr[iv(r), iv(r2)] <- t(t(C * U[, r2]) * U[, r])
	}
	H + cr
}


# posterior mean and sd of rho under a Gaussian approximation of the expected log likelihood about rho_c
# with precision `prec` and the samplers' prior (R/rho_prior.R), on a grid of theta = asin(rho) with Jacobian cos(theta)
.ae_pair_rho_post <- function(rho_c, prec, power = .lame_rho_prior_power(), ngrid = 4001L) {
	th <- (seq_len(ngrid) - 0.5) / ngrid * pi - pi / 2
	r <- sin(th)
	lw <- -0.5 * prec * (r - rho_c)^2 + .lame_rho_prior_log(r, power) + log(cos(th))
	w <- exp(lw - max(lw))
	m <- sum(w * r) / sum(w)
	list(mean = m, sd = sqrt(sum(w * (r - m)^2) / sum(w)))
}

# update of rho: the posterior mean, under the samplers' prior, of the quadratic model peaked at
# rho + score / curvature (a mode would sit on the boundary where the pairs say little); with non-positive
# curvature, a quarter step in the direction of the score
.ae_pair_rho_next <- function(pair, rho, step_max = 0.2, bound = 0.98) {
	prec <- -pair$rc
	if (!(is.finite(prec) && prec > 0)) {
		step <- sign(pair$rs) * step_max / 4
		return(max(min(rho + step, bound), -bound))
	}
	post <- .ae_pair_rho_post(rho + pair$rs / prec, prec)
	max(min(post$mean, rho + step_max, bound), rho - step_max, -bound)
}

# the posterior sd of rho under the same model
.ae_pair_rho_sd <- function(pair, rho) {
	prec <- -pair$rc
	if (!(is.finite(prec) && prec > 0)) return(NA_real_)
	.ae_pair_rho_post(rho + pair$rs / prec, prec)$sd
}

# start of rho: the pairs' maximum likelihood estimate at the marginal predictors, kept within (-0.9, 0.9)
.ae_pair_rho_start <- function(cx, S) {
	if (!is.null(cx$pair_rho0) && is.finite(cx$pair_rho0)) return(max(min(cx$pair_rho0, 0.9), -0.9))
	est <- tryCatch(.ae_binary_pair_rho(cx$y, .ae_mf_predictor(cx, S)$eta, "probit"), error = function(e) NULL)
	if (is.null(est) || !is.finite(est$rho)) {
		.ae_note_fallback("pair_rho_start")
		return(0)
	}
	max(min(est$rho, 0.9), -0.9)
}
