# dyadic correlation of the binary probit EB fit of a directed network (dcor = TRUE in ame() and lame()):
# the two latent errors of a dyad are bivariate normal with correlation rho, so a pair of outcomes has
#   P(y_ij, y_ji | eta_ij, eta_ji, rho) = Phi2(s_ij eta_ij, s_ji eta_ji; s_ij s_ji rho),  s = 2 y - 1,
# and a dyad with one missing cell has the probit probability of the other. the mean-field fit
# (R/ame_als_meanfield.R) uses the expected pair log probability, its score and 2 x 2 information by tensor
# Gauss-Hermite quadrature; rho takes the posterior mean of a quadratic model under the samplers' prior.
# compiled kernels are in src/ae_kernels.cpp; the _r functions are the R references (options(lame.als_engine = "r"))


.ae_pair_cache <- new.env(parent = emptyenv())

# 64-point Gauss-Legendre rule on (-1, 1)
.ae_pair_gl <- function() {
	hit <- .ae_pair_cache[["gl"]]
	if (!is.null(hit)) return(hit)
	n <- 64L
	i <- seq_len(n - 1L)
	J <- matrix(0, n, n)
	J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- i / sqrt(4 * i^2 - 1)
	e <- eigen(J, symmetric = TRUE)
	out <- list(x = e$values, w = 2 * e$vectors[1L, ]^2)
	.ae_pair_cache[["gl"]] <- out
	out
}

# log of the bivariate normal cdf Phi2(a, b; r). pbivnorm is inaccurate below a probability of 1e-6, so small
# probabilities are integrated in the log domain by Gauss-Legendre over (a - 14, a), a the smaller limit
.ae_pair_logphi2_r <- function(a, b, r) {
	p <- pbivnorm::pbivnorm(a, b, r)
	lp <- log(pmax(p, 1e-300))
	bad <- which(!is.finite(p) | p < 1e-6)
	if (length(bad)) {
		gl <- .ae_pair_gl()
		# the integrand's mass sits below the smaller limit, so the window is placed under it
		ab <- pmin(a[bad], b[bad]); bb <- pmax(a[bad], b[bad]); rb <- r[bad]
		w <- 14
		xs <- ab - outer(rep(1, length(ab)), w * (1 - (gl$x + 1) / 2))
		lw <- matrix(log(gl$w * w / 2), length(ab), length(gl$x), byrow = TRUE)
		lf <- stats::dnorm(xs, log = TRUE) + stats::pnorm((bb - rb * xs) / sqrt(1 - rb^2), log.p = TRUE) + lw
		mx <- apply(lf, 1L, max)
		lp[bad] <- mx + log(rowSums(exp(lf - mx)))
	}
	lp
}

# compiled kernel ae_pair_logphi2_cpp: the same scheme and cutoffs in C++
.ae_pair_logphi2 <- function(a, b, r) {
	if (identical(.ae_engine(), "r")) return(.ae_pair_logphi2_r(a, b, rep_len(r, length(a))))
	ae_pair_logphi2_cpp(a, b, r, .ae_threads())
}

# log P, P = Phi2(a, b; r), and its derivatives: gradient (La, Lb), Hessian (Laa, Lbb, Lab),
# score of r (Lr = phi2 / P) and its derivative (Lrr)
.ae_pair_logp_r <- function(a, b, r) {
	kap <- sqrt(1 - r^2)
	lP <- .ae_pair_logphi2_r(a, b, r)
	ua <- (a - r * b) / kap; ub <- (b - r * a) / kap
	la <- stats::dnorm(a, log = TRUE); lb <- stats::dnorm(b, log = TRUE)
	lna <- stats::pnorm(ua, log.p = TRUE); lnb <- stats::pnorm(ub, log.p = TRUE)
	lda <- stats::dnorm(ua, log = TRUE); ldb <- stats::dnorm(ub, log = TRUE)
	ga <- exp(la + lnb - lP); gb <- exp(lb + lna - lP)
	hab <- exp(la + ldb - log(kap) - lP)
	paa <- -a * ga - (r / kap) * exp(la + ldb - lP)
	pbb <- -b * gb - (r / kap) * exp(lb + lda - lP)
	nn <- a^2 - 2 * r * a * b + b^2
	dl2 <- (r * (1 - r^2) + a * b * (1 - r^2) - r * nn) / (1 - r^2)^2
	list(
		L = lP,
		La = ga,
		Lb = gb,
		Laa = paa - ga^2,
		Lbb = pbb - gb^2,
		Lab = hab - ga * gb,
		Lr = hab,
		Lrr = hab * dl2 - hab^2
	)
}

.ae_pair_logp <- function(a, b, r) {
	if (identical(.ae_engine(), "r")) return(.ae_pair_logp_r(a, b, rep_len(r, length(a))))
	ae_pair_logp_cpp(a, b, r, .ae_threads())
}


# expected log probability of a pair, its score psi and information W in the two predictors, and the score
# and curvature of rho, for (eta1, eta2) bivariate normal (m1, m2, v1, v2, cv), by a tensor Gauss-Hermite
# rule of K nodes a side; s1, s2 the signs 2 y - 1
.ae_pair_quad_r <- function(m1, m2, v1, v2, cv, s1, s2, rho, K = 7L, chunk = 8192L) {
	n <- length(m1)
	if (n <= chunk) return(.ae_pair_quad_block(m1, m2, v1, v2, cv, s1, s2, rho, K))
	# blocks of pairs keep the per-node temporaries in cache; the result is identical
	idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
	res <- lapply(idx, function(i) .ae_pair_quad_block(m1[i], m2[i], v1[i], v2[i], cv[i], s1[i], s2[i], rho, K))
	nm <- names(res[[1L]])
	stats::setNames(lapply(nm, function(k) unlist(lapply(res, `[[`, k), use.names = FALSE)), nm)
}

.ae_pair_quad_block <- function(m1, m2, v1, v2, cv, s1, s2, rho, K) {
	nodes <- .ae_mf_gh(K)
	n <- length(m1)
	sd1 <- sqrt(pmax(v1, 0))
	l21 <- ifelse(sd1 > 0, cv / pmax(sd1, 1e-300), 0)
	l22 <- sqrt(pmax(v2 - l21^2, 0))
	r <- s1 * s2 * rho
	ss <- s1 * s2
	ell <- psi1 <- psi2 <- w11 <- w22 <- w12 <- rs <- rc <- numeric(n)
	for (k in seq_len(K)) {
		e1 <- m1 + sd1 * nodes$z[k]
		for (l in seq_len(K)) {
			wt <- nodes$w[k] * nodes$w[l]
			e2 <- m2 + l21 * nodes$z[k] + l22 * nodes$z[l]
			d <- .ae_pair_logp_r(s1 * e1, s2 * e2, r)
			ell <- ell + wt * d$L
			psi1 <- psi1 + wt * s1 * d$La
			psi2 <- psi2 + wt * s2 * d$Lb
			w11 <- w11 - wt * d$Laa
			w22 <- w22 - wt * d$Lbb
			w12 <- w12 - wt * ss * d$Lab
			rs <- rs + wt * ss * d$Lr
			rc <- rc + wt * d$Lrr
		}
	}
	list(ell = ell, psi1 = psi1, psi2 = psi2, w11 = w11, w22 = w22, w12 = w12, rs = rs, rc = rc)
}

# the same by the compiled kernel (mode 0 of ae_pair_quad_cpp), one row per pair; `chunk` is ignored
.ae_pair_quad <- function(m1, m2, v1, v2, cv, s1, s2, rho, K = 7L, chunk = 8192L) {
	if (identical(.ae_engine(), "r")) return(.ae_pair_quad_r(m1, m2, v1, v2, cv, s1, s2, rho, K, chunk))
	ae_pair_quad_cpp(m1, m2, v1, v2, cv, s1, s2, rho, as.integer(K), 0L, FALSE, FALSE, .ae_threads())
}


# the four outcome probabilities P(s1, s2) = Phi2(s1 a, s2 b; s1 s2 r) as an n x 4 matrix, from one bivariate
# probability and the margins; any below 1e-5 is recomputed directly, since the subtraction loses it
.ae_pair_probs_r <- function(a, b, r) {
	pa <- stats::pnorm(a)
	pb <- stats::pnorm(b)
	p11 <- exp(.ae_pair_logphi2_r(a, b, r))
	P <- cbind(p11, pa - p11, pb - p11, 1 - pa - pb + p11)
	sg1 <- c(1, 1, -1, -1)
	sg2 <- c(1, -1, 1, -1)
	for (j in 2:4) {
		bad <- which(!(P[, j] >= 1e-5))
		if (length(bad)) {
			P[bad, j] <- exp(.ae_pair_logphi2_r(sg1[j] * a[bad], sg2[j] * b[bad], sg1[j] * sg2[j] * r[bad]))
		}
	}
	P
}

# per-node terms of the pair problem at (a, b) with correlation r and observed signs (s1, s2): log probability,
# score, Fisher information (I11, I22, I12) over the four outcomes, with `deriv` its derivatives in a and b and
# with `rho_deriv` those in r. probabilities are floored at 1e-150 and non-finite derivatives set to zero
.ae_pair_node_fisher_r <- function(a, b, r, s1, s2, deriv = TRUE, rho_deriv = FALSE) {
	n <- length(a)
	kap <- sqrt(1 - r^2)
	ta <- (a - r * b) / kap
	tb <- (b - r * a) / kap
	fa <- stats::dnorm(a)
	fb <- stats::dnorm(b)
	P <- .ae_pair_probs_r(a, b, r)
	Pc <- pmax(P, 1e-150)
	sg1 <- c(1, 1, -1, -1)
	sg2 <- c(1, -1, 1, -1)
	ss <- sg1 * sg2
	# A_s = Phi(s2 tb), B_s = Phi(s1 ta)
	ptb <- stats::pnorm(tb)
	ntb <- stats::pnorm(-tb)
	pta <- stats::pnorm(ta)
	nta <- stats::pnorm(-ta)
	A <- cbind(ptb, ntb, ptb, ntb)
	B <- cbind(pta, pta, nta, nta)
	F11 <- rowSums(A^2 / Pc)
	F22 <- rowSums(B^2 / Pc)
	F12 <- rowSums(rep(ss, each = n) * A * B / Pc)
	I11 <- fa^2 * F11
	I22 <- fb^2 * F22
	I12 <- fa * fb * F12
	# the observed pair: its column of the matrix
	k <- 1L + 2L * (s1 < 0) + (s2 < 0)
	ix <- cbind(seq_len(n), k)
	Po <- Pc[ix]
	psi1 <- s1 * fa * A[ix] / Po
	psi2 <- s2 * fb * B[ix] / Po
	out <- list(lo = log(Po), psi1 = psi1, psi2 = psi2, I11 = I11, I22 = I22, I12 = I12)
	if (!deriv) return(out)
	ft <- stats::dnorm(tb)
	fta <- stats::dnorm(ta)
	# derivatives of A_s, B_s, P_s in a and b: dt_b/da = -r/kap, dt_b/db = 1/kap, dt_a/da = 1/kap, dt_a/db = -r/kap,
	# dP_s/da = s1 phi(a) A_s, dP_s/db = s2 phi(b) B_s
	sg1m <- rep(sg1, each = n)
	sg2m <- rep(sg2, each = n)
	Aa <- (-r / kap * ft) * sg2m
	Ab <- (ft / kap) * sg2m
	Ba <- (fta / kap) * sg1m
	Bb <- (-r / kap * fta) * sg1m
	Pa <- (fa * sg1m) * A
	Pb <- (fb * sg2m) * B
	P2 <- Pc^2
	F11a <- rowSums(2 * A * Aa / Pc - A^2 * Pa / P2)
	F11b <- rowSums(2 * A * Ab / Pc - A^2 * Pb / P2)
	F22a <- rowSums(2 * B * Ba / Pc - B^2 * Pa / P2)
	F22b <- rowSums(2 * B * Bb / Pc - B^2 * Pb / P2)
	ssm <- rep(ss, each = n)
	F12a <- rowSums(ssm * ((Aa * B + A * Ba) / Pc - A * B * Pa / P2))
	F12b <- rowSums(ssm * ((Ab * B + A * Bb) / Pc - A * B * Pb / P2))
	out$I11a <- -2 * a * fa^2 * F11 + fa^2 * F11a
	out$I11b <- fa^2 * F11b
	out$I22a <- fb^2 * F22a
	out$I22b <- -2 * b * fb^2 * F22 + fb^2 * F22b
	out$I12a <- -a * fa * fb * F12 + fa * fb * F12a
	out$I12b <- -b * fa * fb * F12 + fa * fb * F12b
	for (f in c("I11a", "I11b", "I22a", "I22b", "I12a", "I12b")) out[[f]][!is.finite(out[[f]])] <- 0
	if (!rho_deriv) return(out)
	# the probability of an outcome has P_r = s phi2(a, b; r) (Plackett, s = s1 s2 of the outcome, phi2 the bivariate
	# normal density, phi(a) phi(tb) / kap) and the mixed derivatives P_ar = s phi(a) phi(tb) (r b - a) / kap^3,
	# P_br = s phi(b) phi(ta) (r a - b) / kap^3 (dt_b/dr = (r b - a) / kap^3, dt_a/dr = (r a - b) / kap^3)
	phi2 <- fa * ft / kap
	Pr <- matrix(ssm * phi2, n, 4L)
	Par <- matrix(ssm * (fa * ft * (r * b - a) / kap^3), n, 4L)
	Pbr <- matrix(ssm * (fb * fta * (r * a - b) / kap^3), n, 4L)
	out$I11r <- rowSums(2 * Pa * Par / Pc - Pa^2 * Pr / P2)
	out$I22r <- rowSums(2 * Pb * Pbr / Pc - Pb^2 * Pr / P2)
	out$I12r <- rowSums((Par * Pb + Pa * Pbr) / Pc - Pa * Pb * Pr / P2)
	Po_r <- Pr[ix]
	out$s0 <- Po_r / Po
	out$psi1r <- Par[ix] / Po - Pa[ix] * Po_r / Po^2
	out$psi2r <- Pbr[ix] / Po - Pb[ix] * Po_r / Po^2
	Lr <- phi2 / Po
	nn <- a^2 - 2 * r * a * b + b^2
	dl2 <- s1 * s2 * (r * (1 - r^2) + a * b * (1 - r^2) - r * nn) / (1 - r^2)^2
	out$rc0 <- Lr * dl2 - Lr^2
	for (f in c("I11r", "I22r", "I12r", "s0", "psi1r", "psi2r", "rc0")) out[[f]][!is.finite(out[[f]])] <- 0
	out
}

# expected Fisher information problem of a pair under q, (eta_1, eta_2) bivariate normal, by a tensor
# Gauss-Hermite rule of K[i] nodes a side (K may differ by pair); returns ell, the expected score and information,
# with `deriv` their derivatives in the two means and with `rho_deriv` (needs `deriv`) those in the correlation
.ae_pair_quad_fisher_r <- function(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv = FALSE, rho_deriv = FALSE) {
	n <- length(m1)
	rho_deriv <- deriv && rho_deriv
	v1 <- rep_len(v1, n); v2 <- rep_len(v2, n); cv <- rep_len(cv, n); s1 <- rep_len(s1, n); s2 <- rep_len(s2, n)
	rho <- rep_len(rho, n); K <- rep_len(as.integer(K), n)
	nm <- c(
		"ell",
		"psi1",
		"psi2",
		"I11",
		"I22",
		"I12",
		if (deriv) c("I11a", "I11b", "I22a", "I22b", "I12a", "I12b"),
		if (rho_deriv) c("s0", "rc0", "psi1r", "psi2r", "I11r", "I22r", "I12r")
	)
	out <- stats::setNames(replicate(length(nm), numeric(n), simplify = FALSE), nm)
	for (kk in sort(unique(K))) {
		ix <- which(K == kk)
		nd <- .ae_mf_gh(kk)
		sd1 <- sqrt(pmax(v1[ix], 0))
		l21 <- ifelse(sd1 > 0, cv[ix] / pmax(sd1, 1e-300), 0)
		l22 <- sqrt(pmax(v2[ix] - l21^2, 0))
		acc <- stats::setNames(replicate(length(nm), numeric(length(ix)), simplify = FALSE), nm)
		for (k in seq_len(kk)) {
			a <- m1[ix] + sd1 * nd$z[k]
			for (l in seq_len(kk)) {
				wt <- nd$w[k] * nd$w[l]
				b <- m2[ix] + l21 * nd$z[k] + l22 * nd$z[l]
				nd1 <- .ae_pair_node_fisher_r(a, b, rho[ix], s1[ix], s2[ix], deriv, rho_deriv)
				for (f in nm) acc[[f]] <- acc[[f]] + wt * nd1[[switch(f, ell = "lo", f)]]
			}
		}
		for (f in nm) out[[f]][ix] <- acc[[f]]
	}
	out
}

# the same by the compiled kernel (mode 1 of ae_pair_quad_cpp): vectorised over the pairs, K a vector, one thread
# unless options(lame.als_threads) says more
.ae_pair_quad_fisher <- function(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv = FALSE, rho_deriv = FALSE) {
	if (identical(.ae_engine(), "r")) return(
		.ae_pair_quad_fisher_r(m1, m2, v1, v2, cv, s1, s2, rho, K, deriv, rho_deriv)
	)
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
# sum cw x_ij x_ji' and the (theta, actor) blocks. cw is symmetric in (i, j) within a period
.ae_pair_cross_r <- function(cx, cw) {
	nr <- cx$nr; Tt <- cx$Tt; p <- cx$p
	cw[!is.finite(cw)] <- 0
	Csum <- rowSums(cw, dims = 2L)
	gmc <- .ae_eb_hprime_sums(cw, cx$X, nr, nr, Tt, p)
	pa <- p + 1L
	x1 <- function(t, k) {
		if (k == 1L) return(matrix(1, nr, nr))
		x <- cx$X[, , k - 1L, t]
		x[!is.finite(x)] <- 0
		x
	}
	xx <- matrix(0, pa, pa)
	Mk <- lapply(seq_len(pa), function(k) matrix(0, nr, nr))
	for (t in seq_len(Tt)) {
		ct <- cw[, , t]
		for (k in seq_len(pa)) {
			xk <- x1(t, k) * ct
			Mk[[k]] <- Mk[[k]] + xk
			for (l in k:pa) {
				v <- sum(xk * t(x1(t, l)))
				xx[k, l] <- xx[k, l] + v
				if (l != k) xx[l, k] <- xx[l, k] + v
			}
		}
	}
	list(
		Csum = Csum,
		gmc = gmc,
		xx = xx,
		Mk = Mk,
		d = rowSums(Csum),
		th_a = t(vapply(Mk, colSums, numeric(nr))),
		th_b = t(vapply(Mk, rowSums, numeric(nr)))
	)
}

# the same with the sums over time and dyads from the compiled kernel (ae_pair_cross_cpp)
.ae_pair_cross <- function(cx, cw) {
	if (identical(.ae_engine(), "r")) return(.ae_pair_cross_r(cx, cw))
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
	if (is.null(est) || !is.finite(est$rho)) return(0)
	max(min(est$rho, 0.9), -0.9)
}
