# dyadic correlation of a directed count network at the empirical-Bayes Poisson fit (ame_als_pln.R):
# the probability of a pair of Poisson-lognormal cells with bivariate normal noise, by adaptive
# Gauss-Hermite quadrature centred at the pair mode with the Laplace Cholesky factor as axes

# nodes per side of the rule, by the noise standard deviation
.ae_pln2_kin <- function(s2max) {
	if (s2max <= 0.1) 8L else if (s2max <= 0.3) 10L else if (s2max <= 0.6) 12L else if (s2max <= 1.2) 16L else 20L
}

# R reference for vectors of pairs: counts y1, y2, means m1, m2 and latent covariance (s11, s22, s12);
# returns the log probability L, the score (L1, L2) and the information (w11, w22, w12)
.ae_pln2_node_r <- function(y1, y2, m1, m2, s11, s22, s12, K = .ae_pln2_kin(max(s11, s22)), chunk = 4000L) {
	n <- length(y1)
	s11 <- rep_len(s11, n); s22 <- rep_len(s22, n); s12 <- rep_len(s12, n)
	out <- vector("list", ceiling(n / chunk))
	for (k in seq_along(out)) {
		ix <- ((k - 1L) * chunk + 1L):min(k * chunk, n)
		out[[k]] <- .ae_pln2_block(y1[ix], y2[ix], m1[ix], m2[ix], s11[ix], s22[ix], s12[ix], K)
	}
	if (!length(out)) return(
		list(L = numeric(), L1 = numeric(), L2 = numeric(), w11 = numeric(), w22 = numeric(), w12 = numeric())
	)
	if (length(out) == 1L) return(out[[1L]])
	nm <- names(out[[1L]])
	setNames(lapply(nm, function(a) unlist(lapply(out, `[[`, a), use.names = FALSE)), nm)
}

.ae_pln2_block <- function(y1, y2, m1, m2, s11, s22, s12, K) {
	n <- length(y1)
	m1 <- pmin(pmax(m1, -700), 700); m2 <- pmin(pmax(m2, -700), 700)
	s11 <- pmin(pmax(s11, 1e-8), 1e6); s22 <- pmin(pmax(s22, 1e-8), 1e6)
	# hold the correlation away from one so the covariance stays invertible
	s12 <- pmax(pmin(s12, 0.999999 * sqrt(s11 * s22)), -0.999999 * sqrt(s11 * s22))
	det <- s11 * s22 - s12 * s12
	p11 <- s22 / det; p22 <- s11 / det; p12 <- -s12 / det
	# start from the single-cell modes: e = zeta - m, exp(zeta) + zeta / s = y + m / s
	mode1 <- function(y, m, s) {
		zeta <- pmax(m, ifelse(y > 0, log(y), -Inf)); u <- y + m / s
		act <- seq_len(n)
		for (it in seq_len(100L)) {
			if (!length(act)) break
			zt <- zeta[act]; ez <- exp(zt)
			step <- (ez + zt / s[act] - u[act]) / (ez + 1 / s[act])
			zt <- zt - step; zeta[act] <- zt
			act <- act[!(abs(step) <= 1e-13 * pmax(1, abs(zt)))]
		}
		zeta - m
	}
	e1 <- mode1(y1, m1, s11); e2 <- mode1(y2, m2, s22)
	phi <- function(e1, e2) y1 * e1 - exp(m1 + e1) + y2 * e2 - exp(m2 + e2) - 0.5 * (p11 * e1 * e1 + 2 * p12 * e1 * e2 + p22 * e2 * e2)
	act <- seq_len(n)
	for (it in seq_len(100L)) {
		if (!length(act)) break
		a <- act
		X1 <- exp(m1[a] + e1[a]); X2 <- exp(m2[a] + e2[a])
		g1 <- y1[a] - X1 - (p11[a] * e1[a] + p12[a] * e2[a]); g2 <- y2[a] - X2 - (p12[a] * e1[a] + p22[a] * e2[a])
		h11 <- X1 + p11[a]; h22 <- X2 + p22[a]; h12 <- p12[a]
		dt <- h11 * h22 - h12 * h12
		d1 <- (h22 * g1 - h12 * g2) / dt; d2 <- (h11 * g2 - h12 * g1) / dt
		# damp the step until the objective does not fall
		t <- rep(1, length(a)); ok <- rep(FALSE, length(a))
		n1 <- e1[a]; n2 <- e2[a]
		for (h in 0:30) {
			todo <- which(!ok)
			if (!length(todo)) break
			c1 <- e1[a][todo] + t[todo] * d1[todo]; c2 <- e2[a][todo] + t[todo] * d2[todo]
			aa <- a[todo]
			fc <- y1[aa] * c1 - exp(m1[aa] + c1) + y2[aa] * c2 - exp(m2[aa] + c2) - 0.5 * (p11[aa] * c1 * c1 + 2 * p12[aa] * c1 * c2 + p22[aa] * c2 * c2)
			fo <- y1[aa] * e1[aa] - exp(m1[aa] + e1[aa]) + y2[aa] * e2[aa] - exp(m2[aa] + e2[aa]) - 0.5 * (p11[aa] * e1[aa] * e1[aa] + 2 * p12[aa] * e1[aa] * e2[aa] + p22[aa] * e2[aa] * e2[aa])
			good <- is.finite(fc) & fc >= fo - 1e-12 * (1 + abs(fo))
			n1[todo[good]] <- c1[good]; n2[todo[good]] <- c2[good]; ok[todo[good]] <- TRUE
			t[todo[!good]] <- t[todo[!good]] / 2
		}
		e1[a] <- n1; e2[a] <- n2
		act <- a[!(pmax(abs(t * d1), abs(t * d2)) <= 1e-12 * pmax(1, abs(n1), abs(n2)))]
	}
	X1 <- exp(m1 + e1); X2 <- exp(m2 + e2)
	h11 <- X1 + p11; h22 <- X2 + p22; h12 <- p12
	dh <- h11 * h22 - h12 * h12
	a11 <- h22 / dh; a22 <- h11 / dh; a12 <- -h12 / dh # Laplace covariance
	c11 <- sqrt(a11); c21 <- a12 / c11; c22 <- sqrt(pmax(a22 - c21 * c21, 1e-300))
	rule <- .ae_pln_rule(K)
	z <- rule$z
	zi <- rep(z, each = K); zj <- rep(z, times = K) # product grid, i slowest
	lwg <- rep(rule$lw, each = K) + rep(rule$lw, times = K)
	d1m <- outer(c11, zi) # n x K^2
	d2m <- outer(c21, zi) + outer(c22, zj)
	r <- function(t) {
		out <- expm1(t) - t - 0.5 * t * t
		sm <- abs(t) < 0.1
		if (any(sm)) { t1 <- t[sm]; out[sm] <- t1 * t1 * t1 / 6 * (1 + t1 / 4 * (1 + t1 / 5 * (1 + t1 / 6 * (1 + t1 / 7 * (1 + t1 / 8 * (1 + t1 / 9 * (1 + t1 / 10 * (1 + t1 / 11)))))))) }
		out
	}
	lw <- matrix(lwg, n, K * K, byrow = TRUE) - X1 * r(d1m) - X2 * r(d2m)
	mx <- lw[cbind(seq_len(n), max.col(lw, ties.method = "first"))]
	pw <- exp(lw - mx); tot <- rowSums(pw); pw <- pw / tot
	x1 <- expm1(d1m); x2 <- expm1(d2m)
	x1[pw == 0] <- 0; x2[pw == 0] <- 0
	E1 <- rowSums(pw * d1m); E2 <- rowSums(pw * d2m)
	q1 <- d1m - E1; q2 <- d2m - E2
	ex1 <- rowSums(pw * x1); ex2 <- rowSums(pw * x2)
	u1 <- x1 - ex1; u2 <- x2 - ex2
	# m_kj = Cov(e_k, X_j) = X_j Cov(d_k, expm1(d_j))
	M11 <- X1 * rowSums(pw * q1 * u1); M12 <- X2 * rowSums(pw * q1 * u2)
	M21 <- X1 * rowSums(pw * q2 * u1); M22 <- X2 * rowSums(pw * q2 * u2)
	Ee1 <- e1 + E1; Ee2 <- e2 + E2
	w11 <- p11 * M11 + p12 * M21; w12a <- p11 * M12 + p12 * M22
	w12b <- p12 * M11 + p22 * M21; w22 <- p12 * M12 + p22 * M22
	gs <- y1 * e1 - X1 + y2 * e2 - X2 - 0.5 * (p11 * e1 * e1 + 2 * p12 * e1 * e2 + p22 * e2 * e2)
	L <- y1 * m1 + y2 * m2 - lgamma(y1 + 1) - lgamma(y2 + 1) - 0.5 * log(det / (a11 * a22 - a12 * a12)) + gs + mx + log(tot)
	list(
		L = L,
		L1 = p11 * Ee1 + p12 * Ee2,
		L2 = p12 * Ee1 + p22 * Ee2,
		w11 = w11,
		w22 = w22,
		w12 = 0.5 * (w12a + w12b)
	)
}

# compiled version (ae_pln2_node_cpp); pairs are vectors, the covariance may be a scalar
.ae_pln2_node <- function(y1, y2, m1, m2, s11, s22, s12, K = .ae_pln2_kin(max(s11, s22))) {
	if (identical(.ae_engine(), "r")) return(.ae_pln2_node_r(y1, y2, m1, m2, s11, s22, s12, K))
	rule <- .ae_pln_rule(K)
	ae_pln2_node_cpp(
		as.numeric(y1),
		as.numeric(y2),
		lgamma(as.numeric(y1) + 1),
		lgamma(as.numeric(y2) + 1),
		as.numeric(m1),
		as.numeric(m2),
		as.numeric(s11),
		as.numeric(s22),
		as.numeric(s12),
		rule$z,
		rule$lw,
		.ae_threads()
	)
}

# dyadic correlation and pair information at a fitted state `S` with cell means `ev$eta`, variances
# `ev$V` and noise variance `ve`: rho is the posterior mean of the plug-in pair profile under the
# samplers' prior, and the information integrates the pair probability over the dyad's predictors
.ae_pln_pair_stage <- function(cx, S, ev, ve) {
	pc <- .ae_pair_context(cx$y)
	if (pc$n_pairs < 3L) return(NULL)
	cvv <- .ae_pair_cov(list(nr = cx$nr, pair = pc), S)
	V3 <- array(ev$V, dim(cx$y))
	m1 <- ev$eta[pc$one]; m2 <- ev$eta[pc$two]
	v1 <- V3[pc$one]; v2 <- V3[pc$two]
	Kp <- .ae_pln2_kin(ve)
	ll <- function(rho) sum(.ae_pln2_node(pc$y1, pc$y2, m1, m2, ve, ve, rho * ve, Kp)$L)
	# coarse grid in theta = asin(rho), where the samplers' prior has density cos(theta) (1 - rho^2)^power
	th <- seq(-1.4, 1.4, length.out = 13L)
	lg <- vapply(sin(th), ll, 0)
	rho_c <- sin(th[which.max(lg)])
	bound <- 0.99
	window <- function(r) max(min(0.05, bound - abs(r)), 0.005)
	for (it in seq_len(5L)) {
		h <- window(rho_c)
		l0 <- ll(rho_c); lp <- ll(rho_c + h); lm <- ll(rho_c - h)
		g <- (lp - lm) / (2 * h)
		curv <- (lp - 2 * l0 + lm) / h^2
		step <- if (is.finite(curv) && curv < 0) -g / curv else sign(g) * 0.1
		step <- max(min(step, 0.2), -0.2)
		rho_new <- max(min(rho_c + step, bound), -bound)
		if (abs(rho_new - rho_c) < 1e-4) break
		rho_c <- rho_new
	}
	h <- window(rho_c)
	interior <- bound - abs(rho_c) >= 0.02
	prec <- NA_real_
	if (interior) {
		l0 <- ll(rho_c); lp <- ll(rho_c + h); lm <- ll(rho_c - h)
		prec <- -(lp - 2 * l0 + lm) / h^2
	}
	post <- if (interior && is.finite(prec) && prec > 0 && 1 / sqrt(prec) < 0.3) {
		# interior peak: a Gaussian about it, under the prior
		.ae_pair_rho_post(rho_c, prec)
	} else {
		# flat profile or edge peak: the posterior on a fine grid, the profile interpolated in the angle
		fine <- seq(-1.4, 1.4, length.out = 2801L)
		lpf <- stats::splinefun(th, lg, method = "fmm")(fine)
		rf <- sin(fine)
		lw <- lpf + .lame_rho_prior_log(rf, .lame_rho_prior_power()) + log(cos(fine))
		wt <- exp(lw - max(lw)); mu <- sum(wt * rf) / sum(wt)
		list(mean = mu, sd = sqrt(sum(wt * (rf - mu)^2) / sum(wt)))
	}
	rho <- max(min(post$mean, 0.98), -0.98)
	Kw <- .ae_pln2_kin(max(v1, v2) + ve)
	nd <- .ae_pln2_node(pc$y1, pc$y2, m1, m2, v1 + ve, v2 + ve, cvv + rho * ve, Kw)
	w <- cw <- array(0, dim(cx$y))
	w[pc$one] <- pmax(nd$w11, 0); w[pc$two] <- pmax(nd$w22, 0)
	cw[pc$one] <- nd$w12; cw[pc$two] <- nd$w12
	list(
		rho = rho,
		sd = post$sd,
		peak = rho_c,
		prec = prec,
		n_pairs = pc$n_pairs,
		quad = Kw,
		w = w,
		cw = cw,
		one = pc$one,
		two = pc$two
	)
}
