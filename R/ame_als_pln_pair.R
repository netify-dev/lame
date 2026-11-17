# dyadic correlation of a directed count network at the empirical-Bayes Poisson fit (ame_als_pln.R):
# the probability of a pair of Poisson-lognormal cells with bivariate normal noise, by adaptive
# Gauss-Hermite quadrature centred at the pair mode with the Laplace Cholesky factor as axes

# nodes per side of the rule, by the noise standard deviation
.ae_pln2_kin <- function(s2max) {
	if (s2max <= 0.1) 8L else if (s2max <= 0.3) 10L else if (s2max <= 0.6) 12L else if (s2max <= 1.2) 16L else 20L
}

# log probability L, score (L1, L2) and information (w11, w22, w12) of pairs of counts y1, y2 with means
# m1, m2 and latent covariance (s11, s22, s12), by the compiled kernel (ae_pln2_node_cpp); pairs are
# vectors, the covariance may be a scalar
.ae_pln2_node <- function(y1, y2, m1, m2, s11, s22, s12, K = .ae_pln2_kin(max(s11, s22))) {
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
