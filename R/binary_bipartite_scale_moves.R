# moves along ridges of the bipartite sampler that the conditional sweep crosses slowly
# (the stacked unipartite target and the random-walk driver are in R/dynamic_ab_hyper_move.R)

# additive-effect ridge under dynamic_ab: va and a persistent AR(1) deviation trade off on short
# panels; with the levels integrated out the totals are N(0, V), V = va 11' + AR(1) covariance,
# so the hyperparameters take a Metropolis step and the sampler then redraws the levels

# N x N covariance of a stationary AR(1) path
.lame_bip_ab_cov <- function(sigma2, rho, N) {
	sigma2 / (1 - rho * rho) * rho^abs(outer(seq_len(N), seq_len(N), "-"))
}

# log density of the totals of one margin, up to a constant (W = a'a, df = n - 1 when centred)
.lame_bip_ab_side_loglik <- function(v, Sigma, W, df) {
	ch <- tryCatch(chol(Sigma + v), error = function(e) NULL)
	if (is.null(ch)) return(-Inf)
	-df * sum(log(diag(ch))) - 0.5 * sum(chol2inv(ch) * W)
}

# prior settings and sizes the target needs
.lame_bip_ab_context <- function(a_mat, b_mat, rvar, cvar, centred, prior, stat_var_cap = Inf) {
	cap <- if (is.finite(stat_var_cap) && stat_var_cap > 0) stat_var_cap else Inf
	list(
		N = ncol(a_mat),
		nA = nrow(a_mat),
		nB = nrow(b_mat),
		rvar = isTRUE(rvar),
		cvar = isTRUE(cvar),
		centred = isTRUE(centred),
		eta = prior$etaab,
		S0a = prior$Sab0[1L, 1L],
		S0b = prior$Sab0[2L, 2L],
		shape = prior$sigma_ab_shape %||% 2,
		scale = prior$sigma_ab_scale %||% 1,
		rho_mean = prior$rho_ab_mean %||% 0.8,
		rho_sd = prior$rho_ab_sd %||% 0.15,
		cap = cap,
		bound = 0.99
	)
}

# free coordinates: log va, log vb, log sigma^2 and atanh-scale persistence (rho = bound tanh(z))
.lame_bip_ab_unpack <- function(th, cx, va, vb) {
	k <- 0L
	if (cx$rvar) { k <- k + 1L; va <- exp(th[k]) }
	if (cx$cvar) { k <- k + 1L; vb <- exp(th[k]) }
	list(va = va, vb = vb, s2 = exp(th[k + 1L]), rho = cx$bound * tanh(th[k + 2L]))
}

# log target in those coordinates, Jacobians included
.lame_bip_ab_logpost <- function(th, cx, Wa, Wb) {
	pa <- .lame_bip_ab_unpack(th, cx, 1, 1)
	if (!all(is.finite(c(th, pa$s2)))) return(-Inf)
	lp <- 0
	if (cx$rvar) lp <- lp - 0.5 * cx$eta * log(pa$va) - cx$eta * cx$S0a / (2 * pa$va)
	if (cx$cvar) lp <- lp - 0.5 * cx$eta * log(pa$vb) - cx$eta * cx$S0b / (2 * pa$vb)
	lp <- lp - cx$shape * log(pa$s2) - cx$scale / pa$s2
	z <- th[length(th)]
	lp <- lp - 0.5 * ((pa$rho - cx$rho_mean) / cx$rho_sd)^2 + log1p(-tanh(z)^2)
	if (is.finite(cx$cap) && pa$s2 > cx$cap * (1 - pa$rho * pa$rho)) return(-Inf)
	Sigma <- .lame_bip_ab_cov(pa$s2, pa$rho, cx$N)
	df <- as.integer(cx$centred)
	if (cx$rvar) lp <- lp + .lame_bip_ab_side_loglik(pa$va, Sigma, Wa, cx$nA - df)
	if (cx$cvar) lp <- lp + .lame_bip_ab_side_loglik(pa$vb, Sigma, Wb, cx$nB - df)
	lp
}

# proposal state as one numeric vector so a checkpoint can keep it: count, global scale,
# proposals, accepts, running mean (d) and covariance (d x d)
.lame_bip_ab_state_new <- function(th) {
	d <- length(th)
	c(10, 2.38^2 / d, 0, 0, th, diag(0.05^2, d))
}
.lame_bip_ab_state_get <- function(par, d) {
	list(
		n = par[1L],
		lambda = par[2L],
		att = par[3L],
		acc = par[4L],
		mean = par[4L + seq_len(d)],
		S = matrix(par[4L + d + seq_len(d * d)], d)
	)
}
.lame_bip_ab_state_put <- function(st) {
	as.numeric(c(st$n, st$lambda, st$att, st$acc, st$mean, st$S))
}

# one visit: `tries` adaptive random-walk Metropolis proposals (.lame_ab_rw_step),
# adapting during burn-in only; returns the new values and the state
.lame_bip_ab_collapsed_step <- function(a_mat, b_mat, va, vb, sigma, rho, cx, state = NULL, tries = 8L, adapt = FALSE) {
	Wa <- if (cx$rvar) crossprod(a_mat) else NULL
	Wb <- if (cx$cvar) crossprod(b_mat) else NULL
	th <- c(
		if (cx$rvar) log(va),
		if (cx$cvar) log(vb),
		log(sigma * sigma),
		atanh(max(min(rho / cx$bound, 1 - 1e-12), -1 + 1e-12))
	)
	out <- .lame_ab_rw_step(th, function(x) .lame_bip_ab_logpost(x, cx, Wa, Wb), state, tries, adapt)
	if (!out$ok) return(list(va = va, vb = vb, sigma = sigma, rho = rho, state = out$state))
	pa <- .lame_bip_ab_unpack(out$th, cx, va, vb)
	list(va = pa$va, vb = pa$vb, sigma = sqrt(pa$s2), rho = pa$rho, state = out$state)
}

# scale ridge of the static positions: (U D1, D1^-1 G D2^-1, V D2) leaves U G V' unchanged,
# so only the priors fix the split; each column scale c has an exact generalised inverse
# Gaussian conditional in x = c^2, truncated by rejection under a bound on max |G|

# generalised inverse Gaussian draw, density proportional to x^(lambda - 1) exp(-(chi / x + psi x) / 2),
# optionally truncated to x >= lower, by ratio of uniforms on log x (log-concave for every lambda)
# returns NA when the density is degenerate or the truncation is never met
.lame_rgig <- function(lambda, chi, psi, lower = 0) {
	if (!all(is.finite(c(lambda, chi, psi))) || chi <= 0 || psi <= 0) return(NA_real_)
	m <- log((lambda + sqrt(lambda * lambda + chi * psi)) / psi)
	logf <- function(u) lambda * u - 0.5 * (chi * exp(-u) + psi * exp(u))
	dlogf <- function(u) lambda + 0.5 * chi * exp(-u) - 0.5 * psi * exp(u)
	lf_m <- logf(m)
	root <- function(u) 1 + (u - m) * dlogf(u) / 2
	s <- 1 / sqrt(0.5 * (chi * exp(-m) + psi * exp(m)))
	hi <- s
	k <- 0L
	while (root(m + hi) > 0 && (k <- k + 1L) < 60L) hi <- 2 * hi
	lo <- s
	k <- 0L
	while (root(m - lo) > 0 && (k <- k + 1L) < 60L) lo <- 2 * lo
	u_hi <- stats::uniroot(root, c(m, m + hi), tol = 1e-12)$root
	u_lo <- stats::uniroot(root, c(m - lo, m), tol = 1e-12)$root
	v_hi <- (u_hi - m) * exp(0.5 * (logf(u_hi) - lf_m))
	v_lo <- (u_lo - m) * exp(0.5 * (logf(u_lo) - lf_m))
	for (k in seq_len(5000L)) {
		v1 <- stats::runif(1L)
		v2 <- stats::runif(1L, v_lo, v_hi)
		u <- m + v2 / v1
		if (2 * log(v1) <= logf(u) - lf_m && exp(u) >= lower) return(exp(u))
	}
	NA_real_
}

# one pass over the orbit coordinates: columns of U against rows of G, then columns of V
# against columns of G; `g_max` bounds max |G|, `g_var` is the prior variance of an entry of G
.lame_bip_uv_g_scale_step <- function(U, V, G, g_max = Inf, g_var = 1) {
	nA <- nrow(U)
	nB <- nrow(V)
	RA <- ncol(U)
	RB <- ncol(V)
	drawn <- 0L
	bound <- function(g) if (is.finite(g_max) && g_max > 0) (max(abs(g)) / g_max)^2 else 0
	for (r in seq_len(RA)) {
		x <- .lame_rgig((nA - RB) / 2, sum(G[r, ]^2) / g_var, sum(U[, r]^2), bound(G[r, ]))
		if (is.finite(x) && x > 0) {
			cc <- sqrt(x)
			U[, r] <- cc * U[, r]
			G[r, ] <- G[r, ] / cc
			drawn <- drawn + 1L
		}
	}
	for (s in seq_len(RB)) {
		x <- .lame_rgig((nB - RA) / 2, sum(G[, s]^2) / g_var, sum(V[, s]^2), bound(G[, s]))
		if (is.finite(x) && x > 0) {
			cc <- sqrt(x)
			V[, s] <- cc * V[, s]
			G[, s] <- G[, s] / cc
			drawn <- drawn + 1L
		}
	}
	list(U = U, V = V, G = G, drawn = drawn)
}
