# collapsed update of the dynamic additive-effect hyperparameters for every
# network mode (the one-margin target, proposal state and position-scale ridge
# draw live in R/binary_bipartite_scale_moves.R).
#
# directed unipartite fits with both margins use a stacked target: with the
# levels (l_i^a, l_i^b) ~ N(0, Sab) integrated out, the totals
# x_i = (a_i1, ..., a_iT, b_i1, ..., b_iT) are independent N(0, V2),
#   V2 = Sab (x) 11' + I2 (x) Sigma,   Sigma_ts = sigma^2 / (1 - rho^2) rho^{|t - s|},
# with n - 1 degrees of freedom when the totals are centred over the actors.
# the determinant lemma and Woodbury identity (q = 1' Sigma^{-1} 1,
# w = Sigma^{-1} 1) reduce V2 to small matrices:
#   |V2| = |Sigma|^2 |I2 + q Sab|,
#   tr(V2^{-1} X'X) = tr(Sigma^{-1} (a'a + b'b)) - sum_kl M_kl (X_k w)'(X_l w),
#   M = Sab (I2 + q Sab)^{-1}.
# the priors match the sampler's conditional draws, and the levels are drawn
# next given the totals, so the pair is one Gibbs block

# one visit of an adaptive random-walk Metropolis sampler on the coordinates th
# of `logpost`: `tries` proposals from a Gaussian whose covariance (the running
# covariance of the chain times a global scale) adapts during burn-in only and
# is frozen afterwards. `state` is the numeric vector of .lame_bip_ab_state_new()
# (NULL at the start). returns the new coordinates and state, and `ok = FALSE`
# when the current point has no finite density (nothing moves then)
.lame_ab_rw_step <- function(th, logpost, state = NULL, tries = 8L, adapt = FALSE) {
	d <- length(th)
	if (is.null(state) || length(state) != 4L + d + d * d) state <- .lame_bip_ab_state_new(th)
	st <- .lame_bip_ab_state_get(state, d)
	lp <- logpost(th)
	if (!is.finite(lp)) return(list(th = th, state = state, ok = FALSE))
	for (k in seq_len(tries)) {
		L <- tryCatch(t(chol(st$lambda * (st$S + diag(1e-8, d)))), error = function(e) diag(0.02, d))
		prop <- th + drop(L %*% stats::rnorm(d))
		lpp <- logpost(prop)
		st$att <- st$att + 1
		accept <- is.finite(lpp) && log(stats::runif(1L)) < lpp - lp
		if (accept) {
			th <- prop
			lp <- lpp
			st$acc <- st$acc + 1
		}
		if (adapt) {
			st$lambda <- min(max(exp(log(st$lambda) + (as.numeric(accept) - 0.25) / sqrt(st$att)), 1e-4), 50)
			n <- st$n
			delta <- th - st$mean
			st$S <- st$S * (n - 1) / n + tcrossprod(delta) / (n + 1)
			st$mean <- st$mean + delta / (n + 1)
			st$n <- n + 1
		}
	}
	list(th = th, state = .lame_bip_ab_state_put(st), ok = TRUE)
}

# the cross products of the stacked totals that the target needs
.lame_uni_ab_cross <- function(a_mat, b_mat) {
	waa <- crossprod(a_mat)
	wbb <- crossprod(b_mat)
	list(sum = waa + wbb, aa = waa, bb = wbb, ab = crossprod(a_mat, b_mat))
}

# log density of the stacked totals up to a constant: Sab the 2 x 2 level
# covariance, Sigma the T x T AR(1) covariance, W the cross products, df the
# number of independent rows (n, or n - 1 when the totals are centred)
.lame_uni_ab_loglik <- function(Sab, Sigma, W, df) {
	ch <- tryCatch(chol(Sigma), error = function(e) NULL)
	if (is.null(ch)) return(-Inf)
	Si <- chol2inv(ch)
	w <- rowSums(Si)
	q <- sum(w)
	det_sab <- Sab[1L, 1L] * Sab[2L, 2L] - Sab[1L, 2L]^2
	D <- 1 + q * (Sab[1L, 1L] + Sab[2L, 2L]) + q * q * det_sab
	if (!is.finite(D) || D <= 0) return(-Inf)
	m11 <- (Sab[1L, 1L] + q * det_sab) / D
	m22 <- (Sab[2L, 2L] + q * det_sab) / D
	m12 <- Sab[1L, 2L] / D
	s_aa <- sum(w * (W$aa %*% w))
	s_bb <- sum(w * (W$bb %*% w))
	s_ab <- sum(w * (W$ab %*% w))
	-2 * df * sum(log(diag(ch))) - 0.5 * df * log(D) -
		0.5 * (sum(Si * W$sum) - (m11 * s_aa + m22 * s_bb + 2 * m12 * s_ab))
}

# the settings of the update for a given fit: one-margin settings (the bipartite
# context, which also serves a symmetric network and a directed network with one
# disabled margin) plus the flags of the stacked case. `rvar` is the row switch
# (the actor switch of a symmetric network), `cvar` the column switch
.lame_ab_hyper_context <- function(a_mat, b_mat, rvar, cvar, symmetric, bip, centred, prior, stat_var_cap = Inf) {
	rv <- isTRUE(rvar)
	cv <- isTRUE(cvar) && !isTRUE(symmetric)
	cx <- .lame_bip_ab_context(a_mat, b_mat, rv, cv, centred, prior, stat_var_cap)
	cx$symmetric <- isTRUE(symmetric)
	cx$stacked <- !isTRUE(bip) && rv && cv
	if (cx$stacked) cx$S0 <- matrix(as.numeric(prior$Sab0), 2L)
	cx
}

# the free coordinates of the stacked case: log va, log vb, the atanh of the
# correlation of the levels, log sigma^2 and the atanh-scale persistence
.lame_uni_ab_unpack <- function(th, cx) {
	va <- exp(th[1L])
	vb <- exp(th[2L])
	r <- tanh(th[3L])
	cab <- r * sqrt(va * vb)
	list(
		va = va,
		vb = vb,
		r = r,
		Sab = matrix(c(va, cab, cab, vb), 2L),
		s2 = exp(th[4L]),
		rho = cx$bound * tanh(th[5L])
	)
}

# log of the stacked target in those coordinates (Jacobians included): the
# inverse Wishart density of Sab in (log va, log vb, atanh r) is proportional to
# (va vb)^(-eta / 2) (1 - r^2)^(-(eta + 1) / 2) exp(-eta tr(Sab0 Sab^-1) / 2)
.lame_uni_ab_logpost <- function(th, cx, W) {
	if (!all(is.finite(th)) || abs(th[3L]) > 15) return(-Inf)
	pa <- .lame_uni_ab_unpack(th, cx)
	if (!all(is.finite(c(pa$va, pa$vb, pa$s2)))) return(-Inf)
	one_m_r2 <- 1 / cosh(th[3L])^2
	cab <- pa$Sab[1L, 2L]
	S0 <- cx$S0
	tr_ratio <- (S0[1L, 1L] * pa$vb - 2 * S0[1L, 2L] * cab + S0[2L, 2L] * pa$va) / (pa$va * pa$vb * one_m_r2)
	lp <- -0.5 * cx$eta * (th[1L] + th[2L]) - 0.5 * (cx$eta + 1) * log(one_m_r2) - 0.5 * cx$eta * tr_ratio
	lp <- lp - cx$shape * log(pa$s2) - cx$scale / pa$s2
	z <- th[5L]
	lp <- lp - 0.5 * ((pa$rho - cx$rho_mean) / cx$rho_sd)^2 + log1p(-tanh(z)^2)
	if (is.finite(cx$cap) && pa$s2 > cx$cap * (1 - pa$rho * pa$rho)) return(-Inf)
	Sigma <- .lame_bip_ab_cov(pa$s2, pa$rho, cx$N)
	lp + .lame_uni_ab_loglik(pa$Sab, Sigma, W, cx$nA - as.integer(cx$centred))
}

# one visit of the update for any mode. Sab is the 2 x 2 level covariance of the
# sampler; the function returns it with the drawn entries updated (the diagonal
# of the margins that are on, both diagonal entries of a symmetric network, and
# the covariance of a stacked directed fit), the innovation scale, the
# persistence and the proposal state. `stat_var_cap` is the bound on the
# stationary variance at the current scale of the latent response
.lame_ab_hyper_step <- function(
	a_mat,
	b_mat,
	Sab,
	sigma,
	rho,
	cx,
	state = NULL,
	tries = 8L,
	adapt = FALSE,
	stat_var_cap = cx$cap
) {
	cx$cap <- if (is.finite(stat_var_cap) && stat_var_cap > 0) stat_var_cap else Inf
	if (!isTRUE(cx$stacked)) {
		hm <- .lame_bip_ab_collapsed_step(
			a_mat,
			b_mat,
			Sab[1L, 1L],
			Sab[2L, 2L],
			sigma,
			rho,
			cx,
			state,
			tries,
			adapt
		)
		Sab[1L, 1L] <- hm$va
		Sab[2L, 2L] <- if (isTRUE(cx$symmetric)) hm$va else hm$vb
		return(list(Sab = Sab, sigma = hm$sigma, rho = hm$rho, state = hm$state))
	}
	va <- Sab[1L, 1L]
	vb <- Sab[2L, 2L]
	r <- Sab[1L, 2L] / sqrt(va * vb)
	th <- c(
		log(va),
		log(vb),
		atanh(max(min(r, 1 - 1e-12), -1 + 1e-12)),
		log(sigma * sigma),
		atanh(max(min(rho / cx$bound, 1 - 1e-12), -1 + 1e-12))
	)
	W <- .lame_uni_ab_cross(a_mat, b_mat)
	out <- .lame_ab_rw_step(th, function(x) .lame_uni_ab_logpost(x, cx, W), state, tries, adapt)
	if (!out$ok) return(list(Sab = Sab, sigma = sigma, rho = rho, state = out$state))
	pa <- .lame_uni_ab_unpack(out$th, cx)
	list(Sab = (pa$Sab + t(pa$Sab)) / 2, sigma = sqrt(pa$s2), rho = pa$rho, state = out$state)
}
