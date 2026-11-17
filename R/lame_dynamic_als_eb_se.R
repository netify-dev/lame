# standard errors of the dynamic empirical-Bayes report: the additive effects
# (levels, paths and their time averages), the variance components, the
# hyperparameters of the paths and the cells of the latent product, all from the
# Laplace posterior the fit already factorised

# posterior covariance times several right-hand sides at once: RS the Tt period
# blocks (one column per functional), RZ the border's
.lde_solve_mat <- function(fac, RS, RZ) {
	Tt <- fac$Tt
	yh <- .lde_xx_solve(fac$Si, RS, fac$b)
	rhs <- RZ
	for (t in seq_len(Tt)) rhs <- rhs - crossprod(fac$C[[t]], yh[[t]])
	z <- fac$Szi %*% rhs
	s <- lapply(seq_len(Tt), function(t) yh[[t]] - fac$Y[[t]] %*% z)
	list(s = s, z = z)
}

# sds of the paths (q x Tt), of their presence-weighted time averages and of the
# levels of the q active effects, each centred over the actors of its margin as
# reported (within each period for the paths); `wts` holds the averaging weights
# (q x Tt) and `marg` the margin of each effect
.lde_effect_sd <- function(ctx, fac, wts, marg) {
	q <- ctx$q; Tt <- ctx$Tt
	qs <- nrow(fac$Si[[1L]]); nz <- nrow(fac$Szi)
	post <- .lde_posterior(fac, cond = FALSE)
	# centring contrasts within each margin
	Cm <- diag(q)
	for (g in unique(marg)) {
		k <- which(marg == g)
		Cm[k, k] <- Cm[k, k] - 1 / length(k)
	}
	# one column per path state (q x Tt), time average (q) and level (q)
	nf <- q * Tt + 2L * q
	RS <- lapply(seq_len(Tt), function(t) {
		m <- matrix(0, qs, nf)
		m[seq_len(q), (t - 1L) * q + seq_len(q)] <- Cm
		m[seq_len(q), q * Tt + seq_len(q)] <- t(t(Cm) * wts[, t])
		m
	})
	RZ <- matrix(0, nz, nf)
	RZ[seq_len(q), q * Tt + q + seq_len(q)] <- Cm
	w <- .lde_solve_mat(fac, RS, RZ)
	v <- colSums(RZ * w$z)
	for (t in seq_len(Tt)) v <- v + colSums(RS[[t]] * w$s[[t]])
	v <- sqrt(pmax(v, 0))
	list(
		path = matrix(v[seq_len(q * Tt)], q, Tt),
		bar = v[q * Tt + seq_len(q)],
		level = v[q * Tt + q + seq_len(q)],
		post = post
	)
}

# variance of the cells u_i' v_j of a Gaussian (u_i, v_j) with means U, V, blocks
# Suu[r, s, i], Svv[r, s, j] and cross covariances Suv[r, s, i, j] (NULL: none)
.lde_cell_var <- function(U, V, Suu, Svv, Suv = NULL) {
	R <- ncol(U)
	out <- matrix(0, nrow(U), nrow(V))
	for (r in seq_len(R)) for (s in seq_len(R)) {
		a <- Suu[r, s, ]; b <- Svv[r, s, ]
		out <- out + outer(a, V[, r] * V[, s]) + outer(U[, r] * U[, s], b) + outer(a, b)
		if (!is.null(Suv)) {
			cuv <- Suv[r, s, , ]
			cvu <- Suv[s, r, , ]
			out <- out + 2 * (U[, s] %o% V[, r]) * cuv + cuv * cvu
		}
	}
	pmax(out, 0)
}

# variance terms of the grand mean of a directed or undirected static product folded
# into the intercept: its gradient g in the factor rows (as .ae_eb_fold_grad()), g' S g
# and its covariance with each reported coefficient and with each period's intercept.
# `pred` holds what the coefficients' means given the factor rows contribute (law of
# total variance): with B the regression of the reported coefficients on the factor rows
# under the joint covariance, w_k = B_k' + g for the intercept and B_k' otherwise, and
# pred = w' S w, a quadratic form that cannot be negative
.lde_fold_var <- function(ctx, jt) {
	R <- ctx$R; na <- ctx$na
	g <- if (ctx$sym) {
		rep(2 * ctx$sgn * colMeans(jt$U) / na, each = na)
	} else {
		c(rep(colMeans(jt$V) / na, each = na), rep(colMeans(jt$U) / ctx$nb, each = ctx$nb))
	}
	Sz <- jt$fac$Szi
	ip <- ctx$q + seq_len(ctx$mth)
	Sff <- Sz[jt$ifi, jt$ifi, drop = FALSE]
	ctf <- as.vector(Sz[ip, jt$ifi, drop = FALSE] %*% g)
	vg <- sum(g * (Sff %*% g))
	A <- .lde_avg_map(ctx)
	v_path <- vapply(seq_len(ctx$Tt), function(t) vg + 2 * ctf[ctx$thidx[[t]][1L]], 0)
	v_path <- ifelse(is.finite(v_path), v_path, vg)
	# the reported coefficients' covariance with the factor rows, and their regression on them
	Cf <- A %*% Sz[ip, jt$ifi, drop = FALSE]
	B <- tryCatch(t(solve(Sff, t(Cf))), error = function(e) Cf %*% .lda_safe_pinv(Sff))
	Wt <- t(B)
	Wt[, 1L] <- Wt[, 1L] + g
	pred <- crossprod(Wt, Sff %*% Wt)
	pred <- (pred + t(pred)) / 2
	list(vg = vg, cov = as.vector(A %*% ctf), v_path = v_path, pred = if (all(is.finite(pred))) pred)
}

# a reported coefficient covariance with the folded grand mean in the intercept. a
# marginal intercept variance adds g' S g and twice the covariance with g' f; a normal
# fit's intercept carries its variance given the factors (.lde_centre_report()), to
# which the variance of the mean given the factors, folded mean included, is added (the
# law of total variance), so the sum cannot be negative; its covariances are those given
# the factors (`Vc`) plus the same term
.lde_fold_vcov <- function(V, fold, cond = FALSE, Vc = NULL) {
	if (isTRUE(cond) && !is.null(fold$pred) && !is.null(Vc)) {
		P <- fold$pred
		V[1L, 1L] <- Vc[1L, 1L] + max(P[1L, 1L], 0)
		V[1L, -1L] <- Vc[1L, -1L] + P[1L, -1L]
		V[-1L, 1L] <- Vc[-1L, 1L] + P[-1L, 1L]
		return(V)
	}
	v11 <- V[1L, 1L] + fold$vg + 2 * fold$cov[1L]
	if (!is.finite(v11) || v11 <= 0) {
		# not reached by a marginal intercept variance, which already holds the covariance
		V[1L, 1L] <- V[1L, 1L] + fold$vg
		return(V)
	}
	V[1L, 1L] <- v11
	V[1L, -1L] <- V[1L, -1L] + fold$cov[-1L]
	V[-1L, 1L] <- V[-1L, 1L] + fold$cov[-1L]
	V
}

# cell sds of the static product from the joint Laplace factorisation of
# .lde_fac_joint(): levels, coefficients and factor rows together
.lde_product_sd_fac <- function(ctx, jt) {
	R <- ctx$R; na <- ctx$na; nb <- ctx$nb
	Sf <- jt$fac$Szi[jt$ifi, jt$ifi, drop = FALSE]
	if (ctx$sym) {
		# F and V = F S share their unknowns
		sg <- ctx$sgn
		iu <- lapply(seq_len(R), function(r) (r - 1L) * na + seq_len(na))
		Suu <- array(0, c(R, R, na))
		Suv <- array(0, c(R, R, na, na))
		for (r in seq_len(R)) for (s in seq_len(R)) {
			blk <- Sf[iu[[r]], iu[[s]], drop = FALSE]
			Suu[r, s, ] <- diag(blk)
			Suv[r, s, , ] <- blk * sg[s]
		}
		Svv <- Suu
		for (r in seq_len(R)) for (s in seq_len(R)) Svv[r, s, ] <- Suu[r, s, ] * sg[r] * sg[s]
		return(sqrt(.lde_cell_var(jt$U, jt$V, Suu, Svv, Suv)))
	}
	ix <- .lde_widx(ctx)
	Suu <- array(0, c(R, R, na)); Svv <- array(0, c(R, R, nb)); Suv <- array(0, c(R, R, na, nb))
	for (r in seq_len(R)) for (s in seq_len(R)) {
		Suu[r, s, ] <- Sf[cbind(ix$u[[r]], ix$u[[s]])]
		Svv[r, s, ] <- Sf[cbind(ix$v[[r]], ix$v[[s]])]
		Suv[r, s, , ] <- Sf[ix$u[[r]], ix$v[[s]], drop = FALSE]
	}
	sqrt(.lde_cell_var(jt$U, jt$V, Suu, Svv, Suv))
}

# sds of the reported static positions from the joint Laplace factorisation of
# .lde_fac_joint(), formed as the static route forms them (.ae_eb_report_se()): directed
# and bipartite rows with the directions that leave the product unchanged projected out
# (.ae_eb_gauge_dirs(), .ae_eb_factor_cov(); the rows are reported uncentred) and turned
# by the rotation of the report (.lde_fac_report()); undirected eigenvectors and
# eigenvalues of F S F' by first-order perturbation (.ae_eb_eigen_sd()). actors without
# data keep the prior sd of a row
.lde_position_sd_fac <- function(ctx, st, jt, ev = NULL) {
	R <- ctx$R; na <- ctx$na; nb <- ctx$nb
	Sf <- jt$fac$Szi[jt$ifi, jt$ifi, drop = FALSE]
	Sf <- (Sf + t(Sf)) / 2
	prior <- st$hyp$fac
	if (ctx$sym) {
		es <- .ae_eb_eigen_sd(jt$U, ctx$sgn, Sf, diag(na))
		# the variational row covariances (cavity ones under the cavity rule) for the gamma interval
		Sig <- st$pos$qu$sig[, , , 1L, drop = TRUE]
		if (identical(ctx$fscale, "cavity") && !is.null(ev$lap)) {
			cs <- tryCatch(.lde_fac_cavity(ctx, st, ev$lap)$S$Sig, error = function(e) NULL)
			if (!is.null(cs) && all(is.finite(cs))) Sig <- cs
		}
		Sig <- array(Sig, c(R, R, na))
		return(list(U = es$U, L = es$L, L_Sig = Sig, prior_u = sqrt(rep_len(prior$tau, R))))
	}
	gd <- .ae_eb_gauge_dirs(jt$U, jt$V, rep(1, R), FALSE)
	Sc <- .ae_eb_factor_cov(Sf, gd, c(rep(na, R), rep(nb, R)), centre = FALSE)
	rp <- .lde_fac_report(ctx, st$pos)
	Qr <- rp$rotation
	turn <- function(C) sqrt(pmax(diag(crossprod(Qr, C %*% Qr)), 0))
	ix <- .lde_widx(ctx)
	side <- function(idx, n) {
		o <- matrix(0, n, R)
		for (i in seq_len(n)) {
			k <- vapply(idx, function(x) x[i], 0L)
			o[i, ] <- turn(Sc[k, k, drop = FALSE])
		}
		o
	}
	pb <- .ae_mf_prior_blocks(prior, R)
	Cp <- tryCatch(solve(rbind(cbind(pb$uu, pb$uv), cbind(t(pb$uv), pb$vv))), error = function(e) NULL)
	pu <- pv <- NULL
	if (!is.null(Cp)) {
		pu <- turn(Cp[seq_len(R), seq_len(R), drop = FALSE])
		pv <- turn(Cp[R + seq_len(R), R + seq_len(R), drop = FALSE])
	}
	# the variational row sds in the same gauge, for a weak product (.lde_widen_weak_positions()):
	# under the cavity rule the cavity covariances, as the static route reports them
	# (fit$eb$factor_posterior); the rotation depends on the means only
	rv <- rp
	S_v <- NULL
	if (identical(ctx$fscale, "cavity") && !is.null(ev$lap)) {
		cv <- tryCatch(.lde_fac_cavity(ctx, st, ev$lap), error = function(e) NULL)
		rc <- if (!is.null(cv)) tryCatch(.lde_fac_report(ctx, cv$pos), error = function(e) NULL)
		if (!is.null(rc) && all(is.finite(rc$U_sd)) && all(is.finite(rc$V_sd))) { rv <- rc; S_v <- cv$S; W_v <- cv$C }
	}
	if (is.null(S_v) && !is.null(ev$wk$w)) { S_v <- .lde_fac_state(ctx, st$pos, st$hyp$fac); W_v <- ev$wk$w }
	# with the spread each row takes from the other side, as the static route (.ae_eb_weak_row_var()),
	# in the reported gauge (the rotation is orthogonal, so the gauge projection commutes with it)
	if (!is.null(S_v) && !is.null(S_v$Sig) && !is.null(S_v$Om)) {
		turn_cov <- function(A) array(vapply(seq_len(dim(A)[3L]), function(i) crossprod(Qr, matrix(A[, , i], R, R) %*% Qr), matrix(0, R, R)), c(R, R, dim(A)[3L]))
		vw <- tryCatch(.ae_eb_weak_row_var(S_v$U %*% Qr, S_v[["V"]] %*% Qr, turn_cov(S_v$Sig), turn_cov(S_v$Om), W_v, centre = FALSE), error = function(e) NULL)
		if (length(vw) == (na + nb) * R && all(is.finite(vw))) {
			rv$U_sd <- matrix(sqrt(pmax(vw[seq_len(na * R)], 0)), na, R)
			rv$V_sd <- matrix(sqrt(pmax(vw[na * R + seq_len(nb * R)], 0)), nb, R)
		}
	}
	list(U = side(ix$u, na), V = side(ix$v, nb), prior_u = pu, prior_v = pv, U_var = rv$U_sd, V_var = rv$V_sd)
}

# variance each eigenvalue of an undirected static product takes from the factor-prior scale
# (fit$eb$factor_scale$L_scale, from the refits of .lde_fac_scale_uncertainty()), as in the static
# route; NULL without it
.lde_eigen_scale_var <- function(out) {
	ls <- out$eb$factor_scale$L_scale
	if (length(ls) == length(out$L) && all(is.finite(ls))) ls^2
}

# a column of an undirected binary or count fit keeps an undetermined sign up to this multiple of the
# product threshold (.lde_glm_product_threshold): on simulated d12 panels (n 30, 4 periods) the share of
# columns fitted with the wrong sign was 20% at 1-2 times the threshold, 2.6% above 2 and 1.2% above 2.5;
# a normal fit had none above the threshold (0 of 82)
.lde_glm_sign_margin <- 2.5

# a weak static product: the position sds of a directed or bipartite fit are at least the
# variational row sds, the rule of the static route (.ae_eb_widen_weak_positions()), read
# against the signal statistic the fit reports (.lde_withhold_product())
.lde_widen_weak_positions <- function(out) {
	vs <- out$eb$position_sd_variational
	if (is.null(vs) || is.null(out$U_sd)) return(out)
	# an undirected one: eigenvalue sign and eigenvector direction not determined, as in the static
	# route (.ae_eb_weak_symmetric(), positions and rows uncentred)
	if (isTRUE(out$symmetric)) {
		margin <- if (identical(out$family, "normal")) 1 else .lde_glm_sign_margin
		return(.ae_eb_weak_symmetric(out, out$eb$information, vs$L_Sig, centred = FALSE, L_extra = .lde_eigen_scale_var(out), margin = margin, cells = FALSE))
	}
	.ae_eb_widen_weak_positions(out, out$eb$information, vs$U, vs$V, align = TRUE)
}

# variance of the row means u_it' vbar_t (and column means ubar_t' v_jt) of the moving
# product under the mean-field moments, one column per period
.lde_rowmean_var <- function(ctx, pos) {
	Tt <- ctx$Tt; R <- ctx$R
	na <- dim(pos$qu$mean)[1L]; nb <- dim(pos$qv$mean)[1L]
	va <- matrix(0, na, Tt); vb <- matrix(0, nb, Tt)
	for (t in seq_len(Tt)) {
		U <- matrix(pos$qu$mean[, , t], na, R); V <- matrix(pos$qv$mean[, , t], nb, R)
		Su <- array(pos$qu$sig[, , , t], c(R, R, na)); Sv <- array(pos$qv$sig[, , , t], c(R, R, nb))
		ub <- colMeans(U); vbar <- colMeans(V)
		Svb <- apply(Sv, c(1L, 2L), sum) / nb^2
		Sub <- apply(Su, c(1L, 2L), sum) / na^2
		va[, t] <- vapply(seq_len(na), function(i) sum(vbar * (Su[, , i] %*% vbar)) + sum(U[i, ] * (Svb %*% U[i, ])), 0)
		vb[, t] <- vapply(seq_len(nb), function(j) sum(ub * (Sv[, , j] %*% ub)) + sum(V[j, ] * (Sub %*% V[j, ])), 0)
	}
	list(a = va, b = vb)
}

# cell sds of the moving product from the mean-field moments of every period
.lde_product_sd_dyn <- function(ctx, pos) {
	Tt <- ctx$Tt; R <- ctx$R
	na <- dim(pos$qu$mean)[1L]; nb <- dim(pos$qv$mean)[1L]
	out <- array(0, c(na, nb, Tt))
	for (t in seq_len(Tt)) {
		U <- matrix(pos$qu$mean[, , t], na, R)
		V <- matrix(pos$qv$mean[, , t], nb, R)
		Suu <- array(pos$qu$sig[, , , t], c(R, R, na))
		Svv <- array(pos$qv$sig[, , , t], c(R, R, nb))
		out[, , t] <- sqrt(.lde_cell_var(U, V, Suu, Svv))
	}
	out
}

# posterior variances of the level covariance (inverse-Wishart, or inverse-gamma
# per margin, at the expected second moments of the levels), plus the variance
# that the uncertainty of the levels adds to those moments; NA where a margin is off
.lde_sab_var <- function(ctx, lev, Cl) {
	pr <- ctx$prior; na <- ctx$na; nb <- ctx$nb
	ia <- seq_len(na)
	ig_var <- function(idx, k) {
		m <- lev[idx]; C <- Cl[idx, idx, drop = FALSE]
		nu <- length(idx) + pr$eta0
		psi <- pr$eta0 * pr$Sab0[k, k] + sum(m^2 + diag(C))
		d <- nu - 2
		v_ig <- 2 * psi^2 / (d^2 * max(d - 2, 1))
		v_s <- 2 * sum(C * C) + 4 * sum(m * (C %*% m))
		v_ig + v_s / d^2
	}
	out <- c(va = NA_real_, cab = NA_real_, vb = NA_real_)
	if (ctx$sym) {
		out[["va"]] <- ig_var(ia, 1L)
	} else if (ctx$kind == "bipartite") {
		out[["va"]] <- ig_var(ia, 1L)
		out[["vb"]] <- ig_var(na + seq_len(nb), 2L)
	} else {
		ib <- na + seq_len(nb)
		ma <- lev[ia]; mb <- lev[ib]
		Caa <- Cl[ia, ia, drop = FALSE]; Cbb <- Cl[ib, ib, drop = FALSE]; Cab <- Cl[ia, ib, drop = FALSE]
		S2 <- matrix(c(
			sum(ma^2) + sum(diag(Caa)),
			sum(ma * mb) + sum(diag(Cab)),
			sum(ma * mb) + sum(diag(Cab)),
			sum(mb^2) + sum(diag(Cbb))
		), 2L, 2L)
		psi <- pr$eta0 * pr$Sab0 + S2
		nu <- na + pr$eta0
		d <- nu - 3
		den <- max(d - 2, 1)
		v11 <- 2 * psi[1L, 1L]^2 / (d^2 * den)
		v22 <- 2 * psi[2L, 2L]^2 / (d^2 * den)
		# the hyperprior is diagonal and adds nothing to the covariance S2_ab / d, whose spread is
		# that of the levels' cross moment: the wishart variance (S2_ab^2 + S2_aa S2_bb) / n at the
		# levels' own second moments, not the inverse-wishart one at psi, which carries the
		# hyperprior's scale on the variances
		v12 <- (S2[1L, 2L]^2 + S2[1L, 1L] * S2[2L, 2L]) / (na * d^2)
		s11 <- 2 * sum(Caa * Caa) + 4 * sum(ma * (Caa %*% ma))
		s22 <- 2 * sum(Cbb * Cbb) + 4 * sum(mb * (Cbb %*% mb))
		s12 <- sum(Cab * Cab) + sum(Caa * Cbb) + sum(mb * (Caa %*% mb)) + sum(ma * (Cbb %*% ma)) +
			2 * sum(ma * (Cab %*% mb))
		out[] <- c(v11 + s11 / d^2, v12 + s12 / d^2, v22 + s22 / d^2)
	}
	if (isTRUE(ctx$fzm[1L])) out[["va"]] <- NA_real_
	if (length(ctx$fzm) > 1L && isTRUE(ctx$fzm[2L]) && !ctx$sym) out[["vb"]] <- NA_real_
	if (any(ctx$fzm)) out[["cab"]] <- NA_real_
	out
}

# hyperparameters: the EM estimates take their covariance from Louis' identity,
# the curvature of the expected complete-data log posterior less the variance of
# its score, which is linear in the expected statistics S of the states:
# (I_c - G Var(S) G')^-1 with G the derivative of the score in S

# numerical gradient and hessian of f at x
.lde_num_grad <- function(f, x, h = 1e-4) {
	vapply(seq_along(x), function(k) {
		e <- numeric(length(x)); e[k] <- h
		(f(x + e) - f(x - e)) / (2 * h)
	}, 0)
}
.lde_num_hess <- function(f, x, h = 1e-3) {
	k <- length(x)
	H <- matrix(0, k, k)
	for (i in seq_len(k)) for (j in seq_len(i)) {
		ei <- numeric(k); ej <- numeric(k); ei[i] <- h; ej[j] <- h
		H[i, j] <- H[j, i] <- (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h^2)
	}
	H
}

# covariance of the hyperparameters `th` under the log posterior lp(th, S), given
# the expected statistics S and their posterior covariance VS
.lde_louis_cov <- function(lp, th, S, VS) {
	f0 <- function(x) lp(x, S)
	Ic <- -.lde_num_hess(f0, th)
	G <- vapply(seq_along(S), function(j) {
		hs <- 1e-6 * max(1, abs(S[[j]]))
		e <- numeric(length(S)); e[j] <- hs
		(.lde_num_grad(function(x) lp(x, S + e), th) - .lde_num_grad(function(x) lp(x, S - e), th)) / (2 * hs)
	}, numeric(length(th)))
	G <- matrix(G, length(th))
	Io <- Ic - G %*% VS %*% t(G)
	Io <- (Io + t(Io)) / 2
	ok <- function(M) all(is.finite(M)) && all(eigen(M, symmetric = TRUE, only.values = TRUE)$values > 0)
	if (ok(Io)) return(solve(Io))
	if (ok(Ic)) return(solve(Ic))
	matrix(NA_real_, length(th), length(th))
}

# covariance of the quadratic forms m' A_k m of a Gaussian vector (mean m, covariance Sg)
.lde_qf_cov <- function(m, Sg, mats) {
	K <- length(mats)
	AS <- lapply(mats, function(A) A %*% Sg)
	Am <- lapply(mats, function(A) as.vector(A %*% m))
	out <- matrix(0, K, K)
	for (k in seq_len(K)) for (l in seq_len(k)) {
		out[k, l] <- out[l, k] <- 2 * sum(AS[[k]] * t(AS[[l]])) + 4 * sum(Am[[k]] * as.vector(Sg %*% Am[[l]]))
	}
	out
}

# sd of the innovation sd sqrt(s2) from a covariance in (rho, log s2)
.lde_ar_sigma_sd <- function(V, s2) {
	sqrt(s2) * sqrt(max(V[2L, 2L], 0)) / 2
}

# sds of the innovation scales of the moving coefficients (one per block); the persistence and
# innovation scale of the additive paths, the persistences of the coefficient blocks and the AR(1)
# law of moving positions are smoothing settings chosen by the fit, kept on the fit and not reported: the EM
# values of the additive and position laws are biased at a few periods (intervals around them
# covered 0.48 to 0.79 of simulated panels) and a block persistence stays near its prior; the
# intervals of sigma_beta covered 0.96
.lde_hyper_se <- function(ctx, st, post, fit) {
	out <- numeric(0)
	hyp <- st$hyp
	lay <- ctx$beta
	if (isTRUE(lay$any)) {
		bp <- ctx$betaprior
		share <- if (isTRUE(ctx$pool$sigma)) hyp$beta$share else bp$scale
		for (g in seq_len(lay$ng)) {
			b <- .lde_beta_group_stats(ctx, st, post, g)
			rw <- lay$kind_g[g] %in% c("rw1", "rw2")
			nm <- lay$groups[g]
			s2 <- hyp$beta$s2[g]
			if (rw) {
				# the scale alone, at the statistic A of the path
				lp <- function(th, S) -(b$cnt / 2 + bp$shape) * th[1L] - (S[[1L]] / 2 + share) / exp(th[1L])
				V <- .lde_louis_cov(lp, log(s2), b$A1, b$VA1)
				sd_g <- sqrt(s2) * sqrt(max(V[1L, 1L], 0)) / 2
			} else {
				# the scale jointly with the persistence, whose own error is not reported
				lp <- function(th, S) {
					r <- th[1L]
					if (abs(r) >= 1) return(NA_real_)
					a <- S[[1L]] - 2 * r * S[[2L]] + r^2 * S[[3L]] + (1 - r^2) * S[[4L]]
					-(b$cnt / 2 + bp$shape) * th[2L] - (a / 2 + share) / exp(th[2L]) + b$n_stat / 2 * log1p(-r^2) -
						(r - bp$rho_mean)^2 / (2 * bp$rho_sd^2)
				}
				V <- .lde_louis_cov(lp, c(hyp$beta$rho[g], log(s2)), b$S, b$VS)
				sd_g <- .lde_ar_sigma_sd(V, s2)
				if (!is.finite(sd_g) || !(sd_g > 0)) {
					# a persistence at its boundary has no joint curvature: the scale alone at the fitted
					# persistence, as for a random walk
					r0 <- hyp$beta$rho[g]
					V1 <- .lde_louis_cov(function(th, S) lp(c(r0, th[1L]), S), log(s2), b$S, b$VS)
					sd_g <- sqrt(s2) * sqrt(max(V1[1L, 1L], 0)) / 2
				}
			}
			# last resort, the curvature of the complete-data posterior of the log scale at its mode
			if (!is.finite(sd_g) || !(sd_g > 0)) sd_g <- sqrt(s2) / (2 * sqrt(b$cnt / 2 + bp$shape))
			out[[paste0("sigma_beta.", nm)]] <- sd_g
		}
	}
	out
}

# statistics of block g of the moving coefficients (as .lde_update_beta_hyper()):
# S = (A0, A1, A2, E0) of the deviations from the long-run means and their exact
# posterior covariance VS; for random walks A1 = A(1) (or the rw2 form) and VA1
.lde_beta_group_stats <- function(ctx, st, post, g) {
	lay <- ctx$beta; Tt <- ctx$Tt; q <- ctx$q
	K2 <- if (identical(lay$kind, "rw2")) .lame_beta_change_precision("rw2", Tt)
	members <- which(lay$grp == g)
	ns <- length(lay$pos_all(members[1L]))
	# deviations of every member's states from its long-run mean: a linear map of theta
	idx <- unlist(lapply(members, lay$pos_all))
	M <- matrix(0, length(idx), ctx$mth)
	M[cbind(seq_along(idx), idx)] <- 1
	n_stat <- 0L
	wts <- numeric(0)
	for (k in seq_along(members)) {
		i <- members[k]
		rows <- (k - 1L) * ns + seq_len(ns)
		if (isTRUE(lay$ar1[i])) { M[rows, lay$mupos[i]] <- -1; n_stat <- n_stat + 1L }
		wts <- c(wts, rep(1 / lay$lam[i], ns))
	}
	m <- as.vector(M %*% st$theta)
	Cz <- post$cov_zz[q + seq_len(ctx$mth), q + seq_len(ctx$mth), drop = FALSE]
	Sg <- M %*% Cz %*% t(M)
	blk <- function(A) {
		out <- matrix(0, length(idx), length(idx))
		for (k in seq_along(members)) {
			rows <- (k - 1L) * ns + seq_len(ns)
			out[rows, rows] <- A * wts[rows[1L]]
		}
		out
	}
	cnt <- length(members) * (if (!is.null(K2)) Tt - 1L else Tt) + n_stat
	if (!is.null(K2)) {
		Aq <- blk(K2)
		return(list(cnt = cnt, n_stat = n_stat, A1 = sum(m * (Aq %*% m)) + sum(Aq * Sg),
			VA1 = .lde_qf_cov(m, Sg, list(Aq))))
	}
	sel <- function(v) diag(v, ns)
	L <- matrix(0, ns, ns); for (t in 2:ns) L[t, t - 1L] <- L[t - 1L, t] <- 0.5
	mats <- list(blk(sel(c(0, rep(1, ns - 1L)))), blk(L), blk(sel(c(rep(1, ns - 1L), 0))), blk(sel(c(1, rep(0, ns - 1L)))))
	S <- vapply(mats, function(A) sum(m * (A %*% m)) + sum(A * Sg), 0)
	VS <- .lde_qf_cov(m, Sg, mats)
	if (lay$kind_g[g] %in% c("rw1", "rw2")) {
		a1 <- c(1, -2, 1, 0)
		return(list(cnt = cnt, n_stat = n_stat, A1 = sum(a1 * S), VA1 = matrix(sum(a1 * (VS %*% a1)), 1L, 1L)))
	}
	list(cnt = cnt, n_stat = n_stat, S = S, VS = VS)
}

# variance of a noise variance ve (normal residual, count noise) by Louis' identity
# in log ve: complete-data information (n + nu0) / 2 over `n` cells, less the
# variance of the score, v_sum / (4 ve^2) with v_sum the variance of the sum of squares
.lde_ve_var <- function(ve, n, nu0, v_sum) {
	ic <- (n + nu0) / 2
	io <- ic - v_sum / (4 * ve^2)
	io <- if (is.finite(io) && io > 0) io else ic
	ve^2 / io
}

# attach the standard errors to the fitted object; a failure leaves the fit as it
# is and notes what is missing
.lde_attach_se <- function(out, ctx, st, ev, fit, prep, fac_joint = NULL) {
	res <- tryCatch(.lde_report_se(out, ctx, st, ev, fit, prep, fac_joint), error = function(e) e)
	if (inherits(res, "error")) {
		.ae_note_failure("dynamic_se", "the standard errors of the effects, variance components and product", res)
		return(out)
	}
	out$a_sd <- res$a_sd; out$b_sd <- res$b_sd
	out$a_dynamic_sd <- res$a_dynamic_sd; out$b_dynamic_sd <- res$b_dynamic_sd
	out$a_level_sd <- res$a_level_sd; out$b_level_sd <- res$b_level_sd
	if (!is.null(res$UVPM_sd)) out$UVPM_sd <- res$UVPM_sd
	if (!is.null(res$U_sd) && is.null(out$U_sd)) { out$U_sd <- res$U_sd; out$V_sd <- res$V_sd }
	# eigenvalues of an undirected static product, reported uncentred: the mean direction of the
	# rows is set by the prior, along which an eigenvalue is quadratic in the rows and the
	# linearised sd misses its spread, so every column takes the static route's gamma interval
	# (.ae_eb_eigen_interval(), the linearised sd where the variational covariances are missing)
	if (!is.null(res$L_sd) && is.null(out$L_sd)) {
		out <- .ae_eb_eigen_interval(out, res$L_sd, NULL, res$L_Sig, centred = FALSE, every = TRUE, L_extra = .lde_eigen_scale_var(out))
	}
	if (!is.null(res$pos_var)) out$eb$position_sd_variational <- res$pos_var
	out$eb$vc_se <- res$vc_se
	out$eb$hyper_se <- res$hyper_se
	out
}

.lde_report_se <- function(out, ctx, st, ev, fit, prep, fac_joint) {
	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt
	na <- ctx$na; nb <- ctx$nb; q <- ctx$q
	fac <- if (!is.null(fac_joint)) fac_joint$fac else ev$fac
	# time-averaging weights of the reported a and b: the periods an actor is present
	pres <- function(P, idx, n) {
		if (is.null(P) || !all(dim(P) == c(n, Tt))) return(matrix(1 / Tt, length(idx), Tt))
		w <- matrix(as.numeric(P[idx, , drop = FALSE]), length(idx), Tt)
		w[rowSums(w) == 0, ] <- 1
		w / rowSums(w)
	}
	wa <- pres(prep$row_presence, ctx$ia, nr)
	wts <- if (ctx$sym) wa else rbind(wa, pres(prep$col_presence, ctx$ib, nc))
	marg <- if (ctx$sym) rep(1L, na) else c(rep(1L, na), rep(2L, nb))
	es <- .lde_effect_sd(ctx, fac, wts, marg)
	# moving positions: the effects are conditional on them, and the row means of the
	# product u_i' vbar_t trade against a_i; their variance is added
	if (identical(ctx$fmode, "dynamic") && !is.null(st$pos)) {
		vr <- .lde_rowmean_var(ctx, st$pos)
		add <- if (ctx$sym) vr$a else rbind(vr$a, vr$b)
		es$path <- sqrt(es$path^2 + add)
		es$bar <- sqrt(es$bar^2 + rowSums(add * wts))
		es$level <- sqrt(es$level^2 + rowMeans(add))
	}
	sab <- st$hyp$sab
	# actors without data keep their prior sd
	emb <- function(v, idx, n, prior_sd, names) {
		o <- rep(prior_sd, n); o[idx] <- v; stats::setNames(o, names)
	}
	embm <- function(M, idx, n, prior_sd, names) {
		o <- matrix(prior_sd, n, Tt, dimnames = list(names, prep$time_names)); o[idx, ] <- M; o
	}
	ra <- seq_len(na); rb <- if (ctx$sym) ra else na + seq_len(nb)
	pa <- sqrt(max(sab[1L, 1L], 0)); pb <- sqrt(max(sab[2L, 2L], 0))
	res <- list(
		a_sd = emb(es$bar[ra], ctx$ia, nr, pa, prep$row_names),
		b_sd = emb(es$bar[rb], ctx$ib, nc, pb, prep$col_names),
		a_dynamic_sd = embm(es$path[ra, , drop = FALSE], ctx$ia, nr, pa, prep$row_names),
		b_dynamic_sd = embm(es$path[rb, , drop = FALSE], ctx$ib, nc, pb, prep$col_names),
		a_level_sd = emb(es$level[ra], ctx$ia, nr, pa, prep$row_names),
		b_level_sd = emb(es$level[rb], ctx$ib, nc, pb, prep$col_names)
	)
	if (!ctx$sym && length(ctx$fzm) > 1L && isTRUE(ctx$fzm[2L])) { res$b_sd[] <- 0; res$b_dynamic_sd[] <- 0; res$b_level_sd[] <- 0 }
	if (isTRUE(ctx$fzm[1L])) { res$a_sd[] <- 0; res$a_dynamic_sd[] <- 0; res$a_level_sd[] <- 0 }
	# variance components
	Cl <- fac$Szi[seq_len(q), seq_len(q), drop = FALSE]
	vs <- .lde_sab_var(ctx, st$lev, Cl)
	vc_se <- c(va = NA_real_, cab = NA_real_, vb = NA_real_, rho = NA_real_, ve = NA_real_)
	vc_se[c("va", "cab", "vb")] <- sqrt(pmax(vs[c("va", "cab", "vb")], 0))
	if (ctx$sym) vc_se[["vb"]] <- NA_real_
	pr_sd <- out$eb$pair$sd
	if (!is.null(pr_sd) && is.finite(pr_sd)) vc_se[["rho"]] <- pr_sd
	if (isTRUE(ctx$nrm) || isTRUE(ctx$pln)) vc_se[["ve"]] <- sqrt(.lde_ve_vc(ctx, st, ev))
	res$vc_se <- vc_se
	res$hyper_se <- .lde_hyper_se(ctx, st, es$post, fit)
	# the latent product
	if (ctx$R > 0L && identical(ctx$fmode, "fac") && !is.null(fac_joint)) {
		sdc <- .lde_product_sd_fac(ctx, fac_joint)
		hv <- .lde_held_cell_var(ctx, fac_joint, ev)
		if (hv > 0) sdc <- sqrt(sdc^2 + hv)
		P <- matrix(0, nr, nc, dimnames = list(prep$row_names, prep$col_names))
		P[ctx$ia, ctx$ib] <- sdc
		res$UVPM_sd <- P
		# the positions in the gauge they are reported in, on the full actor set
		if (is.matrix(out$U) && is.matrix(out[["V"]])) {
			ps <- .lde_position_sd_fac(ctx, st, fac_joint, ev)
			R <- ctx$R
			embp <- function(M, Pout, idx, prior_sd) {
				o <- matrix(0, nrow(Pout), ncol(Pout), dimnames = dimnames(Pout))
				if (!is.null(prior_sd)) o[, seq_len(R)] <- rep(prior_sd, each = nrow(o))
				o[idx, seq_len(R)] <- M
				o
			}
			res$U_sd <- embp(ps$U, out$U, ctx$ia, ps$prior_u)
			res$V_sd <- if (ctx$sym) res$U_sd else embp(ps[["V"]], out[["V"]], ctx$ib, ps$prior_v)
			if (!ctx$sym) {
				res$pos_var <- list(U = embp(ps$U_var, out$U, ctx$ia, ps$prior_u), V = embp(ps$V_var, out[["V"]], ctx$ib, ps$prior_v))
			} else if (!is.null(ps$L_Sig)) res$pos_var <- list(L_Sig = ps$L_Sig)
			if (!is.null(ps$L)) { res$L_sd <- stats::setNames(ps$L, names(out$L)); res$L_Sig <- ps$L_Sig }
		}
	} else if (ctx$R > 0L && identical(ctx$fmode, "dynamic") && !is.null(st$pos)) {
		sdc <- .lde_product_sd_dyn(ctx, st$pos)
		P <- array(0, c(nr, nc, Tt), dimnames = list(prep$row_names, prep$col_names, prep$time_names))
		P[ctx$ia, ctx$ib, ] <- sdc
		res$UVPM_sd <- P
	}
	res
}

# sampling variance of the noise variance of a count or normal panel: the observed
# information by the Louis identity, and for counts at least the supplemented-EM value,
# whose rate is that of ve's own EM image (.lde_ve_rate_own()); the rate read from the
# tail of the iteration carries the slowest hyperparameter's and serves only without it
.lde_ve_vc <- function(ctx, st, ev, sens = NULL) {
	vp <- ctx$prior$ve_prior
	obs <- ctx$obs
	vcell <- if (!is.null(ev$lap$vcell)) ev$lap$vcell else array(0, dim(obs))
	ve <- st$hyp$ve
	n_obs <- ctx$nu * sum(obs)
	v_sum <- if (isTRUE(ctx$nrm)) {
		r <- ctx$y - ev$eta
		sum((2 * vcell^2 + 4 * r^2 * vcell)[obs], na.rm = TRUE)
	} else .lde_pln_noise_var_sum(ctx, ev, ve)
	v_ve <- .lde_ve_var(ve, n_obs, vp[["nu0"]], ctx$nu^2 * v_sum)
	# supplemented EM: the rate of ve's own EM image is its fraction of missing
	# information, that of the predictor included; the tail rate serves without it
	r_ve <- st$ve_rate
	if (isTRUE(ctx$pln)) {
		if (is.null(sens)) sens <- tryCatch(.lde_ve_sens(ctx, st, ev), error = function(e) NULL)
		r_own <- tryCatch(.lde_ve_rate_own(ctx, st, ev, sens), error = function(e) NA_real_)
		if (isTRUE(is.finite(r_own) && r_own > 0 && r_own < 1)) r_ve <- r_own
	}
	if (isTRUE(ctx$pln) && isTRUE(is.finite(r_ve)) && r_ve > 0 && r_ve < 1) {
		v_ve <- max(v_ve, ve^2 / ((n_obs + vp[["nu0"]]) / 2) / (1 - r_ve))
	}
	v_ve
}

# the coefficients of a count panel are conditional on ve, while the counts pin the mean
# of the latent predictor plus half the noise variance: an error in ve moves the intercept
# by about minus half of it. its variance enters by the delta method, g g' var(ve), with g
# the derivative of the reported coefficients along the mode, H^-1 times the change of the
# stacked cell scores with ve. returns the derivatives of the unknowns (z), of the reported
# coefficients (g) and of the coefficient paths (g_path), or NULL
.lde_ve_sens <- function(ctx, st, ev, moving_icpt = FALSE) {
	ve <- st$hyp$ve
	if (!isTRUE(ctx$pln) || is.null(ev$fac) || is.null(ve) || !is.finite(ve) || ve <= 0) return(NULL)
	Tt <- ctx$Tt; q <- ctx$q; m <- ctx$mth; na <- ctx$na; nb <- ctx$nb
	score_at <- function(v) {
		hyp <- st$hyp
		hyp$ve <- v
		wk <- if (!is.null(st$Vpos)) .lde_work_noise(ctx, ev$eta, st$Vpos, hyp, 1L) else .lde_work_plain(ctx, ev$eta, hyp)
		wk$score
	}
	h <- 1e-3 * ve
	ds <- (score_at(ve + h) - score_at(ve - h)) / (2 * h)
	ds[!ctx$obs | !is.finite(ds)] <- 0
	rs <- vector("list", Tt)
	rth <- numeric(m)
	for (t in seq_len(Tt)) {
		dt <- matrix(ds[, , t], na, nb)
		rs[[t]] <- .lde_stack_sums(dt, ctx)
		idx <- ctx$thidx[[t]]
		rth[idx] <- rth[idx] + ctx$nu * as.vector(crossprod(ctx$Xc[[t]], as.vector(dt)))
	}
	sol <- .lde_solve(ev$fac, rs, c(numeric(q), rth))
	z <- sol$z[q + seq_len(m)]
	g <- as.vector(.lde_avg_map(ctx) %*% z)
	g_path <- .lde_coef_path(ctx, z)
	# the reported intercept carries the mean of the enabled margins' paths
	fzm <- ctx$fzm
	on_a <- !isTRUE(fzm[1L]); on_b <- !ctx$sym && !isTRUE(fzm[2L])
	cs <- numeric(q)
	if (on_a) cs[seq_len(na)] <- if (ctx$sym) 2 / na else 1 / na
	if (on_b) cs[na + seq_len(nb)] <- 1 / nb
	m_t <- vapply(seq_len(Tt), function(t) sum(cs * sol$s[[t]]), 0)
	# (only an intercept carries them: without one no coefficient takes the margins' means)
	if (identical(ctx$names[1L], "intercept")) {
		g[1L] <- g[1L] + mean(m_t)
		g_path[1L, ] <- g_path[1L, ] + if (isTRUE(moving_icpt)) m_t else mean(m_t)
	}
	if (!all(is.finite(c(g, z, g_path)))) return(NULL)
	names(g) <- ctx$names
	list(z = z, g = g, g_path = g_path, s = sol$s)
}

# fraction of missing information of ve (supplemented EM): the derivative of its EM image
# with the additive predictor following ve along the mode (`sens`, .lde_ve_sens()); the
# other hyperparameters and the positions are held. NA when it cannot be formed
.lde_ve_rate_own <- function(ctx, st, ev, sens) {
	ve <- st$hyp$ve
	if (is.null(sens) || !isTRUE(ctx$pln)) return(NA_real_)
	q <- ctx$q
	sd_st <- list(x = do.call(cbind, sens$s), theta = sens$z, O = NULL)
	d_eta <- .lde_eta(ctx, sd_st)
	vcell <- if (!is.null(ev$lap$vcell)) ev$lap$vcell else array(0, dim(ctx$obs))
	if (!is.null(st$Vpos)) vcell <- vcell + st$Vpos
	img <- function(v) {
		hyp <- st$hyp
		hyp$ve <- v
		.lde_ve_image_at(ctx, ev$eta + (v - ve) * d_eta, vcell, ev$lap$cvpair, NULL, hyp)
	}
	h <- 1e-3 * ve
	r <- (img(ve + h) - img(ve - h)) / (2 * h)
	if (is.finite(r)) r else NA_real_
}

# variance of the sum of the squared count noises over the observed cells, from
# each cell's posterior mean and variance of its noise at the fitted predictor
.lde_pln_noise_var_sum <- function(ctx, ev, ve) {
	nd <- .lde_pln_terms(ctx$y, ev$eta, ve, ctx$lgy)
	m <- nd$me; v <- pmax(nd$ee - nd$me^2, 0)
	sum((2 * v^2 + 4 * m^2 * v)[ctx$obs], na.rm = TRUE)
}

# tables of the summary

# estimate, standard error and 95 percent interval of the coefficients
.lde_coef_table <- function(est, se) {
	if (is.null(se)) return(cbind(Estimate = est))
	se <- se[names(est)]
	tab <- cbind(Estimate = est, Std.Error = se, `CI 2.5%` = est - 1.96 * se, `CI 97.5%` = est + 1.96 * se)
	rownames(tab) <- names(est)
	tab
}

# the same for every coefficient and period of the coefficient paths
.lde_path_table <- function(path, se_path) {
	if (is.null(path)) return(NULL)
	nm <- as.vector(outer(rownames(path), colnames(path), function(a, b) paste0(a, "[", b, "]")))
	est <- as.vector(path)
	if (is.null(se_path)) return(matrix(est, ncol = 1L, dimnames = list(nm, "Estimate")))
	se <- as.vector(se_path[rownames(path), colnames(path), drop = FALSE])
	out <- cbind(Estimate = est, Std.Error = se, `CI 2.5%` = est - 1.96 * se, `CI 97.5%` = est + 1.96 * se)
	rownames(out) <- nm
	out
}

# variance components with standard errors; intervals on the log scale for
# variances, on the Fisher z scale for correlations, score-type for the covariance
# (.ae_eb_cab_score_ci)
.lde_vc_table <- function(object) {
	se <- object$eb$vc_se
	if (is.null(se)) return(NULL)
	est <- object$VC
	# components the model does not estimate (no dyad pairs, a fixed latent scale,
	# a margin switched off) have no row
	nm <- names(est)[is.finite(est) & names(est) %in% names(se)[is.finite(se)]]
	if (!length(nm)) return(NULL)
	se <- se[nm]
	lo <- hi <- rep(NA_real_, length(nm))
	for (k in seq_along(nm)) {
		e <- est[[nm[k]]]; s <- se[[k]]
		if (!is.finite(s)) next
		if (nm[k] %in% c("va", "vb", "ve") && e > 0) {
			lo[k] <- e * exp(-1.96 * s / e); hi[k] <- e * exp(1.96 * s / e)
		} else if (identical(nm[k], "rho") && abs(e) < 1) {
			z <- atanh(e); sz <- s / (1 - e^2)
			lo[k] <- tanh(z - 1.96 * sz); hi[k] <- tanh(z + 1.96 * sz)
		} else if (identical(nm[k], "cab") && all(c("va", "vb") %in% names(est)) &&
				isTRUE(est[["va"]] > 0 && est[["vb"]] > 0)) {
			sc <- sqrt(est[["va"]] * est[["vb"]])
			ci <- .ae_eb_cab_score_ci(e, s, sc^2, 1.96)
			lo[k] <- max(ci[1L], -sc); hi[k] <- min(ci[2L], sc)
		} else {
			lo[k] <- e - 1.96 * s; hi[k] <- e + 1.96 * s
		}
	}
	out <- cbind(Estimate = unname(est[nm]), Std.Error = unname(se), `CI 2.5%` = lo, `CI 97.5%` = hi)
	rownames(out) <- nm
	out
}

# the innovation scales of the moving coefficients, with intervals on the log scale
.lde_hyper_table <- function(object) {
	se <- object$eb$hyper_se
	if (is.null(se) || !length(se)) return(NULL)
	est <- .lde_hyper_estimates(object)
	nm <- intersect(names(se), names(est))
	if (!length(nm)) return(NULL)
	e <- est[nm]; s <- se[nm]
	lo <- e * exp(-1.96 * s / e)
	hi <- e * exp(1.96 * s / e)
	out <- cbind(Estimate = unname(e), Std.Error = unname(s), `CI 2.5%` = lo, `CI 97.5%` = hi)
	rownames(out) <- nm
	out
}

# reported values of the hyperparameters with an error, named as in fit$eb$hyper_se
.lde_hyper_estimates <- function(object) {
	if (is.null(object$sigma_beta)) return(numeric(0))
	stats::setNames(object$sigma_beta, paste0("sigma_beta.", names(object$sigma_beta)))
}
