# static latent factors (dynamic_uv = FALSE, R > 0) of the dynamic binary EB estimator of
# lame_dynamic_als_eb.R, taken as the Gaussian variational posterior of the static estimator
# (ame_als_meanfield.R) so both routes solve one problem: each actor's factors are a one-period block of
# the positions machinery (lame_dynamic_als_eb_pos.R), with the static sweeps, prior update and gauge.
# the prior scale follows options(lame.eb_scale) (ctx$fscale):
#   "cavity" (default): updated from the rows' cavity second moments (.lde_fac_scale_update()), with
#     the sampler's priors: "joint" (directed, ame()/lame() inverse-Wishart, R/uv_prior.R), "diag"
#     (bipartite, receiver at unit scale) or "columns" (undirected, one scale per column)
#   "plugin": the variational update from the mean-field second moments, with a closed-form gauge
# an undirected network keeps its start's sign pattern (ctx$sgn) and mirrors the receiver side


# prior and its constants from the start's factors
.lde_fac_prior_init <- function(ctx, U0, V0) {
	R <- ctx$R
	suv0 <- ctx$suv0
	cavity <- identical(ctx$fscale, "cavity")
	if (ctx$sym) {
		if (cavity) return(list(kind = "columns", tau = pmax(colMeans(U0^2), 1e-3), suv0 = suv0, kappa0 = 2 * R + 2))
		return(list(kind = "scalar", tau = max(mean(colMeans(U0^2)), 0.05), suv0 = suv0, kappa0 = 2 * R + 2))
	}
	if (identical(ctx$kind, "directed")) {
		sc <- pmax(c(colMeans(U0^2), colMeans(V0^2)), 0.05)
		uvp <- ctx$uv_static
		df <- if (is.null(uvp)) .lame_uv_prior_df(R) else uvp$df
		psi <- if (is.null(uvp)) df * suv0 else uvp$psi
		pr <- list(
			kind = "joint",
			Suv = diag(sc, 2L * R),
			Lam = diag(1 / sc, 2L * R),
			psi = psi * diag(2L * R),
			suv0 = suv0,
			kappa0 = df
		)
		# the cavity rule starts at the isotropic, uncorrelated prior of the rows' mean second moment
		if (cavity) pr <- .ae_mf_prior_iso(pr, R, max((sum(U0^2) + sum(V0^2)) / (2 * nrow(U0) * R), 1e-3), 0)
		return(pr)
	}
	# bipartite cavity rule: receiver factor at the sampler's unit scale, sender columns at their
	# second moments
	if (cavity) return(
		list(kind = "diag", su = pmax(colMeans(U0^2), 1e-3), sv = rep(1, R), suv0 = suv0, kappa0 = 2 * R + 2)
	)
	list(
		kind = "diag",
		su = pmax(colMeans(U0^2), 0.05),
		sv = pmax(colMeans(V0^2), 0.05),
		suv0 = suv0,
		kappa0 = 2 * R + 2
	)
}

# the prior's scales as a vector (their relative move measures stationarity)
.lde_fac_scales <- function(prior) {
	switch(prior$kind, joint = diag(prior$Suv), scalar = prior$tau, columns = prior$tau, diag = c(prior$su, prior$sv))
}

# make a mixed or extrapolated prior consistent (symmetric covariance, precision its inverse); an
# invalid prior is an error, which callers trying a state treat as a rejection
.lde_fac_fix <- function(prior) {
	if (identical(prior$kind, "joint")) {
		prior$Suv <- (prior$Suv + t(prior$Suv)) / 2
		if (any(!is.finite(prior$Suv)) || any(eigen(prior$Suv, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
			stop("invalid factor prior")
		}
		prior$Lam <- .ae_safe_inverse(prior$Suv)
	} else if (any(!(.lde_fac_scales(prior) > 0))) {
		stop("invalid factor prior")
	}
	prior
}

# mix two priors of one kind: scales geometrically, a joint covariance linearly
.lde_fac_mix <- function(p0, p1, a) {
	out <- p1
	if (identical(p1$kind, "joint")) {
		out$Suv <- (1 - a) * p0$Suv + a * p1$Suv
		out$Lam <- .ae_safe_inverse(out$Suv)
	} else if (p1$kind %in% c("scalar", "columns")) {
		out$tau <- exp((1 - a) * log(p0$tau) + a * log(p1$tau))
	} else {
		out$su <- exp((1 - a) * log(p0$su) + a * log(p1$su))
		out$sv <- exp((1 - a) * log(p0$sv) + a * log(p1$sv))
	}
	out
}


# static estimator's sweep state from the state's factor arrays (one period)
.lde_fac_state <- function(ctx, pos, prior) {
	R <- ctx$R; na <- ctx$na; nb <- ctx$nb
	S <- list(prior = prior)
	if (ctx$sym) {
		S$F <- matrix(pos$qu$mean[, , 1L], na, R)
		S$sg <- ctx$sgn
		S$Sig <- array(pos$qu$sig[, , , 1L], c(R, R, na))
	} else {
		S$U <- matrix(pos$qu$mean[, , 1L], na, R)
		S[["V"]] <- matrix(pos$qv$mean[, , 1L], nb, R)
		S$Sig <- array(pos$qu$sig[, , , 1L], c(R, R, na))
		S$Om <- array(pos$qv$sig[, , , 1L], c(R, R, nb))
	}
	S
}

# factor arrays of the state from the static estimator's state
.lde_fac_pos <- function(ctx, S) {
	R <- ctx$R
	side <- function(M, Sg) {
		n <- nrow(M)
		list(mean = array(M, c(n, R, 1L)), sig = array(Sg, c(R, R, n, 1L)), cl = array(0, c(R, R, n, 1L)))
	}
	if (ctx$sym) {
		qu <- side(S$F, S$Sig)
		return(list(qu = qu, qv = .lde_pos_mirror(qu, ctx$sgn)))
	}
	list(qu = side(S$U, S$Sig), qv = side(S[["V"]], S$Om))
}

# start: the rank-R weighted least-squares factors of the residual of the fit without factors
# (.lde_pos_svd(..., mm = TRUE)), each shrunk by `scale`, with their implied covariances and prior.
# a much smaller start can fall into the null fixed point on weak panels, a much larger one can diverge
.lde_fac_start <- function(ctx, fit0, scale = if (is.null(ctx$fac_start)) 0.5 else ctx$fac_start, jitter = 0) {
	R <- ctx$R; na <- ctx$na
	sv <- .lde_pos_svd(ctx, fit0, scale, jitter, mm = TRUE)
	U0 <- sv$U0; V0 <- sv$V0
	if (identical(ctx$fscale, "cavity") && identical(ctx$kind, "bipartite")) {
		# the receiver factor takes the sampler's unit scale; U D and V D^-1 leave the product unchanged
		dv <- sqrt(pmax(colMeans(V0^2), 1e-8))
		U0 <- sweep(U0, 2L, dv, "*")
		V0 <- sweep(V0, 2L, dv, "/")
	}
	prior <- .lde_fac_prior_init(ctx, U0, V0)
	W0 <- fit0$ev$wk$w
	W0[!ctx$obs | !is.finite(W0)] <- 0
	C0 <- rowSums(W0, dims = 2L)
	pb <- .ae_mf_prior_blocks(prior, R)
	S <- list(prior = prior)
	if (ctx$sym) {
		Cn <- C0
		diag(Cn) <- 0
		Fs <- sweep(U0, 2L, sv$sgn, "*")
		S$F <- U0
		S$sg <- sv$sgn
		S$Sig <- .ae_mf_row_solve(Cn %*% .ae_mf_outer_rows(Fs), pb$uu, matrix(0, na, R))$Sig
	} else {
		cv <- .ae_mf_cov0(U0, V0, C0, pb$uu, pb$vv)
		S$U <- U0; S[["V"]] <- V0; S$Sig <- cv$Sig; S$Om <- cv$Om
	}
	if (ctx$sym) ctx$sgn <- sv$sgn
	pp <- .lde_fac_product(ctx, S)
	list(pos = .lde_fac_pos(ctx, S), prior = prior, sgn = sv$sgn, O = pp$O, Vpos = pp$Vpos)
}


# product of the factors (the same in every period) and its variance under q
.lde_fac_product <- function(ctx, S) {
	na <- ctx$na; nb <- ctx$nb; Tt <- ctx$Tt
	if (ctx$sym) {
		P <- tcrossprod(sweep(S$F, 2L, S$sg, "*"), S$F)
		V <- .ae_mf_prod_var_sym(S$F, S$sg, S$Sig)
	} else {
		P <- tcrossprod(S$U, S[["V"]])
		V <- .ae_mf_prod_var(S$U, S[["V"]], S$Sig, S$Om)
	}
	list(O = array(P, c(na, nb, Tt)), Vpos = array(V, c(na, nb, Tt)))
}

# static-estimator sweeps at the incoming working problem (C_ij = sum_t w_ijt, Rbar_ij the weighted
# mean working residual after the additive part); returns the factor arrays, the product, its variance
# and the objective on the ordered-cell deviance scale. the plug-in rule also applies the gauge and the
# variational prior update; the cavity rule leaves the prior to .lde_fac_scale_update() and compares
# states at their own prior (.ae_mf_factor_term_fixed())
.lde_fac_step <- function(ctx, st, wk0, eta_add, inner = 3L) {
	R <- ctx$R
	pr <- st$hyp$fac
	w <- wk0$w
	w[!ctx$obs | !is.finite(w)] <- 0
	r <- wk0$z - eta_add
	r[!ctx$obs | !is.finite(r)] <- 0
	C <- rowSums(w, dims = 2L)
	Rbar <- matrix(0, ctx$na, ctx$nb)
	up <- C > 0
	Rbar[up] <- rowSums(w * r, dims = 2L)[up] / C[up]
	S <- .lde_fac_state(ctx, st$pos, pr)
	if (ctx$sym) {
		sw <- .ae_mf_sweeps_sym(Rbar, C, S$F, S$sg, S$Sig, pr, inner)
		S$F <- sw$F; S$Sig <- sw$Sig
	} else {
		sw <- .ae_mf_sweeps(Rbar, C, S$U, S[["V"]], S$Sig, S$Om, pr, inner)
		S$U <- sw$U; S[["V"]] <- sw[["V"]]; S$Sig <- sw$Sig; S$Om <- sw$Om
	}
	if (identical(ctx$fscale, "cavity")) {
		pp <- .lde_fac_product(ctx, S)
		return(list(
			pos = .lde_fac_pos(ctx, S),
			O = pp$O,
			Vpos = pp$Vpos,
			prior = S$prior,
			obj = -2 * .ae_mf_factor_term_fixed(S, R),
			change = NA_real_
		))
	}
	cxf <- list(R = R, sym = isTRUE(ctx$sym))
	S <- .ae_mf_gauge(cxf, S)
	old <- .lde_fac_scales(S$prior)
	S$prior <- .ae_mf_prior_update(S, R)
	change <- max(abs(log(.lde_fac_scales(S$prior) / old)))
	pp <- .lde_fac_product(ctx, S)
	list(
		pos = .lde_fac_pos(ctx, S),
		O = pp$O,
		Vpos = pp$Vpos,
		prior = S$prior,
		obj = -2 * .ae_mf_factor_term(S, R),
		change = change
	)
}

# cavity-rule update of the prior scale: the root of log target(s) - log s, where the target is the
# plug-in scale from the rows' cavity second moments (.ae_mf_scale_target()), found by a few secant
# Newton steps (plain steps bounded by 10 smax in the log); solving it within each evaluation keeps
# the update a function of the state alone, as the squared extrapolation needs. the directed
# correlation takes a plain step of at most smax; returns the prior and the size of the move
.lde_fac_scale_update <- function(ctx, pos, prior, C, smax = 0.05, max_eval = 6L) {
	R <- ctx$R
	cx <- list(R = R, sym = isTRUE(ctx$sym))
	S <- .lde_fac_state(ctx, pos, prior)
	clip <- function(x, m) pmax(pmin(x, m), -m)
	direct <- identical(prior$kind, "joint")
	key <- if (isTRUE(ctx$sym)) "tau" else if (direct) "s" else "su"
	# scales as a vector of logs and the prior they give
	cur <- switch(key, tau = prior$tau, su = prior$su, s = prior$Suv[1L, 1L])
	rho0 <- if (direct) prior$Suv[1L, R + 1L] / prior$Suv[1L, 1L]
	with_scale <- function(x, rho) {
		pr <- prior
		if (key == "tau") pr$tau <- exp(x) else if (key == "su") pr$su <- exp(x) else pr <- .ae_mf_prior_iso(pr, R, exp(x), rho)
		pr
	}
	resid <- function(x, rho) {
		S2 <- S
		S2$prior <- with_scale(x, rho)
		tg <- .ae_mf_scale_target(cx, S2, C)
		list(r = log(tg[[key]]) - x, rho = if (direct) tg$rho)
	}
	x0 <- log(cur)
	rho <- rho0
	e0 <- resid(x0, rho)
	r0 <- e0$r
	x <- x0 + clip(r0, 10 * smax)
	if (direct) rho <- rho + clip(e0$rho - rho, smax)
	for (j in seq_len(max_eval - 1L)) {
		e1 <- resid(x, rho)
		r1 <- e1$r
		m <- (r1 - r0) / (x - x0)
		step <- ifelse(is.finite(m) & m < -1e-3 & m > -2, -r1 / m, r1)
		x0 <- x
		r0 <- r1
		x <- x + clip(step, 10 * smax)
		if (direct) rho <- rho + clip(e1$rho - rho, smax)
		if (max(abs(step)) < 1e-9) break
	}
	prior <- with_scale(x, rho)
	change <- max(abs(x - log(cur)), if (direct) abs(rho - rho0))
	list(prior = prior, change = change)
}


# sender and receiver factor means in the principal-axes gauge of the product (one period)
.lde_fac_report <- function(ctx, pos) {
	ctx1 <- ctx
	ctx1$Tt <- 1L
	rp <- .lde_pos_report(ctx1, pos)
	list(
		U = matrix(rp$U[, , 1L], ctx$na, ctx$R),
		V = matrix(rp$V[, , 1L], ctx$nb, ctx$R),
		U_sd = matrix(rp$U_sd[, , 1L], ctx$na, ctx$R),
		V_sd = matrix(rp$V_sd[, , 1L], ctx$nb, ctx$R),
		rotation = rp$rotation
	)
}


# information blocks of one period's factors at the cell weights Wt (layout of .lde_wblocks()), with
# E[v v'] = vhat vhat' + Om under the variational posterior; with a dyadic correlation `Ct` holds the
# cross weight of each dyad's two cells, coupling parameters across them (as .ae_pair_extend())
.lde_fac_wblocks <- function(ctx, Wt, Ut, Vt, Sig, Om, Ct = NULL) {
	b <- .lde_wblocks_dir(ctx, Wt, Ut, Vt)
	ix <- .lde_widx(ctx)
	R <- ctx$R
	for (r in seq_len(R)) for (s in seq_len(R)) {
		iu <- ix$u[[r]]; ius <- ix$u[[s]]; iv <- ix$v[[r]]; ivs <- ix$v[[s]]
		b$Kww[cbind(iu, ius)] <- b$Kww[cbind(iu, ius)] + as.vector(Wt %*% Om[r, s, ])
		b$Kww[cbind(iv, ivs)] <- b$Kww[cbind(iv, ivs)] + as.vector(crossprod(Wt, Sig[r, s, ]))
	}
	if (!is.null(Ct) && !ctx$sym) {
		na <- ctx$na; nb <- ctx$nb
		ia <- seq_len(na); ib <- na + seq_len(nb)
		for (r in seq_len(R)) {
			iu <- ix$u[[r]]; iv <- ix$v[[r]]
			# (a_i, u_j r) = c_ij v_i r, (a_i, v_i r) = sum_j c_ij u_j r, (b_j, u_j r) = sum_i c_ij v_i r,
			# (b_j, v_i r) = c_ij u_j r
			b$Kxw[ia, iu] <- b$Kxw[ia, iu] + Ct * Vt[, r]
			b$Kxw[cbind(ia, iv)] <- b$Kxw[cbind(ia, iv)] + as.vector(Ct %*% Ut[, r])
			b$Kxw[cbind(ib, iu)] <- b$Kxw[cbind(ib, iu)] + as.vector(Ct %*% Vt[, r])
			b$Kxw[ib, iv] <- b$Kxw[ib, iv] + Ct * Ut[, r]
			for (s in seq_len(R)) {
				ius <- ix$u[[s]]; ivs <- ix$v[[s]]
				# (u_i r, u_j s) = c_ij v_j r v_i s, (v_i r, v_j s) = c_ij u_j r u_i s, (u_i r, v_i s) = sum_j
				# c_ij v_j r u_j s
				b$Kww[iu, ius] <- b$Kww[iu, ius] + t(t(Ct * Vt[, s]) * Vt[, r])
				b$Kww[iv, ivs] <- b$Kww[iv, ivs] + t(t(Ct * Ut[, s]) * Ut[, r])
				cr <- as.vector(Ct %*% (Vt[, r] * Ut[, s]))
				b$Kww[cbind(iu, ivs)] <- b$Kww[cbind(iu, ivs)] + cr
				b$Kww[cbind(ivs, iu)] <- b$Kww[cbind(ivs, iu)] + cr
			}
		}
	}
	if (!ctx$sym) return(b)
	na <- ctx$na
	M <- .lde_tie(ctx)
	Kxw <- b$Kxw[seq_len(na), , drop = FALSE] + b$Kxw[na + seq_len(na), , drop = FALSE]
	list(Kxw = ctx$nu * (Kxw %*% M), Kww = ctx$nu * crossprod(M, b$Kww %*% M))
}

# factor-row prior precision in the stacked layout (vec(U), vec(V)), with the sender-receiver
# coupling of a directed network
.lde_fac_ridge <- function(ctx, prior) {
	R <- ctx$R; na <- ctx$na; nb <- ctx$nb
	if (ctx$sym) return(diag(rep_len(1 / prior$tau, R) %x% rep(1, na), na * R))
	nf <- R * (na + nb)
	H <- matrix(0, nf, nf)
	ix <- .lde_widx(ctx)
	pb <- .ae_mf_prior_blocks(prior, R)
	for (r in seq_len(R)) for (s in seq_len(R)) {
		H[cbind(ix$u[[r]], ix$u[[s]])] <- pb$uu[r, s]
		H[cbind(ix$v[[r]], ix$v[[s]])] <- pb$vv[r, s]
		if (identical(prior$kind, "joint")) {
			# sender and receiver rows of one actor are coupled
			H[cbind(ix$u[[r]], ix$v[[s]])] <- pb$uv[r, s]
			H[cbind(ix$v[[s]], ix$u[[r]])] <- pb$uv[r, s]
		}
	}
	H
}

# cavity quantities of a fitted iterate: cavity cell weights at the factor and additive predictor
# variances, the implied cavity covariances of the factor rows (.ae_mf_cavity_cov()), the product
# variance and the cavity weights at it, which carry the coefficient information (as .ae_mf_refine());
# `lap` is the iterate's Laplace record (`vcell`, and `cvpair` under a dyadic correlation)
.lde_fac_cavity <- function(ctx, st, lap) {
	cx <- list(R = ctx$R, sym = isTRUE(ctx$sym))
	sta <- st
	sta$O <- NULL
	m <- .lde_eta(ctx, sta) + st$O
	vcell <- lap$vcell
	# cavity covariances of the factor rows use independent cells; cell pairs only enter the
	# coefficient information
	ctx0 <- ctx
	ctx0$pair_info <- NULL
	nu <- .lde_cavity_w(ctx0, m, st$Vpos + vcell, lap, st$hyp)$w
	C <- rowSums(nu, dims = 2L)
	S <- .lde_fac_state(ctx, st$pos, st$hyp$fac)
	rest <- .ae_mf_cavity_cov(cx, S, C)
	S$Sig <- rest$Sig
	if (!isTRUE(ctx$sym)) S$Om <- rest$Om
	pp <- .lde_fac_product(ctx, S)
	list(
		S = S,
		pos = .lde_fac_pos(ctx, S),
		C = C,
		Vpos = pp$Vpos,
		w_info = .lde_cavity_w(ctx, m, pp$Vpos + vcell, lap, st$hyp)
	)
}

# coefficient covariance with the factors integrated out, from the information of (levels, coefficients,
# factors) with the paths eliminated by the block-tridiagonal solver, formed once at the fitted iterate
# (as .ae_eb_extend()); `extra` adds var J J' for the uncertainty of the prior scale
.lde_fac_joint <- function(ctx, st, ev, extra = NULL) {
	na <- ctx$na; nb <- ctx$nb; Tt <- ctx$Tt; q <- ctx$q; m <- ctx$mth; p <- ctx$pall; R <- ctx$R
	cavity <- identical(ctx$fscale, "cavity")
	if (cavity) {
		cav <- .lde_fac_cavity(ctx, st, ev$lap)
		sta <- st
		sta$O <- NULL
		wk_i <- ev$wk
		wk_i$w <- cav$w_info$w
		# cross weights of the cell pairs under a dyadic correlation, NULL otherwise
		wk_i$cw <- cav$w_info$cw
		wk_i$z <- ev$wk$z - st$O
		sm <- .lde_system(ctx, sta, wk_i)
		Wall <- cav$w_info$w
		Call <- cav$w_info$cw
		S <- cav$S
	} else {
		sm <- ev$sysm
		Wall <- ev$wk$w
		Call <- ev$wk$cw
		S <- .lde_fac_state(ctx, st$pos, st$hyp$fac)
	}
	U <- if (ctx$sym) S$F else S$U
	V <- if (ctx$sym) sweep(S$F, 2L, S$sg, "*") else S[["V"]]
	Sig <- S$Sig
	Om <- if (ctx$sym) {
		sgm <- outer(ctx$sgn, ctx$sgn)
		array(vapply(seq_len(na), function(i) sgm * Sig[, , i], matrix(0, R, R)), c(R, R, na))
	} else S$Om
	nf <- (if (ctx$sym) na else na + nb) * R
	Hff <- .lde_fac_ridge(ctx, st$hyp$fac)
	Hthf <- matrix(0, m, nf)
	C <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		Wt <- Wall[, , t]
		Wt[!ctx$obs[, , t] | !is.finite(Wt)] <- 0
		Ct <- if (is.null(Call)) NULL else Call[, , t]
		if (!is.null(Ct)) Ct[!ctx$obs[, , t] | !is.finite(Ct)] <- 0
		wb <- .lde_fac_wblocks(ctx, Wt, U, V, Sig, Om, Ct)
		Xt <- ctx$Xc[[t]]
		XW <- Xt * as.vector(Wt)
		idx <- ctx$thidx[[t]]
		for (k in seq_len(p)) Hthf[idx[k], ] <- Hthf[idx[k], ] + .lde_wsum(ctx, matrix(XW[, k], na, nb), U, V)
		# the cross information of the pairs of cells: (theta_k, u_j r) = sum_i c_ij x_ij,k v_i r,
		# (theta_k, v_i r) = sum_j c_ij x_ij,k u_j r
		if (!is.null(Ct)) {
			for (k in seq_len(p)) Hthf[idx[k], ] <- Hthf[idx[k], ] + .lde_wsum(
				ctx,
				t(matrix(Xt[, k], na, nb) * Ct),
				U,
				V
			)
		}
		Hff <- Hff + wb$Kww
		C[[t]] <- cbind(sm$C[[t]], wb$Kxw)
	}
	ip <- q + seq_len(m)
	ifi <- q + m + seq_len(nf)
	Hz <- matrix(0, q + m + nf, q + m + nf)
	Hz[seq_len(q + m), seq_len(q + m)] <- sm$Hzz
	Hz[ip, ifi] <- Hthf
	Hz[ifi, ip] <- t(Hthf)
	Hz[ifi, ifi] <- Hff
	fj <- .lde_factor(sm$A, sm$b, C, Hz, q, q, unit = .lde_cunit(ctx))
	Vcov <- fj$Szi[ip, ip, drop = FALSE]
	if (!is.null(extra) && length(extra$J) == m && is.finite(extra$var) && extra$var > 0) {
		Vcov <- Vcov + extra$var * tcrossprod(extra$J)
	}
	Vcov
}

# coefficient uncertainty from the unknown factor-prior scale (cavity rule, as .ae_mf_scale_uncertainty()):
# a normal in the log scale multiplier with precision c_data + c_prior, and the coefficient response
# J = d theta / d log f, both from two refits with the scales held at exp(+-h); returns J, `var` and `sd_log`
.lde_fac_scale_uncertainty <- function(ctx, st, h = 0.15, tol = 1e-5, max_iter = 200L) {
	cx <- list(R = ctx$R, sym = isTRUE(ctx$sym))
	ctx1 <- ctx
	ctx1$hold_scale <- TRUE
	# the dyadic correlation is held at the fit's value; its uncertainty is a separate matter
	ctx1$rho_hold <- TRUE
	one <- function(sgn) {
		s1 <- st
		s1$hyp$fac <- .lde_fac_fix(.ae_mf_prior_shift(st$hyp$fac, exp(sgn * h)))
		r1 <- .lde_fit_core_dyn(ctx1, s1, max_iter = max_iter, tol = tol)
		cav <- .lde_fac_cavity(ctx1, r1$st, r1$ev$lap)
		list(
			theta = r1$st$theta,
			g = .ae_mf_scale_score(cx, .lde_fac_state(ctx1, r1$st$pos, r1$st$hyp$fac), cav$C),
			ok = isTRUE(r1$converged)
		)
	}
	a <- one(1)
	b <- one(-1)
	c_data <- -(a$g - b$g) / (2 * h)
	c_prior <- .ae_mf_scale_prior_curvature(cx, .lde_fac_state(ctx, st$pos, st$hyp$fac))
	cc <- max(c_data, 0) + c_prior
	list(
		J = (a$theta - b$theta) / (2 * h),
		c = cc,
		c_data = c_data,
		c_prior = c_prior,
		var = if (is.finite(cc) && cc > 0) 1 / cc else NA_real_,
		sd_log = if (is.finite(cc) && cc > 0) 1 / sqrt(cc) else NA_real_,
		ok = a$ok && b$ok
	)
}

# coefficient covariance marginal over the static factors, plus the scale uncertainty under the
# cavity rule unless options(lame.eb_scale_uncertainty = FALSE)
.lde_fac_vcov <- function(ctx, st, ev, ctx_info = ctx) {
	extra <- NULL
	if (identical(ctx$fscale, "cavity") && .ae_mf_uncertainty_on()) {
		su <- tryCatch(.lde_fac_scale_uncertainty(ctx, st), error = function(e) NULL)
		if (!is.null(su) && is.finite(su$var) && all(is.finite(su$J))) extra <- su
	}
	list(vcov = .lde_fac_joint(ctx_info, st, ev, extra), scale = extra)
}
