# empirical-Bayes ALS estimator (bias_correct = "eb"). the actor effects are random,
# (a_i, b_i) ~ N(0, Sab) under the lame() hyperprior (Sab0, eta0), and node
# covariates enter as fixed effects gamma. (mu, beta, gamma, a, b) are solved
# jointly at each IRLS weight, Sab takes an EM update whose fixed point is the
# inverse-Wishart posterior mean (shift = -(d + 1)), and every cell score loses
# the Laplace correction kappa_c = 0.5 (h'_c v_c - h_c (u_i + u_j)) so the
# coefficients are not shrunk towards zero. a symmetric fit halves the ordered
# cells (nu = 1/2), and mu is the population intercept (effects are not re-centred).

# BCD sweeps per reweighting step of an empirical-Bayes fit
.ae_eb_sweeps <- 5L

# the additive-effects prior lame() uses for binary data (probit-moment scale and
# degrees of freedom from the first slice), so the penalty is the MCMC prior
.ae_eb_hyperprior <- function(Yarr, bip) {
	nact <- max(dim(Yarr)[1:2])
	fallback <- list(Sab0 = diag(2), eta0 = round(4 + 3 * nact / 100))
	Y1 <- matrix(Yarr[, , 1L], dim(Yarr)[1L], dim(Yarr)[2L])
	ydist <- table(Y1)
	if (!length(ydist)) return(fallback)
	ymode <- as.numeric(names(ydist)[ydist == max(ydist)])[1L]
	YB <- 1 * (Y1 != ymode)
	ybar <- mean(YB, na.rm = TRUE)
	if (!is.finite(ybar) || ybar <= 0 || ybar >= 1) return(fallback)
	E <- (YB - ybar) / stats::dnorm(stats::qnorm(ybar))
	if (!bip && nrow(E) == ncol(E)) diag(E) <- 0
	a_tmp <- rowMeans(E, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
	b_tmp <- colMeans(E, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
	vscale <- mean(c(stats::var(a_tmp), stats::var(b_tmp)))
	if (!is.finite(vscale) || vscale <= 0) vscale <- 1
	PHAT <- stats::pnorm(stats::qnorm(ybar) + outer(a_tmp, b_tmp, "+"))
	vdfmlt <- 0.25 / mean(PHAT * (1 - PHAT))
	if (!is.finite(vdfmlt) || vdfmlt < 1) vdfmlt <- 1
	list(Sab0 = diag(2) * vscale, eta0 = round(4 * vdfmlt))
}

# the state of the actor-effect prior carried through the IRLS loop
.ae_eb_prior <- function(prep, symmetric = FALSE, family = "binary", s20 = NULL) {
	kind <- if (symmetric) "symmetric" else if (prep$bip) "bipartite" else "directed"
	hp <- if (identical(family, "poisson")) .ae_eb_hyperprior_pois(prep$Yarr) else
	if (identical(family, "normal")) .ae_eb_hyperprior_norm(prep$Yarr) else .ae_eb_hyperprior(prep$Yarr, prep$bip)
	cells <- rowSums(!is.na(prep$Yarr), dims = 2L)
	cnt_r <- rowSums(cells); cnt_c <- colSums(cells)
	act_a <- cnt_r > 0; act_b <- cnt_c > 0
	if (kind != "bipartite") act_a <- act_b <- act_a | act_b
	centre <- function(W, active) {
		if (is.null(W) || !ncol(W)) return(list(W = NULL, mean = numeric(0)))
		m <- colMeans(W[active, , drop = FALSE])
		list(W = sweep(W, 2L, m), mean = m)
	}
	cr <- centre(prep$W_row, cnt_r > 0)
	cc <- centre(prep$W_col, cnt_c > 0)
	out <- list(
		kind = kind,
		Sab0 = hp$Sab0,
		eta0 = hp$eta0,
		shift = -((if (kind == "directed") 2 else 1) + 1),
		nu = if (symmetric) 0.5 else 1,
		nr = prep$nr,
		nc = prep$nc,
		q = if (symmetric) prep$nr else prep$nr + prep$nc,
		act_a = act_a,
		act_b = act_b,
		Wr = cr$W,
		Wc = cc$W,
		wbar_r = cr$mean,
		wbar_c = cc$mean,
		sab = hp$Sab0,
		gamma = NULL,
		sab_change = Inf,
		hprime = NULL
	)
	# a count network has Poisson-lognormal cell noise (ame_als_pln.R) whose variance
	# has the sampler's inverse-gamma prior on the scale s20
	if (identical(family, "poisson")) {
		out$family <- "poisson"
		out$ve_prior <- c(nu0 = .ae_pln_ve_prior[["nu0"]], s20 = if (is.null(s20)) .ae_pln_ve_prior[["s20"]] else s20)
	}
	# a normal network has Gaussian cell noise (ame_als_gauss.R) with the sampler's
	# variance prior on the scale s20
	if (identical(family, "normal")) {
		out$family <- "normal"
		out$ve_prior <- .ae_gauss_ve_prior(prep$Yarr, s20)
	}
	out
}

# prior-mean design of the stacked effects: q x pg, or NULL
.ae_eb_mean_design <- function(ab) {
	pr <- if (is.null(ab$Wr)) 0L else ncol(ab$Wr)
	pc <- if (is.null(ab$Wc)) 0L else ncol(ab$Wc)
	if (pr + pc == 0L) return(NULL)
	if (identical(ab$kind, "symmetric")) return(cbind(ab$Wr, ab$Wc))
	M <- matrix(0, ab$q, pr + pc)
	if (pr) M[seq_len(ab$nr), seq_len(pr)] <- ab$Wr
	if (pc) M[ab$nr + seq_len(ab$nc), pr + seq_len(pc)] <- ab$Wc
	M
}

# prior precision of the stacked effects: diagonal, or 2 x 2 diagonal blocks for
# a directed fit; the entries are kept in attr `ae_prior` to add them by index
.ae_eb_precision <- function(sab, ab) {
	if (identical(ab$kind, "symmetric")) {
		q <- ab$q
		out <- matrix(0, q, q)
		d <- 1 / sab[1L, 1L]
		out[cbind(seq_len(q), seq_len(q))] <- d
		attr(out, "ae_prior") <- list(d = rep(d, q))
		return(out)
	}
	nr <- ab$nr; nc <- ab$nc
	if (identical(ab$kind, "bipartite")) {
		d <- c(rep(1 / sab[1L, 1L], nr), rep(1 / sab[2L, 2L], nc))
		q <- nr + nc
		out <- matrix(0, q, q)
		out[cbind(seq_len(q), seq_len(q))] <- d
		attr(out, "ae_prior") <- list(d = d)
		return(out)
	}
	P <- solve(sab)
	ia <- seq_len(nr); ib <- nr + ia
	out <- matrix(0, 2L * nr, 2L * nr)
	out[cbind(ia, ia)] <- P[1L, 1L]
	out[cbind(ib, ib)] <- P[2L, 2L]
	out[cbind(ia, ib)] <- P[1L, 2L]
	out[cbind(ib, ia)] <- P[2L, 1L]
	attr(out, "ae_prior") <- list(
		d = c(rep(P[1L, 1L], nr), rep(P[2L, 2L], nr)),
		nr = nr,
		o12 = P[1L, 2L],
		o21 = P[2L, 1L]
	)
	out
}

# inverse of a penalised information matrix: the Cholesky inverse, or the
# minimum-norm pseudo-inverse when it is numerically singular
.ae_eb_inverse <- function(H) {
	ch <- tryCatch(chol(H), error = function(e) NULL)
	if (!is.null(ch)) {
		d <- diag(ch)
		if (min(d)^2 > 1e-12 * max(d)^2) return(chol2inv(ch))
	}
	eg <- eigen((H + t(H)) / 2, symmetric = TRUE)
	keep <- eg$values > max(eg$values) * 1e-10
	eg$vectors[, keep, drop = FALSE] %*% (t(eg$vectors[, keep, drop = FALSE]) / eg$values[keep])
}

# inverses of `H` and of its actor-effect block `H[s, s]` through the Schur
# complement of the block; two separate inverses when either factor is singular
.ae_eb_inverse_pair <- function(H, s) {
	D <- H[s, s, drop = FALSE]
	chD <- tryCatch(chol(D), error = function(e) NULL)
	r <- seq_len(nrow(H))[-s]
	apart <- function() list(S = .ae_eb_inverse(H), Sc = .ae_eb_inverse(D))
	if (is.null(chD) || !length(r)) return(apart())
	dD <- diag(chD)
	if (!(min(dD)^2 > 1e-12 * max(dD)^2)) return(apart())
	Dinv <- chol2inv(chD)
	B <- H[r, s, drop = FALSE]
	BD <- B %*% Dinv
	Sg <- H[r, r, drop = FALSE] - tcrossprod(BD, B)
	chS <- tryCatch(chol((Sg + t(Sg)) / 2), error = function(e) NULL)
	if (is.null(chS)) return(apart())
	d <- c(dD, diag(chS))
	if (!(min(d)^2 > 1e-12 * max(d)^2)) return(apart())
	Srr <- chol2inv(chS)
	Srs <- -(Srr %*% BD)
	S <- matrix(0, nrow(H), ncol(H))
	S[s, s] <- Dinv + crossprod(backsolve(chS, BD, transpose = TRUE))
	S[r, r] <- Srr
	S[r, s] <- Srs
	S[s, r] <- t(Srs)
	list(S = S, Sc = Dinv)
}

# Fisher information of the stacked actor effects at the IRLS weights (nu halves
# the doubly counted cells of a symmetric network)
.ae_eb_info_s <- function(wsum, ab) {
	rs <- rowSums(wsum); cs <- colSums(wsum)
	if (identical(ab$kind, "symmetric")) {
		return(ab$nu * (diag(rs + cs, ab$nr) + wsum + t(wsum)))
	}
	rbind(cbind(diag(rs, ab$nr), wsum), cbind(t(wsum), diag(cs, ab$nc)))
}

# likelihood part of the normal equations of (mu, beta | actor effects | gamma)
# at given IRLS weights, from the BCD dyad-level sums `st`
.ae_eb_context <- function(ab, st, wsum, p) {
	pa <- p + 1L
	M <- .ae_eb_mean_design(ab)
	pg <- if (is.null(M)) 0L else ncol(M)
	q <- ab$q
	m <- pa + q + pg
	i_th <- seq_len(pa); i_s <- pa + seq_len(q); i_g <- pa + q + seq_len(pg)
	C <- if (identical(ab$kind, "symmetric")) ab$nu * (st$gr + st$gc) else rbind(st$gr, st$gc)
	H <- matrix(0, m, m)
	H[i_th, i_th] <- ab$nu * st$xtx
	H[i_s, i_th] <- C; H[i_th, i_s] <- t(C)
	H[i_s, i_s] <- .ae_eb_info_s(wsum, ab)
	list(
		H = H,
		pa = pa,
		pg = pg,
		m = m,
		M = M,
		wsum = wsum,
		gm = st$gm,
		i_th = i_th,
		i_s = i_s,
		i_g = i_g,
		nr = ab$nr,
		nc = ab$nc,
		q = q
	)
}

# add the actor-effect penalty (s - M gamma)' Lam (s - M gamma) to the Hessian
.ae_eb_penalise <- function(ctx, Lam) {
	H <- ctx$H
	s <- ctx$i_s
	pr <- attr(Lam, "ae_prior")
	if (is.null(pr)) {
		H[s, s] <- H[s, s] + Lam
	} else {
		# the precision is diagonal or has diagonal blocks, so its entries are added by index
		dd <- cbind(s, s)
		H[dd] <- H[dd] + pr$d
		if (!is.null(pr$nr)) {
			ia <- s[seq_len(pr$nr)]; ib <- s[pr$nr + seq_len(pr$nr)]
			iab <- cbind(ia, ib); iba <- cbind(ib, ia)
			H[iab] <- H[iab] + pr$o12
			H[iba] <- H[iba] + pr$o21
		}
	}
	if (ctx$pg > 0L) {
		LM <- Lam %*% ctx$M
		H[s, ctx$i_g] <- -LM
		H[ctx$i_g, s] <- t(-LM)
		H[ctx$i_g, ctx$i_g] <- crossprod(ctx$M, LM)
	}
	H
}

# d h / d eta of the Fisher weight h(eta) of one cell (probit: h' = h (lam0 - lam1 - 2 eta);
# logit: h = p (1 - p)); `terms` holds the normal log density and tail probabilities
.ae_eb_hprime <- function(eta, link, terms = NULL) {
	# a count cell (Poisson-lognormal): h' is among its terms
	if (identical(link, "log")) return(terms$hp)
	# the information of a Gaussian cell does not depend on the predictor
	if (identical(link, "identity")) return(array(0, .ae_dims(eta)))
	if (identical(link, "probit")) {
		if (is.null(terms)) terms <- .ae_probit_terms(eta, NULL)
		if (!identical(.ae_engine(), "r")) {
			out <- ae_probit_hprime_cpp(eta, terms$ld, terms$lp1, terms$lp0, .ae_threads())
			dim(out) <- dim(eta)
			return(out)
		}
		lam1 <- exp(terms$ld - terms$lp1)
		lam0 <- exp(terms$ld - terms$lp0)
		return(lam1 * lam0 * (lam0 - lam1 - 2 * eta))
	}
	p <- stats::plogis(eta)
	p * (1 - p) * (1 - 2 * p)
}

# sum_t h' [1, x] for every dyad (n_row * n_col x (p + 1)); h' can be negative, so
# the R reference sums its positive and negative parts apart
.ae_eb_hprime_sums_r <- function(hprime, X, nr, nc, Tt, p) {
	hp <- hprime
	hp[!is.finite(hp)] <- 0
	pos <- ae_bcd_stats_cpp(pmax(hp, 0), hp, X, nr, nc, Tt, p)$gm
	neg <- ae_bcd_stats_cpp(pmax(-hp, 0), hp, X, nr, nc, Tt, p)$gm
	pos - neg
}
.ae_eb_hprime_sums <- function(hprime, X, nr, nc, Tt, p) {
	if (identical(.ae_engine(), "r")) return(.ae_eb_hprime_sums_r(hprime, X, nr, nc, Tt, p))
	hp <- hprime
	hp[!is.finite(hp)] <- 0
	ae_weighted_sums_cpp(hp, X, nr, nc, Tt, p, .ae_threads())
}

# Laplace cell correction at an iterate: the (mu, beta) and gamma score corrections
# and their dyad sum kd. `Hss` is K_w + Lam, `hpg` the time sums of h' [1, x]
.ae_eb_laplace <- function(ctx, ab, Hss, st, hpg, Lam, Sc = NULL) {
	nr <- ctx$nr; nc <- ctx$nc
	# covariance of the actor effects given the coefficients, the inverse of `Hss`
	if (is.null(Sc)) Sc <- .ae_eb_inverse(Hss)
	hp <- matrix(hpg[, 1L], nr, nc)
	if (identical(ab$kind, "symmetric")) {
		dg <- diag(Sc)
		v <- outer(dg, dg, "+") + 2 * Sc
		u <- as.vector(Sc %*% rowSums(hp * v))
		cell <- crossprod(hpg, as.vector(v)) -
			as.vector(crossprod(st$gr + st$gc, u))
		kd <- 0.5 * (hp * v - ctx$wsum * outer(u, u, "+"))
	} else {
		ia <- seq_len(nr); ib <- nr + seq_len(nc)
		dg <- diag(Sc)
		v <- outer(dg[ia], dg[ib], "+") + 2 * Sc[ia, ib, drop = FALSE]
		u <- as.vector(Sc %*% c(rowSums(hp * v), colSums(hp * v)))
		cell <- crossprod(hpg, as.vector(v)) -
			as.vector(crossprod(st$gr, u[ia]) + crossprod(st$gc, u[ib]))
		kd <- 0.5 * (hp * v - ctx$wsum * outer(u[ia], u[ib], "+"))
	}
	# the caller applies nu (a symmetric network counts each dyad once)
	list(
		corr_th = 0.5 * as.vector(cell),
		corr_g = if (ctx$pg > 0L) 0.5 * as.vector(crossprod(ctx$M, Lam %*% u)) else numeric(0),
		kd = kd
	)
}

# posterior covariance of the residual effects alpha = s - M gamma, marginal over
# the coefficients and gamma (flat prior), from the inverse `S` of the information
.ae_eb_css <- function(ctx, S) {
	s <- ctx$i_s
	Css <- S[s, s, drop = FALSE]
	if (ctx$pg > 0L) {
		MS <- ctx$M %*% S[ctx$i_g, s, drop = FALSE]
		Css <- Css - MS - t(MS) + ctx$M %*% S[ctx$i_g, ctx$i_g, drop = FALSE] %*% t(ctx$M)
	}
	Css
}

# one penalised linear solve plus the EM update of Sab
#   Omat : the multiplicative term (NULL at R = 0)
#   st   : the BCD dyad-level sums (sz, szr, szc, xtz, gm, gr, gc)
#   hpg  : sum_t h' [1, x] per dyad, or NULL for the plain penalised solve
.ae_eb_step <- function(ctx, ab, Omat, st, sab, hpg = NULL) {
	nr <- ctx$nr; nc <- ctx$nc; pg <- ctx$pg
	nu <- ab$nu; sym <- identical(ab$kind, "symmetric")
	Lam <- .ae_eb_precision(sab, ab)
	H <- .ae_eb_penalise(ctx, Lam)
	# right-hand side: the working response less the multiplicative term
	g0 <- st$xtz; rs <- st$szr; cs <- st$szc
	if (!is.null(Omat)) {
		g0 <- g0 - ae_bcd_gtc_cpp(numeric(nr), numeric(nc), Omat, st$gm)
		wo <- ctx$wsum * Omat
		rs <- rs - rowSums(wo); cs <- cs - colSums(wo)
	}
	rhs_th <- nu * g0
	rhs_s <- if (sym) nu * (rs + cs) else c(rs, cs)
	s <- ctx$i_s
	lap <- NULL
	rhs <- c(rhs_th, rhs_s, numeric(pg))
	if (!is.null(hpg)) {
		# Laplace-marginal equations for (mu, beta, gamma), corrected at the incoming iterate
		inv <- .ae_eb_inverse_pair(H, s)
		S <- inv$S
		lap <- .ae_eb_laplace(ctx, ab, H[s, s, drop = FALSE], st, hpg, Lam, Sc = inv$Sc)
		rhs <- rhs - c(nu * lap$corr_th, numeric(ctx$q), lap$corr_g)
	} else {
		S <- .ae_eb_inverse(H)
	}
	x <- as.vector(S %*% rhs)
	gam <- x[ctx$i_g]
	alpha <- x[s] - if (pg > 0L) as.vector(ctx$M %*% gam) else 0
	a <- x[s][seq_len(nr)]
	b <- if (sym) a else x[s][nr + seq_len(nc)]
	list(
		theta = x[ctx$i_th],
		a = a,
		b = b,
		gamma = gam,
		alpha = alpha,
		kd = lap$kd,
		sab_new = .ae_eb_update_sab(ab, alpha, .ae_eb_css(ctx, S))
	)
}

# Laplace correction and EM image of Sab at a fitted iterate, from its own weights
# `w` and predictor `eta`, so the equations are checked at the iterate itself
.ae_eb_at_iterate <- function(prep, fit, w, eta, terms, link, ab) {
	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt; p <- prep$p_dyad
	observed <- is.finite(prep$Yarr)
	Warr <- w
	Warr[!observed | !is.finite(Warr)] <- 0
	st <- ae_bcd_stats_cpp(Warr, Warr, prep$X, nr, nc, Tt, p, .ae_threads())
	wsum <- rowSums(Warr, dims = 2L)
	ctx <- .ae_eb_context(ab, st, wsum, p)
	Lam <- .ae_eb_precision(fit$sab, ab)
	H <- .ae_eb_penalise(ctx, Lam)
	if (identical(link, "identity")) {
		# a Gaussian cell: the Laplace approximation is exact and there is no correction
		inv <- list(S = .ae_eb_inverse(H))
		lap <- list(corr_th = numeric(ctx$pa), corr_g = numeric(ctx$pg), kd = matrix(0, nr, nc))
	} else {
		hprime <- .ae_eb_hprime(eta, link, terms)
		hprime[!observed | !is.finite(hprime)] <- 0
		hpg <- .ae_eb_hprime_sums(hprime, prep$X, nr, nc, Tt, p)
		inv <- .ae_eb_inverse_pair(H, ctx$i_s)
		lap <- .ae_eb_laplace(ctx, ab, H[ctx$i_s, ctx$i_s, drop = FALSE], st, hpg, Lam, Sc = inv$Sc)
	}
	al <- if (identical(ab$kind, "symmetric")) fit$a else c(fit$a, fit$b)
	if (ctx$pg > 0L && length(fit$gamma) == ctx$pg) {
		al <- al - as.vector(ctx$M %*% fit$gamma)
	}
	list(
		corr_th = lap$corr_th,
		corr_g = lap$corr_g,
		kd = lap$kd,
		sab_new = .ae_eb_update_sab(ab, al, .ae_eb_css(ctx, inv$S))
	)
}

# EM update of the additive-effect covariance from the posterior second moments
.ae_eb_update_sab <- function(ab, alpha, Css) {
	nr <- ab$nr; nc <- ab$nc
	al_a <- alpha[seq_len(nr)]
	ia <- ab$act_a
	ca <- diag(Css)[seq_len(nr)]
	if (identical(ab$kind, "symmetric")) {
		va <- (ab$eta0 * ab$Sab0[1L, 1L] + sum(al_a[ia]^2 + ca[ia])) /
			(sum(ia) + ab$eta0 + ab$shift)
		va <- max(va, .ae_sab_floor(ab))
		return(diag(va, 2L))
	}
	al_b <- alpha[nr + seq_len(nc)]
	ib <- ab$act_b
	cb <- diag(Css)[nr + seq_len(nc)]
	if (identical(ab$kind, "bipartite")) {
		va <- (ab$eta0 * ab$Sab0[1L, 1L] + sum(al_a[ia]^2 + ca[ia])) / (sum(ia) + ab$eta0 + ab$shift)
		vb <- (ab$eta0 * ab$Sab0[2L, 2L] + sum(al_b[ib]^2 + cb[ib])) / (sum(ib) + ab$eta0 + ab$shift)
		fl <- .ae_sab_floor(ab)
		return(diag(c(max(va, fl), max(vb, fl))))
	}
	cab <- Css[cbind(seq_len(nr), nr + seq_len(nc))]
	act <- ia
	S2 <- matrix(
		c(
			sum(al_a[act]^2 + ca[act]),
			sum(al_a[act] * al_b[act] + cab[act]),
			sum(al_a[act] * al_b[act] + cab[act]),
			sum(al_b[act]^2 + cb[act])
		),
		2L,
		2L
	)
	new <- (ab$eta0 * ab$Sab0 + S2) / (sum(act) + ab$eta0 + ab$shift)
	new <- (new + t(new)) / 2
	# keep the covariance positive definite
	d <- sqrt(pmax(diag(new), .ae_sab_floor(ab)))
	r <- max(min(new[1L, 2L] / (d[1L] * d[2L]), 0.99), -0.99)
	new[1L, 2L] <- new[2L, 1L] <- r * d[1L] * d[2L]
	diag(new) <- d^2
	new
}

# actor-effect part of the penalised objective at `sab` (deviance scale); with
# `wsum` it adds the log determinant of the actor-effect information
.ae_eb_penalty <- function(f, sab, ab, wsum = NULL) {
	nr <- ab$nr; nc <- ab$nc
	sym <- identical(ab$kind, "symmetric")
	al <- if (sym) f$a else c(f$a, f$b)
	M <- .ae_eb_mean_design(ab)
	if (!is.null(M) && length(f$gamma) == ncol(M)) al <- al - as.vector(M %*% f$gamma)
	ld <- if (is.null(wsum)) 0 else {
		K <- .ae_eb_precision(sab, ab) + .ae_eb_info_s(wsum, ab)
		2 * sum(log(diag(chol(K))))
	}
	al_a <- al[seq_len(nr)]
	if (sym) {
		quad <- sum(al_a[ab$act_a]^2) / sab[1L, 1L] + sum(ab$act_a) * log(sab[1L, 1L])
	} else if (identical(ab$kind, "bipartite")) {
		al_b <- al[nr + seq_len(nc)]
		quad <- sum(al_a[ab$act_a]^2) / sab[1L, 1L] + sum(al_b[ab$act_b]^2) / sab[2L, 2L] +
			sum(ab$act_a) * log(sab[1L, 1L]) + sum(ab$act_b) * log(sab[2L, 2L])
	} else {
		al_b <- al[nr + seq_len(nc)]
		A <- cbind(al_a, al_b)[ab$act_a, , drop = FALSE]
		quad <- sum((A %*% solve(sab)) * A) + nrow(A) * log(det(sab))
	}
	(quad + ld) / ab$nu
}

# hyperprior part of the objective at `sab`, c log|Sab| + tr(eta0 Sab0 Sab^-1) with
# c = eta0 + shift; needed to compare fits that converged at different covariances
.ae_eb_hyper <- function(sab, ab) {
	cc <- ab$eta0 + ab$shift
	h <- if (identical(ab$kind, "symmetric")) {
		cc * log(sab[1L, 1L]) + ab$eta0 * ab$Sab0[1L, 1L] / sab[1L, 1L]
	} else if (identical(ab$kind, "bipartite")) {
		sum(cc * log(diag(sab)) + ab$eta0 * diag(ab$Sab0) / diag(sab))
	} else {
		cc * log(det(sab)) + ab$eta0 * sum(diag(ab$Sab0 %*% solve(sab)))
	}
	h / ab$nu
}

# penalised Fisher information of (mu, beta, actor effects, gamma) extended by the
# factors with ridge 1 / sig2 (symmetric: F = U |L|^{1/2} with the signs of L). for
# a variational fit `vb` the factors are posterior means, a row's curvature adds
# the covariance of the other side, and the ridge is the prior precision
.ae_eb_extend <- function(ctx, ab, H0, fit, vb = NULL) {
	nr <- ctx$nr; nc <- ctx$nc; pa <- ctx$pa; m0 <- ctx$m
	wsum <- ctx$wsum; gm <- ctx$gm
	Mk <- lapply(seq_len(pa), function(k) matrix(gm[, k], nr, nc))
	i_s <- ctx$i_s
	if (identical(ab$kind, "symmetric")) {
		R <- ncol(fit$U)
		Fm <- if (is.null(vb)) sweep(fit$U, 2L, sqrt(abs(fit$L)), "*") else vb$S$F
		sg <- if (is.null(vb)) ifelse(fit$L < 0, -1, 1) else vb$S$sg
		H <- matrix(0, m0 + R * nr, m0 + R * nr)
		H[seq_len(m0), seq_len(m0)] <- H0
		jf <- function(r) m0 + (r - 1L) * nr + seq_len(nr)
		for (r in seq_len(R)) {
			f_r <- Fm[, r]
			# (a, F_r): s_r [diag(wsum f_r) + wsum * f_r]
			blk <- sg[r] * (diag(as.vector(wsum %*% f_r), nr) + wsum * f_r)
			H[i_s, jf(r)] <- blk
			for (k in seq_len(pa)) {
				H[ctx$i_th[k], jf(r)] <- sg[r] * as.vector(Mk[[k]] %*% f_r)
			}
			for (r2 in seq_len(R)) {
				f_2 <- Fm[, r2]
				H[jf(r), jf(r2)] <- sg[r] * sg[r2] *
					(diag(as.vector(wsum %*% (f_r * f_2)), nr) + wsum * outer(f_2, f_r))
				if (!is.null(vb)) {
					H[cbind(jf(r), jf(r2))] <- H[cbind(jf(r), jf(r2))] +
						sg[r] * sg[r2] * as.vector(wsum %*% vb$S$Sig[r, r2, ])
				}
			}
			H[cbind(jf(r), jf(r))] <- H[cbind(jf(r), jf(r))] + 1 / fit$sig2_u[r]
		}
		H[lower.tri(H)] <- t(H)[lower.tri(H)]
		return(H)
	}
	U <- fit$U; V <- fit[["V"]]
	R <- ncol(U)
	mf <- R * (nr + nc)
	H <- matrix(0, m0 + mf, m0 + mf)
	H[seq_len(m0), seq_len(m0)] <- H0
	iu <- function(r) m0 + (r - 1L) * nr + seq_len(nr)
	iv <- function(r) m0 + R * nr + (r - 1L) * nc + seq_len(nc)
	ia <- i_s[seq_len(nr)]; ib <- i_s[nr + seq_len(nc)]
	pb <- if (!is.null(vb)) .ae_mf_prior_blocks(vb$S$prior, R)
	for (r in seq_len(R)) {
		wv <- as.vector(wsum %*% V[, r]); wu <- as.vector(crossprod(wsum, U[, r]))
		# (a, U), (a, V), (b, U), (b, V)
		H[cbind(ia, iu(r))] <- wv
		H[ia, iv(r)] <- wsum * U[, r]
		H[ib, iu(r)] <- t(wsum) * V[, r]
		H[cbind(ib, iv(r))] <- wu
		# (theta_k, U), (theta_k, V)
		for (k in seq_len(pa)) {
			H[ctx$i_th[k], iu(r)] <- as.vector(Mk[[k]] %*% V[, r])
			H[ctx$i_th[k], iv(r)] <- as.vector(crossprod(Mk[[k]], U[, r]))
		}
		for (r2 in seq_len(R)) {
			H[cbind(iu(r), iu(r2))] <- as.vector(wsum %*% (V[, r] * V[, r2]))
			H[cbind(iv(r), iv(r2))] <- as.vector(crossprod(wsum, U[, r] * U[, r2]))
			H[iu(r), iv(r2)] <- wsum * outer(U[, r2], V[, r])
			if (!is.null(vb)) {
				H[cbind(iu(r), iu(r2))] <- H[cbind(iu(r), iu(r2))] +
					as.vector(wsum %*% vb$S$Om[r, r2, ]) + pb$uu[r, r2]
				H[cbind(iv(r), iv(r2))] <- H[cbind(iv(r), iv(r2))] +
					as.vector(crossprod(wsum, vb$S$Sig[r, r2, ])) + pb$vv[r, r2]
				if (identical(vb$S$prior$kind, "joint")) {
					H[cbind(iu(r), iv(r2))] <- H[cbind(iu(r), iv(r2))] + pb$uv[r, r2]
				}
			}
		}
		if (is.null(vb)) {
			H[cbind(iu(r), iu(r))] <- H[cbind(iu(r), iu(r))] + 1 / fit$sig2_u[r]
			H[cbind(iv(r), iv(r))] <- H[cbind(iv(r), iv(r))] + 1 / fit$sig2_v[r]
		}
	}
	# the cross blocks were filled in the upper triangle only
	H[lower.tri(H)] <- t(H)[lower.tri(H)]
	# the pair information of a dyadic correlation couples the factors of a dyad's actors
	if (!is.null(ctx$cross)) H <- .ae_pair_extend(H, ctx, ab, fit)
	H
}

# the joint Laplace posterior of the effects and factors: its covariance, the inverse
# of the full penalised information, is formed once at the fitted iterate

# the layout of the stacked parameter vector built by .ae_eb_extend
.ae_eb_layout <- function(ctx, ab, R) {
	nr <- ctx$nr; nc <- ctx$nc; m0 <- ctx$m
	sym <- identical(ab$kind, "symmetric")
	iu <- function(r) m0 + (r - 1L) * nr + seq_len(nr)
	iv <- if (sym) iu else function(r) m0 + R * nr + (r - 1L) * nc + seq_len(nc)
	list(
		sym = sym,
		R = R,
		nr = nr,
		nc = nc,
		m0 = m0,
		iu = iu,
		iv = iv,
		size = m0 + (if (sym) R * nr else R * (nr + nc)),
		ia = ctx$i_s[seq_len(nr)],
		ib = if (sym) ctx$i_s[seq_len(nr)] else ctx$i_s[nr + seq_len(nc)]
	)
}

# the penalised information of everything at a fitted iterate and, unless `inverse`
# is FALSE, its inverse; `w` holds the iterate's Fisher weights
.ae_eb_joint <- function(prep, ab, fit, w, R, inverse = TRUE, cross = NULL) {
	nr <- prep$nr; nc <- prep$nc; p <- prep$p_dyad; Tt <- prep$Tt
	Warr <- w; Warr[!is.finite(Warr)] <- 0
	st <- ae_bcd_stats_cpp(Warr, Warr, prep$X, nr, nc, Tt, p, .ae_threads())
	wsum <- rowSums(Warr, dims = 2L)
	ctx <- .ae_eb_context(ab, st, wsum, p)
	# the cross weights of the pairs of a fit with a dyadic correlation
	if (!is.null(cross)) {
		ctx <- .ae_pair_add_cross(ctx, .ae_pair_cross(list(nr = nr, Tt = Tt, p = p, X = prep[["X"]]), cross))
	}
	H0 <- .ae_eb_penalise(ctx, .ae_eb_precision(fit$sab, ab))
	use_f <- R > 0L && !is.null(fit$U)
	Hf <- if (use_f) .ae_eb_extend(ctx, ab, H0, fit) else H0
	lay <- .ae_eb_layout(ctx, ab, if (use_f) R else 0L)
	list(
		ctx = ctx,
		st = st,
		wsum = wsum,
		H0 = H0,
		Hf = Hf,
		lay = lay,
		S = if (inverse) .ae_eb_inverse(Hf)
	)
}

# predictive probability: the link integrated over the Gaussian Laplace posterior
# of the predictor, probit Phi(m / sqrt(1 + v)) or logit plogis(m / sqrt(1 + pi v / 8)).
# m is the mode shifted by -S q (q the gradient of half the log determinant) plus the
# factor-row covariance; v is the delta-method variance under the joint covariance

# the minimal prep, prior, fit and Fisher weights, rebuilt from a fitted object
.ae_eb_state <- function(object) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	sym <- isTRUE(object$symmetric)
	prep <- .ae_prep_from_canonical(
		object$Y,
		object[["X"]],
		object$W_row,
		object$W_col,
		object$mode,
		sym,
		object$row_names,
		object$col_names,
		object$time_names
	)
	ab <- .ae_eb_prior(prep, sym, object$family, s20 = object$eb$ve_prior[["s20"]])
	eb <- object$eb
	R <- object$R
	fit <- list(
		mu = eb$mu_int,
		beta = unname(object$beta[seq_len(prep$p_dyad)]),
		a = unname(eb$a_raw),
		b = unname(eb$b_raw),
		gamma = unname(eb$gamma),
		sab = eb$sab,
		sig2_u = eb$sig2_u,
		sig2_v = eb$sig2_v
	)
	if (object$family %in% c("poisson", "normal")) fit$ve <- eb$ve
	if (R > 0L) {
		fit$U <- unname(object$U); fit[["V"]] <- unname(object[["V"]]); fit$L <- object$L
	}
	list(prep = prep, ab = ab, fit = fit, R = R, w = object$obs_weights, cross = object$eb$pair$cross)
}

# Laplace state of a fitted object: the joint information at the fitted iterate and
# the time sums of h' at the fitted predictor
.ae_eb_laplace_state <- function(object, inverse = TRUE) {
	st <- .ae_eb_state(object)
	prep <- st$prep; ab <- st$ab; fit <- st$fit; R <- st$R
	eta0 <- simplify2array(object$EZ)
	if (length(dim(eta0)) == 2L) eta0 <- array(eta0, c(dim(eta0), 1L))
	obs <- is.finite(prep$Yarr) & is.finite(eta0)
	link <- .ae_default(object$link, "probit")
	terms <- if (identical(link, "log")) .ae_pln_terms(prep$Yarr, eta0, fit$ve) else
	if (identical(link, "identity")) NULL else
	if (identical(link, "probit")) .ae_probit_terms(eta0, NULL)
	hp3 <- .ae_eb_hprime(eta0, link, terms)
	hp3[!obs | !is.finite(hp3)] <- 0
	list(
		prep = prep,
		ab = ab,
		fit = fit,
		R = R,
		eta0 = eta0,
		obs = obs,
		link = link,
		hp = rowSums(hp3, dims = 2L),
		jp = .ae_eb_joint(prep, ab, fit, st$w, R, inverse, cross = st$cross)
	)
}

# the parameters a cell's predictor depends on (its slots: the sender and
# receiver effect and the R sender and R receiver factor entries), where they
# sit in the stacked vector, and the gradient of the predictor in each
.ae_eb_slots <- function(lp) {
	fit <- lp$fit; R <- lp$R; lay <- lp$jp$lay
	sym <- identical(lp$ab$kind, "symmetric")
	# the factors as the product uses them: F = U |L|^(1/2) with signs when symmetric
	sgn <- if (R > 0L && sym) ifelse(fit$L < 0, -1, 1) else rep(1, R)
	Fu <- Fv <- NULL
	if (R > 0L) {
		Fu <- if (sym) sweep(fit$U, 2L, sqrt(abs(fit$L)), "*") else fit$U
		Fv <- if (sym) Fu else fit[["V"]]
	}
	k <- 2L + 2L * R
	list(
		k = k,
		sgn = sgn,
		Fu = Fu,
		Fv = Fv,
		pos = function(I, J) {
			cols <- list(lay$ia[I], lay$ib[J])
			for (r in seq_len(R)) cols[[2L + r]] <- lay$iu(r)[I]
			for (r in seq_len(R)) cols[[2L + R + r]] <- lay$iv(r)[J]
			do.call(cbind, cols)
		},
		grad = function(I, J, Uc, Vc) {
			g <- matrix(1, length(I), k)
			for (r in seq_len(R)) {
				g[, 2L + r] <- sgn[r] * Vc[J, r]
				g[, 2L + R + r] <- sgn[r] * Uc[I, r]
			}
			g
		}
	)
}

# gradient q of half the log determinant of the effect and factor information
# (whose inverse is `Sz`), summed over observed dyads; `posz` maps stacked
# positions to rows of `Sz`
.ae_eb_logdet_grad <- function(lp, sl, Sz, posz, chunk = 20000L) {
	R <- lp$R; jp <- lp$jp; hp <- lp$hp; wsum <- jp$wsum
	nr <- lp$prep$nr; k <- sl$k
	q <- numeric(nrow(Sz))
	obs_d <- which(wsum > 0 | hp != 0)
	for (start in seq(1L, length(obs_d), by = chunk)) {
		idx <- obs_d[start:min(start + chunk - 1L, length(obs_d))]
		I <- (idx - 1L) %% nr + 1L; J <- (idx - 1L) %/% nr + 1L
		P <- sl$pos(I, J)
		P[] <- posz[P]
		g <- if (R > 0L) sl$grad(I, J, sl$Fu, sl$Fv) else matrix(1, length(I), 2L)
		sg <- matrix(0, length(I), k)
		for (a in seq_len(k)) for (b in seq_len(k)) sg[, a] <- sg[, a] + Sz[cbind(P[, a], P[, b])] * g[, b]
		v <- rowSums(sg * g)
		hv <- 0.5 * hp[idx] * v
		for (a in seq_len(k)) {
			contrib <- hv * g[, a]
			if (R > 0L && a > 2L) {
				pa_ <- if (a <= 2L + R) a + R else a - R # the partner slot
				r <- if (a <= 2L + R) a - 2L else a - 2L - R
				contrib <- contrib + wsum[idx] * sl$sgn[r] * sg[, pa_]
			}
			# the ordered cells of an undirected network count each dyad twice
			ag <- rowsum(lp$ab$nu * contrib, P[, a], reorder = FALSE)
			q[as.integer(rownames(ag))] <- q[as.integer(rownames(ag))] + ag[, 1L]
		}
	}
	q
}

# the posterior mean and variance of the linear predictor of every cell of
# the design `X` (the fit's own when NULL), as arrays [n_row, n_col, T]
.ae_eb_predictive <- function(object, X = NULL, displace = TRUE) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	# a variational fit has its own posterior (ame_als_meanfield_fit.R)
	if (displace && identical(object$eb$estimator, "variational")) {
		return(.ae_mf_predictive(object, X))
	}
	lp <- .ae_eb_laplace_state(object)
	prep <- lp$prep; ab <- lp$ab; fit <- lp$fit; R <- lp$R; jp <- lp$jp
	nr <- prep$nr; nc <- prep$nc; p <- prep$p_dyad; Tt <- prep$Tt
	sym <- identical(ab$kind, "symmetric")
	if (is.null(X)) X <- prep[["X"]]
	ctx <- jp$ctx; lay <- jp$lay; Sf <- jp$S
	pa <- ctx$pa
	iz <- setdiff(seq_len(nrow(jp$Hf)), c(ctx$i_th, ctx$i_g))
	posz <- integer(nrow(jp$Hf)); posz[iz] <- seq_along(iz)
	# covariance of the effects and factors given the coefficients
	it <- c(ctx$i_th, ctx$i_g)
	Sz <- Sf[iz, iz, drop = FALSE] -
		Sf[iz, it, drop = FALSE] %*% solve(Sf[it, it, drop = FALSE], Sf[it, iz, drop = FALSE])
	sl <- .ae_eb_slots(lp)
	k <- sl$k; sgn <- sl$sgn; Fu <- sl$Fu; Fv <- sl$Fv
	chunk <- 20000L
	# shift of the integrated parameters to their posterior mean: -Sz q
	shift <- if (displace) {
		-as.vector(Sz %*% .ae_eb_logdet_grad(lp, sl, Sz, posz, chunk))
	} else numeric(nrow(Sz))
	# the posterior means of the parameters
	get <- function(pos) shift[posz[pos]]
	a_t <- fit$a + get(lay$ia); b_t <- fit$b + get(lay$ib)
	if (R > 0L) {
		Ut <- Fu + vapply(seq_len(R), function(r) get(lay$iu(r)), numeric(nr))
		Vt <- if (sym) Ut else Fv + vapply(seq_len(R), function(r) get(lay$iv(r)), numeric(nc))
	}
	# the mean of the predictor apart from the coefficients' part: additive
	# effects, product and the covariance of the two factor rows
	core <- outer(a_t, b_t, "+")
	if (R > 0L) {
		for (r in seq_len(R)) {
			core <- core + sgn[r] * tcrossprod(Ut[, r], Vt[, r]) +
				sgn[r] * Sf[lay$iu(r), lay$iv(r), drop = FALSE]
		}
	}
	# variance: the coefficients (with their covariance with the effects and
	# factors of the dyad) and the dyad's own parameters, linearised at the mean
	Sth <- Sf[ctx$i_th, ctx$i_th, drop = FALSE]
	var_core <- matrix(0, nr, nc)
	cth <- array(0, c(nr, nc, pa))
	all_d <- seq_len(nr * nc)
	for (start in seq(1L, length(all_d), by = chunk)) {
		idx <- all_d[start:min(start + chunk - 1L, length(all_d))]
		I <- (idx - 1L) %% nr + 1L; J <- (idx - 1L) %/% nr + 1L
		P <- sl$pos(I, J)
		g <- if (R > 0L) sl$grad(I, J, Ut, Vt) else matrix(1, length(I), 2L)
		vz <- numeric(length(I))
		for (a in seq_len(k)) for (b in seq_len(k)) vz <- vz + g[, a] * g[, b] * Sf[cbind(P[, a], P[, b])]
		var_core[idx] <- vz
		for (m in seq_len(pa)) {
			cm <- numeric(length(I))
			for (a in seq_len(k)) cm <- cm + g[, a] * Sf[cbind(ctx$i_th[m], P[, a])]
			cth[cbind(I, J, m)] <- cm
		}
	}
	mean_a <- var_a <- array(NA_real_, c(nr, nc, Tt))
	theta <- c(fit$mu, fit$beta)
	for (t in seq_len(Tt)) {
		Xt <- cbind(rep(1, nr * nc), if (p > 0L) matrix(X[, , , t], nr * nc, p))
		mean_a[, , t] <- matrix(Xt %*% theta, nr, nc) + core
		cx <- matrix(rowSums(Xt * matrix(cth, nr * nc, pa)), nr, nc)
		var_a[, , t] <- matrix(rowSums((Xt %*% Sth) * Xt), nr, nc) + 2 * cx + var_core
	}
	list(mean = mean_a, var = pmax(var_a, 0))
}

# predictive probability from the posterior mean and variance of the predictor
.ae_eb_integrate <- function(m, v, link, ve = 0) {
	# a count: the lognormal mean exp(m + (v + ve) / 2)
	if (identical(link, "log")) return(exp(m + (v + ve) / 2))
	if (identical(link, "logit")) stats::plogis(m / sqrt(1 + pi * v / 8))
	else stats::pnorm(m / sqrt(1 + v))
}

# Laplace approximation to the marginal likelihood: effects and factors integrated,
# coefficients and scales at their estimates. `value` is minus twice its log on the
# deviance scale; the hyperpriors are left out since their constants depend on the rank.
# `ld_joint` and `ld_actor` are the log determinants of the joint and actor information
.ae_eb_evidence_core <- function(prep, ab, fit, R, eta, link, ld_joint, ld_actor) {
	y <- prep$Yarr
	if (identical(link, "log")) {
		dev <- .ae_dev_pln(.ae_pln_terms(y, eta, fit$ve))
	} else if (identical(link, "identity")) {
		dev <- .ae_dev_gauss(.ae_gauss_terms(y, eta, fit$ve))
	} else {
		terms <- if (identical(link, "probit")) .ae_probit_terms(eta, which(y == 1))
		dev <- .ae_irls_deviance(y, eta, "binary", link, terms)
	}
	quad <- .ae_eb_penalty(fit, fit$sab, ab)
	nr <- prep$nr; nc <- prep$nc
	fac <- 0
	if (R > 0L) {
		if (identical(ab$kind, "symmetric")) {
			s2 <- max(fit$sig2_u[1L], 1e-12)
			fac <- 2 * sum(abs(fit$L)) / s2 + 2 * nr * R * log(s2)
		} else {
			su <- pmax(fit$sig2_u, 1e-12); sv <- pmax(fit$sig2_v, 1e-12)
			fac <- sum(colSums(fit$U^2) / su) + nr * sum(log(su)) +
				sum(colSums(fit[["V"]]^2) / sv) + nc * sum(log(sv))
		}
	}
	value <- dev + quad + fac + ld_joint / ab$nu
	list(
		value = value,
		log_marginal = -0.5 * ab$nu * value,
		deviance = dev,
		effects = quad,
		factors = fac,
		logdet_effects = ld_actor / ab$nu,
		logdet_joint = ld_joint / ab$nu
	)
}

# the evidence of a fitted object (the value stored in fit$eb$log_marginal)
.ae_eb_evidence <- function(object) {
	if (identical(object$eb$estimator, "variational")) return(object$eb$laplace)
	lp <- .ae_eb_laplace_state(object, inverse = FALSE)
	jp <- lp$jp; ctx <- jp$ctx
	iz <- setdiff(seq_len(nrow(jp$Hf)), c(ctx$i_th, ctx$i_g))
	.ae_eb_evidence_core(
		lp$prep,
		lp$ab,
		lp$fit,
		lp$R,
		lp$eta0,
		lp$link,
		2 * sum(log(diag(chol(jp$Hf[iz, iz, drop = FALSE])))),
		2 * sum(log(diag(chol(jp$H0[ctx$i_s, ctx$i_s, drop = FALSE]))))
	)
}

# inverse penalised Fisher information at the final IRLS weights `Wobs`: covariance
# of the reported coefficients, joint with the factors and conditional on them
.ae_eb_vcov <- function(
	prep,
	fit,
	Wobs,
	ab,
	R,
	eta = NULL,
	link = NULL,
	vb = NULL,
	extra = NULL,
	evidence_only = FALSE
) {
	nr <- prep$nr; nc <- prep$nc; p <- prep$p_dyad; Tt <- prep$Tt
	Warr <- Wobs; Warr[!is.finite(Warr)] <- 0
	# dyad-level sums of the weights and the design (the weights stand in for the response)
	st <- ae_bcd_stats_cpp(Warr, Warr, prep$X, nr, nc, Tt, p, .ae_threads())
	wsum <- rowSums(Warr, dims = 2L)
	ctx <- .ae_eb_context(ab, st, wsum, p)
	# the cross weights of the pairs of a fit with a dyadic correlation
	if (!is.null(vb$cross)) {
		ctx <- .ae_pair_add_cross(ctx, .ae_pair_cross(list(nr = nr, Tt = Tt, p = p, X = prep[["X"]]), vb$cross))
	}
	H0 <- .ae_eb_penalise(ctx, .ae_eb_precision(fit$sab, ab))
	m0 <- ctx$m
	pa <- ctx$pa; pg <- ctx$pg
	# the evidence alone needs the log determinants, not the inverse (.ae_gauss_evidence)
	if (isTRUE(evidence_only)) {
		Hf <- if (R > 0L && !is.null(fit$U)) .ae_eb_extend(ctx, ab, H0, fit, vb) else H0
		iz <- setdiff(seq_len(nrow(Hf)), c(ctx$i_th, ctx$i_g))
		evidence <- tryCatch(
			.ae_eb_evidence_core(
				prep,
				ab,
				fit,
				R,
				eta,
				link,
				2 * sum(log(diag(chol(Hf[iz, iz, drop = FALSE])))),
				2 * sum(log(diag(chol(H0[ctx$i_s, ctx$i_s, drop = FALSE]))))
			),
			error = function(e) NULL
		)
		return(list(vcov = NULL, vcov_cond = NULL, evidence = evidence, elbo = NULL))
	}
	cond <- .ae_eb_inverse(H0)
	joint <- cond
	Hf <- H0
	if (R > 0L && !is.null(fit$U)) {
		Hf <- .ae_eb_extend(ctx, ab, H0, fit, vb)
		joint <- .ae_eb_inverse(Hf)[seq_len(m0), seq_len(m0), drop = FALSE]
	}
	# the evidence, from the information already formed (a variational fit brings its bound)
	evidence <- NULL
	if (!is.null(vb)) {
		evidence <- vb$evidence
	} else if (!is.null(eta) && !is.null(link)) {
		iz <- setdiff(seq_len(nrow(Hf)), c(ctx$i_th, ctx$i_g))
		evidence <- tryCatch(
			.ae_eb_evidence_core(
				prep,
				ab,
				fit,
				R,
				eta,
				link,
				2 * sum(log(diag(chol(Hf[iz, iz, drop = FALSE])))),
				2 * sum(log(diag(chol(H0[ctx$i_s, ctx$i_s, drop = FALSE]))))
			),
			error = function(e) NULL
		)
	}
	idx <- c(ctx$i_th, ctx$i_g)
	# add the uncertainty of a hyperparameter the coefficients respond to
	# (J = d coefficients / d log scale)
	if (!is.null(extra) && length(extra$J) == length(idx) && is.finite(extra$var) && extra$var > 0) {
		joint[idx, idx] <- joint[idx, idx] + extra$var * tcrossprod(extra$J)
	}
	# reported parametrisation: the intercept refers to uncentred node covariates
	A <- diag(length(idx))
	if (pg > 0L) {
		wbar <- c(ab$wbar_r, ab$wbar_c)
		# a symmetric dyad carries the node part of both of its members
		A[1L, pa + seq_len(pg)] <- -(if (identical(ab$kind, "symmetric")) 2 else 1) * wbar
	}
	nm <- c(
		"intercept",
		prep$x_names_dyad,
		if (!is.null(prep$W_row)) colnames(prep$W_row),
		if (!is.null(prep$W_col)) colnames(prep$W_col)
	)
	wrap <- function(S) {
		V <- A %*% S[idx, idx, drop = FALSE] %*% t(A)
		V <- (V + t(V)) / 2
		dimnames(V) <- list(nm, nm)
		V
	}
	list(
		vcov = wrap(joint),
		vcov_cond = wrap(cond),
		evidence = evidence,
		elbo = if (!is.null(vb)) vb$elbo
	)
}
