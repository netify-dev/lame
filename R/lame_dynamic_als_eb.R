# empirical-Bayes estimator for dynamic binary ALS fits (bias_correct = "eb"),
# the dynamic counterpart of ame_als_eb.R. Poisson panels use the
# Poisson-lognormal cell (`ctx$pln`, lame_dynamic_als_eb_pln.R) and normal
# panels the exact Gaussian cell (`ctx$nrm`, lame_dynamic_als_eb_nrm.R).
#
# the additive effects are random AR(1) paths around actor levels, as in lame():
#     a_it = alpha_i + d_it,  d_it = rho d_i,t-1 + e_it,  e_it ~ N(0, sigma^2),
# with (alpha_i, gamma_i) ~ N(0, Sab). paths, levels and coefficients solve one
# block-tridiagonal system with a border (.lde_factor()); Sab, sigma and rho take
# EM updates from posterior second moments under the hyperprior of lame();
# coefficients solve the Laplace-marginal estimating equations, with the Laplace
# covariance as their standard errors. a fit that does not converge is refused
# with a classed error (.lde_refuse_unfinished()).

# the AR(1) deviation law

# diagonal of Q, the vector Q 1 and the sum 1' Q 1 for a stationary AR(1) of
# length Tt >= 2 with persistence rho
.lde_ar_coefs <- function(rho, Tt) {
	cdiag <- rep(1 + rho^2, Tt)
	cdiag[1L] <- 1
	cdiag[Tt] <- 1
	f <- rep((1 - rho)^2, Tt)
	f[1L] <- 1 - rho
	f[Tt] <- 1 - rho
	list(c = cdiag, f = f, g = sum(f))
}

# the product Q d for a vector or the columns of a matrix d (Tt rows)
.lde_q_times <- function(d, rho) {
	d <- as.matrix(d)
	Tt <- nrow(d)
	out <- d
	if (Tt > 1L) {
		out[] <- 0
		out[1L, ] <- d[1L, ] - rho * d[2L, ]
		out[Tt, ] <- d[Tt, ] - rho * d[Tt - 1L, ]
		if (Tt > 2L) {
			mid <- 2:(Tt - 1L)
			out[mid, ] <- (1 + rho^2) * d[mid, , drop = FALSE] -
				rho * (d[mid - 1L, , drop = FALSE] + d[mid + 1L, , drop = FALSE])
		}
	}
	out
}

# the structured solver

# Cholesky with a growing ridge on failure; the blocks are positive definite,
# so a failure is rounding in a near-degenerate block
.lde_chol <- function(S, unit = 1) {
	ch <- tryCatch(chol(S), error = function(e) NULL)
	if (!is.null(ch)) return(ch)
	S <- (S + t(S)) / 2
	# ridge relative to the largest diagonal entry, floored at the unit cell information
	ridge <- 1e-10 * max(abs(diag(S)), unit)
	for (k in 1:12) {
		ch <- tryCatch(chol(S + diag(ridge, nrow(S))), error = function(e) NULL)
		if (!is.null(ch)) return(ch)
		ridge <- ridge * 100
	}
	cli::cli_abort(
		"A dynamic empirical-Bayes information block is not positive definite.",
		class = "lame_als_numerical"
	)
}

# information of a unit-variance cell: 1, or 1 / var(y) for a normal response
.lde_cunit <- function(ctx) {
	v <- ctx$prior$vunit
	if (is.numeric(v) && length(v) == 1L && is.finite(v) && v > 0) 1 / v else 1
}

# forward elimination of the block-tridiagonal period information with a border:
#   A  : the Tt diagonal blocks (prior included)
#   b  : the AR(1) coupling of neighbouring periods, -diag(b)
#   C  : the Tt blocks coupling each period with the border (levels,
#        coefficients, static factors)
#   Hzz: the information of the border
# Schur complements S_1 = A_1, S_t = A_t - (b b') * S_{t-1}^-1 (elementwise),
# and S_z = Hzz - sum_t C_t' Y_t with Y = H_ss^-1 C
.lde_forward <- function(A, b, unit = 1) {
	Tt <- length(A)
	bb <- tcrossprod(b)
	Si <- vector("list", Tt)
	logdet <- 0
	prev <- NULL
	for (t in seq_len(Tt)) {
		S <- A[[t]]
		if (t > 1L) S <- S - bb * prev
		ch <- .lde_chol(S, unit)
		logdet <- logdet + 2 * sum(log(diag(ch)))
		prev <- chol2inv(ch)
		Si[[t]] <- prev
	}
	list(Si = Si, logdet = logdet)
}

# solve H_ss Y = C for a list C of Tt blocks (dim(s_t) x r each)
.lde_xx_solve <- function(Si, C, b) {
	Tt <- length(Si)
	w <- vector("list", Tt)
	w[[1L]] <- C[[1L]]
	if (Tt > 1L) for (t in 2:Tt) w[[t]] <- C[[t]] + b * (Si[[t - 1L]] %*% w[[t - 1L]])
	Y <- vector("list", Tt)
	Y[[Tt]] <- Si[[Tt]] %*% w[[Tt]]
	if (Tt > 1L) for (t in (Tt - 1L):1L) {
		Y[[t]] <- Si[[t]] %*% (w[[t]] + b * Y[[t + 1L]])
	}
	Y
}

# factorise the information matrix: `qx` path rows lead each period block and
# `nl` level columns lead the border; `cdiag` holds the diagonal path-level
# coupling of each C_t, which saves those products
.lde_factor <- function(A, b, C, Hzz, qx, nl, cdiag = NULL, unit = 1) {
	Tt <- length(A)
	fw <- .lde_forward(A, b, unit)
	Y <- .lde_xx_solve(fw$Si, C, b)
	Sz <- Hzz
	if (is.null(cdiag)) {
		for (t in seq_len(Tt)) Sz <- Sz - crossprod(C[[t]], Y[[t]])
	} else {
		il <- seq_len(nl)
		ir <- seq_len(ncol(Hzz))[-il]
		for (t in seq_len(Tt)) {
			Sz[il, ] <- Sz[il, , drop = FALSE] - cdiag[[t]] * Y[[t]]
			if (length(ir)) Sz[ir, ] <- Sz[ir, , drop = FALSE] - crossprod(C[[t]][, ir, drop = FALSE], Y[[t]])
		}
	}
	Sz <- (Sz + t(Sz)) / 2
	chz <- .lde_chol(Sz, unit)
	Szi <- chol2inv(chz)
	il <- seq_len(nl)
	chl <- .lde_chol(Sz[il, il, drop = FALSE], unit)
	list(
		Tt = Tt,
		qx = qx,
		nl = nl,
		r = ncol(Hzz),
		b = b,
		unit = unit,
		Si = fw$Si,
		C = C,
		Y = Y,
		Sz = Sz,
		Szi = Szi,
		chz = chz,
		chl = chl,
		logdet_xx = fw$logdet,
		logdet_z = 2 * sum(log(diag(chz))),
		logdet_l = 2 * sum(log(diag(chl)))
	)
}

# posterior mode: rs the Tt period right-hand sides, rz the border's
.lde_solve <- function(fac, rs, rz) {
	Tt <- fac$Tt
	Cy <- lapply(rs, function(v) matrix(v, ncol = 1L))
	yh <- .lde_xx_solve(fac$Si, Cy, fac$b)
	rhs <- rz
	for (t in seq_len(Tt)) rhs <- rhs - as.vector(crossprod(fac$C[[t]], yh[[t]]))
	z <- as.vector(fac$Szi %*% rhs)
	s <- lapply(seq_len(Tt), function(t) as.vector(yh[[t]]) - as.vector(fac$Y[[t]] %*% z))
	list(s = s, z = z)
}

# solve for the periods and the levels with the other border unknowns held at
# zero and right-hand side (rs, 0); used by the Laplace correction
.lde_solve_cond <- function(fac, rs) {
	Tt <- fac$Tt
	il <- seq_len(fac$nl)
	Cy <- lapply(rs, function(v) matrix(v, ncol = 1L))
	yh <- .lde_xx_solve(fac$Si, Cy, fac$b)
	rhs <- numeric(fac$nl)
	for (t in seq_len(Tt)) rhs <- rhs - as.vector(crossprod(fac$C[[t]][, il, drop = FALSE], yh[[t]]))
	chl <- if (is.null(fac$chl)) .lde_chol(
		fac$Sz[il, il, drop = FALSE],
		if (is.null(fac$unit)) 1 else fac$unit
	) else fac$chl
	zl <- as.vector(backsolve(chl, backsolve(chl, rhs, transpose = TRUE)))
	s <- lapply(seq_len(Tt), function(t) {
		as.vector(yh[[t]]) - as.vector(fac$Y[[t]][, il, drop = FALSE] %*% zl)
	})
	list(s = s, lev = zl)
}

# diagonal blocks of the inverse period information (Omega_tt) and the diagonals
# of its first off-diagonal blocks, by backward recursion with Z = diag(b) Si_t
.lde_omega <- function(fac, lag = TRUE) {
	Tt <- fac$Tt; b <- fac$b
	Om <- vector("list", Tt)
	Ol <- vector("list", max(Tt - 1L, 0L))
	Om[[Tt]] <- fac$Si[[Tt]]
	if (Tt > 1L) for (t in (Tt - 1L):1L) {
		Z <- b * fac$Si[[t]]
		if (lag) Ol[[t]] <- colSums(Z * Om[[t + 1L]])
		Om[[t]] <- fac$Si[[t]] + crossprod(Z, Om[[t + 1L]] %*% Z)
	}
	list(diag = Om, lag = Ol)
}

# posterior moments: marginal variances and lag-one covariances of each period's
# unknowns, the path-level and border covariances (skipped with moments = FALSE),
# and with `cond` the covariance of the paths with the levels integrated out
.lde_posterior <- function(fac, cond = TRUE, full = FALSE, moments = TRUE) {
	Tt <- fac$Tt; qx <- fac$qx
	ix <- seq_len(qx)
	il <- seq_len(fac$nl)
	om <- .lde_omega(fac, lag = moments)
	qs <- nrow(om$diag[[1L]])
	out <- list()
	if (moments) {
		# X_t' = R^-T Y_t' with Sz = R'R, so that Y_t Szi Y_s' = X_t X_s'
		Xt <- lapply(seq_len(Tt), function(t) backsolve(fac$chz, t(fac$Y[[t]]), transpose = TRUE))
		out$var_s <- matrix(0, qs, Tt)
		out$cov_xl <- matrix(0, qx, Tt)
		out$cov_lag <- matrix(0, qs, max(Tt - 1L, 0L))
		for (t in seq_len(Tt)) {
			out$var_s[, t] <- diag(om$diag[[t]]) + colSums(Xt[[t]]^2)
			out$cov_xl[, t] <- -rowSums(fac$Y[[t]][ix, , drop = FALSE] * t(fac$Szi[, il[ix], drop = FALSE]))
		}
		if (Tt > 1L) for (t in seq_len(Tt - 1L)) {
			out$cov_lag[, t] <- om$lag[[t]] + colSums(Xt[[t]] * Xt[[t + 1L]])
		}
		out$var_x <- out$var_s[ix, , drop = FALSE]
		out$cov_zz <- fac$Szi
	}
	out$om <- om
	if (cond) {
		# Y_l Sl^-1 Y_l' = X X' with X' = R_l^-T Y_l'
		chl <- fac$chl
		out$cond <- lapply(seq_len(Tt), function(t) {
			Xl <- backsolve(chl, t(fac$Y[[t]][ix, il, drop = FALSE]), transpose = TRUE)
			om$diag[[t]][ix, ix, drop = FALSE] + crossprod(Xl)
		})
		if (isTRUE(full)) {
			out$cond_full <- lapply(seq_len(Tt), function(t) {
				Xl <- backsolve(chl, t(fac$Y[[t]][, il, drop = FALSE]), transpose = TRUE)
				om$diag[[t]] + crossprod(Xl)
			})
		}
	}
	out
}

# a dense reference (tests)

# the full information matrix built densely, in the order (s_1, ..., s_Tt, z)
.lde_dense_matrix <- function(A, b, C, Hzz) {
	Tt <- length(A); qs <- nrow(A[[1L]]); r <- ncol(Hzz)
	n <- Tt * qs + r
	H <- matrix(0, n, n)
	blk <- function(t) (t - 1L) * qs + seq_len(qs)
	zi <- Tt * qs + seq_len(r)
	for (t in seq_len(Tt)) {
		H[blk(t), blk(t)] <- A[[t]]
		if (t < Tt) {
			H[blk(t), blk(t + 1L)] <- diag(-b, qs)
			H[blk(t + 1L), blk(t)] <- diag(-b, qs)
		}
		H[blk(t), zi] <- C[[t]]
		H[zi, blk(t)] <- t(C[[t]])
	}
	H[zi, zi] <- Hzz
	H
}

# hyperprior and state

# the binary additive-effects prior of lame() (probit-moment scale and degrees
# of freedom from the first slice); kept separate from the static estimator
.lde_hyperprior <- function(Yarr, bip) {
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

# the unit of a variance: var(y) for a normal panel, 1 for families with a
# fixed latent scale; variance floors are fractions of it
.lde_vunit <- function(prep, nrm) {
	if (!isTRUE(nrm)) return(1)
	v <- stats::var(c(prep$Yarr), na.rm = TRUE)
	if (is.finite(v) && v > 0) v else 1
}

# hyperparameter priors of the dynamic additive effects as lame() sets them;
# `prior` overrides the defaults
.lde_prior <- function(prep, kind, prior = list(), rvar = TRUE, cvar = TRUE, pln = FALSE, nrm = FALSE) {
	prior <- if (is.list(prior)) prior else list()
	hp <- if (isTRUE(pln)) .lde_hyperprior_pois(prep$Yarr) else if (isTRUE(nrm)) .lde_hyperprior_nrm(prep$Yarr) else .lde_hyperprior(prep$Yarr, prep$bip)
	vunit <- .lde_vunit(prep, nrm)
	Sab0 <- prior$Sab0
	if (is.null(Sab0)) Sab0 <- hp$Sab0
	Sab0 <- as.matrix(Sab0)
	eta0 <- prior$etaab
	if (is.null(eta0)) eta0 <- prior$eta0
	if (is.null(eta0)) eta0 <- hp$eta0
	rho_mean <- .lda_default(prior$rho_ab_mean, 0.8)
	rho_sd <- .lda_default(prior$rho_ab_sd, 0.15)
	shape <- .lda_default(prior$sigma_ab_shape, 2)
	scale <- prior$sigma_ab_scale
	if (is.null(scale)) {
		# pooled level variance times (1 - rho^2), so the stationary deviation
		# variance matches the level variance a priori
		mult <- if (is.finite(eta0) && eta0 > 3) eta0 / (eta0 - 3) else 1
		v <- if (identical(kind, "symmetric")) {
			if (isTRUE(rvar)) Sab0[1L, 1L]
		} else c(if (isTRUE(rvar)) Sab0[1L, 1L], if (isTRUE(cvar)) Sab0[2L, 2L])
		v <- v[is.finite(v) & v > 0]
		scale <- if (!length(v)) vunit else max(mult * mean(v) * (1 - min(rho_mean^2, 0.9801)), 1e-6 * vunit)
	}
	out <- list(
		Sab0 = Sab0,
		eta0 = eta0,
		rho_mean = rho_mean,
		rho_sd = rho_sd,
		shape = shape,
		scale = scale,
		rho_lo = -0.99,
		rho_hi = 0.99,
		vunit = vunit,
		# shift making the EM fixed point of Sab the inverse-Wishart posterior mean
		shift = -((if (identical(kind, "directed")) 2 else 1) + 1)
	)
	# the sampler's prior of the count-cell noise variance: one degree of freedom on prior$s20
	if (isTRUE(pln)) out$ve_prior <- c(nu0 = 1, s20 = if (is.null(prior$s20)) 1 else .lame_validate_s20(prior$s20))
	# a normal residual variance: the same inverse-gamma on the sampler's default scale
	if (isTRUE(nrm)) out$ve_prior <- c(nu0 = 1, s20 = .lame_residual_prior_scale(prep$Yarr, "normal", prior$s20))
	out
}

# the working problem of a binary cell

# probit score, unfloored Fisher weight, its first and second derivatives and
# the log likelihood of complete cells
.lde_probit_raw <- function(yy, ee) {
	lp1 <- stats::pnorm(ee, log.p = TRUE)
	lp0 <- stats::pnorm(ee, lower.tail = FALSE, log.p = TRUE)
	ld <- stats::dnorm(ee, log = TRUE)
	lam1 <- exp(ld - lp1)
	lam0 <- exp(ld - lp0)
	hprime <- lam1 * lam0 * (lam0 - lam1 - 2 * ee)
	mm <- lam0 - lam1 - 2 * ee
	list(
		score = yy * lam1 - (1 - yy) * lam0,
		wraw = lam1 * lam0,
		hprime = hprime,
		hpp = hprime * mm + lam1 * lam0 * (lam0^2 + lam1^2 + ee * (lam1 - lam0) - 2),
		ell = ifelse(yy == 1, lp1, lp0)
	)
}

# score, floored Fisher weight w, working response z = eta + s / w and the
# derivatives h' and h'' of the weight (for the Laplace correction) of observed
# cells, zero elsewhere. the working shift is capped at delta_cap, with
# w = |score| / delta_cap there so w * shift stays the score; constants match
# .ae_irls_score_step()
.lde_work <- function(y, eta, link, delta_cap = 30, wfloor = 1e-8) {
	obs <- is.finite(y)
	yy <- y; yy[!obs] <- 0
	ee <- eta; ee[!obs] <- 0
	if (identical(link, "probit")) {
		pr <- .lde_probit_raw(yy, ee)
		score <- pr$score
		w <- pmax(pr$wraw, wfloor)
		hprime <- pr$hprime
		hpp <- pr$hpp
		dev <- -2 * pr$ell
	} else {
		p <- stats::plogis(ee)
		qq <- stats::plogis(-ee)
		score <- yy * qq - (1 - yy) * p
		w <- pmax(p * qq, wfloor)
		hprime <- p * qq * (1 - 2 * p)
		hpp <- p * qq * ((qq - p)^2 - 2 * p * qq)
		dev <- 2 * (log1p(exp(-abs(ee))) + pmax(-(2 * yy - 1) * ee, 0))
	}
	delta <- score / w
	big <- abs(delta) > delta_cap
	if (any(big)) {
		delta[big] <- sign(delta[big]) * delta_cap
		w[big] <- abs(score[big]) / delta_cap
	}
	z <- ee + delta
	z[!obs] <- 0
	w[!obs] <- 0
	score[!obs] <- 0
	hprime[!obs | !is.finite(hprime)] <- 0
	hpp[!obs | !is.finite(hpp)] <- 0
	dev[!obs] <- 0
	list(z = z, w = w, score = score, hprime = hprime, hpp = hpp, dev = sum(dev))
}

# the problem

# prior precision of the stacked levels (alpha, gamma): a directed pair is
# correlated, the two sides of a bipartite panel are not, a symmetric network
# has one effect per actor
.lde_level_precision <- function(sab, kind, na, nb) {
	if (identical(kind, "symmetric")) return(diag(1 / sab[1L, 1L], na))
	if (identical(kind, "bipartite")) {
		return(diag(c(rep(1 / sab[1L, 1L], na), rep(1 / sab[2L, 2L], nb))))
	}
	P <- solve(sab)
	rbind(
		cbind(P[1L, 1L] * diag(na), P[1L, 2L] * diag(na)),
		cbind(P[2L, 1L] * diag(na), P[2L, 2L] * diag(na))
	)
}

# the design of the active actors' cells: intercept, dyadic covariates and node
# covariates broadcast per period; a symmetric network adds both members' values
.lde_design <- function(prep, symmetric, ia, ib) {
	Tt <- prep$Tt
	na <- length(ia); nb <- length(ib)
	lay <- .lda_coef_layout(prep, include_nodes = TRUE)
	nm <- lay$names
	pall <- length(nm)
	X <- array(0, dim = c(na, nb, pall, Tt))
	X[, , 1L, ] <- 1
	pd <- prep$p_dyad
	if (pd > 0L) {
		for (k in seq_len(pd)) {
			xk <- prep$X[ia, ib, k, , drop = FALSE]
			xk[!is.finite(xk)] <- 0
			X[, , 1L + k, ] <- xk
		}
	}
	node_arr <- function(W, W_arr) {
		if (is.null(W) || !ncol(W)) return(NULL)
		if (!is.null(W_arr)) {
			A <- W_arr
		} else {
			A <- array(rep(W, Tt), dim = c(nrow(W), ncol(W), Tt))
		}
		A[!is.finite(A)] <- 0
		A
	}
	Wr <- node_arr(prep$W_row, prep$W_row_arr)
	Wc <- node_arr(prep$W_col, prep$W_col_arr)
	col <- 1L + pd
	one_a <- rep(1, na); one_b <- rep(1, nb)
	if (isTRUE(symmetric)) {
		for (W in list(Wr, Wc)) if (!is.null(W)) for (k in seq_len(dim(W)[2L])) {
			col <- col + 1L
			for (t in seq_len(Tt)) {
				X[, , col, t] <- outer(W[ia, k, t], one_b) + outer(one_a, W[ib, k, t])
			}
		}
	} else {
		if (!is.null(Wr)) for (k in seq_len(dim(Wr)[2L])) {
			col <- col + 1L
			for (t in seq_len(Tt)) X[, , col, t] <- outer(Wr[ia, k, t], one_b)
		}
		if (!is.null(Wc)) for (k in seq_len(dim(Wc)[2L])) {
			col <- col + 1L
			for (t in seq_len(Tt)) X[, , col, t] <- outer(one_a, Wc[ib, k, t])
		}
	}
	list(X = X, names = nm, block = lay$block)
}

# layout of the coefficient unknowns: static coefficients first, then for each
# moving coefficient its states of periods 0 to Tt (1 to Tt under rw2), grouped
# by block. `thidx[[t]]` gives each coefficient's position in period t and
# `pos_all(i)` all states of moving coefficient i
.lde_beta_layout <- function(names, block, mask, kind, Tt) {
	pall <- length(names)
	if (is.null(mask) || !any(mask)) {
		return(list(
			any = FALSE,
			mth = pall,
			ps = pall,
			pd = 0L,
			static_idx = seq_len(pall),
			dyn_idx = integer(0),
			thidx = rep(list(seq_len(pall)), Tt),
			pos0 = integer(0),
			groups = character(0),
			grp = integer(0),
			ng = 0L,
			kind = NULL
		))
	}
	static_idx <- which(!mask); dyn_idx <- which(mask)
	ps <- length(static_idx); pd <- length(dyn_idx)
	first <- if (identical(kind, "rw2")) 1L else 0L
	nst <- Tt + 1L - first
	posd <- function(i, t) ps + (i - 1L) * nst + (t - first) + 1L
	pos_all <- function(i) ps + (i - 1L) * nst + seq_len(nst)
	thidx <- lapply(seq_len(Tt), function(t) {
		v <- integer(pall)
		v[static_idx] <- seq_len(ps)
		v[dyn_idx] <- vapply(seq_len(pd), function(i) posd(i, t), integer(1))
		v
	})
	groups <- unique(block[dyn_idx])
	list(
		any = TRUE,
		mth = ps + pd * nst,
		ps = ps,
		pd = pd,
		static_idx = static_idx,
		dyn_idx = dyn_idx,
		thidx = thidx,
		pos0 = if (first == 0L) vapply(seq_len(pd), function(i) posd(i, 0L), integer(1)) else integer(0),
		posd = posd,
		pos_all = pos_all,
		groups = groups,
		grp = match(block[dyn_idx], groups),
		ng = length(groups),
		kind = kind
	)
}

# precision of the coefficient unknowns: the g-prior of the static ones and the
# AR(1) prior of lame() on moving states (or the sampler's rw2 precision)
.lde_theta_prior <- function(ctx, hyp) {
	lay <- ctx$beta
	P <- ctx$Lstat
	if (!lay$any) return(P)
	Tt <- ctx$Tt
	K2 <- if (identical(lay$kind, "rw2")) .rw2_precision(Tt)
	for (i in seq_len(lay$pd)) {
		g <- lay$grp[i]
		rho <- hyp$beta$rho[g]
		tau <- 1 / (hyp$beta$s2[g] * lay$lam[i])
		pos <- lay$pos_all(i)
		if (!is.null(K2)) {
			P[pos, pos] <- P[pos, pos] + K2 * tau
			next
		}
		d <- c(1 / lay$v0[i] + rho^2 * tau, rep((1 + rho^2) * tau, max(Tt - 1L, 0L)), tau)
		P[cbind(pos, pos)] <- P[cbind(pos, pos)] + d
		P[cbind(pos[-1L], pos[-length(pos)])] <- P[cbind(pos[-1L], pos[-length(pos)])] - rho * tau
		P[cbind(pos[-length(pos)], pos[-1L])] <- P[cbind(pos[-length(pos)], pos[-1L])] - rho * tau
	}
	P
}

# the context of a fit: everything fixed across iterations.
#   coef_prior : "mcmc" places the Zellner g-prior of lame() on the coefficients; "flat" none
#   R          : latent rank; fmode "fac" (static factors) or "dynamic" (moving positions)
#   beta_mask  : the moving design columns; beta_kind "ar1", "rw1" or "rw2"
.lde_context <- function(
	prep,
	symmetric,
	link = "probit",
	prior = list(),
	coef_prior = c("mcmc", "flat"),
	centre = TRUE,
	R = 0L,
	beta_mask = NULL,
	beta_kind = "ar1",
	dynamic_ab = TRUE,
	rvar = TRUE,
	cvar = TRUE,
	beta_pool = "none",
	dynamic_uv = FALSE,
	dynamic_uv_kind = "ar1",
	dcor = FALSE
) {
	coef_prior <- match.arg(coef_prior)
	bip <- isTRUE(prep$bip)
	# count panels use the Poisson-lognormal cell and normal panels the Gaussian
	# cell, each with the noise variance as a hyperparameter
	pln <- identical(link, "log")
	nrm <- identical(link, "identity")
	vehyp <- pln || nrm
	# a constant normal response has no variation to model
	if (nrm) {
		vy <- stats::var(c(prep$Yarr), na.rm = TRUE)
		if (!is.finite(vy) || vy <= 0) {
			.lde_refuse(
				"a constant normal response", "The observed responses are all equal (or there are fewer than two), so there is no variation to model and the residual variance is not identified."
			)
		}
	}
	kind <- if (isTRUE(symmetric)) "symmetric" else if (bip) "bipartite" else "directed"
	Tt <- prep$Tt
	Y <- prep$Yarr
	if (kind == "symmetric") {
		for (t in seq_len(Tt)) Y[, , t] <- .ae_symmetrize(Y[, , t])
	}
	obs <- is.finite(Y)
	row_any <- apply(obs, 1L, any)
	col_any <- apply(obs, 2L, any)
	if (bip) {
		ia <- which(row_any); ib <- which(col_any)
	} else {
		ia <- ib <- which(row_any | col_any)
	}
	na <- length(ia); nb <- length(ib)
	sym <- identical(kind, "symmetric")
	q <- if (sym) na else na + nb
	# a disabled margin (rvar / cvar) is a point mass at zero, as in the sampler: a
	# huge prior precision and no EM update. `fzm` flags the margins, `fz` the effects
	fzm <- if (sym) rep(!isTRUE(rvar), 2L) else c(!isTRUE(rvar), !isTRUE(cvar))
	fz <- if (sym) rep(fzm[1L], q) else c(rep(fzm[1L], na), rep(fzm[2L], nb))
	R <- as.integer(R)
	# "dynamic": moving positions; "fac": static factors with the variational
	# posterior of the static estimator (lame_dynamic_als_eb_fac.R)
	fmode <- if (R == 0L) "none" else if (isTRUE(dynamic_uv)) "dynamic" else "fac"
	des <- .lde_design(prep, sym, ia, ib)
	pall <- dim(des$X)[3L]
	if (identical(beta_kind, "rw2") && !is.null(beta_mask) && any(beta_mask) && Tt < 3L) {
		cli::cli_abort("{.arg dynamic_beta_kind = \"rw2\"} requires at least 3 time periods.")
	}
	lay <- .lde_beta_layout(des$names, des$block, beta_mask, beta_kind, Tt)
	# pooled moving-coefficient hyperparameters, as dynamic_beta_pool in the
	# sampler: a shared gamma scale for the innovations, a logit-normal law for the
	# persistences (ar1 only)
	pool <- list(
		sigma = lay$any && lay$ng >= 2L && beta_pool %in% c("sigma", "both"),
		rho = lay$any && lay$ng >= 2L && beta_pool %in% c("rho", "both") &&
			identical(lay$kind, "ar1")
	)
	yred <- Y[ia, ib, , drop = FALSE]
	obsr <- is.finite(yred)
	# dyadic correlation: a directed probit or normal network pairs the two cells
	# of a dyad in the iteration (lame_dynamic_als_eb_pair.R); a count panel
	# estimates it as a stage at the fitted state
	pair_stage <- pln && isTRUE(dcor) && identical(kind, "directed")
	dcor <- isTRUE(dcor) && identical(kind, "directed") && (identical(link, "probit") || nrm)
	lgy <- if (pln) { g <- yred; g[!obsr] <- 0; lgamma(g + 1) }
	# cells-by-covariates matrix of every period (NA design entries are zero)
	Xc <- lapply(seq_len(Tt), function(t) matrix(des$X[, , , t], na * nb, pall))
	centre_v <- numeric(pall)
	# no centring with moving coefficients, whose prior refers to the raw column
	if (isTRUE(centre) && pall > 1L && !lay$any) {
		cells <- do.call(rbind, lapply(seq_len(Tt), function(t) Xc[[t]][as.vector(obsr[, , t]), , drop = FALSE]))
		centre_v[-1L] <- colMeans(cells)[-1L]
		Xc <- lapply(Xc, function(M) sweep(M, 2L, centre_v, "-"))
	}
	# the Zellner g-prior of lame(): precision sum_t X_t' X_t / (T g) on the
	# observed cells, g = n^2 (times var(y) for a normal response)
	XtX <- matrix(0, pall, pall)
	ms <- numeric(pall); nobs <- 0
	for (t in seq_len(Tt)) {
		o <- as.vector(obsr[, , t])
		XtX <- XtX + crossprod(Xc[[t]][o, , drop = FALSE])
		nobs <- nobs + sum(o)
	}
	ms <- diag(XtX) / max(nobs, 1)
	unit_var <- if (nrm) .lame_unit_var(prep$Yarr, "normal") else 1
	gscale <- unit_var
	Lth <- if (identical(coef_prior, "mcmc")) XtX / (Tt * prep$nr^2 * gscale) else matrix(0, pall, pall)
	Lstat <- matrix(0, lay$mth, lay$mth)
	Lstat[seq_len(lay$ps), seq_len(lay$ps)] <- Lth[lay$static_idx, lay$static_idx, drop = FALSE]
	if (lay$any) {
		# initial-state variance per unit of the covariate: the response variance
		# over the mean square of the column, unless prior$beta0_var is given
		v0 <- if (is.null(prior$beta0_var)) .lame_beta0_var_default(
			prep$Yarr,
			if (nrm) "normal" else "binary"
		) else prior$beta0_var
		lay$v0 <- if (is.null(prior$beta0_var)) v0 / ifelse(
			ms[lay$dyn_idx] > 0,
			ms[lay$dyn_idx],
			1
		) else rep(v0, length(lay$dyn_idx))
		# innovation multiplier: one over the mean square of the column
		lay$lam <- 1 / ifelse(ms[lay$dyn_idx] > 0, ms[lay$dyn_idx], 1)
	}
	shape_b <- if (is.null(prior$sigma_beta_shape)) 2 else prior$sigma_beta_shape
	# the sampler's default innovation scale (.lame_sigma_beta_prior_scale())
	scale_b <- if (!is.null(prior$sigma_beta_scale)) prior$sigma_beta_scale else if (pln) {
		.lame_sigma_beta_prior_scale(prep$Yarr, "poisson", shape_b)
	} else if (nrm) {
		.lame_sigma_beta_prior_scale(prep$Yarr, "normal", shape_b)
	} else max((shape_b - 1) * 0.04, 1e-8)
	bp <- list(
		rho_mean = if (is.null(prior$rho_beta_mean)) 0.8 else prior$rho_beta_mean,
		rho_sd = if (is.null(prior$rho_beta_sd)) 0.15 else prior$rho_beta_sd,
		shape = shape_b,
		scale = scale_b,
		rho_lo = if (is.null(prior$rho_beta_lower)) 0 else prior$rho_beta_lower,
		rho_hi = if (is.null(prior$rho_beta_upper)) 0.999 else prior$rho_beta_upper,
		init_s2 = (if (is.null(prior$sigma_beta_init)) 0.25 * sqrt(unit_var) else prior$sigma_beta_init)^2,
		vunit = .lde_vunit(prep, nrm)
	)
	list(
		kind = kind,
		sym = sym,
		bip = bip,
		nu = if (sym) 0.5 else 1,
		static_ab = !isTRUE(dynamic_ab),
		static_s2 = 1e-4 * .lde_vunit(prep, nrm),
		nr = prep$nr,
		nc = prep$nc,
		na = na,
		nb = nb,
		q = q,
		Tt = Tt,
		ia = ia,
		ib = ib,
		y = yred,
		obs = obsr,
		Xc = Xc,
		centre = centre_v,
		pall = pall,
		mth = lay$mth,
		beta = lay,
		thidx = lay$thidx,
		names = des$names,
		block = des$block,
		link = link,
		pln = pln,
		nrm = nrm,
		vehyp = vehyp,
		lgy = lgy,
		prior = .lde_prior(prep, kind, prior, rvar = !fzm[1L], cvar = !fzm[2L], pln = pln, nrm = nrm),
		fzm = fzm,
		fz = fz,
		pm = ifelse(fz, 1e8, 1),
		Lth = Lth,
		Lstat = Lstat,
		coef_prior = coef_prior,
		betaprior = bp,
		pool = pool,
		beta_pool = beta_pool,
		R = R,
		fmode = fmode,
		uvp = if (fmode == "dynamic") .lde_pos_prior(
			prior,
			gaussian_bip_ar = nrm && bip && identical(dynamic_uv_kind, "ar1"),
			# positions carry the response unit here: there is no G to carry it
			scale_default = if (nrm && (bip || sym)) .lame_unit_sd(prep$Yarr, "normal") else .lame_sigma_uv_scale_default(
				prep$Yarr,
				if (nrm) "normal" else "binary",
				bip = bip,
				symmetric = sym
			)
		),
		# moving-position transition: Gaussian AR(1) ("ar1") or Student-t with
		# t_nu degrees of freedom ("t", default 4)
		uv_kind = if (fmode == "dynamic") dynamic_uv_kind else "ar1",
		t_nu = if (is.null(prior$t_nu)) 4 else prior$t_nu,
		# scale of the factor covariance prior (R/uv_prior.R)
		suv0 = .lame_uv_prior_scale(prep$Yarr, if (nrm) "normal" else if (pln) "poisson" else "binary"),
		# shared prior of the static factors of a directed panel
		uv_static = if (fmode == "fac" && identical(kind, "directed")) .lame_uv_prior(
			prep$Yarr,
			if (nrm) "normal" else if (pln) "poisson" else "binary",
			R
		),
		dcor = dcor,
		pair_stage = pair_stage,
		pc = if (dcor || pair_stage) .lde_pair_context(yred),
		# rule for the scale of the static factors' prior (options(lame.eb_scale)):
		# "cavity" (default) or "plugin"; hold_scale keeps it fixed
		fscale = if (fmode == "fac") .ae_mf_scale_rule(R) else "plugin",
		hold_scale = FALSE,
		row_names = prep$row_names,
		col_names = prep$col_names,
		time_names = prep$time_names
	)
}

# the level covariance with disabled margins at a point mass, as in the sampler
.lde_freeze_sab <- function(ctx, sab) {
	for (k in 1:2) if (ctx$fzm[k]) {
		sab[k, ] <- 0; sab[, k] <- 0; sab[k, k] <- 1e-8 * ctx$prior$vunit
	}
	sab
}

# initial state: pooled intercept, no effects, hyperparameters at prior means
.lde_init <- function(ctx) {
	ybar <- mean(ctx$y[ctx$obs])
	ybar <- min(max(ybar, 1e-3), 1 - 1e-3)
	theta <- numeric(ctx$mth)
	b0 <- if (isTRUE(ctx$pln)) .lde_pln_start_intercept(ctx) else if (isTRUE(ctx$nrm)) .lde_nrm_start_intercept(ctx)
	else if (identical(ctx$link, "probit")) stats::qnorm(ybar) else stats::qlogis(ybar)
	lay <- ctx$beta
	if (!lay$any) {
		theta[1L] <- b0
	} else {
		# the intercept is the first design column: static or moving
		if (1L %in% lay$static_idx) theta[match(1L, lay$static_idx)] <- b0
		else theta[lay$pos_all(match(1L, lay$dyn_idx))] <- b0
	}
	pr <- ctx$prior
	sab <- .lde_freeze_sab(ctx, if (!identical(ctx$kind, "directed")) diag(diag(pr$Sab0), 2L) else pr$Sab0)
	hyp <- list(sab = sab, s2 = pr$scale, rho = pr$rho_mean)
	# static additive effects: deviation paths held at a negligible scale (sd 0.01)
	# so one layout serves both
	if (isTRUE(ctx$static_ab)) hyp <- list(sab = sab, s2 = ctx$static_s2, rho = 0)
	# the dyadic correlation starts at zero
	if (isTRUE(ctx$dcor)) hyp$rho_d <- 0
	# noise variance starts: the static count estimator's, or half var(y) for normal
	if (isTRUE(ctx$pln)) hyp$ve <- .ae_pln_ve_start
	if (isTRUE(ctx$nrm)) hyp$ve <- .lde_nrm_start_ve(ctx)
	if (lay$any) {
		bp <- ctx$betaprior
		hyp$beta <- list(
			rho = rep(if (lay$kind %in% c("rw1", "rw2")) 1 else bp$rho_mean, lay$ng),
			s2 = rep(bp$scale / max(bp$shape - 1, 1), lay$ng)
		)
		# pooled hyperparameters start where the sampler starts them
		if (ctx$pool$sigma) hyp$beta$share <- bp$scale
		if (ctx$pool$rho) { hyp$beta$mu <- 0; hyp$beta$tau2 <- 1 }
	}
	list(
		theta = theta,
		x = matrix(0, ctx$q, ctx$Tt),
		lev = numeric(ctx$q),
		O = NULL,
		U = NULL,
		V = NULL,
		hyp = hyp
	)
}

# the information of the static factors (for the marginal covariance)

# index of the (actor, column) entries of the stacked position vector
# (vec(U), vec(V)), column-major within each factor
.lde_widx <- function(ctx) {
	na <- ctx$na; nb <- ctx$nb; R <- ctx$R
	list(
		u = lapply(seq_len(R), function(r) (r - 1L) * na + seq_len(na)),
		v = lapply(seq_len(R), function(r) na * R + (r - 1L) * nb + seq_len(nb))
	)
}

# Gauss-Newton information blocks of one period's static factors in the directed
# layout: with the additive paths (Kxw) and among themselves (Kww); used by
# .lde_fac_joint()
.lde_wblocks_dir <- function(ctx, Wt, Ut, Vt) {
	na <- ctx$na; nb <- ctx$nb; R <- ctx$R; q <- na + nb
	qw <- (na + nb) * R
	ix <- .lde_widx(ctx)
	Kxw <- matrix(0, q, qw)
	Kww <- matrix(0, qw, qw)
	ia <- seq_len(na); ib <- na + seq_len(nb)
	for (r in seq_len(R)) {
		iu <- ix$u[[r]]; iv <- ix$v[[r]]
		# (a_i, u_i r), (a_i, v_j r), (b_j, u_i r), (b_j, v_j r)
		Kxw[cbind(ia, iu)] <- as.vector(Wt %*% Vt[, r])
		Kxw[ia, iv] <- Wt * Ut[, r]
		Kxw[ib, iu] <- t(Wt * rep(Vt[, r], each = na))
		Kxw[cbind(ib, iv)] <- as.vector(crossprod(Wt, Ut[, r]))
		for (s in seq_len(R)) {
			ius <- ix$u[[s]]; ivs <- ix$v[[s]]
			Kww[cbind(iu, ius)] <- as.vector(Wt %*% (Vt[, r] * Vt[, s]))
			Kww[cbind(iv, ivs)] <- as.vector(crossprod(Wt, Ut[, r] * Ut[, s]))
			# (u_i r, v_j s) = w_ij v_j r u_i s
			blk <- Wt * outer(Ut[, s], Vt[, r])
			Kww[iu, ivs] <- blk
			Kww[ivs, iu] <- t(blk)
		}
	}
	list(Kxw = Kxw, Kww = Kww)
}

# map from the free positions F of a symmetric network to (U, V) = (F, F S),
# S the fixed diagonal sign pattern
.lde_tie <- function(ctx) {
	rbind(diag(ctx$na * ctx$R), kronecker(diag(ctx$sgn, ctx$R), diag(ctx$na)))
}

# (M V, M' U) as one stacked position vector
.lde_wsum <- function(ctx, M, Ut, Vt) {
	v <- c(as.vector(M %*% Vt), as.vector(crossprod(M, Ut)))
	if (ctx$sym) v <- ctx$nu * as.vector(crossprod(.lde_tie(ctx), v))
	v
}

.lde_eta <- function(ctx, st) {
	na <- ctx$na; nb <- ctx$nb; Tt <- ctx$Tt
	out <- array(0, dim = c(na, nb, Tt))
	a <- st$x[seq_len(na), , drop = FALSE]
	b <- if (ctx$sym) a else st$x[na + seq_len(nb), , drop = FALSE]
	for (t in seq_len(Tt)) {
		out[, , t] <- matrix(ctx$Xc[[t]] %*% st$theta[ctx$thidx[[t]]], na, nb) +
			outer(a[, t], b[, t], "+")
	}
	if (!is.null(st$O)) out <- out + st$O
	out
}

# the penalised normal equations

# likelihood information of the stacked paths of one period
.lde_kblock <- function(Wt, ctx, Ct = NULL) {
	na <- ctx$na; nb <- ctx$nb
	rs <- rowSums(Wt); cs <- colSums(Wt)
	if (ctx$sym) {
		K <- ctx$nu * (Wt + t(Wt))
		K[cbind(seq_len(na), seq_len(na))] <- K[cbind(seq_len(na), seq_len(na))] +
			ctx$nu * (rs + cs)
		return(K)
	}
	K <- matrix(0, ctx$q, ctx$q)
	ia <- seq_len(na); ib <- na + seq_len(nb)
	K[ia, ib] <- Wt
	K[ib, ia] <- t(Wt)
	K[cbind(ia, ia)] <- rs
	K[cbind(ib, ib)] <- cs
	if (!is.null(Ct)) {
		# cross information of the two cells of a dyad under the dyadic correlation
		K[ia, ia] <- K[ia, ia] + Ct
		K[ib, ib] <- K[ib, ib] + Ct
		d <- rowSums(Ct)
		K[cbind(ia, ib)] <- K[cbind(ia, ib)] + d
		K[cbind(ib, ia)] <- K[cbind(ib, ia)] + d
	}
	K
}

# row and column sums of a cell matrix, in the stacked effect order
.lde_stack_sums <- function(M, ctx) {
	if (ctx$sym) return(ctx$nu * (rowSums(M) + colSums(M)))
	c(rowSums(M), colSums(M))
}

# information blocks, right-hand side and score pieces at the IRLS weights `wk`;
# the multiplicative term enters in Gauss-Newton form, so the working response
# is z + O0 with O0 the term at the iterate
.lde_system <- function(ctx, st, wk) {
	Tt <- ctx$Tt; q <- ctx$q; p <- ctx$pall; m <- ctx$mth; nu <- ctx$nu
	na <- ctx$na; nb <- ctx$nb
	s <- 1 / st$hyp$s2; rho <- st$hyp$rho
	ar <- .lde_ar_coefs(rho, Tt)
	iq <- seq_len(q)
	ip <- q + seq_len(m)
	r_b <- q + m
	K <- vector("list", Tt); G <- vector("list", Tt)
	A <- vector("list", Tt); C <- vector("list", Tt)
	Ax <- vector("list", Tt); Cx <- vector("list", Tt)
	cdiag <- vector("list", Tt)
	rs <- vector("list", Tt)
	rx <- matrix(0, q, Tt)
	Hth <- matrix(0, m, m)
	rth <- numeric(m)
	# dyadic correlation: pair cross weights enter the information and the
	# right-hand side, with the partner's predictor (lame_dynamic_als_eb_pair.R)
	cw <- wk$cw
	eta_a <- if (!is.null(cw)) .lde_eta(ctx, st)
	for (t in seq_len(Tt)) {
		Wt <- wk$w[, , t]
		Ct <- if (is.null(cw)) NULL else cw[, , t]
		K[[t]] <- .lde_kblock(Wt, ctx, Ct)
		Xt <- ctx$Xc[[t]]
		wv <- as.vector(Wt)
		XW <- Xt * wv
		idx <- ctx$thidx[[t]]
		Hth[idx, idx] <- Hth[idx, idx] + nu * crossprod(Xt, XW)
		Gt <- matrix(0, m, q)
		for (k in seq_len(p)) Gt[idx[k], ] <- .lde_stack_sums(matrix(XW[, k], na, nb), ctx)
		if (!is.null(Ct)) {
			# pair cross terms between the coefficients and the partner's effects
			Mk <- lapply(seq_len(p), function(k) matrix(Xt[, k], na, nb) * Ct)
			for (k in seq_len(p)) Gt[idx[k], ] <- Gt[idx[k], ] + c(colSums(Mk[[k]]), rowSums(Mk[[k]]))
			Hx <- matrix(0, p, p)
			for (k in seq_len(p)) for (l in seq_len(p)) Hx[k, l] <- sum(Mk[[k]] * t(matrix(Xt[, l], na, nb)))
			Hth[idx, idx] <- Hth[idx, idx] + nu * Hx
		}
		G[[t]] <- Gt
		Ot <- if (is.null(st$O)) 0 else st$O[, , t]
		res <- Wt * (wk$z[, , t] + Ot)
		if (!is.null(Ct)) res <- res + Ct * t(eta_a[, , t])
		rx[, t] <- .lde_stack_sums(res, ctx)
		rth[idx] <- rth[idx] + nu * as.vector(crossprod(Xt, as.vector(res)))
		At <- K[[t]]
		At[cbind(iq, iq)] <- At[cbind(iq, iq)] + s * ctx$pm * ar$c[t]
		Ax[[t]] <- At
		Cxt <- matrix(0, q, q)
		cdiag[[t]] <- rep_len(-s * ctx$pm * ar$f[t], q)
		Cxt[cbind(iq, iq)] <- cdiag[[t]]
		Cx[[t]] <- Cxt
		A[[t]] <- At
		Ct <- matrix(0, q, q + m)
		Ct[, iq] <- Cxt
		Ct[, ip] <- t(Gt)
		C[[t]] <- Ct
		rs[[t]] <- rx[, t]
	}
	Pth <- .lde_theta_prior(ctx, st$hyp)
	Hth <- Hth + Pth
	Lam <- .lde_level_precision(st$hyp$sab, ctx$kind, na, nb)
	Hzz <- matrix(0, r_b, r_b)
	Hzz[iq, iq] <- Lam
	Hzz[cbind(iq, iq)] <- Hzz[cbind(iq, iq)] + s * ctx$pm * ar$g
	Hzz[ip, ip] <- Hth
	b <- s * rho * ctx$pm
	Lamx <- matrix(0, q, q)
	Lamx[] <- Lam
	Hzx <- Lamx
	Hzx[cbind(iq, iq)] <- Hzx[cbind(iq, iq)] + s * ctx$pm * ar$g
	list(
		K = K,
		G = G,
		Hth = Hth,
		Pth = Pth,
		rx = rx,
		rth = rth,
		A = A,
		C = C,
		rs = rs,
		b = b,
		Hzz = Hzz,
		Lam = Lam,
		s = s,
		rho = rho,
		Ax = Ax,
		Cx = Cx,
		bx = s * rho * ctx$pm,
		Hzx = Hzx,
		cdiag = cdiag,
		ip = ip,
		iq = iq
	)
}

# Laplace-marginal correction of the coefficient equations: the derivative of
# 0.5 log|K + Lam| through the cell predictor and the mode of the paths. per cell
#     kappa = 0.5 (h' v - h (u_i + u_j)),
# v the posterior variance of the cell's additive predictor, h the weight, h' its
# derivative and u = (K + Lam)^-1 E'(h' v). `facc` factorises the additive system
# (paths first, levels leading the border); `cond` holds the conditional
# covariance blocks of the paths
.lde_laplace <- function(ctx, facc, cond, wk, st) {
	Tt <- ctx$Tt; na <- ctx$na; nb <- ctx$nb
	vq <- vector("list", Tt)
	vcell <- array(0, dim = c(na, nb, Tt))
	vdiag <- vector("list", Tt)
	# dyadic correlation: a pair's information is a 2 x 2 matrix, so h' v becomes
	#   T_m = I11,m v1 + I22,m v2 + 2 I12,m cv,   m = 1, 2 (the cell and its partner),
	# with cv the covariance of the two predictors (lame_dynamic_als_eb_pair.R)
	pr <- wk$pair
	Tcell <- array(0, dim = c(na, nb, Tt))
	cvpair <- if (!is.null(pr)) numeric(length(pr$one))
	for (t in seq_len(Tt)) {
		C <- cond[[t]]
		dg <- diag(C)
		vdiag[[t]] <- dg
		if (ctx$sym) {
			v <- outer(dg, dg, "+") + 2 * C
		} else {
			ia <- seq_len(na); ib <- na + seq_len(nb)
			v <- outer(dg[ia], dg[ib], "+") + 2 * C[ia, ib, drop = FALSE]
		}
		vcell[, , t] <- v
		Th <- wk$hprime[, , t] * v
		if (!is.null(pr) && (!is.null(pr$I11a) || isTRUE(ctx$nrm))) {
			ip <- which(pr$t == t)
			if (length(ip)) {
				i <- pr$i[ip]; j <- pr$j[ip]
				cv <- C[cbind(i, j)] + C[cbind(i, na + i)] + C[cbind(na + j, j)] + C[cbind(na + j, na + i)]
				cvpair[ip] <- cv
				# a Gaussian pair has constant information and needs no correction
				if (!is.null(pr$I11a)) {
					v1 <- v[cbind(i, j)]
					v2 <- v[cbind(j, i)]
					Th[cbind(i, j)] <- pr$I11a[ip] * v1 + pr$I22a[ip] * v2 + 2 * pr$I12a[ip] * cv
					Th[cbind(j, i)] <- pr$I11b[ip] * v1 + pr$I22b[ip] * v2 + 2 * pr$I12b[ip] * cv
				}
			}
		}
		Tcell[, , t] <- Th
		vq[[t]] <- .lde_stack_sums(Th, ctx)
	}
	us <- .lde_solve_cond(facc, vq)$s
	kappa <- array(0, dim = c(na, nb, Tt))
	ucell <- array(0, dim = c(na, nb, Tt))
	for (t in seq_len(Tt)) {
		ua <- us[[t]][seq_len(na)]
		ub <- if (ctx$sym) ua else us[[t]][na + seq_len(nb)]
		uc <- outer(ua, ub, "+")
		ucell[, , t] <- uc
		kt <- 0.5 * (Tcell[, , t] - wk$w[, , t] * uc)
		if (!is.null(wk$cw)) kt <- kt - 0.5 * wk$cw[, , t] * t(uc)
		kappa[, , t] <- kt
	}
	kappa[!ctx$obs] <- 0
	# Laplace-marginal score of the dyadic correlation and its curvature `rc`
	# (lame_dynamic_als_eb_pair.R)
	rho_score <- NULL
	if (!is.null(pr) && !is.null(pr$s0)) {
		Tr <- sum(pr$I11r * vcell[pr$one] + pr$I22r * vcell[pr$two] + 2 * pr$I12r * cvpair)
		shift <- sum(ucell[pr$one] * pr$psi1r) + sum(ucell[pr$two] * pr$psi2r)
		rho_score <- list(rs = sum(pr$s0) - 0.5 * Tr - 0.5 * shift, rc = sum(pr$rc0))
	}
	# derivative of the correction in the cell predictor: the loop gain of the
	# fixed-point map
	kprime <- array(0, dim = c(na, nb, Tt))
	for (t in seq_len(Tt)) {
		ua <- us[[t]][seq_len(na)]
		ub <- if (ctx$sym) ua else us[[t]][na + seq_len(nb)]
		kprime[, , t] <- 0.5 * (wk$hpp[, , t] * vcell[, , t] - wk$hprime[, , t] * outer(ua, ub, "+"))
	}
	kprime[!ctx$obs] <- 0
	corr <- numeric(ctx$mth)
	for (t in seq_len(Tt)) {
		idx <- ctx$thidx[[t]]
		corr[idx] <- corr[idx] + ctx$nu * as.vector(crossprod(ctx$Xc[[t]], as.vector(kappa[, , t])))
	}
	list(
		kappa = kappa,
		kprime = kprime,
		corr = corr,
		vcell = vcell,
		vdiag = vdiag,
		cvpair = cvpair,
		Tcell = Tcell,
		ucell = ucell,
		rho_score = rho_score
	)
}

# the hyperparameter EM updates

# expected sums of squares and lag-one products of paths with means d, variances
# Vd (both nser x Tt) and lag-one covariances Cd (nser x Tt-1)
.lde_ar_stats_core <- function(d, Vd, Cd, w = NULL) {
	nser <- nrow(d); Tt <- ncol(d)
	E1 <- d^2 + Vd
	E2 <- d[, -1L, drop = FALSE] * d[, -Tt, drop = FALSE] + Cd
	# scale-mixture weights multiply the transition sums, not the first period
	wt <- if (is.null(w)) 1 else w[, -1L, drop = FALSE]
	list(
		S1 = sum(E1[, 1L]),
		Sn = sum(wt * E1[, -1L, drop = FALSE]),
		Sl = sum(wt * E1[, -Tt, drop = FALSE]),
		Sc = sum(wt * E2),
		nser = nser,
		Tt = Tt
	)
}

# sufficient statistics of the deviations d = x - l of the additive paths from
# their marginal posterior moments
.lde_ar_stats <- function(x, lev, post, free = seq_len(nrow(x))) {
	q <- nrow(x); Tt <- ncol(x)
	d <- x - lev
	vl <- diag(post$cov_zz)[seq_len(q)]
	Vd <- post$var_s[seq_len(q), , drop = FALSE] + vl - 2 * post$cov_xl
	Cd <- post$cov_lag[seq_len(q), , drop = FALSE] - post$cov_xl[, -1L, drop = FALSE] -
		post$cov_xl[, -Tt, drop = FALSE] + vl
	.lde_ar_stats_core(d[free, , drop = FALSE], Vd[free, , drop = FALSE], Cd[free, , drop = FALSE])
}

# EM update of (sigma^2, rho): rho maximises the profiled expected log-posterior,
# then sigma^2 is the inverse-gamma posterior mean
.lde_update_ar <- function(S, pr) {
	cnt <- S$nser * S$Tt
	A <- function(rho) (1 - rho^2) * S$S1 + S$Sn - 2 * rho * S$Sc + rho^2 * S$Sl
	nstar <- cnt / 2 + pr$shape + 1
	prof <- function(rho) {
		a <- A(rho)
		if (!is.finite(a) || a < 0) return(-Inf)
		-nstar * log(a / 2 + pr$scale) + 0.5 * S$nser * log(1 - rho^2) -
			(rho - pr$rho_mean)^2 / (2 * pr$rho_sd^2)
	}
	grid <- seq(pr$rho_lo, pr$rho_hi, length.out = 81L)
	pv <- vapply(grid, prof, numeric(1))
	k <- which.max(pv)
	lo <- grid[max(k - 1L, 1L)]; hi <- grid[min(k + 1L, length(grid))]
	opt <- stats::optimize(prof, c(lo, hi), maximum = TRUE, tol = 1e-10)
	rho <- opt$maximum
	s2 <- (A(rho) / 2 + pr$scale) / (cnt / 2 + pr$shape - 1)
	list(rho = rho, s2 = max(s2, 1e-10 * pr$vunit))
}

# EM update of the level covariance (inverse-Wishart or inverse-gamma posterior
# mean from the levels' second moments)
.lde_update_sab <- function(ctx, lev, post) {
	na <- ctx$na; nb <- ctx$nb; pr <- ctx$prior
	cl <- post$cov_zz[seq_len(ctx$q), seq_len(ctx$q), drop = FALSE]
	ia <- seq_len(na)
	if (ctx$sym) {
		va <- (pr$eta0 * pr$Sab0[1L, 1L] + sum(lev[ia]^2 + diag(cl)[ia])) /
			(na + pr$eta0 + pr$shift)
		return(.lde_freeze_sab(ctx, diag(max(va, 1e-8 * pr$vunit), 2L)))
	}
	ib <- na + seq_len(nb)
	if (ctx$kind == "bipartite") {
		va <- (pr$eta0 * pr$Sab0[1L, 1L] + sum(lev[ia]^2 + diag(cl)[ia])) / (na + pr$eta0 + pr$shift)
		vb <- (pr$eta0 * pr$Sab0[2L, 2L] + sum(lev[ib]^2 + diag(cl)[ib])) / (nb + pr$eta0 + pr$shift)
		return(.lde_freeze_sab(ctx, diag(c(max(va, 1e-8 * pr$vunit), max(vb, 1e-8 * pr$vunit)))))
	}
	cab <- cl[cbind(ia, ib)]
	S2 <- matrix(
		c(
			sum(lev[ia]^2 + diag(cl)[ia]),
			sum(lev[ia] * lev[ib] + cab),
			sum(lev[ia] * lev[ib] + cab),
			sum(lev[ib]^2 + diag(cl)[ib])
		),
		2L,
		2L
	)
	new <- (pr$eta0 * pr$Sab0 + S2) / (na + pr$eta0 + pr$shift)
	new <- (new + t(new)) / 2
	d <- sqrt(pmax(diag(new), 1e-8 * pr$vunit))
	r <- max(min(new[1L, 2L] / (d[1L] * d[2L]), 0.99), -0.99)
	new[1L, 2L] <- new[2L, 1L] <- r * d[1L] * d[2L]
	diag(new) <- d^2
	.lde_freeze_sab(ctx, new)
}

# EM image of all hyperparameters of state `st` under posterior moments `post`
.lde_em <- function(ctx, st, post) {
	free <- which(!ctx$fz)
	ar <- if (isTRUE(ctx$static_ab) || !length(free)) list(s2 = st$hyp$s2, rho = st$hyp$rho) else
		.lde_update_ar(.lde_ar_stats(st$x, st$lev, post, free), ctx$prior)
	out <- list(sab = .lde_update_sab(ctx, st$lev, post), s2 = ar$s2, rho = ar$rho)
	# the count noise variance takes its EM image in the Laplace stage
	# (.lde_pln_ee_sum()); carry the incoming value
	if (!is.null(st$hyp$ve)) out$ve <- st$hyp$ve
	if (ctx$beta$any) out$beta <- .lde_update_beta_hyper(ctx, st, post)
	# the dyadic correlation takes a scoring step elsewhere (.lde_pair_rho_step())
	if (!is.null(st$hyp$rho_d)) out$rho_d <- st$hyp$rho_d
	out
}

# EM update of the innovation scale and persistence of the moving coefficients,
# one pair per block, from their states' second moments (scale only for random
# walks). with dynamic_beta_pool, as in the sampler: "sigma" shares a gamma
# distributed inverse-gamma scale; "rho" puts a normal law on the logit of the
# scaled persistences (mu ~ N(0, 10^2), tau^2 ~ IG(2, 0.5))
.lde_update_beta_hyper <- function(ctx, st, post) {
	lay <- ctx$beta; bp <- ctx$betaprior; Tt <- ctx$Tt
	q <- ctx$q
	h <- st$hyp$beta
	rho <- h$rho; s2 <- h$s2
	pool_s <- isTRUE(ctx$pool$sigma); pool_r <- isTRUE(ctx$pool$rho)
	share <- if (pool_s) h$share else bp$scale
	kind <- lay$kind
	K2 <- if (identical(kind, "rw2")) .rw2_precision(Tt)
	logit <- function(r) stats::qlogis((r - bp$rho_lo) / (bp$rho_hi - bp$rho_lo))
	Ag <- numeric(lay$ng); cg <- numeric(lay$ng)
	for (g in seq_len(lay$ng)) {
		members <- which(lay$grp == g)
		A0 <- 0; A1 <- 0; A2 <- 0; Aq <- 0
		for (i in members) {
			pos <- lay$pos_all(i)
			m <- st$theta[pos]
			Cv <- post$cov_zz[q + pos, q + pos, drop = FALSE]
			w <- 1 / lay$lam[i]
			if (!is.null(K2)) {
				# expected second-difference quadratic form (with the sampler's ridge)
				Aq <- Aq + w * (sum(m * (K2 %*% m)) + sum(K2 * Cv))
			} else {
				E2 <- m^2 + diag(Cv)
				lagc <- m[-1L] * m[-(Tt + 1L)] + Cv[cbind(2:(Tt + 1L), 1:Tt)]
				A0 <- A0 + w * sum(E2[-1L]) # sum_t E b_t^2, t = 1..Tt
				A1 <- A1 + w * sum(lagc) # sum_t E b_t b_{t-1}
				A2 <- A2 + w * sum(E2[-(Tt + 1L)]) # sum_t E b_{t-1}^2
			}
		}
		cnt <- length(members) * Tt
		A <- if (!is.null(K2)) function(r) Aq else function(r) A0 - 2 * r * A1 + r^2 * A2
		if (kind %in% c("rw1", "rw2")) {
			r_new <- 1
		} else {
			nstar <- cnt / 2 + bp$shape + 1
			prof <- function(r) {
				a <- A(r)
				if (!is.finite(a) || a < 0) return(-Inf)
				prior_term <- if (pool_r) {
					-(logit(r) - h$mu)^2 / (2 * h$tau2)
				} else -(r - bp$rho_mean)^2 / (2 * bp$rho_sd^2)
				v <- -nstar * log(a / 2 + share) + prior_term
				if (is.finite(v)) v else -Inf
			}
			grid <- seq(bp$rho_lo, bp$rho_hi, length.out = 81L)
			k <- which.max(vapply(grid, prof, numeric(1)))
			lo <- grid[max(k - 1L, 1L)]; hi <- grid[min(k + 1L, length(grid))]
			r_new <- stats::optimize(prof, c(lo, hi), maximum = TRUE, tol = 1e-10)$maximum
		}
		rho[g] <- r_new
		s2[g] <- max((A(r_new) / 2 + share) / (cnt / 2 + bp$shape - 1), 1e-10 * bp$vunit)
		Ag[g] <- A(r_new); cg[g] <- cnt
	}
	out <- list(rho = rho, s2 = s2)
	# mean-field update of the shared scale
	if (pool_s) {
		out$share <- (bp$shape + lay$ng * bp$shape) /
			(bp$scale + sum((bp$shape + cg / 2) / (share + Ag / 2)))
	}
	if (pool_r) {
		th <- logit(rho)
		prec <- 1 / 100 + lay$ng / h$tau2
		out$mu <- (sum(th) / h$tau2) / prec
		out$tau2 <- (0.5 + sum((th - out$mu)^2) / 2) / (2 + lay$ng / 2 - 1)
	}
	out
}

# largest relative move between two hyperparameter sets
.lde_hyper_gap <- function(ctx, h0, h1) {
	nd <- if (ctx$sym) 1L else 2L
	g <- max(
		abs(log(diag(h1$sab) / diag(h0$sab))[seq_len(nd)]),
		abs(log(h1$s2 / h0$s2)),
		abs(h1$rho - h0$rho)
	)
	if (identical(ctx$kind, "directed")) {
		g <- max(g, abs(h1$sab[1L, 2L] - h0$sab[1L, 2L]) / sqrt(prod(diag(h0$sab))))
	}
	if (!is.null(h0$ve) && !is.null(h1$ve)) g <- max(g, abs(log(h1$ve / h0$ve)))
	if (!is.null(h0$beta)) {
		g <- max(g, abs(log(h1$beta$s2 / h0$beta$s2)), abs(h1$beta$rho - h0$beta$rho))
		if (!is.null(h0$beta$share)) g <- max(g, abs(log(h1$beta$share / h0$beta$share)))
		if (!is.null(h0$beta$mu)) g <- max(g, abs(h1$beta$mu - h0$beta$mu), abs(log(h1$beta$tau2 / h0$beta$tau2)))
	}
	g
}

# one evaluation of a state

# penalised objective of a state on the deviance scale, Laplace-marginal in the
# paths and levels, with the hyperprior terms whose stationary points are the EM
# updates
.lde_objective <- function(ctx, st, wk, facc) {
	hyp <- st$hyp; pr <- ctx$prior; Tt <- ctx$Tt; q <- ctx$q
	s <- 1 / hyp$s2
	d <- st$x - st$lev
	Qd <- .lde_q_times(t(d), hyp$rho)
	pen_dev <- s * sum(t(Qd) * d)
	Lam <- .lde_level_precision(hyp$sab, ctx$kind, ctx$na, ctx$nb)
	pen_lev <- sum(st$lev * (Lam %*% st$lev))
	Pth <- .lde_theta_prior(ctx, hyp)
	pen_th <- sum(st$theta * (Pth %*% st$theta))
	if (ctx$beta$any) {
		bp <- ctx$betaprior; lay <- ctx$beta; hb <- hyp$beta
		for (i in seq_len(lay$pd)) {
			g <- lay$grp[i]
			v0_term <- if (identical(lay$kind, "rw2")) 0 else log(lay$v0[i])
			pen_th <- pen_th + Tt * log(hb$s2[g] * lay$lam[i]) + v0_term
		}
		share <- if (ctx$pool$sigma) hb$share else bp$scale
		rho_term <- if (!identical(lay$kind, "ar1")) 0 else if (ctx$pool$rho) {
			(stats::qlogis((hb$rho - bp$rho_lo) / (bp$rho_hi - bp$rho_lo)) - hb$mu)^2 / hb$tau2
		} else (hb$rho - bp$rho_mean)^2 / bp$rho_sd^2
		pen_th <- pen_th + sum(2 * (bp$shape - 1) * log(hb$s2) + 2 * share / hb$s2 + rho_term)
		# hyperprior of the pooled hyperparameters
		if (ctx$pool$sigma) {
			pen_th <- pen_th + 2 * bp$scale * hb$share - 2 * (bp$shape + lay$ng * bp$shape) * log(hb$share)
		}
		if (ctx$pool$rho) {
			pen_th <- pen_th + hb$mu^2 / 100 + lay$ng * log(hb$tau2) + 2 * (2 - 1) * log(hb$tau2) + 2 * 0.5 / hb$tau2
		}
	}
	nser <- sum(!ctx$fz)
	logdet_prior <- nser * (Tt * log(hyp$s2) - log(1 - hyp$rho^2))
	sab <- hyp$sab
	fm <- !ctx$fzm
	logdet_lev <- if (ctx$sym) {
		if (fm[1L]) ctx$na * log(sab[1L, 1L]) else 0
	} else if (ctx$kind == "bipartite" || !all(fm)) {
		(if (fm[1L]) ctx$na * log(sab[1L, 1L]) else 0) + (if (fm[2L]) ctx$nb * log(sab[2L, 2L]) else 0)
	} else ctx$na * log(det(sab))
	hyper_sab <- if (ctx$sym) {
		if (fm[1L]) (pr$eta0 + pr$shift) * log(sab[1L, 1L]) + pr$eta0 * pr$Sab0[1L, 1L] / sab[1L, 1L] else 0
	} else if (ctx$kind == "bipartite" || !all(fm)) {
		k <- which(fm)
		sum((pr$eta0 + pr$shift) * log(diag(sab)[k]) + pr$eta0 * diag(pr$Sab0)[k] / diag(sab)[k])
	} else {
		(pr$eta0 + pr$shift) * log(det(sab)) + pr$eta0 * sum(diag(pr$Sab0 %*% solve(sab)))
	}
	hyper_ar <- 2 * (pr$shape - 1) * log(hyp$s2) + 2 * pr$scale / hyp$s2 +
		(hyp$rho - pr$rho_mean)^2 / pr$rho_sd^2
	# a symmetric network counts each dyad once (nu), so its penalties and log
	# determinants are divided by nu
	if (isTRUE(ctx$static_ab)) {
		# the stiff deviation block contributes a constant: q Tt log(s)
		pen_dev <- 0; logdet_prior <- 0; hyper_ar <- 0
		ldx <- facc$logdet_xx + facc$logdet_l - q * Tt * log(s)
	} else {
		ldx <- facc$logdet_xx + facc$logdet_l
	}
	# the sampler's prior of the dyadic correlation (R/rho_prior.R)
	hyper_rho <- if (isTRUE(ctx$dcor)) -2 * .lame_rho_prior_log(hyp$rho_d) else 0
	# the sampler's inverse-gamma prior of the count noise variance
	hyper_ve <- if (!is.null(hyp$ve)) .lde_pln_ve_hyper(ctx, hyp$ve) else 0
	wk$dev + (pen_th + pen_dev + pen_lev + logdet_prior + logdet_lev +
			ldx + hyper_sab + hyper_ar + hyper_rho + hyper_ve) / ctx$nu
}

# penalised deviance of a state without log determinants: what a step of the
# reweighting loop must not raise grossly
.lde_pdev <- function(ctx, st, hyp = st$hyp, dev = NULL) {
	if (is.null(dev)) {
		eta <- .lde_eta(ctx, st)
		dev <- .lde_dev_plain(ctx, eta, hyp)
	}
	st$hyp <- hyp
	s <- 1 / hyp$s2
	d <- st$x - st$lev
	Qd <- .lde_q_times(t(d), hyp$rho)
	Lam <- .lde_level_precision(hyp$sab, ctx$kind, ctx$na, ctx$nb)
	pen <- sum(st$lev * (Lam %*% st$lev)) + sum(st$theta * (.lde_theta_prior(ctx, hyp) %*% st$theta))
	if (!isTRUE(ctx$static_ab)) pen <- pen + s * sum(t(Qd) * d)
	dev + pen / ctx$nu
}

# stationarity certificate of a state: the scaled score of every block and the
# hyperparameter move under one EM step from the state
.lde_certificate <- function(ctx, st, wk, sysm, lap, em_self, tol) {
	Tt <- ctx$Tt; q <- ctx$q; nu <- ctx$nu
	hyp <- st$hyp; s <- sysm$s; rho <- hyp$rho
	ar <- .lde_ar_coefs(rho, Tt)
	g_th <- numeric(ctx$mth)
	for (t in seq_len(Tt)) {
		idx <- ctx$thidx[[t]]
		g_th[idx] <- g_th[idx] + nu * as.vector(crossprod(
			ctx$Xc[[t]],
			as.vector(wk$score[, , t] - lap$kappa[, , t])
		))
	}
	g_th <- g_th - as.vector(sysm$Pth %*% st$theta)
	i_th <- diag(sysm$Hth)
	d <- st$x - st$lev
	Qd <- t(.lde_q_times(t(d), rho))
	g_x <- matrix(0, q, Tt)
	i_x <- matrix(0, q, Tt)
	for (t in seq_len(Tt)) {
		g_x[, t] <- .lde_stack_sums(wk$score[, , t], ctx) - s * ctx$pm * Qd[, t]
		i_x[, t] <- diag(sysm$K[[t]]) + s * ctx$pm * ar$c[t]
	}
	g_l <- -as.vector(sysm$Lam %*% st$lev) + s * ctx$pm * rowSums(Qd)
	i_l <- diag(sysm$Lam) + s * ctx$pm * ar$g
	sc <- function(g, i) sqrt(2) * abs(g) / sqrt(pmax(i, .Machine$double.xmin))
	by_block <- c(
		theta = max(sc(g_th, i_th)),
		paths = max(sc(g_x, i_x)),
		levels = max(sc(g_l, i_l))
	)
	gap <- if (is.null(em_self)) 0 else .lde_hyper_gap(ctx, hyp, em_self)
	# the dyadic correlation is held to a tenth of the tolerance
	if (!is.null(em_self) && !is.null(hyp$rho_d)) gap <- max(
		gap,
		abs(em_self$rho_d - hyp$rho_d) * sqrt(tol) / (10 * tol)
	)
	thr <- sqrt(tol)
	mx <- max(by_block)
	list(
		pass = all(is.finite(c(by_block, gap))) && mx <= thr && gap <= thr,
		max_score = mx,
		scale_change = gap,
		by_block = by_block,
		threshold = thr,
		merit = max(mx, gap)
	)
}

# the candidate state from a solve of the corrected equations
.lde_candidate <- function(ctx, st, sol, em) {
	q <- ctx$q; p <- ctx$mth
	cst <- st
	cst$x <- do.call(cbind, lapply(sol$s, function(v) v[seq_len(q)]))
	cst$lev <- sol$z[seq_len(q)]
	cst$theta <- sol$z[q + seq_len(p)]
	if (!is.null(em)) cst$hyp <- em
	cst
}

# step factorisation with a modified border information: the positive part of
# kappa' added to the coefficient weights damps the feedback of the Laplace
# correction, which otherwise makes the map cycle; the fixed point and the
# posterior moments are unchanged
.lde_step_factor <- function(ctx, st, fac, lap) {
	Tt <- ctx$Tt; q <- ctx$q; m <- ctx$mth; nu <- ctx$nu
	r <- fac$r
	ip <- q + seq_len(m)
	dH <- matrix(0, r, r)
	for (t in seq_len(Tt)) {
		kp <- pmax(lap$kprime[, , t], 0)
		Xt <- ctx$Xc[[t]]
		idx <- ctx$thidx[[t]]
		dH[ip[idx], ip[idx]] <- dH[ip[idx], ip[idx]] + nu * crossprod(Xt, Xt * as.vector(kp))
	}
	Sz <- fac$Sz + dH
	fac$Sz <- Sz
	fac$Szi <- chol2inv(.lde_chol((Sz + t(Sz)) / 2, if (is.null(fac$unit)) 1 else fac$unit))
	# the stored factors belong to the unmodified block
	fac$chz <- NULL
	fac$chl <- NULL
	fac$dH <- dH
	fac
}

# evaluate a state: working problem, factorised system, Laplace correction,
# certificate and the next state of the iteration
.lde_eval <- function(ctx, st, tol = 5e-5, candidate = TRUE, move_ve = TRUE) {
	q <- ctx$q
	eta <- .lde_eta(ctx, st)
	ve_gap <- NULL
	if (isTRUE(ctx$vehyp) && move_ve) {
		# count or normal panel: move the noise variance to its EM image first, so
		# the step solves the equations of the state it returns; the size of the
		# move is the EM residual the certificate reads
		ve_new <- .lde_ve_image_at(ctx, eta, st$vcell, st$cvpair, st$dshift, st$hyp)
		ve_gap <- abs(log(ve_new / st$hyp$ve))
		st$hyp$ve <- ve_new
	}
	wk <- .lde_work_plain(ctx, eta, st$hyp)
	sysm <- .lde_system(ctx, st, wk)
	fac <- .lde_factor(sysm$A, sysm$b, sysm$C, sysm$Hzz, q, q, sysm$cdiag, unit = .lde_cunit(ctx))
	post <- .lde_posterior(fac, cond = TRUE)
	lap <- .lde_laplace(ctx, fac, post$cond, wk, st)
	em_self <- .lde_em(ctx, st, post)
	# EM residual of a state evaluated where it stands
	dshift <- if (isTRUE(ctx$pln)) .lde_pln_dshift(ctx, st, eta, lap, fac)
	if (isTRUE(ctx$vehyp) && is.null(ve_gap)) ve_gap <- abs(
		log(.lde_ve_image_at(ctx, eta, lap$vcell, lap$cvpair, dshift, st$hyp) / st$hyp$ve)
	)
	# the dyadic correlation takes its scoring step at the state
	rstep <- .lde_pair_rho_step(ctx, st$hyp, lap, eta)
	rho_new <- rstep$rho
	if (!is.null(rho_new)) em_self$rho_d <- rho_new
	cert <- .lde_certificate(ctx, st, wk, sysm, lap, em_self, tol)
	if (isTRUE(ctx$vehyp)) cert <- .lde_cert_with_gap(cert, ve_gap)
	ev <- list(
		st = st,
		eta = eta,
		wk = wk,
		sysm = sysm,
		fac = fac,
		post = post,
		lap = lap,
		em_self = em_self,
		cert = cert,
		rho_step = rstep,
		objective = .lde_objective(ctx, st, wk, fac)
	)
	if (isTRUE(ctx$vehyp)) ev$ve_gap <- ve_gap
	if (candidate) {
		rz <- c(numeric(q), sysm$rth - lap$corr)
		facs <- .lde_step_factor(ctx, st, fac, lap)
		# the step solves the modified equations in absolute form
		z0 <- c(st$lev, st$theta)
		rz <- rz + as.vector(facs$dH %*% z0)
		sol <- .lde_solve(facs, sysm$rs, rz)
		cst <- .lde_candidate(ctx, st, sol, NULL)
		cst$hyp <- .lde_em(ctx, cst, post)
		if (!is.null(rho_new)) cst$hyp$rho_d <- rho_new
		if (isTRUE(ctx$vehyp)) { cst$vcell <- lap$vcell; cst$cvpair <- lap$cvpair; cst$dshift <- dshift }
		ev$cand <- cst
	}
	ev
}

# blend two states (a damped step); positions are blended by .lde_blend_dyn()
.lde_blend <- function(ctx, s0, s1, a) {
	out <- s1
	mix <- function(x0, x1) if (is.null(x0) || is.null(x1)) x1 else (1 - a) * x0 + a * x1
	out$theta <- mix(s0$theta, s1$theta)
	out$x <- mix(s0$x, s1$x)
	out$lev <- mix(s0$lev, s1$lev)
	h0 <- s0$hyp; h1 <- s1$hyp
	hy <- list(
		sab = (1 - a) * h0$sab + a * h1$sab,
		s2 = exp((1 - a) * log(h0$s2) + a * log(h1$s2)),
		rho = (1 - a) * h0$rho + a * h1$rho
	)
	if (!is.null(h1$rho_d)) hy$rho_d <- (1 - a) * h0$rho_d + a * h1$rho_d
	if (!is.null(h1$ve) && !is.null(h0$ve)) hy$ve <- exp((1 - a) * log(h0$ve) + a * log(h1$ve))
	if (!is.null(h1$beta)) {
		hy$beta <- if (!is.null(h0$beta)) {
			b <- list(
				s2 = exp((1 - a) * log(h0$beta$s2) + a * log(h1$beta$s2)),
				rho = (1 - a) * h0$beta$rho + a * h1$beta$rho
			)
			if (!is.null(h1$beta$share)) b$share <- exp((1 - a) * log(h0$beta$share) + a * log(h1$beta$share))
			if (!is.null(h1$beta$mu)) {
				b$mu <- (1 - a) * h0$beta$mu + a * h1$beta$mu
				b$tau2 <- exp((1 - a) * log(h0$beta$tau2) + a * log(h1$beta$tau2))
			}
			b
		} else h1$beta
	}
	out$hyp <- hy
	out
}

# gauge-free coordinates of a state, watched for oscillation of the iteration
.lde_hvec <- function(ctx, st) {
	h <- st$hyp
	v <- c(st$theta, st$lev, log(diag(h$sab)), log(h$s2), h$rho)
	if (!is.null(h$rho_d)) v <- c(v, h$rho_d)
	if (!is.null(h$ve)) v <- c(v, log(h$ve))
	if (identical(ctx$kind, "directed")) v <- c(v, h$sab[1L, 2L] / sqrt(prod(diag(h$sab))))
	if (!is.null(h$beta)) {
		v <- c(v, log(h$beta$s2), h$beta$rho)
		if (!is.null(h$beta$share)) v <- c(v, log(h$beta$share))
		if (!is.null(h$beta$mu)) v <- c(v, h$beta$mu, log(h$beta$tau2))
	}
	if (!is.null(h$uv)) v <- c(v, log(h$uv$s2), h$uv$rho)
	if (!is.null(h$fac)) v <- c(v, log(.lde_fac_scales(h$fac)))
	# moving-position means, so a 2-cycle of the position sweeps is seen
	if (!is.null(st$pos)) v <- c(v, as.vector(st$pos$qu$mean), as.vector(st$pos$qv$mean))
	# the mixing weights of a scale-mixture transition
	if (!is.null(st$pos$qu$lam)) v <- c(
		v,
		log(as.vector(st$pos$qu$lam[, -1L])),
		if (!ctx$sym) log(as.vector(st$pos$qv$lam[, -1L]))
	)
	v
}

# the reweighting loop: each iteration evaluates the state and moves a fraction
# omega towards the next state of the fixed-point map. a step that raises the
# penalised deviance grossly is damped at once. omega is halved while successive
# moves reverse (cosine below -0.5) and grown by 1.5 while they agree (cosine
# above 0.5), with a floor of 1/64
.lde_fit_core <- function(ctx, st0, max_iter = 200L, tol = 5e-5, verbose = FALSE) {
	st <- st0
	ev <- .lde_eval(ctx, st, tol, move_ve = !isTRUE(st$seeded))
	st <- ev$st
	trace <- data.frame(
		iter = 0L,
		merit = ev$cert$merit,
		score = ev$cert$max_score,
		gap = ev$cert$scale_change,
		objective = ev$objective,
		step = NA_real_
	)
	converged <- isTRUE(ev$cert$pass)
	pln <- isTRUE(ctx$vehyp)
	thr <- sqrt(tol)
	if (pln) {
		# a count fit is also certified by the distance of its noise variance to
		# the limit of its iteration (.lde_tail())
		r0 <- if (isTRUE(st$seeded) && isTRUE(is.finite(st$ve_rate))) min(st$ve_rate, 0.9999) else NA_real_
		converged <- converged && isTRUE(is.finite(r0)) && isTRUE(ev$ve_gap * r0 / (1 - r0) <= thr / 2)
	}
	termination <- if (converged) "converged" else "iteration limit"
	it <- 0L
	pdev_cur <- .lde_pdev(ctx, st)
	# refuse a step that raises the penalised deviance grossly at the incoming
	# hyperparameters
	gross <- function(cand, p_ref) {
		p <- .lde_pdev(ctx, cand, st$hyp)
		is.finite(p) && p <= p_ref + max(2, 0.01 * abs(p_ref))
	}
	omega <- 1
	v_old <- .lde_hvec(ctx, st)
	move_old <- NULL
	best <- NULL
	# log noise variance at accepted iterates, whose tail is extrapolated to its
	# limit (.lde_tail())
	ve_hist <- if (pln) log(st$hyp$ve) else numeric(0)
	ve_rate <- NA_real_
	jump_from <- 0L
	while (!converged && it < max_iter) {
		it <- it + 1L
		g <- ev$cand
		step <- omega
		cand <- if (omega < 1) .lde_blend(ctx, st, g, omega) else g
		if (!gross(cand, pdev_cur)) {
			accepted <- FALSE
			for (k in seq_len(14L)) {
				step <- omega * 2^(-k)
				bl <- .lde_blend(ctx, st, g, step)
				if (gross(bl, pdev_cur)) { cand <- bl; accepted <- TRUE; break }
			}
			if (!accepted) {
				termination <- "no improving step"
				break
			}
		}
		evc <- tryCatch(.lde_eval(ctx, cand, tol), error = function(e) NULL)
		if (is.null(evc) || !is.finite(evc$cert$merit) || !is.finite(evc$objective)) {
			termination <- "no improving step"
			break
		}
		st <- evc$st
		ev <- evc
		# the deviance of a count panel's state is the evaluation's own
		pdev_cur <- .lde_pdev(ctx, st, dev = if (isTRUE(ctx$vehyp)) ev$wk$dev)
		if (is.null(best) || ev$cert$merit < best$ev$cert$merit) best <- list(st = st, ev = ev)
		# relaxation: compare this move with the previous one
		v_new <- .lde_hvec(ctx, st)
		move <- v_new - v_old
		v_old <- v_new
		if (!is.null(move_old)) {
			nn <- sqrt(sum(move^2)) * sqrt(sum(move_old^2))
			if (is.finite(nn) && nn > 0) {
				cs <- sum(move * move_old) / nn
				if (cs < -0.5) omega <- max(omega / 2, 1 / 64)
				else if (cs > 0.5) omega <- min(1, omega * 1.5)
			}
		}
		move_old <- move
		trace <- rbind(
			trace,
			data.frame(
				iter = it,
				merit = ev$cert$merit,
				score = ev$cert$max_score,
				gap = ev$cert$scale_change,
				objective = ev$objective,
				step = step
			)
		)
		if (verbose) {
			iter_line <- sprintf(
				"iter %3d  merit %.3e  score %.3e  gap %.3e  obj %.4f  pdev %.3f  step %.3g",
				it,
				ev$cert$merit,
				ev$cert$max_score,
				ev$cert$scale_change,
				ev$objective,
				pdev_cur,
				step
			)
			cli::cli_verbatim(iter_line)
		}
		if (!pln) {
			if (isTRUE(ev$cert$pass)) { converged <- TRUE; termination <- "converged" }
		} else {
			ve_hist <- c(ve_hist, log(st$hyp$ve))
			tl <- .lde_tail(ve_hist)
			ve_rate <- tl$rate
			if (isTRUE(ev$cert$pass) && .lde_tail_done(tl, thr)) {
				converged <- TRUE
				termination <- "converged"
			} else if (it >= jump_from && !isTRUE(ctx$nrm) && tl$mono && isTRUE(tl$rate > 0.9) && tl$dist > thr / 4 && !isFALSE(
				getOption("lame.als_aitken", TRUE)
			)) {
				# slow monotone tail: jump the variance towards the limit of its moves (at
				# most a factor e, or half of it) and keep the state if its EM residual is
				# smaller than the move it follows; normal panels never jump
				jump <- sign(tl$move) * min(tl$dist, 1)
				took <- FALSE
				for (frac in c(1, 0.5)) {
					stj <- st
					stj$hyp$ve <- st$hyp$ve * exp(frac * jump)
					evj <- tryCatch(.lde_eval(ctx, stj, tol), error = function(e) NULL)
					if (!is.null(evj) && is.finite(evj$cert$merit) && is.finite(evj$objective) && isTRUE(evj$ve_gap < abs(tl$move))) {
						st <- evj$st
						ev <- evj
						pdev_cur <- .lde_pdev(ctx, st, dev = ev$wk$dev)
						v_old <- .lde_hvec(ctx, st)
						move_old <- NULL
						ve_hist <- log(st$hyp$ve)
						took <- TRUE
						break
					}
				}
				if (!took) jump_from <- it + 5L
			}
		}
	}
	if (pln) {
		st$ve_rate <- ve_rate
		ev$st$ve_rate <- ve_rate
	}
	# not certified: return the state with the smallest residual seen
	if (!converged && !is.null(best) && best$ev$cert$merit < ev$cert$merit) { st <- best$st; ev <- best$ev }
	list(
		st = st,
		ev = ev,
		converged = converged,
		iterations = it,
		termination = termination,
		trace = trace
	)
}

.lde_fit <- function(
	ctx,
	st0 = .lde_init(ctx),
	max_iter = 200L,
	tol = 5e-5,
	verbose = FALSE,
	jitter = 0,
	progress = FALSE
) {
	if (identical(ctx$fmode, "dynamic")) {
		# moving positions: fit without them, start the mean-field iteration from the
		# spectral start and search the position scale by profile (lame_dynamic_als_eb_pos.R)
		if (is.null(st0$pos)) {
			ctx0 <- ctx
			ctx0$R <- 0L; ctx0$fmode <- "none"
			# the start ignores the dyadic correlation, which starts at zero afterwards
			ctx0$dcor <- FALSE; ctx0$pc <- NULL
			fit0 <- .lde_fit_core(ctx0, st0, max_iter = max_iter, tol = max(tol, 1e-3))
			sp <- .lde_pos_start(ctx, fit0, jitter = jitter)
			if (ctx$sym) ctx$sgn <- sp$sgn
			st0 <- fit0$st
			st0$pos <- sp$pos
			st0$hyp$uv <- sp$hp
			st0$O <- .lde_pos_product(sp$pos, ctx$Tt)
			st0$Vpos <- .lde_pos_prodvar(sp$pos, ctx$R, ctx$Tt)
		}
		if (ctx$sym && !is.null(st0$sgn)) ctx$sgn <- st0$sgn
		if (isTRUE(st0$seeded)) {
			# a restart that is stationary at its fitted scale is returned as it is;
			# any other restart searches the scale afresh
			chk <- .lde_fit_core_dyn(ctx, st0, max_iter = 0L, tol = tol, settle = FALSE)
			if (chk$converged) {
				chk$ctx <- ctx
				if (is.data.frame(st0$search)) chk$profile <- list(grid = st0$search)
				return(chk)
			}
		}
		pr <- if (isTRUE(ctx$dcor) && !isTRUE(getOption("lame.dcor_full_search", FALSE))) {
			.lde_pos_profile_pair(ctx, st0, tol = tol, max_iter = max_iter, verbose = progress)
		} else {
			.lde_pos_profile(ctx, st0, tol = tol, max_iter = max_iter, verbose = progress)
		}
		out <- pr$fit
		out$profile <- pr
		out$ctx <- ctx
		return(out)
	}
	if (identical(ctx$fmode, "fac")) {
		# static factors: fit without them, then the variational iteration from the
		# spectral start (lame_dynamic_als_eb_fac.R), with prior scales set by the
		# static estimator's rule
		if (is.null(st0$pos)) {
			ctx0 <- ctx
			ctx0$R <- 0L; ctx0$fmode <- "none"
			ctx0$dcor <- FALSE; ctx0$pc <- NULL
			fit0 <- .lde_fit_core(ctx0, st0, max_iter = max_iter, tol = max(tol, 1e-3))
			sp <- .lde_fac_start(ctx, fit0, jitter = jitter)
			if (ctx$sym) ctx$sgn <- sp$sgn
			st0 <- fit0$st
			st0$pos <- sp$pos
			st0$hyp$fac <- sp$prior
			st0$O <- sp$O
			st0$Vpos <- sp$Vpos
		}
		if (ctx$sym && !is.null(st0$sgn)) ctx$sgn <- st0$sgn
		if (isTRUE(st0$seeded)) {
			chk <- .lde_fit_core_dyn(ctx, st0, max_iter = 0L, tol = tol, settle = FALSE)
			if (chk$converged) {
				chk$ctx <- ctx
				return(chk)
			}
		}
		out <- .lde_fit_core_dyn(ctx, st0, max_iter = max_iter, tol = tol)
		out$ctx <- ctx
		return(out)
	}
	out <- .lde_fit_core(ctx, st0, max_iter = max_iter, tol = tol, verbose = verbose)
	out$ctx <- ctx
	out
}

# refusals

# refuse a specification the estimator does not cover; bias_correct = "none"
# still fits the uncorrected penalised model
.lde_refuse <- function(what, why, call = parent.frame(), .envir = parent.frame()) {
	cli::cli_abort(
		c(
			paste0("The empirical-Bayes estimator of a dynamic ALS fit is not available for ", what, "."),
			"x" = why,
			"i" = "{.code bias_correct = \"none\"} fits the uncorrected penalised model (finite-sample bias of order 10 to 50 percent on small or sparse networks); {.code method = \"mcmc\"} fits the calibrated model."
		),
		class = c("lame_dynamic_eb_unsupported", "lame_als_eb"),
		call = call,
		.envir = .envir
	)
}

.lde_check_supported <- function(
	prep,
	symmetric,
	R_eff,
	R_row_eff,
	R_col_eff,
	dynamic_ab,
	dynamic_uv,
	dynamic_uv_kind,
	beta_dyn,
	dynamic_beta_kind,
	lambda_ab,
	lambda_beta,
	lambda_uv,
	prior,
	link
) {
	if (!is.null(lambda_ab) || !is.null(lambda_beta) || !is.null(lambda_uv) ||
			!is.null(prior$lambda_ab_als) || !is.null(prior$lambda_beta_als) ||
			!is.null(prior$lambda_uv_als)) {
		cli::cli_abort(
			c(
				"{.arg lambda_ab}, {.arg lambda_beta} and {.arg lambda_uv} fix the smoothness penalties of the uncorrected penalised fit.",
				"i" = "The empirical-Bayes estimator estimates the innovation scales and persistences of the paths by EM under the {.fn lame} hyperpriors ({.code prior$sigma_ab_shape}, {.code prior$sigma_ab_scale}, {.code prior$rho_ab_mean}, {.code prior$rho_ab_sd}); use {.code bias_correct = \"none\"} to fix the penalties."
			),
			class = c("lame_dynamic_eb_unsupported", "lame_als_eb")
		)
	}
	if (isTRUE(dynamic_uv) && !dynamic_uv_kind %in% c("ar1", "t")) {
		.lde_refuse(
			"{.code dynamic_uv_kind = \"{dynamic_uv_kind}\"} latent positions",
			"Only the Gaussian AR(1) and the Student-t transitions have a posterior-moment (variational) update validated against the sampler."
		)
	}
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") && isTRUE(prep$bip)) {
		.lde_refuse(
			"{.code dynamic_uv_kind = \"t\"} latent positions of a bipartite network",
			"The interaction U G V' of a bipartite network leaves the size of the positions to G, and the transition weights then read innovations that are not on the scale of the prior; the sampler does not fit it either."
		)
	}
	if (isTRUE(dynamic_uv) && !link %in% c("probit", "log", "identity")) {
		.lde_refuse(
			"moving latent positions with a {.val {link}} link",
			"The sampler's model has a probit link, and the cell problem under the position noise is built for it only."
		)
	}
	if (isTRUE(beta_dyn$any) && !dynamic_beta_kind %in% c("ar1", "rw1", "rw2")) {
		.lde_refuse(
			"{.code dynamic_beta_kind = \"{dynamic_beta_kind}\"}",
			"Only the {.val ar1}, {.val rw1} and {.val rw2} kinds have a posterior-moment (EM) update of the state scale (and persistence) validated against the sampler."
		)
	}
	invisible(TRUE)
}

# the fitted object

# coefficient path of every period (p x Tt); the intercept refers to the
# uncentred covariates
.lde_coef_path <- function(ctx, theta) {
	p <- ctx$pall; Tt <- ctx$Tt
	cp <- vapply(seq_len(Tt), function(t) theta[ctx$thidx[[t]]], numeric(p))
	cp <- matrix(cp, p, Tt)
	if (any(ctx$centre != 0)) cp[1L, ] <- cp[1L, ] - sum(ctx$centre * cp[, 1L])
	dimnames(cp) <- list(ctx$names, ctx$time_names)
	cp
}

# map from the coefficient unknowns to the reported time-averaged coefficients
.lde_avg_map <- function(ctx) {
	p <- ctx$pall; m <- ctx$mth; Tt <- ctx$Tt
	A <- matrix(0, p, m)
	for (t in seq_len(Tt)) {
		idx <- ctx$thidx[[t]]
		A[cbind(seq_len(p), idx)] <- A[cbind(seq_len(p), idx)] + 1 / Tt
	}
	if (any(ctx$centre != 0)) {
		# intercept of the uncentred covariates (static centring only)
		A[1L, ] <- A[1L, ] - colSums(ctx$centre[-1L] * A[-1L, , drop = FALSE])
	}
	A
}

.lde_report_theta <- function(ctx, theta) {
	th <- as.vector(.lde_avg_map(ctx) %*% theta)
	names(th) <- ctx$names
	th
}

# covariance of the reported coefficients from the covariance of the unknowns
.lde_report_vcov <- function(ctx, V) {
	A <- .lde_avg_map(ctx)
	out <- A %*% V %*% t(A)
	out <- (out + t(out)) / 2
	dimnames(out) <- list(ctx$names, ctx$names)
	out
}

# assemble the fit object: the fields of a dynamic ALS fit (see .lda_as_fit)
# with the empirical-Bayes record
.lde_as_fit <- function(
	prep,
	ctx,
	fit,
	call,
	call_args,
	mode,
	symmetric,
	link,
	bias_requested,
	tol,
	seed,
	R_row = NULL,
	R_col = NULL
) {
	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt
	na <- ctx$na; nb <- ctx$nb; q <- ctx$q; R <- ctx$R
	st <- fit$st; ev <- fit$ev
	sym <- ctx$sym
	# count panels: the noise variance is a fitted hyperparameter and the fitted
	# mean is exp(eta + ve / 2)
	pln <- isTRUE(ctx$pln)
	nrm <- isTRUE(ctx$nrm)
	ve_hat <- if (pln || nrm) st$hyp$ve else 1
	# a disabled margin is a point mass at zero
	if (any(ctx$fz)) { st$x[ctx$fz, ] <- 0; st$lev[ctx$fz] <- 0 }
	th <- .lde_report_theta(ctx, st$theta)
	p <- ctx$pall
	lay <- ctx$beta
	cpath <- .lde_coef_path(ctx, st$theta)
	Vc <- ev$fac$Szi[q + seq_len(ctx$mth), q + seq_len(ctx$mth), drop = FALSE]
	vc <- .lde_report_vcov(ctx, Vc)
	vcc <- vc
	fac_scale <- NULL
	# directed count panel: the dyadic correlation enters the coefficient
	# covariance only (lame_dynamic_als_eb_plnpair.R)
	pair <- fit$pair
	ctx_p <- ctx
	if (!is.null(pair)) ctx_p$pair_info <- list(rho = pair$rho, cvv = pair$cvv)
	if (!is.null(pair) && !identical(ctx$fmode, "fac")) {
		Vc <- .lde_pln_pair_vcov(ctx, st, ev, pair)
		vc <- .lde_report_vcov(ctx, Vc)
	}
	if (identical(ctx$fmode, "fac")) {
		# static factors: coefficient covariance marginal over the factors, with the
		# conditional one beside it
		fv <- .lde_fac_vcov(ctx, st, ev, ctx_p)
		Vc <- fv$vcov
		vc <- .lde_report_vcov(ctx, Vc)
		fac_scale <- fv$scale
	}
	embed_paths <- function(x) {
		a <- matrix(0, nr, Tt, dimnames = list(prep$row_names, prep$time_names))
		b <- matrix(0, nc, Tt, dimnames = list(prep$col_names, prep$time_names))
		a[ctx$ia, ] <- x[seq_len(na), , drop = FALSE]
		b[ctx$ib, ] <- if (sym) x[seq_len(na), , drop = FALSE] else x[na + seq_len(nb), , drop = FALSE]
		list(a = a, b = b)
	}
	embed_levels <- function(l) {
		a <- structure(numeric(nr), names = prep$row_names)
		b <- structure(numeric(nc), names = prep$col_names)
		a[ctx$ia] <- l[seq_len(na)]
		b[ctx$ib] <- if (sym) l[seq_len(na)] else l[na + seq_len(nb)]
		list(a = a, b = b)
	}
	# count panels report paths and levels at the first-order Laplace posterior
	# mean, mode - u / 2 (lame_dynamic_als_eb_pln.R); the state stays the mode
	x_rep <- st$x
	lev_rep <- st$lev
	if (pln && !is.null(ev$lap$Tcell)) {
		ms <- .lde_pln_mean_shift(ctx, ev)
		x_rep <- st$x - 0.5 * ms$x
		lev_rep <- st$lev - 0.5 * ms$lev
		if (any(ctx$fz)) { x_rep[ctx$fz, ] <- 0; lev_rep[ctx$fz] <- 0 }
	}
	ab <- embed_paths(x_rep)
	lv <- embed_levels(lev_rep)
	# effects reported at mean zero over the actors, with the matching intercept
	cen <- .lde_centre_report(ctx, ev, ab, lv, th, cpath, vc, vcc, st = st)
	ab <- cen$ab; lv <- cen$lv; th <- cen$th; cpath <- cen$cpath; vc <- cen$vc; vcc <- cen$vcc
	eb_centre <- cen$centre
	a_mat <- ab$a; b_mat <- ab$b
	# the multiplicative term on the full actor set
	Oarr <- array(0, dim = c(nr, nc, Tt), dimnames = list(prep$row_names, prep$col_names, prep$time_names))
	if (!is.null(st$O)) Oarr[ctx$ia, ctx$ib, ] <- st$O
	# remove the constant moved into the intercept; the predictor is unchanged
	Oarr <- Oarr - cen$o_shift
	# bipartite ranks may differ: U G V' has rank min(R_row, R_col), the larger
	# factor carries zero columns and G is the rectangular identity
	Rr <- if (prep$bip && !is.null(R_row)) R_row else R
	Rc <- if (prep$bip && !is.null(R_col)) R_col else R
	rn_dim <- if (Rr > 0L) paste0(if (prep$bip) "row_dim" else "dim", seq_len(Rr)) else character(0)
	cn_dim <- if (Rc > 0L) paste0(if (prep$bip) "col_dim" else "dim", seq_len(Rc)) else character(0)
	Uout <- Vout <- Lout <- U_sd <- V_sd <- lam_u <- lam_v <- NULL
	if (ctx$fmode == "fac" && sym) {
		# undirected networks report the eigendecomposition of F S F'
		eg <- eigen((st$O[, , 1L] + t(st$O[, , 1L])) / 2, symmetric = TRUE)
		idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
		Uout <- matrix(0, nr, R, dimnames = list(prep$row_names, rn_dim))
		Uout[ctx$ia, ] <- eg$vectors[, idx, drop = FALSE]
		Vout <- Uout
		Lout <- eg$values[idx]
	} else if (ctx$fmode == "fac") {
		# static factors: variational means in the principal-axis gauge of the product
		rp <- .lde_fac_report(ctx, st$pos)
		Uout <- matrix(0, nr, Rr, dimnames = list(prep$row_names, rn_dim))
		Vout <- matrix(0, nc, Rc, dimnames = list(prep$col_names, cn_dim))
		Uout[ctx$ia, seq_len(R)] <- rp$U; Vout[ctx$ib, seq_len(R)] <- rp$V
	} else if (ctx$fmode == "dynamic") {
		# moving positions: means in the gauge of the time-averaged product, with sds
		rp <- .lde_pos_report(ctx, st$pos)
		Uout <- array(0, dim = c(nr, Rr, Tt), dimnames = list(prep$row_names, rn_dim, prep$time_names))
		Vout <- array(0, dim = c(nc, Rc, Tt), dimnames = list(prep$col_names, cn_dim, prep$time_names))
		# unobserved actors keep the prior sd sqrt(S); bipartite padding columns are zero
		U_sd <- array(0, dim = dim(Uout), dimnames = dimnames(Uout))
		V_sd <- array(0, dim = dim(Vout), dimnames = dimnames(Vout))
		U_sd[, seq_len(R), ] <- sqrt(st$hyp$uv$S); V_sd[, seq_len(R), ] <- sqrt(st$hyp$uv$S)
		Uout[ctx$ia, seq_len(R), ] <- rp$U; Vout[ctx$ib, seq_len(R), ] <- rp$V
		U_sd[ctx$ia, seq_len(R), ] <- rp$U_sd; V_sd[ctx$ib, seq_len(R), ] <- rp$V_sd
		if (sym) Lout <- stats::setNames(ctx$sgn, rn_dim)
		if (identical(ctx$uv_kind, "t")) {
			lam_u <- matrix(1, nr, Tt, dimnames = list(prep$row_names, prep$time_names))
			lam_v <- matrix(1, nc, Tt, dimnames = list(prep$col_names, prep$time_names))
			lam_u[ctx$ia, ] <- st$pos$qu$lam
			lam_v[ctx$ib, ] <- st$pos$qv$lam
		}
	}
	Omat <- if (ctx$fmode == "none") matrix(
		0,
		nr,
		nc,
		dimnames = list(prep$row_names, prep$col_names)
	) else apply(Oarr, c(1, 2), mean)
	# linear predictor on the full actor set; unobserved actors keep the prior mean zero
	des <- .lde_design(prep, sym, seq_len(nr), seq_len(nc))
	EZ <- vector("list", Tt)
	fitted_list <- vector("list", Tt)
	resid_list <- vector("list", Tt)
	W <- array(0, dim = c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		xt <- matrix(des$X[, , , t], nr * nc, p) %*% cpath[, t]
		ez <- matrix(xt, nr, nc) + outer(a_mat[, t], b_mat[, t], "+") + Oarr[, , t]
		W[ctx$ia, ctx$ib, t] <- ev$wk$w[, , t]
		# a normal fit's residual is read from the data, not the working response
		zz <- if (nrm) ctx$y[, , t] else ev$wk$z[, , t]; zz[!ctx$obs[, , t]] <- NA_real_
		zt <- matrix(NA_real_, nr, nc)
		zt[ctx$ia, ctx$ib] <- zz
		fit_t <- if (nrm) ez else if (pln) exp(ez + ve_hat / 2) else if (identical(link, "probit")) stats::pnorm(ez) else stats::plogis(ez)
		if (!prep$bip && nr == nc) {
			diag(ez) <- NA_real_
			diag(fit_t) <- NA_real_
		}
		dimnames(ez) <- list(prep$row_names, prep$col_names)
		dimnames(fit_t) <- list(prep$row_names, prep$col_names)
		EZ[[t]] <- ez
		fitted_list[[t]] <- fit_t
		# a count fit's working residual is taken at the mode, where it was formed
		if (pln) {
			et <- matrix(NA_real_, nr, nc)
			et[ctx$ia, ctx$ib] <- ev$eta[, , t]
			r <- zt - et
		} else r <- zt - ez
		dimnames(r) <- list(prep$row_names, prep$col_names)
		resid_list[[t]] <- r
	}
	names(EZ) <- names(fitted_list) <- names(resid_list) <- prep$time_names
	dimnames(W) <- list(prep$row_names, prep$col_names, prep$time_names)
	a_bar <- .lda_presence_row_means(a_mat, prep$row_presence)
	b_bar <- .lda_presence_row_means(b_mat, prep$col_presence)
	names(a_bar) <- prep$row_names; names(b_bar) <- prep$col_names
	hyp <- st$hyp
	sab <- hyp$sab
	for (k in which(ctx$fzm)) { sab[k, ] <- 0; sab[, k] <- 0 }
	# no path law with static additive effects or with every margin off
	no_ab <- isTRUE(ctx$static_ab) || all(ctx$fz)
	rho_ab <- hyp$rho
	sigma_ab <- sqrt(hyp$s2)
	va <- sab[1L, 1L]
	vb <- if (sym) NA_real_ else sab[2L, 2L]
	cab <- if (sym || prep$bip) NA_real_ else sab[1L, 2L]
	# weighted dispersion of the working response (about 1 under a correct
	# specification), a diagnostic beside ve = 1
	obs_w <- ctx$obs & is.finite(ev$wk$z)
	ve_working <- sum((ev$wk$w * (ev$wk$z - ev$eta)^2)[obs_w]) / sum(ev$wk$w[obs_w])
	# descriptive dyadic correlation of the working residuals, pooled and by period
	rho <- 0
	rho_path <- NULL
	if (!prep$bip && !sym) {
		rho_path <- stats::setNames(rep(NA_real_, Tt), prep$time_names)
		e_ij <- e_ji <- numeric(0)
		for (t in seq_len(Tt)) {
			rt <- resid_list[[t]]
			ut <- upper.tri(rt)
			e_t <- rt[ut]; et_t <- t(rt)[ut]
			ok_t <- is.finite(e_t) & is.finite(et_t)
			if (sum(ok_t) > 2L) {
				rho_path[t] <- suppressWarnings(stats::cor(e_t[ok_t], et_t[ok_t]))
				if (!is.finite(rho_path[t])) rho_path[t] <- NA_real_
			}
			e_ij <- c(e_ij, e_t); e_ji <- c(e_ji, et_t)
		}
		ok <- is.finite(e_ij) & is.finite(e_ji)
		if (sum(ok) > 2L) {
			rho <- suppressWarnings(stats::cor(e_ij[ok], e_ji[ok]))
			if (!is.finite(rho)) rho <- 0
		}
	}
	# a normal fit with dcor reports the estimated pair correlation instead
	if (nrm && isTRUE(ctx$dcor)) rho <- unname(hyp$rho_d)
	VC <- c(va = va, cab = cab, vb = vb, rho = rho, ve = ve_hat)
	# effective degrees of freedom of the mean model, read by the Pearson
	# overdispersion test
	df_mean <- if ((pln || nrm) && !is.null(ev$lap$vcell)) sum((ev$wk$w * ev$lap$vcell)[ctx$obs]) + p else NA_real_
	# empirical variances of the level estimates, descriptive only
	act_a <- seq_len(nr) %in% ctx$ia; act_b <- seq_len(nc) %in% ctx$ib
	va_raw <- if (sum(act_a) > 1L) stats::var(a_bar[act_a]) else NA_real_
	vb_raw <- if (sym) NA_real_ else if (sum(act_b) > 1L) stats::var(b_bar[act_b]) else NA_real_
	coef_names <- ctx$names
	coef_path <- cpath
	dimnames(coef_path) <- list(coef_names, prep$time_names)
	BETA <- if (lay$any) {
		array(coef_path, dim = c(1L, p, Tt), dimnames = list("point", coef_names, prep$time_names))
	} else {
		matrix(th, nrow = 1L, dimnames = list("point", coef_names))
	}
	beta_mask <- stats::setNames(seq_len(p) %in% lay$dyn_idx, coef_names)
	beta_groups <- stats::setNames(rep("", p), coef_names)
	if (lay$any) beta_groups[lay$dyn_idx] <- lay$groups[lay$grp]
	hb <- st$hyp$beta
	se_path <- NULL
	if (lay$any) {
		dv <- sqrt(pmax(diag(Vc), 0))
		se_path <- vapply(seq_len(Tt), function(t) dv[ctx$thidx[[t]]], numeric(p))
		dimnames(se_path) <- list(coef_names, prep$time_names)
		if (!is.null(cen$d_var_path)) se_path[1L, ] <- sqrt(se_path[1L, ]^2 + cen$d_var_path)
	}
	trace <- fit$trace
	# objective after each reweighting step (the start is not a step)
	obj_hist <- if (length(trace$objective) > 1L) trace$objective[-1L] else trace$objective
	obj_diff <- diff(trace$objective)
	# deviance of the reported fit, at the mean position offset when positions are fitted
	dev_final <- if (!is.null(ev$dev)) ev$dev
	else if (ctx$fmode %in% c("dynamic", "fac")) .lde_dev_plain(ctx, ev$eta, st$hyp)
	else ev$wk$dev
	dyn_uv <- identical(ctx$fmode, "dynamic")
	uvh <- if (dyn_uv) st$hyp$uv
	out <- list(
		call = call,
		family = if (nrm) "normal" else if (pln) "poisson" else "binary",
		mode = mode,
		symmetric = symmetric,
		R = if (prep$bip) max(Rr, Rc) else R,
		R_row = if (prep$bip) Rr else NULL,
		R_col = if (prep$bip) Rc else NULL,
		longitudinal = TRUE,
		mu = th["intercept"],
		beta = th[-1L],
		coef_path = coef_path,
		beta_dynamic_per_t = coef_path,
		beta_dynamic_mask = beta_mask,
		beta_dynamic_groups = beta_groups,
		BETA = BETA,
		a = a_bar,
		b = b_bar,
		a_dynamic = a_mat,
		b_dynamic = b_mat,
		a_level = lv$a,
		b_level = lv$b,
		APM = a_bar,
		BPM = b_bar,
		U = Uout,
		V = Vout,
		L = Lout,
		G = if (prep$bip && R > 0L) {
			g <- .lda_rect_identity(Rr, Rc)
			dimnames(g) <- list(rn_dim, cn_dim)
			g
		} else NULL,
		Omat = Omat,
		Oarr = Oarr,
		UVPM = if (dyn_uv) Oarr else Omat,
		coefficients = th,
		VC = VC,
		s2 = ve_hat,
		ve_working = ve_working,
		vc_raw = c(va = va_raw, vb = vb_raw),
		vc_period = c(va = NA_real_, vb = NA_real_),
		vc_noise = c(va = NA_real_, vb = NA_real_),
		df_mean = df_mean,
		pois_latent_var = if (pln) ve_hat,
		EZ = EZ,
		fitted = fitted_list,
		YPM = fitted_list,
		residuals = resid_list,
		iterations = fit$iterations,
		converged = fit$converged,
		deviance = dev_final,
		objective = ev$objective,
		objective_label = if (dyn_uv) "empirical-Bayes profile objective (Laplace-marginal deviance, position KL, hyperprior)" else "empirical-Bayes penalized deviance (Laplace-marginal)",
		dev_history = obj_hist,
		objective_history = obj_hist,
		penalized_deviance = obj_hist,
		n_time = Tt,
		dims = list(n_row = nr, n_col = nc, p = p),
		row_names = prep$row_names,
		col_names = prep$col_names,
		x_names = coef_names[-1L],
		time_names = prep$time_names,
		Y = prep$Yarr,
		X = prep$X,
		W_row = prep$W_row,
		W_col = prep$W_col,
		W_row_arr = prep$W_row_arr,
		W_col_arr = prep$W_col_arr,
		changing_composition = isTRUE(prep$changing_composition),
		row_presence = prep$row_presence,
		col_presence = prep$col_presence,
		actor_presence = prep$actor_presence,
		dynamic_uv = dyn_uv,
		dynamic_uv_kind = if (dyn_uv) ctx$uv_kind else NULL,
		dynamic_ab = !no_ab,
		dynamic_beta_pool = ctx$beta_pool,
		dynamic_beta = lay$any,
		dynamic_rho = FALSE,
		rho_path = rho_path,
		residual_rho_path = rho_path,
		rho_path_role = if (is.null(rho_path)) NULL else "residual_diagnostic",
		rho_uv = if (dyn_uv) uvh$rho else NULL,
		sigma_uv = if (dyn_uv) sqrt(uvh$s2) else NULL,
		U_sd = U_sd,
		V_sd = V_sd,
		# expected mixing weights of a Student-t transition; small weights mark heavy steps
		lambda_u = lam_u,
		lambda_v = lam_v,
		rho_ab = if (no_ab) NULL else rho_ab,
		sigma_ab = if (no_ab) NULL else sigma_ab,
		rho_beta = if (lay$any) stats::setNames(hb$rho, lay$groups) else NULL,
		sigma_beta = if (lay$any) stats::setNames(sqrt(hb$s2), lay$groups) else NULL,
		lambda_uv_als = NULL,
		lambda_ab_als = NULL,
		lambda_beta_als = NULL,
		lowrank_method = if (dyn_uv) "empirical-Bayes mean-field Gaussian positions" else if (R > 0L) "empirical-Bayes joint Gauss-Newton" else NULL,
		non_normal_method = if (nrm) "none" else "irls",
		link = if (nrm) NULL else link,
		obs_weights = W,
		irls_trace = obj_hist,
		stability = NULL,
		dynamic = list(
			uv = dyn_uv,
			uv_kind = if (dyn_uv) ctx$uv_kind else NULL,
			ab = !no_ab,
			ab_kind = if (no_ab) NULL else "ar1_eb",
			beta = lay$any,
			beta_kind = if (lay$any) lay$kind else NULL,
			G = FALSE,
			snap_model_estimated = FALSE,
			t_model_estimated = FALSE,
			mode = mode
		),
		diagnostics = list(
			converged = fit$converged,
			objective = obj_hist,
			convergence_trace = trace,
			objective_label = if (dyn_uv) "empirical-Bayes profile objective (Laplace-marginal deviance, position KL, hyperprior)" else "empirical-Bayes penalized deviance (Laplace-marginal)",
			objective_scope = "reweighting steps; the empirical-Bayes equations carry the Laplace correction of the incoming iterate, so a step need not lower it",
			reweight_steps = NULL,
			penalized_deviance = obj_hist,
			reported_objective = ev$objective,
			reported_objective_label = if (dyn_uv) "empirical-Bayes profile objective" else "empirical-Bayes penalized deviance",
			objective_increases = 0L,
			objective_rises = if (length(obj_diff)) sum(
				obj_diff > 1e-8 * (abs(utils::head(trace$objective, -1L)) + 1)
			) else 0L,
			max_objective_increase = 0,
			stationarity = ev$cert[c("pass", "max_score", "threshold", "by_block", "scale_change")],
			termination = fit$termination,
			convergence_tolerance = tol,
			message = if (fit$converged) "Converged" else paste0(
				"Stopped (",
				fit$termination,
				") before the convergence criteria were reached."
			)
		),
		eb = list(
			sab = sab,
			Sab0 = ctx$prior$Sab0,
			eta0 = ctx$prior$eta0,
			rho_ab = if (no_ab) NULL else rho_ab,
			sigma_ab = if (no_ab) NULL else sigma_ab,
			sigma2_ab = if (no_ab) NULL else hyp$s2,
			rho_ab_prior = c(mean = ctx$prior$rho_mean, sd = ctx$prior$rho_sd),
			sigma_ab_prior = c(shape = ctx$prior$shape, scale = ctx$prior$scale),
			levels = list(a = lv$a, b = lv$b),
			# the fitted prior of the static factors' rows
			factor_prior = if (ctx$fmode == "fac") { fp <- hyp$fac; attr(fp, "acc") <- NULL; fp },
			# the rule for the factor prior scale and the sd of its log under the cavity rule
			factor_scale = if (ctx$fmode == "fac") list(
				rule = ctx$fscale,
				s = mean(.lde_fac_scales(hyp$fac)),
				sd_log = if (is.null(fac_scale)) NA_real_ else fac_scale$sd_log
			),
			rho_uv = if (dyn_uv) uvh$rho else NULL,
			sigma_uv = if (dyn_uv) sqrt(uvh$s2) else NULL,
			uv_prior = if (dyn_uv) c(
				rho_mean = ctx$uvp$rho_mean,
				rho_sd = ctx$uvp$rho_sd,
				shape = ctx$uvp$shape,
				scale = ctx$uvp$scale
			),
			uv_variance = if (dyn_uv) mean(st$Vpos[ctx$obs]),
			uv_scale = if (dyn_uv) uvh$S,
			uv_kind = if (dyn_uv) ctx$uv_kind,
			t_nu = if (dyn_uv && identical(ctx$uv_kind, "t")) ctx$t_nu,
			restart = .lde_restart_record(ctx, st, if (dyn_uv) fit$profile$grid, isTRUE(fit$converged)),
			uv_search = if (dyn_uv && !is.null(fit$profile)) fit$profile$grid,
			effects = c(sender = !ctx$fzm[1L], receiver = !ctx$fzm[2L]),
			beta_pool = if (lay$any) list(
				mode = ctx$beta_pool,
				sigma_share = hb$share,
				rho_mu = hb$mu,
				rho_tau = if (is.null(hb$tau2)) NULL else sqrt(hb$tau2)
			),
			gamma = NULL,
			vcov = vc,
			vcov_cond = vcc,
			centre = eb_centre,
			predictive = if (!pln && !nrm) .lde_predictive_record(
				ctx,
				st,
				ev,
				Vc,
				nr,
				nc,
				Tt,
				p,
				st$hyp,
				rho_ab,
				sigma_ab
			),
			se = sqrt(diag(vc)),
			se_path = se_path,
			vcov_path = if (lay$any) Vc else NULL,
			beta_rho = if (lay$any) stats::setNames(hb$rho, lay$groups) else NULL,
			beta_sigma = if (lay$any) stats::setNames(sqrt(hb$s2), lay$groups) else NULL,
			se_method = if (dyn_uv) "laplace (penalised Fisher information of the paths, levels and coefficients, at the mean-field mean of the positions with the position noise in the weights)" else if (R > 0L) "laplace (penalised Fisher information of the paths, levels, coefficients and positions)" else "laplace (penalised Fisher information of the paths, levels and coefficients)",
			coef_prior = ctx$coef_prior,
			hyper_gap = ev$cert$scale_change
		),
		meta = list(
			sampler = "dynamic_als",
			algorithm = if (nrm) "dynamic empirical-Bayes smoother of the Gaussian model (AR(1) deviation paths around actor levels, exact posterior given the hyperparameters, EM hyperparameters including the residual variance and the dyadic correlation)" else if (pln) "dynamic empirical-Bayes IRLS of the Poisson-lognormal model (AR(1) deviation paths around actor levels, EM hyperparameters including the variance of the noise, Laplace-marginal coefficients)" else "dynamic empirical-Bayes IRLS (AR(1) deviation paths around actor levels, EM hyperparameters, Laplace-marginal coefficients)",
			uncertainty_available = FALSE,
			approximation_note = paste0(
				"Dynamic empirical-Bayes point estimate: penalised iteratively reweighted least squares (",
				if (nrm) "Gaussian cells with the residual variance as a hyperparameter" else if (pln) "Poisson-lognormal cells with the noise variance as a hyperparameter" else paste0("binary, ", link, " link"),
				") with Gaussian priors on the actor levels and zero-mean AR(1) deviation paths, ",
				"the level covariance, innovation scale and persistence estimated by EM from posterior second moments, ",
				"and Laplace-marginal coefficients. For posterior summaries use lame() with method = \"mcmc\"."
			)
		)
	)
	# the noise variance of a count or normal fit and its prior
	if (pln || nrm) out$eb[c("ve", "ve_prior")] <- list(ve_hat, ctx$prior$ve_prior)
	# dyadic correlation of a directed probit fit (lame_dynamic_als_eb_pair.R), with
	# its score, curvature and sd; .ae_update_response_rho() reads it into VC["rho"]
	if (isTRUE(ctx$dcor)) {
		rs <- if (is.null(ev$rho_step)) 0 else ev$rho_step$rs
		rc <- if (is.null(ev$rho_step)) 0 else ev$rho_step$rc
		rho_d <- unname(hyp$rho_d)
		out$eb$pair <- list(
			rho = rho_d,
			score = rs,
			curvature = -rc,
			sd = if (ctx$pc$n_pairs > 0L) .ae_pair_rho_sd(list(rs = rs, rc = rc), rho_d) else NA_real_,
			n_pairs = ctx$pc$n_pairs,
			quad = if (is.null(ev$wk$pair$K)) NA_integer_ else max(ev$wk$pair$K)
		)
	}
	# dyadic correlation of a directed count panel at the fitted state
	# (lame_dynamic_als_eb_plnpair.R), read into VC["rho"] when overdispersed
	if (!is.null(pair)) {
		out$eb$pair <- list(
			rho = pair$rho,
			score = 0,
			curvature = pair$prec,
			sd = pair$sd,
			n_pairs = pair$n_pairs,
			quad = NA_integer_,
			method = pair$method
		)
	}
	class(out) <- c("lame_dynamic_als", "lame_als", "ame_als")
	out <- .ae_update_response_rho(out, link)
	# a normal fit's dyadic correlation comes from the exact bivariate-normal pair likelihood
	if (nrm && isTRUE(ctx$dcor) && ctx$pc$n_pairs > 0L) {
		rho_d <- unname(hyp$rho_d)
		out$rho_estimation <- list(
			rho = rho_d,
			identified = TRUE,
			n_pairs = ctx$pc$n_pairs,
			boundary = abs(rho_d) >= 0.98,
			sd = out$eb$pair$sd,
			method = "exact bivariate-normal pair likelihood (EM)"
		)
	}
	if (isTRUE(ctx$dcor) && !nrm && is.list(out$rho_estimation)) out$rho_estimation$method <- "Laplace-marginal pair likelihood"
	out$bias_correct <- list(
		method = "eb",
		applied = FALSE,
		requested = bias_requested,
		default = is.null(bias_requested),
		dynamic = TRUE
	)
	# information about the latent product, withheld below the threshold
	out <- .lde_withhold_product(out)
	out
}

# withhold the latent product of a dynamic fit whose signal per loading (the
# weighted squared double-centred product over the noise variance and the
# loadings, as in .ae_gauss_information) is below its threshold, as the static
# estimator does (.ae_product_empty, fit$eb$product_withheld)
.lde_withhold_product <- function(out) {
	R <- out$R
	P <- out$UVPM
	nrm <- identical(out$family, "normal")
	ve <- if (nrm) out$VC[["ve"]] else 1
	if (is.null(R) || R < 1L || is.null(P) || !is.finite(ve) || ve <= 0) return(out)
	Pa <- if (length(dim(P)) == 3L) P else array(P, c(dim(P), out$n_time))
	# cell information: 1 / ve for a normal response, the working Fisher weight otherwise
	Wa <- if (nrm) NULL else out$obs_weights
	signal <- sum(vapply(seq_len(dim(Pa)[3L]), function(t) {
		x <- Pa[, , t]
		x <- x - rowMeans(x) - rep(colMeans(x), each = nrow(x)) + mean(x)
		if (is.null(Wa)) sum(x^2) else sum(Wa[, , t] * x^2, na.rm = TRUE)
	}, 0))
	moving <- isTRUE(out$dynamic_uv)
	rows <- if (isTRUE(out$symmetric)) nrow(Pa) else nrow(Pa) + ncol(Pa)
	loadings <- (if (identical(out$mode, "bipartite")) max(
		.ae_default(out$R_row, R),
		.ae_default(out$R_col, R)
	) else R) * rows *
		(if (moving) dim(Pa)[3L] else 1L)
	thr <- if (!nrm) .lde_glm_product_threshold else if (moving) .lde_moving_product_threshold else .ae_product_threshold
	info <- c(signal = signal, loadings = loadings, per_loading = signal / (ve * loadings), threshold = thr)
	if (!nrm) info <- c(info, signal_per_loading = signal / loadings)
	out$eb$information <- info
	out <- .lde_withhold_intercept(out, moving)
	if (isFALSE(
		getOption("lame.eb_withhold_product", TRUE)
	) || !is.finite(info[["per_loading"]]) || info[["per_loading"]] >= thr) {
		return(out)
	}
	em <- .ae_product_empty(out)
	out <- em$fit
	out$eb$product_withheld <- list(
		reason = "the latent product is not identified by the data at this rank",
		per_loading = info[["per_loading"]],
		threshold = thr,
		loadings = loadings,
		rank = R,
		family = out$family,
		moving = moving,
		log_marginal = NULL,
		rank_chosen = NA_integer_,
		means = em$means
	)
	out
}

# signal per loading of a period below which the moving positions of a dynamic
# normal fit are withheld, calibrated against the products of lame()
.lde_moving_product_threshold <- 1.3

# signal per loading below which the latent product of a dynamic binary or count
# fit is withheld (static factors and moving positions), calibrated against lame()
.lde_glm_product_threshold <- 0.7

# signal per loading from which the intercept of a dynamic binary fit with static
# factors is withheld: the intercept compensates the shrinkage of the variational
# posterior and sits below the sampler's
.lde_binary_intercept_threshold <- 1.5

# withhold the intercept of a dynamic binary fit with static factors above the
# threshold; fit$eb$intercept_withheld keeps the values internal computations use
.lde_withhold_intercept <- function(out, moving) {
	info <- out$eb$information
	if (isFALSE(getOption("lame.eb_withhold_intercept", TRUE)) || !identical(out$family, "binary") || isTRUE(moving) ||
			is.null(info) || !"signal_per_loading" %in% names(info) || !isTRUE(
				is.finite(info[["signal_per_loading"]])
			) ||
			info[["signal_per_loading"]] < .lde_binary_intercept_threshold || !"intercept" %in% names(out$coefficients)) {
		return(out)
	}
	top <- c("mu", "coef_path", "beta_dynamic_per_t", "BETA", "coefficients")
	ebs <- c("vcov", "vcov_cond", "se", "se_path", "vcov_path")
	means <- list(top = out[intersect(top, names(out))], eb = out$eb[intersect(ebs, names(out$eb))])
	for (k in names(means$top)) out[[k]] <- .lde_intercept_na(out[[k]])
	for (k in names(means$eb)) out$eb[[k]] <- .lde_intercept_na(out$eb[[k]])
	out$eb$intercept_withheld <- list(
		reason = "the intercept is not reproduced by the variational posterior of the static factors at this signal",
		per_loading = info[["signal_per_loading"]],
		threshold = .lde_binary_intercept_threshold,
		family = out$family,
		means = means
	)
	out
}

# a vector, matrix or array with the intercept's elements, rows and columns set to NA
.lde_intercept_na <- function(x) {
	if (is.null(x)) return(x)
	if (is.matrix(x)) {
		r <- which(rownames(x) == "intercept"); cc <- which(colnames(x) == "intercept")
		if (length(r)) x[r, ] <- NA_real_
		if (length(cc)) x[, cc] <- NA_real_
	} else if (!is.null(names(x))) {
		x[names(x) == "intercept"] <- NA_real_
	}
	x
}

# restore the withheld intercept values for internal computations
.lde_intercept_restore <- function(object) {
	iw <- object$eb$intercept_withheld
	if (is.null(iw) || is.null(iw$means)) return(object)
	for (k in names(iw$means$top)) object[[k]] <- iw$means$top[[k]]
	for (k in names(iw$means$eb)) object$eb[[k]] <- iw$means$eb[[k]]
	object
}

# the sentence print() and summary() give for a withheld intercept
.lde_intercept_withheld_text <- function(iw) {
	paste0(
		"The intercept is withheld (NA, with its standard error and covariances; the coefficient path of a moving intercept too): at ",
		signif(iw$per_loading, 2),
		" units of signal per loading (the threshold is ",
		iw$threshold,
		") the variational posterior of the static factors ",
		"shrinks the linear predictor and the intercept compensates, so that it sits below the sampler's. Fitted values, predictions and residuals ",
		"carry it; the slopes are reported as usual ({.code fit$eb$intercept_withheld})."
	)
}

# the entry point

# warm start from a previous fit's trajectories: the coefficients and the active
# actors' paths, when they fit this problem
.lde_start_state <- function(ctx, start) {
	st <- .lde_init(ctx)
	if (!is.list(start)) return(st)
	cf <- start$coefficients
	if (!ctx$beta$any && is.numeric(cf) && !is.null(names(cf)) && all(ctx$names %in% names(cf)) &&
			all(is.finite(cf[ctx$names]))) {
		th <- unname(cf[ctx$names])
		th[1L] <- th[1L] + sum(ctx$centre[-1L] * th[-1L])
		st$theta <- th
	}
	# undo the reporting centring to return to the iteration's parametrisation
	ce <- start$eb$centre
	if (!(is.numeric(ce) && length(ce) %in% 2:3 && all(is.finite(ce)))) ce <- c(a = 0, b = 0, o = 0)
	if (length(ce) == 2L) ce <- c(ce, o = 0)
	ap <- if (is.null(start$a_path)) start$a_dynamic else start$a_path
	bp <- if (is.null(start$b_path)) start$b_dynamic else start$b_path
	if (is.matrix(ap) && is.matrix(bp) && ncol(ap) == ctx$Tt && ncol(bp) == ctx$Tt &&
			nrow(ap) == ctx$nr && nrow(bp) == ctx$nc && all(is.finite(ap)) && all(is.finite(bp))) {
		st$x <- rbind(
			ap[ctx$ia, , drop = FALSE] + ce[["a"]],
			if (ctx$sym) NULL else bp[ctx$ib, , drop = FALSE] + ce[["b"]]
		)
		if (!ctx$beta$any) st$theta[1L] <- st$theta[1L] - sum(ce)
		st$lev <- rowMeans(st$x)
		# the levels of an empirical-Bayes fit, when present
		la <- start$a_level; lb <- start$b_level
		if (is.numeric(la) && is.numeric(lb) && length(la) == ctx$nr && length(lb) == ctx$nc &&
				all(is.finite(la)) && all(is.finite(lb))) {
			st$lev <- c(la[ctx$ia] + ce[["a"]], if (ctx$sym) NULL else lb[ctx$ib] + ce[["b"]])
		}
		# variance components of an empirical-Bayes fit of the same kind
		eb <- start$eb
		if (!isTRUE(ctx$static_ab) && is.list(eb) && is.matrix(eb$sab) && all(dim(eb$sab) == 2L) &&
				all(is.finite(eb$sab)) && is.finite(eb$rho_ab) && abs(eb$rho_ab) < 0.99 &&
				is.finite(eb$sigma_ab) && eb$sigma_ab > 0 && all(
					eigen(eb$sab, only.values = TRUE)$values > 0
				)) {
			st$hyp$sab <- if (identical(ctx$kind, "directed")) eb$sab else diag(diag(eb$sab), 2L)
			st$hyp$rho <- eb$rho_ab
			st$hyp$s2 <- eb$sigma_ab^2
		}
		# noise variance of an empirical-Bayes count fit
		if (isTRUE(ctx$vehyp) && is.list(eb) && is.numeric(eb$ve) && length(eb$ve) == 1L && is.finite(eb$ve) && eb$ve > 0) st$hyp$ve <- eb$ve
	}
	.lde_start_restart(ctx, st, start)
}

# signature of the problem a restart record belongs to
.lde_signature <- function(ctx) {
	lay <- ctx$beta
	list(
		kind = ctx$kind,
		na = ctx$na,
		nb = ctx$nb,
		Tt = ctx$Tt,
		R = ctx$R,
		fmode = ctx$fmode,
		uv_kind = ctx$uv_kind,
		t_nu = if (identical(ctx$uv_kind, "t")) ctx$t_nu,
		fscale = if (identical(ctx$fmode, "fac")) ctx$fscale,
		dcor = if (isTRUE(ctx$dcor)) TRUE,
		names = ctx$names,
		mth = ctx$mth,
		dyn_idx = lay$dyn_idx,
		beta_kind = lay$kind,
		beta_pool = ctx$beta_pool,
		link = ctx$link,
		coef_prior = ctx$coef_prior,
		static_ab = isTRUE(ctx$static_ab),
		fz = ctx$fz,
		ia = ctx$ia,
		ib = ctx$ib
	)
}

# the part of the final state the reported quantities do not determine exactly
# (coefficient states, hyperparameters, the latent block in the iteration's gauge);
# a converged refit of the same problem starts at its fixed point
.lde_restart_record <- function(ctx, st, search = NULL, converged = TRUE) {
	rec <- list(signature = .lde_signature(ctx), theta = st$theta, hyp = st$hyp, converged = converged)
	# count fits keep the mode and the cell variances the next noise update reads
	if (isTRUE(ctx$vehyp)) {
		rec[c("x", "lev", "vcell", "cvpair", "dshift")] <- list(st$x, st$lev, st$vcell, st$cvpair, st$dshift)
		# and the tail rate of the noise variance iteration (.lde_tail())
		rec["ve_rate"] <- list(st$ve_rate)
	}
	if (ctx$fmode %in% c("dynamic", "fac")) {
		slim <- function(qs) list(mean = qs$mean, sig = qs$sig, cl = qs$cl, lam = qs$lam)
		rec$pos <- list(qu = slim(st$pos$qu), qv = slim(st$pos$qv))
		rec$sgn <- ctx$sgn
		# the position-scale search that led here (a diagnostic)
		rec$search <- search
	}
	rec
}

# seed the state from a restart record of the same problem; the record wins over
# the reported fields it covers (drop `eb$restart` to edit them)
.lde_start_restart <- function(ctx, st, start) {
	rs <- if (is.list(start$eb)) start$eb$restart
	if (!is.list(rs) || !identical(rs$signature, .lde_signature(ctx))) return(st)
	fin <- function(x) is.numeric(x) && all(is.finite(x))
	if (fin(rs$theta) && length(rs$theta) == ctx$mth) st$theta <- rs$theta
	hy <- rs$hyp
	if (is.list(hy) && is.matrix(hy$sab) && all(dim(hy$sab) == 2L) && fin(hy$sab) && fin(hy$s2) && fin(hy$rho) &&
			(!isTRUE(ctx$vehyp) || (fin(hy$ve) && hy$ve > 0))) {
		st$hyp <- hy
		if (isTRUE(ctx$vehyp) && fin(rs$x) && identical(dim(rs$x), dim(st$x)) && fin(rs$lev) && length(rs$lev) == length(st$lev)) {
			st$x <- rs$x
			st$lev <- rs$lev
			if (fin(rs$vcell) && identical(dim(rs$vcell), dim(ctx$y))) st$vcell <- rs$vcell
			if (fin(rs$cvpair)) st$cvpair <- rs$cvpair
			if (fin(rs$dshift) && identical(dim(rs$dshift), dim(ctx$y))) st$dshift <- rs$dshift
			if (fin(rs$ve_rate) && length(rs$ve_rate) == 1L) st$ve_rate <- rs$ve_rate
			# only a converged record counts as stationary
			if (!ctx$fmode %in% c("dynamic", "fac")) st$seeded <- isTRUE(rs$converged)
		}
	}
	if (ctx$fmode %in% c("dynamic", "fac") && is.list(rs$pos) && is.list(
		if (ctx$fmode == "dynamic") st$hyp$uv else rs$hyp$fac
	)) {
		Tp <- if (ctx$fmode == "fac") 1L else ctx$Tt
		ok <- function(qs, n) is.list(qs) && identical(dim(qs$mean), c(n, ctx$R, Tp)) &&
			identical(dim(qs$sig), c(ctx$R, ctx$R, n, Tp)) && fin(qs$mean) && fin(qs$sig) && fin(qs$cl) &&
			(if (identical(ctx$uv_kind, "t")) is.matrix(qs$lam) && identical(dim(qs$lam), c(n, ctx$Tt)) && fin(qs$lam) && all(qs$lam > 0) else is.null(qs$lam))
		if (ok(rs$pos$qu, ctx$na) && ok(rs$pos$qv, ctx$nb)) {
			st$pos <- rs$pos
			if (ctx$fmode == "fac") {
				st$hyp$fac <- rs$hyp$fac
				if (ctx$sym) st$sgn <- rs$sgn
				if (ctx$sym) ctx$sgn <- rs$sgn
				pp <- .lde_fac_product(ctx, .lde_fac_state(ctx, st$pos, st$hyp$fac))
				st$O <- pp$O
				st$Vpos <- pp$Vpos
			} else {
				st$O <- .lde_pos_product(st$pos, ctx$Tt)
				st$Vpos <- .lde_pos_prodvar(st$pos, ctx$R, ctx$Tt)
			}
			# only a converged record counts as stationary
			st$seeded <- isTRUE(rs$converged)
			st$search <- rs$search
			if (ctx$sym) st$sgn <- rs$sgn
		}
	}
	st
}

# refuse an unfinished fit (no convergence, or a position-scale search at the edge
# or unsettled); the condition carries the last fit object in `fit`
.lde_refuse_unfinished <- function(out, fit, ctx, max_iter, edge, unsettled, call = parent.frame(2L)) {
	st <- out$diagnostics$stationarity
	score <- st$max_score
	thr <- st$threshold
	blk <- if (length(st$by_block)) names(st$by_block)[which.max(st$by_block)] else NA_character_
	tr <- out$diagnostics$convergence_trace
	nt <- if (is.data.frame(tr)) nrow(tr) else 0L
	fall <- if (nt > 10L) tr$merit[nt - 10L] / tr$merit[nt] else NA_real_
	pr <- fit$profile
	if (edge) {
		rg <- pr$range
		lower <- is.finite(pr$S) && length(rg) == 2L && pr$S < sqrt(rg[1L] * rg[2L])
		cli::cli_abort(
			c(
				"The panel does not identify the stationary variance of the moving positions.",
				"x" = if (lower) {
					"The marginal objective is lowest at the smallest scale searched, {signif(rg[1L], 3)} (the range searched: {signif(rg[1L], 3)} to {signif(rg[2L], 3)}): at rank {ctx$R} the panel carries no variation that the positions can resolve."
				} else {
					"The marginal objective keeps falling up to the largest scale searched, {signif(rg[2L], 3)} (the range searched: {signif(rg[1L], 3)} to {signif(rg[2L], 3)}): at rank {ctx$R} the positions reproduce the panel exactly."
				},
				"i" = "Lower {.arg R}, use {.code bias_correct = \"none\"}, or fit with {.code method = \"mcmc\"}. The fit at the end of the range is in the {.field fit} field of the condition."
			),
			class = c("lame_dynamic_uv_scale_edge", "lame_als_eb"),
			call = call,
			fit = out
		)
	}
	if (unsettled) {
		# the unsettled count over the five evaluations nearest the best scale
		n_near <- pr$n_near
		n_bad <- pr$n_unsettled
		if (is.null(n_near) || is.null(n_bad)) {
			g <- pr$grid
			near <- if (is.data.frame(g) && nrow(g) && is.finite(pr$S)) order(abs(log(g$S) - log(pr$S)))[seq_len(min(5L, nrow(g)))] else integer(0)
			n_near <- length(near)
			n_bad <- sum(!g$converged[near])
		}
		# an unsettled search failed at one evaluation at least
		n_bad <- max(n_bad, 1L)
		cli::cli_abort(
			c(
				"The iteration did not settle at {n_bad} of the {n_near} position scale{?s} nearest the minimum of the search within {max_iter} passes.",
				"x" = "The profile of the marginal objective over the scale of the positions is not determined: the positions reproduce the panel at these scales, or the panel is too sparse to identify them.",
				"i" = "Lower {.arg R}, use {.code bias_correct = \"none\"}, or fit with {.code method = \"mcmc\"}. The last iterate is in the {.field fit} field of the condition."
			),
			class = c("lame_dynamic_uv_scale_unsettled", "lame_als_convergence", "lame_als_eb"),
			call = call,
			fit = out
		)
	}
	cli::cli_abort(
		c(
			"{.fn lame_dynamic_als} did not reach the empirical-Bayes convergence criterion: it stopped after {fit$iterations} pass{?es} ({fit$termination}).",
			"x" = "The largest scaled residual of the empirical-Bayes equations is {signif(score, 3)} against the threshold {signif(thr, 3)}, in the block {.field {blk}}{if (is.finite(fall)) paste0('; the merit fell by a factor of ', signif(fall, 3), ' over the last ten passes') else ''}.",
			"i" = "The last iterate is in the {.field fit} field of the condition. {.arg max_iter} is {max_iter}; {.code fit$diagnostics$convergence_trace} holds the passes."
		),
		class = c("lame_als_convergence", "lame_als_eb"),
		call = call,
		fit = out
	)
}

# the dynamic empirical-Bayes fit of lame_dynamic_als(): context, reweighting
# loop and fitted object
.lde_dynamic_fit <- function(
	prep,
	family,
	mode,
	symmetric,
	link,
	prior,
	max_iter,
	tol,
	verbose,
	seed,
	call,
	call_args,
	bias_requested,
	start = NULL,
	start_jitter = 0,
	stability = "none",
	coef_prior = "mcmc",
	R = 0L,
	beta_mask = NULL,
	beta_kind = "ar1",
	dynamic_ab = TRUE,
	rvar = TRUE,
	cvar = TRUE,
	R_row = NULL,
	R_col = NULL,
	beta_pool = "none",
	dynamic_uv = FALSE,
	dynamic_uv_kind = "ar1",
	dcor = FALSE
) {
	# the pair quadrature cache belongs to one fit
	.lde_pair_reset()
	on.exit(.lde_pair_reset(), add = TRUE)
	ctx <- .lde_context(
		prep,
		symmetric,
		link,
		prior,
		coef_prior = coef_prior,
		R = R,
		beta_mask = beta_mask,
		beta_kind = beta_kind,
		dynamic_ab = dynamic_ab,
		rvar = rvar,
		cvar = cvar,
		beta_pool = beta_pool,
		dynamic_uv = dynamic_uv,
		dynamic_uv_kind = dynamic_uv_kind,
		dcor = dcor
	)
	if (verbose) {
		cli::cli_h3("Dynamic empirical-Bayes ALS fit")
		cli::cli_text("Network: {prep$nr} x {prep$nc}, {prep$Tt} time slices; {ctx$q} actor effects with AR(1) deviation paths")
		if (isTRUE(ctx$dcor)) cli::cli_text("Dyadic correlation: the {ctx$pc$n_pairs} dyad-periods observed in both directions are bivariate {if (isTRUE(ctx$nrm)) 'normal' else 'probit'} pairs")
	}
	st0 <- .lde_start_state(ctx, start)
	if (start_jitter > 0) {
		set.seed(seed)
		st0$theta <- st0$theta + stats::rnorm(length(st0$theta), 0, start_jitter / 5)
		st0$x <- st0$x + matrix(stats::rnorm(length(st0$x), 0, start_jitter), nrow(st0$x))
	}
	fit <- .lde_fit(
		ctx,
		st0,
		max_iter = max_iter,
		tol = tol,
		verbose = FALSE,
		jitter = start_jitter,
		progress = verbose
	)
	ctx <- fit$ctx
	# a fit is finished when it converged and any position-scale search found an
	# interior minimiser; anything else is refused below
	edge <- isTRUE(fit$profile$at_edge)
	unsettled <- isTRUE(fit$profile$unsettled)
	finished <- isTRUE(fit$converged) && !edge && !unsettled
	if (verbose) {
		if (finished) {
			cli::cli_alert_success("Dynamic empirical-Bayes ALS converged in {fit$iterations} iteration{?s}.")
		} else {
			cli::cli_alert_danger("Dynamic empirical-Bayes ALS stopped ({fit$termination}) without reaching its convergence criterion.")
		}
	}
	# directed count panel: the dyadic correlation at the fitted state
	if (finished && isTRUE(ctx$pair_stage)) {
		fit$pair <- tryCatch(.lde_pln_pair_stage(ctx, fit), error = function(e) {
			cli::cli_warn(
				c(
					"The dyadic correlation of the counts could not be computed; the cells of a dyad are treated as independent.",
					"x" = conditionMessage(e)
				),
				class = "lame_als_pair"
			)
			NULL
		})
	}
	out <- .lde_as_fit(
		prep,
		ctx,
		fit,
		call,
		call_args,
		mode,
		symmetric,
		link,
		bias_requested,
		tol,
		seed,
		R_row = R_row,
		R_col = R_col
	)
	out$als_refit <- list(
		fitter = "lame_dynamic_als",
		args = call_args[setdiff(names(call_args), c("Y", "Xdyad", "call"))],
		seed = seed
	)
	if (!finished) .lde_refuse_unfinished(out, fit, ctx, max_iter, edge, unsettled)
	if (!identical(stability, "none")) {
		out$stability <- .lda_dynamic_stability(
			out,
			starts = stability,
			seed = seed,
			call_args = call_args
		)
	}
	out
}

# printing

# the lines a dynamic empirical-Bayes fit adds to print() and summary()
.lde_print_eb <- function(x, digits = 4) {
	if (!is.null(x$eb$product_withheld)) cli::cli_text(.ae_product_withheld_text(x$eb$product_withheld))
	if (!is.null(x$eb$intercept_withheld)) cli::cli_text(.lde_intercept_withheld_text(x$eb$intercept_withheld))
	eb <- if (is.list(x$eb)) x$eb else x$diagnostics$eb
	if (is.null(x$bias_correct) || !isTRUE(x$bias_correct$dynamic) || is.null(eb)) {
		return(invisible(NULL))
	}
	if (!is.null(eb$rho_ab)) {
		cli::cli_text("Dynamic effects: persistence rho_ab = {round(eb$rho_ab, digits)}, innovation scale sigma_ab = {round(eb$sigma_ab, digits)} (EM estimates under the {.fn lame} hyperpriors); level covariance in {.code fit$eb$sab}.")
	} else if (!any(eb$effects)) {
		cli::cli_text("Additive effects: none (every margin disabled).")
	} else {
		cli::cli_text("Additive effects: static levels (EM estimate of their covariance under the {.fn lame} hyperprior in {.code fit$eb$sab}).")
	}
	if (!is.null(eb$pair)) {
		pr <- eb$pair
		if (identical(x$family, "normal")) {
			cli::cli_text("Dyadic correlation: rho = {round(pr$rho, digits)}{if (is.finite(pr$sd)) paste0(' (posterior sd ', round(pr$sd, digits), ')') else ''}, one for all periods: the residuals of the two cells of each of the {pr$n_pairs} dyad-periods observed in both directions are one bivariate normal observation, under the {.fn lame} prior of the correlation.")
		} else if (identical(x$family, "poisson")) {
			cli::cli_text("Dyadic correlation: rho = {round(pr$rho, digits)}{if (is.finite(pr$sd)) paste0(' (posterior sd ', round(pr$sd, digits), ')') else ''}, one for all periods: the correlation of the noises of the two cells of each of the {pr$n_pairs} dyad-periods observed in both directions, a bivariate Poisson-lognormal observation, estimated at the fitted state of independent cells under the {.fn lame} prior of the correlation.")
		} else {
			cli::cli_text("Dyadic correlation: rho = {round(pr$rho, digits)}{if (is.finite(pr$sd)) paste0(' (posterior sd ', round(pr$sd, digits), ')') else ''}, one for all periods: the two cells of each of the {pr$n_pairs} dyad-periods observed in both directions are one bivariate probit observation, under the {.fn lame} prior of the correlation.")
		}
	}
	if (!is.null(eb$rho_uv)) {
		cli::cli_text("Moving latent positions (mean-field Gaussian posterior{if (identical(eb$uv_kind, 't')) paste0(' with Student-t transitions, nu = ', eb$t_nu)}): persistence rho_uv = {round(eb$rho_uv, digits)}, innovation scale sigma_uv = {round(eb$sigma_uv, digits)} (stationary variance {round(eb$uv_scale, digits)} chosen by the profile of the marginal objective, persistence by EM, under the {.fn lame} hyperpriors).")
	}
	if (!is.null(eb$factor_scale)) {
		fs <- eb$factor_scale
		cli::cli_text("Static latent factors (mean-field Gaussian posterior, the static estimator's): the scale of their prior is set by the {fs$rule} rule{if (identical(fs$rule, 'cavity')) ' (second moments at the cavity covariances of the rows)' else ' (second moments at the mean-field covariances)'}.")
	}
	if (!is.null(eb$se)) {
		cli::cli_text("Laplace standard errors (coefficients, marginal over the paths, levels{if (!is.null(eb$factor_scale)) ' and factors' else ''}{if (!is.null(eb$factor_scale) && is.finite(eb$factor_scale$sd_log)) ', with the uncertainty of the factor scale' else ''}, conditional on the variance components):")
		cli::cli_verbatim(utils::capture.output(print(round(eb$se, digits))))
	}
	invisible(NULL)
}
