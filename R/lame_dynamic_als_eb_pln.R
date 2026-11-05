# dynamic empirical-Bayes estimator for count panels (lame_dynamic_als() with
# family = "Poisson"): the estimator of lame_dynamic_als_eb.R with the
# Poisson-lognormal cell of ame_als_pln.R in place of the probit cell.
# the cell integrates out the noise e ~ N(0, ve); its derivatives are cumulants of
# the noise posterior given the count (src/ae_pln.h). ve is an extra hyperparameter
# updated by its EM image (Gauss-Seidel, before each working problem), and its slow
# contraction is certified and extrapolated from the tail of its moves.
# the linear predictor is on the latent scale: fitted() is exp(EZ + ve / 2)

# cell likelihood

# nodes of the inner rule, chosen once per stage from 1.5 ve so the objective
# does not change as ve settles
.lde_pln_kin <- function(ve) .ae_pln_kin(1.5 * ve)

# R reference for pln_node() of src/ae_pln.h with the fourth cumulant, in chunks to bound memory
.lde_pln_node_block <- function(y, eta, ve, K) {
	n <- length(y)
	eta <- pmin(pmax(eta, -700), 700)
	ve <- pmin(pmax(ve, 1e-8), 1e6)
	iv <- 1 / ve
	ly <- ifelse(y > 0, log(y), -Inf)
	zeta <- pmax(eta, ly)
	u <- y + eta / ve
	act <- seq_len(n)
	for (it in seq_len(100L)) {
		if (!length(act)) break
		zt <- zeta[act]
		ez <- exp(zt)
		step <- (ez + zt / ve[act] - u[act]) / (ez + iv[act])
		zt <- zt - step
		zeta[act] <- zt
		act <- act[!(abs(step) <= 1e-13 * pmax(1, abs(zt)))]
	}
	Xs <- exp(zeta)
	s <- 1 / sqrt(Xs + iv)
	estar <- zeta - eta
	rule <- .ae_pln_rule(K)
	Z <- matrix(rule$z, n, K, byrow = TRUE)
	tt <- s * Z
	r <- expm1(tt) - tt - 0.5 * tt * tt
	sm <- abs(tt) < 0.1
	if (any(sm)) {
		t1 <- tt[sm]
		r[sm] <- t1 * t1 * t1 / 6 * (1 + t1 / 4 * (1 + t1 / 5 * (1 + t1 / 6 * (1 + t1 / 7 * (1 + t1 / 8 *
				(1 + t1 / 9 * (1 + t1 / 10 * (1 + t1 / 11))))))))
	}
	lw <- matrix(rule$lw, n, K, byrow = TRUE) - Xs * r
	ex <- expm1(tt)
	mx <- lw[cbind(seq_len(n), max.col(lw, ties.method = "first"))]
	p <- exp(lw - mx)
	tot <- rowSums(p)
	p <- p / tot
	ex[p == 0] <- 0
	ez1 <- rowSums(p * Z)
	d <- Z - ez1
	vz <- rowSums(p * d * d)
	k3 <- rowSums(p * d * d * d)
	k4 <- rowSums(p * d * d * d * d) - 3 * vz * vz
	eex <- rowSums(p * ex)
	cze <- rowSums(p * d * (ex - eex))
	me <- estar + s * ez1
	gs <- y * estar - Xs - estar * estar / (2 * ve)
	list(
		l = y * eta - lgamma(y + 1) - 0.5 * log(2 * pi * ve) + gs + log(s) + 0.918938533204672741780329736406 +
			mx + log(tot),
		s = me / ve,
		w = Xs * s * cze / ve,
		hp = -(s * s * s) * k3 / (ve * ve * ve),
		hpp = -(s * s * s * s) * k4 / (ve * ve * ve * ve),
		ee = s * s * vz + me * me,
		me = me
	)
}

.lde_pln_node_r <- function(y, eta, ve, K = .lde_pln_kin(max(ve)), chunk = 20000L) {
	n <- length(y)
	ve <- rep_len(ve, n)
	out <- vector("list", ceiling(n / chunk))
	for (k in seq_along(out)) {
		ix <- ((k - 1L) * chunk + 1L):min(k * chunk, n)
		out[[k]] <- .lde_pln_node_block(y[ix], eta[ix], ve[ix], K)
	}
	if (!length(out)) return(list(
		l = numeric(),
		s = numeric(),
		w = numeric(),
		hp = numeric(),
		hpp = numeric(),
		ee = numeric(),
		me = numeric()
	))
	if (length(out) == 1L) return(out[[1L]])
	nm <- names(out[[1L]])
	setNames(lapply(nm, function(a) unlist(lapply(out, `[[`, a), use.names = FALSE)), nm)
}

# cell terms on the observed entries of (y, eta): log likelihood, score, information,
# its two derivatives, and the noise's second moment and mean; `lgy` is lgamma(y + 1).
# compiled kernel, or the R reference with options(lame.als_engine = "r")
.lde_pln_terms <- function(y, eta, ve, lgy = NULL, K = .lde_pln_kin(ve)) {
	d <- .ae_dims(y)
	if (identical(.ae_engine(), "r")) {
		ix <- which(is.finite(y) & is.finite(eta))
		nd <- .lde_pln_node_r(y[ix], eta[ix], ve, K)
		out <- lapply(nd, function(v) { o <- array(0, d); o[ix] <- v; o })
		return(out)
	}
	if (is.null(lgy)) lgy <- lgamma(pmax(as.numeric(y), 0, na.rm = TRUE) + 1)
	rule <- .ae_pln_rule(K)
	out <- lde_pln_cells_cpp(
		as.numeric(y),
		as.numeric(eta),
		as.numeric(ve),
		as.numeric(lgy),
		rule$z,
		rule$lw,
		.ae_threads()
	)
	for (nm in names(out)) dim(out[[nm]]) <- d
	out
}

# working problem of the count cells at `eta` and `ve` (the fields of .lde_work()),
# with the response increment capped and the weight floored as in the binary problem;
# `ee` is the plug-in statistic of the EM image of ve
.lde_work_pln <- function(ctx, eta, ve, delta_cap = 30, wfloor = 1e-8) {
	y <- ctx$y
	obs <- ctx$obs
	nd <- .lde_pln_terms(y, eta, ve, ctx$lgy, .lde_pln_kin(ve))
	score <- nd$s
	w <- pmax(nd$w, wfloor)
	delta <- score / w
	big <- abs(delta) > delta_cap
	if (any(big)) {
		delta[big] <- sign(delta[big]) * delta_cap
		w[big] <- abs(score[big]) / delta_cap
	}
	z <- eta + delta
	hprime <- nd$hp
	hpp <- nd$hpp
	z[!obs] <- 0
	w[!obs] <- 0
	score[!obs] <- 0
	hprime[!obs | !is.finite(hprime)] <- 0
	hpp[!obs | !is.finite(hpp)] <- 0
	ee <- nd$ee
	ee[!obs] <- 0
	list(z = z, w = w, score = score, hprime = hprime, hpp = hpp, dev = -2 * sum(nd$l[obs]), ee = ee)
}

# deviance of the count cells at `eta`
.lde_dev_pln <- function(ctx, eta, ve) {
	nd <- .lde_pln_terms(ctx$y, eta, ve, ctx$lgy, .lde_pln_kin(ve))
	-2 * sum(nd$l[ctx$obs])
}

# noise variance

# statistic of the EM image of ve: expected second moment of the noise summed over
# observed cells, over N(eta, vcell) and the noise posterior; `vcell` NULL is the plug-in
.lde_pln_ee_sum <- function(ctx, eta, vcell, ve) {
	if (is.null(vcell) || !any(vcell > 0)) {
		nd <- .lde_pln_terms(ctx$y, eta, ve, ctx$lgy, .lde_pln_kin(ve))
		return(sum(nd$ee[ctx$obs]))
	}
	# a short rule over the predictor and half the cell's noise rule are accurate enough for this sum
	vmax <- max(vcell[ctx$obs])
	cl <- .ae_pln_cells(ctx$y, eta, vcell, ve, if (vmax <= 1.5) 3L else 5L, max(12L, ceiling(.lde_pln_kin(ve) / 2)))
	sum(cl$ee[ctx$obs])
}

# EM image of ve: posterior mean of the sampler's inverse-gamma conditional (an
# undirected network counts each dyad once)
.lde_pln_ve_image <- function(ctx, ee_sum) {
	pr <- ctx$prior$ve_prior
	max(
		(pr[["nu0"]] * pr[["s20"]] + ctx$nu * ee_sum) / (ctx$nu * sum(ctx$obs) + pr[["nu0"]] - 2),
		1e-8 * ctx$prior$vunit
	)
}

# hyperprior of ve on the -2 log density scale: the sampler's scaled inverse chi-square
.lde_pln_ve_hyper <- function(ctx, ve) {
	pr <- ctx$prior$ve_prior
	-2 * (-0.5 * (pr[["nu0"]] - 2) * log(ve) - 0.5 * pr[["nu0"]] * pr[["s20"]] / ve)
}

# priors of a count panel

# additive-effects hyperprior for counts: variance of the centred log row and column
# rates, with round(4 + 3 n / 100) degrees of freedom
.lde_hyperprior_pois <- function(Yarr) {
	list(Sab0 = diag(2) * .ame_pois_effect_scale(Yarr), eta0 = round(4 + 3 * dim(Yarr)[1L] / 100))
}

# starting intercept: log mean count less a quarter of the starting noise variance
.lde_pln_start_intercept <- function(ctx) {
	ybar <- mean(ctx$y[ctx$obs])
	if (!is.finite(ybar) || ybar <= 0) ybar <- 1
	log(max(ybar, 1e-3)) - .ae_pln_ve_start / 2
}


# cell expected over a Gaussian predictor

# outer rule by the largest predictor variance: the static rule up to 3, then 21 nodes
# up to 10 and 41 beyond
.lde_pln_outer_rule <- function(vmax) .ae_pln_outer(
	if (!is.finite(vmax) || vmax <= 3) 9L else if (vmax <= 10) 21L else 41L
)

# R reference of lde_pln_expect_cpp(): the cell terms averaged over the predictor
# m + sqrt(v) z by an outer Gauss-Hermite rule; `lm` is the log expected probability of the count
.lde_pln_expect_r <- function(y, m, v, ve, Kin, K2) {
	d <- .ae_dims(y)
	ix <- which(is.finite(y) & is.finite(m))
	mm <- m[ix]
	vv <- pmax(rep_len(v, length(y))[ix], 0)
	yy <- y[ix]
	rules <- .lde_pln_outer_rule(suppressWarnings(max(vv, 0)))
	rule_of <- findInterval(vv, rules$vcut, left.open = TRUE) + 1L
	n <- length(ix)
	acc <- list(es = numeric(n), ew = numeric(n), eh = numeric(n), ehh = numeric(n), el = numeric(n), ee = numeric(n))
	for (r in seq_along(rules$z)) {
		sel <- which(rule_of == r)
		if (!length(sel)) next
		sdv <- sqrt(vv[sel])
		for (k in seq_along(rules$z[[r]])) {
			nd <- .lde_pln_node_r(yy[sel], mm[sel] + sdv * rules$z[[r]][k], ve, Kin)
			wk <- rules$w[[r]][k]
			acc$es[sel] <- acc$es[sel] + wk * nd$s
			acc$ew[sel] <- acc$ew[sel] + wk * nd$w
			acc$eh[sel] <- acc$eh[sel] + wk * nd$hp
			acc$ehh[sel] <- acc$ehh[sel] + wk * nd$hpp
			acc$el[sel] <- acc$el[sel] + wk * nd$l
			acc$ee[sel] <- acc$ee[sel] + wk * nd$ee
		}
	}
	acc$lm <- .lde_pln_node_r(yy, mm, ve + vv, K2)$l
	lapply(acc, function(a) { o <- array(0, d); o[ix] <- a; o })
}

# expected cells by the compiled kernel (or the R reference with options(lame.als_engine = "r"))
.lde_pln_expect <- function(y, m, v, ve, lgy = NULL) {
	d <- .ae_dims(y)
	v <- array(pmax(rep_len(as.numeric(v), length(y)), 0), d)
	Kin <- .lde_pln_kin(ve)
	vmax <- suppressWarnings(max(v[is.finite(y) & is.finite(m)], 0))
	K2 <- .ae_pln_kin(1.5 * (ve + vmax))
	if (identical(.ae_engine(), "r")) return(.lde_pln_expect_r(y, m, v, ve, Kin, K2))
	if (is.null(lgy)) lgy <- lgamma(pmax(as.numeric(y), 0, na.rm = TRUE) + 1)
	outer <- .lde_pln_outer_rule(vmax)
	rule <- .ae_pln_rule(Kin)
	rule2 <- .ae_pln_rule(K2)
	out <- lde_pln_expect_cpp(
		as.numeric(y),
		as.numeric(lgy),
		as.numeric(m),
		as.numeric(v),
		as.numeric(ve),
		outer$z,
		outer$w,
		outer$vcut,
		rule$z,
		rule$lw,
		rule2$z,
		rule2$lw,
		.ae_threads()
	)
	for (nm in names(out)) dim(out[[nm]]) <- d
	out
}

# working problem of the count cells whose predictor is N(eta, vnoise): the fields of
# .lde_work_gh(), the marginal deviance and the expected second moment of the noise
.lde_work_gh_pln <- function(ctx, eta, vnoise, ve, delta_cap = 30, wfloor = 1e-8) {
	obs <- ctx$obs
	cl <- .lde_pln_expect(ctx$y, eta, vnoise, ve, ctx$lgy)
	score <- cl$es
	w <- pmax(cl$ew, wfloor)
	delta <- score / w
	big <- abs(delta) > delta_cap
	if (any(big)) {
		delta[big] <- sign(delta[big]) * delta_cap
		w[big] <- abs(score[big]) / delta_cap
	}
	z <- eta + delta
	hprime <- cl$eh
	hpp <- cl$ehh
	z[!obs] <- 0
	w[!obs] <- 0
	score[!obs] <- 0
	hprime[!obs | !is.finite(hprime)] <- 0
	hpp[!obs | !is.finite(hpp)] <- 0
	ee <- cl$ee
	ee[!obs] <- 0
	list(
		z = z,
		w = w,
		score = score,
		hprime = hprime,
		hpp = hpp,
		dev = -2 * sum(cl$el[obs]),
		dev_marg = -2 * sum(cl$lm[obs]),
		ee = ee
	)
}

# cavity weights of the count cells (.ae_pln_info_w())
.lde_cavity_w_pln <- function(ctx, eta, V3, hyp) {
	# pairs of a directed panel at the fitted correlation (lame_dynamic_als_eb_plnpair.R),
	# or the independent cells of the iteration
	if (!is.null(ctx$pair_info)) return(
		.lde_cavity_w_pln_pair(ctx, eta, V3, ctx$pair_info$cvv, hyp$ve, ctx$pair_info$rho)
	)
	list(w = .ae_pln_info_w(ctx$y, eta, V3, hyp$ve), cw = NULL)
}

# posterior mean of the paths

# first-order shift of the posterior mean of the additive paths and levels from their
# mode, E[x] - mode = -u / 2 with u = H^-1 E'(h' v) (the vector .lde_laplace() forms);
# the sampler reports posterior means, and a skewed posterior's mean is not its mode
.lde_pln_mean_shift <- function(ctx, ev) {
	vq <- lapply(seq_len(ctx$Tt), function(t) .lde_stack_sums(ev$lap$Tcell[, , t], ctx))
	sol <- .lde_solve_cond(ev$fac, vq)
	list(x = do.call(cbind, sol$s), lev = sol$lev)
}

# shift of the cells' predictors from the mode to the first-order posterior mean; the
# EM image of ve is taken about this centre because its fixed point is sensitive to it
.lde_pln_dshift <- function(ctx, st, eta, lap, fac) {
	if (is.null(lap$Tcell)) return(NULL)
	ms <- .lde_pln_mean_shift(ctx, list(lap = lap, fac = fac))
	xs <- st$x - 0.5 * ms$x
	ls <- st$lev - 0.5 * ms$lev
	if (any(ctx$fz)) { xs[ctx$fz, ] <- 0; ls[ctx$fz] <- 0 }
	stw <- st
	stw$x <- xs
	stw$lev <- ls
	.lde_eta(ctx, stw) - eta
}

# certificate with the residual of the noise variance's EM equation among the hyperparameters' moves
.lde_cert_with_gap <- function(cert, ve_gap) {
	if (is.null(ve_gap)) return(cert)
	cert$scale_change <- max(cert$scale_change, ve_gap)
	cert$merit <- max(cert$max_score, cert$scale_change)
	cert$pass <- isTRUE(cert$pass) && is.finite(ve_gap) && ve_gap <= cert$threshold
	cert
}

# tail of a hyperparameter's iteration

# distance of a slow geometric EM tail to its limit, from the last four moves of
# coordinate `h`; returns the distance, rate, last move, sign agreement and whether the tail is geometric
.lde_tail <- function(h, n_inc = 4L, cap = 0.9999) {
	n <- length(h)
	unknown <- list(dist = Inf, rate = NA_real_, move = NA_real_, mono = FALSE, stable = FALSE)
	if (n < n_inc + 1L) return(unknown)
	inc <- diff(h[(n - n_inc):n])
	if (!all(is.finite(inc))) return(unknown)
	a <- abs(inc)
	last <- a[n_inc]
	if (last <= 1e-14) return(list(dist = 0, rate = 0, move = inc[n_inc], mono = TRUE, stable = TRUE))
	ratios <- a[-1L] / pmax(a[-n_inc], 1e-300)
	r <- min(max(ratios), cap)
	list(
		dist = last * r / (1 - r),
		rate = r,
		move = inc[n_inc],
		mono = all(sign(inc) == sign(inc[n_inc])),
		stable = (max(ratios) - min(ratios)) <= 0.1 * (1 - r)
	)
}

# whether a scalar coordinate's tail is at its limit within half the threshold
.lde_tail_done <- function(tl, thr) {
	isTRUE(tl$dist <= thr / 2) && (isTRUE(tl$stable) || isTRUE(tl$dist <= thr / 20))
}

# the same for the columns of `H`, one row per accepted iterate
.lde_tail_ok <- function(H, thr) {
	if (is.null(H) || nrow(H) < 5L) return(FALSE)
	all(vapply(seq_len(ncol(H)), function(j) .lde_tail_done(.lde_tail(H[, j]), thr), TRUE))
}

# scalar hyperparameters whose tails a count fit with latent factors watches: log ve
# and the log scales of the factor prior
.lde_tail_scalars <- function(ctx, st) {
	v <- log(st$hyp$ve)
	if (identical(ctx$fmode, "fac") && !is.null(st$hyp$fac)) v <- c(v, log(.lde_fac_scales(st$hyp$fac)))
	v
}
