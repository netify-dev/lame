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

# cell terms on the observed entries of (y, eta): log likelihood, score, information,
# its two derivatives, and the noise's second moment and mean; `lgy` is lgamma(y + 1).
# compiled kernel (lde_pln_cells_cpp)
.lde_pln_terms <- function(y, eta, ve, lgy = NULL, K = .lde_pln_kin(ve)) {
	d <- .ae_dims(y)
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

# expected cells by the compiled kernel (lde_pln_expect_cpp)
.lde_pln_expect <- function(y, m, v, ve, lgy = NULL) {
	d <- .ae_dims(y)
	v <- array(pmax(rep_len(as.numeric(v), length(y)), 0), d)
	Kin <- .lde_pln_kin(ve)
	vmax <- suppressWarnings(max(v[is.finite(y) & is.finite(m)], 0))
	K2 <- .ae_pln_kin(1.5 * (ve + vmax))
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
