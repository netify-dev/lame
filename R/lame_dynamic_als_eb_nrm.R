# dynamic empirical-Bayes estimator for normal panels (lame_dynamic_als(), family = "normal"): the engine of
# lame_dynamic_als_eb.R with a Gaussian cell, whose working problem is exact (weight 1 / ve, no Laplace
# correction); it adds the EM step of ve and, for directed panels with dcor, an exact bivariate-normal pair

# additive-effects hyperprior lame() uses for a normal response: Sab0 from the mean variance of the row and
# column means of the centred response, with round(4 + 3 n / 100) degrees of freedom (n row actors)
.lde_hyperprior_nrm <- function(Yarr) {
	Ec <- Yarr - mean(Yarr, na.rm = TRUE)
	a_tmp <- apply(Ec, 1L, mean, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
	b_tmp <- apply(Ec, 2L, mean, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
	vscale <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
	if (!is.finite(vscale) || vscale <= 0) vscale <- stats::var(c(Yarr), na.rm = TRUE)
	if (!is.finite(vscale) || vscale <= 0) vscale <- 1
	list(Sab0 = diag(2) * vscale, eta0 = round(4 + 3 * dim(Yarr)[1L] / 100))
}

# starting intercept: the mean response
.lde_nrm_start_intercept <- function(ctx) {
	ybar <- mean(ctx$y[ctx$obs])
	if (!is.finite(ybar)) ybar <- 0
	ybar
}

# starting residual variance: half the response variance
.lde_nrm_start_ve <- function(ctx) {
	v <- stats::var(ctx$y[ctx$obs])
	if (!is.finite(v) || v <= 0) v <- 1
	v / 2
}

# 2 x 2 precision of a pair of cells: diagonal weight w and cross weight c
.lde_nrm_pair_info <- function(ve, rho) {
	d <- ve * (1 - rho^2)
	list(w = 1 / d, c = -rho / d)
}

# working problem of the cells at `eta` (the fields of .lde_work()) with residual variance hyp$ve and,
# for directed panels with dcor, hyp$rho_d; with `vnoise` (each predictor's variance under q) it also
# returns the expected and marginal deviances, which carry the normalising constants
.lde_work_nrm <- function(ctx, eta, hyp, vnoise = NULL) {
	ve <- hyp$ve
	obs <- ctx$obs
	y <- ctx$y
	r <- y - eta
	r[!obs] <- 0
	vn <- if (is.null(vnoise)) NULL else { v <- vnoise; v[!obs] <- 0; v }
	w <- array(0, dim(y))
	w[obs] <- 1 / ve
	score <- r / ve
	z <- y
	z[!obs] <- 0
	zero <- array(0, dim(y))
	nobs <- sum(obs)
	# independent cells: -2 log N(y; eta, ve) = (r^2 + v) / ve + log(2 pi ve), marginally N(eta, ve + v)
	dev <- sum((r^2 + (if (is.null(vn)) 0 else vn))[obs]) / ve + nobs * log(2 * pi * ve)
	dev_marg <- if (is.null(vn)) NULL else sum((r^2 / (ve + vn))[obs]) + sum(log(2 * pi * (ve + vn[obs])))
	out <- list(z = z, w = w, score = score, hprime = zero, hpp = zero, dev = dev, dev_marg = dev_marg)
	if (!isTRUE(ctx$dcor) || is.null(ctx$pc) || !ctx$pc$n_pairs) return(out)
	# replace the independent terms of paired cells with the bivariate normal's
	pc <- ctx$pc
	rho <- hyp$rho_d
	one <- pc$one
	two <- pc$two
	inf <- .lde_nrm_pair_info(ve, rho)
	r1 <- r[one]
	r2 <- r[two]
	psi1 <- inf$w * r1 + inf$c * r2
	psi2 <- inf$w * r2 + inf$c * r1
	out$w[c(one, two)] <- inf$w
	out$score[one] <- psi1
	out$score[two] <- psi2
	out$z[one] <- eta[one] + psi1 / inf$w
	out$z[two] <- eta[two] + psi2 / inf$w
	cw <- array(0, dim(y))
	cw[c(one, two)] <- inf$c
	out$cw <- cw
	v1 <- if (is.null(vn)) 0 else vn[one]
	v2 <- if (is.null(vn)) 0 else vn[two]
	# deviance of the paired cells, with the predictors' variances and the normalising constants
	q_ind <- (r1^2 + v1) / ve + (r2^2 + v2) / ve + 2 * log(2 * pi * ve)
	q_pair <- ((r1^2 + v1 + r2^2 + v2) - 2 * rho * r1 * r2) / (ve * (1 - rho^2)) + 2 * log(2 * pi * ve) + log(1 - rho^2)
	out$dev <- out$dev - sum(q_ind) + sum(q_pair)
	if (!is.null(dev_marg_cells <- out$dev_marg)) {
		# the marginal law of a pair: N((eta1, eta2), ve [1 rho; rho 1] + diag(v1, v2))
		s11 <- ve + v1
		s22 <- ve + v2
		s12 <- ve * rho
		det <- s11 * s22 - s12^2
		m_ind <- r1^2 / (ve + v1) + r2^2 / (ve + v2) + log(2 * pi * (ve + v1)) + log(2 * pi * (ve + v2))
		m_pair <- (s22 * r1^2 - 2 * s12 * r1 * r2 + s11 * r2^2) / det + log((2 * pi)^2 * det)
		out$dev_marg <- dev_marg_cells - sum(m_ind) + sum(m_pair)
	}
	out$pair <- c(pc[c("one", "two", "i", "j", "t")])
	out
}

# deviance of the cells at `eta`
.lde_dev_nrm <- function(ctx, eta, hyp) .lde_work_nrm(ctx, eta, hyp)$dev

# cavity weights: 1 / (ve + V) per cell, and for a pair the precision of the bivariate marginal law
# with the covariance `cvpair` of its two predictors
.lde_cavity_w_nrm <- function(ctx, eta, V3, lap, hyp) {
	ve <- hyp$ve
	obs <- ctx$obs
	w <- array(0, dim(ctx$y))
	w[obs] <- 1 / (ve + V3[obs])
	if (!isTRUE(ctx$dcor) || is.null(ctx$pc) || !ctx$pc$n_pairs || is.null(lap$cvpair)) return(list(w = w, cw = NULL))
	pc <- ctx$pc
	rho <- hyp$rho_d
	one <- pc$one
	two <- pc$two
	s11 <- ve + V3[one]
	s22 <- ve + V3[two]
	s12 <- ve * rho + lap$cvpair
	det <- s11 * s22 - s12^2
	w[one] <- s22 / det
	w[two] <- s11 / det
	cw <- array(0, dim(ctx$y))
	cw[one] <- -s12 / det
	cw[two] <- -s12 / det
	list(w = w, cw = cw)
}

# expected residual sum of squares under the posterior of the predictors (the statistic of the EM image
# of ve); pairs use the bivariate precision, and `vcell` NULL gives the plug-in at the point predictor
.lde_nrm_q_sum <- function(ctx, eta, vcell, cvpair, rho_d) {
	obs <- ctx$obs
	r <- ctx$y - eta
	r[!obs] <- 0
	v <- if (is.null(vcell)) 0 else vcell
	e2 <- r^2 + v
	e2[!obs] <- 0
	pc <- ctx$pc
	if (!isTRUE(ctx$dcor) || is.null(pc) || !pc$n_pairs || is.null(rho_d)) return(sum(e2))
	one <- pc$one
	two <- pc$two
	cv <- if (is.null(cvpair)) 0 else cvpair
	pair <- (e2[one] + e2[two] - 2 * rho_d * (r[one] * r[two] + cv)) / (1 - rho_d^2)
	sum(e2) - sum(e2[one]) - sum(e2[two]) + sum(pair)
}

# EM step of the dyadic correlation under the samplers' (1 - rho^2)^k prior: the best real root in (-1, 1)
# of -(M - 2k) ve rho^3 + S12 rho^2 + ((M - 2k) ve - S11) rho + S12 = 0, with the score `rs` and
# curvature `rc` of the expected log likelihood there; NULL without a pair
.lde_nrm_rho_step <- function(ctx, hyp, eta, lap, vextra = NULL) {
	pc <- ctx$pc
	if (!isTRUE(ctx$dcor) || is.null(pc) || !pc$n_pairs || !is.null(ctx$rho_hold)) return(NULL)
	obs <- ctx$obs
	r <- ctx$y - eta
	r[!obs] <- 0
	v <- lap$vcell
	if (!is.null(vextra)) v <- v + vextra
	one <- pc$one
	two <- pc$two
	cv <- if (is.null(lap$cvpair)) 0 else lap$cvpair
	S11 <- sum(r[one]^2 + v[one] + r[two]^2 + v[two])
	S12 <- sum(r[one] * r[two] + cv)
	M <- pc$n_pairs
	ve <- hyp$ve
	k <- .lame_rho_prior_power()
	lp <- function(rho) -((M - 2 * k) / 2) * log(1 - rho^2) - (S11 - 2 * rho * S12) / (2 * ve * (1 - rho^2))
	a <- (M - 2 * k) * ve
	rt <- polyroot(c(S12, a - S11, S12, -a))
	rt <- Re(rt[abs(Im(rt)) < 1e-8 * max(1, abs(rt))])
	rt <- rt[abs(rt) < 0.999]
	cand <- c(rt, -0.98, 0.98)
	val <- vapply(cand, lp, numeric(1))
	new <- cand[which.max(val)]
	# curvature of the expected log likelihood without the prior, by central difference of its derivative
	gr <- function(rho) M * rho / (1 - rho^2) + (S12 * (1 + rho^2) - rho * S11) / (ve * (1 - rho^2)^2)
	h <- 1e-5
	hi <- min(new + h, 0.9999)
	lo <- max(new - h, -0.9999)
	list(rho = new, rs = gr(new), rc = (gr(hi) - gr(lo)) / (hi - lo))
}

# EM image of ve at predictors `eta` with posterior variances `vcell` (and pair covariances `cvpair`);
# count panels use the statistic of lame_dynamic_als_eb_pln.R at the shifted predictor `eta + dshift`
.lde_ve_image_at <- function(ctx, eta, vcell, cvpair, dshift, hyp) {
	if (isTRUE(ctx$pln)) return(
		.lde_pln_ve_image(ctx, .lde_pln_ee_sum(ctx, if (is.null(dshift)) eta else eta + dshift, vcell, hyp$ve))
	)
	.lde_pln_ve_image(ctx, .lde_nrm_q_sum(ctx, eta, vcell, cvpair, hyp$rho_d))
}
