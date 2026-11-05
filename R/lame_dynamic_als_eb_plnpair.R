# dyadic correlation of a directed count panel at the fitted state of the dynamic Poisson EB
# estimator. the two cells of a dyad have bivariate normal noise with variance ve and correlation
# rho, so a pair of counts has a Poisson-lognormal pair probability (see ame_als_pln_pair.R).
# effects, factors and ve come from the independent-cell fit; this file adds rho and the pair
# information used for the coefficient covariance

# covariance of the two predictors of each pair (i < j, period t) from the additive effects:
# cov(a_i + b_j, a_j + b_i); the factors of two actors are independent under q
.lde_pair_cov <- function(ctx, cond) {
	pc <- ctx$pc
	na <- ctx$na
	cv <- numeric(pc$n_pairs)
	for (t in unique(pc$t)) {
		ip <- which(pc$t == t)
		C <- cond[[t]]
		i <- pc$i[ip]
		j <- pc$j[ip]
		cv[ip] <- C[cbind(i, j)] + C[cbind(i, na + i)] + C[cbind(na + j, j)] + C[cbind(na + j, na + i)]
	}
	cv
}

# posterior mean and sd of the dyadic correlation under the sampler prior, from the plug-in
# pair likelihood: grid in asin(rho), Newton refinement, then a Gaussian or fine-grid posterior
.lde_pln_pair_rho <- function(y1, y2, m1, m2, ve) {
	Kp <- .ae_pln2_kin(ve)
	ll <- function(rho) sum(.ae_pln2_node(y1, y2, m1, m2, ve, ve, rho * ve, Kp)$L)
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
		.ae_pair_rho_post(rho_c, prec)
	} else {
		fine <- seq(-1.4, 1.4, length.out = 2801L)
		lpf <- stats::splinefun(th, lg, method = "fmm")(fine)
		rf <- sin(fine)
		lw <- lpf + .lame_rho_prior_log(rf, .lame_rho_prior_power()) + log(cos(fine))
		wt <- exp(lw - max(lw))
		mu <- sum(wt * rf) / sum(wt)
		list(mean = mu, sd = sqrt(sum(wt * (rf - mu)^2) / sum(wt)))
	}
	list(rho = max(min(post$mean, 0.98), -0.98), sd = post$sd, peak = rho_c, prec = prec)
}

# cell weights under dyadic correlation `rho`: the pair cells get the observed information of
# the pair probability integrated over the predictors (`V3` their variances, `cvv` their
# covariance), the other cells the single-cell cavity weight
.lde_cavity_w_pln_pair <- function(ctx, eta, V3, cvv, ve, rho) {
	pc <- ctx$pc
	w <- .ae_pln_info_w(ctx$y, eta, V3, ve)
	cw <- array(0, dim(ctx$y))
	if (pc$n_pairs > 0L) {
		m1 <- eta[pc$one]; m2 <- eta[pc$two]
		v1 <- V3[pc$one]; v2 <- V3[pc$two]
		Kw <- .ae_pln2_kin(max(v1, v2) + ve)
		nd <- .ae_pln2_node(ctx$y[pc$one], ctx$y[pc$two], m1, m2, v1 + ve, v2 + ve, cvv + rho * ve, Kw)
		w[pc$one] <- pmax(nd$w11, 0)
		w[pc$two] <- pmax(nd$w22, 0)
		cw[pc$one] <- nd$w12
		cw[pc$two] <- nd$w12
	}
	list(w = w, cw = cw)
}

# dyadic correlation of a fitted state with the information of its pairs; NULL when
# fewer than three dyads are observed in both directions
.lde_pln_pair_stage <- function(ctx, fit) {
	pc <- ctx$pc
	if (is.null(pc) || pc$n_pairs < 3L) return(NULL)
	st <- fit$st
	ev <- fit$ev
	ve <- st$hyp$ve
	eta <- ev$eta
	rr <- .lde_pln_pair_rho(ctx$y[pc$one], ctx$y[pc$two], eta[pc$one], eta[pc$two], ve)
	c(rr, list(n_pairs = pc$n_pairs, cvv = .lde_pair_cov(ctx, ev$post$cond), method = "integrated pair likelihood"))
}

# coefficient covariance of a fit without static factors, from the system at the pair
# weights; `pr` is the result of .lde_pln_pair_stage()
.lde_pln_pair_vcov <- function(ctx, st, ev, pr) {
	q <- ctx$q
	m <- ctx$mth
	sta <- st
	sta$O <- NULL
	Vp <- if (is.null(st$Vpos)) 0 else st$Vpos
	cw <- .lde_cavity_w_pln_pair(ctx, ev$eta, Vp + ev$lap$vcell, pr$cvv, st$hyp$ve, pr$rho)
	wk_i <- ev$wk
	wk_i$w <- cw$w
	wk_i$cw <- cw$cw
	wk_i$z <- ev$wk$z - (if (is.null(st$O)) 0 else st$O)
	sm <- .lde_system(ctx, sta, wk_i)
	fj <- .lde_factor(sm$A, sm$b, sm$C, sm$Hzz, q, q, sm$cdiag)
	fj$Szi[q + seq_len(m), q + seq_len(m), drop = FALSE]
}
