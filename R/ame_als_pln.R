# Poisson-lognormal cell of the empirical-bayes count estimator: y | z ~ Poisson(exp(z)), z = eta + e, e ~ N(0, ve)
# the noise is integrated out by an adaptive Gauss-Hermite rule at the posterior mode (derivation in src/ae_pln.h)

.ae_dims <- function(x) if (is.null(dim(x))) length(x) else dim(x)

# number of inner-rule nodes by the noise variance, chosen once per stage so the objective
# does not change with the rule as ve settles
.ae_pln_kin <- function(ve) {
	if (ve <= 0.1) 14L else if (ve <= 0.3) 18L else if (ve <= 0.5) 22L else if (ve <= 1) 28L else if (ve <= 1.5) 36L else if (ve <= 2) 44L else if (ve <= 3) 56L else 64L
}

# outer rules over a Gaussian predictor N(m, v): nine nodes above variance 0.3, fewer where the
# predictor is nearly known; `K` = 9 selects by variance, any other number fixes the rule
.ae_pln_outer <- function(K = 9L) {
	if (!identical(as.integer(K), 9L)) {
		g <- .ae_mf_gh(K)
		return(list(z = list(g$z), w = list(g$w), vcut = numeric()))
	}
	gs <- lapply(c(3L, 5L, 7L, 9L), .ae_mf_gh)
	list(z = lapply(gs, `[[`, "z"), w = lapply(gs, `[[`, "w"), vcut = c(0.004, 0.04, 0.3))
}

.ae_pln_rule_cache <- new.env(parent = emptyenv())

# nodes and log weights of the standard normal Gauss-Hermite rule of the inner integral
.ae_pln_rule <- function(K = 24L) {
	key <- as.character(K)
	hit <- .ae_pln_rule_cache[[key]]
	if (!is.null(hit)) return(hit)
	g <- .ae_mf_gh(K)
	out <- list(z = g$z, lw = log(g$w))
	.ae_pln_rule_cache[[key]] <- out
	out
}

# reference implementation of the cell (as pln_node() in src/ae_pln.h): log likelihood, score,
# information, its derivative, and the noise posterior's second moment and mean; processed in chunks
.ae_pln_node_r <- function(y, eta, ve, K = .ae_pln_kin(max(ve)), chunk = 20000L) {
	n <- length(y)
	ve <- rep_len(ve, n)
	out <- vector("list", ceiling(n / chunk))
	for (k in seq_along(out)) {
		ix <- ((k - 1L) * chunk + 1L):min(k * chunk, n)
		out[[k]] <- .ae_pln_node_block(y[ix], eta[ix], ve[ix], K)
	}
	if (!length(out)) return(list(
		l = numeric(),
		s = numeric(),
		w = numeric(),
		hp = numeric(),
		ee = numeric(),
		me = numeric()
	))
	if (length(out) == 1L) return(out[[1L]])
	nm <- names(out[[1L]])
	setNames(lapply(nm, function(a) unlist(lapply(out, `[[`, a), use.names = FALSE)), nm)
}

.ae_pln_node_block <- function(y, eta, ve, K) {
	n <- length(y)
	eta <- pmin(pmax(eta, -700), 700)
	ve <- pmin(pmax(ve, 1e-8), 1e6)
	iv <- 1 / ve
	ly <- ifelse(y > 0, log(y), -Inf)
	# posterior mode of the noise by a monotone Newton iteration from the right of the root
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
	# expm1(t) - t - t^2 / 2, by its series where the subtraction loses every digit
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
	# a node whose weight underflowed carries nothing
	ex[p == 0] <- 0
	ez1 <- rowSums(p * Z)
	d <- Z - ez1
	vz <- rowSums(p * d * d)
	k3 <- rowSums(p * d * d * d)
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
		ee = s * s * vz + me * me,
		me = me
	)
}

# the same by the compiled kernel (ae_pln_node_cpp); zero off the observed cells and where
# the predictor is not finite
.ae_pln_terms_r <- function(y, eta, ve, K = .ae_pln_kin(ve)) {
	ix <- which(is.finite(y) & is.finite(eta))
	nd <- .ae_pln_node_r(y[ix], eta[ix], ve, K)
	fill <- function(v) {
		o <- array(0, .ae_dims(y))
		o[ix] <- v
		o
	}
	list(
		l = fill(nd$l),
		s = fill(nd$s),
		w = fill(nd$w),
		hp = fill(nd$hp),
		ee = fill(nd$ee),
		me = fill(nd$me),
		ve = ve
	)
}

.ae_pln_terms <- function(y, eta, ve, K = .ae_pln_kin(ve)) {
	if (identical(.ae_engine(), "r")) return(.ae_pln_terms_r(y, eta, ve, K))
	rule <- .ae_pln_rule(K)
	out <- ae_pln_node_cpp(
		as.numeric(y),
		as.numeric(eta),
		as.numeric(ve),
		lgamma(as.numeric(y) + 1),
		rule$z,
		rule$lw,
		.ae_threads()
	)
	d <- .ae_dims(y)
	for (nm in names(out)) dim(out[[nm]]) <- d
	out$ve <- ve
	out
}

# expected working problem of cells with predictor N(m, v), as .ae_mf_cells_r() for a binary cell,
# plus `ee`, the expected second moment of the noise; the response increment is capped at 30
.ae_pln_cells_r <- function(y, m, v, ve, K = 9L, Kin = .ae_pln_kin(ve)) {
	obs <- is.finite(y) & is.finite(m)
	ix <- which(obs)
	mm <- m[ix]
	vv <- pmax(rep_len(v, length(y))[ix], 0)
	yy <- y[ix]
	rules <- .ae_pln_outer(K)
	rule_of <- findInterval(vv, rules$vcut, left.open = TRUE) + 1L
	es <- ew <- el <- ee <- numeric(length(ix))
	for (r in seq_along(rules$z)) {
		sel <- which(rule_of == r)
		if (!length(sel)) next
		sdv <- sqrt(vv[sel])
		for (k in seq_along(rules$z[[r]])) {
			nd <- .ae_pln_node_r(yy[sel], mm[sel] + sdv * rules$z[[r]][k], ve, Kin)
			es[sel] <- es[sel] + rules$w[[r]][k] * nd$s
			ew[sel] <- ew[sel] + rules$w[[r]][k] * nd$w
			el[sel] <- el[sel] + rules$w[[r]][k] * nd$l
			ee[sel] <- ee[sel] + rules$w[[r]][k] * nd$ee
		}
	}
	delta <- es / pmax(ew, 1e-12)
	big <- abs(delta) > 30
	if (any(big)) {
		delta[big] <- sign(delta[big]) * 30
		ew[big] <- abs(es[big]) / 30
	}
	z <- w <- ell <- eev <- array(0, .ae_dims(y))
	z[ix] <- mm + delta
	w[ix] <- ew
	ell[ix] <- el
	eev[ix] <- ee
	list(z = z, w = w, ell = ell, ee = eev)
}

.ae_pln_cells <- function(y, m, v, ve, K = 9L, Kin = .ae_pln_kin(ve)) {
	if (identical(.ae_engine(), "r")) return(.ae_pln_cells_r(y, m, v, ve, K, Kin))
	rule <- .ae_pln_rule(Kin)
	outer <- .ae_pln_outer(K)
	out <- ae_pln_cells_cpp(
		as.numeric(y),
		lgamma(as.numeric(y) + 1),
		as.numeric(m),
		as.numeric(v),
		as.numeric(ve),
		outer$z,
		outer$w,
		outer$vcut,
		rule$z,
		rule$lw,
		.ae_threads()
	)
	d <- .ae_dims(y)
	for (nm in names(out)) dim(out[[nm]]) <- d
	out
}

# cavity weight of the cells for the covariance of the factor rows (see .ae_mf_info_w_r): the
# information of a Poisson-lognormal cell with variance v + ve at mean m
.ae_pln_info_w_r <- function(y, m, v, ve, Kin = .ae_pln_kin(max(ve + pmax(v, 0), na.rm = TRUE))) {
	ix <- which(is.finite(y) & is.finite(m))
	w <- .ae_pln_node_r(y[ix], m[ix], ve + pmax(rep_len(v, length(y))[ix], 0), Kin)$w
	w[!is.finite(w) | w < 0] <- 0
	out <- array(0, .ae_dims(y))
	out[ix] <- w
	out
}

.ae_pln_info_w <- function(y, m, v, ve, Kin = .ae_pln_kin(max(ve + pmax(v, 0), na.rm = TRUE))) {
	if (identical(.ae_engine(), "r")) return(.ae_pln_info_w_r(y, m, v, ve, Kin))
	rule <- .ae_pln_rule(Kin)
	out <- ae_pln_info_w_cpp(
		as.numeric(y),
		lgamma(as.numeric(y) + 1),
		as.numeric(m),
		as.numeric(v),
		as.numeric(ve),
		rule$z,
		rule$lw,
		.ae_threads()
	)
	dim(out) <- .ae_dims(y)
	out
}

# variance of the noise

# prior of ve used by rs2_fc() for counts: scaled inverse chi-square, nu0 = 1 and s20 = 1
.ae_pln_ve_prior <- c(nu0 = 1, s20 = 1)

# EM image of ve, (nu0 s20 + SS) / (N + nu0 - 2), the posterior mean of its inverse-gamma full
# conditional; an undirected network counts its dyads once (`nu` = 1/2)
.ae_pln_ve_image <- function(ee_sum, n_cells, nu = 1, pr = .ae_pln_ve_prior, floor = 1e-8) {
	if (is.null(pr)) pr <- .ae_pln_ve_prior
	max((pr[["nu0"]] * pr[["s20"]] + nu * ee_sum) / (nu * n_cells + pr[["nu0"]] - 2), floor)
}

# hyperprior of ve on the log-likelihood scale, stationary at the EM image
.ae_pln_ve_logprior <- function(ve, pr = .ae_pln_ve_prior) {
	if (is.null(pr)) pr <- .ae_pln_ve_prior
	-0.5 * (pr[["nu0"]] - 2) * log(ve) - 0.5 * pr[["nu0"]] * pr[["s20"]] / ve
}

# Aitken delta-squared extrapolation of the slowly contracting EM image of log ve over `lag` steps,
# damped and bounded by `max_step`; NA until three increments agree in ratio
.ae_pln_aitken <- function(hist, lag = 5L, max_step = 0.5, damp = 0.5, agree = 0.05) {
	if (isFALSE(getOption("lame.als_aitken", TRUE))) return(NA_real_)
	n <- length(hist)
	if (n < 3L * lag + 1L) return(NA_real_)
	x3 <- hist[n]; x2 <- hist[n - lag]; x1 <- hist[n - 2L * lag]; x0 <- hist[n - 3L * lag]
	d2 <- x3 - x2; d1 <- x2 - x1; d0 <- x1 - x0
	if (!all(is.finite(c(d0, d1, d2))) || abs(d0) < 1e-8 || abs(d1) < 1e-8 || abs(d2) < 1e-8) return(NA_real_)
	r1 <- d2 / d1; r0 <- d1 / d0
	if (!(r1 > 0.3 && r1 < 0.97 && r0 > 0.3 && r0 < 0.97 && abs(r1 - r0) < agree)) return(NA_real_)
	max(min(damp * d2 * r1 / (1 - r1), max_step), -max_step)
}

# start of ve for a reweighting loop that estimates it by the EM step
.ae_pln_ve_start <- 0.5

# variance of the noise by the Laplace marginal likelihood

# the reported ve maximises the Laplace marginal likelihood of the posterior mode (type-II ML),
# found by continuation on log ve from a cold solve at a small ve, where weak factor columns stay alive,
# then a bracket and parabolic search; the variational stage holds ve at this value
.ae_pln_ve_lo <- 0.05
.ae_pln_mode_init <- function(fit, ve) {
	init <- list(
		mu = fit$mu,
		beta = fit$beta,
		a = fit$a,
		b = fit$b,
		U = fit$U,
		V = fit[["V"]],
		L = fit$L,
		sab = fit$sab,
		gamma = fit$gamma,
		ve = ve
	)
	init
}
# run(ve, from) solves the mode at the value ve held, cold when `from` is NULL and warm from the mode `from` otherwise
.ae_pln_ve_search <- function(run, ve0, from0 = NULL, local = FALSE, max_eval = 24L, ve_prior = .ae_pln_ve_prior) {
	if (is.null(ve_prior)) ve_prior <- .ae_pln_ve_prior
	xs <- numeric(0); Js <- numeric(0); irs <- list(); pass <- logical(0); its <- integer(0)
	objective <- function(ir, ve) {
		lm <- ir$fit$eb_vcov$evidence$log_marginal
		if (is.null(lm) || !is.finite(lm)) return(-Inf)
		lm + .ae_pln_ve_logprior(ve, ve_prior)
	}
	# solved point nearest to x (any point, if none has solved)
	nearest <- function(x) {
		cand <- if (any(pass)) which(pass) else seq_along(xs)
		irs[[cand[which.min(abs(xs[cand] - x))]]]
	}
	ev <- function(x, from) {
		x <- as.numeric(x)
		ir <- run(exp(x), from)
		# a warm solve that misses the certificate is repeated from the nearest solved point
		if (!isTRUE(ir$stationarity$pass) && any(pass)) {
			alt <- nearest(x)
			if (!identical(alt, from)) {
				ir2 <- run(exp(x), alt)
				if (isTRUE(ir2$stationarity$pass)) ir <- ir2
			}
		}
		xs <<- c(xs, x); Js <<- c(Js, objective(ir, exp(x))); irs[[length(irs) + 1L]] <<- ir
		pass <<- c(pass, isTRUE(ir$stationarity$pass)); its <<- c(its, as.integer(ir$irls_iters))
		length(xs)
	}
	h <- if (local) 0.2 else 0.4
	x0 <- log(ve0)
	ev(x0, from0)
	# bracket: step up while the objective rises, down if it fell at the first step
	j <- ev(x0 + h, irs[[1L]])
	grow <- function(dir) {
		step <- h
		repeat {
			if (length(xs) >= max_eval) break
			step <- step * 1.5
			edge <- if (dir > 0) which.max(xs) else which.min(xs)
			j <- ev(xs[edge] + dir * step, irs[[edge]])
			if (!(Js[j] > max(Js[-j]))) break
		}
	}
	if (Js[j] > Js[1L]) {
		grow(+1)
	} else {
		j <- ev(x0 - h, irs[[1L]])
		if (Js[j] > Js[1L]) grow(-1)
	}
	# parabolic (or golden-section) steps inside the bracket, each solved warm, until the parabola
	# through points within 0.5 of each other predicts a gain below 0.01 in the log likelihood
	for (k in seq_len(8L)) {
		if (length(xs) >= max_eval) break
		o <- order(xs); b <- which.max(Js[o])
		if (b == 1L || b == length(o) || !is.finite(Js[o][b])) break
		x3 <- xs[o][(b - 1L):(b + 1L)]; J3 <- Js[o][(b - 1L):(b + 1L)]
		if (!all(is.finite(J3))) break
		gl <- x3[2L] - x3[1L]; gr <- x3[3L] - x3[2L]
		d1 <- (J3[2L] - J3[1L]) / gl; d2 <- (J3[3L] - J3[2L]) / gr
		cc <- (d2 - d1) / (x3[3L] - x3[1L]) # the curvature / 2
		xv <- NA_real_; gain <- Inf
		if (is.finite(cc) && cc < 0 && max(gl, gr) <= 4 * min(gl, gr)) {
			xv <- (x3[1L] + x3[2L]) / 2 - d1 / (2 * cc)
			gain <- J3[1L] + d1 * (xv - x3[1L]) + cc * (xv - x3[1L]) * (xv - x3[2L]) - J3[2L]
			if (!is.finite(xv) || xv <= x3[1L] || xv >= x3[3L]) xv <- NA_real_
		}
		if (is.na(xv)) {
			xv <- if (gr >= gl) x3[2L] + 0.382 * gr else x3[2L] - 0.382 * gl
			gain <- Inf
		} else if ((gain < 0.01 || abs(xv - x3[2L]) < 0.01) && max(gl, gr) <= 0.5) break
		ev(xv, nearest(xv))
	}
	b <- which.max(Js)
	# return the maximiser with a solved mode: a solve stopped at its iteration limit is continued,
	# one that stalled is repeated cold
	through <- function(ir) {
		for (k in seq_len(4L)) {
			if (isTRUE(ir$stationarity$pass) || !identical(ir$stationarity$termination, "iteration limit")) break
			ir <- run(exp(xs[b]), ir)
		}
		ir
	}
	if (!isTRUE(irs[[b]]$stationarity$pass)) {
		res <- through(irs[[b]])
		if (!isTRUE(res$stationarity$pass)) res <- through(run(exp(xs[b]), NULL))
		if (isTRUE(res$stationarity$pass)) {
			irs[[b]] <- res; Js[b] <- objective(res, exp(xs[b])); pass[b] <- TRUE
		}
	}
	list(
		ir = irs[[b]],
		ve = as.numeric(exp(xs[b])),
		objective = as.numeric(Js[b]),
		curve = data.frame(
			ve = as.numeric(exp(xs)),
			objective = as.numeric(Js),
			solved = pass,
			iterations = its
		)[order(xs), ],
		evaluations = length(xs)
	)
}

# prior of the actor effects of a count network

# additive-effects prior of ame() and lame() for counts: scale from .ame_pois_effect_scale,
# df = round(4 + 3 n / 100) with n the (row) actors
.ae_eb_hyperprior_pois <- function(Yarr) {
	nact <- dim(Yarr)[1L]
	list(Sab0 = diag(2) * .ame_pois_effect_scale(Yarr), eta0 = round(4 + 3 * nact / 100))
}

# starting working response: log count plus the samplers' pseudo-count (the mean count, at most one)
.ae_pln_working_start <- function(Yarr) {
	ybar <- mean(Yarr, na.rm = TRUE)
	if (!is.finite(ybar) || ybar <= 0) ybar <- 1
	log(Yarr + min(ybar, 1))
}

# working problem of the mode

# Gaussian working response z = eta + s / w and weight w, with the increment capped at 30
# and the weight raised to |s| / 30 where the cap binds (.ae_irls_score_step)
.ae_irls_pln <- function(eta, terms, delta_cap = 30) {
	.ae_irls_score_step(eta, terms$s, terms$w, delta_cap)
}

# minus twice the log likelihood of the observed cells
.ae_dev_pln <- function(terms) -2 * sum(terms$l)
