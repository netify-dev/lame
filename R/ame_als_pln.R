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

# log likelihood, score, information, its derivative and the noise posterior's moments of each cell (as
# pln_node() in src/ae_pln.h) by the compiled kernel (ae_pln_node_cpp); zero off the observed cells and where
# the predictor is not finite
.ae_pln_terms <- function(y, eta, ve, K = .ae_pln_kin(ve)) {
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

# expected working problem of cells with predictor N(m, v), as .ae_mf_cells() for a binary cell, plus `ee`, the
# expected second moment of the noise; the response increment is capped at 30
.ae_pln_cells <- function(y, m, v, ve, K = 9L, Kin = .ae_pln_kin(ve)) {
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

# cavity weight of the cells for the covariance of the factor rows (see .ae_mf_info_w): the information of a
# Poisson-lognormal cell with variance v + ve at mean m
.ae_pln_info_w <- function(y, m, v, ve, Kin = .ae_pln_kin(max(ve + pmax(v, 0), na.rm = TRUE))) {
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

# prior of ve for counts: scaled inverse chi-square on s20. a scale the user gives carries one
# degree of freedom, as rs2_fc() does; without one the prior has none (nu0 = 0), so no scale
# pulls ve, which the likelihood of a count network barely informs
.ae_pln_ve_prior <- c(nu0 = 1, s20 = 1)
.ae_pln_ve_flat <- c(nu0 = 0, s20 = 1)
.ae_pln_ve_prior_of <- function(s20 = NULL, nu0 = NULL) {
	if (is.null(s20)) return(.ae_pln_ve_flat)
	c(nu0 = if (is.null(nu0)) .ae_pln_ve_prior[["nu0"]] else as.numeric(nu0), s20 = as.numeric(s20))
}

# EM image of ve, (nu0 s20 + SS) / (N + nu0 - 2), the posterior mean of its inverse-gamma full
# conditional; an undirected network counts its dyads once (`nu` = 1/2). a network of two or
# fewer cells has no such mean, and the image is then the ML one, SS / N
.ae_pln_ve_image <- function(ee_sum, n_cells, nu = 1, pr = .ae_pln_ve_flat, floor = 1e-8) {
	if (is.null(pr)) pr <- .ae_pln_ve_flat
	den <- nu * n_cells + pr[["nu0"]] - 2
	if (!(den > 0)) den <- max(nu * n_cells, 1)
	max((pr[["nu0"]] * pr[["s20"]] + nu * ee_sum) / den, floor)
}

# hyperprior of ve on the log-likelihood scale, stationary at the EM image
.ae_pln_ve_logprior <- function(ve, pr = .ae_pln_ve_flat) {
	if (is.null(pr)) pr <- .ae_pln_ve_flat
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
# smallest ve searched: the log density of the prior already holds the maximum off zero, and
# lognormal noise of this variance is not distinguishable from none in any count network
.ae_pln_ve_min <- 1e-6
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
.ae_pln_ve_search <- function(run, ve0, from0 = NULL, local = FALSE, max_eval = 24L, ve_prior = .ae_pln_ve_flat) {
	if (is.null(ve_prior)) ve_prior <- .ae_pln_ve_flat
	xs <- numeric(0); Js <- numeric(0); irs <- list(); pass <- logical(0); its <- integer(0)
	base <- numeric(0); wl <- list(); tl <- list()
	# the Laplace evidence plus the log density of the prior, and the second-order part of the
	# variance of the product cells the evidence omits (.ae_pln_second_order)
	objective <- function(ir, ve) {
		ev <- ir$fit$eb_vcov$evidence
		lm <- ev$log_marginal
		so <- ev$second_order
		ok <- is.list(so) && all(is.finite(so$w)) && all(is.finite(so$t))
		list(
			base = if (is.null(lm) || !is.finite(lm)) -Inf else lm + .ae_pln_ve_logprior(ve, ve_prior),
			w = if (ok) so$w else NULL,
			t = if (ok) so$t else NULL
		)
	}
	# each point is scored with that variance held at its value at the best point, -sum w t_best,
	# so the scores differ only through the cell information w: the maximum is then the EM fixed
	# point whose second moment of the noise carries the full posterior variance of the product
	rescore <- function() {
		b <- which.max(if (length(Js) == length(base)) Js else base)
		for (pass_k in 1:2) {
			tb <- tl[[b]]
			if (is.null(tb)) {
				Js <<- base
				return(invisible())
			}
			Js <<- base - vapply(seq_along(base), function(k) if (is.null(wl[[k]])) 0 else sum(wl[[k]] * tb), 0)
			b <- which.max(Js)
		}
	}
	record <- function(k, ir, ve) {
		o <- objective(ir, ve)
		base[k] <<- o$base
		wl[k] <<- list(o$w)
		tl[k] <<- list(o$t)
		rescore()
	}
	# solved point nearest to x (any point, if none has solved)
	nearest <- function(x) {
		cand <- if (any(pass)) which(pass) else seq_along(xs)
		irs[[cand[which.min(abs(xs[cand] - x))]]]
	}
	ev <- function(x, from) {
		x <- max(as.numeric(x), log(.ae_pln_ve_min))
		ir <- run(exp(x), from)
		# a warm solve that misses the certificate is repeated from the nearest solved point
		if (!isTRUE(ir$stationarity$pass) && any(pass)) {
			alt <- nearest(x)
			if (!identical(alt, from)) {
				ir2 <- run(exp(x), alt)
				if (isTRUE(ir2$stationarity$pass)) ir <- ir2
			}
		}
		xs <<- c(xs, x); irs[[length(irs) + 1L]] <<- ir
		record(length(xs), ir, exp(x))
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
			irs[[b]] <- res; record(b, res, exp(xs[b])); pass[b] <- TRUE
			b <- which.max(Js)
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

# additive-effects prior of a count fit: the scale of ame() and lame() (.ame_pois_effect_scale) with
# the weight of one actor. the log row and column rates that scale is read from carry the latent
# product and the count noise, so it is several times the effect variance under a strong product,
# and the sampler's weight of round(4 + 3 n / 100) actors pulls the estimate towards it
.ae_eb_hyperprior_pois <- function(Yarr) {
	list(Sab0 = diag(2) * .ame_pois_effect_scale(Yarr), eta0 = 1)
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
