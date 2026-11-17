# the warning of the sampler for a binary response with values other than 0/1
.ae_warn_thresholded <- function() {
	cli::cli_warn(
		c(
			"{.arg family} = {.val binary} but {.arg Y} contains values other than 0/1.",
			"i" = "{.arg Y} will be thresholded to {.code 1 * (Y > 0)}; if you meant counts, use {.val poisson}, or {.val normal} as appropriate (ordinal outcomes take {.fn ame} or {.fn lame} with {.code method = \"mcmc\"})."
		),
		class = "lame_als_thresholded"
	)
}

# the response as given (matrix, array or list of matrices): a binary one with values other than
# 0/1 is thresholded at zero, as ame() and lame() do it
.ae_binary_input <- function(Y, family) {
	if (!identical(family, "binary")) return(Y)
	vals <- if (is.list(Y)) unlist(lapply(Y, as.vector), use.names = FALSE) else as.vector(Y)
	if (!(is.numeric(vals) || is.logical(vals)) || all(vals[is.finite(vals)] %in% c(0, 1))) return(Y)
	.ae_warn_thresholded()
	cut0 <- function(y) { y[] <- 1 * (y > 0); y }
	if (is.list(Y)) lapply(Y, cut0) else cut0(Y)
}

# values that do not belong to the stated observation family, on the prepared array: a binary
# response is thresholded as above, a negative count is refused. returns the response
.ae_check_als_response <- function(y, family) {
	observed_y <- y[is.finite(y)]
	if (family == "binary" && !all(observed_y %in% c(0, 1))) {
		.ae_warn_thresholded()
		y[] <- 1 * (y > 0)
	}
	if (family == "poisson" && any(observed_y < 0)) {
		cli::cli_abort(
			"family = {.val poisson} requires non-negative counts.",
			class = "lame_als_data"
		)
	}
	y
}

# iteratively reweighted least squares (IRLS) for the binary (logit / probit)
# and Poisson (log) families of the ALS estimator: a Gaussian working response
# and weight are rebuilt from the current linear predictor and the weighted BCD
# core is refit, giving an approximate GLM fit

# linear predictor clamped to [-eta_cap, eta_cap] (no clamp at the default Inf)
.ae_irls_cap <- function(eta, eta_cap) {
	if (is.finite(eta_cap)) pmin(pmax(eta, -eta_cap), eta_cap) else eta
}

# each map returns the Gaussian working response z = eta + delta and weight w

# Poisson, log link:  z = eta + (y - mu)/mu ,  w = mu = exp(eta)
.ae_irls_poisson <- function(
	y,
	eta,
	eta_cap = Inf,
	wfloor = 1e-6,
	delta_cap = 30
) {
	mu <- exp(.ae_irls_cap(eta, eta_cap))
	.ae_irls_score_step(eta, y - mu, pmax(mu, wfloor), delta_cap)
}

# working response from the score s and Fisher weight w: delta = s / w. where
# delta is capped at `delta_cap` the weight becomes |s| / delta_cap, so w * delta
# still equals the score and the fixed point solves the score equations
.ae_irls_score_step <- function(eta, s, w, delta_cap) {
	delta <- s / w
	big <- !is.na(delta) & abs(delta) > delta_cap
	if (any(big)) {
		delta[big] <- sign(delta[big]) * delta_cap
		w[big] <- abs(s[big]) / delta_cap
	}
	list(z = eta + delta, w = w)
}

# binary, logit link:  s = y - p ,  Fisher weight p (1 - p)
.ae_irls_binary_logit <- function(
	y,
	eta,
	eta_cap = Inf,
	wfloor = 1e-8,
	delta_cap = 30
) {
	ec <- .ae_irls_cap(eta, eta_cap)
	p <- stats::plogis(ec)
	q <- stats::plogis(-ec)
	.ae_irls_score_step(
		eta,
		y * q - (1 - y) * p,
		pmax(p * q, wfloor),
		delta_cap
	)
}

# binary, probit link:  s = y phi / Phi - (1 - y) phi / (1 - Phi) ,  Fisher
# weight phi^2 / (Phi (1 - Phi)), from log-scale Mills ratios for exact tails
.ae_irls_binary_probit <- function(
	y,
	eta,
	eta_cap = Inf,
	wfloor = 1e-8,
	delta_cap = 30,
	terms = NULL
) {
	if (is.null(terms)) terms <- .ae_probit_terms(.ae_irls_cap(eta, eta_cap), NULL)
	lam1 <- exp(terms$ld - terms$lp1)
	lam0 <- exp(terms$ld - terms$lp0)
	.ae_irls_score_step(
		eta,
		y * lam1 - (1 - y) * lam0,
		pmax(lam1 * lam0, wfloor),
		delta_cap
	)
}

# normal log density and both log tail probabilities at a probit predictor, computed once per iterate by the
# compiled kernel (ae_probit_terms_cpp); `outcome` is the log probability of each cell's observed outcome and
# `one` the cells with y = 1 (NULL: not needed)
.ae_probit_terms <- function(eta, one) {
	tm <- ae_probit_terms_cpp(eta, .ae_threads())
	d <- dim(eta)
	if (!is.null(d)) {
		dim(tm$ld) <- d; dim(tm$lp1) <- d; dim(tm$lp0) <- d
	}
	outcome <- NULL
	if (!is.null(one)) {
		outcome <- tm$lp0
		outcome[one] <- tm$lp1[one]
	}
	tm$outcome <- outcome
	tm
}

# family deviances (the IRLS ranking objective)

.ae_dev_poisson <- function(y, mu) {
	mu <- pmax(mu, 1e-12)
	term <- ifelse(y == 0 | !is.finite(y), 0, y * log(y / mu))
	2 * sum((term - (y - mu))[is.finite(y)])
}

.ae_dev_binary <- function(y, p) {
	p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
	ok <- is.finite(y)
	-2 * sum((y * log(p) + (1 - y) * log(1 - p))[ok])
}

# minus twice the observed-cell log likelihood of gaussian or poisson-lognormal cell terms
.ae_dev_terms <- function(terms) -2 * sum(terms$l)

# deviance from the predictor, so tails keep their exact log probabilities and
# the Poisson form avoids cancellation between large counts
.ae_irls_deviance <- function(y, eta, family, link, terms = NULL) {
	observed <- is.finite(y)
	if (!is.null(terms) && family == "binary" && link == "probit") {
		return(-2 * sum(terms$outcome[observed]))
	}
	y <- y[observed]
	eta <- eta[observed]
	if (family == "poisson") {
		term <- exp(eta)
		positive <- y > 0
		delta <- eta[positive] - log(y[positive])
		term[positive] <- y[positive] * (expm1(delta) - delta)
		return(2 * sum(term))
	}
	if (link == "probit") {
		return(-2 * sum(stats::pnorm((2 * y - 1) * eta, log.p = TRUE)))
	}
	2 * sum(log1p(exp(-abs(eta))) + pmax(-(2 * y - 1) * eta, 0))
}

# deviance of a reported static fit (fit$deviance): binomial or Poisson deviance
# of the observed cells at the reported fitted values for an IRLS fit, otherwise
# the residual sum of squares at the reported EZ. the optimiser's own objective
# is kept as fit$objective
.ae_fit_deviance <- function(
	Yarr,
	fitted,
	resid,
	family,
	link,
	eta = NULL,
	latent_var = 0
) {
	if (is.null(link) || !(family %in% c("binary", "poisson"))) {
		return(sum(vapply(resid, function(r) sum(r[is.finite(r)]^2), numeric(1))))
	}
	d <- 0
	for (t in seq_along(fitted)) {
		y <- Yarr[, , t]
		m <- fitted[[t]]
		ok <- is.finite(y) & is.finite(m)
		if (!is.null(eta)) {
			predictor <- eta[[t]]
			if (family == "poisson") predictor <- predictor + latent_var / 2
			d <- d + .ae_irls_deviance(y, predictor, family, link)
		} else {
			d <- d + if (identical(family, "binary")) .ae_dev_binary(y[ok], m[ok]) else .ae_dev_poisson(y[ok], m[ok])
		}
	}
	d
}

# per-coordinate scale of the factor-block ridge on the IRLS path: 1 for binary
# (probit latent scale), the moment-scaled row/column-mean variance of log1p(y)
# for Poisson. the EB fit of a directed unipartite network takes the prior of
# ame() and lame() (R/uv_prior.R) instead
.ae_uv_prior_scale <- function(y, family) {
	# a normal response is in its own unit (R/normal_scale.R)
	if (identical(family, "normal")) return(.lame_unit_sd(y, "normal"))
	if (family != "poisson") return(1)
	Ymom <- log1p(pmax(y, 0))
	Ec <- Ymom - mean(Ymom, na.rm = TRUE)
	a_tmp <- apply(Ec, 1, mean, na.rm = TRUE); a_tmp[is.na(a_tmp)] <- 0
	b_tmp <- apply(Ec, 2, mean, na.rm = TRUE); b_tmp[is.na(b_tmp)] <- 0
	vscale <- mean(c(stats::var(a_tmp), stats::var(b_tmp)), na.rm = TRUE)
	if (!is.finite(vscale) || vscale <= 0) {
		vscale <- max(stats::var(c(Ymom), na.rm = TRUE), 1e-12)
	}
	vscale
}

# BCD-fit-shaped state built from a warm start, with factor scales seeded from
# the incoming factors so its penalized deviance is comparable to the candidates
.ae_irls_warm_state <- function(init, nr, nc, R, symmetric, uv_prior) {
	U <- init$U; V <- init[["V"]]; L <- init$L
	Omat <- matrix(0, nr, nc)
	if (R > 0 && !is.null(U)) {
		if (symmetric && !is.null(L)) {
			Omat <- U %*% (L * t(U))
		} else if (!is.null(V)) {
			Omat <- tcrossprod(U, V)
		}
	}
	sig2_u <- sig2_v <- NULL
	if (R > 0 && !is.null(uv_prior) && !is.null(U)) {
		if (symmetric) {
			s2 <- (uv_prior$kappa0 * uv_prior$suv0 + 2 * sum(abs(L))) /
				(uv_prior$kappa0 + 2 * nr * R)
			sig2_u <- sig2_v <- rep(s2, R)
		} else {
			sig2_u <- (uv_prior$suv0 + colSums(U^2)) / (nr + uv_prior$kappa0)
			sig2_v <- (uv_prior$suv0 + colSums(V^2)) / (nc + uv_prior$kappa0)
		}
	}
	out <- list(
		mu = init$mu,
		beta = init$beta,
		a = init$a,
		b = init$b,
		U = U,
		V = V,
		L = L,
		Omat = Omat,
		sig2_u = sig2_u,
		sig2_v = sig2_v,
		sse = NA_real_,
		iterations = 0L,
		converged = FALSE,
		add_converged = TRUE,
		dev_history = numeric(0)
	)
	# an EB state keeps its actor-effect covariance, node-covariate coefficients
	# and, for counts, the noise variance
	for (nm in c("sab", "gamma", "ve")) if (!is.null(init[[nm]])) out[[nm]] <- init[[nm]]
	out
}

# damped step that stays a rank-R state, with product and penalty evaluated at
# the blended factors
.ae_irls_blend <- function(f0, f1, step, symmetric) {
	out <- f1
	for (nm in c(
		"mu", "beta", "a", "b", "U", "V", "L", "Omat", "sig2_u", "sig2_v",
		"sab", "gamma", "ve"
	)) {
		if (!is.null(f0[[nm]]) && !is.null(f1[[nm]]) && length(f0[[nm]]) == length(f1[[nm]])) {
			out[[nm]] <- (1 - step) * f0[[nm]] + step * f1[[nm]]
		}
	}
	if (!is.null(out$U)) {
		if (symmetric) {
			eg <- eigen((out$Omat + t(out$Omat)) / 2, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(ncol(out$U))]
			out$U <- eg$vectors[, idx, drop = FALSE]
			out[["V"]] <- out$U
			out$L <- eg$values[idx]
			out$Omat <- out$U %*% (out$L * t(out$U))
		} else {
			out$Omat <- tcrossprod(out$U, out[["V"]])
		}
	}
	out
}

# IRLS outer loop
#   prep       canonical prep list (.ae_prepare / .ae_prep_from_canonical)
#   init       optional warm-start state; when NULL the loop is seeded by the
#              one-shot transform fit
#   warm       start the reweighting at `init`
#   ab_prior   the EB actor-effect prior (.ae_eb_prior)
#   warm_retry rerun a warm start that fails the certificate from the spectral start
#   ve_fixed   count network under EB: noise variance held at this value (used by
#              .ae_pln_ve_search)
#   vcov_mode  "full" (covariance and evidence), "coef" (the same without the
#              standard errors of the effects and product, for a mode a
#              variational stage refines), "evidence" (log determinant only) or "none"
# returns list(fit, Z, W, dev_trace, irls_iters, best_obj, data_deviance, eta,
# stationarity)
.ae_irls_run <- function(
	prep,
	family,
	link,
	R,
	symmetric,
	max_iter,
	tol,
	lowrank_method = "mm",
	irls_maxit = NULL,
	irls_tol = NULL,
	init = NULL,
	linear_solver = "eigen",
	warm = FALSE,
	ab_prior,
	ve_fixed = NULL,
	warm_retry = TRUE,
	vcov_mode = c("full", "coef", "evidence", "none")
) {
	vcov_mode <- match.arg(vcov_mode)
	if (is.null(irls_maxit)) irls_maxit <- max_iter
	if (is.null(irls_tol)) irls_tol <- 1e-9

	Tt <- prep$Tt; nr <- prep$nr; nc <- prep$nc
	y <- .ae_check_als_response(prep$Yarr, family)

	unobserved <- !is.finite(y)

	# dyadic design as a cell-by-covariate matrix, cells in array order
	xmat <- if (prep$p_dyad > 0) {
		matrix(aperm(prep$X, c(1L, 2L, 4L, 3L)), nr * nc * Tt, prep$p_dyad)
	}

	# linear predictor array from a BCD fit
	eta_of <- function(f) {
		E <- array(f$mu + outer(f$a, f$b, "+") + f$Omat, dim = c(nr, nc, Tt))
		if (!is.null(xmat)) E <- E + as.vector(xmat %*% f$beta)
		E
	}

	# a count network under EB is the Poisson-lognormal model (ame_als_pln.R),
	# with noise variance ve as a hyperparameter
	pln <- identical(family, "poisson")
	# a normal network under EB has the Gaussian cell of ame_als_gauss.R: working
	# response y and constant weight 1 / ve
	gauss <- identical(family, "normal")
	has_ve <- pln || gauss
	work_fun <- switch(
		paste(family, link),
		"poisson log" = function(eta, terms = NULL) .ae_irls_pln(eta, terms),
		"binary logit" = function(eta, terms = NULL) .ae_irls_binary_logit(y, eta),
		"binary probit" = function(eta, terms = NULL)
			.ae_irls_binary_probit(y, eta, terms = terms),
		"normal identity" = function(eta, terms = NULL) .ae_irls_gauss(y, terms),
		cli::cli_abort("IRLS is not available for family {.val {family}}.")
	)
	dev_fun <- function(eta, terms = NULL)
		if (pln || gauss) .ae_dev_terms(terms) else .ae_irls_deviance(y, eta, family, link, terms)
	# probit terms of an iterate serve its working response, deviance and certificate
	probit <- identical(family, "binary") && identical(link, "probit")
	one_cells <- if (probit) which(y == 1)
	terms_of <- function(eta, f = NULL) {
		if (pln) .ae_pln_terms(y, eta, f$ve) else if (gauss) .ae_gauss_terms(y, eta, f$ve) else if (probit) .ae_probit_terms(eta, one_cells)
	}

	# MAP prior for the factor block (see ame_als_lowrank.R): the unpenalized
	# rank-R GLM fit does not exist under quasi-separation, so the ridge keeps the
	# factors finite
	kappa0 <- 2 * R + 2
	suv0 <- .ae_uv_prior_scale(y, family)
	# the EB fit of a directed unipartite network uses the uv prior of ame() and
	# lame() (R/uv_prior.R); other penalised fits keep the stabilising ridge
	if (R > 0 && !symmetric && !isTRUE(prep$bip)) {
		uvp <- .lame_uv_prior(y, family, R)
		kappa0 <- uvp$df
		suv0 <- uvp$psi
	}
	uv_prior <- if (R > 0) list(suv0 = suv0, kappa0 = kappa0) else NULL

	# an EB step runs a few BCD sweeps, then refreshes weights, Laplace correction
	# and EM scale images, so the loop does not stall on stale weights
	bcd_iter <- min(max_iter, .ae_eb_sweeps)
	# EM state a reweighted solve starts from: incoming actor-effect covariance and
	# sums of h' at its predictor
	eb_state <- function(f, eta, terms) {
		st <- ab_prior
		if (!is.null(f$sab)) st$sab <- f$sab
		st$gamma <- f$gamma
		# the information of a Gaussian cell does not depend on the predictor
		st$hprime <- if (gauss) NULL else .ae_eb_hprime(eta, link, terms)
		if (!gauss) st$hprime[unobserved] <- 0
		st
	}
	# EM image of the noise variance at an iterate (the E-step)
	ve_next <- function(terms) {
		if (!is.null(ve_fixed)) return(ve_fixed)
		.ae_pln_ve_image(sum(terms$ee), sum(!unobserved), ab_prior$nu, ab_prior$ve_prior, .ae_ve_floor(ab_prior))
	}
	# the EM image of a weakly informed variance contracts slowly, so accepted
	# iterates are extrapolated to their geometric limit (Aitken, .ae_pln_aitken)
	ve_hist <- numeric(0)
	# the EB objective is the Laplace-marginal one: actor-effect penalty and log
	# determinant of their information at the iterate's Fisher weights
	pdev_of <- function(f, scales = f, dev = NULL, work = NULL) {
		if (is.null(dev)) dev <- dev_fun(eta_of(f))
		dev <- dev + .ae_eb_penalty(
			f,
			scales$sab,
			ab_prior,
			wsum = rowSums(work$w, dims = 2L, na.rm = TRUE)
		)
		if (is.null(uv_prior) || is.null(f$U) || is.null(f$sig2_u)) {
			return(dev)
		}
		if (symmetric) {
			# EB scale update (.ae_lowrank_em_mm) on the scale of the ordered cells
			s2 <- max(scales$sig2_u[1], 1e-12)
			dev + (2 * suv0 + 2 * sum(abs(f$L))) / s2 +
				(2 * nr * R + 2 * kappa0) * log(s2)
		} else {
			su <- pmax(scales$sig2_u, 1e-12)
			sv <- pmax(scales$sig2_v, 1e-12)
			dev + sum((suv0 + colSums(f$U^2)) / su) +
				(nr + kappa0) * sum(log(su)) +
				sum((suv0 + colSums(f[["V"]]^2)) / sv) +
				(nc + kappa0) * sum(log(sv))
		}
	}

	# certificate of an iterate and its scalar merit (the largest thresholded
	# score), which accept a step
	cert_of <- function(f, eta, work, terms) {
		.ae_irls_stationarity(
			prep,
			f,
			family,
			link,
			R,
			symmetric,
			tol,
			eta = eta,
			xmat = xmat,
			working = work,
			terms = terms,
			ab_prior = ab_prior,
			ve_fixed = !is.null(ve_fixed)
		)
	}
	merit_of <- function(cert) {
		if (!all(is.finite(c(cert$max_score, cert$scale_change)))) return(Inf)
		max(cert$max_score, cert$scale_change, cert$unused_rank_score)
	}

	# seed: warm start, or the one-shot transform fit with the same penalty as
	# the reweighted fits
	Z0 <- if (pln) .ae_pln_working_start(y) else .ae_working_response(y, family)
	if (symmetric) for (t in seq_len(Tt)) Z0[, , t] <- .ae_symmetrize(Z0[, , t])
	warm <- isTRUE(warm) && !is.null(init)
	# a normal network starts at the fixed-effects fit and the variance of its
	# residuals (.ae_gauss_seed)
	ve_seed <- NULL
	if (gauss) {
		gs <- .ae_gauss_seed(prep, Z0, R, symmetric, max_iter, tol, init, lowrank_method, linear_solver, ve_fixed)
		init <- gs$init
		ve_seed <- gs$ve
	}
	fit <- if (warm) {
		.ae_irls_warm_state(init, nr, nc, R, symmetric, uv_prior)
	} else {
		.ae_bcd_fit(
			Z0,
			prep$X,
			R,
			symmetric,
			bcd_iter,
			tol,
			init = init,
			lowrank_method = lowrank_method,
			weights = if (gauss) array(1 / ve_seed, dim(y)),
			linear_solver = linear_solver,
			uv_prior = uv_prior,
			ab_prior = ab_prior
		)
	}
	if (has_ve) fit$ve <- if (!is.null(ve_fixed)) ve_fixed else if (gauss) ve_seed else if (!is.null(init$ve)) init$ve else .ae_pln_ve_start

	# predictor, deviance and working response of the current iterate, computed
	# once per accepted iterate
	eta_fit <- eta_of(fit)
	terms_fit <- terms_of(eta_fit, fit)
	dev_fit <- dev_fun(eta_fit, terms_fit)
	work_fit <- work_fun(eta_fit, terms_fit)
	best_fit <- fit
	best_dev <- pdev_of(fit, dev = dev_fit, work = work_fit)
	merit_fit <- merit_of(cert_of(fit, eta_fit, work_fit, terms_fit))
	merit_best <- merit_fit
	dev_trace <- best_dev
	n_accepted <- 0L
	conditional_gain <- numeric()
	termination <- "iteration limit"
	inner_tol <- min(tol, 1e-8)
	inner_refinements <- 0L
	certificate_now <- NULL

	for (k in seq_len(irls_maxit)) {
		# with a noise variance, move it to its EM image first and form the
		# reweighted problem at that value
		fit_w <- fit; terms_w <- terms_fit; work_w <- work_fit
		if (has_ve) {
			ve_new <- ve_next(terms_fit)
			fit_w$ve <- ve_new
			terms_w <- terms_of(eta_fit, fit_w)
			work_w <- work_fun(eta_fit, terms_w)
		}
		zk <- work_w$z; wt <- work_w$w
		zk[unobserved] <- NA_real_ # keep the observed pattern
		wt[unobserved] <- NA_real_
		if (symmetric) for (t in seq_len(Tt)) {
			zk[, , t] <- .ae_symmetrize(zk[, , t])
			wt[, , t] <- .ae_symmetrize(wt[, , t])
		}
		cand <- .ae_bcd_fit(
			zk,
			prep$X,
			R,
			symmetric,
			bcd_iter,
			inner_tol,
			init = list(
				mu = fit$mu,
				beta = fit$beta,
				a = fit$a,
				b = fit$b,
				U = fit$U,
				V = fit[["V"]],
				L = fit$L
			),
			lowrank_method = lowrank_method,
			weights = wt,
			linear_solver = linear_solver,
			uv_prior = if (is.null(uv_prior)) NULL else
				c(
					uv_prior,
					list(
						sig2_u = fit$sig2_u,
						sig2_v = fit$sig2_v
					)
				),
			ab_prior = eb_state(fit, eta_fit, terms_w)
		)
		if (has_ve) cand$ve <- ve_new
		eta_cand <- eta_of(cand)
		terms_cand <- terms_of(eta_cand, cand)
		dev_cand <- dev_fun(eta_cand, terms_cand)
		# the candidate's Fisher weights enter its objective
		work_cand <- work_fun(eta_cand, terms_cand)
		dev_k <- pdev_of(cand, dev = dev_cand, work = work_cand)
		prev <- pdev_of(fit, cand, dev = dev_fit, work = work_fit)
		# the reweighted equations carry the previous Laplace correction, so a
		# step is accepted when it lowers the scaled residual of the equations
		cert_cand <- cert_of(cand, eta_cand, work_cand, terms_cand)
		merit_cand <- merit_of(cert_cand)
		improves <- function(dev_s, merit_s) {
			# allow the residual to rise by 5 percent, never beyond twice the best
			is.finite(dev_s) && merit_s < min(1.05 * merit_fit, 2 * merit_best)
		}
		anchor <- fit
		# a full step can overshoot: halve toward the current iterate and keep
		# the first improving point; stop when no halved step improves
		if (!improves(dev_k, merit_cand)) {
			for (step in 2^(-seq_len(20L))) {
				damped <- .ae_irls_blend(anchor, cand, step, symmetric)
				eta_damped <- eta_of(damped)
				terms_damped <- terms_of(eta_damped, damped)
				dev_damped <- dev_fun(eta_damped, terms_damped)
				work_damped <- work_fun(eta_damped, terms_damped)
				dev_s <- pdev_of(damped, dev = dev_damped, work = work_damped)
				cert_s <- cert_of(damped, eta_damped, work_damped, terms_damped)
				merit_s <- merit_of(cert_s)
				if (improves(dev_s, merit_s)) {
					cand <- damped; dev_k <- dev_s
					eta_cand <- eta_damped; dev_cand <- dev_damped
					terms_cand <- terms_damped; work_cand <- work_damped
					cert_cand <- cert_s; merit_cand <- merit_s
					break
				}
			}
		}
		dev_trace <- c(dev_trace, dev_k)
		if (!improves(dev_k, merit_cand)) {
			termination <- "no improving conditional step"
			break
		}
		conditional_gain <- c(conditional_gain, prev - dev_k)
		best_dev <- dev_k
		best_fit <- cand
		n_accepted <- n_accepted + 1L
		fit <- cand
		eta_fit <- eta_cand
		terms_fit <- terms_cand
		dev_fit <- dev_cand
		work_fit <- work_cand
		certificate_now <- cert_cand
		merit_fit <- merit_cand
		merit_best <- min(merit_best, merit_cand)
		if (isTRUE(certificate_now$pass)) {
			termination <- "converged"
			break
		}
		# the noise variance converges slowly: jump along its geometric tail
		if (has_ve && is.null(ve_fixed)) {
			ve_hist <- c(ve_hist, log(fit$ve))
			jump <- .ae_pln_aitken(ve_hist)
			if (is.finite(jump)) {
				fit$ve <- best_fit$ve <- exp(log(fit$ve) + jump)
				terms_fit <- terms_of(eta_fit, fit)
				dev_fit <- dev_fun(eta_fit, terms_fit)
				work_fit <- work_fun(eta_fit, terms_fit)
				certificate_now <- cert_of(fit, eta_fit, work_fit, terms_fit)
				merit_fit <- merit_best <- merit_of(certificate_now)
				ve_hist <- numeric(0)
			}
		}
		# a two-sweep inner plateau can hide a remaining score: tighten the inner
		# tolerance
		if (cand$iterations <= 2L && inner_tol > 64 * .Machine$double.eps) {
			inner_tol <- max(64 * .Machine$double.eps, inner_tol / 100)
			inner_refinements <- inner_refinements + 1L
		}
	}

	# the last certificate evaluated is that of the best (current) iterate
	certificate <- if (!is.null(certificate_now)) certificate_now else
		.ae_irls_stationarity(
			prep,
			best_fit,
			family,
			link,
			R,
			symmetric,
			tol,
			eta = eta_fit,
			xmat = xmat,
			working = work_fit,
			terms = terms_fit,
			ab_prior = ab_prior
		)
	certificate$outer_iterations <- n_accepted
	certificate$last_inner_iterations <- best_fit$iterations
	certificate$termination <- if (isTRUE(certificate$pass)) "converged" else termination
	certificate$conditional_gain <- conditional_gain
	certificate$inner_tolerance <- inner_tol
	certificate$inner_refinements <- inner_refinements
	certificate$score <- certificate$gradient <- certificate$information <- NULL
	best_fit$converged <- isTRUE(certificate$pass)
	if (warm && warm_retry && !isTRUE(certificate$pass)) {
		# a warm half-fit can sit in an inferior factor basin: retry from the
		# spectral start
		retry <- .ae_irls_run(
			prep,
			family,
			link,
			R,
			symmetric,
			max_iter,
			tol,
			lowrank_method = lowrank_method,
			irls_maxit = irls_maxit,
			irls_tol = irls_tol,
			linear_solver = linear_solver,
			ab_prior = ab_prior,
			ve_fixed = ve_fixed,
			vcov_mode = vcov_mode
		)
		if (isTRUE(retry$stationarity$pass)) {
			retry$stationarity$warm_retry <- list(
				attempted = TRUE,
				start = "cold",
				original = certificate
			)
			return(retry)
		}
	}

	# working response and weights at the final predictor
	eta <- eta_fit
	if (family == "binary") {
		# a strict separating direction in the unpenalized linear part means no
		# finite optimum exists
		linear_eta <- eta - array(rep(best_fit$Omat, Tt), dim(y))
		observed <- is.finite(y)
		if (any(observed) && all(is.finite(linear_eta[observed])) &&
				all((2 * y[observed] - 1) * linear_eta[observed] > 0)) {
			cli::cli_abort(
				c(
					"Binary ALS has no finite estimate: the unpenalized linear predictor completely separates the observed outcomes.",
					"i" = "Use {.fn ame} or {.fn lame} (MCMC) with proper coefficient and actor-effect priors."
				),
				class = "lame_als_separation"
			)
		}
	}
	working <- work_fit
	working$z[unobserved] <- NA_real_
	working$w[unobserved] <- NA_real_
	if (symmetric) for (t in seq_len(Tt)) {
		working$z[, , t] <- .ae_symmetrize(working$z[, , t])
		working$w[, , t] <- .ae_symmetrize(working$w[, , t])
	}
	best_fit$sse <- sum(working$w * (working$z - eta)^2, na.rm = TRUE)
	# cell scores at the mode summed over periods: the matrix the product of a count fit is fitted to
	# (.ae_pln_bulk_charge)
	if (identical(family, "poisson")) {
		best_fit$score_sum <- rowSums(working$w * (working$z - eta), dims = 2L, na.rm = TRUE)
	}
	# Laplace covariance of the regression coefficients at the final weights
	best_fit$ab_prior <- ab_prior
	best_fit$eb_vcov <- if (identical(vcov_mode, "none")) list() else
		tryCatch(
			.ae_eb_vcov(
				prep,
				best_fit,
				working$w,
				ab_prior,
				R,
				eta = eta,
				link = link,
				evidence_only = identical(vcov_mode, "evidence"),
				report = identical(vcov_mode, "full")
			),
			error = function(e) structure(list(), message = conditionMessage(e))
		)
	list(
		fit = best_fit,
		Z = working$z,
		W = working$w,
		dev_trace = dev_trace,
		irls_iters = n_accepted,
		best_obj = best_dev,
		data_deviance = dev_fit,
		eta = eta,
		stationarity = certificate
	)
}

# mark the effects of dropped boundary actors as not estimated: a / b and the
# predictor and fitted values of their cells are NA, and fit$not_estimated
# records the actors
.ae_mark_not_estimated <- function(out, db) {
	if (is.null(db) || (!length(db$row) && !length(db$col))) return(out)
	sym <- isTRUE(out$symmetric)
	if (length(db$row)) out$a[db$row] <- NA_real_
	if (length(db$col)) out$b[db$col] <- NA_real_
	blank <- function(M) {
		if (is.null(M)) return(M)
		if (length(db$row)) M[db$row, ] <- NA_real_
		if (length(db$col)) M[, db$col] <- NA_real_
		M
	}
	for (nm in c("EZ", "fitted", "YPM")) {
		if (is.list(out[[nm]])) out[[nm]] <- lapply(out[[nm]], blank)
	}
	lab <- function(idx, val, names, role) {
		if (!length(idx)) return(NULL)
		data.frame(
			actor = names[idx],
			side = role,
			outcome = ifelse(
				is.na(val),
				"no observed cell",
				ifelse(val > 0, "only ties", "no tie")
			),
			stringsAsFactors = FALSE
		)
	}
	role_r <- if (sym) "actor" else if (identical(out$mode, "bipartite")) "row" else "sender"
	role_c <- if (identical(out$mode, "bipartite")) "column" else "receiver"
	out$not_estimated <- rbind(
		lab(db$row, db$row_val, out$row_names, role_r),
		if (!sym) lab(db$col, db$col_val, out$col_names, role_c)
	)
	out
}

# the message naming the actors whose effect is not estimated
.ae_inform_not_estimated <- function(ne) {
	if (is.null(ne) || !nrow(ne)) return(invisible(NULL))
	few <- function(x) paste(utils::head(x, 5L), collapse = ", ")
	txt <- paste0(ne$actor, " (", ne$side, ", ", ne$outcome, ")")
	boundary <- any(ne$outcome != "no observed cell")
	empty <- any(ne$outcome == "no observed cell")
	if (!empty) {
		cli::cli_inform(c(
			"i" = "{nrow(ne)} actor effect{?s} not estimated: {few(txt)}{if (nrow(ne) > 5L) ', ...' else ''}.",
			" " = "An actor with a single outcome value on a side has its effect there at the boundary of the probit / logit fit; its cells on that side are left out, its effect is NA in {.code a} / {.code b}, and {.code va} / {.code vb} describe the estimated effects. The list is in {.code fit$not_estimated}."
		))
		return(invisible(NULL))
	}
	# an actor whose covariates are missing in every period has no observed cell
	cli::cli_inform(c(
		"i" = "{nrow(ne)} actor effect{?s} not estimated: {few(txt)}{if (nrow(ne) > 5L) ', ...' else ''}.",
		" " = if (boundary) "An actor with a single outcome value on a side has its effect there at the boundary of the probit / logit fit, and its cells on that side are left out.",
		" " = "An actor with no observed cell on a side (a covariate missing in every period leaves none) has no effect there to estimate. The effect is NA in {.code a} / {.code b}, {.code va} / {.code vb} describe the estimated effects, and the list is in {.code fit$not_estimated}."
	))
	invisible(NULL)
}

# an effect vector with its not-estimated entries at zero (for a warm start)
.ae_na_zero <- function(x) {
	if (is.null(x)) return(x)
	x[is.na(x)] <- 0
	x
}

# dropped-actor record of a refit: actors this data puts at the boundary plus
# those the original fit left out
.ae_refit_boundary <- function(db, ne, row_names, col_names, symmetric) {
	if (is.null(ne) || !nrow(ne)) return(db)
	row_side <- ne$side %in% c("sender", "row", "actor")
	r <- match(ne$actor[row_side], row_names)
	cc <- if (symmetric) r else match(ne$actor[!row_side], col_names)
	val_of <- function(o) {
		ifelse(o == "only ties", 1, ifelse(o == "no tie", 0, NA_real_))
	}
	rv <- val_of(ne$outcome[row_side])
	cv <- if (symmetric) rv else val_of(ne$outcome[!row_side])
	keep_r <- !is.na(r); keep_c <- !is.na(cc)
	old <- list(
		row = r[keep_r],
		col = cc[keep_c],
		row_val = rv[keep_r],
		col_val = cv[keep_c]
	)
	if (is.null(db)) return(old)
	add_r <- !old$row %in% db$row
	add_c <- !old$col %in% db$col
	list(
		row = c(db$row, old$row[add_r]),
		col = c(db$col, old$col[add_c]),
		row_val = c(db$row_val, old$row_val[add_r]),
		col_val = c(db$col_val, old$col_val[add_c])
	)
}
