# ame_als_irls.R
#
# iteratively reweighted least squares (IRLS) for the non-normal families of
# the fast AME estimator. With `non_normal_method = "irls"` the estimator
# repeatedly rebuilds a family-specific Gaussian working response and weight
# from the current linear predictor and refits the weighted BCD core, instead
# of fitting a single fixed transform. The result is a fast *approximate GLM*
# AME fit -- closer to the calibrated MCMC estimator than the one-shot
# `log(y+1)` / rank-normal transforms, at a few times their cost.
#
# supported: binary (logit / probit link) and poisson (log link). Ordinal stays
# on the rank-normal transform (cumulative-link IRLS is left to MCMC).

# --- family working-response / weight maps (elementwise, NA-propagating) -----

# each map clips the link inputs (p, phi) and the working-response increment so
# separated / sparse data cannot produce an infinite or NaN working response;
# the increment cap is wide enough never to bite in ordinary data.

# poisson, log link:  z = eta + (y - mu)/mu ,  w = mu = exp(eta)
.ae_irls_poisson <- function(y, eta, eta_cap = 20, wfloor = 1e-6,
                             delta_cap = 30) {
	mu <- exp(pmin(pmax(eta, -eta_cap), eta_cap))
	delta <- pmin(pmax((y - mu) / pmax(mu, wfloor), -delta_cap), delta_cap)
	list(z = eta + delta, w = pmax(mu, wfloor))
}

# binary, logit link:  z = eta + (y - p)/(p(1-p)) ,  w = p(1-p)
.ae_irls_binary_logit <- function(y, eta, eta_cap = 30, pfloor = 1e-6,
                                  wfloor = 1e-8, delta_cap = 30) {
	p <- stats::plogis(pmin(pmax(eta, -eta_cap), eta_cap))
	p <- pmin(pmax(p, pfloor), 1 - pfloor)
	v <- p * (1 - p)
	delta <- pmin(pmax((y - p) / v, -delta_cap), delta_cap)
	list(z = eta + delta, w = pmax(v, wfloor))
}

# binary, probit link:  z = eta + (y - p)/phi ,  w = phi^2 / (p(1-p))
.ae_irls_binary_probit <- function(y, eta, eta_cap = 8, pfloor = 1e-6,
                                   phifloor = 1e-8, wfloor = 1e-8,
                                   delta_cap = 30) {
	ec  <- pmin(pmax(eta, -eta_cap), eta_cap)
	p   <- pmin(pmax(stats::pnorm(ec), pfloor), 1 - pfloor)
	phi <- pmax(stats::dnorm(ec), phifloor)
	delta <- pmin(pmax((y - p) / phi, -delta_cap), delta_cap)
	list(z = eta + delta, w = pmax(phi^2 / (p * (1 - p)), wfloor))
}

# --- family deviances (the IRLS ranking objective) --------------------------

.ae_dev_poisson <- function(y, mu) {
	mu <- pmax(mu, 1e-12)
	term <- ifelse(y == 0 | !is.finite(y), 0, y * log(y / mu))
	2 * sum((term - (y - mu))[is.finite(y)])
}

.ae_dev_binary <- function(y, p) {
	p  <- pmin(pmax(p, 1e-12), 1 - 1e-12)
	ok <- is.finite(y)
	-2 * sum((y * log(p) + (1 - y) * log(1 - p))[ok])
}

# --------------------------------------------------------------------------
# IRLS outer loop
# --------------------------------------------------------------------------

# prep      : canonical prep list (from .ae_prepare / .ae_prep_from_canonical)
# init      : optional warm-start state (mu, beta, a, b, U, V, L); when NULL the
#             loop is seeded by the one-shot transform fit
# returns   : list(fit = best BCD fit, Z = its working response,
#                  dev_trace, irls_iters)
.ae_irls_run <- function(prep, family, link, R, symmetric, max_iter, tol,
                         lowrank_method = "mm", irls_maxit = 5L,
                         irls_tol = 1e-4, init = NULL,
                         linear_solver = "eigen") {

	Tt <- prep$Tt; nr <- prep$nr; nc <- prep$nc
	y  <- prep$Yarr

	# linear predictor array from a BCD fit
	eta_of <- function(f) {
		ab <- outer(f$a, f$b, "+")
		E  <- array(0, dim = c(nr, nc, Tt))
		for (t in seq_len(Tt)) {
			xb <- matrix(0, nr, nc)
			if (prep$p_dyad > 0) for (k in seq_len(prep$p_dyad)) {
				xb <- xb + f$beta[k] * prep$X[, , k, t]
			}
			E[, , t] <- f$mu + xb + ab + f$Omat
		}
		E
	}

	work_fun <- switch(paste(family, link),
		"poisson log"   = function(eta) .ae_irls_poisson(y, eta),
		"binary logit"  = function(eta) .ae_irls_binary_logit(y, eta),
		"binary probit" = function(eta) .ae_irls_binary_probit(y, eta),
		cli::cli_abort("IRLS is not available for family {.val {family}}."))
	dev_fun <- switch(family,
		poisson = function(eta) .ae_dev_poisson(y,
			exp(pmin(pmax(eta, -20), 20))),
		binary  = function(eta) .ae_dev_binary(y,
			if (link == "logit") stats::plogis(eta) else stats::pnorm(eta)))

	# --- seed: warm start, or the one-shot transform fit ---
	Z0 <- .ae_working_response(y, family)
	if (symmetric) for (t in seq_len(Tt)) Z0[, , t] <- .ae_symmetrize(Z0[, , t])
	fit <- .ae_bcd_fit(Z0, prep$X, R, symmetric, max_iter, tol,
	                   init = init, lowrank_method = lowrank_method,
	                   linear_solver = linear_solver)

	best_fit <- fit
	best_dev <- dev_fun(eta_of(fit))
	best_Z   <- Z0
	best_W   <- NULL                        # the seed fit is unit-weighted
	dev_trace <- best_dev

	for (k in seq_len(irls_maxit)) {
		eta <- eta_of(fit)
		wk  <- work_fun(eta)
		zk  <- wk$z; wt <- wk$w
		zk[!is.finite(y)] <- NA_real_           # keep the observed pattern
		wt[!is.finite(y)] <- NA_real_
		if (symmetric) for (t in seq_len(Tt)) {
			zk[, , t] <- .ae_symmetrize(zk[, , t])
			wt[, , t] <- .ae_symmetrize(wt[, , t])
		}
		cand <- .ae_bcd_fit(zk, prep$X, R, symmetric, max_iter, tol,
		                    init = list(mu = fit$mu, beta = fit$beta,
		                                a = fit$a, b = fit$b,
		                                U = fit$U, V = fit$V, L = fit$L),
		                    lowrank_method = lowrank_method, weights = wt,
		                    linear_solver = linear_solver)
		dev_k <- dev_fun(eta_of(cand))
		prev  <- dev_trace[length(dev_trace)]
		dev_trace <- c(dev_trace, dev_k)
		# accept only an improving iterate; on an overshoot stop and keep the
		# best so far (rather than carrying the overshot state forward)
		if (!is.finite(dev_k) || dev_k > best_dev) break
		best_dev <- dev_k; best_fit <- cand; best_Z <- zk; best_W <- wt
		fit <- cand
		if (is.finite(prev) &&
		    abs(prev - dev_k) / max(1, abs(prev)) < irls_tol) break
	}

	list(fit = best_fit, Z = best_Z, W = best_W, dev_trace = dev_trace,
	     irls_iters = length(dev_trace) - 1L)
}
