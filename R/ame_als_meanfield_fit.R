# the variational refinement stage of the empirical-Bayes fit: it turns a converged
# mode (.ae_irls_run) into the reported fit, with the predictive moments, evidence and
# coefficient covariance of such a fit (the engine is in ame_als_meanfield.R).
# fit$eb$log_marginal, which compares ranks, is the Laplace evidence of the starting
# mode, not the lower bound, whose slack grows with the rank

# refine the converged mode `ir` (from .ae_irls_run) to the variational posterior of
# its effects and factors; with `start` (a state from .ae_mf_object_state) the
# iteration starts there, as in a refit on new data
#
# with `defer_uncertainty` the scale part of the coefficient covariance is left to
# .ae_mf_finish (`pending`), so fits that are never reported skip it
.ae_mf_refine <- function(
	prep,
	ir,
	ab,
	R,
	symmetric,
	link,
	max_iter,
	tol,
	verbose = FALSE,
	start = NULL,
	pair = NULL,
	ve_fixed = FALSE,
	defer_uncertainty = FALSE
) {
	scale <- .ae_mf_scale_rule(R)
	fit_with <- function(stat) .ae_mf_fit(
		prep,
		ir$fit,
		ab,
		R,
		symmetric,
		link,
		ir$W,
		max_iter = .ae_mf_budget(max_iter, scale),
		tol = .ae_mf_tol(tol),
		verbose = verbose,
		start = start,
		scale = scale,
		pair = pair,
		ve_fixed = ve_fixed,
		scale_stat = stat
	)
	res <- fit_with(NULL)
	# the posterior-mean scale update is not a stationary point of the objective the damped
	# steps compare, so when no step improves the fit is repeated with the plug-in update
	if (!isTRUE(res$converged) && identical(res$termination, "no improving step") && R > 0L && !is.null(res$scale) &&
			identical(res$scale$stat, "mean")) {
		res <- fit_with("plugin")
	}
	cx <- res$cx
	S <- res$S
	# the fit reports the cavity covariances (.ae_mf_cavity_cov); the state of the updates,
	# its weights and its lower bound are kept so a refit and the bound start there
	Sr <- S
	evr <- res$ev
	if (!is.null(res$scale)) {
		rest <- .ae_mf_cavity_cov(cx, S, res$scale$C)
		Sr$Sig <- rest$Sig; if (!cx$sym) Sr$Om <- rest$Om
		# the predictor at the cavity covariances (with the dyad covariance under a dyadic correlation)
		pr <- .ae_mf_predictor(cx, Sr)
		evr <- list(eta = pr$eta, V = pr$Vt, ve = S$ve)
		if (!is.null(cx$pair)) evr$pair <- list(cv = .ae_pair_cov(cx, Sr), rho = S$rho)
	}
	f <- ir$fit
	f$mu <- S$mu; f$beta <- S$beta; f$a <- S$a; f$b <- S$b; f$gamma <- S$gamma
	f$sab <- S$sab
	f$sab_change <- 0
	if (cx$has_ve) f$ve <- S$ve
	if (R > 0L) {
		f$Omat <- res$ev$Pm
		if (cx$sym) {
			rep <- .ae_mf_sym_rep(S$F, S$sg)
			f$U <- rep$U; f[["V"]] <- rep$U; f$L <- rep$L
			f$sig2_u <- f$sig2_v <- rep_len(S$prior$tau, R)
		} else {
			f$U <- S$U; f[["V"]] <- S[["V"]]; f$L <- NULL
			if (identical(S$prior$kind, "joint")) {
				f$sig2_u <- diag(S$prior$Suv)[seq_len(R)]
				f$sig2_v <- diag(S$prior$Suv)[R + seq_len(R)]
			} else {
				f$sig2_u <- S$prior$su; f$sig2_v <- S$prior$sv
			}
		}
	}
	W <- res$ev$wk$w; W[!cx$obs] <- NA_real_
	Z <- res$ev$wk$z; Z[!cx$obs] <- NA_real_
	eta <- res$ev$eta
	elbo <- res$elbo
	# converged when both the variational iteration and its starting mode are; a refit
	# has no mode of its own
	mode_converged <- if (is.null(start)) isTRUE(ir$fit$converged) else TRUE
	# events per loading: below about three the product is weakly identified and the
	# rank is not determined by the data
	loadings <- R * (if (cx$sym) cx$nr else cx$nr + cx$nc)
	events <- cx$nu * (if (cx$pln) sum(cx$y > 0, na.rm = TRUE) else sum(cx$y == 1, na.rm = TRUE))
	record <- list(
		mode_converged = mode_converged,
		information = if (R > 0L && !cx$gauss) .ae_glm_information(cx, res$ev$Pm, W, loadings, events) else
		if (R > 0L) .ae_gauss_information(cx, res$ev$Pm, S$ve, loadings),
		converged = res$converged,
		termination = res$termination,
		iterations = res$iterations,
		residual = res$residual,
		gain = res$gain,
		step = res$step,
		tolerance = res$tol,
		objective = res$merit,
		elbo = elbo$elbo,
		ell = elbo$ell,
		kl_actors = elbo$kl_actors,
		kl_factors = elbo$kl_factors,
		scale = if (!is.null(res$scale)) list(
			rule = res$scale$rule,
			s = res$scale$s,
			rho = res$scale$rho,
			sd_log = NA_real_
		) else list(rule = "plugin"),
		start = if (is.null(start)) "posterior mode" else "previous variational fit",
		mode_log_marginal = if (!is.null(ir$fit$eb_vcov$evidence)) ir$fit$eb_vcov$evidence$log_marginal else NA_real_
	)
	# the part of the uncertainty of the coefficients that comes from the scale of the factor prior
	extra <- NULL
	pending <- NULL
	if (!is.null(res$scale) && identical(res$scale$rule, "cavity") && isTRUE(defer_uncertainty)) {
		pending <- list(prep = prep, ab = ab, R = R, symmetric = cx$sym, link = link, res = res)
	} else if (!is.null(res$scale) && identical(res$scale$rule, "cavity") && .ae_mf_uncertainty_on()) {
		su <- tryCatch(.ae_mf_scale_uncertainty(prep, ab, R, cx$sym, link, res), error = function(e) NULL)
		if (!is.null(su) && is.finite(su$sd_log)) {
			extra <- list(J = su$J, var = su$sd_log^2)
			record$scale$sd_log <- su$sd_log
		}
	}
	if (!is.null(res$pair)) record$pair <- res$pair
	posterior <- list(
		kind = cx$ab$kind,
		Sig = Sr$Sig,
		Om = Sr$Om,
		F = S$F,
		sg = S$sg,
		prior = S$prior,
		Sig_mf = if (!is.null(res$scale)) S$Sig,
		Om_mf = if (!is.null(res$scale)) S$Om
	)
	# expected cross weights of the pairs, for the joint information of the predictor
	f$eb_mf <- list(
		record = record,
		posterior = posterior,
		cross = if (!is.null(cx$pair)) res$ev$pair$cw
	)
	f$converged <- isTRUE(res$converged) && mode_converged
	# coefficient covariance under the variational posterior; the rank evidence stays
	# the Laplace value of the starting mode
	mode_evidence <- ir$fit$eb_vcov$evidence
	# a count network with a dyadic correlation (ame_als_pln_pair.R): the pair correlation
	# and information at the fitted state change the standard errors, not the effects
	pair_fit <- NULL
	if (cx$pln && !is.null(pair)) {
		pair_fit <- tryCatch(.ae_pln_pair_stage(cx, Sr, evr, S$ve), error = function(e) {
			cli::cli_warn(
				c(
					"The dyadic correlation of the counts could not be computed; the cells of a dyad are treated as independent.",
					"x" = conditionMessage(e)
				),
				class = "lame_als_pair"
			)
			NULL
		})
		if (!is.null(pair_fit)) {
			record$pair <- list(
				rho = pair_fit$rho,
				score = 0,
				curvature = pair_fit$prec,
				sd = pair_fit$sd,
				n_pairs = pair_fit$n_pairs,
				quad = pair_fit$quad,
				method = "integrated pair likelihood"
			)
			f$eb_mf$record <- record
			W[pair_fit$one] <- pair_fit$w[pair_fit$one]; W[pair_fit$two] <- pair_fit$w[pair_fit$two]
		}
	}
	# under the cavity rule the coefficient information uses the cavity weight; the
	# plug-in rule and a fit without factors keep the expected weight
	cross <- NULL
	# a Gaussian cell's information does not depend on the predictor, so the expected
	# weight 1 / ve is exact there
	Winfo <- if (!is.null(res$scale) && !cx$gauss) {
		wi <- if (!is.null(cx$pair)) {
			ci <- .ae_mf_info_pair(cx, evr)
			cross <- ci$cw
			ci$w
		} else .ae_mf_info_w(cx$y, evr$eta, array(evr$V, dim(cx$y)), link, ve = S$ve)
		wi[!cx$obs] <- NA_real_
		wi
	} else {
		if (!is.null(cx$pair)) cross <- res$ev$pair$cw
		W
	}
	if (!is.null(pair_fit)) {
		Winfo[pair_fit$one] <- pair_fit$w[pair_fit$one]; Winfo[pair_fit$two] <- pair_fit$w[pair_fit$two]
		cross <- pair_fit$cw
		f$eb_mf$cross <- cross
	}
	vb_cov <- list(S = Sr, elbo = elbo, evidence = mode_evidence, cross = cross)
	f$eb_vcov <- tryCatch(
		.ae_eb_vcov(prep, f, Winfo, ab, R, eta = eta, link = link, vb = vb_cov, extra = extra),
		error = function(e) structure(list(), message = conditionMessage(e))
	)
	if (!is.null(pending)) pending$vcov <- list(W = Winfo, eta = eta, vb = vb_cov)
	list(
		fit = f,
		W = W,
		Z = Z,
		eta = eta,
		record = record,
		converged = res$converged,
		objective = -2 * elbo$elbo / ab$nu,
		pending = pending
	)
}

# add the scale uncertainty a refinement with `defer_uncertainty` left out to the
# reported fit `ir` (the same computation and state as in .ae_mf_refine)
.ae_mf_finish <- function(ir) {
	pd <- ir$mf_pending
	ir$mf_pending <- NULL
	if (is.null(pd) || !.ae_mf_uncertainty_on()) return(ir)
	su <- tryCatch(
		.ae_mf_scale_uncertainty(pd$prep, pd$ab, pd$R, pd$symmetric, pd$link, pd$res),
		error = function(e) NULL
	)
	if (is.null(su) || !is.finite(su$sd_log)) return(ir)
	f <- ir$fit
	f$eb_mf$record$scale$sd_log <- su$sd_log
	# a covariance that could not be formed stays as it is
	if (!is.null(f$eb_vcov$vcov)) {
		vc <- tryCatch(
			.ae_eb_vcov(
				pd$prep,
				f,
				pd$vcov$W,
				pd$ab,
				pd$R,
				eta = pd$vcov$eta,
				link = pd$link,
				vb = pd$vcov$vb,
				extra = list(J = su$J, var = su$sd_log^2)
			),
			error = function(e) NULL
		)
		if (!is.null(vc)) f$eb_vcov$vcov <- vc$vcov
	}
	ir$fit <- f
	ir
}

# uncertainty of the coefficients from the scale of the factor prior: the log scale
# has precision c = c_data + c_prior, so the covariance gains J J' / c with
# J = d coefficients / d log scale, both from two fits at the scale times exp(+-h)
.ae_mf_scale_uncertainty <- function(prep, ab, R, symmetric, link, res, h = 0.15, tol = 1e-4, max_iter = 200L) {
	S <- res$S
	one <- function(sgn) {
		S1 <- S
		S1$prior <- .ae_mf_prior_shift(S$prior, exp(sgn * h))
		r1 <- .ae_mf_fit(
			prep,
			NULL,
			ab,
			R,
			symmetric,
			link,
			NULL,
			max_iter = max_iter,
			tol = tol,
			start = S1,
			scale = "fixed",
			ve_fixed = isTRUE(res$cx$ve_fixed),
			pair = if (!is.null(res$cx$pair)) list(rho = S$rho, quad = res$cx$pair_quad, update = FALSE)
		)
		C1 <- .ae_mf_info_sum(r1$cx, r1$ev)
		list(theta = c(r1$S$mu, r1$S$beta, r1$S$gamma), g = .ae_mf_scale_score(r1$cx, r1$S, C1), ok = r1$converged)
	}
	a <- one(1); b <- one(-1)
	c_data <- -(a$g - b$g) / (2 * h)
	c_prior <- .ae_mf_scale_prior_curvature(res$cx, S)
	cc <- max(c_data, 0) + c_prior
	J <- (a$theta - b$theta) / (2 * h)
	list(
		J = J,
		c = cc,
		c_data = c_data,
		c_prior = c_prior,
		sd_log = if (is.finite(cc) && cc > 0) 1 / sqrt(cc) else NA_real_,
		ok = a$ok && b$ok
	)
}

# the refit of a variational fit on new data: the iteration starts at the fit's
# own variational state (the previous fit is the warm start, so no mode is
# computed), on the same actors and the same rank
.ae_mf_refit <- function(
	prep,
	object,
	ab,
	link,
	max_iter,
	tol,
	verbose = FALSE,
	mode_args = NULL,
	ve_hold = FALSE
) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	ms <- .ae_mf_object_state(object)
	pln <- identical(object$family, "poisson")
	mode_evidence <- NULL
	steps <- 0L
	if (pln) {
		# a count network: the noise variance of the new data maximises the Laplace marginal
		# likelihood of its mode, by a local search about the fit's value (.ae_pln_ve_search)
		lrm <- .ae_default(object$lowrank_method, "mm")
		lsv <- .ae_default(object$linear_solver, "eigen")
		run <- function(ve, from) {
			.ae_irls_run(
				prep,
				object$family,
				link,
				object$R,
				isTRUE(object$symmetric),
				max_iter,
				tol,
				lowrank_method = lrm,
				init = .ae_pln_mode_init(from$fit, ve),
				linear_solver = lsv,
				warm = TRUE,
				ab_prior = ab,
				ve_fixed = ve,
				warm_retry = FALSE
			)
		}
		from0 <- list(
			fit = list(
				mu = ms$S$mu,
				beta = ms$S$beta,
				a = ms$S$a,
				b = ms$S$b,
				gamma = ms$S$gamma,
				sab = ms$S$sab,
				U = object$U,
				V = object[["V"]],
				L = object$L
			)
		)
		if (isTRUE(ve_hold)) {
			# a converged fit's own data: hold its noise variance so the refit does not move the
			# estimates; the mode at that value gives the evidence
			r0 <- run(ms$S$ve, from0)
			mode_evidence <- r0$fit$eb_vcov$evidence
			steps <- as.integer(r0$irls_iters)
		} else {
			search <- .ae_pln_ve_search(run, ms$S$ve, from0 = from0, local = TRUE, ve_prior = ab$ve_prior)
			ms$S$ve <- search$ve
			mode_evidence <- search$ir$fit$eb_vcov$evidence
			steps <- as.integer(sum(search$curve$iterations))
		}
	} else if (!is.null(mode_args) && !identical(object$family, "normal")) {
		# the marginal likelihood of the new data (Laplace at its mode) when the caller wants it
		mode_evidence <- tryCatch(
			.ae_irls_run(
				prep,
				object$family,
				link,
				object$R,
				isTRUE(object$symmetric),
				200L,
				1e-6,
				lowrank_method = mode_args$lowrank_method,
				linear_solver = mode_args$linear_solver,
				ab_prior = ab
			)$fit$eb_vcov$evidence,
			error = function(e) NULL
		)
	}
	fit0 <- list(
		mu = ms$S$mu,
		beta = ms$S$beta,
		a = ms$S$a,
		b = ms$S$b,
		gamma = ms$S$gamma,
		sab = ms$S$sab,
		ab_prior = ab,
		Omat = matrix(0, prep$nr, prep$nc),
		U = object$U,
		V = object[["V"]],
		L = object$L,
		sig2_u = object$eb$sig2_u,
		sig2_v = object$eb$sig2_v,
		iterations = object$iterations,
		converged = TRUE,
		sse = NA_real_,
		dev_history = NULL,
		eb_vcov = list(evidence = mode_evidence)
	)
	ir <- list(fit = fit0, W = object$obs_weights)
	# the refits run at least 50 iterations: the variational iteration may need more steps on new data
	vb <- .ae_mf_refine(
		prep,
		ir,
		ab,
		object$R,
		isTRUE(object$symmetric),
		link,
		max(max_iter, 50L),
		tol,
		verbose,
		start = ms$S,
		pair = if (!is.null(object$eb$pair)) list(rho = object$eb$pair$rho, quad = object$eb$pair$quad),
		ve_fixed = pln
	)
	# the reweighting steps this refit ran (none for a binary fit, which computes no mode)
	vb$irls_iters <- steps
	# a normal fit: the evidence of the new data is that of the mode at the variance of the noise the refit settled on
	if (!is.null(mode_args) && identical(object$family, "normal")) {
		ir2 <- .ae_gauss_evidence(
			prep,
			list(fit = vb$fit),
			ab,
			object$R,
			isTRUE(object$symmetric),
			link,
			max(max_iter, 50L),
			tol,
			mode_args$lowrank_method,
			mode_args$linear_solver
		)
		vb$fit <- ir2$fit
	}
	vb
}

# the variational stage of one converged mode `ir`, in the form the multistart
# selection and the assembly read; `defer_uncertainty` leaves the scale
# uncertainty for .ae_mf_finish (`mf_pending`)
.ae_mf_stage <- function(
	prep,
	ir,
	ab,
	R,
	symmetric,
	link,
	max_iter,
	tol,
	pair = NULL,
	ve_fixed = FALSE,
	defer_uncertainty = FALSE
) {
	vb <- .ae_mf_refine(
		prep,
		ir,
		ab,
		R,
		symmetric,
		link,
		max_iter,
		tol,
		pair = pair,
		ve_fixed = ve_fixed,
		defer_uncertainty = defer_uncertainty
	)
	ir$fit <- vb$fit
	ir$W <- vb$W; ir$Z <- vb$Z; ir$eta <- vb$eta
	ir$best_obj <- vb$objective
	ir$selection_converged <- vb$converged && isTRUE(ir$stationarity$pass)
	ir$mf_pending <- vb$pending
	ir
}

# the iteration budget and tolerance of the refinement for the caller's
# `max_iter` and `tol`: the variational iteration is a damped fixed point and
# needs more iterations than a reweighting step
.ae_mf_budget <- function(max_iter, scale = "plugin") {
	(if (identical(scale, "cavity")) 5L else 2L) * as.integer(max_iter)
}

# rule for the scale of the factor prior: "cavity" (default) updates it from the cavity
# second moments of the rows, "plugin" (options(lame.eb_scale)) from the mean-field
# ones, which understate them on weak data; a fit without factors has no scale
.ae_mf_scale_rule <- function(R) {
	opt <- getOption("lame.eb_scale", "cavity")
	if (!is.character(opt) || length(opt) != 1L || !opt %in% c("cavity", "plugin")) {
		cli::cli_abort("{.code options(lame.eb_scale)} must be {.val cavity} or {.val plugin}.")
	}
	if (R == 0L || identical(opt, "plugin")) "plugin" else "cavity"
}
# whether the standard errors carry the uncertainty of the scale (lame.eb_scale_uncertainty, default TRUE)
.ae_mf_uncertainty_on <- function() {
	opt <- getOption("lame.eb_scale_uncertainty", TRUE)
	if (!is.logical(opt) || length(opt) != 1L || is.na(opt)) {
		cli::cli_abort("{.code options(lame.eb_scale_uncertainty)} must be {.val TRUE} or {.val FALSE}.")
	}
	opt
}
.ae_mf_tol <- function(tol) max(tol, 1e-8)

# refuse a refinement that stopped before its fixed point (class lame_als_convergence),
# reporting its state and the remedy
.ae_mf_refuse <- function(record, R = NULL, requested = NULL) {
	why <- switch(
		record$termination,
		"iteration limit" = paste0("it reached the limit of ", record$iterations, " iterations"),
		"no improving step" = paste0("no damped step improved the objective after ", record$iterations, " iterations"),
		record$termination
	)
	rk <- if (is.null(R)) "this rank" else paste0("rank ", R)
	# a fit without factors: nearly collinear covariates are the usual cause
	if (!is.null(R) && R == 0L && (is.null(requested) || requested <= R)) {
		cli::cli_abort(
			c(
				"The variational refinement of the coefficients and effects did not converge: {why}.",
				"x" = "The largest step a full update would still take is {signif(record$residual, 2)} (tolerance {signif(record$tolerance, 2)}); the data do not determine the coefficients and effects.",
				"i" = "Nearly collinear covariates are the usual cause: drop or combine them. Otherwise increase {.arg max_iter}, or use {.fn ame} / {.fn lame} (MCMC)."
			),
			class = "lame_als_convergence"
		)
	}
	# a fit whose columns were retired down to the rank it is refused at (R/ame_als_retire.R)
	left <- if (!is.null(requested) && !is.null(R) && requested > R) paste0(
		" (rank ",
		requested,
		" was requested and ",
		requested - R,
		" column",
		if (requested - R > 1L) "s were" else " was",
		" retired)"
	) else ""
	cli::cli_abort(
		c(
			"The variational refinement of the latent factors did not converge: {why}.",
			"x" = "The largest step a full update would still take is {signif(record$residual, 2)} (tolerance {signif(record$tolerance, 2)}); the data do not determine the product at {rk}{left}.",
			"i" = "Fit a smaller rank, increase {.arg max_iter}, or use {.fn ame} / {.fn lame} (MCMC) at this rank."
		),
		class = "lame_als_convergence"
	)
}

# the eigen representation of an undirected factor: F S F' = U diag(L) U' with
# orthonormal columns U and signed eigenvalues L, by decreasing |L|
.ae_mf_sym_rep <- function(Fm, sg) {
	R <- ncol(Fm)
	sv <- svd(Fm, nu = R, nv = R)
	N <- crossprod(sv$v, sg * sv$v)
	N <- (sv$d * N) * rep(sv$d, each = R)
	eg <- eigen((N + t(N)) / 2, symmetric = TRUE)
	ord <- order(abs(eg$values), decreasing = TRUE)
	list(U = sv$u %*% eg$vectors[, ord, drop = FALSE], L = eg$values[ord])
}

# the variational state of a fitted object (the inverse of the packing in
# .ae_mf_refine): the means and covariances of the factors, the effects and
# their covariance, recomputed from the weights the fit stores
.ae_mf_object_state <- function(object) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	st <- .ae_eb_state(object)
	ab <- st$ab; fit <- st$fit; R <- st$R
	post <- object$eb$factor_posterior
	S <- list(
		mu = fit$mu,
		beta = as.numeric(fit$beta),
		a = as.numeric(fit$a),
		b = as.numeric(fit$b),
		gamma = as.numeric(fit$gamma),
		sab = fit$sab
	)
	if (object$family %in% c("poisson", "normal")) S$ve <- fit$ve
	if (!is.null(object$eb$pair)) S$rho <- object$eb$pair$rho
	W <- st$w; W[!is.finite(W)] <- 0
	info_s <- .ae_eb_info_s(rowSums(W, dims = 2L), ab)
	# the cross weights of the pairs of a fit with a dyadic correlation couple the effects of the two actors of a dyad
	if (!is.null(st$cross)) {
		Cx <- rowSums(st$cross, dims = 2L)
		nr <- ab$nr
		ia <- seq_len(nr); ib <- nr + seq_len(nr)
		info_s[ia, ia] <- info_s[ia, ia] + Cx
		info_s[ib, ib] <- info_s[ib, ib] + Cx
		info_s[cbind(ia, ib)] <- info_s[cbind(ia, ib)] + rowSums(Cx)
		info_s[cbind(ib, ia)] <- info_s[cbind(ib, ia)] + rowSums(Cx)
	}
	S$Sc <- .ae_eb_inverse(info_s + .ae_eb_precision(S$sab, ab))
	S$vact <- .ae_mf_vact(S$Sc, ab)
	if (R > 0L) {
		S$prior <- post$prior
		# the state of the updates: the mean-field covariances when the fit reports cavity ones
		S$Sig <- if (!is.null(post$Sig_mf)) post$Sig_mf else post$Sig
		if (identical(ab$kind, "symmetric")) {
			S$F <- post$F; S$sg <- post$sg
		} else {
			S$U <- fit$U; S$V <- fit[["V"]]; S$Om <- if (!is.null(post$Om_mf)) post$Om_mf else post$Om
		}
	}
	list(S = S, prep = st$prep, ab = ab, R = R, w = st$w)
}

# the evidence lower bound of a variational fit, recomputed from the object at
# its stored state (the value stored in fit$eb$variational$elbo)
.ae_mf_elbo_object <- function(object) {
	ms <- .ae_mf_object_state(object)
	link <- .ae_default(object$link, "probit")
	cx <- .ae_mf_context(
		ms$prep,
		ms$ab,
		ms$R,
		link,
		pair = if (!is.null(object$eb$pair)) list(rho = object$eb$pair$rho, quad = object$eb$pair$quad, update = FALSE)
	)
	ev <- .ae_mf_evaluate(cx, ms$S)
	.ae_mf_elbo(cx, ms$S, ev)
}

# posterior mean and variance of the linear predictor of every cell of `X` (the fit's
# own when NULL) under the variational posterior, as [n_row, n_col, T] arrays: the
# Gaussian additive part at the expected weights plus the product's mean and variance
.ae_mf_predictive <- function(object, X = NULL) {
	# restore a withheld product (.ae_withhold_product) to its variational means
	object <- .ae_product_restore(object)
	add <- object
	add$R <- 0L; add$U <- NULL; add[["V"]] <- NULL; add$L <- NULL
	pr <- .ae_eb_predictive(add, X, displace = FALSE)
	if (object$R > 0L) {
		post <- object$eb$factor_posterior
		pm <- as.matrix(object$UVPM)
		vp <- if (identical(post$kind, "symmetric")) {
			.ae_mf_prod_var_sym(post$F, post$sg, post$Sig)
		} else {
			.ae_mf_prod_var(object$U, object[["V"]], post$Sig, post$Om)
		}
		Tt <- dim(pr$mean)[3L]
		pr$mean <- pr$mean + array(rep(pm, Tt), dim(pr$mean))
		pr$var <- pr$var + array(rep(vp, Tt), dim(pr$var))
	}
	pr
}
