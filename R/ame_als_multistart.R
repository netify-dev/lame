# multistart ALS: the low-rank objective is non-convex for R > 0, so the BCD runs from
# the deterministic cold start plus random low-rank starts and keeps the best one (lowest
# SSE, IRLS objective, or Laplace evidence for an empirical-Bayes fit with a latent
# rank). random starts use `seed` and restore the caller's RNG state

# the n_starts initialisations: the deterministic cold start (NULL) plus n_starts - 1
# random low-rank starts sized from the working response
.ae_multistart_starts <- function(Z, R, symmetric, n_starts, seed) {
	starts <- vector("list", n_starts)
	starts[1] <- list(NULL) # deterministic cold start
	if (n_starts <= 1L) return(starts)

	d <- dim(Z); nr <- d[1]; nc <- d[2]
	z_sd <- stats::sd(Z[is.finite(Z)])
	if (!is.finite(z_sd) || z_sd == 0) z_sd <- 1
	# random factors give u_i' v_j an sd of about 0.5 sd(z), since sd(u' v) = sqrt(R) sc^2
	sc <- sqrt(0.5 * z_sd / sqrt(R))

	# save and restore the caller's RNG stream
	had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	old_seed <- if (had_seed) {
		get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
	} else {
		NULL
	}
	on.exit({
		if (had_seed) {
			assign(".Random.seed", old_seed, envir = .GlobalEnv)
		} else if (exists(
			".Random.seed",
			envir = .GlobalEnv,
			inherits = FALSE
		)) {
			rm(".Random.seed", envir = .GlobalEnv)
		}
	}, add = TRUE)
	set.seed(seed)

	for (s in 2:n_starts) {
		Us <- matrix(stats::rnorm(nr * R, 0, sc), nr, R)
		if (symmetric) {
			starts[[s]] <- list(U = Us, V = Us, L = rep(sc, R))
		} else {
			Vs <- matrix(stats::rnorm(nc * R, 0, sc), nc, R)
			starts[[s]] <- list(U = Us, V = Vs)
		}
	}
	starts
}

# the best finite criterion, preferring converged runs
.ae_multistart_choose <- function(values, converged) {
	eligible <- is.finite(values) & converged
	pool <- which(if (any(eligible)) eligible else is.finite(values))
	if (!length(pool)) {
		cli::cli_abort(
			"No ALS start returned a finite selection criterion.",
			class = "lame_als_numerical"
		)
	}
	pool[which.min(values[pool])]
}

# the value an IRLS start is compared by: the objective its reweighting loop minimises,
# plus the hyperprior of Sab for an empirical-Bayes fit (the raw deviance always
# favours the least-shrunk optimum)
.ae_irls_objective <- function(ir, ab_prior = NULL) {
	obj <- ir$best_obj
	if (!is.null(ab_prior) && is.finite(obj)) {
		obj <- obj + tryCatch(
			.ae_eb_hyper(ir$fit$sab, ab_prior),
			error = function(e) NA_real_
		)
	}
	obj
}

# Laplace marginal likelihood of an empirical-Bayes start (minus twice its log,
# ame_als_eb.R), or NA; it integrates the factors, so starts at one optimum agree
.ae_irls_evidence <- function(ir) {
	ev <- ir$fit$eb_vcov$evidence$value
	if (is.null(ev) || !is.finite(ev)) NA_real_ else ev
}

# deviance units above which a start's distance from the kept start (.ae_irls_distance)
# counts as a materially different optimum
.ae_multistart_material <- 1

# distance between the fitted predictors of two IRLS starts on the deviance scale: the
# squared predictor difference weighted by the mean Fisher weights, independent of the
# factor scales and their representation
.ae_irls_distance <- function(ir, ref) {
	sum((ir$W + ref$W) / 2 * (ir$eta - ref$eta)^2, na.rm = TRUE)
}

# the start an IRLS multistart keeps, the value each start is compared by, each start's
# distance from the kept one, and which converged starts reached a different optimum
.ae_irls_select <- function(irs, ab_prior = NULL, tol = .ae_multistart_material) {
	value <- vapply(irs, .ae_irls_objective, numeric(1), ab_prior = ab_prior)
	# a start followed by a variational refinement converged when the refinement did
	converged <- vapply(irs, function(z) {
		if (!is.null(z$selection_converged)) isTRUE(z$selection_converged)
		else isTRUE(z$stationarity$pass)
	}, logical(1))
	criterion <- if (is.null(ab_prior)) "penalised deviance" else "empirical-Bayes objective"
	# an empirical-Bayes start with a latent rank is rated by its Laplace marginal
	# likelihood when every converged start has one
	if (!is.null(ab_prior) && any(converged) &&
			all(vapply(irs, function(z) !is.null(z$fit$U), logical(1)))) {
		ev <- vapply(irs, .ae_irls_evidence, numeric(1))
		if (all(is.finite(ev[converged]))) {
			value <- ev
			criterion <- "Laplace marginal likelihood"
		}
	}
	best <- .ae_multistart_choose(value, converged)
	distance <- vapply(irs, .ae_irls_distance, numeric(1), ref = irs[[best]])
	different <- converged & is.finite(distance) & distance > tol
	list(
		value = value,
		converged = converged,
		best = best,
		excess = value - value[best],
		distance = distance,
		different = different,
		material = any(different),
		criterion = criterion
	)
}

# relative Frobenius distance of a start's product from the kept start's (a vanishing
# kept product is measured against sqrt(eps)); NA without a product
.ae_product_distance <- function(m, ref) {
	if (is.null(m) || is.null(ref)) return(NA_real_)
	sqrt(sum((m - ref)^2)) / max(sqrt(sum(ref^2)), sqrt(.Machine$double.eps))
}

# the coefficients each start reports and the distance of its product from the kept
# start's; `assemble(j)` builds start j as the kept fit is (.ae_assemble), and a start
# that cannot be assembled gets a row of NA
.ae_multistart_summary <- function(assemble, n_starts, best) {
	grab <- function(j) {
		a <- tryCatch(suppressWarnings(assemble(j)), error = function(e) NULL)
		if (is.null(a)) NULL else list(coefficients = a$coefficients, product = a$UVPM)
	}
	ref <- grab(best)
	if (is.null(ref)) return(NULL)
	nm <- names(ref$coefficients)
	coefficients <- matrix(NA_real_, n_starts, length(nm), dimnames = list(NULL, nm))
	distance <- rep(NA_real_, n_starts)
	for (j in seq_len(n_starts)) {
		s <- if (j == best) ref else grab(j)
		if (is.null(s) || length(s$coefficients) != length(nm)) next
		coefficients[j, ] <- s$coefficients
		distance[j] <- if (j == best) {
			if (is.null(ref$product)) NA_real_ else 0
		} else {
			.ae_product_distance(s$product, ref$product)
		}
	}
	list(coefficients = coefficients, product_distance = distance)
}

.ae_multistart_fit <- function(
	Z,
	X,
	R,
	symmetric,
	max_iter,
	tol,
	lowrank_method,
	linear_solver,
	n_starts,
	seed,
	certify = TRUE,
	assemble = NULL,
	unit_rules = FALSE
) {

	bcd <- function(init) {
		.ae_bcd_fit(
			Z,
			X,
			R,
			symmetric,
			max_iter,
			tol,
			init = init,
			lowrank_method = lowrank_method,
			linear_solver = linear_solver,
			certify = certify,
			unit_rules = unit_rules
		)
	}

	starts <- .ae_multistart_starts(Z, R, symmetric, n_starts, seed)
	fits <- lapply(starts, bcd)

	sse <- vapply(fits, function(f) f$sse, numeric(1))
	converged <- vapply(
		fits,
		function(f) isTRUE(f$converged) &&
			!identical(f$add_converged, FALSE),
		logical(1)
	)
	best <- .ae_multistart_choose(sse, converged)
	out <- fits[[best]]
	out$multistart_sse <- sse
	out$multistart_selected <- best
	out$multistart_converged <- converged
	out$multistart_criterion <- "working SSE"
	out$multistart_gaussian_diagnostics <- lapply(fits, `[[`, "gaussian_diagnostics")
	if (!is.null(assemble)) {
		summ <- .ae_multistart_summary(function(j) assemble(fits[[j]]), length(fits), best)
		out$multistart_coefficients <- summ$coefficients
		out$multistart_product_distance <- summ$product_distance
	}
	out
}
