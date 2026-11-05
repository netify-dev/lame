# split dyads into two halves by a circulant graph that ignores the outcome;
# each dyad's history stays together and both halves keep every actor
.lda_count_split <- function(y) {
	n <- dim(y)[1L]
	if (length(dim(y)) != 3L || dim(y)[2L] != n || n < 4L) {
		cli::cli_abort(
			"The dynamic count split needs at least four actors.",
			class = "lame_jk_split"
		)
	}
	gcd <- function(a, b) { while (b != 0L) { z <- b; b <- a %% b; a <- z }; a }
	stride <- 17L
	while (gcd(stride, n) != 1L) stride <- stride + 1L
	index <- ((seq_len(n) - 1L) * stride) %% n
	distance <- abs(outer(index, index, "-"))
	distance <- pmin(distance, n - distance)
	part <- distance > 0 & distance <= floor(n / 4)
	mask <- array(rep(part, dim(y)[3L]), dim(y))
	obs <- is.finite(y)
	first <- obs & mask
	second <- obs & !mask
	fraction <- c(sum(first), sum(second)) / sum(obs)
	if (any(!is.finite(fraction)) || any(fraction <= 0)) {
		cli::cli_abort(
			"Both count jackknife halves must contain observed dyads.",
			class = "lame_jk_split"
		)
	}
	list(mask = mask, fraction = fraction, stride = stride)
}

# jackknife-correct the variance and covariance, then convert to a correlation
# projected onto the admissible range; raw moments and boundary flags are returned
.lda_count_moments <- function(full, halves, fraction) {
	v <- unname(full$VC["ve"])
	hv <- vapply(halves, function(h) unname(h$VC["ve"]), numeric(1))
	# latent covariance from the unrestricted count moment, not a clipped correlation
	cov_of <- function(est, rho, ve) {
		cv <- if (is.list(est)) est$covariance else NULL
		if (!is.null(cv) && is.finite(cv)) return(cv)
		if (ve == 0 || !is.finite(rho)) 0 else rho * ve
	}
	fc <- cov_of(full$rho_estimation, unname(full$VC["rho"]), v)
	hc <- vapply(seq_along(halves), function(h) {
		cov_of(
			halves[[h]]$rho_estimation,
			unname(halves[[h]]$VC["rho"]),
			hv[h]
		)
	}, numeric(1))
	raw_v <- 2 * v - sum(fraction * hv)
	raw_c <- 2 * fc - sum(fraction * hc)
	new_v <- max(0, raw_v)
	new_c <- if (is.finite(raw_c)) min(new_v, max(-new_v, raw_c)) else NA_real_
	rho <- if (new_v > 0) new_c / new_v else NA_real_
	path <- full$rho_path
	raw_path <- if (is.null(path)) NULL else {
		path_cov <- function(f, ve) {
			vapply(seq_along(path), function(t) {
				est <- if (t <= length(f$rho_path_estimation)) f$rho_path_estimation[[t]]
				cov_of(est, unname(f$rho_path[t]), ve)
			}, numeric(1))
		}
		value <- path_cov(full, v)
		for (h in seq_along(halves)) {
			cp <- path_cov(halves[[h]], hv[h])
			if (h == 1L) half_path <- fraction[h] * cp else half_path <- half_path + fraction[h] * cp
		}
		2 * value - half_path
	}
	path <- if (is.null(raw_path)) NULL else if (new_v == 0) {
		setNames(rep(NA_real_, length(raw_path)), names(full$rho_path))
	} else pmax(-1, pmin(1, raw_path / new_v))
	list(
		variance = new_v,
		covariance = new_c,
		rho = rho,
		rho_path = path,
		unrestricted = list(variance = raw_v, covariance = raw_c, covariance_path = raw_path),
		boundary = raw_v <= 0 || (is.finite(raw_c) && abs(raw_c) >= new_v),
		path_boundary = if (is.null(raw_path)) NULL else abs(raw_path) >= new_v
	)
}

.lda_apply_count_moments <- function(object, moments, details = list()) {
	old_v <- unname(object$VC["ve"])
	new_v <- moments$variance
	shift <- (old_v - new_v) / 2
	object$raw <- object[c(
		"coefficients", "coef_path", "VC", "rho_path",
		"rho_estimation", "rho_path_estimation", "deviance", "residual_rho_path",
		"rho_path_role", "dynamic_rho"
	)]
	object$mu <- object$mu + shift
	object$coefficients["intercept"] <- object$coefficients["intercept"] + shift
	object$coef_path["intercept", ] <- object$coef_path["intercept", ] + shift
	object$beta_dynamic_per_t <- object$coef_path
	if (length(dim(object$BETA)) == 3L) {
		object$BETA[, "intercept", ] <- object$BETA[, "intercept", ] + shift
	} else object$BETA[, "intercept"] <- object$BETA[, "intercept"] + shift
	object$EZ <- lapply(object$EZ, function(x) x + shift)
	object$VC["ve"] <- new_v
	object$VC["rho"] <- moments$rho
	object$s2 <- new_v
	object$pois_latent_var <- new_v
	object$rho_path <- moments$rho_path
	object$residual_rho_path <- object$rho_path
	object$rho_path_role <- if (is.null(object$rho_path)) NULL else "residual_diagnostic"
	object$dynamic_rho <- FALSE
	# without detectable overdispersion the correlation is not identified
	flat <- new_v > 0 && !is.null(object$overdispersion) &&
		!isTRUE(object$overdispersion$overdispersed)
	if (flat) {
		object$VC["rho"] <- NA_real_
		if (!is.null(object$rho_path)) object$rho_path[] <- NA_real_
		object$residual_rho_path <- object$rho_path
	}
	meta <- function(rho, old, boundary) {
		if (flat) rho <- NA_real_
		list(
			rho = rho,
			identified = new_v > 0 && is.finite(rho),
			n_pairs = old$n_pairs,
			boundary = isTRUE(boundary) && !flat,
			method = "split-dyad latent covariance jackknife",
			reason = if (new_v == 0) "zero latent variance" else if (flat) "latent variance not distinguishable from zero" else if (!is.finite(rho)) "no reciprocal pairs" else NULL
		)
	}
	object$rho_estimation <- meta(moments$rho, object$rho_estimation, moments$boundary)
	if (!is.null(object$rho_path_estimation)) {
		object$rho_path_estimation <- lapply(seq_along(object$rho_path), function(t) {
			meta(object$rho_path[t], object$raw$rho_path_estimation[[t]], moments$path_boundary[t])
		})
	}
	object$bias_correct <- c(
		list(
			method = "jackknife",
			applied = TRUE,
			target = "latent count variance and reciprocal covariance",
			split = "reciprocal dyad histories",
			moments = moments
		),
		details
	)
	object
}

.lda_count_jackknife <- function(
	object,
	call_args,
	max_iter,
	tol,
	fitter = lame_dynamic_als
) {
	if (!isTRUE(object$converged) || !isTRUE(object$diagnostics$stationarity$pass)) {
		cli::cli_abort(
			"The full count fit must converge before bias correction.",
			class = "lame_jk_nonconvergence"
		)
	}
	split <- .lda_count_split(object$Y)
	halves <- attempts <- vector("list", 2L)
	summaries <- function(fit, warnings, budget) {
		if (inherits(fit, "error")) return(list(
			error = conditionMessage(fit),
			warnings = warnings,
			max_iter = budget
		))
		list(
			VC = fit$VC,
			rho_path = fit$rho_path,
			converged = fit$converged,
			diagnostics = fit$diagnostics,
			objective = fit$objective,
			iterations = fit$iterations,
			warnings = warnings,
			max_iter = budget
		)
	}
	for (h in 1:2) {
		y <- object$Y
		y[if (h == 1L) !split$mask else split$mask] <- NA_real_
		data <- .ae_dynamic_refit_data(object, y, object[["X"]])
		args <- call_args
		args$Y <- data$Y
		args$Xdyad <- data$Xdyad
		args$bootstrap <- 0L
		args$stability <- "none"
		args$verbose <- FALSE
		args$start <- .ae_static_seed(object)
		args$tol <- tol
		if ("bias_correct" %in% names(formals(fitter))) args$bias_correct <- "none"
		budget <- as.integer(max_iter)
		for (attempt in 1:2) {
			args$max_iter <- budget
			warnings <- character()
			half <- tryCatch(
				withCallingHandlers(
					do.call(fitter, args, quote = TRUE),
					warning = function(w) {
						warnings <<- c(warnings, conditionMessage(w))
						invokeRestart("muffleWarning")
					}
				),
				error = identity
			)
			attempts[[h]][[attempt]] <- summaries(half, warnings, budget)
			# a half is healthy when it converged and is stationary; its warnings
			# repeat the full fit's and are only kept in the attempt record
			healthy <- !inherits(half, "error") && isTRUE(half$converged) &&
				isTRUE(half$diagnostics$stationarity$pass)
			# retry a half that stopped short from its last iterate
			stopped <- inherits(half, "lame_als_convergence") && !is.null(half$fit)
			if (healthy || (inherits(half, "error") && !stopped)) break
			args$start <- .ae_static_seed(if (stopped) half$fit else half)
			budget <- as.integer(5 * max_iter)
		}
		if (!healthy) {
			cli::cli_abort(
				c(
					"The dynamic count jackknife half {h} did not pass its convergence checks.",
					"i" = "Inspect the half-fit attempts attached to this error; increase the iteration budget or use {.code bias_correct = \"none\"} to inspect the uncorrected fit."
				),
				class = "lame_jk_nonconvergence",
				half = h,
				attempts = attempts
			)
		}
		halves[[h]] <- half
	}
	moments <- .lda_count_moments(object, halves, split$fraction)
	.lda_apply_count_moments(
		object,
		moments,
		list(
			fraction = split$fraction,
			stride = split$stride,
			converged = vapply(halves, function(h) isTRUE(h$converged), logical(1)),
			attempts = attempts
		)
	)
}
