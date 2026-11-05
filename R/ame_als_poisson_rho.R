# latent Poisson reciprocity from the cross-dyad identity
# cov(y_ij / mu_ij - 1, y_ji / mu_ji - 1) = exp(rho * ve) - 1, with inverse
# marginal-variance weights. rho is NA unless ve is distinguishable from zero;
# the raw moment and covariance are kept for the count jackknife
.ae_poisson_pair_rho <- function(y, eta, variance, overdispersed = TRUE) {
	context <- .lame_binary_dynamic_collapsed_rho_context(y)
	keep <- context$observed1 & context$observed2 &
		is.finite(eta[context$one]) & is.finite(eta[context$two])
	n_pairs <- sum(keep)
	if (!is.finite(variance) || variance < 0) {
		cli::cli_abort(
			"Latent count variance must be finite and nonnegative.",
			class = "lame_als_numerical"
		)
	}
	if (n_pairs == 0L) {
		return(list(
			rho = NA_real_,
			identified = FALSE,
			n_pairs = n_pairs,
			boundary = FALSE,
			method = "weighted count covariance",
			reason = "no reciprocal pairs"
		))
	}
	if (variance > log(.Machine$double.xmax)) {
		cli::cli_abort(
			"The fitted count covariance exceeds the available numeric range.",
			class = "lame_als_numerical"
		)
	}
	one <- context$one[keep]
	two <- context$two[keep]
	log_phi <- if (variance == 0) -Inf else log(expm1(variance))
	terms <- function(index) {
		linear <- eta[index]
		x <- log_phi + linear
		log_denominator <- pmax(x, 0) + log1p(exp(-abs(x)))
		weight <- exp(linear - log_denominator)
		residual <- y[index] * exp(-log_denominator) - weight
		list(weight = weight, residual = residual)
	}
	a <- terms(one)
	b <- terms(two)
	# common scaling cancels from the ratio and prevents overflow
	scale <- max(abs(a$residual), abs(b$residual), a$weight, b$weight)
	if (!is.finite(scale) || scale == 0) {
		cli::cli_abort(
			"The fitted count covariance cannot be evaluated at this scale.",
			class = "lame_als_numerical"
		)
	}
	numerator <- sum((a$residual / scale) * (b$residual / scale))
	denominator <- sum((a$weight / scale) * (b$weight / scale))
	if (denominator == 0) {
		cli::cli_abort(
			"The reciprocal count information underflowed at the fitted means.",
			class = "lame_als_numerical"
		)
	}
	moment <- numerator / denominator
	covariance <- if (moment > -1) log1p(moment) else NA_real_
	if (variance == 0 || !isTRUE(overdispersed)) {
		return(list(
			rho = NA_real_,
			identified = FALSE,
			n_pairs = n_pairs,
			boundary = FALSE,
			moment = moment,
			covariance = covariance,
			method = "weighted count covariance",
			reason = if (variance == 0) "zero latent variance" else
				"latent variance not distinguishable from zero"
		))
	}
	lower <- expm1(-variance)
	upper <- expm1(variance)
	rho <- if (moment <= lower) -1 else if (moment >= upper) 1 else log1p(moment) / variance
	list(
		rho = rho,
		identified = TRUE,
		n_pairs = n_pairs,
		boundary = abs(rho) == 1,
		moment = moment,
		covariance = covariance,
		method = "weighted count covariance"
	)
}

# one-sided test of the latent count variance against zero: under
# equidispersion the pearson statistic has mean n_obs - df and variance
# about sum(2 + 1 / mu) * (n_obs - df) / n_obs
.ae_pois_overdispersion <- function(Yarr, eta, df_mean, level = 0.05) {
	y <- as.vector(Yarr)
	eta <- as.vector(eta)
	ok <- is.finite(y) & is.finite(eta)
	y <- y[ok]
	mu <- exp(eta[ok])
	n_obs <- length(y)
	df_mean <- .ae_default(df_mean, 0)
	resid_df <- max(n_obs - df_mean, 1)
	if (n_obs < 3L || any(mu <= 0)) {
		return(list(
			statistic = NA_real_,
			p_value = NA_real_,
			overdispersed = FALSE,
			level = level
		))
	}
	pearson <- sum((y - mu)^2 / mu)
	null_var <- sum(2 + 1 / mu) * resid_df / n_obs
	z <- (pearson - resid_df) / sqrt(null_var)
	p <- stats::pnorm(z, lower.tail = FALSE)
	list(
		statistic = z,
		p_value = p,
		overdispersed = is.finite(p) && p < level,
		level = level
	)
}

# refresh reciprocity summaries on the latent scale used for simulation
.ae_update_response_rho <- function(object, link) {
	if (isTRUE(object$symmetric) || object$mode == "bipartite" || is.null(link)) {
		return(object)
	}
	if (object$family == "binary") {
		object <- .ae_update_binary_rho(object, link)
		estimate <- function(y, eta) .ae_binary_pair_rho(y, eta, link)
	} else if (object$family == "poisson" && link == "log") {
		object$rho_working <- unname(object$VC["rho"])
		variance <- unname(object$VC["ve"])
		shift <- .ae_default(object$pois_latent_var, 0) / 2
		eta <- array(unlist(object$EZ), dim(object$Y))
		object$overdispersion <- .ae_pois_overdispersion(
			object$Y,
			eta + shift,
			object$df_mean
		)
		over <- variance > 0 && isTRUE(object$overdispersion$overdispersed)
		estimate <- function(y, eta) .ae_poisson_pair_rho(
			y,
			eta + shift,
			variance,
			overdispersed = over
		)
		object$rho_estimation <- estimate(object$Y, eta)
		object$VC["rho"] <- object$rho_estimation$rho
		# an EB count fit with a pair stage reports its pair-model correlation
		pr <- object$eb$pair
		if (!is.null(pr) && over) {
			object$rho_estimation <- c(
				list(
					rho = pr$rho,
					identified = pr$n_pairs > 0L,
					n_pairs = pr$n_pairs,
					boundary = abs(pr$rho) >= 0.98,
					sd = pr$sd,
					method = pr$method
				),
				object$rho_estimation[c("moment", "covariance")]
			)
			object$VC["rho"] <- pr$rho
		}
	} else if (object$family == "normal" && identical(link, "identity")) {
		# a normal EB fit with a pair stage reports the posterior mean of rho from its pair model
		pr <- object$eb$pair
		if (!is.null(pr)) {
			object$rho_working <- unname(object$VC["rho"])
			object$VC["rho"] <- pr$rho
			object$rho_estimation <- list(
				rho = pr$rho,
				identified = pr$n_pairs > 0L,
				n_pairs = pr$n_pairs,
				boundary = abs(pr$rho) >= 0.98,
				sd = pr$sd,
				method = "bivariate normal pair likelihood"
			)
		}
		return(object)
	} else return(object)
	if (!is.null(object$rho_path)) {
		object$rho_path_working <- object$rho_path
		eta <- array(unlist(object$EZ), dim(object$Y))
		path <- lapply(seq_len(dim(object$Y)[3L]), function(t) {
			estimate(object$Y[, , t, drop = FALSE], eta[, , t, drop = FALSE])
		})
		object$rho_path <- setNames(
			vapply(path, `[[`, numeric(1), "rho"),
			names(object$rho_path_working)
		)
		object$rho_path_estimation <- path
		object$residual_rho_path <- object$rho_path
		object$rho_path_role <- "residual_diagnostic"
		object$dynamic_rho <- FALSE
	}
	object
}

# whether simulation needs a reciprocity the fit did not identify; a Poisson
# fit with no detectable latent variance simulates independent errors instead
.ae_rho_required <- function(object) {
	object$family %in% c("binary", "poisson") &&
		!is.null(object$rho_estimation) && !isTRUE(object$rho_estimation$identified) &&
		!(object$family == "poisson" && (isTRUE(object$VC["ve"] == 0) ||
				identical(
					object$rho_estimation$reason,
					"latent variance not distinguishable from zero"
				)))
}

# text shown by print and summary for a reciprocity that is not estimated
.ae_rho_not_estimated_text <- function(object) {
	est <- object$rho_estimation
	if (is.null(est) || isTRUE(est$identified)) return(NULL)
	reason <- est$reason %||% "the reciprocal pairs do not identify it"
	txt <- sprintf("rho: not estimated (%s", reason)
	od <- object$overdispersion
	if (identical(est$reason, "latent variance not distinguishable from zero") &&
			!is.null(od) && is.finite(od$p_value)) {
		txt <- sprintf(
			"%s: ve = %.3g, overdispersion test p = %.3g",
			txt,
			unname(object$VC["ve"]),
			od$p_value
		)
	}
	paste0(txt, ")")
}
