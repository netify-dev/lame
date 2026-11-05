# stable log probability of an interval under the standard normal
.ae_normal_interval_log <- function(lower, upper) {
	value <- rep(-Inf, length(lower))
	valid <- upper > lower
	left <- valid & upper <= 0
	right <- valid & lower >= 0
	middle <- valid & !left & !right
	high <- stats::pnorm(upper[left], log.p = TRUE)
	low <- stats::pnorm(lower[left], log.p = TRUE)
	value[left] <- high + log(-expm1(low - high))
	high <- stats::pnorm(lower[right], lower.tail = FALSE, log.p = TRUE)
	low <- stats::pnorm(upper[right], lower.tail = FALSE, log.p = TRUE)
	value[right] <- high + log(-expm1(low - high))
	value[middle] <- log(stats::pnorm(upper[middle]) - stats::pnorm(lower[middle]))
	value
}

# log bivariate normal cdf of .lame_binary_dynamic_bvn_log() computed by the compiled kernel
# (pbivnorm's scheme), with the same quadrature fallback below 1e-10
.ae_bvn_log <- function(a, b, rho) {
	if (identical(.ae_engine(), "r")) return(.lame_binary_dynamic_bvn_log(a, b, rho))
	rho <- rep(rho, length.out = length(a))
	stopifnot(all(is.finite(a)), all(is.finite(b)), all(abs(rho) < 1))
	if (!length(a)) return(numeric())
	p <- ae_pbivnorm_cpp(a, b, rho, .ae_threads())
	ans <- rep(NA_real_, length(p))
	good <- is.finite(p) & p > 1e-10 & p <= 1
	ans[good] <- log(p[good])
	bad <- which(!good)
	for (k in bad) ans[k] <- .lame_binary_dynamic_bvn_log_quad(a[k], b[k], rho[k])
	ans
}

# reciprocal binary dependence at the fitted marginal probabilities; the logit link maps
# probabilities to normal quantiles, and exact endpoint probabilities allow rho = -1 or 1
.ae_binary_pair_rho <- function(y, eta, link) {
	if (link == "logit") {
		eta <- stats::qnorm(stats::plogis(eta, log.p = TRUE), log.p = TRUE)
	}
	context <- .lame_binary_dynamic_collapsed_rho_context(y)
	keep <- context$observed1 & context$observed2 &
		is.finite(eta[context$one]) & is.finite(eta[context$two])
	n_pairs <- sum(keep)
	if (n_pairs == 0L) {
		return(list(
			rho = NA_real_, identified = FALSE, n_pairs = 0L,
			boundary = FALSE, log_likelihood = 0
		))
	}
	a <- context$sign1[keep] * eta[context$one[keep]]
	b <- context$sign2[keep] * eta[context$two[keep]]
	sign_pair <- context$sign1[keep] * context$sign2[keep]
	log_likelihood <- function(rho) {
		if (abs(rho) < 1) {
			return(sum(.ae_bvn_log(a, b, rho * sign_pair)))
		}
		same <- rho * sign_pair == 1
		value <- rep(-Inf, n_pairs)
		value[same] <- stats::pnorm(pmin(a[same], b[same]), log.p = TRUE)
		value[!same] <- .ae_normal_interval_log(-b[!same], a[!same])
		sum(value)
	}
	margin <- sqrt(.Machine$double.eps)
	optimum <- stats::optimize(
		function(rho) -log_likelihood(rho),
		c(-1 + margin, 1 - margin),
		tol = 1e-7
	)
	candidates <- c(-1, optimum$minimum, 0, 1)
	values <- vapply(candidates, log_likelihood, numeric(1))
	if (all(is.finite(values)) &&
			diff(range(values)) <= 64 * .Machine$double.eps * max(1, abs(values))) {
		return(list(
			rho = NA_real_,
			identified = FALSE,
			n_pairs = n_pairs,
			boundary = FALSE,
			log_likelihood = max(values)
		))
	}
	best <- which.max(values)
	if (!is.finite(values[best])) {
		cli::cli_abort(
			"The reciprocal binary likelihood is not finite at the fitted marginal probabilities.",
			class = "lame_als_numerical"
		)
	}
	list(
		rho = candidates[best],
		identified = TRUE,
		n_pairs = n_pairs,
		boundary = abs(candidates[best]) == 1,
		log_likelihood = values[best]
	)
}

# refresh dependence after a full fit or a split-sample correction so the likelihood
# and simulations use the reported marginal predictor
.ae_update_binary_rho <- function(object, link) {
	if (object$family != "binary" || isTRUE(object$symmetric) ||
			object$mode == "bipartite" || is.null(link)) return(object)
	one <- unlist(lapply(object$residuals, function(x) x[upper.tri(x)]))
	two <- unlist(lapply(object$residuals, function(x) t(x)[upper.tri(x)]))
	keep <- is.finite(one) & is.finite(two)
	object$rho_working <- if (sum(keep) > 2L) {
		suppressWarnings(stats::cor(one[keep], two[keep]))
	} else NA_real_
	# a fit with a dyadic correlation (ame_als_pair.R) reports its model's correlation
	pr <- object$eb$pair
	if (!is.null(pr)) {
		object$VC["rho"] <- pr$rho
		object$rho_estimation <- list(
			rho = pr$rho,
			identified = pr$n_pairs > 0L,
			n_pairs = pr$n_pairs,
			boundary = abs(pr$rho) >= 0.98,
			sd = pr$sd,
			method = "variational pair likelihood"
		)
		return(object)
	}
	eta <- array(unlist(object$EZ), dim(object$Y))
	estimate <- .ae_binary_pair_rho(object$Y, eta, link)
	object$VC["rho"] <- estimate$rho
	object$rho_estimation <- estimate
	object
}
