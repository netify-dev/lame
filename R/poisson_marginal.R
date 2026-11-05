# integrate the Poisson-lognormal mass on a local normal scale; log concavity
# bounds both omitted tails relative to a positive lower bound on the mass
.lame_pois_marginal_loglik <- function(y_r, eta_obs, ve_s) {
	if (length(y_r) == 0L) return(numeric(0))
	if (length(ve_s) != 1L || !is.finite(ve_s) || ve_s <= 0)
		cli::cli_abort("Poisson marginal log likelihood requires a finite positive observation variance.")
	if (length(y_r) != length(eta_obs) || any(!is.finite(eta_obs)) ||
			any(!is.finite(y_r) | y_r < 0 | y_r != floor(y_r)))
		cli::cli_abort("Poisson log likelihood requires finite predictors and nonnegative integer counts.")
	# integration runs to a fixed accuracy criterion
	sd_prior <- sqrt(ve_s)
	t_mode <- (log(y_r + pmax(eta_obs, 0) / ve_s + 1e-8) - eta_obs) / sd_prior
	converged <- FALSE
	for (iteration in seq_len(200L)) {
		log_vlambda <- log(ve_s) + eta_obs + sd_prior * t_mode
		w <- stats::plogis(-log_vlambda)
		step <- t_mode * w + stats::plogis(log_vlambda) / sd_prior - sd_prior * y_r * w
		log_scale <- -.5 * (pmax(log_vlambda, 0) + log1p(exp(-abs(log_vlambda))))
		if (any(!is.finite(step)) || any(!is.finite(log_scale)))
			cli::cli_abort("Poisson likelihood mode cannot be resolved numerically.")
		if (all(abs(step) <= 1e-7 * exp(log_scale))) { converged <- TRUE; break }
		t_mode <- t_mode - step
	}
	if (!converged) cli::cli_abort("Poisson likelihood mode did not converge.")
	z_mode <- eta_obs + sd_prior * t_mode
	scale <- exp(log_scale)
	if (any(!is.finite(scale) | scale <= 0)) cli::cli_abort("Poisson integration scale cannot be resolved numerically.")
	remainder <- function(x) {
		out <- expm1(x) - x
		small <- abs(x) < 1e-3
		z <- x[small]
		out[small] <- z * z * (.5 + z * (1 / 6 + z * (1 / 24 + z * (1 / 120 + z * (1 / 720 + z / 5040)))))
		out
	}
	vapply(seq_along(y_r), function(i) {
		y <- y_r[i]; tm <- t_mode[i]; zm <- z_mode[i]; sc <- scale[i]
		lambda <- exp(zm)
		log_pois <- .poisson_log_mass(y, zm)
		peak <- log_pois + stats::dnorm(tm, log = TRUE)
		if (!is.finite(peak) || !is.finite(lambda))
			cli::cli_abort("Poisson marginal log mass cannot be resolved numerically.")
		log_relative <- function(x) {
			dt <- sc * x; dz <- sd_prior * dt
			value <- numeric(length(x))
			near <- abs(dz) <= .5
			value[near] <- (y - lambda) * dz[near] - lambda * remainder(dz[near])
			if (any(!near)) value[!near] <- .poisson_log_mass(rep(y, sum(!near)), zm + dz[!near]) - log_pois
			value - tm * dt - .5 * dt * dt
		}
		derivative <- function(x) {
			dt <- sc * x; z <- zm + sd_prior * dt
			(sd_prior * (y - exp(z)) - tm - dt) * sc
		}
		# concavity gives f(x) >= min(f(-1), f(1)) for -1 <= x <= 1
		log_lower <- log(2) + min(log_relative(c(-1, 1)))
		if (!is.finite(log_lower)) cli::cli_abort("Poisson integral lower bound cannot be resolved numerically.")
		limits <- lapply(c(-1, 1), function(sign) {
			points <- sign * c(1, 8)
			for (j in seq_len(1024L)) {
				x <- utils::tail(points, 1L)
				g <- derivative(x)
				bound <- log_relative(x) - log(abs(g))
				if (!is.na(bound) && sign * g < 0 && bound < log_lower + log(1e-12)) return(points)
				if (!is.finite(x * 2)) break
				points <- c(points, x * 2)
			}
			cli::cli_abort("Poisson integral tails cannot be bounded numerically.")
		})
		cuts <- sort(c(limits[[1L]], 0, limits[[2L]]))
		integral <- error <- 0
		for (j in seq_len(length(cuts) - 1L)) {
			part <- stats::integrate(
				function(x) exp(log_relative(x)),
				cuts[j],
				cuts[j + 1L],
				rel.tol = 1e-10,
				abs.tol = 1e-12 / length(cuts),
				subdivisions = 1000L
			)
			integral <- integral + part$value
			error <- error + part$abs.error
		}
		if (!is.finite(integral) || integral <= 0 || !is.finite(error) || error > 1e-9 * integral)
			cli::cli_abort("Poisson marginal likelihood integration did not achieve its accuracy criterion.")
		peak + log(sc) + log(integral)
	}, numeric(1))
}
