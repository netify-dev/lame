# interweaving step for the path hyperparameters of a dynamic_beta fit (lame()).
#
# sigma_g and rho_g are drawn from their exact conditionals given the path every
# sweep; this step adds the non-centred view, which mixes where the data say
# little about the path:
#     beta_1 fixed,  beta_t = f_t beta_{t-1} + sigma_g sqrt(Lambda q_t) w_t,
#     w_t ~ N(0, 1)  (f_t = rho^gap, q_t = (1 - rho^(2 gap)) / (1 - rho^2);
#     rw1: f = 1, q = gap; rw2 and matern32: beta = sigma_g sqrt(Lambda) C w)
# given w and the data, (sigma_g, rho_g) has density
#     p(sigma_g, rho_g) x exp(sum_t h_t' beta_t - beta_t' Omega_t beta_t / 2)
# with (Omega_t, h_t) from beta_dynamic_stats_cpp. a symmetric walk on
# (log sigma_g, logit rho_g) rebuilds the path from the same w and is accepted by
# Metropolis; the path draw that follows completes the interweaving. the walk
# scale comes from a fixed ladder, so the step keeps no state.

# part of a path the rw2 prior leaves flat: its projection on constants and
# linear functions of time (`tpos` the time positions)
.lame_beta_rw2_flat <- function(x, tpos) {
	B <- cbind(1, tpos - mean(tpos))
	B %*% solve(crossprod(B), crossprod(B, x))
}

# path of a block under new hyperparameters with its standardised innovations
# held fixed; x is the (periods x coefficients) path, lam the Lambda diagonal,
# gaps the time gaps. a matern32 path scales with sigma; an rw2 path scales its
# non-flat part and keeps the flat part, whose prior enters the ratio
.lame_beta_nc_path <- function(x, kind, sig0, rho0, sig1, rho1, lam, gaps) {
	if (identical(kind, "matern32")) return(x * (sig1 / sig0))
	if (identical(kind, "rw2")) {
		flat <- .lame_beta_rw2_flat(x, c(0, cumsum(gaps)))
		return(flat + (sig1 / sig0) * (x - flat))
	}
	nt <- nrow(x)
	z0 <- .lame_beta_gap_transition(rho0, gaps)
	z1 <- .lame_beta_gap_transition(rho1, gaps)
	sl <- sqrt(lam)
	out <- x
	for (t in 2:nt) {
		w <- (x[t, ] - z0$f[t - 1L] * x[t - 1L, ]) / (sig0 * sl * sqrt(z0$q[t - 1L]))
		out[t, ] <- z1$f[t - 1L] * out[t - 1L, ] + sig1 * sl * sqrt(z1$q[t - 1L]) * w
	}
	out
}

# log density of a path under the likelihood sufficient statistics (up to a constant)
.lame_beta_nc_loglik <- function(B, st) {
	nt <- nrow(B)
	q <- 0
	for (t in seq_len(nt)) {
		b <- B[t, ]
		q <- q + sum(st$h[t, ] * b) - 0.5 * sum(b * (st$Omega[, , t] %*% b))
	}
	q
}

# one sweep of the step over the blocks. `st` is list(Omega, h) from
# beta_dynamic_stats_cpp, `prior_scale` the inverse-gamma scale of sigma^2 (one
# or per block), `rho_pool` NULL or list(mu, tau) of a pooled logistic-normal
# persistence, `movable_rho` whether rho_g is a parameter of the kind (ar1)
.lame_beta_nc_step <- function(
	path,
	st,
	group_id,
	n_groups,
	kind,
	Lambda_diag,
	sigma_g,
	rho_g,
	prior_shape,
	prior_scale,
	rho_mean,
	rho_sd,
	rho_lower,
	rho_upper,
	rho_pool = NULL,
	movable_rho = identical(kind, "ar1"),
	gaps = rep(1, nrow(path) - 1L),
	tries = 6L,
	scales = c(0.03, 0.1, 0.3, 0.8)
) {
	out <- list(path = path, sigma = sigma_g, rho = rho_g, attempts = 0L, accepts = 0L)
	nt <- nrow(path)
	if (nt < 2L || ncol(path) == 0L) return(out)
	prior_scale <- rep_len(prior_scale, n_groups)
	rho_mean <- rep_len(rho_mean, n_groups)
	rho_sd <- rep_len(rho_sd, n_groups)
	ll_cur <- .lame_beta_nc_loglik(path, st)
	if (!is.finite(ll_cur)) return(out)
	rng <- rho_upper - rho_lower
	log_prior <- function(g, sig, rho) {
		lp <- -2 * prior_shape * log(sig) - prior_scale[g] / sig^2
		if (movable_rho) {
			u <- (rho - rho_lower) / rng
			if (!is.finite(u) || u <= 0 || u >= 1) return(-Inf)
			lp <- lp + if (is.null(rho_pool)) {
				stats::dnorm(rho, rho_mean[g], rho_sd[g], log = TRUE) + log(u) + log1p(-u)
			} else stats::dnorm(stats::qlogis(u), rho_pool$mu, rho_pool$tau, log = TRUE)
		}
		lp
	}
	for (g in seq_len(n_groups)) {
		idx <- which(group_id == g)
		if (!length(idx)) next
		lam <- Lambda_diag[idx]
		sig0 <- sigma_g[g]
		rho0 <- rho_g[g]
		# the flat part of an rw2 path has prior N(0, sigma^2 Lambda / ridge) in
		# each of its two dimensions, which depends on sigma
		flat_ss <- NULL
		if (identical(kind, "rw2")) {
			if (nt < 3L) next
			fl <- .lame_beta_rw2_flat(out$path[, idx, drop = FALSE], c(0, cumsum(gaps)))
			flat_ss <- colSums(fl^2) / lam
		}
		ridge <- eval(formals(.rw2_precision)$jitter)
		log_prior_g <- function(sig, rho) {
			lp <- log_prior(g, sig, rho)
			if (!is.null(flat_ss)) lp <- lp - 2 * length(idx) * log(sig) - ridge * sum(flat_ss) / (2 * sig^2)
			lp
		}
		lp_cur <- log_prior_g(sig0, rho0)
		if (!is.finite(lp_cur)) next
		for (k in seq_len(tries)) {
			sc <- scales[sample.int(length(scales), 1L)]
			# move the scale, the persistence or both (state-independent, so the
			# proposal stays symmetric)
			mv <- if (movable_rho) sample.int(3L, 1L) else 1L
			sig1 <- if (mv != 2L) sig0 * exp(sc * stats::rnorm(1L)) else sig0
			rho1 <- rho0
			if (mv >= 2L) {
				th <- stats::qlogis((rho0 - rho_lower) / rng) + 3 * sc * stats::rnorm(1L)
				rho1 <- rho_lower + rng * stats::plogis(th)
			}
			lp_new <- log_prior_g(sig1, rho1)
			if (!is.finite(lp_new)) next
			x1 <- tryCatch(
				.lame_beta_nc_path(out$path[, idx, drop = FALSE], kind, sig0, rho0, sig1, rho1, lam, gaps),
				error = function(e) NULL
			)
			out$attempts <- out$attempts + 1L
			if (is.null(x1) || !all(is.finite(x1))) next
			B1 <- out$path
			B1[, idx] <- x1
			ll_new <- .lame_beta_nc_loglik(B1, st)
			if (is.finite(ll_new) && log(stats::runif(1L)) < ll_new - ll_cur + lp_new - lp_cur) {
				out$path <- B1
				ll_cur <- ll_new
				lp_cur <- lp_new
				sig0 <- sig1
				rho0 <- rho1
				out$accepts <- out$accepts + 1L
			}
		}
		out$sigma[g] <- sig0
		out$rho[g] <- rho0
	}
	out
}
