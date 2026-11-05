# weighted innovation sums on the centered sender and receiver spaces.
.lame_robust_actor_ss <- function(x, rho) {
	n <- nrow(x) %/% 2L
	nt <- ncol(x)
	q <- x
	q[, 1L] <- sqrt(1 - rho^2) * x[, 1L]
	if (nt > 1L) q[, -1L] <- x[, -1L, drop = FALSE] - rho * x[, -nt, drop = FALSE]
	rbind(
		colSums(q[seq_len(n), , drop = FALSE]^2),
		colSums(q[n + seq_len(n), , drop = FALSE]^2)
	)
}
.lame_robust_actor_q <- function(x, rho, tau) sum(tau * .lame_robust_actor_ss(x, rho))
.lame_robust_actor_rho_log <- function(rho, x, variance, dimension, tau) {
	if (!is.finite(rho) || abs(rho) >= .995) return(-Inf)
	dimension / 2 * log1p(-rho^2) - .lame_robust_actor_q(x, rho, tau) / (2 * variance) +
		dnorm(rho, .8, .3, log = TRUE)
}
.lame_robust_actor_rho_draw <- function(rho, x, variance, dimension, tau) {
	for (k in 1:2) {
		candidate <- rnorm(1, rho, .08)
		if (log(runif(1)) < .lame_robust_actor_rho_log(candidate, x, variance, dimension, tau) -
				.lame_robust_actor_rho_log(rho, x, variance, dimension, tau)) rho <- candidate
	}
	rho
}
# backend of the robust continuous-network model: actor paths and levels sum
# to zero, and coefficients and variances have proper priors. errors are
# Gaussian with unit precision or a Gamma(nu/2, nu/2) precision per directed
# pair shared across periods
.lame_robust_ar_q <- function(x, rho) {
	x <- as.matrix(x)
	tt <- ncol(x)
	(1 - rho^2) * sum(x[, 1]^2) + if (tt > 1L) sum((x[, -1, drop = FALSE] - rho * x[, -tt, drop = FALSE])^2) else 0
}
.lame_robust_rho_log <- function(rho, x, variance, dimension) {
	if (!is.finite(rho) || abs(rho) >= .995) {
		return(-Inf)
	}
	dimension / 2 * log1p(-rho^2) - .lame_robust_ar_q(x, rho) / (2 * variance) + dnorm(rho, .8, .3, log = TRUE)
}
.lame_robust_rho_draw <- function(rho, x, variance, dimension) {
	for (k in 1:2) {
		proposal <- rnorm(1, rho, .08)
		if (log(runif(1)) < .lame_robust_rho_log(
			proposal,
			x,
			variance,
			dimension
		) - .lame_robust_rho_log(rho, x, variance, dimension)) rho <- proposal
	}
	rho
}
.lame_robust_ig <- function(shape, rate) 1 / rgamma(1, shape = shape, rate = rate)
.lame_robust_ar_prior <- function(d, tt, rho, variance, center = FALSE) {
	x <- matrix(0, d, tt)
	x[, 1] <- rnorm(d, sd = sqrt(variance / (1 - rho^2)))
	if (tt > 1L) for (t in 2:tt) x[, t] <- rho * x[, t - 1] + rnorm(d, sd = sqrt(variance))
	if (center) x <- sweep(x, 2, colMeans(x))
	x
}
.lame_robust_prior <- function(Y, p) {
	v <- var(as.numeric(Y), na.rm = TRUE)
	if (!is.finite(v) || v <= 0) v <- 1
	list(
		ve_shape = 2,
		ve_rate = v,
		ar_shape = 2,
		mu_rate = .1 * v,
		ab_rate = .1 * v,
		uv_rate = .1 * sqrt(v),
		level_df = 4,
		level_scale = diag(v, 2),
		beta_precision = rep(.01 / v, p),
		nu = 4,
		actor_df = Inf,
		observation_df = Inf,
		training_variance = v
	)
}
.lame_robust_components <- function(s, D) {
	n <- D$n
	tt <- D$tt
	xb <- array(if (D$p) drop(D$X %*% s$beta) else 0, c(n, n, tt))
	mu <- array(rep(s$mu, each = n * n), c(n, n, tt))
	a <- b <- uv <- array(0, c(n, n, tt))
	for (t in 1:tt) {
		a[, , t] <- s$A[, t]
		b[, , t] <- matrix(rep(s$B[, t], each = n), n, n)
		if (D$rank) uv[, , t] <- tcrossprod(matrix(s$U[, , t], n, D$rank), matrix(s$V[, , t], n, D$rank))
	}
	list(xb = xb, mu = mu, a = a, b = b, uv = uv, eta = xb + mu + a + b + uv)
}
.lame_robust_initialize <- function(D, prior) {
	n <- D$n
	tt <- D$tt
	rr <- D$rank
	v <- prior$training_variance
	Sab <- solve(rWishart(1, prior$level_df, solve(prior$level_scale))[, , 1])
	lev <- matrix(rnorm(n * 2), n, 2) %*% chol(Sab)
	lev <- sweep(lev, 2, colMeans(lev))
	s <- list(
		beta = rnorm(D$p, 0, .2 * sqrt(v)),
		mu = rep(mean(D$Y[D$observed]), tt) + rnorm(tt, 0, .02 * sqrt(v)),
		A = .lame_robust_ar_prior(n, tt, .6, .03 * v, TRUE) + lev[, 1],
		B = .lame_robust_ar_prior(n, tt, .6, .03 * v, TRUE) + lev[, 2],
		a0 = lev[, 1],
		b0 = lev[, 2],
		U = array(rnorm(n * rr * tt, sd = sqrt(.1 * sqrt(v))), c(n, rr, tt)),
		V = array(rnorm(n * rr * tt, sd = sqrt(.1 * sqrt(v))), c(n, rr, tt)),
		ve = v,
		Sab = Sab,
		rho_mu = .6,
		rho_ab = .6,
		rho_uv = .6,
		var_mu = .1 * v,
		var_ab = .1 * v,
		var_uv = .1 * sqrt(v),
		lambda = matrix(1, n, n)
	)
	if (is.finite(prior$actor_df)) s$tau_ab <- matrix(1, 2L, tt)
	if (is.finite(prior$observation_df)) s$omega <- array(1, c(n, n, tt))
	if (D$p) names(s$beta) <- colnames(D$X)
	s
}
.lame_robust_step <- function(s, D, prior, kernel, update_hyper = TRUE) {
	n <- D$n
	tt <- D$tt
	weight <- D$observed * array(rep(s$lambda, tt), dim(D$Y))
	if (is.finite(prior$observation_df)) weight <- weight * s$omega
	comp <- .lame_robust_components(s, D)
	if (D$p) {
		s$beta <- as.numeric(
			kernel$.lame_robust_weighted_beta(D$Y - (comp$eta - comp$xb), weight, D$X, s$ve, prior$beta_precision)
		)
		comp$eta <- comp$eta - comp$xb
		comp$xb <- array(drop(D$X %*% s$beta), dim(D$Y))
		comp$eta <- comp$eta + comp$xb
	}
	s$mu <- as.numeric(kernel$.lame_robust_weighted_mu(D$Y - (comp$eta - comp$mu), weight, s$ve, s$rho_mu, s$var_mu))
	comp$eta <- comp$eta - comp$mu
	comp$mu <- array(rep(s$mu, each = n * n), dim(D$Y))
	comp$eta <- comp$eta + comp$mu
	z <- kernel$.lame_robust_weighted_actor(
		D$Y - (comp$eta - comp$a),
		weight,
		s$ve,
		s$rho_ab,
		s$var_ab,
		s$Sab,
		s$b0,
		TRUE,
		if (is.finite(prior$actor_df)) s$tau_ab[1L, ] else NULL
	)
	s$A <- z$path
	s$a0 <- as.numeric(z$level)
	comp$eta <- comp$eta - comp$a
	for (t in 1:tt) comp$a[, , t] <- s$A[, t]
	comp$eta <- comp$eta + comp$a
	z <- kernel$.lame_robust_weighted_actor(
		D$Y - (comp$eta - comp$b),
		weight,
		s$ve,
		s$rho_ab,
		s$var_ab,
		s$Sab,
		s$a0,
		FALSE,
		if (is.finite(prior$actor_df)) s$tau_ab[2L, ] else NULL
	)
	s$B <- z$path
	s$b0 <- as.numeric(z$level)
	comp$eta <- comp$eta - comp$b
	for (t in 1:tt) comp$b[, , t] <- matrix(rep(s$B[, t], each = n), n, n)
	comp$eta <- comp$eta + comp$b
	if (D$rank) {
		z <- kernel$.lame_robust_weighted_uv(D$Y - (comp$eta - comp$uv), weight, s$U, s$V, s$ve, s$rho_uv, s$var_uv)
		if (D$rank == 2L) {
			z <- .lame_uv_orientation_ar1(z, s$rho_uv, sqrt(s$var_uv), align = FALSE)
		}
		s$U <- z$U
		s$V <- z$V
		comp$eta <- comp$eta - comp$uv
		for (t in 1:tt) comp$uv[, , t] <- tcrossprod(matrix(s$U[, , t], n, D$rank), matrix(s$V[, , t], n, D$rank))
		comp$eta <- comp$eta + comp$uv
	}
	if (update_hyper) {
		s <- .lame_robust_variance_update(s, D$Y - comp$eta, D, prior)
		x <- matrix(s$mu, 1, tt)
		s$var_mu <- .lame_robust_ig(prior$ar_shape + tt / 2, prior$mu_rate + .lame_robust_ar_q(x, s$rho_mu) / 2)
		s$rho_mu <- .lame_robust_rho_draw(s$rho_mu, x, s$var_mu, 1)
		x <- rbind(s$A - s$a0, s$B - s$b0)
		d <- 2 * (n - 1)
		if (is.finite(prior$actor_df)) {
			s$var_ab <- .lame_robust_ig(
				prior$ar_shape + d * tt / 2,
				prior$ab_rate + .lame_robust_actor_q(x, s$rho_ab, s$tau_ab) / 2
			)
			s$rho_ab <- .lame_robust_actor_rho_draw(s$rho_ab, x, s$var_ab, d, s$tau_ab)
			ss <- .lame_robust_actor_ss(x, s$rho_ab)
			s$tau_ab <- matrix(
				rgamma(
					2 * tt,
					(prior$actor_df + n - 1) / 2,
					(prior$actor_df + ss / s$var_ab) / 2
				),
				2,
				tt
			)
		} else {
			s$var_ab <- .lame_robust_ig(prior$ar_shape + d * tt / 2, prior$ab_rate + .lame_robust_ar_q(x, s$rho_ab) / 2)
			s$rho_ab <- .lame_robust_rho_draw(s$rho_ab, x, s$var_ab, d)
		}
		lev <- cbind(s$a0, s$b0)
		s$Sab <- solve(rWishart(1, prior$level_df + n - 1, solve(prior$level_scale + crossprod(lev)))[, , 1])
		if (D$rank) {
			x <- rbind(matrix(s$U, n * D$rank, tt), matrix(s$V, n * D$rank, tt))
			d <- 2 * n * D$rank
			s$var_uv <- .lame_robust_ig(prior$ar_shape + d * tt / 2, prior$uv_rate + .lame_robust_ar_q(x, s$rho_uv) / 2)
			s$rho_uv <- .lame_robust_rho_draw(s$rho_uv, x, s$var_uv, d)
		}
	}
	s
}
.lame_robust_data <- function(Y, X, rank = 2L, persistent = TRUE) {
	stopifnot(is.list(Y), length(Y) > 0, all(vapply(Y, is.matrix, logical(1))))
	n <- nrow(Y[[1]])
	tt <- length(Y)
	yy <- array(unlist(Y), c(n, n, tt), dimnames = list(rownames(Y[[1]]), colnames(Y[[1]]), names(Y)))
	obs <- is.finite(yy)
	for (t in 1:tt) obs[, , t][row(yy[, , t]) == col(yy[, , t])] <- FALSE
	yy[!obs] <- 0
	p <- if (is.null(X)) 0L else dim(X[[1]])[3]
	design <- if (p) do.call(
		rbind,
		lapply(X, function(z) matrix(z, n * n, p, dimnames = list(NULL, dimnames(z)[[3]])))
	) else matrix(numeric(), n * n * tt, 0)
	design[!is.finite(design)] <- 0
	list(
		Y = yy,
		observed = obs,
		n = n,
		tt = tt,
		p = p,
		X = design,
		rank = rank,
		persistent = persistent,
		counts = apply(obs, c(1, 2), sum),
		pairs = row(yy[, , 1]) != col(yy[, , 1]),
		actor_names = rownames(Y[[1]]),
		years = names(Y)
	)
}
.lame_robust_fit <- function(
	D,
	kernel,
	seed,
	burn = 5000L,
	nscan = 20000L,
	thin = 20L,
	checkpoint = NULL,
	prior_override = list(),
	verbose = FALSE
) {
	# a NULL seed continues the caller's stream; set.seed(NULL) would reseed it
	if (!is.null(seed)) set.seed(seed)
	prior <- .lame_robust_prior(ifelse(D$observed, D$Y, NA_real_), D$p)
	prior <- utils::modifyList(prior, prior_override)
	s <- .lame_robust_initialize(D, prior)
	ns <- nscan %/% thin
	n <- D$n
	tt <- D$tt
	rr <- D$rank
	draws <- list(
		beta = matrix(NA_real_, ns, D$p),
		mu = matrix(NA_real_, ns, tt),
		A = array(NA_real_, c(n, tt, ns)),
		B = array(NA_real_, c(n, tt, ns)),
		a0 = matrix(NA_real_, n, ns),
		b0 = matrix(NA_real_, n, ns),
		U = array(NA_real_, c(n, rr, tt, ns)),
		V = array(NA_real_, c(n, rr, tt, ns)),
		lambda = array(NA_real_, c(n, n, ns)),
		hyper = matrix(
			NA_real_,
			ns,
			10,
			dimnames = list(
				NULL,
				c("ve", "rho_mu", "var_mu", "rho_ab", "var_ab", "rho_uv", "var_uv", "va", "vb", "cab")
			)
		)
	)
	if (is.finite(prior$actor_df)) draws$tau_ab <- array(
		NA_real_,
		c(2L, tt, ns),
		dimnames = list(c("sender", "receiver"), D$years, NULL)
	)
	if (is.finite(prior$observation_df)) draws$omega <- array(
		NA_real_,
		c(n, n, tt, ns),
		dimnames = list(D$actor_names, D$actor_names, D$years, NULL)
	)
	started <- Sys.time()
	for (iter in seq_len(burn + nscan)) {
		s <- .lame_robust_step(s, D, prior, kernel)
		if (iter > burn && (iter - burn) %% thin == 0L) {
			j <- (iter - burn) %/% thin
			if (is.finite(prior$actor_df)) draws$tau_ab[, , j] <- s$tau_ab
			if (is.finite(prior$observation_df)) draws$omega[, , , j] <- s$omega
			draws$beta[j, ] <- s$beta
			draws$mu[j, ] <- s$mu
			draws$A[, , j] <- s$A
			draws$B[, , j] <- s$B
			draws$a0[, j] <- s$a0
			draws$b0[, j] <- s$b0
			draws$U[, , , j] <- s$U
			draws$V[, , , j] <- s$V
			draws$lambda[, , j] <- s$lambda
			draws$hyper[j, ] <- c(
				s$ve,
				s$rho_mu,
				s$var_mu,
				s$rho_ab,
				s$var_ab,
				s$rho_uv,
				s$var_uv,
				s$Sab[1, 1],
				s$Sab[2, 2],
				s$Sab[1, 2]
			)
		}
		if (iter %% 5000L == 0L) {
			if (isTRUE(verbose)) {
				n_total <- burn + nscan
				cli::cli_text("Iteration {iter} of {n_total}")
				utils::flush.console()
			}
			if (!is.null(checkpoint)) saveRDS(
				list(
					state = s,
					iteration = iter,
					rng = .Random.seed,
					draws = draws,
					data = D,
					prior = prior,
					seed = seed
				),
				checkpoint
			)
		}
	}
	list(
		scope = "Persistent directed-pair variance model with proper priors and centered actor effects.",
		draws = draws,
		data = D,
		prior = prior,
		seed = seed,
		burn = burn,
		nscan = nscan,
		thin = thin,
		started = started,
		finished = Sys.time(),
		final_state = s
	)
}

.lame_robust_variance_update <- function(s, resid, D, prior) {
	if (is.finite(prior$observation_df %||% Inf)) {
		rr <- D$observed * resid^2
		ss <- apply(s$omega * rr, c(1, 2), sum)
		s$ve <- .lame_robust_ig(
			prior$ve_shape + sum(D$observed) / 2,
			prior$ve_rate + sum(s$lambda * ss) / 2
		)
		if (D$persistent) s$lambda[D$pairs] <- rgamma(
			sum(D$pairs),
			(prior$nu + D$counts[D$pairs]) / 2,
			(prior$nu + ss[D$pairs] / s$ve) / 2
		)
		precision <- array(rep(s$lambda, D$tt), dim(rr))
		s$omega <- array(
			rgamma(
				length(rr),
				(prior$observation_df + D$observed) / 2,
				(prior$observation_df + precision * rr / s$ve) / 2
			),
			dim(rr)
		)
		for (t in seq_len(D$tt)) s$omega[, , t][!D$pairs] <- 1
		return(s)
	}
	ss <- apply(D$observed * resid^2, c(1, 2), sum)
	s$ve <- .lame_robust_ig(prior$ve_shape + sum(D$observed) / 2, prior$ve_rate + sum(s$lambda * ss) / 2)
	if (D$persistent) s$lambda[D$pairs] <- rgamma(
		sum(D$pairs),
		shape = (prior$nu + D$counts[D$pairs]) / 2,
		rate = (prior$nu + ss[D$pairs] / s$ve) / 2
	)
	s
}
