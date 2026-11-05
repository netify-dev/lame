# Gaussian heat bath along affine translations that leave the likelihood
# unchanged: U shifts with V held fixed while the receiver deviations and the
# intercept compensate (roles reversed for the second update)
.lame_binary_dynamic_translation_basis <- function(V) {
	rank <- dim(V)[2]
	nt <- dim(V)[3]
	mu <- apply(V, c(2, 3), mean)
	dim(mu) <- c(rank, nt)
	C <- matrix(0, nt, rank * nt + 1L)
	for (t in seq_len(nt)) C[t, (t - 1L) * rank + seq_len(rank)] <- mu[, t]
	C[, ncol(C)] <- -1
	z <- svd(C, nu = 0, nv = ncol(C))
	tolerance <- max(dim(C)) * .Machine$double.eps * max(z$d)
	keep <- sum(z$d > tolerance)
	z$v[, seq.int(keep + 1L, ncol(C)), drop = FALSE]
}

.lame_binary_dynamic_translation_whiten <- function(x, rho, sigma) {
	nt <- ncol(x)
	z <- x
	z[, 1] <- sqrt(1 - rho^2) * x[, 1]
	if (nt > 1L)
		z[, -1] <- x[, -1, drop = FALSE] - rho * x[, -nt, drop = FALSE]
	z / sigma
}

.lame_binary_dynamic_translation_apply <- function(U, V, margin, beta, basis, q) {
	n <- dim(U)[1]
	rank <- dim(U)[2]
	nt <- dim(U)[3]
	shift <- drop(basis %*% q)
	ct <- matrix(shift[seq_len(rank * nt)], rank, nt)
	for (t in seq_len(nt)) {
		U[, , t] <- sweep(matrix(U[, , t], n, rank), 2, ct[, t], "+")
		vt <- matrix(V[, , t], n, rank)
		margin[, t] <- margin[, t] - drop(sweep(vt, 2, colMeans(vt), "-") %*% ct[, t])
	}
	beta[1] <- beta[1] - shift[length(shift)]
	list(U = U, margin = margin, beta = beta)
}

# translation step with the affected actor levels integrated out, then redrawn;
# given the opposite levels the margin has covariance K + v 11'
.lame_binary_dynamic_translation_collapsed_system <- function(
	U,
	V,
	margin,
	other_level,
	beta,
	iV0,
	rho_uv,
	sigma_uv,
	rho_ab,
	sigma_ab,
	Sab,
	margin_index
) {
	n <- dim(U)[1]
	rank <- dim(U)[2]
	nt <- dim(U)[3]
	k <- margin_index
	o <- 3L - k
	variance <- Sab[k, k] - Sab[k, o]^2 / Sab[o, o]
	mean <- Sab[k, o] / Sab[o, o] * other_level
	stopifnot(is.finite(variance), variance > 0)
	K <- sigma_ab^2 / (1 - rho_ab^2) * rho_ab^abs(outer(seq_len(nt), seq_len(nt), "-"))
	C <- K + variance
	white <- solve(chol(C))
	basis <- .lame_binary_dynamic_translation_basis(V)
	d <- ncol(basis)
	old <- c(
		.lame_binary_dynamic_translation_whiten(matrix(U, ncol = nt), rho_uv, sigma_uv),
		(margin - mean) %*% white
	)
	directions <- matrix(0, length(old), d)
	du <- array(0, c(n, rank, nt))
	dm <- matrix(0, n, nt)
	for (j in seq_len(d)) {
		ct <- matrix(basis[seq_len(rank * nt), j], rank, nt)
		for (t in seq_len(nt)) {
			du[, , t] <- matrix(rep(ct[, t], each = n), n, rank)
			vt <- matrix(V[, , t], n, rank)
			dm[, t] <- -drop(sweep(vt, 2, colMeans(vt), "-") %*% ct[, t])
		}
		directions[, j] <- c(
			.lame_binary_dynamic_translation_whiten(
				matrix(du, ncol = nt),
				rho_uv,
				sigma_uv
			),
			dm %*% white
		)
	}
	delta <- basis[nrow(basis), ]
	list(
		basis = basis,
		precision = crossprod(directions) + iV0[1, 1] * tcrossprod(delta),
		score = drop(crossprod(directions, old)) - delta * sum(iV0[1, ] * beta),
		level_variance = variance,
		level_mean = mean,
		K = K,
		C = C
	)
}

.lame_binary_dynamic_translation_collapsed_step <- function(
	U,
	V,
	margin,
	level,
	other_level,
	beta,
	iV0,
	rho_uv,
	sigma_uv,
	rho_ab,
	sigma_ab,
	Sab,
	margin_index,
	max_abs = Inf
) {
	z <- .lame_binary_dynamic_translation_collapsed_system(
		U,
		V,
		margin,
		other_level,
		beta,
		iV0,
		rho_uv,
		sigma_uv,
		rho_ab,
		sigma_ab,
		Sab,
		margin_index
	)
	L <- chol(z$precision)
	q <- -backsolve(L, forwardsolve(t(L), z$score)) + backsolve(L, rnorm(ncol(L)))
	ans <- .lame_binary_dynamic_translation_apply(U, V, margin, beta, z$basis, q)
	rejected <- max(abs(ans$U)) > max_abs
	if (rejected)
		ans <- list(U = U, margin = margin, beta = beta)
	ki <- drop(solve(z$K, rep(1, ncol(margin))))
	variance <- 1 / (1 / z$level_variance + sum(ki))
	mean <- variance * (z$level_mean / z$level_variance + drop(ans$margin %*% ki))
	draw <- mean + sqrt(variance) * rnorm(nrow(margin))
	ans$level <- draw - mean(draw)
	ans$bound_rejected <- rejected
	ans
}

# non-centered scale move on the factors and their innovation sd together; the
# ratio keeps the likelihood and the inverse-gamma variance prior
.lame_binary_dynamic_uv_scale_step <- function(
	eta,
	uv,
	offset,
	sign,
	sigma,
	shape,
	scale,
	variance_max = Inf,
	factor_max = Inf,
	step = 0.08,
	tries = 8L
) {
	factor <- 1
	accepted <- 0L
	loglik <- function(f) sum(pnorm(sign * (eta + (f * f - 1) * uv + offset), log.p = TRUE))
	old <- loglik(1)
	v <- sigma^2
	for (k in seq_len(tries)) {
		proposal <- factor * exp(rnorm(1, 0, step))
		vn <- v * proposal^2
		if (!is.finite(proposal) || proposal <= 0 || proposal > factor_max || vn > variance_max)
			next
		newer <- loglik(proposal)
		ratio <- newer - old - shape * log(vn / (v * factor^2)) - scale / vn + scale / (v * factor^2)
		if (is.finite(ratio) && log(runif(1)) < ratio) {
			factor <- proposal
			old <- newer
			accepted <- accepted + 1L
		}
	}
	list(factor = factor, sigma = sigma * factor, accepts = accepted, attempts = tries)
}

# partially collapsed reciprocity update: one standard normal factor per dyad
# carries the pair correlation, the independent residuals are integrated out
.lame_binary_dynamic_rho_factor_draw <- function(Z, EZ, rho) {
	n <- dim(Z)[1]
	nt <- dim(Z)[3]
	upper <- upper.tri(matrix(0, n, n))
	raw <- array(0, dim(Z))
	a <- abs(rho)
	sgn <- if (rho < 0)
		-1
	else 1
	for (t in seq_len(nt)) {
		e <- Z[, , t] - EZ[, , t]
		mu <- sqrt(a) * (e + sgn * t(e))[upper] / (1 + a)
		f <- matrix(0, n, n)
		f[upper] <- rnorm(length(mu), mu, sqrt((1 - a) / (1 + a)))
		raw[, , t] <- f + t(f)
	}
	raw
}

.lame_binary_dynamic_rho_offset <- function(raw, rho) {
	out <- sqrt(abs(rho)) * raw
	if (rho < 0)
		for (t in seq_len(dim(out)[3])) {
			x <- out[, , t]
			x[lower.tri(x)] <- -x[lower.tri(x)]
			out[, , t] <- x
		}
	out
}

.lame_binary_dynamic_rho_z_draw <- function(Z, EZ, off, rho, idx1, idx0, idx_na, idx_diag) {
	stopifnot(is.finite(rho), abs(rho) < 1)
	sd_cell <- sqrt(1 - abs(rho))
	z <- Z
	if (length(idx1)) {
		mu <- EZ[idx1] + off[idx1]
		lt <- pnorm(-mu / sd_cell, lower.tail = FALSE, log.p = TRUE)
		z[idx1] <- mu + sd_cell * qnorm(
			log(runif(length(mu))) + lt,
			lower.tail = FALSE,
			log.p = TRUE
		)
	}
	if (length(idx0)) {
		mu <- EZ[idx0] + off[idx0]
		lt <- pnorm(-mu / sd_cell, log.p = TRUE)
		z[idx0] <- mu + sd_cell * qnorm(log(runif(length(mu))) + lt, log.p = TRUE)
	}
	if (length(idx_na))
		z[idx_na] <- rnorm(length(idx_na), EZ[idx_na] + off[idx_na], sd_cell)
	if (length(idx_diag))
		z[idx_diag] <- rnorm(length(idx_diag), EZ[idx_diag], sqrt(1 + rho))
	if (any(!is.finite(z)))
		stop("Nonfinite conditional binary response in the non-centred AR(1) move.")
	z
}

# HMC on the binary pair likelihood with both Gaussian auxiliaries integrated
# out and whitened priors; the caller redraws z after an accepted trajectory
.lame_binary_dynamic_pair_hmc_context <- function(Y, Xlist, rank) {
	n <- dim(Y)[1]
	nt <- dim(Y)[3]
	k <- n - 1L
	p <- dim(Xlist[[1]])[3]
	obs <- is.finite(Y)
	for (t in seq_len(nt)) diag(obs[, , t]) <- FALSE
	dims <- list(
		beta = p,
		level = c(k, 2L),
		adev = c(k, nt),
		bdev = c(k, nt),
		U = c(
			n *
				rank,
			nt
		),
		V = c(n * rank, nt)
	)
	sizes <- vapply(dims, prod, numeric(1))
	ends <- cumsum(sizes)
	indices <- lapply(seq_along(sizes), function(i) seq.int(ends[i] - sizes[i] + 1, ends[i]))
	names(indices) <- names(dims)
	rate <- max(0.001, min(0.999, mean(Y[obs])))
	list(
		n = n,
		nt = nt,
		k = k,
		p = p,
		rank = rank,
		pairs = .lame_binary_dynamic_collapsed_rho_context(Y),
		H = qr.Q(qr(stats::contr.helmert(n))),
		X = lapply(
			Xlist,
			function(x) matrix(
				x,
				n *
					n,
				p
			)
		),
		obs = obs,
		sign = 2 * Y[obs] - 1,
		dims = dims,
		indices = indices,
		size = sum(sizes),
		weight = dnorm(qnorm(rate))^2 / (rate * (1 - rate))
	)
}

.lame_binary_dynamic_pair_hmc_hyper <- function(
	cx,
	iV0,
	Sab,
	rho_ab,
	sigma_ab,
	rho_uv,
	sigma_uv,
	rho,
	max_abs = Inf
) {
	temporal <- function(r, s) chol(
		s^2 / (1 - r^2) * r^abs(outer(
			seq_len(cx$nt),
			seq_len(cx$nt),
			"-"
		))
	)
	list(
		beta_root = solve(chol(iV0)),
		level_root = chol(Sab),
		ab_root = temporal(
			rho_ab,
			sigma_ab
		),
		uv_root = temporal(rho_uv, sigma_uv),
		uv_variance = sigma_uv^2 / (1 - rho_uv^2),
		rho = rho,
		sd = 1,
		max_abs = max_abs
	)
}

.lame_binary_dynamic_pair_hmc_unpack <- function(q, cx) {
	z <- lapply(seq_along(cx$dims), function(j) {
		d <- cx$dims[[j]]
		v <- q[cx$indices[[j]]]
		if (length(d) > 1)
			dim(v) <- d
		v
	})
	names(z) <- names(cx$dims)
	z
}

.lame_binary_dynamic_pair_hmc_pack <- function(z) unlist(z, use.names = FALSE)

.lame_binary_dynamic_pair_hmc_decode <- function(q, cx, h) {
	z <- .lame_binary_dynamic_pair_hmc_unpack(q, cx)
	list(
		beta = drop(h$beta_root %*% z$beta),
		level = cx$H %*% (z$level %*% h$level_root),
		adev = cx$H %*% (z$adev %*% h$ab_root),
		bdev = cx$H %*% (z$bdev %*% h$ab_root),
		U = array(
			z$U %*%
				h$uv_root,
			c(cx$n, cx$rank, cx$nt)
		),
		V = array(
			z$V %*% h$uv_root,
			c(
				cx$n,
				cx$rank,
				cx$nt
			)
		)
	)
}

.lame_binary_dynamic_pair_hmc_encode <- function(beta, level, adev, bdev, U, V, cx, h) {
	.lame_binary_dynamic_pair_hmc_pack(list(
		solve(h$beta_root, beta),
		crossprod(cx$H, level) %*%
			solve(h$level_root),
		crossprod(cx$H, adev) %*% solve(h$ab_root),
		crossprod(
			cx$H,
			bdev
		) %*% solve(h$ab_root),
		matrix(U, ncol = cx$nt) %*% solve(h$uv_root),
		matrix(
			V,
			ncol = cx$nt
		) %*% solve(h$uv_root)
	))
}

.lame_binary_dynamic_pair_hmc_evaluate <- function(q, cx, h, offset) {
	z <- .lame_binary_dynamic_pair_hmc_decode(q, cx, h)
	if (max(abs(z$U), abs(z$V)) > h$max_abs)
		return(list(logp = -Inf, gradient = NULL, bound = TRUE))
	eta <- array(0, c(cx$n, cx$n, cx$nt))
	gbeta <- numeric(cx$p)
	ga <- gb <- matrix(0, cx$n, cx$nt)
	gu <- gv <- array(0, dim(z$U))
	for (t in seq_len(cx$nt)) eta[, , t] <- matrix(cx$X[[t]] %*% z$beta, cx$n, cx$n) + outer(
		z$level[,
			1] + z$adev[, t],
		z$level[, 2] + z$bdev[, t],
		"+"
	) + tcrossprod(
		matrix(
			z$U[, , t],
			cx$n,
			cx$rank
		),
		matrix(z$V[, , t], cx$n, cx$rank)
	)
	pair <- .lame_binary_dynamic_pair_likelihood_gradient(eta, h$rho, cx$pairs)
	lp <- pair$logp
	ge <- pair$gradient
	for (t in seq_len(cx$nt)) {
		e <- ge[, , t]
		gbeta <- gbeta + drop(crossprod(cx$X[[t]], c(e)))
		ga[, t] <- rowSums(e)
		gb[, t] <- colSums(e)
		gu[, , t] <- e %*% matrix(z$V[, , t], cx$n, cx$rank)
		gv[, , t] <- t(e) %*% matrix(z$U[, , t], cx$n, cx$rank)
	}
	gradient <- .lame_binary_dynamic_pair_hmc_pack(list(
		drop(crossprod(h$beta_root, gbeta)),
		crossprod(cx$H, cbind(rowSums(ga), rowSums(gb))) %*% t(h$level_root),
		crossprod(
			cx$H,
			ga
		) %*% t(h$ab_root),
		crossprod(cx$H, gb) %*% t(h$ab_root),
		matrix(gu, ncol = cx$nt) %*%
			t(h$uv_root),
		matrix(gv, ncol = cx$nt) %*% t(h$uv_root)
	)) - q
	list(logp = sum(lp) - sum(q * q) / 2, gradient = gradient, eta = eta, bound = FALSE)
}

# diagonal mass from the data and fixed hyperparameters, not the position
.lame_binary_dynamic_pair_hmc_mass <- function(cx, h) {
	weight <- cx$weight / (h$sd * h$sd)
	xx <- Reduce(
		"+",
		lapply(
			seq_len(cx$nt),
			function(t) crossprod(
				cx$X[[t]][c(
					cx$obs[, ,
						t]
				), , drop = FALSE]
			)
		)
	)
	mb <- 1 + diag(crossprod(h$beta_root, xx %*% h$beta_root)) * weight
	ml <- matrix(
		rep(1 + cx$nt * (cx$n - 1) * weight * rowSums(h$level_root^2), each = cx$k),
		cx$k,
		2
	)
	ma <- matrix(
		rep(1 + (cx$n - 1) * weight * rowSums(h$ab_root^2), each = cx$k),
		cx$k,
		cx$nt
	)
	mu <- matrix(
		rep(
			1 + (cx$n - 1) * weight * h$uv_variance * rowSums(h$uv_root^2),
			each = cx$n *
				cx$rank
		),
		cx$n * cx$rank,
		cx$nt
	)
	.lame_binary_dynamic_pair_hmc_pack(list(mb, ml, ma, ma, mu, mu))
}

.lame_binary_dynamic_pair_likelihood_gradient <- function(eta, rho, cx) {
	both <- cx$observed1 & cx$observed2
	ge <- array(0, dim(eta))
	lp <- 0
	if (any(both)) {
		i <- cx$one[both]
		j <- cx$two[both]
		s1 <- cx$sign1[both]
		s2 <- cx$sign2[both]
		a <- s1 * eta[i]
		b <- s2 * eta[j]
		r <- rho * s1 * s2
		sd <- sqrt(1 - r * r)
		logprob <- .lame_binary_dynamic_bvn_log(a, b, r)
		lp <- sum(logprob)
		ge[i] <- s1 * exp(dnorm(a, log = TRUE) + pnorm((b - r * a) / sd, log.p = TRUE) - logprob)
		ge[j] <- s2 * exp(dnorm(b, log = TRUE) + pnorm((a - r * b) / sd, log.p = TRUE) - logprob)
	}
	single <- xor(cx$observed1, cx$observed2)
	if (any(single)) {
		i <- ifelse(cx$observed1[single], cx$one[single], cx$two[single])
		sign <- ifelse(cx$observed1[single], cx$sign1[single], cx$sign2[single])
		v <- sign * eta[i]
		logprob <- pnorm(v, log.p = TRUE)
		lp <- lp + sum(logprob)
		ge[i] <- sign * exp(dnorm(v, log = TRUE) - logprob)
	}
	list(logp = lp, gradient = ge)
}

# non-centered AR(1) persistence move: standardized innovations and stationary
# variance w stay fixed while rho and v = w * (1 - rho^2) change together; the
# inverse-gamma prior on v in (rho, w) coordinates contributes -shape * log(v)
.lame_binary_dynamic_nc_ar_innovations <- function(paths, rho, sigma) {
	nt <- ncol(paths)
	z <- paths
	z[, 1] <- paths[, 1] * sqrt(1 - rho * rho) / sigma
	if (nt > 1)
		z[, -1] <- (paths[, -1, drop = FALSE] - rho * paths[, -nt, drop = FALSE]) / sigma
	z
}

.lame_binary_dynamic_nc_ar_paths <- function(z, rho, variance) {
	x <- z
	x[, 1] <- sqrt(variance / (1 - rho * rho)) * z[, 1]
	if (ncol(z) > 1)
		for (t in 2:ncol(z)) x[, t] <- rho * x[, t - 1] + sqrt(variance) * z[, t]
	x
}

.lame_binary_dynamic_nc_ar_step <- function(
	paths,
	rho,
	sigma,
	loglikelihood,
	prior_mean,
	prior_sd,
	shape,
	scale,
	stationary_max = Inf,
	sigma_max = Inf,
	max_abs = Inf,
	step = 0.15,
	tries = 4L
) {
	limit <- 0.99
	evaluate <- function(x) {
		v <- loglikelihood(x)
		list(value = if (is.list(v)) v$loglik else v, detail = v)
	}
	ev <- evaluate(paths)
	# skip the move at the persistence bound or the stationary cap
	w <- sigma * sigma / (1 - rho * rho)
	if (!is.finite(rho) || abs(rho) >= limit || !is.finite(sigma) || sigma <= 0 ||
			!is.finite(w) || w > stationary_max * (1 + 1e-12))
		return(list(
			paths = paths,
			rho = rho,
			sigma = sigma,
			evaluation = ev$detail,
			accepts = 0L,
			attempts = 0L
		))
	innovation <- .lame_binary_dynamic_nc_ar_innovations(paths, rho, sigma)
	accepted <- 0L
	hp <- function(r, v) dnorm(r, prior_mean, prior_sd, log = TRUE) - shape * log(v) - scale / v +
		log1p(-(r / limit)^2)
	old <- ev$value + hp(rho, sigma * sigma)
	for (k in seq_len(tries)) {
		prop <- limit * tanh(atanh(rho / limit) + rnorm(1, 0, step))
		vn <- w * (1 - prop * prop)
		if (!is.finite(prop) || abs(prop) >= limit || vn > sigma_max * sigma_max)
			next
		pn <- .lame_binary_dynamic_nc_ar_paths(innovation, prop, vn)
		if (max(abs(pn)) > max_abs)
			next
		en <- evaluate(pn)
		new <- en$value + hp(prop, vn)
		if (is.finite(new) && log(runif(1)) < new - old) {
			rho <- prop
			sigma <- sqrt(vn)
			paths <- pn
			ev <- en
			old <- new
			accepted <- accepted + 1L
		}
	}
	list(
		paths = paths,
		rho = rho,
		sigma = sigma,
		evaluation = ev$detail,
		accepts = accepted,
		attempts = tries
	)
}

# asymptotic Mills ratio series for z <= -10; 32 terms bound the remainder below 8e-21
.lame_normal_mills_series <- function(z) {
	inverse_square <- (1 / z)^2
	term <- rep(1, length(z)); value <- term
	for (k in seq_len(32L)) {
		term <- -term * (2 * k - 1) * inverse_square
		value <- value + term
	}
	value
}

.lame_normal_mills_ratio <- function(z) {
	answer <- numeric(length(z)); tail <- z <= -10
	answer[tail] <- -z[tail] / .lame_normal_mills_series(z[tail])
	moderate <- !tail & z <= -5
	# Laplace continued fraction with 64 levels for the moderate tail
	x <- -z[moderate]
	denominator <- x
	for (k in 64L:1L) denominator <- x + k / denominator
	answer[moderate] <- denominator
	ordinary <- !tail & !moderate
	answer[ordinary] <- exp(dnorm(z[ordinary], log = TRUE) - pnorm(z[ordinary], log.p = TRUE))
	answer
}

.lame_normal_log_cdf_shift <- function(z, displacement) {
	z <- rep(z, length.out = length(displacement)); shifted <- z + displacement
	answer <- numeric(length(displacement)); tail <- z <= -10 & shifted <= -10
	# expand the shift in the far tail rather than subtracting two log tails
	answer[tail] <- -displacement[tail] * (z[tail] + displacement[tail] / 2) -
		log1p(displacement[tail] / z[tail]) +
		log(.lame_normal_mills_series(shifted[tail])) - log(.lame_normal_mills_series(z[tail]))
	answer[!tail] <- pnorm(shifted[!tail], log.p = TRUE) - pnorm(z[!tail], log.p = TRUE)
	answer
}

.lame_binary_dynamic_bvn_log_quad <- function(a, b, rho) {
	if (rho == 0)
		return(pnorm(a, log.p = TRUE) + pnorm(b, log.p = TRUE))
	pa <- pnorm(a, log.p = TRUE)
	pb <- pnorm(b, log.p = TRUE)
	qa <- pnorm(a, lower.tail = FALSE, log.p = TRUE)
	qb <- pnorm(b, lower.tail = FALSE, log.p = TRUE)
	# return a marginal when the other tail is below machine precision
	if (qa - pb < log(4 * .Machine$double.eps)) return(pb)
	if (qb - pa < log(4 * .Machine$double.eps)) return(pa)
	sd <- sqrt((1 - rho) * (1 + rho))
	logf <- function(x) dnorm(x, log = TRUE) + pnorm((b - rho * x) / sd, log.p = TRUE)
	derivative <- function(x) {
		z <- (b - rho * x) / sd
		-x - rho / sd * .lame_normal_mills_ratio(z)
	}
	da <- derivative(a)
	if (da >= 0) mode <- a else {
		low <- min(a - 10, -abs(b) - 10, -10)
		while (derivative(low) < 0) low <- 2 * low - 10
		mode <- stats::uniroot(derivative, c(low, a), tol = 1e-11)$root
	}
	scale <- if (da >= 0) min(sd, 1 / max(1, da)) else sd
	peak <- logf(mode)
	conditional_mode <- (b - rho * mode) / sd
	logg <- function(u) {
		displacement <- scale * u
		-displacement * (mode + displacement / 2) +
			.lame_normal_log_cdf_shift(conditional_mode, (-rho / sd) * displacement)
	}
	left <- -1
	while (logg(left) > -45) left <- 2 * left
	upper <- max(0, (a - mode) / scale)
	right <- min(upper, 1)
	while (right < upper && logg(right) > -45) right <- min(upper, 2 * right)
	# integrate piecewise over geometric intervals from the peak, so both a
	# narrow transition and a broad tail near |rho| = 1 are resolved
	piecewise <- function(edge) {
		if (edge == 0) return(0)
		cuts <- sort(unique(c(
			0,
			edge,
			sign(edge) * 2^seq.int(0, max(0, ceiling(log2(abs(edge)))))
		)))
		cuts <- cuts[cuts >= min(0, edge) & cuts <= max(0, edge)]
		sum(vapply(seq_len(length(cuts) - 1L), function(k) {
			stats::integrate(
				function(u) exp(logg(u)),
				cuts[k],
				cuts[k + 1L],
				rel.tol = 1e-11,
				abs.tol = 1e-13,
				subdivisions = 300L
			)$value
		}, numeric(1)))
	}
	ans <- piecewise(left) + piecewise(right)
	stopifnot(is.finite(ans), ans > 0)
	answer <- peak + log(scale) + log(ans)
	stopifnot(answer <= 1e-8)
	min(answer, 0)
}

.lame_binary_dynamic_bvn_log <- function(a, b, rho) {
	rho <- rep(rho, length.out = length(a))
	stopifnot(all(is.finite(a)), all(is.finite(b)), all(abs(rho) < 1))
	if (!length(a))
		return(numeric())
	p <- pbivnorm::pbivnorm(a, b, rho)
	ans <- rep(NA_real_, length(p))
	good <- is.finite(p) & p > 1e-10 & p <= 1
	ans[good] <- log(p[good])
	bad <- which(!good)
	if (length(bad))
		for (k in bad) ans[k] <- .lame_binary_dynamic_bvn_log_quad(a[k], b[k], rho[k])
	ans
}

.lame_binary_dynamic_collapsed_rho_context <- function(Y) {
	n <- dim(Y)[1]
	nt <- dim(Y)[3]
	ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	one <- unlist(lapply(
		seq_len(nt),
		function(t) ij[, 1] + (ij[, 2] - 1L) * n + (t - 1L) *
			n * n
	))
	two <- unlist(lapply(
		seq_len(nt),
		function(t) ij[, 2] + (ij[, 1] - 1L) * n + (t - 1L) *
			n * n
	))
	list(
		one = one,
		two = two,
		observed1 = is.finite(Y[one]),
		observed2 = is.finite(Y[two]),
		sign1 = 2 * Y[one] - 1,
		sign2 = 2 * Y[two] - 1,
		diagonal = unlist(lapply(
			seq_len(nt),
			function(t) seq_len(n) + (seq_len(n) - 1L) * n + (t - 1L) * n * n
		))
	)
}

.lame_binary_dynamic_collapsed_rho_loglik <- function(rho, eta, cx) {
	both <- cx$observed1 & cx$observed2
	sum(.lame_binary_dynamic_bvn_log(
		cx$sign1[both] * eta[cx$one[both]],
		cx$sign2[both] *
			eta[cx$two[both]],
		rho * cx$sign1[both] * cx$sign2[both]
	))
}

# independence proposals from the prior of rho, so prior and proposal cancel
.lame_binary_dynamic_collapsed_rho_mh <- function(rho, eta, cx, tries = 4L) {
	old <- .lame_binary_dynamic_collapsed_rho_loglik(rho, eta, cx)
	stopifnot(is.finite(old))
	acc <- 0L
	for (k in seq_len(tries)) {
		proposal <- .lame_rho_prior_draw()
		new <- .lame_binary_dynamic_collapsed_rho_loglik(proposal, eta, cx)
		if (is.finite(new) && log(runif(1)) < new - old) {
			rho <- proposal
			old <- new
			acc <- acc + 1L
		}
	}
	list(rho = rho, accepts = acc, attempts = tries)
}

# exact draw of a bivariate normal pair above two thresholds: exponential
# tilting of the first truncated marginal with a tangent envelope from the
# log-concave second survival function, anchored at the first marginal's mode
.lame_binary_dynamic_biv_positive <- function(lower1, lower2, rho) {
	n <- length(lower1)
	r <- rep(rho, length.out = n)
	swap <- lower2 > lower1
	a <- ifelse(swap, lower2, lower1)
	b <- ifelse(swap, lower1, lower2)
	sd <- sqrt(1 - r * r)
	anchor <- pmax(a, 0)
	for (iter in seq_len(40L)) {
		zz <- (r * anchor - b) / sd
		mills <- exp(dnorm(zz, log = TRUE) - pnorm(zz, log.p = TRUE))
		derivative <- -anchor + r / sd * mills
		curvature <- -1 - (r / sd)^2 * pmax(mills * (zz + mills), 0)
		next_anchor <- pmax(a, anchor - derivative / curvature)
		difference <- max(abs(next_anchor - anchor))
		anchor <- next_anchor
		if (difference < 1e-10)
			break
	}
	z <- (r * anchor - b) / sd
	tilt <- r / sd * exp(dnorm(z, log = TRUE) - pnorm(z, log.p = TRUE))
	envelope <- pnorm(z, log.p = TRUE) - tilt * anchor + 1e-10
	tail <- pnorm(a - tilt, lower.tail = FALSE, log.p = TRUE)
	x <- y <- numeric(n)
	pending <- seq_len(n)
	iterations <- 0L
	proposals <- 0L
	while (length(pending)) {
		iterations <- iterations + 1L
		if (iterations > 100000L)
			stop("Truncated-pair rejection failed to finish")
		k <- pending
		xx <- tilt[k] + qnorm(log(runif(length(k))) + tail[k], lower.tail = FALSE, log.p = TRUE)
		log_accept <- pnorm((r[k] * xx - b[k]) / sd[k], log.p = TRUE) - tilt[k] * xx - envelope[k]
		if (any(!is.finite(xx)) || any(log_accept > 1e-07))
			stop("Invalid truncated-pair envelope")
		accept <- log(runif(length(k))) < log_accept
		proposals <- proposals + length(k)
		if (any(accept)) {
			kk <- k[accept]
			xx <- xx[accept]
			x[kk] <- xx
			yy <- r[kk] * xx + sd[kk] * qnorm(
				log(runif(length(kk))) + pnorm(
					(b[kk] - r[kk] *
							xx) / sd[kk],
					lower.tail = FALSE,
					log.p = TRUE
				),
				lower.tail = FALSE,
				log.p = TRUE
			)
			y[kk] <- yy
		}
		pending <- k[!accept]
	}
	ans <- cbind(ifelse(swap, y, x), ifelse(swap, x, y))
	stopifnot(
		all(is.finite(ans)),
		all(ans[, 1] >= lower1 - 1e-10),
		all(
			ans[, 2] >= lower2 -
				1e-10
		)
	)
	attr(ans, "proposals") <- proposals
	attr(ans, "iterations") <- iterations
	ans
}

.lame_binary_dynamic_collapsed_rho_z <- function(eta, Y, rho, cx = .lame_binary_dynamic_collapsed_rho_context(Y)) {
	z <- eta
	i <- cx$one
	j <- cx$two
	both <- cx$observed1 & cx$observed2
	none <- !cx$observed1 & !cx$observed2
	if (any(both)) {
		k <- which(both)
		s1 <- cx$sign1[k]
		s2 <- cx$sign2[k]
		eps <- .lame_binary_dynamic_biv_positive(
			-s1 * eta[i[k]],
			-s2 * eta[j[k]],
			rho *
				s1 * s2
		)
		z[i[k]] <- eta[i[k]] + s1 * eps[, 1]
		z[j[k]] <- eta[j[k]] + s2 * eps[, 2]
	}
	if (any(none)) {
		k <- which(none)
		e <- rnorm(length(k))
		z[i[k]] <- eta[i[k]] + e
		z[j[k]] <- eta[j[k]] + rho * e + sqrt(1 - rho * rho) * rnorm(length(k))
	}
	single <- which(xor(cx$observed1, cx$observed2))
	if (length(single)) {
		observed <- ifelse(cx$observed1[single], i[single], j[single])
		hidden <- ifelse(cx$observed1[single], j[single], i[single])
		sign <- 2 * Y[observed] - 1
		eps <- qnorm(
			log(runif(length(single))) + pnorm(sign * eta[observed], log.p = TRUE),
			lower.tail = FALSE,
			log.p = TRUE
		)
		z[observed] <- eta[observed] + sign * eps
		z[hidden] <- eta[hidden] + rho * sign * eps + sqrt(1 - rho * rho) * rnorm(length(single))
	}
	z[cx$diagonal] <- rnorm(length(cx$diagonal), eta[cx$diagonal], sqrt(1 + rho))
	stopifnot(all(is.finite(z)))
	z
}
