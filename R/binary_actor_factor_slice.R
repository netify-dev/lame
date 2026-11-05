# context for the elliptical slice update of binary actor effects and factors
# binary auxiliaries are integrated out and refreshed after each update (murray, adams and mackay 2010)
.lame_binary_actor_slice_context <- function(
	rank,
	nt,
	Sab,
	margin_index,
	rho_ab,
	sigma_ab,
	rho_uv,
	sigma_uv
) {
	k <- margin_index
	o <- 3L - k
	variance <- Sab[k, k] - Sab[k, o]^2 / Sab[o, o]
	stopifnot(variance > 0, abs(rho_ab) < 1, abs(rho_uv) < 1, sigma_ab > 0, sigma_uv > 0)
	temporal <- function(rho, sigma) sigma^2 / (1 - rho^2) * rho^abs(outer(seq_len(nt), seq_len(nt), "-"))
	p <- 1L + nt + 2L * rank * nt
	ix <- list(
		level = 1L,
		dev = 1L + seq_len(nt),
		first = 1L + nt + seq_len(rank * nt),
		second = 1L + nt + rank * nt + seq_len(rank * nt)
	)
	root <- matrix(0, p, p)
	root[1, 1] <- sqrt(variance / 2)
	root[ix$dev, ix$dev] <- t(chol(temporal(rho_ab, sigma_ab) / 2))
	uv_root <- kronecker(diag(rank), t(chol(temporal(rho_uv, sigma_uv))))
	root[ix$first, ix$first] <- uv_root
	root[ix$second, ix$second] <- uv_root
	list(rank = rank, nt = nt, p = p, ix = ix, root = root, level_coefficient = Sab[k, o] / Sab[o, o])
}
.lame_binary_actor_slice_apply <- function(x, system, position, opposite, margin, level, eta) {
	cx <- system$context
	ix <- cx$ix
	r <- cx$rank
	nt <- cx$nt
	i <- system$pair[1]
	j <- system$pair[2]
	old_position <- position[c(i, j), , , drop = FALSE]
	old_margin <- margin[c(i, j), , drop = FALSE]
	level[i] <- x[ix$level]
	level[j] <- system$levelsum - level[i]
	margin[i, ] <- level[i] + x[ix$dev]
	margin[j, ] <- level[j] + system$devsum - x[ix$dev]
	position[i, , ] <- t(matrix(x[ix$first], nt, r))
	position[j, , ] <- t(matrix(x[ix$second], nt, r))
	for (tt in seq_len(nt)) {
		for (k in 1:2) {
			actor <- system$pair[k]
			delta <- position[actor, , tt] - old_position[k, , tt]
			eta[actor, , tt] <- eta[actor, , tt] + margin[actor, tt] - old_margin[k, tt] + drop(
				matrix(opposite[, , tt], dim(opposite)[1], r) %*% delta
			)
		}
	}
	list(position = position, margin = margin, level = level, eta = eta)
}
# elliptical slice update of one actor's effects and factors under the
# Z-integrated binary likelihood and their Gaussian prior (murray, adams and mackay 2010)
.lame_binary_actor_slice_system <- function(
	position,
	opposite,
	margin,
	level,
	other_level,
	Y,
	eta,
	pair,
	rho,
	context
) {
	n <- dim(position)[1]
	r <- context$rank
	nt <- context$nt
	ix <- context$ix
	p <- context$p
	i <- pair[1]
	j <- pair[2]
	stopifnot(i != j, all(pair >= 1L), all(pair <= n), abs(rho) < 1)
	levelsum <- level[i] + level[j]
	devsum <- margin[i, ] + margin[j, ] - levelsum
	pack <- function(actor) as.vector(t(matrix(position[actor, , , drop = FALSE], r, nt)))
	current <- c(level[i], margin[i, ] - level[i], pack(i), pack(j))
	prior_mean <- c(
		(levelsum + context$level_coefficient * (other_level[i] - other_level[j])) / 2,
		devsum / 2,
		rep(0, 2 * r * nt)
	)
	outside <- setdiff(seq_len(n), pair)
	from <- rep(c(rep(i, length(outside)), rep(j, length(outside)), i), nt)
	to <- rep(c(outside, outside, j), nt)
	tt <- rep(seq_len(nt), each = 2 * n - 3L)
	one <- from + (to - 1L) * n + (tt - 1L) * n * n
	two <- to + (from - 1L) * n + (tt - 1L) * n * n
	keep <- is.finite(Y[one]) | is.finite(Y[two])
	from <- from[keep]
	to <- to[keep]
	tt <- tt[keep]
	one <- one[keep]
	two <- two[keep]
	design <- function(actors, targets, times) {
		D <- matrix(0, length(actors), p)
		for (actor in pair) {
			rows <- which(actors == actor)
			if (!length(rows)) next
			sign <- if (actor == i) 1 else -1
			D[rows, ix$level] <- sign
			D[cbind(rows, ix$dev[times[rows]])] <- sign
			loc <- if (actor == i) ix$first else ix$second
			for (k in seq_len(r)) D[cbind(rows, loc[times[rows] + (k - 1L) * nt])] <- opposite[cbind(
				targets[rows],
				rep(k, length(rows)),
				times[rows]
			)]
		}
		D
	}
	D1 <- design(from, to, tt)
	D2 <- design(to, from, tt)
	list(
		context = context,
		pair = pair,
		prior_mean = prior_mean,
		current = current,
		levelsum = levelsum,
		devsum = devsum,
		D1 = D1,
		D2 = D2,
		offset1 = eta[one] - drop(D1 %*% current),
		offset2 = eta[two] - drop(D2 %*% current),
		observed1 = is.finite(Y[one]),
		observed2 = is.finite(Y[two]),
		sign1 = 2 * Y[one] - 1,
		sign2 = 2 * Y[two] - 1,
		rho = rho
	)
}
.lame_binary_actor_slice_loglik <- function(eta1, eta2, system, bvn) {
	z <- system
	both <- z$observed1 & z$observed2
	first <- z$observed1 & !z$observed2
	second <- z$observed2 & !z$observed1
	sum(bvn(z$sign1[both] * eta1[both], z$sign2[both] * eta2[both], z$rho * z$sign1[both] * z$sign2[both])) +
		sum(pnorm(z$sign1[first] * eta1[first], log.p = TRUE)) + sum(
			pnorm(z$sign2[second] * eta2[second], log.p = TRUE)
		)
}
.lame_binary_actor_slice_step <- function(system, max_abs = Inf, bvn = .lame_binary_dynamic_bvn_log) {
	z <- system
	cx <- z$context
	mu <- z$prior_mean
	center <- z$current - mu
	stopifnot(all(is.finite(z$current)), max(abs(z$current[c(cx$ix$first, cx$ix$second)])) <= max_abs)
	direction <- drop(cx$root %*% rnorm(cx$p))
	m1 <- z$offset1 + drop(z$D1 %*% mu)
	m2 <- z$offset2 + drop(z$D2 %*% mu)
	a1 <- drop(z$D1 %*% center)
	a2 <- drop(z$D2 %*% center)
	b1 <- drop(z$D1 %*% direction)
	b2 <- drop(z$D2 %*% direction)
	current_loglik <- .lame_binary_actor_slice_loglik(m1 + a1, m2 + a2, z, bvn)
	stopifnot(is.finite(current_loglik))
	threshold <- current_loglik + log(runif(1))
	theta <- runif(1, 0, 2 * pi)
	lower <- theta - 2 * pi
	upper <- theta
	evaluations <- 0L
	bounds <- 0L
	repeat{
		x <- mu + center * cos(theta) + direction * sin(theta)
		if (max(abs(x[c(cx$ix$first, cx$ix$second)])) > max_abs) {
			loglik <- -Inf
			bounds <- bounds + 1L
		} else {
			loglik <- .lame_binary_actor_slice_loglik(
				m1 + a1 * cos(theta) + b1 * sin(theta),
				m2 + a2 * cos(theta) + b2 * sin(theta),
				z,
				bvn
			)
			if (!is.finite(loglik)) stop("Nonfinite marginal elliptical-slice likelihood.")
		}
		evaluations <- evaluations + 1L
		if (loglik > threshold) break
		if (theta < 0) lower <- theta else if (theta > 0) upper <- theta else stop("Elliptical-slice bracket collapsed at current state.")
		theta <- runif(1, lower, upper)
	}
	list(x = x, evaluations = evaluations, bounds = bounds, angle = theta, loglik = loglik)
}

.lame_binary_actor_slice_round <- function(
	U,
	V,
	A,
	B,
	alevel,
	blevel,
	eta,
	Y,
	Sab,
	rho_ab,
	sigma_ab,
	rho_uv,
	sigma_uv,
	rho,
	max_abs = Inf
) {
	n <- dim(U)[1L]
	rank <- dim(U)[2L]
	nt <- dim(U)[3L]
	transpose <- function(x) aperm(x, c(2L, 1L, 3L))
	blocks <- evaluations <- bounds <- 0L
	for (side in 1:2) {
		pos <- if (side == 1) U else V
		opp <- if (side == 1) V else U
		mar <- if (side == 1) A else B
		lev <- if (side == 1) alevel else blevel
		other <- if (side == 1) blevel else alevel
		yy <- if (side == 1) Y else transpose(Y)
		ez <- if (side == 1) eta else transpose(eta)
		cx <- .lame_binary_actor_slice_context(rank, nt, Sab, side, rho_ab, sigma_ab, rho_uv, sigma_uv)
		perm <- sample.int(n)
		pairs <- matrix(perm[seq_len(2L * (n %/% 2L))], nrow = 2L)
		for (k in seq_len(ncol(pairs))) {
			ss <- .lame_binary_actor_slice_system(pos, opp, mar, lev, other, yy, ez, pairs[, k], rho, cx)
			z <- .lame_binary_actor_slice_step(ss, max_abs = max_abs)
			state <- .lame_binary_actor_slice_apply(z$x, ss, pos, opp, mar, lev, ez)
			pos <- state$position
			mar <- state$margin
			lev <- state$level
			ez <- state$eta
			blocks <- blocks + 1L
			evaluations <- evaluations + z$evaluations
			bounds <- bounds + z$bounds
		}
		if (side == 1) {
			U <- pos
			A <- mar
			alevel <- lev
			eta <- ez
		} else {
			V <- pos
			B <- mar
			blevel <- lev
			eta <- transpose(ez)
		}
	}
	list(
		U = U,
		V = V,
		A = A,
		B = B,
		alevel = alevel,
		blevel = blevel,
		eta = eta,
		blocks = blocks,
		evaluations = evaluations,
		bound_rejections = bounds
	)
}
.lame_report_actor_slice_health <- function(fit) {
	z <- fit$actor_factor_slice
	if (!is.null(z) && z["numerical_rejections"] > 0) {
		cli::cli_abort(
			c(
				"Joint actor-factor slice rounds failed numerically on {z[['numerical_rejections']]} of {z[['rounds']]} rounds.",
				"i" = "A round that fails is rejected and the chain keeps its state, so it is not a draw from the posterior and no fit is returned."
			),
			class = "lame_sampler_failure"
		)
	}
	invisible(z)
}
