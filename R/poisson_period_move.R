# moves of a Poisson lame() fit that carry the latent array along. on a sparse
# panel the latent slice of a period tracks its intercept, so the pair moves
# only in small conditional steps; these moves shift both together.
#   period / path translation: coefficient k of period t (or of the whole path)
#     rises by delta and the latent values by delta x_kt, leaving every residual
#     unchanged; delta is drawn by an elliptical slice step from the Gaussian path
#     prior against the log-concave Poisson factor (murray, adams and mackay 2010)
#   static coefficients of a dynamic_beta fit: the same translation under the
#     coefficient's g-prior (as in R/latent_translation_move.R)
#   scale move (.lame_ppm_scale): s2 and rho with the residuals carried along
#     their conditional law

# cells that carry a Poisson factor, period by period: positions in the period's
# matrix, counts and log exposure (`ltm` an already built context of
# R/latent_translation_move.R)
.lame_ppm_context <- function(Y, mode, offset = NULL, ltm = NULL) {
	cx <- if (is.null(ltm)) .lame_ltm_context(Y, "poisson", mode, offset = offset) else ltm
	nr <- cx$nr
	nc <- cx$nc
	nt <- cx$N
	tp <- (cx$ok_idx - 1L) %/% (nr * nc) + 1L
	pos <- cx$ok_idx - (tp - 1L) * nr * nc
	per <- lapply(seq_len(nt), function(t) {
		w <- which(tp == t)
		list(pos = pos[w], y = cx$y_ok[w], off = cx$off_ok[w])
	})
	list(nr = nr, nc = nc, N = nt, per = per, ltm = cx, cache = new.env(parent = emptyenv()))
}

# prior of the path of one coefficient as -x' P x / 2 + b' x: P the precision of
# the kind over sigma^2 lambda, b the pull of the first state toward its prior
# mean (ar1 and rw1 only)
.lame_ppm_path_prior <- function(kind, nt, rho, sigma, lam, time_positions, length_scale, m0, v0) {
	if (kind %in% c("ar1", "rw1")) {
		K <- .lame_beta_transition_precision(if (identical(kind, "rw1")) 1 else rho, diff(time_positions))
		P <- K / (sigma^2 * lam)
		P[1L, 1L] <- P[1L, 1L] + 1 / v0
		b <- numeric(nt)
		b[1L] <- m0 / v0
	} else {
		P <- .dynamic_beta_kind_to_prec(kind, nt, length_scale, time_positions) / (sigma^2 * lam)
		b <- numeric(nt)
	}
	list(P = P, b = b)
}

# design averaged over the unordered pair, for a symmetric network whose latent
# value per pair has the mean of the two directed predictors; `M` holds the
# nr x nr slices of one period as columns, `x` is an nr x nr x N array
.lame_ppm_sym_cols <- function(M, nr) {
	if (!ncol(M)) return(M)
	vapply(seq_len(ncol(M)), function(j) {
		m <- matrix(M[, j], nr, nr)
		as.numeric((m + t(m)) / 2)
	}, numeric(nrow(M)))
}
.lame_ppm_sym_array <- function(x) (x + aperm(x, c(2L, 1L, 3L))) / 2

# one visit: translations of every period and of the whole path of each dynamic
# coefficient, in random order. `path` is the N x p path, `Xdyn` the per-period
# design (nr nc x p) and `pr` the path prior list(kind, rho, sigma, lam, m0, v0,
# time_positions, length_scale). the factors exp(z + offset) are kept current as
# shifts are drawn, and Z and EZ move once at the end. a symmetric network moves
# its latent array by the pair-averaged design `Xz` and its predictor by the
# design itself
.lame_ppm_visit <- function(cx, Z, EZ, path, Xdyn, pr, whole = TRUE, Xz = Xdyn) {
	nt <- cx$N
	nrc <- cx$nr * cx$nc
	p <- ncol(path)
	stats <- c(period_calls = 0, period_evals = 0, path_calls = 0, path_evals = 0, stuck = 0)
	if (nt < 2L || p < 1L) return(list(Z = Z, EZ = EZ, path = path, stats = stats))
	ez <- lapply(seq_len(nt), function(t) exp(Z[, , t][cx$per[[t]]$pos] + cx$per[[t]]$off))
	shifted <- matrix(0, nt, p)
	for (k in sample.int(p)) {
		xo <- lapply(seq_len(nt), function(t) Xz[[t]][cx$per[[t]]$pos, k])
		if (!all(vapply(xo, function(v) all(is.finite(v)), NA))) next
		if (!all(
			vapply(seq_len(nt), function(t) all(is.finite(Xz[[t]][, k])) && all(is.finite(Xdyn[[t]][, k])), NA)
		)) next
		pri <- .lame_ppm_path_prior(
			pr$kind,
			nt,
			pr$rho[k],
			pr$sigma[k],
			pr$lam[k],
			pr$time_positions,
			pr$length_scale,
			pr$m0[k],
			pr$v0[k]
		)
		Px <- drop(pri$P %*% path[, k])
		ysx <- vapply(seq_len(nt), function(t) sum(cx$per[[t]]$y * xo[[t]]), 0)
		unit <- all(vapply(xo, function(v) all(v == 1), NA))
		shift <- function(t, d) {
			shifted[t, k] <<- shifted[t, k] + d
			ez[[t]] <<- if (unit) ez[[t]] * exp(d) else ez[[t]] * exp(d * xo[[t]])
		}
		for (t in sample.int(nt)) {
			Q <- pri$P[t, t]
			lin <- Px[t] - pri$b[t]
			if (!is.finite(Q) || Q <= 0 || !is.finite(lin)) next
			ans <- if (!identical(getOption("lame.ltm_engine"), "r")) {
				if (unit) ess_pois_scalar_cpp(ysx[t], sum(ez[[t]]), -lin / Q, 1 / sqrt(Q))
				else ess_pois_cells_cpp(ez[[t]], xo[[t]], ysx[t], -lin / Q, 1 / sqrt(Q))
			} else {
				ll <- if (unit) {
					s_t <- sum(ez[[t]])
					function(d) d * ysx[t] - s_t * expm1(d)
				} else {
					e_t <- ez[[t]]
					x_t <- xo[[t]]
					function(d) d * ysx[t] - ppm_shift_sum_cpp(e_t, x_t, d)
				}
				.lame_ess_step(0, -lin / Q, matrix(1 / sqrt(Q)), ll, 0)
			}
			stats["period_calls"] <- stats["period_calls"] + 1
			stats["period_evals"] <- stats["period_evals"] + ans$evals
			stats["stuck"] <- stats["stuck"] + ans$stuck
			d <- ans$x
			if (d == 0) next
			shift(t, d)
			path[t, k] <- path[t, k] + d
			Px <- Px + d * pri$P[, t]
		}
		if (whole) {
			Q <- sum(pri$P)
			lin <- sum(Px - pri$b)
			if (!is.finite(Q) || Q <= 0 || !is.finite(lin)) next
			ll <- if (unit) {
				s_all <- sum(vapply(ez, sum, 0))
				y_all <- sum(ysx)
				function(d) d * y_all - s_all * expm1(d)
			} else {
				y_all <- sum(ysx)
				function(d) d * y_all - sum(vapply(seq_len(nt), function(t) ppm_shift_sum_cpp(ez[[t]], xo[[t]], d), 0))
			}
			ans <- .lame_ess_step(0, -lin / Q, matrix(1 / sqrt(Q)), ll, 0)
			stats["path_calls"] <- stats["path_calls"] + 1
			stats["path_evals"] <- stats["path_evals"] + ans$evals
			stats["stuck"] <- stats["stuck"] + ans$stuck
			d <- ans$x
			if (d == 0) next
			for (t in seq_len(nt)) shift(t, d)
			path[, k] <- path[, k] + d
		}
	}
	if (any(shifted != 0)) {
		shift_of <- function(X) array(vapply(seq_len(nt), function(t) {
			moved <- which(shifted[t, ] != 0)
			if (length(moved)) drop(X[[t]][, moved, drop = FALSE] %*% shifted[t, moved]) else numeric(nrc)
		}, numeric(nrc)), dim(Z))
		sh <- shift_of(Xz)
		Z <- Z + sh
		EZ <- EZ + if (identical(Xz, Xdyn)) sh else shift_of(Xdyn)
	}
	list(Z = Z, EZ = EZ, path = path, stats = stats)
}

# dyadic variance s2 and correlation rho of the Poisson family, moved with the
# latent residuals along their approximate conditional law.
#
# with the Gaussian approximation of the Poisson factors at e = 0 (b = y - lambda,
# w = lambda) the residuals of a dyad are N(m, V), V = (Sigma^-1 + diag(w))^-1,
# m = V b. a random walk in (log s2, rho) proposes (s2', rho') and maps each dyad by
#     e' = m' + L' L^-1 (e - m)        (L L' = V, lower Cholesky factors)
# a bijection with Jacobian det(L') / det(L), so the step is an exact Metropolis
# step for any m and V. one latent value per cell or per unordered pair uses
# V = 1 / (1 / s2 + w).
# `s2_prior` c(shape0, rate0) of the s2 conditional, `two` whether rho is a
# parameter, `rho_logprior` its log prior. diagonal cells are left to the host's
# refresh unless `diag_in_target` (ame() counts them in the variance conditional)
.lame_ppm_scale <- function(
	cx,
	Z,
	EZ,
	s2,
	rho,
	s2_prior,
	two,
	rho_logprior,
	tries = 2L,
	steps = c(0.1, 0.3, 0.8),
	diag_in_target = FALSE
) {
	lx <- cx$ltm
	mode <- lx$mode
	s2_start <- s2
	rho_start <- rho
	stats <- c(scale_attempts = 0, scale_accepts = 0, scale_rho_attempts = 0, scale_rho_accepts = 0)
	lpi_s <- function(s) -(s2_prior[["shape0"]] + 1) * log(s) - s2_prior[["rate0"]] / s
	cl <- .lame_ppm_cells(cx, two)
	if (two) {
		iu <- cl$iu
		il <- cl$il
		wu <- exp(EZ[iu] + cl$offu)
		wl <- exp(EZ[il] + cl$offl)
		if (length(cl$mu)) wu[cl$mu] <- 0
		if (length(cl$ml)) wl[cl$ml] <- 0
		yu <- cl$yu
		yl <- cl$yl
		bu <- yu - wu
		bl <- yl - wl
		eu <- Z[iu] - EZ[iu]
		el <- Z[il] - EZ[il]
		st <- ppm_dyad_state_cpp(eu, el, wu, wl, yu, yl, s2, rho)
		moved <- FALSE
		for (k in seq_len(tries)) {
			tau <- steps[sample.int(length(steps), 1L)]
			type <- sample.int(3L, 1L)
			s2n <- if (type != 3L) s2 * exp(tau * stats::rnorm(1L)) else s2
			rhon <- if (type != 2L) rho + tau * stats::rnorm(1L) else rho
			if (abs(rhon) >= 1) next
			stats["scale_attempts"] <- stats["scale_attempts"] + 1
			if (type == 3L) stats["scale_rho_attempts"] <- stats["scale_rho_attempts"] + 1
			nw <- ppm_dyad_carry_cpp(eu, el, wu, wl, bu, bl, yu, yl, s2, rho, s2n, rhon)
			la <- (nw$ll - st$ll) + (nw$g - st$g) + (nw$logdet - st$logdet) +
				lpi_s(s2n) - lpi_s(s2) + log(s2n / s2) + rho_logprior(rhon) - rho_logprior(rho)
			if (is.finite(la) && log(stats::runif(1L)) < la) {
				eu <- nw$eu
				el <- nw$el
				st <- nw[c("ll", "g", "logdet")]
				s2 <- s2n
				rho <- rhon
				moved <- TRUE
				stats["scale_accepts"] <- stats["scale_accepts"] + 1
				if (type == 3L) stats["scale_rho_accepts"] <- stats["scale_rho_accepts"] + 1
			}
		}
		if (moved) {
			Z[iu] <- EZ[iu] + eu
			Z[il] <- EZ[il] + el
			# diagonal cells (no data, variance s2 (1 + rho)) scale with the variance,
			# which leaves their density times the Jacobian unchanged
			if (diag_in_target && !is.null(lx$pd)) {
				r <- sqrt(s2 * (1 + rho) / (s2_start * (1 + rho_start)))
				Z[lx$pd] <- EZ[lx$pd] + r * (Z[lx$pd] - EZ[lx$pd])
			}
		}
		return(list(Z = Z, s2 = s2, rho = rho, stats = stats))
	}
	# one latent value per cell (bipartite, directed without rho) or per unordered
	# pair (symmetric), whose mean is the mean of the two directed predictors
	if (identical(mode, "sym")) {
		mu_pair <- (EZ[lx$pu] + EZ[lx$pl]) / 2
		EZ[lx$pu] <- mu_pair
		EZ[lx$pl] <- mu_pair
	}
	idx <- cl$idx
	w <- exp(EZ[idx] + cl$off)
	if (length(cl$mb)) w[cl$mb] <- 0
	y <- cl$y
	b <- y - w
	e <- Z[idx] - EZ[idx]
	st <- ppm_cell_state_cpp(e, w, y, s2)
	moved <- FALSE
	for (k in seq_len(tries)) {
		tau <- steps[sample.int(length(steps), 1L)]
		s2n <- s2 * exp(tau * stats::rnorm(1L))
		stats["scale_attempts"] <- stats["scale_attempts"] + 1
		nw <- ppm_cell_carry_cpp(e, w, b, y, s2, s2n)
		la <- (nw$ll - st$ll) + (nw$g - st$g) + (nw$logdet - st$logdet) + lpi_s(s2n) - lpi_s(s2) + log(s2n / s2)
		if (is.finite(la) && log(stats::runif(1L)) < la) {
			e <- nw$e
			st <- nw[c("ll", "g", "logdet")]
			s2 <- s2n
			moved <- TRUE
			stats["scale_accepts"] <- stats["scale_accepts"] + 1
		}
	}
	if (moved) {
		Z[idx] <- EZ[idx] + e
		if (identical(mode, "sym")) Z[lx$pl] <- Z[lx$pu]
		if (diag_in_target && !is.null(lx$pd)) {
			r <- sqrt(s2 / s2_start)
			Z[lx$pd] <- EZ[lx$pd] + r * (Z[lx$pd] - EZ[lx$pd])
		}
	}
	list(Z = Z, s2 = s2, rho = rho, stats = stats)
}

# cells the scale move reads, found once per fit and cached: the upper and lower
# cell of each dyad when rho is a parameter, otherwise the cells that carry one
# latent value each
.lame_ppm_cells <- function(cx, two) {
	key <- if (two) "dyad" else "cell"
	got <- get0(key, envir = cx$cache, inherits = FALSE)
	if (!is.null(got)) return(got)
	lx <- cx$ltm
	obs <- lx$obs3
	Yv <- lx$Y
	off <- lx$offarr
	if (two) {
		iu <- lx$pu
		il <- lx$pl
		ou <- obs[iu]
		ol <- obs[il]
		got <- list(
			iu = iu,
			il = il,
			offu = off[iu],
			offl = off[il],
			yu = ifelse(ou, Yv[iu], 0),
			yl = ifelse(ol, Yv[il], 0),
			mu = which(!ou),
			ml = which(!ol)
		)
	} else {
		idx <- if (identical(lx$mode, "bip")) seq_along(Yv) else if (identical(lx$mode, "sym")) lx$pu else c(lx$pu, lx$pl)
		ob <- obs[idx]
		got <- list(idx = idx, off = off[idx], y = ifelse(ob, Yv[idx], 0), mb = which(!ob))
	}
	assign(key, got, envir = cx$cache)
	got
}

# static coefficients of a dynamic_beta fit with the latent array (Poisson):
# coefficient j and every latent value by delta x_j rise together, so the
# conditional of delta is its prior times the log-concave Poisson factor, drawn
# by an elliptical slice step (the translation of .lame_ltm_coef with a compiled
# factor sum). `Xz` the pair-averaged slices of a symmetric network
.lame_ppm_static <- function(cx, Z, EZ, beta, Xk, iV0, Xz = Xk) {
	lx <- cx$ltm
	stats <- c(coef_calls = 0, coef_evals = 0, coef_stuck = 0)
	for (k in sample.int(length(beta))) {
		x <- Xk[[k]]
		xz <- Xz[[k]]
		xo <- xz[lx$ok_idx]
		if (!all(is.finite(xo)) || !any(xo != 0)) next
		ez <- exp(Z[lx$ok_idx] + lx$off_ok)
		ysx <- sum(lx$y_ok * xo)
		a1 <- drop(iV0 %*% beta)[k]
		a2 <- iV0[k, k]
		if (!is.finite(a1) || !is.finite(a2) || a2 <= 0) next
		ans <- if (identical(getOption("lame.ltm_engine"), "r")) {
			ll <- function(dl) dl * ysx - ppm_shift_sum_cpp(ez, xo, dl)
			.lame_ess_step(0, -a1 / a2, matrix(1 / sqrt(a2)), ll, 0)
		} else ess_pois_cells_cpp(ez, xo, ysx, -a1 / a2, 1 / sqrt(a2))
		stats["coef_calls"] <- stats["coef_calls"] + 1
		stats["coef_evals"] <- stats["coef_evals"] + ans$evals
		stats["coef_stuck"] <- stats["coef_stuck"] + ans$stuck
		d <- ans$x
		if (d == 0) next
		Z <- Z + d * xz
		EZ <- EZ + d * x
		beta[k] <- beta[k] + d
	}
	list(Z = Z, EZ = EZ, beta = beta, stats = stats)
}
