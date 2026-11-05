# joint translations of the actor-period effects and the latent array in a Poisson dynamic_ab fit,
# since effect draws given Z mix slowly on sparse panels. shifting a_it and its row of Z (a column for a
# receiver) by delta keeps every residual, so delta's conditional is the AR(1) prior times the Poisson factor
# exp(delta S1 - exp(delta) S2) of the row's observed cells. with an intercept the effects of a period move in
# random pairs (delta, -delta). elliptical slice draw (Murray, Adams and MacKay 2010); a symmetric network
# moves an actor's row and column together

# the observed counts of a panel (diagonal out, a symmetric pair read from its upper cell) and the offsets
.lame_pab_context <- function(Y, mode, offset = NULL) {
	d <- dim(Y)
	nr <- d[1L]
	nc <- d[2L]
	N <- if (length(d) >= 3L) d[3L] else 1L
	Y <- array(Y, c(nr, nc, N))
	if (!identical(mode, "bip")) for (t in seq_len(N)) diag(Y[, , t]) <- NA_real_
	if (identical(mode, "sym")) for (t in seq_len(N)) {
		y <- Y[, , t]
		up <- upper.tri(y)
		pair <- y[up]
		alt <- t(y)[up]
		pair[is.na(pair)] <- alt[is.na(pair)]
		m <- matrix(NA_real_, nr, nc)
		m[up] <- pair
		m[lower.tri(m)] <- t(m)[lower.tri(m)]
		Y[, , t] <- m
	}
	list(Y = Y, nr = nr, nc = nc, N = N, mode = mode, offset = if (is.null(offset)) numeric(N) else as.numeric(offset))
}

# elliptical slice steps from N(mu, 1/prec) for the densities exp(delta A - exp(delta) B - exp(-delta) C), one per entry (vectors)
.lame_pab_slice <- function(mu, prec, A, B, C) {
	n <- length(mu)
	lf <- function(d) d * A - exp(d) * B - exp(-d) * C
	nu <- stats::rnorm(n, mu, 1 / sqrt(prec))
	thr <- lf(rep(0, n)) + log(stats::runif(n))
	th <- stats::runif(n, 0, 2 * pi)
	lo <- th - 2 * pi
	hi <- th
	out <- numeric(n)
	todo <- rep(TRUE, n)
	evals <- 0L
	for (it in 1:200) {
		w <- which(todo)
		if (!length(w)) break
		d <- mu[w] - mu[w] * cos(th[w]) + (nu[w] - mu[w]) * sin(th[w])
		ok <- (d * A[w] - exp(d) * B[w] - exp(-d) * C[w]) > thr[w]
		evals <- evals + length(w)
		out[w[ok]] <- d[ok]
		todo[w[ok]] <- FALSE
		bad <- w[!ok]
		if (length(bad)) {
			neg <- th[bad] < 0
			lo[bad[neg]] <- th[bad[neg]]
			hi[bad[!neg]] <- th[bad[!neg]]
			th[bad] <- stats::runif(length(bad), lo[bad], hi[bad])
		}
	}
	out[todo] <- 0
	list(delta = out, evals = evals)
}

# one margin of one period: `Zt` latent slice, `y` counts (NA unobserved), `off` period offset, `eff` the
# effects with `m`, `v` their conditional prior mean and variance, `rows` TRUE for rows; returns shifts and work
.lame_pab_margin <- function(Zt, y, off, eff, m, v, rows, centred) {
	n <- length(eff)
	act <- which(is.finite(eff) & is.finite(m) & is.finite(v) & v > 0)
	shift <- numeric(n)
	if (length(act) < (if (centred) 2L else 1L)) return(list(shift = shift, pairs = 0L, evals = 0L))
	obs <- is.finite(y) & is.finite(Zt)
	y0 <- ifelse(obs, y, 0)
	e0 <- ifelse(obs, exp(Zt + off), 0)
	S1 <- if (rows) rowSums(y0) else colSums(y0)
	S2 <- if (rows) rowSums(e0) else colSums(e0)
	if (centred) {
		p <- act[sample.int(length(act))]
		np <- length(p) %/% 2L
		i <- p[seq_len(np)]
		k <- p[np + seq_len(np)]
		A <- S1[i] - S1[k]
		B <- S2[i]
		C <- S2[k]
		prec <- 1 / v[i] + 1 / v[k]
		mu <- ((m[i] - eff[i]) / v[i] + (eff[k] - m[k]) / v[k]) / prec
		sl <- .lame_pab_slice(mu, prec, A, B, C)
		shift[i] <- sl$delta
		shift[k] <- -sl$delta
		return(list(shift = shift, pairs = np, evals = sl$evals, i = i, k = k, delta = sl$delta))
	}
	i <- act
	sl <- .lame_pab_slice(m[i] - eff[i], 1 / v[i], S1[i], S2[i], numeric(length(i)))
	shift[i] <- sl$delta
	list(shift = shift, pairs = length(i), evals = sl$evals)
}

# symmetric network: an actor's cells are its row and column, so actors move one pair (or one actor) at
# a time from the current latent array; the cell shared by a pair does not move
.lame_pab_sym_period <- function(Zt, Et, y, off, eff, m, v, centred) {
	n <- length(eff)
	act <- which(is.finite(eff) & is.finite(m) & is.finite(v) & v > 0)
	pairs <- 0L
	evals <- 0L
	if (length(act) < (if (centred) 2L else 1L)) return(list(Z = Zt, E = Et, eff = eff, pairs = 0L, evals = 0L))
	# the compiled loop takes the order of the actors from here and consumes the random numbers in the order of the loop below
	if (!identical(getOption("lame.pab_engine"), "r")) {
		return(pab_sym_period_cpp(Zt, Et, y, off, eff, m, v, centred, as.integer(act[sample.int(length(act))])))
	}
	obs <- is.finite(y)
	y0 <- ifelse(obs, y, 0)
	apply_shift <- function(i, d) {
		Zt[i, ] <<- Zt[i, ] + d
		Zt[, i] <<- Zt[, i] + d
		Et[i, ] <<- Et[i, ] + d
		Et[, i] <<- Et[, i] + d
	}
	stats_of <- function(i, k) {
		j <- setdiff(seq_len(n), c(i, k))
		o <- obs[i, j] & is.finite(Zt[i, j])
		c(sum(y0[i, j][o]), sum(exp(Zt[i, j][o] + off)))
	}
	if (centred) {
		p <- act[sample.int(length(act))]
		np <- length(p) %/% 2L
		for (q in seq_len(np)) {
			i <- p[q]
			k <- p[np + q]
			si <- stats_of(i, k)
			sk <- stats_of(k, i)
			prec <- 1 / v[i] + 1 / v[k]
			mu <- ((m[i] - eff[i]) / v[i] + (eff[k] - m[k]) / v[k]) / prec
			sl <- .lame_pab_slice(mu, prec, si[1L] - sk[1L], si[2L], sk[2L])
			d <- sl$delta
			if (d != 0) {
				eff[i] <- eff[i] + d
				eff[k] <- eff[k] - d
				apply_shift(i, d)
				apply_shift(k, -d)
			}
			pairs <- pairs + 1L
			evals <- evals + sl$evals
		}
	} else {
		for (i in act[sample.int(length(act))]) {
			si <- stats_of(i, 0L)
			sl <- .lame_pab_slice(m[i] - eff[i], 1 / v[i], si[1L], si[2L], 0)
			d <- sl$delta
			eff[i] <- eff[i] + d
			apply_shift(i, d)
			pairs <- pairs + 1L
			evals <- evals + sl$evals
		}
	}
	list(Z = Zt, E = Et, eff = eff, pairs = pairs, evals = evals)
}

# one visit: `cx` from .lame_pab_context, `Z`, `EZ` latent array and predictor, `a_mat`, `b_mat` the effects
# (NA for an absent actor-period), `a_lev`, `b_lev` levels, `rho`, `sigma` AR(1) persistence and innovation
# scale, `rv`, `cv` margins in the model, `centred` effects sum to zero per period. symmetric: b_mat = a_mat
.lame_pab_visit <- function(cx, Z, EZ, a_mat, b_mat, a_lev, b_lev, rho, sigma, rv, cv, centred) {
	nr <- cx$nr
	nc <- cx$nc
	N <- cx$N
	sym <- identical(cx$mode, "sym")
	pairs <- 0L
	evals <- 0L
	for (t in seq_len(N)) {
		Zt <- Z[, , t]
		Et <- EZ[, , t]
		if (rv && sym) {
			cond <- .lame_ab_ar1_cond(a_mat, t, N, rho, sigma, level = a_lev)
			r <- .lame_pab_sym_period(
				Zt,
				Et,
				cx$Y[, , t],
				cx$offset[t],
				a_mat[, t],
				cond$m,
				rep_len(cond$v, nr),
				centred
			)
			Zt <- r$Z
			Et <- r$E
			a_mat[, t] <- r$eff
			pairs <- pairs + r$pairs
			evals <- evals + r$evals
		} else if (rv) {
			cond <- .lame_ab_ar1_cond(a_mat, t, N, rho, sigma, level = a_lev)
			r <- .lame_pab_margin(Zt, cx$Y[, , t], cx$offset[t], a_mat[, t], cond$m, rep_len(cond$v, nr), TRUE, centred)
			pairs <- pairs + r$pairs
			evals <- evals + r$evals
			if (any(r$shift != 0)) {
				a_mat[, t] <- a_mat[, t] + r$shift
				sh <- r$shift
				sh[!is.finite(sh)] <- 0
				Zt <- Zt + outer(sh, rep(1, nc))
				Et <- Et + outer(sh, rep(1, nc))
			}
		}
		if (cv && !sym) {
			cond <- .lame_ab_ar1_cond(b_mat, t, N, rho, sigma, level = b_lev)
			r <- .lame_pab_margin(
				Zt,
				cx$Y[, , t],
				cx$offset[t],
				b_mat[, t],
				cond$m,
				rep_len(cond$v, nc),
				FALSE,
				centred
			)
			pairs <- pairs + r$pairs
			evals <- evals + r$evals
			if (any(r$shift != 0)) {
				b_mat[, t] <- b_mat[, t] + r$shift
				sh <- r$shift
				sh[!is.finite(sh)] <- 0
				Zt <- Zt + outer(rep(1, nr), sh)
				Et <- Et + outer(rep(1, nr), sh)
			}
		}
		Z[, , t] <- Zt
		EZ[, , t] <- Et
	}
	if (sym) b_mat <- a_mat
	list(Z = Z, EZ = EZ, a_mat = a_mat, b_mat = b_mat, stats = c(pairs = pairs, evals = evals))
}
