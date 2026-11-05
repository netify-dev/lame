# tempered transition (Neal 1996) of one eigen-direction of a symmetric Poisson fit with dynamic positions.
# the signs of the weights L in U_t diag(L) U_t' make the posterior multimodal, and no move that keeps the
# predictor can change a weight's sign; with the latent array carried (non-centred), the conditional of
# (L_r, u_r) is its prior times the Poisson factor of the observed cells, tempered over a ladder of inverse
# temperatures and sampled by elliptical slice steps (Murray, Adams and MacKay 2010), with the down pass the
# reverse of the up pass so the tempered-transition acceptance is exact. runs every `uv_temper_every` sweeps

# transition data from .lame_ltm_context(): counts `Y` (n x n x N, mirrored, NA where unobserved),
# observed off-diagonal cells `ok` and the period log exposure `off`
.lame_temper_data <- function(cx, n, N) {
	Y <- array(NA_real_, c(n, n, N))
	Y[cx$ok_idx] <- cx$y_ok
	ok <- array(FALSE, c(n, n, N))
	ok[cx$ok_idx] <- TRUE
	for (t in seq_len(N)) {
		y <- Y[, , t]
		o <- ok[, , t]
		y[lower.tri(y)] <- t(y)[lower.tri(y)]
		o[lower.tri(o)] <- t(o)[lower.tri(o)]
		Y[, , t] <- y
		ok[, , t] <- o
	}
	list(Y = Y, ok = ok, off = cx$offset, n = n, N = N)
}

# Poisson log factor of the carried latent array z (n x n x N, symmetric)
.lame_temper_ll <- function(dt, z) {
	s <- 0
	for (t in seq_len(dt$N)) {
		o <- dt$ok[, , t] & upper.tri(dt$ok[, , t])
		zt <- z[, , t][o]
		s <- s + sum(dt$Y[, , t][o] * zt - exp(zt + dt$off[t]))
	}
	s
}

# one tempered transition of direction `r`: `U_cube` the paths (n x R x N), `L` the weights, `Z` the latent
# array (n x n x N), `dt` the data, `rho` the persistence, `l_var` the weight prior variance, `bs` the inverse
# temperatures (decreasing from 1), `cycles` passes per level, `max_abs` the coordinate bound; the compiled
# version (src/uv_tempered_move.cpp) runs unless options(lame.temper_engine = "r")
.lame_uv_sym_temper <- function(U_cube, L, Z, dt, r, rho, l_var, bs, cycles = 2L, max_abs = Inf) {
	n <- dt$n; N <- dt$N
	still <- list(U = U_cube, L = L, Z = Z, accepted = FALSE, logacc = NA_real_)
	if (!(is.finite(rho) && abs(rho) < 1 && is.finite(l_var) && l_var > 0) || !all(is.finite(Z))) return(still)
	root_n <- t(chol(rho^abs(outer(seq_len(N), seq_len(N), "-"))))
	# latent array read by the factor, upper triangle mirrored
	z <- Z
	for (t in seq_len(N)) { m <- z[, , t]; m[lower.tri(m)] <- t(m)[lower.tri(m)]; z[, , t] <- m }
	if (identical(.lame_move_engine("lame.temper_engine"), "r")) {
		o <- .lame_uv_sym_temper_r(U_cube, L, z, dt, r, root_n, l_var, bs, cycles, max_abs)
	} else {
		if (is.null(dt$Yv)) { dt$Yv <- as.numeric(dt$Y); dt$Yv[!is.finite(dt$Yv)] <- 0; dt$okv <- as.integer(dt$ok) }
		o <- uv_temper_transition_cpp(
			as.numeric(z),
			as.numeric(U_cube),
			L,
			r - 1L,
			dt$Yv,
			dt$okv,
			dt$off,
			unname(root_n),
			l_var,
			bs,
			as.integer(cycles),
			max_abs
		)
	}
	if (!is.finite(o$logacc) || !(log(stats::runif(1L)) < o$logacc)) { still$logacc <- o$logacc; return(still) }
	list(
		U = array(o$U, dim(U_cube)),
		L = o$L,
		Z = Z + (array(o$z, dim(Z)) - z),
		accepted = TRUE,
		logacc = o$logacc,
		flipped = sign(o$L[r]) != sign(L[r])
	)
}

# R reference of the transition, with the same steps, sums (long double, as in C++) and random numbers;
# `z` the mirrored latent array, `root_n` the lower Cholesky factor of a path's AR(1) correlation
.lame_uv_sym_temper_r <- function(U_cube, L, z, dt, r, root_n, l_var, bs, cycles, max_abs) {
	n <- dt$n; N <- dt$N
	K <- length(bs) - 1L
	U <- U_cube
	Lr <- L[r]
	okm <- dt$ok
	Y0 <- dt$Y
	Y0[!is.finite(Y0)] <- 0
	E <- array(0, c(n, n, N))
	for (t in seq_len(N)) E[, , t] <- ifelse(okm[, , t], exp(z[, , t] + dt$off[t]), 0)
	up <- upper.tri(matrix(0, n, n))
	ll <- sum(
		unlist(
			lapply(
				seq_len(N),
				function(t) { o <- okm[, , t] & up; (Y0[, , t] * z[, , t] - exp(z[, , t] + dt$off[t]))[o] }
			)
		)
	)
	rail <- is.finite(max_abs)
	actor <- function(a, b) {
		ui <- U[a, r, ]
		uj <- matrix(U[, r, ], n, N)
		diffv <- function(x) {
			du <- x - ui
			dz <- Lr * uj * matrix(rep(du, each = n), n, N)
			em <- expm1(dz)
			dz[a, ] <- 0
			em[a, ] <- 0
			cells <- matrix(okm[a, , ], n, N)
			cells[a, ] <- FALSE
			list(dz = dz, em = em, val = sum((matrix(Y0[a, , ], n, N) * dz - matrix(E[a, , ], n, N) * em)[cells]))
		}
		e <- stats::rnorm(N)
		nu <- vapply(seq_len(N), function(t) { m <- 0; for (q in seq_len(t)) m <- m + root_n[t, q] * e[q]; m }, 0)
		thr <- log(stats::runif(1L))
		th <- stats::runif(1L, 0, 2 * pi)
		lo <- th - 2 * pi
		hi <- th
		ev <- 0L
		repeat {
			xp <- ui * cos(th) + nu * sin(th)
			ev <- ev + 1L
			feasible <- !(rail && any(abs(xp) > max_abs))
			if (feasible) {
				d <- diffv(xp)
				if (is.finite(d$val) && b * d$val > thr) break
			}
			if (th < 0) lo <- th else hi <- th
			if (ev >= 60L) return(invisible(NULL))
			th <- stats::runif(1L, lo, hi)
		}
		for (t in seq_len(N)) {
			for (j in seq_len(n)) {
				if (j == a) next
				z[a, j, t] <<- z[a, j, t] + d$dz[j, t]
				z[j, a, t] <<- z[j, a, t] + d$dz[j, t]
				if (okm[a, j, t]) { f <- 1 + d$em[j, t]; E[a, j, t] <<- E[a, j, t] * f; E[j, a, t] <<- E[j, a, t] * f }
			}
			z[a, a, t] <<- z[a, a, t] + Lr * (xp[t] * xp[t] - ui[t] * ui[t])
			U[a, r, t] <<- xp[t]
		}
		ll <<- ll + d$val
		invisible(NULL)
	}
	weight <- function(b) {
		L0 <- Lr
		base <- array(0, c(n, n, N))
		for (t in seq_len(N)) base[, , t] <- tcrossprod(U[, r, t])
		diffv <- function(v) sum(
			unlist(
				lapply(
					seq_len(N),
					function(t) { o <- okm[, , t] & up; dz <- (v - L0) * base[, , t]; (Y0[, , t] * dz - E[, , t] * expm1(dz))[o] }
				)
			)
		)
		nuL <- sqrt(l_var) * stats::rnorm(1L)
		thr <- log(stats::runif(1L))
		th <- stats::runif(1L, 0, 2 * pi)
		lo <- th - 2 * pi
		hi <- th
		ev <- 0L
		repeat {
			xp <- L0 * cos(th) + nuL * sin(th)
			ev <- ev + 1L
			val <- diffv(xp)
			if (is.finite(val) && b * val > thr) {
				dv <- xp - L0
				for (t in seq_len(N)) {
					dz <- dv * base[, , t]
					z[, , t] <<- z[, , t] + dz
					E[, , t] <<- ifelse(okm[, , t], E[, , t] * (1 + expm1(dz)), E[, , t])
				}
				Lr <<- xp
				ll <<- ll + val
				return(invisible(NULL))
			}
			if (th < 0) lo <- th else hi <- th
			if (ev >= 60L) return(invisible(NULL))
			th <- stats::runif(1L, lo, hi)
		}
	}
	pass <- function(b, rev, with_weight) {
		steps <- c(sample.int(n), if (with_weight) 0L)
		if (rev) steps <- rev(steps)
		for (st in steps) if (st == 0L) weight(b) else actor(st, b)
	}
	logacc <- 0
	for (j in seq_len(K)) {
		logacc <- logacc + (bs[j + 1L] - bs[j]) * ll
		for (cc in seq_len(cycles)) pass(bs[j + 1L], FALSE, cc == 1L)
	}
	for (j in K:1) {
		logacc <- logacc + (bs[j] - bs[j + 1L]) * ll
		for (cc in seq_len(cycles)) pass(bs[j + 1L], TRUE, cc == cycles)
	}
	Lo <- L
	Lo[r] <- Lr
	list(z = z, U = U, L = Lo, logacc = logacc, ll = ll)
}
