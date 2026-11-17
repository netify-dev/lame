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

# one tempered transition of direction `r`: `U_cube` the paths (n x R x N), `L` the weights, `Z` the latent
# array (n x n x N), `dt` the data, `rho` the persistence, `l_var` the weight prior variance, `bs` the inverse
# temperatures (decreasing from 1), `cycles` passes per level, `max_abs` the coordinate bound; the transition
# is compiled (src/uv_tempered_move.cpp)
.lame_uv_sym_temper <- function(U_cube, L, Z, dt, r, rho, l_var, bs, cycles = 2L, max_abs = Inf) {
	n <- dt$n; N <- dt$N
	still <- list(U = U_cube, L = L, Z = Z, accepted = FALSE, logacc = NA_real_)
	if (!(is.finite(rho) && abs(rho) < 1 && is.finite(l_var) && l_var > 0) || !all(is.finite(Z))) return(still)
	root_n <- t(chol(rho^abs(outer(seq_len(N), seq_len(N), "-"))))
	# latent array read by the factor, upper triangle mirrored
	z <- Z
	for (t in seq_len(N)) { m <- z[, , t]; m[lower.tri(m)] <- t(m)[lower.tri(m)]; z[, , t] <- m }
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
