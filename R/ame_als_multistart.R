# ame_als_multistart.r
#
# multi-start optimisation for the fast ame estimator. for r > 0 the
# multiplicative low-rank objective is non-convex, so block coordinate descent
# can land in different local optima depending on the starting factors.
# `multistart = "cheap"` / "full" runs the bcd from several starts -- the
# deterministic cold start plus random low-rank starts -- and keeps the
# lowest-sse fit, with a diagnostic when the starts disagree.
#
# the random starts use `seed` for reproducibility and the global rng state is
# restored afterwards, so the estimator stays deterministic given `seed` and
# does not mutate the caller's rng stream.

.ae_multistart_fit <- function(Z, X, R, symmetric, max_iter, tol,
                               lowrank_method, linear_solver,
                               n_starts, seed) {

	bcd <- function(init) {
		.ae_bcd_fit(Z, X, R, symmetric, max_iter, tol, init = init,
		            lowrank_method = lowrank_method,
		            linear_solver = linear_solver)
	}

	fits <- vector("list", n_starts)
	fits[[1]] <- bcd(NULL)                         # deterministic cold start

	if (n_starts > 1L) {
		d  <- dim(Z); nr <- d[1]; nc <- d[2]
		sdZ <- stats::sd(Z[is.finite(Z)])
		if (!is.finite(sdZ) || sdZ == 0) sdZ <- 1
		# size the random factors so the multiplicative term u_i' v_j has sd
		# about 0.5 sd(z): with iid n(0, sc^2) entries, sd(u' v) = sqrt(r) sc^2,
		# so sc = sqrt(0.5 sd(z) / sqrt(r))
		sc <- sqrt(0.5 * sdZ / sqrt(R))

		# save / restore the caller's rng stream
		had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
		old_seed <- if (had_seed) {
			get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
		} else {
			NULL
		}
		on.exit({
			if (had_seed) {
				assign(".Random.seed", old_seed, envir = .GlobalEnv)
			} else if (exists(".Random.seed", envir = .GlobalEnv,
			                  inherits = FALSE)) {
				rm(".Random.seed", envir = .GlobalEnv)
			}
		}, add = TRUE)
		set.seed(seed)

		for (s in 2:n_starts) {
			Us <- matrix(stats::rnorm(nr * R, 0, sc), nr, R)
			if (symmetric) {
				fits[[s]] <- bcd(list(U = Us, V = Us, L = rep(sc, R)))
			} else {
				Vs <- matrix(stats::rnorm(nc * R, 0, sc), nc, R)
				fits[[s]] <- bcd(list(U = Us, V = Vs))
			}
		}
	}

	sse  <- vapply(fits, function(f) f$sse, numeric(1))
	best <- which.min(sse)
	out  <- fits[[best]]
	out$multistart_sse <- sse
	out$multistart_selected <- best
	out
}
