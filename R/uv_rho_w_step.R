# block draw of the persistence rho and the stationary variance w = sigma^2 / (1 - rho^2) of AR(1)
# coordinate paths. in (rho, w) the cap w <= prior$uv_var_cap bounds a rectangle, so the draws do not
# stall against the curved boundary the cap forms in (rho, sigma). with S1, n1 the sum of squares and
# count of the first period, Stt, Sxy, Syy the sums of u_t^2, u_t u_{t-1}, u_{t-1}^2 over the ntr
# transitions, and (a, b) the inverse-gamma prior of sigma^2:
#   log f(rho, w) = -((n1 + ntr) / 2 + a + 1) log w - (ntr / 2 + a) log(1 - rho^2) - Q(rho) / (2 w) + log pi(rho)
#   Q(rho) = S1 + (Stt - 2 rho Sxy + rho^2 Syy + 2 b) / (1 - rho^2)
# w given rho is an inverse gamma truncated at the cap; rho given w is a slice step (Neal 2003)

# `U_cube`, `V_cube` the coordinate paths (n x r x N; V_cube is ignored when symmetric); `rho`, `sigma`
# the current values; `prior_mean`, `prior_sd` the normal prior of rho; `shape`, `scale` the inverse-gamma
# prior of sigma^2; `sigma_max`, `stat_var_cap` the bounds on sigma and w; `lambda_u`, `lambda_v` (n x N)
# the Student-t precision weights, where a transition without a positive finite weight is left out.
# returns the new rho and sigma
.lame_uv_rho_w_step <- function(
	U_cube,
	V_cube,
	rho,
	sigma,
	prior_mean,
	prior_sd,
	shape,
	scale,
	sigma_max = Inf,
	stat_var_cap = Inf,
	symmetric = FALSE,
	bound = 0.99,
	lambda_u = NULL,
	lambda_v = NULL
) {
	N <- dim(U_cube)[3L]
	if (is.null(N) || N < 2L || !is.finite(rho) || !is.finite(sigma) || sigma <= 0) return(
		list(rho = rho, sigma = sigma)
	)
	S1 <- n1 <- Stt <- Sxy <- Syy <- ntr <- 0
	use_v <- !(symmetric || is.null(V_cube))
	cubes <- if (use_v) list(U_cube, V_cube) else list(U_cube)
	lams <- if (use_v) list(lambda_u, lambda_v) else list(lambda_u)
	for (k in seq_along(cubes)) {
		cu <- cubes[[k]]
		lam <- lams[[k]]
		S1 <- S1 + sum(cu[, , 1L]^2)
		n1 <- n1 + prod(dim(cu)[1:2])
		if (is.null(lam)) {
			cur <- cu[, , -1L, drop = FALSE]
			prev <- cu[, , -N, drop = FALSE]
			Stt <- Stt + sum(cur^2)
			Sxy <- Sxy + sum(cur * prev)
			Syy <- Syy + sum(prev^2)
			ntr <- ntr + length(cur)
		} else {
			for (t in 2:N) {
				w <- lam[, t]
				keep <- is.finite(w) & w > 0
				if (!any(keep)) next
				cur <- matrix(cu[keep, , t], nrow = sum(keep))
				prev <- matrix(cu[keep, , t - 1L], nrow = sum(keep))
				wk <- w[keep]
				Stt <- Stt + sum(wk * rowSums(cur^2))
				Sxy <- Sxy + sum(wk * rowSums(cur * prev))
				Syy <- Syy + sum(wk * rowSums(prev^2))
				ntr <- ntr + length(cur)
			}
		}
	}
	if (!all(is.finite(c(S1, Stt, Sxy, Syy)))) return(list(rho = rho, sigma = sigma))
	rho <- max(-bound, min(bound, rho))
	w <- sigma^2 / (1 - rho * rho)
	lp <- function(r) {
		if (!is.finite(r) || abs(r) > bound) return(-Inf)
		om <- 1 - r * r
		if (is.finite(sigma_max) && w > sigma_max^2 / om) return(-Inf)
		-(ntr / 2 + shape) * log(om) - (S1 + (Stt - 2 * r * Sxy + r * r * Syy + 2 * scale) / om) / (2 * w) +
			dnorm(r, prior_mean, prior_sd, log = TRUE)
	}
	lp0 <- lp(rho)
	if (is.finite(lp0)) {
		# slice step for rho given w: step out from a random start until both ends leave the slice or reach
		# the support bound, then shrink
		logy <- lp0 + log(runif(1))
		step <- 0.02
		lo <- max(rho - step * runif(1), -bound)
		hi <- min(lo + step, bound)
		while (lo > -bound && lp(lo) > logy) lo <- max(lo - step, -bound)
		while (hi < bound && lp(hi) > logy) hi <- min(hi + step, bound)
		repeat {
			cand <- runif(1, lo, hi)
			if (lp(cand) > logy) {
				rho <- cand
				break
			}
			if (cand < rho) lo <- cand else hi <- cand
			if (hi - lo < 1e-12) break
		}
	}
	om <- 1 - rho * rho
	rate <- (S1 + (Stt - 2 * rho * Sxy + rho * rho * Syy + 2 * scale) / om) / 2
	shape_w <- (n1 + ntr) / 2 + shape
	wmax <- min(stat_var_cap, if (is.finite(sigma_max)) sigma_max^2 / om else Inf)
	# draw on the sd scale of w: sqrt(1 / g), g ~ Gamma(shape, rate) cut at 1 / bound^2, in log survival space
	sw <- if (is.finite(wmax)) .lame_rinvgamma_sd_trunc(shape_w, rate, sqrt(wmax), sqrt(w)) else sqrt(
		1 / rgamma(1, shape = shape_w, rate = rate)
	)
	list(rho = rho, sigma = sw * sqrt(om))
}
