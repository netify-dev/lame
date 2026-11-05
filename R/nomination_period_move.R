# period-level move of the intercept path of a cbin / frn dynamic_beta fit (lame()). a period in which nearly
# every row is at its cap leaves its intercept weakly informed, so the Gibbs pair (intercept_t, Z_t) moves slowly,
# and .lame_nom_level_shift cannot help because rows with room in other periods hold it.
# a row below its cap in period t is a pinned cell (i, t) whose beta_t + a_i must not change; periods linked by a
# common pinned row form a component K. with r = |R_K| / n_rows the move shifts beta_t by (1 - r) delta on P_K
# and -r delta elsewhere, a_i by -(1 - r) delta on R_K and r delta elsewhere, and the free latent cells so every
# residual is kept; delta is then a truncated normal, drawn exactly. pinned / free labels depend only on each
# row slice, which the translation keeps, so the step is an exact Gibbs draw along a translation group

# per-period bounds of the row translations: the shifts (lo, hi) of each row of
# each period (matrices rows x periods) that keep the row consistent with its
# observed cells (see .lame_nom_shift_bounds for the reduction over periods)
.lame_nom_shift_bounds_t <- function(Z, odmax, masks) {
	n_r <- dim(Z)[1L]
	N <- length(masks)
	cap <- rep_len(as.numeric(odmax), n_r)
	rr <- seq_len(n_r)
	lo <- matrix(-Inf, n_r, N)
	hi <- matrix(Inf, n_r, N)
	for (t in seq_len(N)) {
		mk <- masks[[t]]
		z <- Z[, , t]
		z1 <- -z; z1[!mk$o1] <- -Inf
		m1 <- -z1[cbind(rr, max.col(z1, ties.method = "first"))]
		lo[, t] <- ifelse(mk$has1, -m1, -Inf)
		z0 <- z; z0[!mk$o0] <- -Inf
		M0 <- z0[cbind(rr, max.col(z0, ties.method = "first"))]
		q0 <- rowSums(mk$elig & (z > M0))
		hi[, t] <- ifelse(mk$has0 & q0 < cap, -M0, Inf)
	}
	list(lo = lo, hi = hi)
}

# the components of the graph of periods and pinned rows: a list of
# list(periods, rows), `free` the rows x periods matrix of cells a translation
# may carry (TRUE) against the pinned ones (FALSE)
.lame_nom_period_components <- function(free) {
	N <- ncol(free)
	lab <- seq_len(N)
	pin <- !free
	multi <- which(rowSums(pin) > 1L)
	repeat {
		changed <- FALSE
		for (i in multi) {
			ps <- which(pin[i, ])
			m <- min(lab[ps])
			if (any(lab[ps] != m)) {
				lab[lab %in% lab[ps]] <- m
				changed <- TRUE
			}
		}
		if (!changed) break
	}
	lapply(split(seq_len(N), lab), function(P) {
		list(periods = P, rows = which(rowSums(pin[, P, drop = FALSE]) > 0L))
	})
}

# one sweep over the components in random order. `ic` the intercept path from .lame_nom_icpt_state (kind "path");
# the other arguments as .lame_nom_level_shift. returns Z, a, ic and the counts
.lame_nom_period_shift <- function(Z, a, b, Sab, odmax, masks, bip, rows_ok, ic) {
	n_r <- dim(Z)[1L]
	N <- dim(Z)[3L]
	none <- list(Z = Z, a = a, ic = ic, calls = 0, comps = 0L)
	if (is.null(ic) || !identical(ic$kind, "path") || N < 2L || length(ic$path) != N) return(none)
	bd <- .lame_nom_shift_bounds_t(Z, odmax, masks)
	# a cell is free when its row is not shrunk, is bounded above by nothing (the row is at
	# its cap or has no observed non-nomination) and satisfies its own bounds; the rest are
	# pinned and stay where they are
	free <- rows_ok & bd$lo < 0 & is.infinite(bd$hi)
	if (!any(free)) return(none)
	comps <- .lame_nom_period_components(free)
	lo <- bd$lo
	if (bip) {
		m <- rep(0, n_r)
		s2 <- max(Sab[1L, 1L], 1e-12)
	} else {
		r_ab <- if (is.finite(Sab[2L, 2L]) && Sab[2L, 2L] > 0) Sab[1L, 2L] / Sab[2L, 2L] else 0
		m <- b * r_ab
		s2 <- max(Sab[1L, 1L] - r_ab * Sab[1L, 2L], 1e-12)
	}
	ord <- if (length(comps) > 1L) sample.int(length(comps)) else 1L
	calls <- 0
	for (k in ord) {
		cm <- comps[[k]]
		P <- cm$periods
		R <- cm$rows
		nR <- length(R)
		if (nR >= n_r) next
		r <- nR / n_r
		inP <- seq_len(N) %in% P
		inR <- seq_len(n_r) %in% R
		ua <- ifelse(inR, -(1 - r), r)
		ut <- ifelse(inP, 1 - r, -r)
		path <- ic$path
		du <- ut[-1L] - ic$rho * ut[-N]
		Q <- sum(ua^2) / s2 + ut[1L]^2 / ic$v0 + sum(du^2) / ic$sd2
		lin <- sum(ua * (a - m)) / s2 + ut[1L] * (path[1L] - ic$m0) / ic$v0 +
			sum(du * (path[-1L] - ic$rho * path[-N])) / ic$sd2
		if (!is.finite(Q) || Q <= 0 || !is.finite(lin)) next
		# the cells translated up (rows outside R_K in the periods of K) bound
		# delta from below by -(weakest nomination); the cells translated down (rows
		# of R_K in the other periods) bound it from above by the weakest nomination
		up <- matrix(FALSE, n_r, N)
		up[!inR, P] <- TRUE
		dn <- matrix(FALSE, n_r, N)
		if (length(P) < N) dn[R, -P] <- TRUE
		lower <- if (any(up)) max(lo[up]) else -Inf
		upper <- if (any(dn)) min(-lo[dn]) else Inf
		sd_d <- 1 / sqrt(Q)
		if (!(lower < upper)) next
		delta <- sd_d * rtnorm_interval_logp((-lin / Q) / sd_d, lower / sd_d, upper / sd_d)
		if (!is.finite(delta)) next
		for (t in which(colSums(up) > 0L)) {
			rows <- which(up[, t])
			Z[rows, , t] <- Z[rows, , t] + delta
			lo[rows, t] <- lo[rows, t] - delta
		}
		for (t in which(colSums(dn) > 0L)) {
			rows <- which(dn[, t])
			Z[rows, , t] <- Z[rows, , t] - delta
			lo[rows, t] <- lo[rows, t] + delta
		}
		a <- a + ua * delta
		ic$path <- path + ut * delta
		calls <- calls + 1
	}
	list(Z = Z, a = a, ic = ic, calls = calls, comps = length(comps))
}
