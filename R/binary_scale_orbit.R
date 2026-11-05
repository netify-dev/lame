# exact draw of the scale of a binary probit predictor with the latent array integrated out
#
# when a covariate (or an actor, or a factor block) separates ties from non-ties the
# posterior is wide along the ray that multiplies the coefficients and the additive
# effects by c > 0, and the augmented sweep crawls along it in steps of 1 / sqrt(cells).
# this move scales the static coefficients and additive effects theta_s (d of them) by c,
# a group move with invariant measure dc / c (Liu and Wu 1999) and density
#     c^(d - 1) exp(-(A c^2 - 2 B c) / 2) prod_cells Phi(s (c e_s + m + off) / sd),
# drawn by a slice step on log c (Neal 2003). a directed network keeps the shared dyad
# factor of binary_actor_collapsed_move.R, and Z is redrawn given it and the scaled state

# the index sets of the observed cells of a panel, for the likelihood and for the refresh of the latent array
# `Y` the outcomes (nr x nc x T), mode "dir", "sym" or "bip"
.lame_bso_context <- function(Y, mode) {
	d <- dim(Y)
	nr <- d[1L]
	nc <- d[2L]
	nt <- if (length(d) >= 3L) d[3L] else 1L
	Y <- array(Y, c(nr, nc, nt))
	cx <- list(mode = mode, nr = nr, nc = nc, nt = nt, Y = Y)
	if (identical(mode, "sym")) {
		# the observed value of an unordered pair is its upper cell's, or its lower cell's when that is missing
		# (the kernel of the latent array reads the upper triangle first)
		iu <- which(upper.tri(matrix(0, nr, nc)))
		ij <- arrayInd(iu, c(nr, nc))
		lo <- ij[, 2L] + (ij[, 1L] - 1L) * nr
		up_idx <- lo_idx <- integer(0)
		sg <- numeric(0)
		for (t in seq_len(nt)) {
			y <- Y[, , t]
			pair <- y[iu]
			alt <- y[lo]
			pair[is.na(pair)] <- alt[is.na(pair)]
			ok <- is.finite(pair)
			off <- (t - 1L) * nr * nc
			up_idx <- c(up_idx, iu[ok] + off)
			lo_idx <- c(lo_idx, lo[ok] + off)
			sg <- c(sg, 2 * pair[ok] - 1)
		}
		cx$up <- up_idx
		cx$lo <- lo_idx
		cx$sg <- sg
		return(cx)
	}
	dg <- array(FALSE, dim(Y))
	if (identical(mode, "dir")) dg[cbind(
		rep(seq_len(nr), nt),
		rep(seq_len(nr), nt),
		rep(seq_len(nt), each = nr)
	)] <- TRUE
	cx$idx <- which(is.finite(Y) & !dg)
	cx$sg <- 2 * Y[cx$idx] - 1
	cx$i1 <- which(Y == 1 & !dg)
	cx$i0 <- which(Y == 0 & !dg)
	cx$ina <- which(is.na(Y) & !dg)
	cx$idg <- which(dg)
	cx$upper <- upper.tri(matrix(0, nr, nc))
	cx
}

# one slice step on u = log c for the log density `lg`, started at u = 0 (Neal 2003:
# stepping out with width `w` and at most `m` steps, then shrinkage)
.lame_bso_slice <- function(lg, w = 0.25, m = 40L) {
	lg0 <- lg(0)
	evals <- 1L
	if (!is.finite(lg0)) return(list(u = 0, evals = evals))
	thr <- lg0 + log(stats::runif(1L))
	left <- -w * stats::runif(1L)
	right <- left + w
	j <- floor(m * stats::runif(1L))
	k <- (m - 1L) - j
	while (j > 0L) {
		evals <- evals + 1L
		if (!(lg(left) > thr)) break
		left <- left - w
		j <- j - 1L
	}
	while (k > 0L) {
		evals <- evals + 1L
		if (!(lg(right) > thr)) break
		right <- right + w
		k <- k - 1L
	}
	repeat {
		u <- stats::runif(1L, left, right)
		val <- lg(u)
		evals <- evals + 1L
		if (is.finite(val) && val > thr) return(list(u = u, evals = evals))
		if (u < 0) left <- u else right <- u
		# the bracket collapses only if the density fails numerically: the state stays
		if (right - left < 1e-12) return(list(u = 0, evals = evals))
	}
}

# one visit. `Z` and `E` are the latent array and the predictor (nr x nc x T), `m` the part
# of the predictor the scale does not carry, `rho` the dyadic correlation (zero unless
# directed), and `A_pri`, `B_pri`, `d_th` the prior terms and number of scaled parameters.
# returns the latent array and the factor c the caller applies to the scaled parameters
.lame_bso_visit <- function(Z, E, m, cx, rho, A_pri, B_pri, d_th, w = 0.25) {
	mode <- cx$mode
	dim3 <- c(cx$nr, cx$nc, cx$nt)
	Z <- array(Z, dim3)
	E <- array(E, dim3)
	still <- list(Z = Z, c = 1, evals = 0L)
	if (!all(is.finite(E)) || !is.finite(A_pri) || !is.finite(B_pri) || !is.finite(d_th) || A_pri <= 0) return(still)
	dirn <- identical(mode, "dir")
	if (dirn && (!is.finite(rho) || abs(rho) >= 1)) return(still)
	off <- if (dirn && rho != 0) .lame_dyad_factor_draw(Z, E, rho, cx$upper) else NULL
	sdc <- if (dirn) sqrt(1 - abs(rho)) else 1
	cells <- if (identical(mode, "sym")) function(a) (a[cx$up] + a[cx$lo]) / 2 else function(a) a[cx$idx]
	e_all <- cells(E)
	if (length(e_all) == 0L) return(still)
	m_all <- if (is.null(m)) numeric(length(e_all)) else cells(array(m, dim3))
	es <- e_all - m_all
	mm <- m_all + if (is.null(off)) 0 else off[cx$idx]
	if (length(mm) == 1L) mm <- rep(mm, length(es))
	if (!all(is.finite(es)) || !all(is.finite(mm))) return(still)
	inv_sd <- 1 / sdc
	lg <- function(u) {
		k <- exp(u)
		d_th * u - A_pri * k^2 / 2 + B_pri * k + bso_loglik_cpp(es, mm, cx$sg, k, inv_sd)
	}
	sl <- .lame_bso_slice(lg, w)
	k <- exp(sl$u)
	if (k == 1) {
		still$evals <- sl$evals
		return(still)
	}
	E1 <- if (is.null(m)) k * E else k * (E - array(m, dim3)) + array(m, dim3)
	Z1 <- if (identical(mode, "sym")) {
		for (t in seq_len(cx$nt)) Z[, , t] <- rZ_bin_sym_fc(Z[, , t], E1[, , t], cx$Y[, , t], rho)
		Z
	} else {
		.lame_binary_z_given_factor(
			Z,
			E1,
			if (is.null(off)) array(0, dim3) else off,
			if (dirn) rho else 0,
			cx$i1,
			cx$i0,
			cx$ina,
			cx$idg
		)
	}
	list(Z = Z1, c = k, evals = sl$evals)
}
