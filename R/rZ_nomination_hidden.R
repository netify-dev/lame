# hidden cells of cbin nomination rows: a row with hidden cells takes a complete-data sweep given the
# labels its latent values imply (.lame_cbin_implied_labels), then a redraw of each hidden cell given
# the rest of the row and the observed cells (.lame_cbin_hidden_moves)

# hidden cells of each row: NA, off the diagonal, in a column observed in the slice and a row with an
# observed cell; `bip` is the fit's mode, so a square bipartite slice keeps its diagonal cells
.lame_nom_hidden_cells <- function(Y, bip) {
	stopifnot(is.logical(bip), length(bip) == 1L, !is.na(bip))
	obs <- !is.na(Y)
	h <- !obs
	if (!bip) {
		diag(h) <- FALSE
		diag(obs) <- FALSE
	}
	h[, colSums(obs) == 0L] <- FALSE
	h[rowSums(obs) == 0L, ] <- FALSE
	h
}

# nominations each row has made (the odobs of the cap check): observed nominations plus hidden cells
# with positive latent values; cells of an actor absent from the slice never count
.lame_nom_room <- function(odobs, hid, Z) {
	odobs + rowSums(hid & (Z > 0))
}

# refuse a nomination cap below a row's observed out-degree; od is one per row or rows x periods,
# odmax a scalar or one per row
.lame_check_odmax_outdegree <- function(od, odmax, actors = NULL) {
	od <- as.matrix(od)
	od_max <- apply(od, 1L, max)
	cap <- rep_len(as.numeric(odmax), length(od_max))
	bad <- which(od_max > cap)
	if (!length(bad)) return(invisible(NULL))
	nm <- if (length(actors) == length(od_max)) as.character(actors) else as.character(seq_along(od_max))
	show <- utils::head(bad, 5L)
	rows <- paste0(nm[show], " (", od_max[show], " nominations, cap ", cap[show], ")")
	if (length(bad) > length(show)) rows <- c(rows, paste0(length(bad) - length(show), " more"))
	cli::cli_abort(
		c(
			"{.arg odmax} is below the observed out-degree of {length(bad)} row{?s}: {rows}.",
			"i" = "A row cannot make more nominations than its cap, so these values of {.arg odmax} contradict {.arg Y}.",
			"i" = "Set {.arg odmax} to at least each row's observed out-degree (the largest is {max(od_max)}), or pass {.code odmax = NULL} to use the largest observed out-degree."
		),
		class = "lame_odmax_below_outdegree"
	)
}

# YO template for the posterior-predictive simY_frn() draws: Y with hidden cells set to zero, so they
# compete for the row's cap as in the samplers; only the diagonal and absent actors' cells stay out
.lame_nom_predictive_yo <- function(Y, bip) {
	YO <- Y
	YO[.lame_nom_hidden_cells(Y, bip)] <- 0
	YO
}

# nomination labels the current latent values imply for the hidden cells: the highest positive ones,
# as many as the cap leaves after the observed nominations (compiled in cbin_implied_labels_cpp)
.lame_cbin_implied_labels <- function(Z, Y, odmax, hid) {
	cbin_implied_labels_cpp(Z, Y, as.numeric(rep_len(odmax, nrow(Y))), hid)
}
# R reference of cbin_implied_labels_cpp, used in the tests
.lame_cbin_implied_labels_r <- function(Z, Y, odmax, hid) {
	Ya <- Y
	cap <- rep_len(odmax, nrow(Y))
	n1 <- rowSums(!is.na(Y) & Y > 0)
	for (i in which(rowSums(hid) > 0L)) {
		h <- which(hid[i, ])
		zh <- Z[i, h]
		lab <- numeric(length(h))
		r <- max(cap[i] - n1[i], 0)
		pos <- which(zh > 0)
		if (r > 0 && length(pos)) {
			top <- pos[order(zh[pos], decreasing = TRUE)][seq_len(min(r, length(pos)))]
			lab[top] <- 1
		}
		Ya[i, h] <- lab
	}
	Ya
}

# log of the standard normal mass on (lo, hi), stable in both tails
.lame_log_norm_mass <- function(lo, hi) {
	out <- rep(-Inf, length(lo))
	ok <- hi > lo
	up <- ok & lo >= 0
	dn <- ok & !up
	if (any(up)) {
		a <- stats::pnorm(lo[up], lower.tail = FALSE, log.p = TRUE)
		b <- stats::pnorm(hi[up], lower.tail = FALSE, log.p = TRUE)
		out[up] <- a + log1p(-exp(b - a))
	}
	if (any(dn)) {
		a <- stats::pnorm(hi[dn], log.p = TRUE)
		b <- stats::pnorm(lo[dn], log.p = TRUE)
		out[dn] <- a + log1p(-exp(b - a))
	}
	out
}

# redraw each hidden cell from its conditional given the rest of its row and the observed cells; the
# row's consistency changes only at 0, M0 and m1, so the support is a union of intervals. one hidden
# cell per row per round, the upper triangle before the lower when rho couples the rows
.lame_cbin_hidden_moves <- function(Z, EZ, rho, Y, odmax, hid, bip) {
	n_r <- nrow(Z)
	# hidden cells of each row in a random order (the same draws as which() and sample.int() per row)
	hl <- split(col(hid)[hid], factor(row(hid)[hid], levels = seq_len(n_r)))
	names(hl) <- NULL
	for (i in which(lengths(hl) > 1L)) hl[[i]] <- hl[[i]][sample.int(length(hl[[i]]))]
	if (!max(lengths(hl))) return(Z)
	# compiled rounds; .lame_cbin_hidden_moves_r is the R reference with the same draws
	cbin_hidden_moves_cpp(Z, EZ, rho, !bip, as.numeric(rep_len(odmax, n_r)), Y, hid, hl)
}
.lame_cbin_hidden_moves_r <- function(Z, EZ, rho, Y, odmax, hid, bip) {
	n_r <- nrow(Z); n_c <- ncol(Z)
	sq <- !bip
	cond_sd <- if (sq) sqrt(1 - rho^2) else 1
	cap <- rep_len(odmax, n_r)
	o1 <- !is.na(Y) & Y > 0
	o0 <- !is.na(Y) & Y == 0
	if (sq) { diag(o1) <- FALSE; diag(o0) <- FALSE }
	n1 <- rowSums(o1)
	r <- cap - n1
	z1 <- Z; z1[!o1] <- Inf
	m1 <- apply(z1, 1L, min)
	z0 <- Z; z0[!o0] <- -Inf
	M0 <- apply(z0, 1L, max)
	P0 <- rowSums(o0 & Z > 0)
	# the hidden cells of each row, in a random order
	hl <- lapply(seq_len(n_r), function(i) {
		h <- which(hid[i, ])
		if (length(h) > 1L) h[sample.int(length(h))] else h
	})
	n_round <- max(lengths(hl))
	if (!n_round) return(Z)
	for (k in seq_len(n_round)) {
		rows <- which(lengths(hl) >= k)
		cols <- vapply(rows, function(i) hl[[i]][k], integer(1))
		blocks <- if (sq) list(rows < cols, rows > cols) else list(rep(TRUE, length(rows)))
		for (bl in blocks) {
			if (!any(bl)) next
			ii <- rows[bl]; jj <- cols[bl]
			pos <- ii + (jj - 1L) * n_r
			x_old <- Z[pos]
			# the other hidden cells of each row above 0, M0 and m1
			cnt <- function(t) rowSums(hid[ii, , drop = FALSE] & (Z[ii, , drop = FALSE] > t)) - (x_old > t)
			a0 <- cnt(0); aM <- cnt(M0[ii]); am <- cnt(m1[ii])
			M <- M0[ii]; m <- m1[ii]; rr <- r[ii]; cp <- cap[ii]; nn <- n1[ii]; p0 <- P0[ii]
			Mp <- pmax(M, 0)
			lo <- cbind(-Inf, 0, Mp, m)
			hi <- cbind(0, Mp, m, Inf)
			feas <- matrix(FALSE, length(ii), 4L)
			for (g in 1:4) {
				gt0 <- g >= 2L
				gtM <- g >= 3L & M > 0
				gtm <- g == 4L
				K <- nn + a0 + p0 + gt0
				ok_le <- K <= cp & p0 == 0
				ok_gt <- K > cp & (am + gtm) <= rr & (M <= 0 | (aM + gtM) >= rr)
				feas[, g] <- (ok_le | ok_gt) & hi[, g] > lo[, g]
			}
			mu <- EZ[pos]
			if (sq && rho != 0) {
				mp <- jj + (ii - 1L) * n_r
				mu <- mu + rho * (Z[mp] - EZ[mp])
			}
			lm <- .lame_log_norm_mass((lo - mu) / cond_sd, (hi - mu) / cond_sd)
			lm <- matrix(lm, length(ii), 4L)
			lm[!feas] <- -Inf
			mx <- apply(lm, 1L, max)
			good <- is.finite(mx)
			if (!any(good)) next
			w <- exp(lm[good, , drop = FALSE] - mx[good])
			w <- w / rowSums(w)
			u <- stats::runif(sum(good))
			g_pick <- 1L + rowSums(u > t(apply(w, 1L, cumsum)))
			g_pick <- pmin(g_pick, 4L)
			sel <- cbind(seq_len(sum(good)), g_pick)
			lb <- lo[good, , drop = FALSE][sel]; ub <- hi[good, , drop = FALSE][sel]
			val <- mu[good] + cond_sd * rtnorm_interval_logp(
				rep(0, sum(good)),
				(lb - mu[good]) / cond_sd,
				(ub - mu[good]) / cond_sd
			)
			bad <- !is.finite(val)
			if (any(bad)) val[bad] <- x_old[good][bad]
			Z[pos[good]] <- val
		}
	}
	Z
}

# frn rows with hidden cells: observed ranks are read for their order only, since a hidden
# nomination could hold any rank; every cell is drawn given the implied labels, then the hidden
# cells are redrawn as for cbin (.lame_cbin_hidden_moves)
.rZ_frn_hidden_sweep <- function(Z, EZ, rho, Y, Ya, odmax, hid, bip) {
	n_r <- nrow(Z); n_c <- ncol(Z)
	sq <- !bip
	cond_sd <- if (sq) sqrt(1 - rho^2) else 1
	cap <- rep_len(odmax, n_r)
	obs <- !is.na(Y)
	if (sq) diag(obs) <- FALSE
	o1 <- obs & Y > 0
	h1 <- hid & !is.na(Ya) & Ya > 0
	non <- (obs & Y == 0) | (hid & !is.na(Ya) & Ya == 0)
	# the order of the observed nominations within each row, 1 the weakest,
	# and the column holding each order position
	k_row <- rowSums(o1)
	k_max <- max(c(0L, k_row))
	P <- matrix(NA_integer_, n_r, max(k_max, 1L))
	cls <- matrix(-1L, n_r, n_c)
	cls[non] <- 0L
	cls[h1] <- -2L
	if (k_max > 0L) {
		idx <- which(o1)
		rr <- row(o1)[idx]; cc <- col(o1)[idx]
		o <- order(rr, Y[idx])
		rr <- rr[o]; cc <- cc[o]
		qpos <- sequence(tabulate(rr, n_r))
		P[cbind(rr, qpos)] <- cc
		cls[cbind(rr, cc)] <- qpos
	}
	full <- k_row + rowSums(h1) >= cap
	if (sq) diag(cls) <- -3L
	rows <- seq_len(n_r)
	row_of <- row(Z)
	held <- function(q) {
		out <- rep(NA_real_, n_r)
		has <- !is.na(P[, q])
		out[has] <- Z[cbind(rows[has], P[has, q])]
		out
	}
	limits <- function(cl) {
		lo <- rep(-Inf, n_r); hi <- rep(Inf, n_r)
		if (cl == 0L || cl == -2L || cl == 1L) {
			zt <- Z; zt[!non] <- -Inf
			maxnon <- apply(zt, 1L, max)
		}
		if (cl == 0L) {
			zt <- Z; zt[!(o1 | h1)] <- Inf
			hi <- ifelse(full, apply(zt, 1L, min), 0)
		} else if (cl == -2L) {
			lo <- pmax(0, maxnon)
		} else if (cl >= 1L) {
			lo <- if (cl == 1L) pmax(0, maxnon) else held(cl - 1L)
			if (cl < k_max) hi <- held(cl + 1L)
		}
		lo[is.na(lo)] <- -Inf
		hi[is.na(hi)] <- Inf
		list(lo = lo, hi = hi)
	}
	tmask <- if (sq) list(U = upper.tri(Z), L = lower.tri(Z)) else list(A = matrix(TRUE, n_r, n_c))
	for (cl in sample(c(-2L, -1L, 0L, seq_len(k_max)))) {
		lim <- limits(cl)
		for (side in sample(names(tmask))) {
			cells <- tmask[[side]] & cls == cl
			if (!any(cells)) next
			r <- row_of[cells]
			mu <- EZ[cells]
			if (sq && rho != 0) mu <- mu + rho * (t(Z)[cells] - t(EZ)[cells])
			vals <- cond_sd * rtnorm_interval_logp(mu / cond_sd, lim$lo[r] / cond_sd, lim$hi[r] / cond_sd)
			bad <- !is.finite(vals)
			if (any(bad)) vals[bad] <- Z[cells][bad]
			Z[cells] <- vals
		}
	}
	if (sq) diag(Z) <- stats::rnorm(n_r, diag(EZ), sqrt(1 + rho))
	Z
}

# joint move of a sender's level a_i and its latent row: shifting both by delta leaves z - ez
# unchanged, so a_i is drawn from its conditional prior restricted to the shifts that keep the row
# consistent with its observed cells (a row at its cap mixes slowly under the Gibbs sweep alone)
.lame_nom_shift_masks <- function(Y, bip) {
	lapply(seq_len(dim(Y)[3L]), function(t) {
		y <- Y[, , t]
		obs <- !is.na(y)
		if (!bip) diag(obs) <- FALSE
		hid <- .lame_nom_hidden_cells(y, bip)
		o1 <- obs & y > 0
		o0 <- obs & y == 0
		list(
			o1 = o1,
			o0 = o0,
			elig = obs | hid,
			has1 = rowSums(o1) > 0L,
			has0 = rowSums(o0) > 0L
		)
	})
}
# the shifts (lo, hi) of each row of the latent array Z (rows x columns x
# periods) that keep the row consistent with its observed cells
.lame_nom_shift_bounds <- function(Z, odmax, masks) {
	n_r <- dim(Z)[1L]
	cap <- rep_len(as.numeric(odmax), n_r)
	lo <- rep(-Inf, n_r)
	hi <- rep(Inf, n_r)
	rr <- seq_len(n_r)
	for (t in seq_along(masks)) {
		mk <- masks[[t]]
		z <- Z[, , t]
		# row minimum of the observed nominations and row maximum of the
		# observed non-nominations (max.col with ties "first" draws nothing)
		z1 <- -z; z1[!mk$o1] <- -Inf
		m1 <- -z1[cbind(rr, max.col(z1, ties.method = "first"))]
		lo <- pmax(lo, ifelse(mk$has1, -m1, -Inf))
		z0 <- z; z0[!mk$o0] <- -Inf
		M0 <- z0[cbind(rr, max.col(z0, ties.method = "first"))]
		# cells that could hold a nomination above the strongest observed
		# non-nomination, row by row (M0 recycles down the columns)
		q0 <- rowSums(mk$elig & (z > M0))
		hi <- pmin(hi, ifelse(mk$has0 & q0 < cap, -M0, Inf))
	}
	list(lo = lo, hi = hi)
}
.lame_nom_row_shift <- function(Z, a, b, Sab, odmax, masks, bip, rows_ok, center = FALSE) {
	n_r <- dim(Z)[1L]
	bd <- .lame_nom_shift_bounds(Z, odmax, masks)
	lo <- bd$lo
	hi <- bd$hi
	if (bip) {
		m <- rep(0, n_r)
		s <- sqrt(max(Sab[1L, 1L], 1e-12))
	} else {
		r_ab <- if (is.finite(Sab[2L, 2L]) && Sab[2L, 2L] > 0) Sab[1L, 2L] / Sab[2L, 2L] else 0
		m <- b * r_ab
		s <- sqrt(max(Sab[1L, 1L] - r_ab * Sab[1L, 2L], 1e-12))
	}
	ii <- which(rows_ok & lo < 0 & hi > 0)
	if (!length(ii)) return(list(Z = Z, a = a))
	if (isTRUE(center)) {
		# pairwise translations keep sum(a) = 0: given the pair's sum the contrast is Gaussian on the
		# intersection of both rows' constraints, and disjoint random pairs form a valid sweep
		if (length(ii) < 2L) return(list(Z = Z, a = a))
		ii <- ii[sample.int(length(ii))]
		for (k in seq.int(1L, length(ii) - 1L, by = 2L)) {
			i <- ii[k]
			j <- ii[k + 1L]
			lower <- max(lo[i], -hi[j])
			upper <- min(hi[i], -lo[j])
			if (lower >= upper) next
			sd_delta <- s / sqrt(2)
			mean_delta <- ((m[i] - a[i]) - (m[j] - a[j])) / 2
			delta <- sd_delta * rtnorm_interval_logp(
				mean_delta / sd_delta,
				lower / sd_delta,
				upper / sd_delta
			)
			if (!is.finite(delta)) next
			a[i] <- a[i] + delta
			a[j] <- a[j] - delta
			Z[i, , ] <- Z[i, , ] + delta
			Z[j, , ] <- Z[j, , ] - delta
		}
		return(list(Z = Z, a = a))
	}
	val <- s * rtnorm_interval_logp(m[ii] / s, (a[ii] + lo[ii]) / s, (a[ii] + hi[ii]) / s)
	ok <- is.finite(val)
	ii <- ii[ok]
	if (!length(ii)) return(list(Z = Z, a = a))
	d <- val[ok] - a[ii]
	Z[ii, , ] <- Z[ii, , ] + d
	a[ii] <- val[ok]
	list(Z = Z, a = a)
}

# effective sample size of the intercept (and va) of a cbin / frn fit with hidden cells, kept in
# fit$hidden_cell_ess as a record of mixing (NULL otherwise); .lame_cap_level_unreported decides reporting
.lame_nom_hidden_ess <- function(fit, has_hidden) {
	fam <- fit[["family"]] %||% ""
	if (!isTRUE(has_hidden) || length(fam) != 1L || !(fam %in% c("cbin", "frn"))) return(NULL)
	ess1 <- function(x) {
		x <- as.numeric(x)
		x <- x[is.finite(x)]
		if (length(x) < 10L || stats::sd(x) < 1e-12) return(NA_real_)
		as.numeric(coda::effectiveSize(coda::mcmc(x)))
	}
	out <- c(intercept = NA_real_, va = NA_real_, draws = NA_real_)
	B <- fit$BETA
	if (!is.null(B) && length(B)) {
		nm <- if (length(dim(B)) == 3L) dimnames(B)[[2L]] else colnames(B)
		if ("intercept" %in% nm) {
			x <- if (length(dim(B)) == 3L) rowMeans(matrix(B[, "intercept", ], nrow = dim(B)[1L])) else B[, "intercept"]
			out[["intercept"]] <- ess1(x)
			out[["draws"]] <- length(x)
		}
	}
	if (is.matrix(fit$VC) && "va" %in% colnames(fit$VC)) {
		out[["va"]] <- ess1(fit$VC[, "va"])
		if (!is.finite(out[["draws"]])) out[["draws"]] <- nrow(fit$VC)
	}
	out
}
