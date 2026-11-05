# joint moves of the static positions and their latent row (or column) in directed cbin, frn and ordinal fits:
# an actor with few events has a wide position that alternating Gibbs draws cross slowly. moving u_i to x shifts
# row i of Z by D_j = (x - u_i)' v_j and keeps every residual, so x is drawn by an elliptical slice step from its
# prior (Suv drawn per visit, R/binary_actor_collapsed_move.R) restricted to rows the data allow; receivers mirror
# this with D_i = u_i' (y - v_j). cbin / frn visits are compiled (src/nomination_position_move.cpp), and
# options(lame.spm_engine = "r") runs the R visit below with the same draws

# whether one cbin / frn row `z` is allowed by the data: nominees above zero and every non-nominee, frn
# nominees increasing in rank order, and non-nominees below zero while the row has room (classes of .lame_nom_col_prep)
.lame_spm_row_ok <- function(z, cls, kind, room) {
	non <- cls == 0L
	nom <- if (identical(kind, "cbin")) cls == 1L else (cls >= 1L | cls == -2L)
	m0 <- if (any(non)) max(z[non]) else -Inf
	if (any(nom)) {
		if (!(min(z[nom]) > max(0, m0))) return(FALSE)
		if (!identical(kind, "cbin")) {
			ch <- cls >= 1L
			if (sum(ch) > 1L) {
				zz <- z[ch][order(cls[ch])]
				if (any(diff(zz) <= 0)) return(FALSE)
			}
		}
	}
	if (isTRUE(room) && m0 >= 0) return(FALSE)
	TRUE
}

# per-slice cell classes, row flags and aggregates for the column limits; NULL when a cell is hidden
# (a missing cell that could hold a nomination), so the host's updates run alone
.lame_spm_nom_slices <- function(Z, Y, family, YL, cap, ref = NULL) {
	N <- dim(Z)[3L]
	sts <- vector("list", N)
	for (t in seq_len(N)) {
		if (any(.lame_nom_hidden_cells(Y[, , t], FALSE))) return(NULL)
		st <- .lame_nom_col_slice(
			Z[, , t],
			Y[, , t],
			family,
			if (identical(family, "frn")) YL[[t]] else NULL,
			cap,
			FALSE
		)
		st$room_row <- if (identical(st$kind, "frn_h")) !st$full else st$room
		sts[[t]] <- st
	}
	sts
}

# the row aggregates of a slice after row i is set to `znew`
.lame_spm_nom_set_row <- function(st, i, znew) {
	st$Z[i, ] <- znew
	cl <- st$cls[i, ]
	st$Zm0[i, ] <- ifelse(cl == 0L, znew, -Inf)
	st$Zmn[i, ] <- ifelse(st$isnom[i, ], -znew, -Inf)
	st$M0[i] <- max(st$Zm0[i, ])
	st$Nq[i] <- max(st$Zmn[i, ])
	invisible(NULL)
}

# prior the positions are drawn under (as in .lame_ltm_uv_visit); NULL when unavailable
.lame_spm_prior <- function(U, V, shrink, uv_prior, Suv) {
	n <- nrow(U)
	R <- ncol(U)
	if (!isTRUE(shrink) && is.null(uv_prior) && is.null(Suv)) return(NULL)
	pri <- tryCatch({
		S <- if (is.null(Suv)) .lame_actor_move_suv(U, V, n, R, shrink, uv_prior) else Suv
		P <- solve(S)
		iu <- seq_len(R)
		iv <- R + seq_len(R)
		cov_u <- solve(P[iu, iu, drop = FALSE])
		cov_v <- solve(P[iv, iv, drop = FALSE])
		cov_u <- (cov_u + t(cov_u)) / 2
		cov_v <- (cov_v + t(cov_v)) / 2
		list(
			root_u = t(chol(cov_u)),
			root_v = t(chol(cov_v)),
			reg_u = -cov_u %*% P[iu, iv, drop = FALSE],
			reg_v = -cov_v %*% P[iv, iu, drop = FALSE]
		)
	}, error = function(e) NULL)
	if (is.null(pri) || any(!is.finite(unlist(pri)))) return(NULL)
	pri
}

# one visit of every sender and receiver position of a directed cbin / frn panel (static positions);
# `uv_prior` list(df, scale) of the factor covariance prior (NULL: unit prior of R/uv_prior.R), `Suv` a fixed covariance
.lame_spm_nom_visit <- function(Z, Y, family, YL, cap, U, V, shrink = TRUE, uv_prior = NULL, Suv = NULL) {
	n <- nrow(U)
	R <- ncol(U)
	N <- dim(Z)[3L]
	stats <- c(spm_calls = 0, spm_evals = 0, spm_stuck = 0)
	Z0 <- Z
	state <- function() list(Z = Z, U = U, V = V, shift = Z - Z0, stats = stats)
	if (R == 0L || n < 3L) return(state())
	pri <- .lame_spm_prior(U, V, shrink, uv_prior, Suv)
	if (is.null(pri)) return(state())
	sts <- .lame_spm_nom_slices(Z, Y, family, YL, cap)
	if (is.null(sts)) return(state())
	# senders: the whole row i shifts by D_j = (x - u_i)' v_j
	for (i in sample.int(n)) {
		cur <- U[i, ]
		ll_row <- function(x) {
			d <- drop(V %*% (x - cur))
			for (st in sts) if (!.lame_spm_row_ok(st$Z[i, ] + d, st$cls[i, ], st$kind, st$room_row[i])) return(-Inf)
			0
		}
		ans <- .lame_ess_step(cur, drop(pri$reg_u %*% V[i, ]), pri$root_u, ll_row, 0)
		stats["spm_calls"] <- stats["spm_calls"] + 1
		stats["spm_evals"] <- stats["spm_evals"] + ans$evals
		stats["spm_stuck"] <- stats["spm_stuck"] + ans$stuck
		dx <- ans$x - cur
		if (!any(dx != 0)) next
		d <- drop(V %*% dx)
		for (t in seq_len(N)) .lame_spm_nom_set_row(sts[[t]], i, sts[[t]]$Z[i, ] + d)
		U[i, ] <- ans$x
	}
	# receivers: column j shifts by D_i = u_i' (y - v_j) within the limits each row sets
	for (j in sample.int(n)) {
		cur <- V[j, ]
		lims <- lapply(sts, .lame_nom_col_limits, j = j)
		zj <- lapply(sts, function(st) st$Z[, j])
		ll_col <- function(y) {
			d <- drop(U %*% (y - cur))
			for (t in seq_len(N)) if (any(zj[[t]] + d <= lims[[t]]$lo | zj[[t]] + d >= lims[[t]]$hi)) return(-Inf)
			0
		}
		ans <- .lame_ess_step(cur, drop(pri$reg_v %*% U[j, ]), pri$root_v, ll_col, 0)
		stats["spm_calls"] <- stats["spm_calls"] + 1
		stats["spm_evals"] <- stats["spm_evals"] + ans$evals
		stats["spm_stuck"] <- stats["spm_stuck"] + ans$stuck
		dx <- ans$x - cur
		if (!any(dx != 0)) next
		d <- drop(U %*% dx)
		for (t in seq_len(N)) .lame_nom_col_apply(sts[[t]], j, d)
		V[j, ] <- ans$x
	}
	for (t in seq_len(N)) Z[, , t] <- sts[[t]]$Z
	state()
}

# the extremes of the moving cells of an ordinal line, per category and period (K x N matrices of the maximum and
# minimum and of their positions within the period's matrix); `ZM`, `WM` the values and codes of the line (cells x N)
.lame_spm_ord_line <- function(cx, ZM, WM, i, row) {
	storage.mode(WM) <- "integer"
	spm_ord_line_cpp(ZM, WM, cx$K, as.integer(i), row, cx$nr)
}

# whether line `ZM` with codes `WM` fits the data given the other cells' extremes (`omax`, `omin`): per period
# the rank-likelihood ordering, or every cell inside its cutpoint interval (`af` = c(-Inf, alpha, Inf))
.lame_spm_ord_ok <- function(cx, ZM, WM, omax, omin, af = NULL, i = 0L, row = TRUE) {
	spm_ord_ok_cpp(
		ZM,
		WM,
		if (is.null(omax)) matrix(0, 1L, 1L) else omax,
		if (is.null(omin)) matrix(0, 1L, 1L) else omin,
		if (is.null(af)) numeric(0) else af
	)
}

# one visit of every sender and receiver position of a directed ordinal panel (static positions);
# `cx` from .lame_ltm_context, `alpha` explicit cutpoints (NULL: rank likelihood), other arguments as .lame_spm_nom_visit
.lame_spm_ord_visit <- function(cx, Z, U, V, shrink = TRUE, uv_prior = NULL, Suv = NULL, alpha = NULL) {
	n <- nrow(U)
	R <- ncol(U)
	N <- cx$N
	stats <- c(spm_calls = 0, spm_evals = 0, spm_stuck = 0)
	Z0 <- Z
	state <- function() list(Z = Z, U = U, V = V, shift = Z - Z0, stats = stats)
	if (R == 0L || n < 3L) return(state())
	pri <- .lame_spm_prior(U, V, shrink, uv_prior, Suv)
	if (is.null(pri)) return(state())
	af <- if (is.null(alpha)) NULL else c(-Inf, alpha, Inf)
	for (row in c(TRUE, FALSE)) {
		cache <- if (is.null(af)) .lame_ltm_ord_cache(cx, Z) else NULL
		for (i in sample.int(n)) {
			if (row) {
				ZM <- matrix(Z[i, , ], n, N)
				WM <- matrix(cx$W[i, , ], n, N)
				cur <- U[i, ]
				opp <- V
				mu <- drop(pri$reg_u %*% V[i, ])
				root <- pri$root_u
			} else {
				ZM <- matrix(Z[, i, ], n, N)
				WM <- matrix(cx$W[, i, ], n, N)
				cur <- V[i, ]
				opp <- U
				mu <- drop(pri$reg_v %*% U[i, ])
				root <- pri$root_v
			}
			if (is.null(af)) {
				bd <- .lame_ltm_ord_bounds_cached(cx, Z, cache, i, row)
				omax <- bd$upd$omax
				omin <- bd$upd$omin
				omax[is.na(omax)] <- -Inf
				omin[is.na(omin)] <- Inf
			} else omax <- omin <- NULL
			ll <- function(x) {
				d <- drop(opp %*% (x - cur))
				if (.lame_spm_ord_ok(cx, ZM + d, WM, omax, omin, af, i, row)) 0 else -Inf
			}
			ans <- .lame_ess_step(cur, mu, root, ll, 0)
			stats["spm_calls"] <- stats["spm_calls"] + 1
			stats["spm_evals"] <- stats["spm_evals"] + ans$evals
			stats["spm_stuck"] <- stats["spm_stuck"] + ans$stuck
			dx <- ans$x - cur
			if (!any(dx != 0)) next
			d <- drop(opp %*% dx)
			ZM1 <- ZM + d
			if (row) Z[i, , ] <- ZM1 else Z[, i, ] <- ZM1
			if (is.null(af)) {
				ln <- .lame_spm_ord_line(cx, ZM1, WM, i, row)
				upd <- bd$upd
				upd[c("mmax", "mmin", "lm", "lmn")] <- ln[c("mmax", "mmin", "lm", "lmn")]
				cache <- .lame_ltm_ord_update(cache, upd, 0)
			}
			if (row) U[i, ] <- ans$x else V[i, ] <- ans$x
		}
	}
	state()
}

# the panel a fit hands to the visit; `ok` is FALSE when the move does not serve it (hidden cbin / frn cells,
# an ordinal panel the translation context refuses), and the host's updates then run alone
.lame_spm_context <- function(Y, family, YL = NULL, cap = NULL, Y_int = NULL) {
	d <- dim(Y)
	Y <- array(Y, c(d[1L], d[2L], if (length(d) >= 3L) d[3L] else 1L))
	ctx <- list(family = family, Y = Y, YL = YL, cap = cap, ok = d[1L] == d[2L] && d[1L] >= 3L)
	if (!ctx$ok) return(ctx)
	if (identical(family, "ordinal")) {
		ctx$cx <- .lame_ltm_context(Y, "ordinal", "dir", Y_int = Y_int)
		ctx$ok <- isTRUE(ctx$cx$ok)
	} else {
		ctx$ok <- !any(
			vapply(seq_len(dim(Y)[3L]), function(t) any(.lame_nom_hidden_cells(Y[, , t], FALSE)), logical(1))
		)
		if (ctx$ok) ctx$prep <- .lame_nom_cells_prep(Y, family, YL, cap, FALSE)
	}
	ctx
}

# the engine of a compiled move: options(<name> = "r") runs the R reference, which gives the same draws from the same seed
.lame_move_engine <- function(name) {
	opt <- getOption(name, "cpp")
	if (!is.character(opt) || length(opt) != 1L || !opt %in% c("cpp", "r")) cli::cli_abort("{.code options({name})} must be {.val cpp} or {.val r}.")
	opt
}

# cell classes of a cbin / frn panel without hidden cells for the compiled visits, formed once per fit
.lame_nom_cells_prep <- function(Y, family, YL, cap, bip) {
	nr <- dim(Y)[1L]; nc <- dim(Y)[2L]; N <- dim(Y)[3L]
	cls <- array(0L, dim(Y))
	room <- matrix(0L, nr, N)
	yl <- NULL
	mr <- 0L
	for (t in seq_len(N)) {
		pr <- .lame_nom_col_prep(
			matrix(0, nr, nc),
			Y[, , t],
			family,
			if (identical(family, "frn")) YL[[t]] else NULL,
			cap,
			bip
		)
		cls[, , t] <- pr$cls
		room[, t] <- as.integer(pr$room)
		if (identical(family, "frn")) {
			if (is.null(yl)) { mr <- pr$max_rank; yl <- array(0L, c(nr, mr, N)) }
			h <- pr$YL
			h[is.na(h)] <- 0L
			yl[, , t] <- h
		}
	}
	list(
		cls = as.integer(cls),
		kind = if (identical(family, "cbin")) 1L else 2L,
		room = as.integer(room),
		YL = if (is.null(yl)) integer(0) else as.integer(yl),
		maxrank = mr
	)
}

# compiled visit of cbin / frn static positions (src/nomination_position_move.cpp), with the same steps and
# random numbers as .lame_spm_nom_visit()
.lame_spm_nom_visit_cpp <- function(ctx, Z, U, V, shrink = TRUE, uv_prior = NULL, Suv = NULL) {
	n <- nrow(U)
	R <- ncol(U)
	stats <- c(spm_calls = 0, spm_evals = 0, spm_stuck = 0)
	Z0 <- Z
	state <- function() list(Z = Z, U = U, V = V, shift = Z - Z0, stats = stats)
	if (R == 0L || n < 3L) return(state())
	pri <- .lame_spm_prior(U, V, shrink, uv_prior, Suv)
	if (is.null(pri)) return(state())
	pp <- ctx$prep
	zz <- as.numeric(Z)
	o <- spm_nom_phase_cpp(
		zz,
		pp$cls,
		pp$kind,
		pp$room,
		pp$YL,
		pp$maxrank,
		n,
		n,
		U,
		V,
		pri$reg_u,
		pri$root_u,
		sample.int(n),
		TRUE
	)
	U[] <- o$X
	stats <- stats + c(o$calls, o$evals, o$stuck)
	o2 <- spm_nom_phase_cpp(
		o$Z,
		pp$cls,
		pp$kind,
		pp$room,
		pp$YL,
		pp$maxrank,
		n,
		n,
		V,
		U,
		pri$reg_v,
		pri$root_v,
		sample.int(n),
		FALSE
	)
	V[] <- o2$X
	stats <- stats + c(o2$calls, o2$evals, o2$stuck)
	Z[] <- o2$Z
	state()
}

# one visit of the panel `ctx`; `Z` the latent array (n x n x N)
.lame_spm_visit <- function(ctx, Z, U, V, shrink = TRUE, uv_prior = NULL, Suv = NULL, alpha = NULL) {
	Z <- array(Z, dim(ctx$Y))
	if (identical(ctx$family, "ordinal")) return(.lame_spm_ord_visit(ctx$cx, Z, U, V, shrink, uv_prior, Suv, alpha))
	if (!is.null(ctx$prep) && identical(.lame_move_engine("lame.spm_engine"), "cpp")) return(
		.lame_spm_nom_visit_cpp(ctx, Z, U, V, shrink, uv_prior, Suv)
	)
	.lame_spm_nom_visit(Z, ctx$Y, ctx$family, ctx$YL, ctx$cap, U, V, shrink, uv_prior, Suv)
}
