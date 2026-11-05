# joint move of the coordinate paths of cbin and frn fits under dynamic_uv and the latent
# rows and columns they bound
#
# a sparse nominator's latent row is bounded on one side only, so the augmented sweep moves
# (position, latent row) by about 1 / sqrt(cells) per sweep and the AR(1) prior slows the
# path further. shifting actor i's path by d and the latent values of its cells by W_t d_t
# (W_t = V_t G' when bipartite, V_t when directed) leaves every residual unchanged, so the
# path's conditional is its AR(1) prior times the order conditions of
# nomination_actor_move.R, drawn by an elliptical slice step. senders, then receivers, in
# random order; the compiled visit (src/nomination_position_move.cpp) and the R visit
# (options(lame.npm_engine = "r")) give the same draws from the same seed

# prior roots of the paths: stationary start N(0, sigma^2 / (1 - rho^2)), innovation sd sigma
.lame_npm_root <- function(N, rho, sigma) {
	ch <- tryCatch(chol(.lame_dpm_cov(N, rho, sigma)), error = function(e) NULL)
	if (is.null(ch)) NULL else t(ch)
}

# the slices of the panel (see .lame_spm_nom_slices()); NULL when a cell is hidden
.lame_npm_slices <- function(Z, Y, family, YL, cap, bip) {
	N <- dim(Z)[3L]
	sts <- vector("list", N)
	for (t in seq_len(N)) {
		if (any(.lame_nom_hidden_cells(Y[, , t], bip))) return(NULL)
		st <- .lame_nom_col_slice(Z[, , t], Y[, , t], family, if (identical(family, "frn")) YL[[t]] else NULL, cap, bip)
		st$room_row <- if (identical(st$kind, "frn_h")) !st$full else st$room
		sts[[t]] <- st
	}
	sts
}

# context of the compiled visit: the cell classes, functions of the data alone (NULL when a cell is hidden)
.lame_npm_context <- function(Y, family, YL, cap, bip) {
	if (any(
		vapply(seq_len(dim(Y)[3L]), function(tt) any(.lame_nom_hidden_cells(Y[, , tt], bip)), logical(1))
	)) return(NULL)
	.lame_nom_cells_prep(Y, family, YL, cap, bip)
}

# the compiled visit: the same steps and random numbers as the R visit below
.lame_npm_visit_cpp <- function(ctx, Z, U_cube, V_cube, G, rho_uv, sigma_uv, max_abs = Inf) {
	nr <- dim(Z)[1L]; nc <- dim(Z)[2L]; N <- dim(Z)[3L]
	Z0 <- Z
	stats <- c(npm_calls = 0, npm_evals = 0, npm_stuck = 0)
	state <- function() list(Z = Z, U = U_cube, V = V_cube, shift = Z - Z0, stats = stats)
	if (!(is.finite(rho_uv) && abs(rho_uv) < 1 && is.finite(sigma_uv) && sigma_uv > 0)) return(state())
	root <- .lame_npm_root(N, rho_uv, sigma_uv)
	if (is.null(root)) return(state())
	RA <- dim(U_cube)[2L]; RB <- dim(V_cube)[2L]
	Gm <- if (is.null(G)) matrix(0, 0L, 0L) else G
	o <- npm_nom_phase_cpp(
		as.numeric(Z),
		ctx$cls,
		ctx$kind,
		ctx$room,
		ctx$YL,
		ctx$maxrank,
		nr,
		nc,
		as.numeric(U_cube),
		as.numeric(V_cube),
		RA,
		RB,
		Gm,
		root,
		sample.int(nr),
		TRUE,
		max_abs
	)
	stats <- stats + c(o$calls, o$evals, o$stuck)
	o2 <- npm_nom_phase_cpp(
		o$Z,
		ctx$cls,
		ctx$kind,
		ctx$room,
		ctx$YL,
		ctx$maxrank,
		nr,
		nc,
		o$U,
		o$V,
		RA,
		RB,
		Gm,
		root,
		sample.int(nc),
		FALSE,
		max_abs
	)
	stats <- stats + c(o2$calls, o2$evals, o2$stuck)
	Z[] <- o2$Z
	U_cube[] <- o2$U
	V_cube[] <- o2$V
	state()
}

# one visit. `Z`, `Y` are the latent and observed arrays (nr x nc x N), `YL` the frn
# rank-holder matrices, `U_cube`, `V_cube` the paths, `G` the bipartite interaction
# matrix (NULL when directed), `max_abs` the coordinate bound, and `ctx` the compiled
# context (NULL runs the R visit)
.lame_npm_visit <- function(Z, Y, family, YL, cap, U_cube, V_cube, G, rho_uv, sigma_uv, max_abs = Inf, ctx = NULL) {
	if (!is.null(ctx) && identical(.lame_move_engine("lame.npm_engine"), "cpp")) return(
		.lame_npm_visit_cpp(ctx, Z, U_cube, V_cube, G, rho_uv, sigma_uv, max_abs)
	)
	nr <- dim(Z)[1L]; nc <- dim(Z)[2L]; N <- dim(Z)[3L]
	bip <- !is.null(G)
	Z0 <- Z
	stats <- c(npm_calls = 0, npm_evals = 0, npm_stuck = 0)
	state <- function() list(Z = Z, U = U_cube, V = V_cube, shift = Z - Z0, stats = stats)
	if (!(is.finite(rho_uv) && abs(rho_uv) < 1 && is.finite(sigma_uv) && sigma_uv > 0)) return(state())
	root <- .lame_npm_root(N, rho_uv, sigma_uv)
	if (is.null(root)) return(state())
	sts <- .lame_npm_slices(Z, Y, family, YL, cap, bip)
	if (is.null(sts)) return(state())
	rail <- is.finite(max_abs) && max_abs > 0
	RA <- dim(U_cube)[2L]; RB <- dim(V_cube)[2L]
	draw_path <- function(cur) {
		# prior draw of a path (R x N): the roots act along the periods, coordinate by coordinate
		nu <- matrix(stats::rnorm(length(cur)), nrow(cur))
		t(root %*% t(nu))
	}
	ess <- function(cur, ll) {
		nu <- draw_path(cur)
		thr <- log(stats::runif(1L))
		th <- stats::runif(1L, 0, 2 * pi)
		lo <- th - 2 * pi
		hi <- th
		ev <- 0L
		repeat {
			xp <- cur * cos(th) + nu * sin(th)
			ev <- ev + 1L
			if (!(rail && max(abs(xp)) > max_abs) && ll(xp)) return(list(x = xp, evals = ev, stuck = FALSE))
			if (th < 0) lo <- th else hi <- th
			if (ev >= 60L) return(list(x = cur, evals = ev, stuck = TRUE))
			th <- stats::runif(1L, lo, hi)
		}
	}
	# senders: row i of every period shifts by W_t d_t, W_t = V_t G' (bipartite) or V_t
	for (i in sample.int(nr)) {
		cur <- matrix(U_cube[i, , ], RA, N)
		W <- lapply(seq_len(N), function(t) { Vt <- matrix(V_cube[, , t], nc, RB); if (bip) Vt %*% t(G) else Vt })
		ll <- function(x) {
			d <- x - cur
			for (t in seq_len(N)) {
				st <- sts[[t]]
				if (!.lame_spm_row_ok(
					st$Z[i, ] + drop(W[[t]] %*% d[, t]),
					st$cls[i, ],
					st$kind,
					st$room_row[i]
				)) return(FALSE)
			}
			TRUE
		}
		a <- ess(cur, ll)
		stats["npm_calls"] <- stats["npm_calls"] + 1; stats["npm_evals"] <- stats["npm_evals"] + a$evals; stats["npm_stuck"] <- stats["npm_stuck"] + a$stuck
		d <- a$x - cur
		if (!any(d != 0)) next
		for (t in seq_len(N)) {
			sh <- drop(W[[t]] %*% d[, t])
			.lame_spm_nom_set_row(sts[[t]], i, sts[[t]]$Z[i, ] + sh)
		}
		U_cube[i, , ] <- a$x
	}
	# receivers: column j of every period shifts by (U_t G) d_t or U_t d_t, inside its row limits
	for (j in sample.int(nc)) {
		cur <- matrix(V_cube[j, , ], RB, N)
		W <- lapply(seq_len(N), function(t) { Ut <- matrix(U_cube[, , t], nr, RA); if (bip) Ut %*% G else Ut })
		lims <- lapply(sts, .lame_nom_col_limits, j = j)
		zj <- lapply(sts, function(st) st$Z[, j])
		ll <- function(x) {
			d <- x - cur
			for (t in seq_len(N)) {
				zn <- zj[[t]] + drop(W[[t]] %*% d[, t])
				if (any(zn <= lims[[t]]$lo | zn >= lims[[t]]$hi)) return(FALSE)
			}
			TRUE
		}
		a <- ess(cur, ll)
		stats["npm_calls"] <- stats["npm_calls"] + 1; stats["npm_evals"] <- stats["npm_evals"] + a$evals; stats["npm_stuck"] <- stats["npm_stuck"] + a$stuck
		d <- a$x - cur
		if (!any(d != 0)) next
		for (t in seq_len(N)) .lame_nom_col_apply(sts[[t]], j, drop(W[[t]] %*% d[, t]))
		V_cube[j, , ] <- a$x
	}
	for (t in seq_len(N)) Z[, , t] <- sts[[t]]$Z
	state()
}
