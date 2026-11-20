# reported positions of a static fit. the latent product is identified and its
# factors only up to a gauge, so the positions factor the posterior mean product
# at the fitted rank k, with columns in decreasing order of its singular values:
#   directed:  U V' with the scale split evenly between U and V;
#   bipartite: U G V' with the scale split evenly and G the k x k identity
#              inside the R_row x R_col shape (columns past k are zero);
#   symmetric: U L U' with U orthonormal and L the k eigenvalues largest in size.
# a column's sign is fixed so that its positive entries (over the unit row and
# column vectors) carry more of its squared length than its negative ones.
# posterior sds come from the saved draws, each factored the same way and matched
# onto the reported positions by a move that leaves its product unchanged: a
# rotation of the balanced factors, or for symmetric fits a relabelling and sign
# change of the columns (L is fixed). the matching works on the factors only,
# never on an n x m product

# reported factors of a posterior mean product P at ranks R_u (rows) and R_v
# (columns); returns U, V, G (bipartite), the vector L (symmetric) and the
# strengths d of the k reported dimensions
.uv_static_factors <- function(P, R_u, R_v = R_u, kind = "directed") {
	if (!is.matrix(P)) P <- as.matrix(P)
	if (!all(is.finite(P))) P[!is.finite(P)] <- 0
	k <- min(R_u, R_v, dim(P))
	sign_of <- function(x) if (sum(x * abs(x)) < 0) -1 else 1
	if (kind == "symmetric") {
		e <- eigen((P + t(P)) / 2, symmetric = TRUE)
		pick <- order(-abs(e$values))[seq_len(k)]
		vec <- e$vectors[, pick, drop = FALSE]
		for (j in seq_len(k)) vec[, j] <- vec[, j] * sign_of(vec[, j])
		rownames(vec) <- rownames(P)
		return(list(U = vec, L = e$values[pick], d = abs(e$values[pick])))
	}
	s <- svd(P, nu = k, nv = k)
	for (j in seq_len(k)) {
		f <- sign_of(c(s$u[, j], s$v[, j]))
		s$u[, j] <- s$u[, j] * f
		s$v[, j] <- s$v[, j] * f
	}
	w <- diag(sqrt(s$d[seq_len(k)]), nrow = k)
	U <- matrix(0, nrow(P), R_u, dimnames = list(rownames(P), NULL))
	V <- matrix(0, ncol(P), R_v, dimnames = list(colnames(P), NULL))
	U[, seq_len(k)] <- s$u %*% w
	V[, seq_len(k)] <- s$v %*% w
	out <- list(U = U, V = V, d = s$d[seq_len(k)])
	if (kind == "bipartite") out[["G"]] <- .lda_rect_identity(R_u, R_v)
	out
}

# posterior sds of the reported positions from the saved draws, and each draw's own
# strengths (d_draws, one row per complete draw, in the order of the reported
# dimensions); NULL when fewer than two complete draws are stored
.uv_static_sd <- function(fit, U, V, L, kind) {
	Ud <- fit$U_samples
	Vd <- fit$V_samples
	Gd <- fit$G_samples
	if (length(dim(Ud)) != 3L || length(dim(Vd)) != 3L || dim(Ud)[3L] != dim(Vd)[3L]) return(NULL)
	S_all <- dim(Ud)[3L]
	if (kind == "bipartite" && (length(dim(Gd)) != 3L || dim(Gd)[3L] < S_all)) return(NULL)
	n_u <- nrow(U)
	n_v <- if (kind == "symmetric") n_u else nrow(V)
	if (dim(Ud)[1L] != n_u || dim(Vd)[1L] != n_v) return(NULL)
	R_u <- dim(Ud)[2L]
	R_v <- dim(Vd)[2L]
	k <- if (kind == "symmetric") length(L) else min(R_u, R_v, n_u, n_v)
	if (k == 0L) return(NULL)
	U_sum <- U_sq <- matrix(0, n_u, k)
	V_sum <- V_sq <- matrix(0, n_v, k)
	Ur <- U[, seq_len(k), drop = FALSE]
	Vr <- if (kind == "symmetric") NULL else V[, seq_len(k), drop = FALSE]
	n_pos <- sum(L > 0)
	col_sign <- ifelse(L < 0, -1, 1)
	rows <- seq_len(n_u)
	d_draws <- matrix(NA_real_, S_all, k)
	S <- 0L
	for (s in seq_len(S_all)) {
		Us <- matrix(Ud[, , s], n_u, R_u)
		Vs <- matrix(Vd[, , s], n_v, R_v)
		G_s <- if (kind == "bipartite") matrix(Gd[, , s], R_u, R_v) else diag(1, R_u)
		# only complete draws enter, so unfilled slots of a short run are left out
		if (!all(is.finite(Us)) || !all(is.finite(Vs)) || !all(is.finite(G_s))) next
		S <- S + 1L
		if (kind == "symmetric") {
			pr <- .dyn_draw_pairs(Us, Vs, NULL, k, "symmetric", n_pos)
			m <- .dyn_match(pr, Ur, Ur, rows, rows, "symmetric", col_sign)
			Uc <- .dyn_place(pr, m, L, n_u, n_u)$U
			# a column takes the draw's eigenvalue of its sign and rank within that sign
			pos <- sort(pr$d[pr$d > 0], decreasing = TRUE)
			neg <- sort(-pr$d[pr$d < 0], decreasing = TRUE)
			d_draws[S, L > 0] <- pos[seq_len(n_pos)]
			d_draws[S, L <= 0] <- neg[seq_len(k - n_pos)]
		} else {
			# the draw's product is (qu A)(qv B)' with qu, qv orthonormal bases of U and V
			qu <- qr.Q(qr(Us))
			qv <- qr.Q(qr(Vs))
			core <- crossprod(qu, Us) %*% G_s %*% crossprod(Vs, qv)
			sv <- svd(core, nu = k, nv = k)
			d_draws[S, ] <- sv$d[seq_len(k)]
			w <- sqrt(pmax(sv$d[seq_len(k)], 0))
			A <- sweep(sv$u, 2L, w, "*")
			B <- sweep(sv$v, 2L, w, "*")
			# rotation of the balanced factors onto the reported ones (orthogonal Procrustes)
			M <- crossprod(A, crossprod(qu, Ur)) + crossprod(B, crossprod(qv, Vr))
			Q <- if (all(is.finite(M))) {
				ms <- svd(M)
				ms$u %*% t(ms$v)
			} else diag(k)
			Uc <- qu %*% (A %*% Q)
			Vc <- qv %*% (B %*% Q)
			V_sum <- V_sum + Vc
			V_sq <- V_sq + Vc^2
		}
		U_sum <- U_sum + Uc
		U_sq <- U_sq + Uc^2
	}
	if (S < 2L) return(NULL)
	sd_of <- function(sq, sm, R) {
		out <- matrix(0, nrow(sq), R)
		out[, seq_len(k)] <- sqrt(pmax((sq - sm^2 / S) / (S - 1L), 0))
		out
	}
	U_sd <- sd_of(U_sq, U_sum, ncol(U))
	dimnames(U_sd) <- dimnames(U)
	d_draws <- d_draws[seq_len(S), , drop = FALSE]
	if (kind == "symmetric") return(list(U_sd = U_sd, d_draws = d_draws))
	V_sd <- sd_of(V_sq, V_sum, ncol(V))
	dimnames(V_sd) <- dimnames(V)
	list(U_sd = U_sd, V_sd = V_sd, d_draws = d_draws)
}

# strength and share of each reported dimension; dimensions j and j + 1 are tied
# when the reported gap d_j - d_(j+1) is below tie_z draw-to-draw sds of that gap
# (the data do not order the two strengths, so the axes within their plane are not
# identified). symmetric dimensions of opposite sign are never tied: their
# eigenvectors are identified. tied_with_next is NA without draws
.uv_strength_table <- function(d, d_draws = NULL, L = NULL, tie_z = 2) {
	k <- length(d)
	if (k == 0L) return(NULL)
	tot <- sum(d^2)
	tied <- rep(NA, k)
	if (!is.null(d_draws) && nrow(d_draws) >= 2L) {
		tied <- rep(FALSE, k)
		for (j in seq_len(k - 1L)) {
			if (!is.null(L) && sign(L[j]) != sign(L[j + 1L])) next
			g_s <- d_draws[, j] - d_draws[, j + 1L]
			g_s <- g_s[is.finite(g_s)]
			if (length(g_s) < 2L) {
				tied[j] <- NA
				next
			}
			tied[j] <- d[j] - d[j + 1L] <= tie_z * stats::sd(g_s)
		}
	}
	data.frame(
		dimension = seq_len(k),
		strength = d,
		share = if (is.finite(tot) && tot > 0) d^2 / tot else rep(NA_real_, k),
		tied_with_next = tied
	)
}

# reported positions and posterior sds of a static ame / lame fit (see the top of
# this file); the fit builders and the accessors share it, so a fit saved before
# the positions were summarised this way reads the same as a new one
.lame_static_positions <- function(fit) {
	U0 <- fit[["U"]]
	if (isTRUE(fit$dynamic_uv) || is.null(U0) || length(dim(U0)) < 2L || dim(U0)[2L] == 0L) {
		return(fit)
	}
	bip <- identical(fit$mode, "bipartite")
	sym <- !bip && isTRUE(fit$symmetric)
	kind <- if (bip) "bipartite" else if (sym) "symmetric" else "directed"
	P <- if (sym) fit[["ULUPM"]] %||% fit[["UVPM"]] else fit[["UVPM"]]
	if (is.null(P) || length(dim(P)) != 2L) return(fit)
	V0 <- fit[["V"]]
	R_u <- dim(U0)[2L]
	R_v <- if (bip) {
		if (is.null(V0) || length(dim(V0)) < 2L) return(fit)
		dim(V0)[2L]
	} else R_u
	# a stored product without actor names takes them from the positions
	if (is.null(rownames(P)) && length(dimnames(U0)[[1L]]) == nrow(P)) rownames(P) <- dimnames(U0)[[1L]]
	V_names <- if (bip) dimnames(V0)[[1L]] else dimnames(U0)[[1L]]
	if (is.null(colnames(P)) && length(V_names) == ncol(P)) colnames(P) <- V_names
	f <- .uv_static_factors(P, R_u, R_v, kind)
	fit[["U"]] <- f$U
	if (sym) {
		# L is the diagonal matrix of the k eigenvalues
		fit[["L"]] <- diag(f$L, nrow = length(f$L))
		# a symmetric fit that carries V keeps it as U L
		if (!is.null(fit[["V"]])) fit[["V"]] <- sweep(f$U, 2L, f$L, "*")
	} else {
		fit[["V"]] <- f[["V"]]
	}
	if (bip) fit[["G"]] <- f[["G"]]
	sds <- .uv_static_sd(fit, f$U, f[["V"]], f$L, kind)
	fit$U_sd <- sds$U_sd
	fit$V_sd <- sds$V_sd
	fit$latent_strength <- .uv_strength_table(f$d, sds$d_draws, if (sym) f$L)
	fit$uv_positions <- "posterior_mean_product"
	fit
}

# reported positions of any ame / lame fit, static or dynamic; a fit already
# summarised, or an ALS fit (which reports its own estimates), is returned unchanged
.lame_positions <- function(fit) {
	if (!inherits(fit, c("ame", "lame")) || inherits(fit, "ame_als") ||
			identical(fit$uv_positions, "posterior_mean_product")) {
		return(fit)
	}
	if (isTRUE(fit$dynamic_uv)) .lame_dynamic_positions(fit) else .lame_static_positions(fit)
}

# strength of each reported latent dimension of a fit with static positions: a data
# frame with columns dimension, strength (singular value of the reported product,
# |eigenvalue| for symmetric fits), share (strength^2 / sum(strength^2)) and
# tied_with_next (see .uv_strength_table); NULL for a fit without static positions
# (R = 0, dynamic positions, a withheld ALS product)
.lame_dim_strength <- function(fit) {
	if (inherits(fit, "ame_als")) return(.als_dim_strength(fit))
	if (!inherits(fit, c("ame", "lame")) || isTRUE(fit$dynamic_uv)) return(NULL)
	.lame_positions(fit)$latent_strength
}

# strengths of the reported factors of a static ALS fit, one row per column of U;
# each column takes the singular (eigen) pair of the product it agrees with most
.als_dim_strength <- function(fit) {
	U <- fit[["U"]]
	if (is.null(U) || length(dim(U)) != 2L || ncol(U) == 0L || !all(is.finite(U))) return(NULL)
	k <- ncol(U)
	rows_u <- seq_len(nrow(U))
	if (isTRUE(fit$symmetric)) {
		L <- fit[["L"]]
		if (is.matrix(L)) L <- diag(L)
		if (length(L) != k || !all(is.finite(L))) return(NULL)
		pr <- .dyn_draw_pairs(U, sweep(U, 2L, L, "*"), NULL, k, "symmetric", sum(L > 0))
		m <- .dyn_match(pr, U, U, rows_u, rows_u, "symmetric", ifelse(L < 0, -1, 1))
		d <- abs(pr$d[m$pair_of])
		# a column held at zero matches no pair: its strength is zero
		d[!is.finite(d)] <- 0
		return(.uv_strength_table(d, L = L))
	}
	V <- fit[["V"]]
	if (is.null(V) || length(dim(V)) != 2L || ncol(V) != k || !all(is.finite(V))) return(NULL)
	G <- fit[["G"]]
	if (is.null(G) || !identical(dim(G), c(k, k))) G <- diag(1, k)
	pr <- .dyn_draw_pairs(U, V, G, k, "directed")
	m <- .dyn_match(pr, U, V, rows_u, seq_len(nrow(V)), "directed")
	d <- pr$d[m$pair_of]
	d[!is.finite(d)] <- 0
	.uv_strength_table(d)
}

# runs of tied dimensions of a strength table: each element holds the dimensions
# j..m that share one space
.lame_tie_runs <- function(st) {
	if (is.null(st) || !any(st$tied_with_next %in% TRUE)) return(list())
	tied <- st$tied_with_next %in% TRUE
	start <- which(tied & !c(FALSE, utils::head(tied, -1L)))
	lapply(start, function(j) {
		m <- j
		while (m <= length(tied) && tied[m]) m <- m + 1L
		seq.int(j, m)
	})
}

# the sentence for tied dimensions of a strength table, NULL when none is tied;
# `runs` limits it to some of the runs of .lame_tie_runs
.lame_tie_sentence <- function(st, runs = .lame_tie_runs(st)) {
	if (length(runs) == 0L) return(NULL)
	parts <- vapply(runs, function(dims) {
		if (length(dims) == 2L) {
			sprintf("Dimensions %d and %d have nearly equal strength, so their two axes are not well separated; read them together as a plane.", dims[1L], dims[2L])
		} else {
			sprintf("Dimensions %s and %d have nearly equal strength, so their axes are not well separated; read them together as one space.",
				paste(utils::head(dims, -1L), collapse = ", "), utils::tail(dims, 1L))
		}
	}, character(1))
	paste(parts, collapse = " ")
}

# the summary line of the dimension shares
.lame_share_line <- function(st, rank_txt = NULL) {
	if (is.null(st) || !any(is.finite(st$share))) return(NULL)
	pct <- paste0(formatC(100 * st$share, format = "f", digits = 0), "%", collapse = ", ")
	head_txt <- if (is.null(rank_txt)) "Latent dimensions" else paste0("Latent dimensions, ", rank_txt)
	paste0(head_txt, ": ", pct, " of the latent product")
}

# fitted rank of the latent product: k = min(R_row, R_col) for a bipartite fit with
# unequal ranks, as "rank k (R_row = a, R_col = b)"; NULL otherwise
.lame_rank_text <- function(fit) {
	if (!identical(fit$mode, "bipartite")) return(NULL)
	U <- fit[["U"]]
	V <- fit[["V"]]
	if (is.null(U) || is.null(V) || length(dim(U)) < 2L || length(dim(V)) < 2L) return(NULL)
	R_u <- dim(U)[2L]
	R_v <- dim(V)[2L]
	if (R_u == R_v) return(NULL)
	sprintf("rank %d (R_row = %d, R_col = %d)", min(R_u, R_v), R_u, R_v)
}

# summary lines of the latent dimensions: their shares (with the fitted rank when
# the bipartite ranks differ) and the statement on tied dimensions
.lame_print_latent_lines <- function(x) {
	line <- .lame_share_line(x$latent_strength, x$rank_text)
	if (is.null(line) && !is.null(x$rank_text)) line <- paste0("Latent dimensions: ", x$rank_text)
	if (!is.null(line)) cli::cli_bullets(c("*" = "{line}"))
	tie <- .lame_tie_sentence(x$latent_strength)
	if (!is.null(tie)) cli::cli_alert_info("{tie}", wrap = TRUE)
}
