# dynamic als/map point estimator for longitudinal ame models
#
# covers normal, binary, and poisson panels with optional named
# changing-composition alignment and optional
# dynamic additive effects, selected dynamic regression coefficients, static
# low-rank uv when r > 0, smooth dynamic uv, bipartite dynamic g, and
# node-covariate coefficients identified by the als orthogonality decomposition.
# dynamic node coefficients use period-specific node values; static node
# coefficients use the same actor-level summaries as static als. snap-shift uv
# uses the separate lame_snap_als() route.

.lda_default <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.lda_max_abs <- function(x) {
	if (is.null(x) || length(x) == 0L) return(0)
	out <- max(abs(x), na.rm = TRUE)
	if (is.finite(out)) out else 0
}

.lda_presence_row_means <- function(mat, presence = NULL) {
	if (is.null(presence) || !all(dim(presence) == dim(mat))) {
		return(rowMeans(mat))
	}
	out <- numeric(nrow(mat))
	for (i in seq_len(nrow(mat))) {
		ok <- presence[i, ] & is.finite(mat[i, ])
		out[i] <- if (any(ok)) mean(mat[i, ok]) else 0
	}
	out
}

.lda_presence_col_vars <- function(mat, presence = NULL) {
	if (is.null(presence) || !all(dim(presence) == dim(mat))) {
		return(apply(mat, 2, stats::var))
	}
	vapply(seq_len(ncol(mat)), function(t) {
		ok <- presence[, t] & is.finite(mat[, t])
		if (sum(ok) > 1L) stats::var(mat[ok, t]) else NA_real_
	}, numeric(1))
}

.lda_presence_col_cov <- function(a_mat, b_mat, row_presence = NULL,
                                  col_presence = NULL) {
	if (is.null(row_presence) || is.null(col_presence) ||
	    !all(dim(row_presence) == dim(a_mat)) ||
	    !all(dim(col_presence) == dim(b_mat)) ||
	    nrow(a_mat) != nrow(b_mat)) {
		return(vapply(seq_len(ncol(a_mat)), function(t) {
			stats::cov(a_mat[, t], b_mat[, t])
		}, numeric(1)))
	}
	vapply(seq_len(ncol(a_mat)), function(t) {
		ok <- row_presence[, t] & col_presence[, t] &
			is.finite(a_mat[, t]) & is.finite(b_mat[, t])
		if (sum(ok) > 1L) stats::cov(a_mat[ok, t], b_mat[ok, t]) else NA_real_
	}, numeric(1))
}

.lda_transition_from_presence <- function(presence) {
	if (is.null(presence) || ncol(presence) < 2L) return(NULL)
	out <- matrix(1, nrow(presence), ncol(presence),
	              dimnames = dimnames(presence))
	p <- presence
	p[is.na(p)] <- FALSE
	out[, -1L] <- as.numeric(p[, -1L, drop = FALSE] &
		p[, -ncol(p), drop = FALSE])
	out
}

.lda_combine_transition_weights <- function(weights, presence_weights) {
	if (is.null(presence_weights)) return(weights)
	if (is.null(weights)) return(presence_weights)
	if (!all(dim(weights) == dim(presence_weights))) return(weights)
	weights * presence_weights
}

.lda_clip_rho <- function(x, default = 0.95) {
	if (is.null(x) || length(x) == 0L || !is.finite(x[1L])) return(default)
	pmin(pmax(as.numeric(x[1L]), 0), 0.999)
}

.lda_path_penalty <- function(Tt, lambda, kind = "rw1", rho = 1) {
	K <- matrix(0, Tt, Tt)
	if (lambda <= 0 || Tt < 2L) return(K)
	if (identical(kind, "rw2")) {
		if (Tt < 3L) {
			cli::cli_abort("{.arg dynamic_beta_kind = \"rw2\"} requires at least 3 time periods.")
		}
		D <- matrix(0, Tt - 2L, Tt)
		for (t in seq_len(Tt - 2L)) {
			D[t, t] <- 1
			D[t, t + 1L] <- -2
			D[t, t + 2L] <- 1
		}
		return(lambda * crossprod(D))
	}
	D <- matrix(0, Tt - 1L, Tt)
	for (t in seq_len(Tt - 1L)) {
		D[t, t] <- if (identical(kind, "ar1")) -rho else -1
		D[t, t + 1L] <- 1
	}
	lambda * crossprod(D)
}

.lda_path_penalty_weighted <- function(Tt, lambda, weights = NULL, rho = 1) {
	K <- matrix(0, Tt, Tt)
	if (lambda <= 0 || Tt < 2L) return(K)
	w <- if (is.null(weights)) rep(1, Tt) else as.numeric(weights)
	if (length(w) != Tt) w <- rep(1, Tt)
	w[!is.finite(w) | w < 0] <- 0
	D <- matrix(0, Tt - 1L, Tt)
	for (t in seq_len(Tt - 1L)) {
		D[t, t] <- -rho
		D[t, t + 1L] <- 1
	}
	lambda * crossprod(sqrt(w[-1L]) * D)
}

.lda_solve_path <- function(obs_prec, rhs, lambda, kind = "rw1",
                            rho = 1, ridge = 1e-8,
                            transition_weights = NULL) {
	Tt <- length(rhs)
	P <- if (!is.null(transition_weights) && !identical(kind, "rw2")) {
		.lda_path_penalty_weighted(
			Tt, lambda, transition_weights,
			rho = if (identical(kind, "ar1")) rho else 1)
	} else {
		.lda_path_penalty(Tt, lambda, kind, rho)
	}
	K <- diag(obs_prec, nrow = Tt) + P
	K <- K + diag(ridge, Tt)
	.ae_safe_solve(K, rhs)
}

.lda_unit_weights <- function(Z) {
	W <- array(0, dim = dim(Z))
	W[is.finite(Z)] <- 1
	W
}

.lda_eta_array <- function(E) {
	array(unlist(E, use.names = FALSE), dim = c(nrow(E[[1L]]), ncol(E[[1L]]), length(E)))
}

.lda_irls_work <- function(Yarr, eta_arr, family, link) {
	work_fun <- switch(paste(family, link),
		"poisson log"   = function(eta) .ae_irls_poisson(Yarr, eta),
		"binary probit" = function(eta) .ae_irls_binary_probit(Yarr, eta),
		"binary logit"  = function(eta) .ae_irls_binary_logit(Yarr, eta),
		cli::cli_abort("IRLS is not available for family {.val {family}} with link {.val {link}}."))
	wk <- work_fun(eta_arr)
	Z <- wk$z
	W <- wk$w
	Z[!is.finite(Yarr)] <- NA_real_
	W[!is.finite(Yarr)] <- NA_real_
	W[!is.finite(W) | W < 0] <- 0
	list(Z = Z, W = W)
}

.lda_o_array <- function(O, Tt) {
	if (is.null(O)) return(NULL)
	d <- dim(O)
	if (length(d) == 3L) return(O)
	array(rep(O, Tt), dim = c(nrow(O), ncol(O), Tt))
}

.lda_o_slice <- function(Oarr, t) {
	if (length(dim(Oarr)) == 3L) Oarr[, , t] else Oarr
}

.lda_symmetric_factor_slice <- function(O, R, row_names = NULL,
                                       dim_names = NULL) {
	O <- (O + t(O)) / 2
	n <- nrow(O)
	if (R <= 0L) {
		U <- V <- matrix(0, n, 0L)
		return(list(U = U, V = V, O = matrix(0, n, n)))
	}
	eg <- eigen(O, symmetric = TRUE)
	rr <- min(R, length(eg$values))
	idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(rr)]
	vals <- eg$values[idx]
	scale <- sqrt(pmax(abs(vals), 0))
	U <- sweep(eg$vectors[, idx, drop = FALSE], 2, scale, "*")
	V <- sweep(U, 2, sign(vals), "*")
	if (rr < R) {
		U <- cbind(U, matrix(0, n, R - rr))
		V <- cbind(V, matrix(0, n, R - rr))
	}
	U <- U[, seq_len(R), drop = FALSE]
	V <- V[, seq_len(R), drop = FALSE]
	if (is.null(dim_names)) dim_names <- paste0("dim", seq_len(R))
	dimnames(U) <- dimnames(V) <- list(row_names, dim_names)
	O_rank <- tcrossprod(U, V)
	dimnames(O_rank) <- list(row_names, row_names)
	list(U = U, V = V, O = O_rank)
}

.lda_symmetric_o_to_uv_cube <- function(Oarr, R, row_names = NULL,
                                        dim_names = NULL,
                                        time_names = NULL) {
	if (length(dim(Oarr)) != 3L) Oarr <- .lda_o_array(Oarr, 1L)
	d <- dim(Oarr)
	if (is.null(row_names)) row_names <- dimnames(Oarr)[[1L]]
	if (is.null(time_names)) time_names <- dimnames(Oarr)[[3L]]
	if (is.null(dim_names)) dim_names <- paste0("dim", seq_len(R))
	U <- array(0, dim = c(d[1L], R, d[3L]),
	           dimnames = list(row_names, dim_names, time_names))
	V <- array(0, dim = c(d[1L], R, d[3L]),
	           dimnames = list(row_names, dim_names, time_names))
	O_rank <- array(0, dim = d,
	                dimnames = list(row_names, row_names, time_names))
	for (t in seq_len(d[3L])) {
		f <- .lda_symmetric_factor_slice(Oarr[, , t], R, row_names,
		                                 dim_names)
		U[, , t] <- f$U
		V[, , t] <- f$V
		O_rank[, , t] <- f$O
	}
	list(U = U, V = V, Oarr = O_rank)
}

.lda_rect_identity <- function(nr, nc) {
	out <- matrix(0, nr, nc)
	if (nr > 0L && nc > 0L) {
		for (k in seq_len(min(nr, nc))) out[k, k] <- 1
	}
	out
}

.lda_coerce_factor <- function(M, n, R, scale = 0.1, seed = NULL) {
	if (R <= 0L) return(matrix(0, n, 0L))
	if (!is.null(seed)) set.seed(seed)
	out <- matrix(0, n, R)
	if (!is.null(M) && length(dim(M)) == 2L && nrow(M) == n) {
		k <- min(ncol(M), R)
		if (k > 0L) out[, seq_len(k)] <- M[, seq_len(k), drop = FALSE]
		if (R > k) out[, (k + 1L):R] <- stats::rnorm(n * (R - k), 0, scale)
	} else {
		out[,] <- stats::rnorm(n * R, 0, scale)
	}
	out
}

.lda_uv_to_o_array <- function(U, V) {
	if (is.null(U) || is.null(V)) return(NULL)
	du <- dim(U)
	dv <- dim(V)
	if (length(du) == 2L) {
		return(array(tcrossprod(U, V), dim = c(du[1L], dv[1L], 1L)))
	}
	Tt <- du[3L]
	out <- array(0, dim = c(du[1L], dv[1L], Tt))
	for (t in seq_len(Tt)) {
		Ut <- matrix(U[, , t], nrow = du[1L], ncol = du[2L])
		Vt <- matrix(V[, , t], nrow = dv[1L], ncol = dv[2L])
		out[, , t] <- tcrossprod(Ut, Vt)
	}
	dimnames(out) <- list(dimnames(U)[[1L]], dimnames(V)[[1L]], dimnames(U)[[3L]])
	out
}

.lda_uvg_to_o_array <- function(U, V, G, Tt = NULL, time_names = NULL) {
	if (is.null(U) || is.null(V) || is.null(G)) return(NULL)
	du <- dim(U)
	dv <- dim(V)
	dg <- dim(G)
	Tt <- max(if (length(du) == 3L) du[3L] else 1L,
	          if (length(dv) == 3L) dv[3L] else 1L,
	          if (length(dg) == 3L) dg[3L] else 1L,
	          .lda_default(Tt, 1L))
	out <- array(0, dim = c(du[1L], dv[1L], Tt))
	for (t in seq_len(Tt)) {
		Ut <- if (length(du) == 3L) {
			matrix(U[, , t], nrow = du[1L], ncol = du[2L])
		} else U
		Vt <- if (length(dv) == 3L) {
			matrix(V[, , t], nrow = dv[1L], ncol = dv[2L])
		} else V
		Gt <- if (length(dg) == 3L) {
			matrix(G[, , t], nrow = dg[1L], ncol = dg[2L])
		} else G
			out[, , t] <- Ut %*% Gt %*% t(Vt)
	}
	time_names <- .lda_default(time_names, if (length(du) == 3L) {
		dimnames(U)[[3L]]
	} else if (length(dv) == 3L) {
		dimnames(V)[[3L]]
	} else if (length(dg) == 3L) {
		dimnames(G)[[3L]]
	} else NULL)
	dimnames(out) <- list(dimnames(U)[[1L]], dimnames(V)[[1L]], time_names)
	out
}

.lda_procrustes_q <- function(A, ref) {
	if (ncol(A) == 1L) {
		sgn <- sign(sum(A * ref, na.rm = TRUE))
		if (!is.finite(sgn) || sgn == 0) sgn <- 1
		return(matrix(sgn, 1L, 1L))
	}
	M <- crossprod(A, ref)
	if (!all(is.finite(M))) return(diag(ncol(A)))
	sv <- tryCatch(svd(M), error = function(e) NULL)
	if (is.null(sv)) return(diag(ncol(A)))
	sv$u %*% t(sv$v)
}

.lda_align_uv_cube <- function(U, V) {
	if (is.null(U) || is.null(V) || length(dim(U)) != 3L || dim(U)[3L] < 2L) {
		return(list(U = U, V = V))
	}
	Tt <- dim(U)[3L]
	for (t in 2:Tt) {
		Ut <- matrix(U[, , t], nrow = dim(U)[1L], ncol = dim(U)[2L])
		Vt <- matrix(V[, , t], nrow = dim(V)[1L], ncol = dim(V)[2L])
		Ur <- matrix(U[, , t - 1L], nrow = dim(U)[1L], ncol = dim(U)[2L])
		Vr <- matrix(V[, , t - 1L], nrow = dim(V)[1L], ncol = dim(V)[2L])
		A <- rbind(Ut, Vt)
		ref <- rbind(Ur, Vr)
		Q <- .lda_procrustes_q(A, ref)
		U[, , t] <- Ut %*% Q
		V[, , t] <- Vt %*% Q
	}
	list(U = U, V = V)
}

.lda_align_uvg_cube <- function(U, V, G) {
	if (is.null(U) || is.null(V) || is.null(G) ||
	    length(dim(U)) != 3L || length(dim(V)) != 3L ||
	    length(dim(G)) != 3L || dim(U)[3L] < 2L) {
		return(list(U = U, V = V, G = G))
	}
	Tt <- dim(U)[3L]
	for (t in 2:Tt) {
		Ut <- matrix(U[, , t], nrow = dim(U)[1L], ncol = dim(U)[2L])
		Vt <- matrix(V[, , t], nrow = dim(V)[1L], ncol = dim(V)[2L])
		Gt <- matrix(G[, , t], nrow = dim(G)[1L], ncol = dim(G)[2L])
		Ur <- matrix(U[, , t - 1L], nrow = dim(U)[1L], ncol = dim(U)[2L])
		Vr <- matrix(V[, , t - 1L], nrow = dim(V)[1L], ncol = dim(V)[2L])
		Qu <- .lda_procrustes_q(Ut, Ur)
		Qv <- .lda_procrustes_q(Vt, Vr)
		U[, , t] <- Ut %*% Qu
		V[, , t] <- Vt %*% Qv
		G[, , t] <- t(Qu) %*% Gt %*% Qv
	}
	list(U = U, V = V, G = G)
}

.lda_coef_index <- function(mask, Tt) {
	p <- length(mask)
	static_idx <- which(!mask)
	dynamic_idx <- which(mask)
	n_static <- length(static_idx)
	n_dynamic <- length(dynamic_idx) * Tt
	static_map <- rep(NA_integer_, p)
	if (n_static > 0L) static_map[static_idx] <- seq_len(n_static)
	dynamic_map <- matrix(NA_integer_, p, Tt)
	pos <- n_static
	for (k in dynamic_idx) {
		dynamic_map[k, ] <- pos + seq_len(Tt)
		pos <- pos + Tt
	}
	list(n = n_static + n_dynamic,
	     static_idx = static_idx,
	     dynamic_idx = dynamic_idx,
	     static_map = static_map,
	     dynamic_map = dynamic_map)
}

.lda_coef_layout <- function(prep, include_nodes = TRUE) {
	nms <- c("intercept", prep$x_names_dyad)
	block <- c("intercept", rep("dyad", prep$p_dyad))
	if (isTRUE(include_nodes)) {
		if (!is.null(prep$W_row) && ncol(prep$W_row) > 0L) {
			nms <- c(nms, colnames(prep$W_row))
			block <- c(block, rep("row", ncol(prep$W_row)))
		}
		if (!is.null(prep$W_col) && ncol(prep$W_col) > 0L) {
			nms <- c(nms, colnames(prep$W_col))
			block <- c(block, rep("col", ncol(prep$W_col)))
		}
	}
	list(names = nms, block = block)
}

.lda_beta_dyn_from_mask <- function(mask, coef_names, coef_block) {
	mask <- as.logical(mask)
	groups <- ifelse(mask, coef_block, "")
	names(mask) <- coef_names
	names(groups) <- coef_names
	group_names <- sort(unique(groups[groups != ""]))
	list(any = any(mask),
	     mask = mask,
	     static_idx = which(!mask),
	     dynamic_idx = which(mask),
	     groups = groups,
	     group_names = group_names,
	     n_groups = length(group_names))
}

.lda_dynamic_node_names <- function(beta_dyn) {
	if (is.null(beta_dyn) || !isTRUE(beta_dyn$any)) return(character(0))
	mask <- beta_dyn$mask
	groups <- beta_dyn$groups
	if (is.null(names(mask)) || is.null(names(groups))) return(character(0))
	names(mask)[mask & groups %in% c("row", "col")]
}

.lda_node_matrix <- function(prep, margin, t = NULL, dynamic = FALSE) {
	W <- if (identical(margin, "row")) prep$W_row else prep$W_col
	W_arr <- if (identical(margin, "row")) prep$W_row_arr else prep$W_col_arr
	if (isTRUE(dynamic) && !is.null(W_arr) && !is.null(t)) {
		out <- W_arr[, , t, drop = FALSE]
		out <- matrix(out, nrow = dim(W_arr)[1L], ncol = dim(W_arr)[2L])
		dimnames(out) <- list(NULL, dimnames(W_arr)[[2L]])
		return(out)
	}
	W
}

.lda_bind_node_matrices <- function(W_row, W_col) {
	if (is.null(W_row)) return(W_col)
	if (is.null(W_col)) return(W_row)
	cbind(W_row, W_col)
}

.lda_bind_node_arrays <- function(A, B) {
	if (is.null(A)) return(B)
	if (is.null(B)) return(A)
	n <- dim(A)[1L]
	p <- dim(A)[2L] + dim(B)[2L]
	Tt <- dim(A)[3L]
	out <- array(0, dim = c(n, p, Tt))
	out[, seq_len(dim(A)[2L]), ] <- A
	out[, dim(A)[2L] + seq_len(dim(B)[2L]), ] <- B
	dimnames(out) <- list(NULL,
	                      c(dimnames(A)[[2L]], dimnames(B)[[2L]]),
	                      dimnames(A)[[3L]])
	out
}

.lda_node_varies_names <- function(prep, margin) {
	W_arr <- if (identical(margin, "row")) prep$W_row_arr else prep$W_col_arr
	if (is.null(W_arr) || dim(W_arr)[3L] < 2L) return(character(0))
	nms <- dimnames(W_arr)[[2L]]
	varies <- logical(dim(W_arr)[2L])
	for (k in seq_along(varies)) {
		vals <- W_arr[, k, , drop = FALSE]
		vals <- matrix(vals, nrow = dim(W_arr)[1L], ncol = dim(W_arr)[3L])
		within <- mean(apply(vals, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
		varies[k] <- is.finite(within) && within > 1e-8
	}
	nms[varies]
}

.lda_design_obs <- function(prep, t, obs, coef_names, symmetric = FALSE) {
	p <- length(coef_names)
	Xt <- matrix(1, nrow = sum(obs), ncol = p)
	colnames(Xt) <- coef_names
	if (sum(obs) == 0L) return(Xt)
	ij <- which(obs, arr.ind = TRUE)
	ii <- ij[, 1L]
	jj <- ij[, 2L]
	W_node <- NULL
	if (isTRUE(symmetric)) {
		W_node <- .lda_bind_node_matrices(
			.lda_node_matrix(prep, "row", t, dynamic = TRUE),
			.lda_node_matrix(prep, "col", t, dynamic = TRUE))
	}
	for (k in seq_along(coef_names)) {
		nm <- coef_names[k]
		if (identical(nm, "intercept")) next
		if (prep$p_dyad > 0L && nm %in% prep$x_names_dyad) {
			Xt[, k] <- prep$X[, , match(nm, prep$x_names_dyad), t][obs]
		} else if (!is.null(W_node) && isTRUE(symmetric) &&
		           nm %in% colnames(W_node)) {
			Xt[, k] <- W_node[ii, nm] + W_node[jj, nm]
		} else {
			W_row_t <- .lda_node_matrix(prep, "row", t, dynamic = TRUE)
			W_col_t <- .lda_node_matrix(prep, "col", t, dynamic = TRUE)
			if (!is.null(W_row_t) && nm %in% colnames(W_row_t)) {
				Xt[, k] <- W_row_t[ii, nm]
			} else if (!is.null(W_col_t) && nm %in% colnames(W_col_t)) {
				Xt[, k] <- W_col_t[jj, nm]
			} else {
				Xt[, k] <- 0
			}
		}
	}
	Xt
}

.lda_update_beta_path <- function(Z, W, prep, coef_path_old, a_mat, b_mat, Oarr,
                                  beta_dyn, lambda_beta, kind, rho_beta,
                                  symmetric = FALSE) {
	Tt <- prep$Tt
	p <- nrow(coef_path_old)
	coef_names <- rownames(coef_path_old)
	mask <- beta_dyn$mask
	idx <- .lda_coef_index(mask, Tt)
	H <- matrix(0, idx$n, idx$n)
	g <- numeric(idx$n)

	for (t in seq_len(Tt)) {
		obs <- is.finite(Z[, , t])
		if (!any(obs)) next
		target <- Z[, , t] - outer(a_mat[, t], b_mat[, t], "+") -
			.lda_o_slice(Oarr, t)
		y <- target[obs]
		wt <- W[, , t][obs]
		keep <- is.finite(wt) & wt > 0
		if (!any(keep)) next
		y <- y[keep]
		wt <- wt[keep]
		Xt <- .lda_design_obs(prep, t, obs, coef_names,
		                      symmetric = symmetric)
		Xt <- Xt[keep, , drop = FALSE]
		Xt_w <- Xt * sqrt(wt)
		XtX <- crossprod(Xt_w)
		Xty <- crossprod(Xt, wt * y)
		param_idx <- integer(p)
		for (k in seq_len(p)) {
			param_idx[k] <- if (mask[k]) idx$dynamic_map[k, t] else idx$static_map[k]
		}
		for (k in seq_len(p)) {
			g[param_idx[k]] <- g[param_idx[k]] + Xty[k]
			for (l in seq_len(p)) {
				H[param_idx[k], param_idx[l]] <-
					H[param_idx[k], param_idx[l]] + XtX[k, l]
			}
		}
	}

	if (length(idx$dynamic_idx) > 0L && lambda_beta > 0) {
		for (k in idx$dynamic_idx) {
			idk <- idx$dynamic_map[k, ]
			H[idk, idk] <- H[idk, idk] +
				.lda_path_penalty(Tt, lambda_beta, kind, rho_beta)
		}
	}
	H <- H + diag(1e-8, idx$n)
	theta <- .ae_safe_solve(H, g)

	out <- matrix(0, p, Tt, dimnames = dimnames(coef_path_old))
	for (k in seq_len(p)) {
		if (mask[k]) {
			out[k, ] <- theta[idx$dynamic_map[k, ]]
		} else {
			out[k, ] <- theta[idx$static_map[k]]
		}
	}
	out
}

.lda_xbeta_array <- function(prep, coef_path, symmetric = FALSE,
                             beta_dyn = NULL) {
	nr <- prep$nr
	nc <- prep$nc
	Tt <- prep$Tt
	p <- nrow(coef_path)
	cn <- rownames(coef_path)
	dyn_node_names <- .lda_dynamic_node_names(beta_dyn)
	out <- array(0, dim = c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		xt <- matrix(coef_path[1L, t], nr, nc)
		if (prep$p_dyad > 0L && p > 1L) {
			kmax <- min(prep$p_dyad, p - 1L)
			for (k in seq_len(kmax)) {
				xt <- xt + coef_path[k + 1L, t] * prep$X[, , k, t]
			}
		}
		if (!is.null(cn) && isTRUE(symmetric)) {
			W_node_static <- .lda_bind_node_matrices(prep$W_row, prep$W_col)
			W_node_dyn <- .lda_bind_node_matrices(
				.lda_node_matrix(prep, "row", t, dynamic = TRUE),
				.lda_node_matrix(prep, "col", t, dynamic = TRUE))
			W_node <- W_node_static
			if (!is.null(W_node)) {
				for (nm in intersect(colnames(W_node), cn)) {
					W_use <- if (nm %in% dyn_node_names &&
					             !is.null(W_node_dyn) &&
					             nm %in% colnames(W_node_dyn)) {
						W_node_dyn
					} else W_node_static
					g <- as.numeric(coef_path[nm, t] * W_use[, nm])
					xt <- xt + outer(g, g, "+")
				}
			}
		} else if (!is.null(cn)) {
			if (!is.null(prep$W_row)) {
				for (nm in intersect(colnames(prep$W_row), cn)) {
					W_row_t <- .lda_node_matrix(
						prep, "row", t,
						dynamic = nm %in% dyn_node_names)
					g <- as.numeric(coef_path[nm, t] * W_row_t[, nm])
					xt <- xt + matrix(g, nr, nc)
				}
			}
			if (!is.null(prep$W_col)) {
				for (nm in intersect(colnames(prep$W_col), cn)) {
					W_col_t <- .lda_node_matrix(
						prep, "col", t,
						dynamic = nm %in% dyn_node_names)
					g <- as.numeric(coef_path[nm, t] * W_col_t[, nm])
					xt <- xt + matrix(g, nr, nc, byrow = TRUE)
				}
			}
		}
		out[, , t] <- xt
	}
	out
}

.lda_update_dynamic_ab <- function(Z, W, Xb, Oarr, a_mat, b_mat,
                                   lambda_ab, rho_ab, kind = "ar1",
                                   transition_a = NULL,
                                   transition_b = NULL) {
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)

	for (i in seq_len(nr)) {
		rhs <- numeric(Tt)
		den <- numeric(Tt)
		for (t in seq_len(Tt)) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			Ot <- .lda_o_slice(Oarr, t)
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			ok_idx <- which(ok)[keep]
			wt <- wt[keep]
			base <- Z[i, ok_idx, t] - Xb[i, ok_idx, t] - Ot[i, ok_idx] -
				b_mat[ok_idx, t]
			rhs[t] <- sum(wt * base)
			den[t] <- sum(wt)
		}
		tw <- if (is.null(transition_a)) NULL else transition_a[i, ]
		a_mat[i, ] <- .lda_solve_path(
			den, rhs, lambda_ab, kind, rho_ab, transition_weights = tw)
	}
	for (j in seq_len(nc)) {
		rhs <- numeric(Tt)
		den <- numeric(Tt)
		for (t in seq_len(Tt)) {
			ok <- mask[, j, t]
			if (!any(ok)) next
			Ot <- .lda_o_slice(Oarr, t)
			wt <- W[ok, j, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			ok_idx <- which(ok)[keep]
			wt <- wt[keep]
			base <- Z[ok_idx, j, t] - Xb[ok_idx, j, t] - Ot[ok_idx, j] -
				a_mat[ok_idx, t]
			rhs[t] <- sum(wt * base)
			den[t] <- sum(wt)
		}
		tw <- if (is.null(transition_b)) NULL else transition_b[j, ]
		b_mat[j, ] <- .lda_solve_path(
			den, rhs, lambda_ab, kind, rho_ab, transition_weights = tw)
	}
	for (t in seq_len(Tt)) {
		a_mat[, t] <- a_mat[, t] - mean(a_mat[, t])
		b_mat[, t] <- b_mat[, t] - mean(b_mat[, t])
	}
	list(a = a_mat, b = b_mat)
}

.lda_update_dynamic_ab_symmetric <- function(Z, W, Xb, Oarr, a_mat,
                                             lambda_ab, rho_ab,
                                             kind = "ar1",
                                             transition_a = NULL) {
	n <- dim(Z)[1L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)

	for (i in seq_len(n)) {
		rhs <- numeric(Tt)
		den <- numeric(Tt)
		for (t in seq_len(Tt)) {
			Ot <- .lda_o_slice(Oarr, t)

			ok_r <- mask[i, , t]
			if (any(ok_r)) {
				wt <- W[i, ok_r, t]
				keep <- is.finite(wt) & wt > 0
				if (any(keep)) {
					idx <- which(ok_r)[keep]
					wt <- wt[keep]
					base <- Z[i, idx, t] - Xb[i, idx, t] - Ot[i, idx] -
						a_mat[idx, t]
					rhs[t] <- rhs[t] + sum(wt * base)
					den[t] <- den[t] + sum(wt)
				}
			}

			ok_c <- mask[, i, t]
			if (any(ok_c)) {
				wt <- W[ok_c, i, t]
				keep <- is.finite(wt) & wt > 0
				if (any(keep)) {
					idx <- which(ok_c)[keep]
					wt <- wt[keep]
					base <- Z[idx, i, t] - Xb[idx, i, t] - Ot[idx, i] -
						a_mat[idx, t]
					rhs[t] <- rhs[t] + sum(wt * base)
					den[t] <- den[t] + sum(wt)
				}
			}
		}
		tw <- if (is.null(transition_a)) NULL else transition_a[i, ]
		a_mat[i, ] <- .lda_solve_path(
			den, rhs, lambda_ab, kind, rho_ab, transition_weights = tw)
	}
	for (t in seq_len(Tt)) a_mat[, t] <- a_mat[, t] - mean(a_mat[, t])
	list(a = a_mat, b = a_mat)
}

.lda_update_static_ab <- function(Z, W, Xb, Oarr, symmetric = FALSE) {
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)
	base <- array(0, dim = dim(Z))
	for (t in seq_len(Tt)) {
		Ot <- .lda_o_slice(Oarr, t)
		base[, , t] <- ifelse(mask[, , t],
			W[, , t] * (Z[, , t] - Xb[, , t] - Ot), 0)
	}
	wsum <- apply(ifelse(mask, W, 0), c(1, 2), sum)
	base_rs <- apply(base, 1, sum)
	base_cs <- apply(base, 2, sum)
	cnt_row <- rowSums(wsum)
	cnt_col <- colSums(wsum)
	ab <- .ae_additive_solve(base_rs, base_cs, wsum, cnt_row, cnt_col,
	                         symmetric = symmetric)
	ab$a <- ab$a - mean(ab$a)
	ab$b <- ab$b - mean(ab$b)
	list(a = matrix(ab$a, nr, Tt), b = matrix(ab$b, nc, Tt))
}

.lda_update_static_lowrank <- function(Z, W, Xb, a_mat, b_mat, Omat, R,
                                       symmetric, lowrank_method) {
	if (R <= 0L) {
		return(list(Omat = matrix(0, dim(Z)[1L], dim(Z)[2L]),
		            U = NULL, V = NULL, L = NULL))
	}
	nr <- dim(Z)[1L]
	nc <- dim(Z)[2L]
	Tt <- dim(Z)[3L]
	mask <- is.finite(Z)
	wsum <- apply(ifelse(mask, W, 0), c(1, 2), sum)
	res_sum <- matrix(0, nr, nc)
	for (t in seq_len(Tt)) {
		Wt <- W[, , t]
		Wt[!mask[, , t] | !is.finite(Wt)] <- 0
		res_t <- Z[, , t] - Xb[, , t] - outer(a_mat[, t], b_mat[, t], "+")
		res_t[!mask[, , t]] <- 0
		res_sum <- res_sum + Wt * res_t
	}
	obs_dyad <- wsum > 0
	Rbar <- matrix(0, nr, nc)
	Rbar[obs_dyad] <- res_sum[obs_dyad] / wsum[obs_dyad]
	.ae_lowrank_update(Omat, Rbar, wsum, obs_dyad, R, symmetric,
	                   lowrank_method)
}

# marginal level prior for the multiplicative terms. the temporal
# difference penalty (lambda_uv / lambda_g) only constrains
# u_{t} - rho * u_{t-1}: the common level of a whole path lies in (or near)
# the null space of that penalty, so with only a 1e-6 numerical ridge a
# binary / poisson IRLS fit can push the factor paths to arbitrarily large
# levels while improving working deviance (quasi-separation runaway).
# .lda_level_prec is the ALS counterpart of the mcmc sampler's per-slice
# N(0, s2 * I) factor prior expressed in the working scale of the weighted
# least-squares blocks: precision 1 equals prior variance 1 on the probit /
# log latent scale for the IRLS families and prior variance ve for the
# weight-1 normal working response. it is kept fixed (not adapted to the
# factor magnitudes) so exact block minimization preserves ALS descent and
# the level prior cannot chase a shrinking factor scale into collapse.
.lda_level_prec <- 1

.lda_update_uv_path_block <- function(target, mask, W, other_cube, margin,
                                      lambda_uv, rho_uv,
                                      transition_weights = NULL,
                                      ridge = 1e-6, level_prec = 0) {
	d <- dim(other_cube)
	R <- d[2L]
	Tt <- d[3L]
	n_path <- if (identical(margin, "row")) dim(target)[1L] else dim(target)[2L]
	out <- array(0, dim = c(n_path, R, Tt))

	for (idx_path in seq_len(n_path)) {
		tw <- if (is.null(transition_weights)) rep(1, Tt) else
			transition_weights[idx_path, ]
		P <- kronecker(.lda_path_penalty_weighted(Tt, lambda_uv, tw, rho_uv),
		               diag(R))
		H <- P
		g <- numeric(R * Tt)
		for (t in seq_len(Tt)) {
			block <- ((t - 1L) * R + 1L):(t * R)
			if (identical(margin, "row")) {
				ok <- mask[idx_path, , t]
				if (!any(ok)) next
				wt <- W[idx_path, ok, t]
				keep <- is.finite(wt) & wt > 0
				if (!any(keep)) next
				Xt <- matrix(other_cube[ok, , t], ncol = R)
				yt <- target[idx_path, ok, t]
			} else {
				ok <- mask[, idx_path, t]
				if (!any(ok)) next
				wt <- W[ok, idx_path, t]
				keep <- is.finite(wt) & wt > 0
				if (!any(keep)) next
				Xt <- matrix(other_cube[ok, , t], ncol = R)
				yt <- target[ok, idx_path, t]
			}
			Xt <- Xt[keep, , drop = FALSE]
			yt <- yt[keep]
			wt <- wt[keep]
			H[block, block] <- H[block, block] + crossprod(Xt * sqrt(wt))
			g[block] <- g[block] + as.numeric(crossprod(Xt, wt * yt))
		}
		# level prior on every time block plus the numerical floor
		H <- H + diag(level_prec + ridge, R * Tt)
		sol <- .ae_safe_solve(H, g)
		sol[!is.finite(sol)] <- 0
		for (t in seq_len(Tt)) {
			block <- ((t - 1L) * R + 1L):(t * R)
			out[idx_path, , t] <- sol[block]
		}
	}
	out
}

.lda_update_dynamic_uv_directed <- function(Z, W, Xb, a_mat, b_mat,
                                            U_cube, V_cube, lambda_uv,
                                            rho_uv, symmetric = FALSE,
                                            transition_u = NULL,
                                            transition_v = NULL,
                                            align = TRUE,
                                            level_prec = 0) {
	Tt <- dim(Z)[3L]
	udn <- dimnames(U_cube)
	vdn <- dimnames(V_cube)
	target <- array(0, dim = dim(Z))
	for (t in seq_len(Tt)) {
		target[, , t] <- Z[, , t] - Xb[, , t] -
			outer(a_mat[, t], b_mat[, t], "+")
	}
	mask <- is.finite(Z)
	U_cube <- .lda_update_uv_path_block(
		target, mask, W, V_cube, margin = "row",
		lambda_uv = lambda_uv, rho_uv = rho_uv,
		transition_weights = transition_u,
		level_prec = level_prec)
	dimnames(U_cube) <- udn
	V_cube <- .lda_update_uv_path_block(
		target, mask, W, U_cube, margin = "col",
		lambda_uv = lambda_uv, rho_uv = rho_uv,
		transition_weights = transition_v,
		level_prec = level_prec)
	dimnames(V_cube) <- vdn
	if (isTRUE(align)) {
		uv <- .lda_align_uv_cube(U_cube, V_cube)
		U_cube <- uv$U
		V_cube <- uv$V
	}
	Oarr <- .lda_uv_to_o_array(U_cube, V_cube)
	if (isTRUE(symmetric)) {
		sf <- .lda_symmetric_o_to_uv_cube(
			Oarr, dim(U_cube)[2L],
			row_names = dimnames(U_cube)[[1L]],
			dim_names = dimnames(U_cube)[[2L]],
			time_names = dimnames(U_cube)[[3L]])
		U_cube <- sf$U
		V_cube <- sf$V
		Oarr <- sf$Oarr
	}
	list(U = U_cube, V = V_cube, Oarr = Oarr)
}

.lda_update_dynamic_uv_bip_g <- function(Z, W, Xb, a_mat, b_mat,
                                         U_cube, V_cube, G,
                                         lambda_uv, rho_uv,
                                         transition_u = NULL,
                                         transition_v = NULL,
                                         align = TRUE,
                                         level_prec = 0) {
	Tt <- dim(Z)[3L]
	RA <- dim(U_cube)[2L]
	RB <- dim(V_cube)[2L]
	dg <- dim(G)
	target <- array(0, dim = dim(Z))
	for (t in seq_len(Tt)) {
		target[, , t] <- Z[, , t] - Xb[, , t] -
			outer(a_mat[, t], b_mat[, t], "+")
	}
	mask <- is.finite(Z)
	V_eff <- array(0, dim = c(dim(V_cube)[1L], RA, Tt),
	               dimnames = list(dimnames(V_cube)[[1L]],
	                               dimnames(U_cube)[[2L]],
	                               dimnames(V_cube)[[3L]]))
	for (t in seq_len(Tt)) {
		Vt <- matrix(V_cube[, , t], nrow = dim(V_cube)[1L], ncol = RB)
		Gt <- if (length(dg) == 3L) {
			matrix(G[, , t], nrow = RA, ncol = RB)
		} else G
		V_eff[, , t] <- Vt %*% t(Gt)
	}
	udn <- dimnames(U_cube)
	vdn <- dimnames(V_cube)
	gdn <- dimnames(G)
	U_cube <- .lda_update_uv_path_block(
		target, mask, W, V_eff, margin = "row",
		lambda_uv = lambda_uv, rho_uv = rho_uv,
		transition_weights = transition_u,
		level_prec = level_prec)
	dimnames(U_cube) <- udn
	U_eff <- array(0, dim = c(dim(U_cube)[1L], RB, Tt),
	               dimnames = list(dimnames(U_cube)[[1L]],
	                               dimnames(V_cube)[[2L]],
	                               dimnames(U_cube)[[3L]]))
	for (t in seq_len(Tt)) {
		Ut <- matrix(U_cube[, , t], nrow = dim(U_cube)[1L], ncol = RA)
		Gt <- if (length(dg) == 3L) {
			matrix(G[, , t], nrow = RA, ncol = RB)
		} else G
		U_eff[, , t] <- Ut %*% Gt
	}
	V_cube <- .lda_update_uv_path_block(
		target, mask, W, U_eff, margin = "col",
		lambda_uv = lambda_uv, rho_uv = rho_uv,
		transition_weights = transition_v,
		level_prec = level_prec)
	dimnames(V_cube) <- vdn
	if (isTRUE(align) && length(dg) == 3L) {
		uvg <- .lda_align_uvg_cube(U_cube, V_cube, G)
		U_cube <- uvg$U
		V_cube <- uvg$V
		G <- uvg$G
	}
	dimnames(G) <- gdn
	list(U = U_cube, V = V_cube, G = G,
	     Oarr = .lda_uvg_to_o_array(U_cube, V_cube, G))
}

.lda_update_static_uv_bip_g <- function(target, mask, W, U, V, G,
                                        ridge = 1e-6, level_prec = 0) {
	dg <- dim(G)
	RA <- ncol(U)
	RB <- ncol(V)
	Tt <- dim(target)[3L]
	udn <- dimnames(U)
	vdn <- dimnames(V)

	for (i in seq_len(nrow(U))) {
		H <- diag(level_prec + ridge, RA)
		g <- numeric(RA)
		for (t in seq_len(Tt)) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			jj <- which(ok)[keep]
			wt <- wt[keep]
			Gt <- if (length(dg) == 3L) matrix(G[, , t], RA, RB) else G
			Xt <- matrix(0, length(jj), RA)
			for (m in seq_along(jj)) Xt[m, ] <- as.numeric(V[jj[m], ] %*% t(Gt))
			Xtw <- Xt * sqrt(wt)
			H <- H + crossprod(Xtw)
			g <- g + as.numeric(crossprod(Xt, wt * target[i, jj, t]))
		}
		U[i, ] <- .ae_safe_solve(H, g)
	}
	dimnames(U) <- udn

	for (j in seq_len(nrow(V))) {
		H <- diag(level_prec + ridge, RB)
		g <- numeric(RB)
		for (t in seq_len(Tt)) {
			ok <- mask[, j, t]
			if (!any(ok)) next
			wt <- W[ok, j, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			ii <- which(ok)[keep]
			wt <- wt[keep]
			Gt <- if (length(dg) == 3L) matrix(G[, , t], RA, RB) else G
			Xt <- matrix(0, length(ii), RB)
			for (m in seq_along(ii)) Xt[m, ] <- as.numeric(U[ii[m], ] %*% Gt)
			Xtw <- Xt * sqrt(wt)
			H <- H + crossprod(Xtw)
			g <- g + as.numeric(crossprod(Xt, wt * target[ii, j, t]))
		}
		V[j, ] <- .ae_safe_solve(H, g)
	}
	dimnames(V) <- vdn
		list(U = U, V = V,
		     Oarr = .lda_uvg_to_o_array(U, V, G, Tt = dim(target)[3L]))
}

.lda_update_static_g_bip <- function(target, mask, W, U, V, G,
                                     ridge = 1e-6) {
	du <- dim(U)
	dv <- dim(V)
	RA <- if (length(du) == 3L) du[2L] else ncol(U)
	RB <- if (length(dv) == 3L) dv[2L] else ncol(V)
	Tt <- dim(target)[3L]
	p <- RA * RB
	H <- diag(ridge, p)
	g <- numeric(p)
	for (t in seq_len(Tt)) {
		Ut <- if (length(du) == 3L) matrix(U[, , t], nrow = du[1L], ncol = RA) else U
		Vt <- if (length(dv) == 3L) matrix(V[, , t], nrow = dv[1L], ncol = RB) else V
		for (i in seq_len(nrow(target))) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			jj <- which(ok)[keep]
			wt <- wt[keep]
			yt <- target[i, jj, t]
			for (m in seq_along(jj)) {
				x <- as.vector(outer(Ut[i, ], Vt[jj[m], ]))
				H <- H + wt[m] * tcrossprod(x)
				g <- g + wt[m] * x * yt[m]
			}
		}
	}
	out <- matrix(.ae_safe_solve(H, g), RA, RB)
	out[!is.finite(out)] <- 0
	dimnames(out) <- dimnames(G)
	out
}

.lda_update_dynamic_g_path <- function(target, mask, W, U, V, G_cube,
                                       lambda_g, rho_g, ridge = 1e-6,
                                       level_prec = 0) {
	du <- dim(U)
	dv <- dim(V)
	dg <- dim(G_cube)
	Tt <- dim(target)[3L]
	RA <- if (length(du) == 3L) du[2L] else ncol(U)
	RB <- if (length(dv) == 3L) dv[2L] else ncol(V)
	p <- RA * RB
	P <- kronecker(.lda_path_penalty(Tt, lambda_g, "ar1", rho_g), diag(p))
	H <- P
	g <- numeric(p * Tt)
	for (t in seq_len(Tt)) {
		block <- ((t - 1L) * p + 1L):(t * p)
		Ut <- if (length(du) == 3L) {
			matrix(U[, , t], nrow = du[1L], ncol = RA)
		} else U
		Vt <- if (length(dv) == 3L) {
			matrix(V[, , t], nrow = dv[1L], ncol = RB)
		} else V
		for (i in seq_len(nrow(target))) {
			ok <- mask[i, , t]
			if (!any(ok)) next
			wt <- W[i, ok, t]
			keep <- is.finite(wt) & wt > 0
			if (!any(keep)) next
			jj <- which(ok)[keep]
			wt <- wt[keep]
			yt <- target[i, jj, t]
			for (m in seq_along(jj)) {
				x <- as.vector(outer(Ut[i, ], Vt[jj[m], ]))
				H[block, block] <- H[block, block] + wt[m] * tcrossprod(x)
				g[block] <- g[block] + wt[m] * x * yt[m]
			}
		}
	}
	# level prior on every time block plus the numerical floor
	H <- H + diag(level_prec + ridge, p * Tt)
	theta <- .ae_safe_solve(H, g)
	theta[!is.finite(theta)] <- 0
	out <- array(0, dim = c(RA, RB, Tt), dimnames = dimnames(G_cube))
	for (t in seq_len(Tt)) {
		block <- ((t - 1L) * p + 1L):(t * p)
		out[, , t] <- matrix(theta[block], RA, RB)
	}
	if (!is.null(dg) && length(dg) == 3L) dimnames(out) <- dimnames(G_cube)
	out
}

.lda_t_transition_weights <- function(cube, rho, nu = 4, scale = NULL) {
	d <- dim(cube)
	dn <- dimnames(cube)
	out <- matrix(1, d[1L], d[3L],
	              dimnames = if (!is.null(dn)) list(dn[[1L]], dn[[3L]]) else NULL)
	if (d[3L] < 2L) return(out)
	innov <- cube[, , -1L, drop = FALSE] -
		rho * cube[, , -d[3L], drop = FALSE]
	ss <- apply(innov^2, c(1, 3), sum)
	if (is.null(scale) || !is.finite(scale) || scale <= 0) {
		scale <- sqrt(stats::median(ss[is.finite(ss) & ss > 0], na.rm = TRUE) / d[2L])
		if (!is.finite(scale) || scale <= 0) scale <- 1
	}
	out[, -1L] <- (nu + d[2L]) / (nu + ss / max(scale^2, 1e-8))
	out[!is.finite(out) | out <= 0] <- 1
	pmin(out, 10)
}

.lda_xdyad_from_fit <- function(fit) {
	if (is.null(fit$X) || length(dim(fit$X)) < 4L || dim(fit$X)[3L] == 0L) {
		return(NULL)
	}
	lapply(seq_len(dim(fit$X)[4L]), function(t) fit$X[, , , t, drop = FALSE])
}

.lda_fit_surface_delta <- function(base, other) {
	if (is.null(base$EZ) || is.null(other$EZ)) return(NA_real_)
	max(vapply(seq_len(base$n_time), function(t) {
		max(abs(base$EZ[[t]] - other$EZ[[t]]), na.rm = TRUE)
	}, numeric(1)), na.rm = TRUE)
}

.lda_dynamic_stability <- function(fit, starts = c("quick", "validation"),
                                   seed = 6886, call_args = list()) {
	starts <- match.arg(starts)
	n_start <- if (identical(starts, "quick")) 2L else 4L
	refits <- vector("list", n_start)
	for (s in seq_len(n_start)) {
		args <- call_args
		args$seed <- seed + 1000L + s
		args$verbose <- FALSE
		args$bootstrap <- 0L
		args$stability <- "none"
		args$start_jitter <- 0.05
		args$call <- NULL
		refits[[s]] <- tryCatch(do.call(lame_dynamic_als, args),
		                        error = function(e) e)
	}
	ok <- vapply(refits, function(x) inherits(x, "lame_dynamic_als"), logical(1))
	valid <- refits[ok]
	surface_delta <- if (length(valid) > 0L) {
		vapply(valid, function(x) .lda_fit_surface_delta(fit, x), numeric(1))
	} else numeric(0)
	coef_delta <- if (length(valid) > 0L) {
		vapply(valid, function(x) max(abs(x$coef_path - fit$coef_path), na.rm = TRUE),
		       numeric(1))
	} else numeric(0)
	list(
		preset = starts,
		n_start = n_start,
		n_success = sum(ok),
		converged = vapply(valid, function(x) isTRUE(x$converged), logical(1)),
		max_surface_delta = if (length(surface_delta)) max(surface_delta, na.rm = TRUE) else NA_real_,
		max_coef_delta = if (length(coef_delta)) max(coef_delta, na.rm = TRUE) else NA_real_,
		refits = valid,
		failures = refits[!ok])
}

.lda_dynamic_bootstrap <- function(fit, R_boot, type, block_length, seed,
                                   call_args = list(), verbose = FALSE) {
	if (R_boot <= 0L) return(NULL)
	if (!identical(type, "parametric")) {
		cli::cli_warn(c(
			"{.arg bootstrap_type = \"block\"} is not yet available for dynamic ALS.",
			"i" = "Using parametric dynamic refits instead."))
		type <- "parametric"
	}
	sims <- simulate(fit, nsim = R_boot, seed = seed)
	obs_mask <- is.finite(fit$Y)
	mask_sim <- function(Ysim) {
		if (is.list(Ysim) && !is.data.frame(Ysim)) {
			out <- Ysim
			for (t in seq_along(out)) {
				out[[t]][!obs_mask[, , t]] <- NA_real_
				dimnames(out[[t]]) <- dimnames(fit$Y)[1:2]
				if (isTRUE(fit$changing_composition)) {
					row_keep <- (fit$row_presence %||% fit$actor_presence)[, t]
					col_keep <- (fit$col_presence %||% fit$actor_presence)[, t]
					out[[t]] <- out[[t]][row_keep, col_keep, drop = FALSE]
				}
			}
			names(out) <- dimnames(fit$Y)[[3L]]
			return(out)
		}
		if (is.array(Ysim) && length(dim(Ysim)) == 3L) {
			Ysim[!obs_mask] <- NA_real_
			dimnames(Ysim) <- dimnames(fit$Y)
			if (isTRUE(fit$changing_composition)) {
				out <- lapply(seq_len(dim(Ysim)[3L]), function(t) {
					row_keep <- (fit$row_presence %||% fit$actor_presence)[, t]
					col_keep <- (fit$col_presence %||% fit$actor_presence)[, t]
					Ysim[row_keep, col_keep, t, drop = FALSE][, , 1L]
				})
				names(out) <- dimnames(fit$Y)[[3L]]
				return(out)
			}
			return(Ysim)
		}
		Ysim
	}
	refits <- vector("list", R_boot)
	for (b in seq_len(R_boot)) {
		args <- call_args
		args$Y <- mask_sim(sims$Y[[b]])
		args$seed <- (if (is.null(seed)) 6886 else seed) + b
		args$verbose <- FALSE
		args$bootstrap <- 0L
		args$stability <- "none"
		args$call <- NULL
		refits[[b]] <- tryCatch(do.call(lame_dynamic_als, args),
		                        error = function(e) e)
		if (isTRUE(verbose) && b %% 10L == 0L) {
			cli::cli_text("Dynamic bootstrap refit {b}/{R_boot}")
		}
	}
	ok <- vapply(refits, function(x) inherits(x, "lame_dynamic_als"), logical(1))
	valid <- refits[ok]
	coef_paths <- if (length(valid) > 0L) {
		array(unlist(lapply(valid, function(x) x$coef_path), use.names = FALSE),
		      dim = c(nrow(fit$coef_path), ncol(fit$coef_path), length(valid)),
		      dimnames = list(rownames(fit$coef_path), colnames(fit$coef_path), NULL))
	} else NULL
	list(
		type = type,
		R = R_boot,
		n_success = sum(ok),
		coef_paths = coef_paths,
		converged = vapply(valid, function(x) isTRUE(x$converged), logical(1)),
		refits = valid,
		failures = refits[!ok])
}

.lda_eta <- function(prep, coef_path, a_mat, b_mat, Oarr, symmetric = FALSE,
                     beta_dyn = NULL) {
	Xb <- .lda_xbeta_array(prep, coef_path, symmetric = symmetric,
	                       beta_dyn = beta_dyn)
	E <- vector("list", prep$Tt)
	for (t in seq_len(prep$Tt)) {
		E[[t]] <- Xb[, , t] + outer(a_mat[, t], b_mat[, t], "+") +
			.lda_o_slice(Oarr, t)
	}
	E
}

.lda_objective <- function(Z, W, prep, coef_path, a_mat, b_mat, Oarr,
	                           beta_dyn, lambda_beta, beta_kind, rho_beta,
	                           dynamic_ab, lambda_ab, rho_ab,
	                           transition_a = NULL, transition_b = NULL,
	                           symmetric = FALSE,
	                           dynamic_uv = FALSE, U = NULL, V = NULL,
                           lambda_uv = 0, rho_uv = 1,
                           transition_u = NULL, transition_v = NULL,
                           dynamic_G = FALSE, G_cube = NULL,
                           lambda_g = 0, rho_g = 1,
                           level_prec_uv = 0, level_prec_g = 0) {
	E <- .lda_eta(prep, coef_path, a_mat, b_mat, Oarr,
	              symmetric = symmetric, beta_dyn = beta_dyn)
	sse <- 0
	for (t in seq_len(prep$Tt)) {
		ok <- is.finite(Z[, , t])
		r <- Z[, , t] - E[[t]]
		sse <- sse + sum(W[, , t][ok] * r[ok]^2)
	}
	pen_ab <- 0
		if (isTRUE(dynamic_ab) && lambda_ab > 0 && prep$Tt > 1L) {
			for (i in seq_len(nrow(a_mat))) {
				d <- a_mat[i, -1L] - rho_ab * a_mat[i, -prep$Tt]
				wa <- if (is.null(transition_a)) rep(1, prep$Tt - 1L) else
					pmax(0, transition_a[i, -1L])
				pen_ab <- pen_ab + lambda_ab * sum(wa * d^2)
			}
			if (!isTRUE(symmetric)) {
				for (j in seq_len(nrow(b_mat))) {
					d <- b_mat[j, -1L] - rho_ab * b_mat[j, -prep$Tt]
					wb <- if (is.null(transition_b)) rep(1, prep$Tt - 1L) else
						pmax(0, transition_b[j, -1L])
					pen_ab <- pen_ab + lambda_ab * sum(wb * d^2)
				}
			}
		}
	pen_beta <- 0
	if (isTRUE(beta_dyn$any) && lambda_beta > 0) {
		P <- .lda_path_penalty(prep$Tt, lambda_beta, beta_kind, rho_beta)
		for (k in which(beta_dyn$mask)) {
			pen_beta <- pen_beta + as.numeric(t(coef_path[k, ]) %*% P %*% coef_path[k, ])
		}
	}
	pen_uv <- 0
	if (isTRUE(dynamic_uv) && lambda_uv > 0 && prep$Tt > 1L &&
	    !is.null(U) && !is.null(V) && length(dim(U)) == 3L) {
		for (i in seq_len(dim(U)[1L])) {
			d <- U[i, , -1L, drop = FALSE] -
				rho_uv * U[i, , -prep$Tt, drop = FALSE]
			wu <- if (is.null(transition_u)) rep(1, prep$Tt - 1L) else
				pmax(0, transition_u[i, -1L])
			pen_uv <- pen_uv + lambda_uv * sum(wu * apply(d^2, 3, sum))
		}
		for (j in seq_len(dim(V)[1L])) {
			d <- V[j, , -1L, drop = FALSE] -
				rho_uv * V[j, , -prep$Tt, drop = FALSE]
			wv <- if (is.null(transition_v)) rep(1, prep$Tt - 1L) else
				pmax(0, transition_v[j, -1L])
			pen_uv <- pen_uv + lambda_uv * sum(wv * apply(d^2, 3, sum))
		}
	}
	pen_g <- 0
	if (isTRUE(dynamic_G) && lambda_g > 0 && prep$Tt > 1L &&
	    !is.null(G_cube) && length(dim(G_cube)) == 3L) {
		for (r in seq_len(dim(G_cube)[1L])) for (c in seq_len(dim(G_cube)[2L])) {
			d <- G_cube[r, c, -1L] - rho_g * G_cube[r, c, -prep$Tt]
			pen_g <- pen_g + lambda_g * sum(d^2)
		}
	}
	# fixed-precision marginal level prior on the multiplicative terms:
	# matches the diag(level_prec) added to the factor block solves, so the
	# monitored objective is the exact penalized least-squares / IRLS
	# working objective the block updates minimize
	pen_level <- 0
	if (level_prec_uv > 0 && !is.null(U) && !is.null(V)) {
		pen_level <- pen_level + level_prec_uv * (sum(U^2) + sum(V^2))
	}
	if (level_prec_g > 0 && !is.null(G_cube) &&
	    length(dim(G_cube)) == 3L) {
		pen_level <- pen_level + level_prec_g * sum(G_cube^2)
	}
	sse + pen_ab + pen_beta + pen_uv + pen_g + pen_level
}

.lda_default_lambdas <- function(prep, lambda_ab, lambda_beta, lambda_uv) {
	mask <- is.finite(prep$Yarr)
	row_deg <- numeric(0)
	for (t in seq_len(prep$Tt)) {
		row_deg <- c(row_deg, rowSums(mask[, , t]), colSums(mask[, , t]))
	}
	deg_pos <- row_deg[row_deg > 0]
	deg_med <- if (length(deg_pos) > 0L) stats::median(deg_pos, na.rm = TRUE) else 1
	if (is.null(lambda_ab)) {
		lambda_ab <- max(1, deg_med / 4)
	}
	if (is.null(lambda_beta)) {
		nobs_t <- vapply(seq_len(prep$Tt), function(t) sum(mask[, , t]), numeric(1))
		nobs_pos <- nobs_t[nobs_t > 0]
		nobs_med <- if (length(nobs_pos) > 0L) stats::median(nobs_pos, na.rm = TRUE) else 1
		lambda_beta <- max(1, sqrt(nobs_med))
	}
	if (is.null(lambda_uv)) {
		lambda_uv <- max(1, deg_med / 4)
	}
	list(ab = lambda_ab, beta = lambda_beta, uv = lambda_uv)
}

.lda_decompose_margin_node_paths <- function(coef_path, e_mat, W, active,
                                             time_names,
                                             intercept_multiplier = 1,
                                             W_arr = NULL,
                                             presence = NULL) {
	if (is.null(W) || ncol(W) == 0L) {
		return(list(coef_path = coef_path, e_mat = e_mat,
		            beta_names = character(0)))
	}
	Tt <- ncol(e_mat)
	n_actor <- nrow(e_mat)
	node_names <- colnames(W)
	dynamic_names <- intersect(node_names, rownames(coef_path))
	static_names <- setdiff(node_names, dynamic_names)
	beta_names <- character(0)

		if (length(dynamic_names) > 0L) {
			for (t in seq_len(Tt)) {
				W_t <- if (!is.null(W_arr)) {
				tmp <- W_arr[, , t, drop = FALSE]
				tmp <- matrix(tmp, nrow = dim(W_arr)[1L],
				              ncol = dim(W_arr)[2L])
				dimnames(tmp) <- list(NULL, dimnames(W_arr)[[2L]])
				tmp
				} else W
				W_dyn <- W_t[, dynamic_names, drop = FALSE]
				active_t <- if (is.matrix(active) && ncol(active) >= t) {
					as.numeric(active[, t])
				} else active
				dec <- .ae_decompose(e_mat[, t], W_dyn, active_t)
				e_mat[, t] <- dec$resid
			coef_path[1L, t] <- coef_path[1L, t] +
				intercept_multiplier * dec$intercept
			coef_path[names(dec$beta), t] <-
				coef_path[names(dec$beta), t] + dec$beta
		}
		beta_names <- c(beta_names, dynamic_names)
	}

		if (length(static_names) > 0L) {
			W_static <- W[, static_names, drop = FALSE]
			active_static <- if (is.matrix(active)) {
				as.numeric(rowSums(active) > 0)
			} else active
			dec <- .ae_decompose(
				.lda_presence_row_means(e_mat, presence),
				W_static, active_static)
		if (length(dec$beta) > 0L) {
			node_level <- as.numeric(dec$intercept + W_static %*% dec$beta)
			e_mat <- e_mat - matrix(node_level, n_actor, Tt)
			coef_path[1L, ] <- coef_path[1L, ] +
				intercept_multiplier * dec$intercept
			node_path <- matrix(rep(dec$beta, Tt), nrow = length(dec$beta),
			                    dimnames = list(names(dec$beta), time_names))
			coef_path <- rbind(coef_path, node_path)
			beta_names <- c(beta_names, names(dec$beta))
		}
	}

	list(coef_path = coef_path, e_mat = e_mat,
	     beta_names = unique(beta_names))
}

.lda_decompose_node_paths <- function(prep, symmetric, coef_path, a_mat, b_mat) {
	has_node <- !is.null(prep$W_row) || !is.null(prep$W_col)
	if (!isTRUE(has_node)) {
		return(list(coef_path = coef_path, a_mat = a_mat, b_mat = b_mat,
		            beta_names = character(0)))
	}

	Tt <- prep$Tt
	obs_count <- apply(!is.na(prep$Yarr), c(1, 2), sum)
	obs_row_time <- matrix(FALSE, prep$nr, Tt)
	obs_col_time <- matrix(FALSE, prep$nc, Tt)
	for (t in seq_len(Tt)) {
		yt <- prep$Yarr[, , t]
		obs_row_time[, t] <- rowSums(!is.na(yt)) > 0
		obs_col_time[, t] <- colSums(!is.na(yt)) > 0
	}
	dimnames(obs_row_time) <- list(prep$row_names, prep$time_names)
	dimnames(obs_col_time) <- list(prep$col_names, prep$time_names)
	if (isTRUE(symmetric)) obs_count <- obs_count + t(obs_count)
	active_row <- as.numeric(rowSums(obs_count) > 0)
	active_col <- as.numeric(colSums(obs_count) > 0)
	active_row_time <- obs_row_time
	active_col_time <- obs_col_time

	if (isTRUE(symmetric)) {
		W_node <- prep$W_row
		if (!is.null(prep$W_col)) {
			W_node <- if (is.null(W_node)) prep$W_col else cbind(W_node, prep$W_col)
			}
			W_node_arr <- .lda_bind_node_arrays(prep$W_row_arr, prep$W_col_arr)
			active_sym_time <- obs_row_time | obs_col_time
			dec <- .lda_decompose_margin_node_paths(
				coef_path, a_mat, W_node, active_sym_time, prep$time_names,
				intercept_multiplier = 2, W_arr = W_node_arr,
				presence = prep$row_presence)
		coef_path <- dec$coef_path
		a_mat <- dec$e_mat
		b_mat <- a_mat
		return(list(coef_path = coef_path, a_mat = a_mat, b_mat = b_mat,
		            beta_names = dec$beta_names))
	}

	dec_r <- .lda_decompose_margin_node_paths(
		coef_path, a_mat, prep$W_row, active_row_time, prep$time_names,
		W_arr = prep$W_row_arr, presence = prep$row_presence)
	coef_path <- dec_r$coef_path
	a_mat <- dec_r$e_mat
	dec_c <- .lda_decompose_margin_node_paths(
		coef_path, b_mat, prep$W_col, active_col_time, prep$time_names,
		W_arr = prep$W_col_arr, presence = prep$col_presence)
	coef_path <- dec_c$coef_path
	b_mat <- dec_c$e_mat
	list(coef_path = coef_path, a_mat = a_mat, b_mat = b_mat,
	     beta_names = unique(c(dec_r$beta_names, dec_c$beta_names)))
}

.lda_as_fit <- function(prep, family, mode, symmetric, R, R_row, R_col,
                        call, coef_path,
                        beta_dyn, dynamic_beta_kind, a_mat, b_mat, Omat, Oarr,
                        U, V, L, G_static, G_cube,
                        objective_trace, converged, iterations,
                        lambda_ab, lambda_beta, lambda_uv,
                        lambda_g, rho_ab, rho_beta, rho_uv, rho_g,
                        dynamic_ab, dynamic_uv, dynamic_uv_kind, dynamic_G,
                        lowrank_method, link = NULL, Zwork = NULL,
                        obs_weights = NULL, irls_trace = NULL,
                        stability = NULL,
                        transition_u = NULL, transition_v = NULL,
                        convergence_trace = NULL, tol = NULL) {
	nr <- prep$nr
	nc <- prep$nc
	Tt <- prep$Tt
	if (is.null(Oarr)) Oarr <- .lda_o_array(Omat, Tt)
	objective_diff <- diff(objective_trace)
	objective_increase_tol <- 1e-10 * (abs(utils::head(objective_trace, -1L)) + 1)
	objective_increases <- if (length(objective_diff) > 0L) {
		sum(objective_diff > objective_increase_tol, na.rm = TRUE)
	} else {
		0L
	}
	max_objective_increase <- if (length(objective_diff) > 0L) {
		max(c(0, objective_diff), na.rm = TRUE)
	} else {
		0
	}
	objective_label <- if (identical(family, "normal")) {
		"penalized least-squares objective"
	} else {
		"penalized IRLS working objective"
	}
	Omat_store <- if (isTRUE(dynamic_uv) || isTRUE(dynamic_G)) {
		apply(Oarr, c(1, 2), mean)
	} else Omat
	dec <- .lda_decompose_node_paths(prep, symmetric, coef_path, a_mat, b_mat)
	coef_path <- dec$coef_path
	a_mat <- dec$a_mat
	b_mat <- dec$b_mat
	EZ <- .lda_eta(prep, coef_path, a_mat, b_mat, Oarr,
	               symmetric = symmetric, beta_dyn = beta_dyn)
	fitted_list <- vector("list", Tt)
	resid_list <- vector("list", Tt)
	Zresid <- if (is.null(Zwork)) prep$Yarr else Zwork
	ve_num <- 0
	ve_den <- 0
	for (t in seq_len(Tt)) {
		fitted_list[[t]] <- .ae_link_inverse(EZ[[t]], family, link)
		if (!prep$bip && nr == nc) {
			diag(EZ[[t]]) <- NA_real_
			diag(fitted_list[[t]]) <- NA_real_
		}
		dimnames(EZ[[t]]) <- list(prep$row_names, prep$col_names)
		dimnames(fitted_list[[t]]) <- list(prep$row_names, prep$col_names)
		ok <- is.finite(Zresid[, , t])
		r <- Zresid[, , t] - EZ[[t]]
		resid_list[[t]] <- r
		ve_num <- ve_num + sum(r[ok]^2)
		ve_den <- ve_den + sum(ok)
	}
	names(EZ) <- prep$time_names
	names(fitted_list) <- prep$time_names
	names(resid_list) <- prep$time_names
	ve <- if (ve_den > 0) ve_num / ve_den else NA_real_

	dimnames(a_mat) <- list(prep$row_names, prep$time_names)
	dimnames(b_mat) <- list(prep$col_names, prep$time_names)
	a_bar <- .lda_presence_row_means(a_mat, prep$row_presence)
	b_bar <- .lda_presence_row_means(b_mat, prep$col_presence)
	va <- if (nr > 1L) mean(.lda_presence_col_vars(a_mat, prep$row_presence), na.rm = TRUE) else NA_real_
	vb <- if (symmetric) NA_real_ else
		if (nc > 1L) mean(.lda_presence_col_vars(b_mat, prep$col_presence), na.rm = TRUE) else NA_real_
	cab <- if (symmetric || prep$bip) NA_real_ else
		mean(.lda_presence_col_cov(a_mat, b_mat, prep$row_presence, prep$col_presence),
		     na.rm = TRUE)
	rho <- 0
	rho_path <- NULL
	if (!prep$bip && !symmetric) {
		rho_path <- rep(NA_real_, Tt)
		names(rho_path) <- prep$time_names
		e_ij <- numeric(0)
		e_ji <- numeric(0)
		for (t in seq_len(Tt)) {
			rt <- resid_list[[t]]
			ut <- upper.tri(rt)
			e_t <- rt[ut]
			et_t <- t(rt)[ut]
			ok_t <- is.finite(e_t) & is.finite(et_t)
			if (sum(ok_t) > 2L) {
				rho_path[t] <- suppressWarnings(stats::cor(e_t[ok_t], et_t[ok_t]))
				if (!is.finite(rho_path[t])) rho_path[t] <- NA_real_
			}
			e_ij <- c(e_ij, e_t)
			e_ji <- c(e_ji, et_t)
		}
		ok <- is.finite(e_ij) & is.finite(e_ji)
		if (sum(ok) > 2L) {
			rho <- suppressWarnings(stats::cor(e_ij[ok], e_ji[ok]))
			if (!is.finite(rho)) rho <- 0
		}
	}
	VC <- c(va = va, cab = cab, vb = vb, rho = rho, ve = ve)

	coef_names <- rownames(coef_path)
	coef_overall <- rowMeans(coef_path)
	names(coef_overall) <- coef_names
	beta_mask <- beta_dyn$mask
	if (length(beta_mask) < length(coef_names)) {
		beta_mask <- c(beta_mask, rep(FALSE, length(coef_names) - length(beta_mask)))
	}
	beta_mask <- beta_mask[seq_along(coef_names)]
	names(beta_mask) <- coef_names
	beta_groups <- beta_dyn$groups
	if (length(beta_groups) < length(coef_names)) {
		beta_groups <- c(beta_groups,
		                 rep("", length(coef_names) - length(beta_groups)))
	}
	beta_groups <- beta_groups[seq_along(coef_names)]
	names(beta_groups) <- coef_names
	BETA <- if (isTRUE(beta_dyn$any)) {
		array(coef_path, dim = c(1L, nrow(coef_path), ncol(coef_path)),
		      dimnames = list("point", coef_names, prep$time_names))
	} else {
		matrix(coef_path[, 1L], nrow = 1L,
		       dimnames = list("point", coef_names))
	}

	out <- list(
		call = call,
		family = family, mode = mode, symmetric = symmetric, R = R,
		R_row = if (prep$bip) R_row else NULL,
		R_col = if (prep$bip) R_col else NULL,
		longitudinal = TRUE,
		mu = coef_overall["intercept"],
		beta = coef_overall[-1L],
		coef_path = coef_path,
		beta_dynamic_per_t = coef_path,
		beta_dynamic_mask = beta_mask,
		beta_dynamic_groups = beta_groups,
		BETA = BETA,
		a = a_bar, b = b_bar,
		a_dynamic = a_mat,
		b_dynamic = b_mat,
		APM = a_bar,
		BPM = b_bar,
		U = U, V = V, L = L,
		G = if (prep$bip && R_row > 0L && R_col > 0L) {
			if (isTRUE(dynamic_G) && !is.null(G_cube)) {
				apply(G_cube, c(1, 2), mean)
			} else G_static
		} else NULL,
		G_cube = if (isTRUE(dynamic_G)) G_cube else NULL,
		G_cube_post_mean = if (isTRUE(dynamic_G)) G_cube else NULL,
		G_cube_post_sd = if (isTRUE(dynamic_G) && !is.null(G_cube)) {
			array(0, dim = dim(G_cube), dimnames = dimnames(G_cube))
		} else NULL,
		Omat = Omat_store,
		Oarr = Oarr,
		UVPM = if (isTRUE(dynamic_uv) || isTRUE(dynamic_G)) Oarr else Omat_store,
		coefficients = coef_overall,
		VC = VC, s2 = ve,
		EZ = EZ, fitted = fitted_list, YPM = fitted_list,
		residuals = resid_list,
		iterations = iterations, converged = converged,
		deviance = objective_trace[length(objective_trace)],
		dev_history = objective_trace,
		objective_history = objective_trace,
		n_time = Tt,
		dims = list(n_row = nr, n_col = nc, p = length(coef_overall)),
		row_names = prep$row_names, col_names = prep$col_names,
			x_names = coef_names[-1L], time_names = prep$time_names,
			Y = prep$Yarr, X = prep$X, W_row = prep$W_row, W_col = prep$W_col,
			W_row_arr = prep$W_row_arr, W_col_arr = prep$W_col_arr,
			changing_composition = isTRUE(prep$changing_composition),
			row_presence = prep$row_presence,
			col_presence = prep$col_presence,
			actor_presence = prep$actor_presence,
			dynamic_uv = isTRUE(dynamic_uv),
		dynamic_uv_kind = if (isTRUE(dynamic_uv)) dynamic_uv_kind else NULL,
		dynamic_ab = isTRUE(dynamic_ab),
		dynamic_beta = isTRUE(beta_dyn$any),
		dynamic_G = isTRUE(dynamic_G),
		dynamic_rho = !is.null(rho_path) && length(rho_path) > 1L,
		rho_path = rho_path,
		rho_uv = if (isTRUE(dynamic_uv)) rho_uv else NULL,
		sigma_uv = NULL,
		rho_G = if (isTRUE(dynamic_G)) rho_g else NULL,
		RHO_G = if (isTRUE(dynamic_G)) rho_g else NULL,
		SIGMA_G2 = if (isTRUE(dynamic_G)) 1 / max(lambda_g, 1e-8) else NULL,
		rho_ab = if (isTRUE(dynamic_ab)) rho_ab else NULL,
		sigma_ab = NULL,
		rho_beta = if (isTRUE(beta_dyn$any)) rho_beta else NULL,
		sigma_beta = NULL,
		lambda_u = if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) transition_u else NULL,
		lambda_v = if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) transition_v else NULL,
		lambda_uv_als = if (isTRUE(dynamic_uv)) lambda_uv else NULL,
		lambda_ab_als = if (isTRUE(dynamic_ab)) lambda_ab else NULL,
		lambda_beta_als = if (isTRUE(beta_dyn$any)) lambda_beta else NULL,
		lambda_G_als = if (isTRUE(dynamic_G)) lambda_g else NULL,
		lowrank_method = lowrank_method,
		non_normal_method = if (identical(family, "normal")) "none" else "irls",
		link = link,
		obs_weights = obs_weights,
		irls_trace = irls_trace,
		stability = stability,
		dynamic = list(
			uv = isTRUE(dynamic_uv),
			uv_kind = if (isTRUE(dynamic_uv)) dynamic_uv_kind else NULL,
			ab = isTRUE(dynamic_ab), ab_kind = "ar1_penalty",
			beta = isTRUE(beta_dyn$any), beta_kind = dynamic_beta_kind,
			G = isTRUE(dynamic_G),
			snap_model_estimated = FALSE,
			t_model_estimated = isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t"),
			mode = mode),
			diagnostics = list(
				converged = converged,
				objective = objective_trace,
				convergence_trace = convergence_trace,
				objective_label = objective_label,
				objective_increases = objective_increases,
				max_objective_increase = max_objective_increase,
				last_objective_delta = if (length(objective_trace) > 1L) {
					utils::tail(diff(objective_trace), 1L)
				} else NA_real_,
			last_relative_objective_delta = if (length(objective_trace) > 1L) {
				abs(utils::tail(diff(objective_trace), 1L)) /
					max(1, abs(objective_trace[length(objective_trace) - 1L]))
			} else NA_real_,
			max_objective_delta = if (length(objective_trace) > 1L)
				max(abs(diff(objective_trace))) else NA_real_,
			last_max_fitted_delta = if (!is.null(convergence_trace) &&
			                            nrow(convergence_trace) > 0L) {
				utils::tail(convergence_trace$max_fitted_delta, 1L)
			} else NA_real_,
			last_max_parameter_delta = if (!is.null(convergence_trace) &&
			                               nrow(convergence_trace) > 0L) {
				utils::tail(convergence_trace$max_parameter_delta, 1L)
			} else NA_real_,
			convergence_tolerance = tol,
			fitted_value_tolerance = if (!is.null(tol)) sqrt(tol) else NA_real_,
			fitted_value_delta_scope = "observed_cells",
			message = if (converged) "Converged" else "Stopped at max_iter before the convergence criteria were reached."
			),
		meta = list(
			sampler = "dynamic_als",
			algorithm = "dynamic penalized least-squares block coordinate descent",
			uncertainty_available = FALSE,
			approximation_note = paste0(
				if (identical(family, "normal")) "Gaussian" else "IRLS",
				" dynamic ALS/MAP point estimate with RW/AR smoothing penalties. ",
				"It provides fitted paths and diagnostics, not posterior draws.")
		)
	)
	class(out) <- c("lame_dynamic_als", "lame_als", "ame_als")
	out
}

lame_dynamic_als <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
                             R = 0, R_row = NULL, R_col = NULL,
                             family = "normal",
                             mode = c("unipartite", "bipartite"),
                             symmetric = FALSE,
                             dynamic_ab = FALSE,
                             dynamic_uv = FALSE,
                             dynamic_uv_kind = c("ar1", "snap", "t"),
                             dynamic_beta = FALSE,
                             dynamic_beta_kind = c("ar1", "rw1", "rw2", "matern32"),
                             dynamic_G = FALSE,
                             prior = list(),
                             lambda_ab = NULL,
                             lambda_beta = NULL,
                             lambda_uv = NULL,
                             max_iter = 200,
                             tol = 5e-5,
                             link = c("probit", "logit"),
                             stability = c("none", "quick", "validation"),
                             lowrank_method = c("hybrid", "mm", "als"),
                             verbose = TRUE,
                             seed = 6886,
                             bootstrap = 0L,
                             bootstrap_type = c("parametric", "block"),
                             bootstrap_block_length = 1L,
                             bootstrap_seed = NULL,
                             call = NULL,
                             start_jitter = 0,
                             ...) {
	mode <- match.arg(mode)
	dynamic_uv_kind <- match.arg(dynamic_uv_kind)
	dynamic_beta_kind <- match.arg(dynamic_beta_kind)
	link <- match.arg(link)
	stability <- match.arg(stability)
	lowrank_method <- match.arg(lowrank_method)
	bootstrap_type <- match.arg(bootstrap_type)
	# seed locally: restore the global rng stream on exit so a downstream
	# random draw is not silently perturbed by having fit a model
	if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		.old_seed <- get(".Random.seed", envir = globalenv())
		on.exit(assign(".Random.seed", .old_seed, envir = globalenv()),
		        add = TRUE)
	}
	amen_aliases <- c(nrm = "normal", bin = "binary", pois = "poisson")
	if (length(family) == 1L && is.character(family) &&
	    family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
			"i" = "Update your script to {.code family = \"{new_family}\"} when convenient."))
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
	    !family %in% c("normal", "binary", "poisson")) {
		cli::cli_abort(c(
			"{.fn lame_dynamic_als} currently supports {.val normal}, {.val binary}, and {.val poisson}.",
			"i" = "Use {.code method = \"mcmc\"} for dynamic {.val {family}} models."))
	}
	if (identical(family, "poisson")) link <- "log"
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "snap")) {
		cli::cli_abort(c(
			"{.fn lame_dynamic_als} supports smooth {.arg dynamic_uv_kind} values {.val ar1} and {.val t}.",
			"i" = "The snap-only ALS route is available for supported normal unipartite and bipartite fits without other dynamic effects.",
			"i" = "Use {.fn lame_snap_als} or {.code lame(..., method = \"als\", dynamic_uv = TRUE, dynamic_uv_kind = \"snap\")} for that route."))
	}
	if (isTRUE(dynamic_G)) {
		if (!identical(mode, "bipartite")) {
			cli::cli_abort(c(
				"{.arg dynamic_G = TRUE} only applies to bipartite networks.",
				"i" = "For unipartite networks the multiplicative term uses {.arg dynamic_uv}."))
		}
	}
	if (identical(dynamic_beta_kind, "matern32")) {
		cli::cli_abort(c(
			"{.arg dynamic_beta_kind = \"matern32\"} is not yet supported by {.fn lame_dynamic_als}.",
			"i" = "Use {.val ar1}, {.val rw1}, or {.val rw2}, or fit with {.code method = \"mcmc\"}."))
	}
	if (!is.numeric(bootstrap) || bootstrap < 0L ||
	    !isTRUE(all.equal(bootstrap, round(bootstrap)))) {
		cli::cli_abort("{.arg bootstrap} must be a non-negative integer.")
	}
	bootstrap <- as.integer(bootstrap)
	if (!is.numeric(start_jitter) || length(start_jitter) != 1L ||
	    !is.finite(start_jitter) || start_jitter < 0) {
		cli::cli_abort("{.arg start_jitter} must be a non-negative finite number.")
	}

	if (is.null(call)) call <- match.call()
	call_args <- list(
		Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
		R = R, R_row = R_row, R_col = R_col,
		family = family, mode = mode, symmetric = symmetric,
		dynamic_ab = dynamic_ab, dynamic_uv = dynamic_uv,
		dynamic_uv_kind = dynamic_uv_kind,
		dynamic_beta = dynamic_beta,
		dynamic_beta_kind = dynamic_beta_kind,
		dynamic_G = dynamic_G, prior = prior,
		lambda_ab = lambda_ab, lambda_beta = lambda_beta,
		lambda_uv = lambda_uv, max_iter = max_iter, tol = tol,
		link = link, stability = "none",
		lowrank_method = lowrank_method,
		bootstrap_type = bootstrap_type,
		bootstrap_block_length = bootstrap_block_length,
		call = call)
	prep <- .ae_prepare(Y, Xdyad, Xrow, Xcol, mode, longitudinal = TRUE,
	                    node_time_policy = "quiet")
	if (mode == "unipartite" && prep$nr != prep$nc) {
		cli::cli_abort("Unipartite dynamic ALS requires square network slices.")
	}
	if (symmetric && prep$bip) {
		cli::cli_warn("{.arg symmetric} is ignored for bipartite networks.")
		symmetric <- FALSE
	}
	if (symmetric && prep$nr != prep$nc) {
		cli::cli_abort("Symmetric dynamic ALS requires a square network.")
	}
	if (symmetric && prep$nr == prep$nc) {
		for (t in seq_len(prep$Tt)) {
			Yt <- prep$Yarr[, , t]
			off <- Yt - t(Yt)
			if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
				cli::cli_abort(c(
					"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
					"i" = "Symmetrize the input or set {.arg symmetric} = FALSE."))
			}
		}
	}
	if (prep$Tt < 2L) {
		cli::cli_abort("{.fn lame_dynamic_als} requires at least 2 time periods.")
	}
	if (length(R) != 1L || anyNA(R) || R < 0 ||
	    !isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	R <- as.integer(R)
	.validate_rank <- function(x, nm) {
		if (length(x) != 1L || anyNA(x) || x < 0 ||
		    !isTRUE(all.equal(x, round(x)))) {
			cli::cli_abort("{.arg {nm}} must be a single non-negative integer.")
		}
		as.integer(x)
	}
	if (prep$bip) {
		R_row_eff <- .validate_rank(.lda_default(R_row, R), "R_row")
		R_col_eff <- .validate_rank(.lda_default(R_col, R), "R_col")
	} else {
		if (!is.null(R_row) || !is.null(R_col)) {
			cli::cli_warn("{.arg R_row} and {.arg R_col} are ignored for unipartite dynamic ALS.")
		}
		R_row_eff <- R
		R_col_eff <- R
	}
	R_eff <- if (prep$bip) max(R_row_eff, R_col_eff) else R
	if ((isTRUE(dynamic_uv) || isTRUE(dynamic_G)) &&
	    (R_row_eff <= 0L || R_col_eff <= 0L)) {
		cli::cli_abort(c(
			"Dynamic multiplicative ALS requires a positive multiplicative latent rank.",
			"i" = if (prep$bip) {
				"Set positive {.arg R_row} and {.arg R_col} values, or set {.arg R} > 0 for both."
			} else {
				"Set {.arg R} to at least 1, or turn off {.arg dynamic_uv} / {.arg dynamic_G}."
			}))
	}
	if (prep$bip) {
		if (R_row_eff >= prep$nr) {
			cli::cli_abort("{.arg R_row} = {R_row_eff} must be smaller than n_row = {prep$nr}.")
		}
		if (R_col_eff >= prep$nc) {
			cli::cli_abort("{.arg R_col} = {R_col_eff} must be smaller than n_col = {prep$nc}.")
		}
	} else {
		rmax <- min(prep$nr, prep$nc)
		if (R >= rmax) {
			cli::cli_abort("{.arg R} = {R} must be smaller than min(n_row, n_col) = {rmax}.")
		}
	}
	if (symmetric && prep$p_dyad > 0L) {
		for (k in seq_len(prep$p_dyad)) for (t in seq_len(prep$Tt)) {
			prep$X[, , k, t] <- .ae_symmetrize(prep$X[, , k, t])
		}
	}
	if (symmetric && lowrank_method != "mm" && R > 0L) {
		cli::cli_warn(c(
			"{.arg lowrank_method} = {.val {lowrank_method}} is not available for symmetric models.",
			"i" = "The symmetric eigenmodel uses the majorise-minimise solver."))
		lowrank_method <- "mm"
	}

	base_layout <- .lda_coef_layout(prep, include_nodes = FALSE)
	full_layout <- .lda_coef_layout(prep, include_nodes = TRUE)
	dynamic_beta_parse <- dynamic_beta
	if (is.logical(dynamic_beta) && length(dynamic_beta) > 1L &&
	    length(dynamic_beta) == length(base_layout$names) &&
	    length(full_layout$names) > length(base_layout$names)) {
		dynamic_beta_parse <- c(dynamic_beta,
		                        rep(FALSE, length(full_layout$names) -
		                             length(base_layout$names)))
	}
	beta_dyn_full <- parse_dynamic_beta(dynamic_beta_parse,
	                                    full_layout$names, full_layout$block,
	                                    intercept = TRUE, T = prep$Tt,
	                                    family = family, mode = mode)
	node_dynamic <- beta_dyn_full$mask &
		full_layout$block %in% c("row", "col")
	if (any(node_dynamic) && "intercept" %in% full_layout$names) {
		beta_dyn_full$mask[match("intercept", full_layout$names)] <- TRUE
		beta_dyn_full <- .lda_beta_dyn_from_mask(
			beta_dyn_full$mask, full_layout$names, full_layout$block)
	}
	base_idx <- which(!full_layout$block %in% c("row", "col"))
	node_idx <- which(full_layout$block %in% c("row", "col"))
	include_idx <- sort(unique(c(base_idx, node_idx[beta_dyn_full$mask[node_idx]])))
	coef_names <- full_layout$names[include_idx]
	coef_block <- full_layout$block[include_idx]
	beta_dyn <- .lda_beta_dyn_from_mask(beta_dyn_full$mask[include_idx],
	                                    coef_names, coef_block)
	dyn_node_names <- .lda_dynamic_node_names(beta_dyn)
	vary_row_static <- setdiff(.lda_node_varies_names(prep, "row"),
	                           dyn_node_names)
	vary_col_static <- setdiff(.lda_node_varies_names(prep, "col"),
	                           dyn_node_names)
	if (length(vary_row_static) > 0L) {
		cli::cli_warn(c(
			"Some row node covariates vary within actor across time slices.",
			"i" = "{.fn lame_dynamic_als} uses period-specific row values only for row coefficients selected by {.arg dynamic_beta}.",
			"i" = "Using per-actor means for static row coefficient{?s}: {.val {vary_row_static}}."))
	}
	if (length(vary_col_static) > 0L) {
		cli::cli_warn(c(
			"Some column node covariates vary within actor across time slices.",
			"i" = "{.fn lame_dynamic_als} uses period-specific column values only for column coefficients selected by {.arg dynamic_beta}.",
			"i" = "Using per-actor means for static column coefficient{?s}: {.val {vary_col_static}}."))
	}
	if (!isTRUE(dynamic_ab) && !isTRUE(beta_dyn$any) &&
	    !isTRUE(dynamic_uv) && !isTRUE(dynamic_G)) {
		return(lame_als(Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
		                R = R_eff, family = family, mode = mode,
		                symmetric = symmetric, verbose = verbose, seed = seed))
	}

	lams <- .lda_default_lambdas(prep,
		lambda_ab = .lda_default(lambda_ab, prior$lambda_ab_als),
		lambda_beta = .lda_default(lambda_beta, prior$lambda_beta_als),
		lambda_uv = .lda_default(lambda_uv, prior$lambda_uv_als))
	lambda_ab <- as.numeric(lams$ab)
	lambda_beta <- as.numeric(lams$beta)
	lambda_uv <- as.numeric(lams$uv)
	lambda_g <- as.numeric(.lda_default(prior$lambda_G_als, lambda_uv))
	rho_ab <- .lda_clip_rho(prior$rho_ab_mean, 0.95)
	rho_beta <- .lda_clip_rho(prior$rho_beta_mean, 0.95)
	rho_uv <- .lda_clip_rho(prior$rho_uv_mean, 0.95)
	rho_g <- .lda_clip_rho(prior$rho_G_mean, rho_uv)

	if (verbose) {
		cli::cli_h3("Dynamic ALS / MAP point fit")
		rank_msg <- if (prep$bip) {
			paste0("R_row = ", R_row_eff, ", R_col = ", R_col_eff)
		} else paste0("R = ", R)
		cli::cli_text("Network: {prep$nr} x {prep$nc}, {prep$Tt} time slices, {rank_msg}")
		cli::cli_text("Family: {.field {family}}{ifelse(identical(family, 'normal'), '', paste0(' / ', link, ' IWLS'))}")
		cli::cli_text("Dynamic UV: {.val {isTRUE(dynamic_uv)}} | Dynamic G: {.val {isTRUE(dynamic_G)}} | Dynamic additive effects: {.val {isTRUE(dynamic_ab)}} | Dynamic beta: {.val {isTRUE(beta_dyn$any)}}")
	}

	fit0 <- .ae_run(Y, Xdyad, Xrow, Xcol, R_eff, family, mode, symmetric,
	                max_iter = min(100L, max_iter), tol = tol,
	                verbose = FALSE, seed = seed, longitudinal = TRUE,
	                call = call, lowrank_method = lowrank_method,
	                non_normal_method = "transform", link = "probit",
	                linear_solver = "eigen", multistart = "none",
	                node_time_policy = "quiet")
	coef0 <- fit0$coefficients[coef_names]
	coef0[is.na(coef0)] <- 0
	coef0[!is.finite(coef0)] <- 0
	coef_path <- matrix(rep(coef0, prep$Tt), nrow = length(coef_names),
	                    dimnames = list(coef_names, prep$time_names))
	a_mat <- matrix(fit0$a, prep$nr, prep$Tt)
	b_mat <- matrix(fit0$b, prep$nc, prep$Tt)
	Omat <- matrix(0, prep$nr, prep$nc)
	U <- NULL
	V <- NULL
	L <- NULL
	G_static <- NULL
	G_cube <- NULL

	if (R_eff > 0L && !is.null(fit0$U) &&
	    (isTRUE(symmetric) || prep$bip || !is.null(fit0$V))) {
		if (prep$bip) {
			row_dim_names <- if (R_row_eff > 0L) {
				paste0("row_dim", seq_len(R_row_eff))
			} else character(0)
			col_dim_names <- if (R_col_eff > 0L) {
				paste0("col_dim", seq_len(R_col_eff))
			} else character(0)
			U <- .lda_coerce_factor(fit0$U, prep$nr, R_row_eff,
			                        seed = seed + 11L)
			V <- .lda_coerce_factor(fit0$V, prep$nc, R_col_eff,
			                        seed = seed + 17L)
			dimnames(U) <- list(prep$row_names, row_dim_names)
			dimnames(V) <- list(prep$col_names, col_dim_names)
			G_static <- .lda_rect_identity(R_row_eff, R_col_eff)
			dimnames(G_static) <- list(row_dim_names, col_dim_names)
			if (isTRUE(dynamic_G)) {
				G_cube <- array(rep(G_static, prep$Tt),
				                dim = c(R_row_eff, R_col_eff, prep$Tt),
				                dimnames = list(row_dim_names, col_dim_names,
				                                prep$time_names))
			}
			if (isTRUE(dynamic_uv)) {
				U <- array(rep(U, prep$Tt),
				           dim = c(prep$nr, R_row_eff, prep$Tt),
				           dimnames = list(prep$row_names, row_dim_names,
				                           prep$time_names))
				V <- array(rep(V, prep$Tt),
				           dim = c(prep$nc, R_col_eff, prep$Tt),
				           dimnames = list(prep$col_names, col_dim_names,
				                           prep$time_names))
					if (isTRUE(dynamic_G)) {
						uvg <- .lda_align_uvg_cube(U, V, G_cube)
						U <- uvg$U
						V <- uvg$V
						G_cube <- uvg$G
						Oarr <- .lda_uvg_to_o_array(U, V, G_cube,
						                            Tt = prep$Tt,
						                            time_names = prep$time_names)
					} else {
						Oarr <- .lda_uvg_to_o_array(U, V, G_static,
						                            Tt = prep$Tt,
						                            time_names = prep$time_names)
					}
				} else {
					Oarr <- .lda_uvg_to_o_array(U, V,
						if (isTRUE(dynamic_G)) G_cube else G_static,
						Tt = prep$Tt, time_names = prep$time_names)
				}
			Omat <- apply(Oarr, c(1, 2), mean)
			} else {
				uv_names <- paste0("dim", seq_len(R))
				if (isTRUE(dynamic_uv) && isTRUE(symmetric) && !is.null(fit0$L)) {
					O0 <- fit0$U %*% (fit0$L * t(fit0$U))
					sf <- .lda_symmetric_factor_slice(O0, R, prep$row_names, uv_names)
					U <- sf$U
					V <- sf$V
					L <- NULL
					Omat <- sf$O
				} else {
					U <- fit0$U
					V <- if (isTRUE(symmetric) && !is.null(fit0$U)) {
						if (!is.null(fit0$L)) {
							sweep(fit0$U, 2, sign(fit0$L), "*")
						} else fit0$U
					} else fit0$V
					L <- fit0$L
					Omat <- if (isTRUE(symmetric) && !is.null(fit0$L)) {
						fit0$U %*% (fit0$L * t(fit0$U))
					} else {
						tcrossprod(U, V)
					}
				}
				if (isTRUE(dynamic_uv)) {
					U <- array(rep(U, prep$Tt), dim = c(prep$nr, R, prep$Tt),
					           dimnames = list(prep$row_names, uv_names, prep$time_names))
					V <- array(rep(V, prep$Tt), dim = c(prep$nc, R, prep$Tt),
					           dimnames = list(prep$col_names, uv_names, prep$time_names))
					uv <- .lda_align_uv_cube(U, V)
					U <- uv$U
					V <- uv$V
					Oarr <- .lda_uv_to_o_array(U, V)
					if (isTRUE(symmetric)) {
						sf <- .lda_symmetric_o_to_uv_cube(
							Oarr, R, row_names = prep$row_names,
							dim_names = uv_names, time_names = prep$time_names)
						U <- sf$U
						V <- sf$V
						Oarr <- sf$Oarr
					}
					Omat <- apply(Oarr, c(1, 2), mean)
					L <- NULL
				} else {
				if (!is.null(U)) dimnames(U) <- list(prep$row_names, uv_names)
				if (!is.null(V)) dimnames(V) <- list(prep$col_names, uv_names)
				Oarr <- .lda_o_array(Omat, prep$Tt)
			}
		}
	} else {
		Oarr <- .lda_o_array(Omat, prep$Tt)
	}
	if (start_jitter > 0) {
		set.seed(seed)
		coef_path <- coef_path +
			matrix(stats::rnorm(length(coef_path), 0, start_jitter / 5),
			       nrow = nrow(coef_path), dimnames = dimnames(coef_path))
		a_mat <- a_mat + matrix(stats::rnorm(length(a_mat), 0, start_jitter),
		                        nrow = nrow(a_mat), dimnames = dimnames(a_mat))
		b_mat <- b_mat + matrix(stats::rnorm(length(b_mat), 0, start_jitter),
		                        nrow = nrow(b_mat), dimnames = dimnames(b_mat))
		if (R_eff > 0L && !is.null(U) && !is.null(V)) {
			U <- U + array(stats::rnorm(length(U), 0, start_jitter),
			               dim = dim(U), dimnames = dimnames(U))
			V <- V + array(stats::rnorm(length(V), 0, start_jitter),
			               dim = dim(V), dimnames = dimnames(V))
			if (isTRUE(dynamic_G) && !is.null(G_cube)) {
				G_cube <- G_cube + array(stats::rnorm(length(G_cube), 0, start_jitter / 5),
				                         dim = dim(G_cube), dimnames = dimnames(G_cube))
				if (length(dim(U)) == 3L) {
					uvg <- .lda_align_uvg_cube(U, V, G_cube)
					U <- uvg$U
					V <- uvg$V
					G_cube <- uvg$G
				}
				Oarr <- .lda_uvg_to_o_array(U, V, G_cube,
				                            Tt = prep$Tt,
				                            time_names = prep$time_names)
				Omat <- apply(Oarr, c(1, 2), mean)
			} else if (prep$bip && !is.null(G_static)) {
				Oarr <- .lda_uvg_to_o_array(U, V, G_static,
				                            Tt = prep$Tt,
				                            time_names = prep$time_names)
				Omat <- apply(Oarr, c(1, 2), mean)
				} else if (length(dim(U)) == 3L) {
					uv <- .lda_align_uv_cube(U, V)
					U <- uv$U
					V <- uv$V
					Oarr <- .lda_uv_to_o_array(U, V)
					if (isTRUE(symmetric)) {
						sf <- .lda_symmetric_o_to_uv_cube(
							Oarr, R, row_names = prep$row_names,
							dim_names = dimnames(U)[[2L]],
							time_names = prep$time_names)
						U <- sf$U
						V <- sf$V
						Oarr <- sf$Oarr
					}
					Omat <- apply(Oarr, c(1, 2), mean)
				} else {
				Omat <- tcrossprod(U, V)
				Oarr <- .lda_o_array(Omat, prep$Tt)
			}
		}
	}

	Z <- if (identical(family, "normal")) prep$Yarr else .ae_working_response(prep$Yarr, family)
	W <- .lda_unit_weights(prep$Yarr)
	if (isTRUE(symmetric)) {
		for (t in seq_len(prep$Tt)) {
			Z[, , t] <- .ae_symmetrize(Z[, , t])
			W[, , t] <- .ae_symmetrize(W[, , t])
		}
	}
	use_irls <- !identical(family, "normal")
	irls_trace <- numeric(0)
	transition_u <- NULL
	transition_v <- NULL
	transition_a <- .lda_transition_from_presence(prep$row_presence)
	transition_b <- .lda_transition_from_presence(prep$col_presence)
	fit_mask <- is.finite(Z) & is.finite(W) & W > 0
	objective_trace <- numeric(0)
	convergence_trace <- data.frame(
		iter = integer(0),
		objective = numeric(0),
		relative_objective_delta = numeric(0),
		max_fitted_delta = numeric(0),
		max_parameter_delta = numeric(0)
	)
	converged <- FALSE
	old_eta <- NULL
	# fixed marginal level prior for the multiplicative terms (see
	# .lda_level_prec): active for every dynamic factor path and for the
	# static bipartite factor blocks that share the same block solver. the
	# unipartite static low-rank path uses the shared static-ALS updater
	# and is regularized there.
	uv_level_active <- R_eff > 0L && (isTRUE(dynamic_uv) || prep$bip)
	g_level_active <- isTRUE(dynamic_G)
	level_prec_uv <- if (uv_level_active) .lda_level_prec else 0
	level_prec_g <- if (g_level_active) .lda_level_prec else 0
	for (iter in seq_len(max_iter)) {
		if (isTRUE(use_irls)) {
			eta_start <- .lda_eta(prep, coef_path, a_mat, b_mat, Oarr,
			                      symmetric = symmetric, beta_dyn = beta_dyn)
			work <- .lda_irls_work(prep$Yarr, .lda_eta_array(eta_start),
			                       family = family, link = link)
			Z <- work$Z
			W <- work$W
			if (isTRUE(symmetric)) {
				for (t in seq_len(prep$Tt)) {
					Z[, , t] <- .ae_symmetrize(Z[, , t])
					W[, , t] <- .ae_symmetrize(W[, , t])
				}
			}
		}
		transition_u <- .lda_transition_from_presence(prep$row_presence)
		transition_v <- .lda_transition_from_presence(prep$col_presence)
		if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t")) {
			transition_u <- .lda_combine_transition_weights(
				.lda_t_transition_weights(U, rho_uv), transition_u)
			transition_v <- .lda_combine_transition_weights(
				.lda_t_transition_weights(V, rho_uv), transition_v)
		}
		coef_old <- coef_path
		a_old <- a_mat
		b_old <- b_mat
		O_old <- Oarr
		G_old <- if (isTRUE(dynamic_G)) G_cube else G_static

		coef_path <- .lda_update_beta_path(
			Z, W, prep, coef_path, a_mat, b_mat, Oarr,
			beta_dyn = beta_dyn, lambda_beta = lambda_beta,
			kind = dynamic_beta_kind, rho_beta = rho_beta,
			symmetric = symmetric)
		Xb <- .lda_xbeta_array(prep, coef_path, symmetric = symmetric,
		                       beta_dyn = beta_dyn)
		if (isTRUE(dynamic_ab)) {
				ab <- if (isTRUE(symmetric)) {
					.lda_update_dynamic_ab_symmetric(
						Z, W, Xb, Oarr, a_mat,
						lambda_ab = lambda_ab, rho_ab = rho_ab,
						kind = "ar1",
						transition_a = transition_a)
				} else {
					.lda_update_dynamic_ab(
						Z, W, Xb, Oarr, a_mat, b_mat,
						lambda_ab = lambda_ab, rho_ab = rho_ab,
						kind = "ar1",
						transition_a = transition_a,
						transition_b = transition_b)
				}
			a_mat <- ab$a
			b_mat <- ab$b
		} else {
			ab <- .lda_update_static_ab(Z, W, Xb, Oarr, symmetric = symmetric)
			a_mat <- ab$a
			b_mat <- ab$b
		}
		target_mult <- array(0, dim = dim(Z))
		for (t in seq_len(prep$Tt)) {
			target_mult[, , t] <- Z[, , t] - Xb[, , t] -
				outer(a_mat[, t], b_mat[, t], "+")
		}
		if (isTRUE(dynamic_G)) {
			if (isTRUE(dynamic_uv)) {
				uv <- .lda_update_dynamic_uv_bip_g(
					Z, W, Xb, a_mat, b_mat, U, V, G_cube,
					lambda_uv = lambda_uv, rho_uv = rho_uv,
					transition_u = transition_u,
					transition_v = transition_v,
					level_prec = level_prec_uv)
				U <- uv$U
				V <- uv$V
				G_cube <- uv$G
			} else if (prep$bip && R_eff > 0L) {
				lr <- .lda_update_static_uv_bip_g(
					target_mult, is.finite(Z), W, U, V, G_cube,
					level_prec = level_prec_uv)
				U <- lr$U
				V <- lr$V
				L <- NULL
			}
			G_cube <- .lda_update_dynamic_g_path(
				target_mult, is.finite(Z), W, U, V, G_cube,
				lambda_g = lambda_g, rho_g = rho_g,
				level_prec = level_prec_g)
			if (length(dim(U)) == 3L) {
				uvg <- .lda_align_uvg_cube(U, V, G_cube)
				U <- uvg$U
				V <- uvg$V
				G_cube <- uvg$G
			}
			Oarr <- .lda_uvg_to_o_array(U, V, G_cube,
			                            Tt = prep$Tt,
			                            time_names = prep$time_names)
			Omat <- apply(Oarr, c(1, 2), mean)
			L <- NULL
		} else if (isTRUE(dynamic_uv)) {
			if (prep$bip) {
				uv <- .lda_update_dynamic_uv_bip_g(
					Z, W, Xb, a_mat, b_mat, U, V, G_static,
					lambda_uv = lambda_uv, rho_uv = rho_uv,
					transition_u = transition_u,
					transition_v = transition_v,
					align = FALSE,
					level_prec = level_prec_uv)
				U <- uv$U
				V <- uv$V
				G_static <- .lda_update_static_g_bip(
					target_mult, is.finite(Z), W, U, V, G_static)
				Oarr <- .lda_uvg_to_o_array(U, V, G_static,
				                            Tt = prep$Tt,
				                            time_names = prep$time_names)
			} else {
				uv <- .lda_update_dynamic_uv_directed(
					Z, W, Xb, a_mat, b_mat, U, V,
					lambda_uv = lambda_uv, rho_uv = rho_uv,
					symmetric = symmetric,
					transition_u = transition_u,
					transition_v = transition_v,
					level_prec = level_prec_uv)
				U <- uv$U
				V <- uv$V
				Oarr <- uv$Oarr
			}
			Omat <- apply(Oarr, c(1, 2), mean)
			L <- NULL
		} else if (R_eff > 0L) {
			if (prep$bip) {
				lr <- .lda_update_static_uv_bip_g(
					target_mult, is.finite(Z), W, U, V, G_static,
					level_prec = level_prec_uv)
				U <- lr$U
				V <- lr$V
				G_static <- .lda_update_static_g_bip(
					target_mult, is.finite(Z), W, U, V, G_static)
				Oarr <- .lda_uvg_to_o_array(U, V, G_static,
				                            Tt = prep$Tt,
				                            time_names = prep$time_names)
				Omat <- apply(Oarr, c(1, 2), mean)
				L <- NULL
			} else {
				lr <- .lda_update_static_lowrank(Z, W, Xb, a_mat, b_mat, Omat, R,
				                                 symmetric = symmetric,
				                                 lowrank_method = lowrank_method)
				Omat <- lr$Omat
				U <- lr$U
				V <- lr$V
				L <- lr$L
				Oarr <- .lda_o_array(Omat, prep$Tt)
			}
		}

			obj <- .lda_objective(Z, W, prep, coef_path, a_mat, b_mat, Oarr,
			                      beta_dyn, lambda_beta, dynamic_beta_kind,
			                      rho_beta, dynamic_ab, lambda_ab, rho_ab,
			                      transition_a = transition_a,
			                      transition_b = transition_b,
			                      symmetric = symmetric,
		                      dynamic_uv = dynamic_uv, U = U, V = V,
		                      lambda_uv = lambda_uv, rho_uv = rho_uv,
		                      transition_u = transition_u,
		                      transition_v = transition_v,
		                      dynamic_G = dynamic_G, G_cube = G_cube,
		                      lambda_g = lambda_g, rho_g = rho_g,
		                      level_prec_uv = level_prec_uv,
		                      level_prec_g = level_prec_g)
		objective_trace <- c(objective_trace, obj)
		if (isTRUE(use_irls)) irls_trace <- c(irls_trace, obj)
		eta_now <- .lda_eta(prep, coef_path, a_mat, b_mat, Oarr,
		                    symmetric = symmetric, beta_dyn = beta_dyn)
		max_eta_delta <- if (is.null(old_eta)) Inf else {
			max(vapply(seq_len(prep$Tt), function(t) {
				ok <- fit_mask[, , t]
				if (!any(ok)) return(0)
				max(abs(eta_now[[t]][ok] - old_eta[[t]][ok]), na.rm = TRUE)
			}, numeric(1)))
		}
		old_eta <- eta_now
		rel_obj <- NA_real_
		param_delta <- NA_real_
		if (iter > 1L) {
			rel_obj <- abs(objective_trace[iter - 1L] - obj) /
				max(1, abs(objective_trace[iter - 1L]))
			param_delta <- max(
				.lda_max_abs(coef_path - coef_old),
				.lda_max_abs(a_mat - a_old),
				.lda_max_abs(b_mat - b_old),
				.lda_max_abs(Oarr - O_old),
				if (isTRUE(dynamic_G) && !is.null(G_cube) && !is.null(G_old))
					.lda_max_abs(G_cube - G_old) else 0,
				if (!isTRUE(dynamic_G) && !is.null(G_static) && !is.null(G_old))
					.lda_max_abs(G_static - G_old) else 0)
		}
		convergence_trace <- rbind(
			convergence_trace,
			data.frame(
				iter = iter,
				objective = obj,
				relative_objective_delta = rel_obj,
				max_fitted_delta = max_eta_delta,
				max_parameter_delta = param_delta
			)
		)
		if (iter > 1L) {
			if (is.finite(rel_obj) && rel_obj < tol &&
			    is.finite(max_eta_delta) && max_eta_delta < sqrt(tol)) {
				converged <- TRUE
				break
			}
		}
	}

	if (verbose) {
		if (converged) {
			cli::cli_alert_success("Dynamic ALS converged in {iter} iteration{?s}.")
		} else {
			cli::cli_alert_warning("Dynamic ALS stopped at max_iter = {max_iter}; returning the last iterate.")
		}
	}
	if (!converged) {
		cli::cli_warn(c(
			"{.fn lame_dynamic_als} did not meet the objective and observed fitted-value convergence criteria.",
			"i" = "Inspect {.code fit$diagnostics$convergence_trace}; increase {.arg max_iter} if the trace is still moving."))
	}
	dynamic_objective_diff <- diff(objective_trace)
	dynamic_objective_increase_tol <- tol * (abs(utils::head(objective_trace, -1L)) + 1)
	dynamic_objective_increases <- if (length(dynamic_objective_diff) > 0L) {
		sum(dynamic_objective_diff > dynamic_objective_increase_tol, na.rm = TRUE)
	} else {
		0L
	}
	if (verbose && dynamic_objective_increases > 0L) {
		cli::cli_warn(c(
			"{.fn lame_dynamic_als} objective increased in {dynamic_objective_increases} step{?s}.",
			"i" = "For binary and Poisson fits this is the IRLS working objective, not a likelihood trace.",
			"i" = "Inspect {.code fit$diagnostics$objective_increases} and the fitted-value trace."))
	}
	if (isTRUE(dynamic_uv) && identical(dynamic_uv_kind, "t") &&
	    !is.null(U) && length(dim(U)) == 3L) {
		transition_u <- .lda_combine_transition_weights(
			.lda_t_transition_weights(U, rho_uv),
			.lda_transition_from_presence(prep$row_presence))
		transition_v <- if (!is.null(V) && length(dim(V)) == 3L) {
			.lda_combine_transition_weights(
				.lda_t_transition_weights(V, rho_uv),
				.lda_transition_from_presence(prep$col_presence))
		} else NULL
	}

	fit_out <- .lda_as_fit(prep, family, mode, symmetric, R_eff,
	                       R_row_eff, R_col_eff, call, coef_path,
	                       beta_dyn, dynamic_beta_kind, a_mat, b_mat, Omat, Oarr,
	                       U, V, L, G_static, G_cube,
	                       objective_trace, converged, iter,
	                       lambda_ab, lambda_beta, lambda_uv, lambda_g,
	                       rho_ab, rho_beta, rho_uv, rho_g,
	                       dynamic_ab, dynamic_uv, dynamic_uv_kind, dynamic_G,
	                       lowrank_method,
	                       link = if (identical(family, "normal")) NULL else link,
	                       Zwork = Z, obs_weights = W, irls_trace = irls_trace,
	                       transition_u = transition_u,
	                       transition_v = transition_v,
	                       convergence_trace = convergence_trace,
	                       tol = tol)
	if (!identical(stability, "none")) {
		fit_out$stability <- .lda_dynamic_stability(
			fit_out, starts = stability, seed = seed, call_args = call_args)
	}
	if (bootstrap > 0L) {
		fit_out$bootstrap_dynamic <- .lda_dynamic_bootstrap(
			fit_out, R_boot = bootstrap, type = bootstrap_type,
			block_length = bootstrap_block_length,
			seed = if (is.null(bootstrap_seed)) seed else bootstrap_seed,
			call_args = call_args, verbose = verbose)
		fit_out$meta$uncertainty_available <- TRUE
		fit_out$meta$uncertainty_source <- "dynamic_bootstrap"
	}
	fit_out
}

#' @export
coef.lame_dynamic_als <- function(object, ...) {
	if (isTRUE(object$dynamic_beta)) object$coef_path else object$coefficients
}

#' @export
predict.lame_dynamic_als <- function(object, newdata = NULL,
                                     type = c("response", "link"), ...) {
	type <- match.arg(type)
	EZ <- object$EZ
	if (!is.null(newdata)) {
		p_dyad <- dim(object$X)[3L]
		if (p_dyad == 0L) {
			cli::cli_abort("This fit has no dyadic covariates; {.arg newdata} cannot be applied.")
		}
		Xn <- newdata
		if (is.matrix(Xn)) Xn <- array(Xn, c(dim(Xn), 1L, 1L))
		if (length(dim(Xn)) == 3L) Xn <- array(Xn, c(dim(Xn), 1L))
		if (!identical(dim(Xn), dim(object$X))) {
			cli::cli_abort(c(
				"{.arg newdata} must match the fitted dyadic design dimensions.",
				"i" = "Expected [{paste(dim(object$X), collapse = ', ')}], got [{paste(dim(Xn), collapse = ', ')}]."))
		}
		EZ <- vector("list", object$n_time)
		dyn_node_names <- character(0)
		if (!is.null(object$beta_dynamic_mask) &&
		    !is.null(object$beta_dynamic_groups)) {
			bm <- object$beta_dynamic_mask
			bg <- object$beta_dynamic_groups
			dyn_node_names <- names(bm)[bm & bg %in% c("row", "col")]
		}
		for (t in seq_len(object$n_time)) {
			if (isTRUE(object$symmetric) && p_dyad > 0L) {
				for (k in seq_len(p_dyad)) Xn[, , k, t] <- .ae_symmetrize(Xn[, , k, t])
			}
			xb_new <- matrix(object$coef_path["intercept", t],
			                 object$dims$n_row, object$dims$n_col)
			for (k in seq_len(p_dyad)) {
				xb_new <- xb_new + object$coef_path[k + 1L, t] * Xn[, , k, t]
			}
			cn <- rownames(object$coef_path)
			if (isTRUE(object$symmetric)) {
				W_node <- .lda_bind_node_matrices(object$W_row, object$W_col)
				W_node_dyn <- .lda_bind_node_matrices(
					.lda_node_matrix(object, "row", t, dynamic = TRUE),
					.lda_node_matrix(object, "col", t, dynamic = TRUE))
				if (!is.null(W_node)) {
					for (nm in intersect(colnames(W_node), cn)) {
						W_use <- if (nm %in% dyn_node_names &&
						             !is.null(W_node_dyn) &&
						             nm %in% colnames(W_node_dyn)) {
							W_node_dyn
						} else W_node
						g <- as.numeric(object$coef_path[nm, t] * W_use[, nm])
						xb_new <- xb_new + outer(g, g, "+")
					}
				}
			} else {
				if (!is.null(object$W_row)) {
					for (nm in intersect(colnames(object$W_row), cn)) {
						W_row_t <- .lda_node_matrix(
							object, "row", t,
							dynamic = nm %in% dyn_node_names)
						g <- as.numeric(object$coef_path[nm, t] * W_row_t[, nm])
						xb_new <- xb_new + matrix(g, object$dims$n_row,
						                          object$dims$n_col)
					}
				}
				if (!is.null(object$W_col)) {
					for (nm in intersect(colnames(object$W_col), cn)) {
						W_col_t <- .lda_node_matrix(
							object, "col", t,
							dynamic = nm %in% dyn_node_names)
						g <- as.numeric(object$coef_path[nm, t] * W_col_t[, nm])
						xb_new <- xb_new + matrix(g, object$dims$n_row,
						                          object$dims$n_col,
						                          byrow = TRUE)
					}
				}
			}
			Ot <- if (!is.null(object$Oarr)) object$Oarr[, , t] else object$Omat
			ez <- xb_new +
				outer(object$a_dynamic[, t], object$b_dynamic[, t], "+") +
				Ot
			if (!identical(object$mode, "bipartite") &&
			    object$dims$n_row == object$dims$n_col) {
				diag(ez) <- NA_real_
			}
			dimnames(ez) <- list(object$row_names, object$col_names)
			EZ[[t]] <- ez
		}
		names(EZ) <- object$time_names
	}
	pred <- if (type == "link") EZ else {
		lapply(EZ, function(ez) .ae_link_inverse(ez, object$family, object$link))
	}
	pred
}

#' @export
print.lame_dynamic_als <- function(x, digits = 4, ...) {
	cli::cli_text("{.strong Dynamic AME fit via penalized ALS/MAP} (fast, MCMC-free)")
	cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}} | R = {x$R}")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slices")
	cli::cli_text("Dynamic UV: {.val {isTRUE(x$dynamic_uv)}} | Dynamic G: {.val {isTRUE(x$dynamic_G)}} | Dynamic additive effects: {.val {isTRUE(x$dynamic_ab)}} | Dynamic beta: {.val {isTRUE(x$dynamic_beta)}}")
	cli::cli_text("Converged: {.val {x$converged}} in {x$iterations} iteration{?s}")
	cli::cli_text("")
	cli::cli_text("{.strong Coefficients:}")
	if (isTRUE(x$dynamic_beta)) {
		print(round(x$coef_path, digits))
	} else {
		print(round(x$coefficients, digits))
	}
	cli::cli_text("")
	cli::cli_text(cli::col_grey("Point estimate only; no posterior draws."))
	invisible(x)
}

#' @export
summary.lame_dynamic_als <- function(object, ...) {
	out <- list(
		call = object$call,
		family = object$family,
		mode = object$mode,
		R = object$R,
		n_time = object$n_time,
		dims = object$dims,
		coefficients = if (isTRUE(object$dynamic_beta)) object$coef_path else
			cbind(Estimate = object$coefficients),
		VC = object$VC,
		dynamic_uv = object$dynamic_uv,
		dynamic_G = object$dynamic_G,
		dynamic_ab = object$dynamic_ab,
		dynamic_beta = object$dynamic_beta,
		converged = object$converged,
		iterations = object$iterations,
		objective = object$deviance,
		diagnostics = object$diagnostics)
	class(out) <- "summary.lame_dynamic_als"
	out
}

#' @export
print.summary.lame_dynamic_als <- function(x, digits = 4, ...) {
	cli::cli_h2("Dynamic ALS summary")
	cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}} | R = {x$R}")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slices")
	cli::cli_text("Dynamic UV: {.val {isTRUE(x$dynamic_uv)}} | Dynamic G: {.val {isTRUE(x$dynamic_G)}} | Dynamic additive effects: {.val {isTRUE(x$dynamic_ab)}} | Dynamic beta: {.val {isTRUE(x$dynamic_beta)}}")
	objective_label <- x$diagnostics$objective_label %||% "objective"
	cli::cli_text("Converged: {.val {x$converged}} ({x$iterations} iteration{?s}); {objective_label} = {round(x$objective, digits)}")
	if (isTRUE((x$diagnostics$objective_increases %||% 0L) > 0L)) {
		cli::cli_text(cli::col_yellow("Objective increased in {x$diagnostics$objective_increases} step{?s}; max increase = {round(x$diagnostics$max_objective_increase, digits)}."))
	}
	cli::cli_rule()
	cli::cli_text("{.strong Coefficients}")
	print(round(x$coefficients, digits))
	cli::cli_text("")
	cli::cli_text("{.strong Variance summaries}")
	print(round(x$VC, digits))
	cli::cli_text(cli::col_grey("Point estimate only; no posterior draws or credible intervals."))
	invisible(x)
}
