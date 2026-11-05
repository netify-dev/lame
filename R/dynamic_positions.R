# reported positions of a dynamic_uv fit. each period's product is identified and
# its factors only up to a gauge, so the positions factor the posterior mean
# product of every period at the fitted rank:
#   directed:  U_t V_t' with the scale split evenly between U_t and V_t;
#   symmetric: U_t diag(L) U_t' with L fixed across periods;
#   bipartite: U_t G V_t' with G diagonal and fixed across periods.
# L and G carry the scale, so the coordinates have entries of order one, as the
# sampler's coordinate paths do.
# periods are matched to the one before by moves that leave every product
# unchanged (a joint rotation of U_t and V_t; a relabelling and sign change of the
# columns when L or G is fixed). posterior sds come from the saved draws, each
# factored the same way and matched onto the reported positions

# rows of a cube slice with finite entries
.dyn_rows <- function(M) {
	which(rowSums(!is.finite(M)) == 0L)
}

# leading singular or eigen pairs of a product slice; a symmetric slice keeps
# n_pos positive and k - n_pos negative eigenvalues
.dyn_pairs <- function(P, k, kind, n_pos = k) {
	if (kind == "symmetric") {
		e <- eigen((P + t(P)) / 2, symmetric = TRUE)
		pos <- which(e$values > 0)
		neg <- which(e$values < 0)
		pick <- c(
			pos[order(-e$values[pos])][seq_len(min(n_pos, length(pos)))],
			neg[order(e$values[neg])][seq_len(min(k - n_pos, length(neg)))]
		)
		vec <- e$vectors[, pick, drop = FALSE]
		return(list(u = vec, v = vec, d = e$values[pick]))
	}
	s <- svd(P, nu = k, nv = k)
	list(u = s$u, v = s$v, d = s$d[seq_len(k)])
}

# leading pairs of the rank-limited product U G V' (symmetric: U V' with V = U L)
# from small decompositions, without forming the n x m product
.dyn_draw_pairs <- function(U, V, G, k, kind, n_pos = k) {
	qu <- qr.Q(qr(U))
	if (kind == "symmetric") {
		core <- crossprod(qu, U) %*% t(crossprod(qu, V))
		pr <- .dyn_pairs(core, k, kind, n_pos)
		pr$u <- qu %*% pr$u
		pr$v <- pr$u
		return(pr)
	}
	qv <- qr.Q(qr(V))
	core <- crossprod(qu, U) %*% G %*% t(crossprod(qv, V))
	pr <- .dyn_pairs(core, k, kind)
	pr$u <- qu %*% pr$u
	pr$v <- qv %*% pr$v
	pr
}

# assigns pairs to columns: each column takes the pair that agrees most with the
# reference column, with the sign of that agreement; a symmetric column only takes
# a pair whose eigenvalue has the column's sign
.dyn_match <- function(pr, ref_u, ref_v, rows_u, rows_v, kind, col_sign = NULL) {
	k <- length(pr$d)
	K <- ncol(ref_u)
	agree <- crossprod(pr$u[rows_u, , drop = FALSE], ref_u[rows_u, , drop = FALSE])
	if (kind != "symmetric") {
		agree <- agree + crossprod(pr$v[rows_v, , drop = FALSE], ref_v[rows_v, , drop = FALSE])
	}
	agree[!is.finite(agree)] <- 0
	allowed <- matrix(TRUE, k, K)
	if (!is.null(col_sign)) allowed <- outer(sign(pr$d), col_sign, "==")
	score <- ifelse(allowed, abs(agree), -Inf)
	pair_of <- rep(NA_integer_, K)
	flip <- rep(1, K)
	for (step in seq_len(min(k, K))) {
		if (!any(is.finite(score))) break
		best <- which(score == max(score), arr.ind = TRUE)[1L, ]
		pair_of[best[2L]] <- best[1L]
		flip[best[2L]] <- if (agree[best[1L], best[2L]] < 0) -1 else 1
		score[best[1L], ] <- -Inf
		score[, best[2L]] <- -Inf
	}
	list(pair_of = pair_of, flip = flip)
}

# columns of U and V from matched pairs, each scaled against its fixed scale;
# side splits that scale between U and V (side * 1 / side keeps the product)
.dyn_place <- function(pr, m, scale, n_u, n_v, side = 1) {
	K <- length(scale)
	U <- matrix(0, n_u, K)
	V <- matrix(0, n_v, K)
	for (j in seq_len(K)) {
		p <- m$pair_of[j]
		if (is.na(p) || scale[j] == 0 || pr$d[p] / scale[j] <= 0) next
		w <- m$flip[j] * sqrt(pr$d[p] / scale[j])
		U[, j] <- pr$u[, p] * w * side
		V[, j] <- pr$v[, p] * w / side
	}
	list(U = U, V = V)
}

# rotation of X onto Y over the given rows (orthogonal Procrustes)
.dyn_rotation <- function(X, Y) {
	M <- crossprod(X, Y)
	if (!all(is.finite(M))) return(diag(ncol(X)))
	s <- svd(M)
	s$u %*% t(s$v)
}

# positions, fixed scale and posterior sds of a dynamic_uv fit (see the top of
# this file); a fit without dynamic positions is returned unchanged
.lame_dynamic_positions <- function(fit) {
	U0 <- fit[["U"]]
	if (!isTRUE(fit$dynamic_uv) || is.null(U0) || length(dim(U0)) != 3L || dim(U0)[2L] == 0L) {
		return(fit)
	}
	bip <- identical(fit$mode, "bipartite")
	sym <- !bip && isTRUE(fit$symmetric)
	kind <- if (bip) "bipartite" else if (sym) "symmetric" else "directed"
	P_all <- if (sym) fit$ULUPM %||% fit$UVPM else fit$UVPM
	V0 <- if (sym) U0 else fit[["V"]]
	if (is.null(P_all) || length(dim(P_all)) != 3L || is.null(V0) || length(dim(V0)) != 3L) {
		return(fit)
	}
	n_u <- dim(U0)[1L]
	n_v <- dim(V0)[1L]
	R_u <- dim(U0)[2L]
	R_v <- dim(V0)[2L]
	Tn <- dim(P_all)[3L]
	k <- min(R_u, R_v, n_u, n_v)
	slice <- function(A, t) matrix(A[, , t], dim(A)[1L], dim(A)[2L])

	# leading pairs of each period's posterior mean product on its present actors
	rows_u <- cols_v <- vector("list", Tn)
	pairs <- vector("list", Tn)
	n_pos <- k
	if (sym) {
		# the signature of L: the usual number of positive leading eigenvalues
		pos_counts <- vapply(seq_len(Tn), function(t) {
			P <- slice(P_all, t)
			r <- .dyn_rows(P)
			if (length(r) == 0L) return(NA_real_)
			ev <- eigen((P[r, r] + t(P[r, r])) / 2, symmetric = TRUE, only.values = TRUE)$values
			sum(ev[order(-abs(ev))][seq_len(min(k, length(ev)))] > 0)
		}, numeric(1))
		n_pos <- if (any(is.finite(pos_counts))) as.integer(round(stats::median(pos_counts, na.rm = TRUE))) else k
	}
	for (t in seq_len(Tn)) {
		P <- slice(P_all, t)
		ru <- which(rowSums(is.finite(P)) > 0L)
		cv <- which(colSums(is.finite(P)) > 0L)
		Pt <- P[ru, cv, drop = FALSE]
		Pt[!is.finite(Pt)] <- 0
		k_t <- min(k, length(ru), length(cv))
		pr <- if (k_t > 0L) {
			.dyn_pairs(Pt, k_t, kind, n_pos)
		} else {
			list(u = matrix(0, length(ru), 0L), v = matrix(0, length(cv), 0L), d = numeric(0))
		}
		u <- matrix(NA_real_, n_u, ncol(pr$u))
		v <- matrix(NA_real_, n_v, ncol(pr$v))
		u[ru, ] <- pr$u
		v[cv, ] <- pr$v
		rows_u[[t]] <- ru
		cols_v[[t]] <- cv
		pairs[[t]] <- list(u = u, v = v, d = pr$d)
	}

	U_new <- array(NA_real_, c(n_u, R_u, Tn), dimnames = dimnames(U0))
	V_new <- array(NA_real_, c(n_v, R_v, Tn), dimnames = dimnames(V0))
	scale <- NULL
	col_sign <- NULL
	if (kind == "directed") {
		for (t in seq_len(Tn)) {
			pr <- pairs[[t]]
			w <- sqrt(pmax(pr$d, 0))
			Ut <- matrix(0, n_u, R_u)
			Vt <- matrix(0, n_v, R_v)
			Ut[, seq_along(w)] <- sweep(pr$u, 2L, w, "*")
			Vt[, seq_along(w)] <- sweep(pr$v, 2L, w, "*")
			if (t > 1L) {
				ru <- intersect(rows_u[[t]], rows_u[[t - 1L]])
				cv <- intersect(cols_v[[t]], cols_v[[t - 1L]])
				Q <- .dyn_rotation(
					rbind(Ut[ru, , drop = FALSE], Vt[cv, , drop = FALSE]),
					rbind(slice(U_new, t - 1L)[ru, , drop = FALSE], slice(V_new, t - 1L)[cv, , drop = FALSE])
				)
				Ut <- Ut %*% Q
				Vt <- Vt %*% Q
			}
			U_new[, , t] <- Ut
			V_new[, , t] <- Vt
		}
	} else {
		# match the unit pairs period by period, then fix each column's scale
		K <- k
		col_sign <- if (sym) c(rep(1, n_pos), rep(-1, K - n_pos)) else NULL
		# size of each period's product, so that L or G holds the scale per actor
		size <- vapply(seq_len(Tn), function(t) {
			if (sym) length(rows_u[[t]]) else sqrt(length(rows_u[[t]]) * length(cols_v[[t]]))
		}, numeric(1))
		side <- if (sym) 1 else (n_u / n_v)^(1 / 4)
		unit_u <- unit_v <- vector("list", Tn)
		matched_d <- matrix(NA_real_, Tn, K)
		for (t in seq_len(Tn)) {
			pr <- pairs[[t]]
			if (t == 1L) {
				ref_u <- matrix(0, n_u, K)
				ref_v <- matrix(0, n_v, K)
				ref_u[, seq_len(ncol(pr$u))] <- pr$u
				ref_v[, seq_len(ncol(pr$v))] <- pr$v
				ref_u[!is.finite(ref_u)] <- 0
				ref_v[!is.finite(ref_v)] <- 0
				ru <- rows_u[[1L]]
				cv <- cols_v[[1L]]
			} else {
				ref_u <- unit_u[[t - 1L]]
				ref_v <- unit_v[[t - 1L]]
				ru <- intersect(rows_u[[t]], rows_u[[t - 1L]])
				cv <- intersect(cols_v[[t]], cols_v[[t - 1L]])
			}
			m <- .dyn_match(pr, ref_u, ref_v, ru, cv, kind, col_sign)
			uu <- matrix(0, n_u, K)
			vv <- matrix(0, n_v, K)
			for (j in seq_len(K)) {
				p <- m$pair_of[j]
				if (is.na(p)) next
				uu[, j] <- m$flip[j] * pr$u[, p]
				vv[, j] <- m$flip[j] * pr$v[, p]
				matched_d[t, j] <- pr$d[p] / max(size[t], 1)
			}
			uu[!is.finite(uu)] <- 0
			vv[!is.finite(vv)] <- 0
			unit_u[[t]] <- uu
			unit_v[[t]] <- vv
			pairs[[t]]$match <- m
		}
		scale <- colMeans(matched_d, na.rm = TRUE)
		scale[!is.finite(scale)] <- 0
		for (t in seq_len(Tn)) {
			pl <- .dyn_place(pairs[[t]], pairs[[t]]$match, scale, n_u, n_v, side)
			Ut <- matrix(0, n_u, R_u)
			Vt <- matrix(0, n_v, R_v)
			Ut[, seq_len(K)] <- pl$U
			Vt[, seq_len(K)] <- pl$V
			U_new[, , t] <- Ut
			V_new[, , t] <- Vt
		}
	}
	# absent actors have no position in a period
	for (t in seq_len(Tn)) {
		U_new[setdiff(seq_len(n_u), rows_u[[t]]), , t] <- NA_real_
		V_new[setdiff(seq_len(n_v), cols_v[[t]]), , t] <- NA_real_
	}

	# posterior sds of the reported coordinates from the saved draws
	U_sd <- V_sd <- NULL
	Ud <- fit$U_draws
	Vd <- fit$V_draws
	Gd <- fit$G_samples
	S <- if (!is.null(Ud) && length(dim(Ud)) == 4L) dim(Ud)[4L] else 0L
	draws_ok <- S > 1L && !is.null(Vd) && identical(dim(Vd)[3:4], dim(Ud)[3:4]) &&
		dim(Ud)[3L] == Tn && (!bip || (length(dim(Gd)) == 3L && dim(Gd)[3L] == S))
	# only complete draws enter, so unfilled slots of a short run are left out
	keep <- if (draws_ok) {
		which(vapply(seq_len(S), function(s) {
			all(is.finite(Ud[, , , s])) && all(is.finite(Vd[, , , s])) && (!bip || all(is.finite(Gd[, , s])))
		}, logical(1)))
	} else {
		integer(0)
	}
	S <- length(keep)
	if (S > 1L) {
		U_sum <- U_sq <- array(0, dim(U_new))
		V_sum <- V_sq <- array(0, dim(V_new))
		for (s in keep) {
			G_s <- if (bip) matrix(Gd[, , s], dim(Gd)[1L], dim(Gd)[2L]) else diag(R_u)
			for (t in seq_len(Tn)) {
				Us <- matrix(Ud[, , t, s], n_u, R_u)
				Vs <- matrix(Vd[, , t, s], n_v, R_v)
				ref_U <- slice(U_new, t)
				ref_V <- slice(V_new, t)
				ru <- rows_u[[t]]
				cv <- cols_v[[t]]
				pr <- .dyn_draw_pairs(Us, Vs, G_s, k, kind, n_pos)
				if (kind == "directed") {
					w <- sqrt(pmax(pr$d, 0))
					Uc <- matrix(0, n_u, R_u)
					Vc <- matrix(0, n_v, R_v)
					Uc[, seq_along(w)] <- sweep(pr$u, 2L, w, "*")
					Vc[, seq_along(w)] <- sweep(pr$v, 2L, w, "*")
					Q <- .dyn_rotation(
						rbind(Uc[ru, , drop = FALSE], Vc[cv, , drop = FALSE]),
						rbind(ref_U[ru, , drop = FALSE], ref_V[cv, , drop = FALSE])
					)
					Uc <- Uc %*% Q
					Vc <- Vc %*% Q
				} else {
					m <- .dyn_match(
						pr,
						ref_U[, seq_len(k), drop = FALSE],
						ref_V[, seq_len(k), drop = FALSE],
						ru,
						cv,
						kind,
						col_sign
					)
					pl <- .dyn_place(pr, m, scale, n_u, n_v, side)
					Uc <- matrix(0, n_u, R_u)
					Vc <- matrix(0, n_v, R_v)
					Uc[, seq_len(k)] <- pl$U
					Vc[, seq_len(k)] <- pl$V
				}
				U_sum[, , t] <- U_sum[, , t] + Uc
				U_sq[, , t] <- U_sq[, , t] + Uc^2
				V_sum[, , t] <- V_sum[, , t] + Vc
				V_sq[, , t] <- V_sq[, , t] + Vc^2
			}
		}
		U_sd <- sqrt(pmax((U_sq - U_sum^2 / S) / (S - 1L), 0))
		V_sd <- sqrt(pmax((V_sq - V_sum^2 / S) / (S - 1L), 0))
		U_sd[!is.finite(U_new)] <- NA_real_
		V_sd[!is.finite(V_new)] <- NA_real_
		dimnames(U_sd) <- dimnames(U_new)
		dimnames(V_sd) <- dimnames(V_new)
	}

	fit[["U"]] <- U_new
	fit$U_sd <- U_sd
	if (sym) {
		fit$L <- scale
		fit$V_sd <- NULL
	} else {
		fit[["V"]] <- V_new
		fit$V_sd <- V_sd
	}
	if (bip) {
		G <- matrix(0, R_u, R_v, dimnames = dimnames(fit[["G"]]))
		G[cbind(seq_len(k), seq_len(k))] <- scale
		fit[["G"]] <- G
	}
	fit$uv_positions <- "posterior_mean_product"
	fit
}
