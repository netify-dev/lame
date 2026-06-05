.lame_snap_default <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.lame_snap_logit <- function(p) {
	p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
	log(p) - log1p(-p)
}

.lame_snap_xbeta <- function(beta, X) {
	d <- dim(X)
	out <- array(0, dim = c(d[1L], d[2L], d[4L]))
	if (d[3L] == 0L || length(beta) == 0L) return(out)
	for (t in seq_len(d[4L])) {
		xt <- matrix(0, d[1L], d[2L])
		for (k in seq_len(d[3L])) xt <- xt + beta[k] * X[, , k, t]
		out[, , t] <- xt
	}
	out
}

.lame_snap_uv_cube <- function(U, V, symmetric = FALSE) {
	d <- dim(U)
	out <- array(0, dim = c(d[1L], dim(V)[1L], d[3L]))
	for (t in seq_len(d[3L])) {
		ot <- U[, , t, drop = FALSE]
		vt <- V[, , t, drop = FALSE]
		dim(ot) <- d[1:2]
		dim(vt) <- dim(V)[1:2]
		out[, , t] <- ot %*% t(vt)
		if (symmetric) out[, , t] <- (out[, , t] + t(out[, , t])) / 2
	}
	out
}

.lame_snap_to_3d <- function(x) {
	if (is.matrix(x)) {
		out <- array(x, dim = c(nrow(x), ncol(x), 1L))
		dimnames(out) <- list(rownames(x), colnames(x), "dyad1")
		return(out)
	}
	if (is.array(x) && length(dim(x)) == 3L) return(x)
	cli::cli_abort("Each {.arg Xdyad} slice must be a matrix or 3D array.")
}

.lame_snap_balance_panel <- function(Y, Xdyad = NULL) {
	no_change <- list(Y = Y, Xdyad = Xdyad, actor_presence = NULL,
	                  balanced_actor_set = FALSE)
	if (!(is.list(Y) && !is.data.frame(Y))) return(no_change)
	Tt <- length(Y)
	if (Tt == 0L) return(no_change)
	Ylist <- lapply(Y, as.matrix)
	dims <- vapply(Ylist, dim, integer(2L))
	rnames <- lapply(Ylist, rownames)
	cnames <- lapply(Ylist, colnames)
	all_named <- all(vapply(seq_len(Tt), function(t) {
		!is.null(rnames[[t]]) && !is.null(cnames[[t]])
	}, logical(1L)))
	same_dims <- all(dims[1L, ] == dims[1L, 1L]) && all(dims[2L, ] == dims[2L, 1L])
	same_names <- all_named && all(vapply(seq_len(Tt), function(t) {
		identical(rnames[[t]], rnames[[1L]]) &&
			identical(cnames[[t]], cnames[[1L]])
	}, logical(1L)))
	if (same_dims && (!all_named || same_names)) return(no_change)
	if (!all_named) {
		cli::cli_abort(c(
			"Changing actor composition in {.fn lame_snap_als} requires named network slices.",
			"i" = "Give every slice row names and column names so actors can be aligned across periods."))
	}
	for (t in seq_len(Tt)) {
		if (anyDuplicated(rnames[[t]]) || anyDuplicated(cnames[[t]])) {
			cli::cli_abort("Network slice {t} has duplicate actor names.")
		}
		if (!setequal(rnames[[t]], cnames[[t]])) {
			cli::cli_abort(c(
				"{.fn lame_snap_als} currently handles changing composition only for unipartite actor sets.",
				"i" = "Slice {t} has different row and column actor names."))
		}
	}
	actors <- unique(unlist(rnames, use.names = FALSE))
	n <- length(actors)
	Ybal <- vector("list", Tt)
	names(Ybal) <- names(Y)
	actor_presence <- matrix(FALSE, n, Tt,
	                         dimnames = list(actors, .lame_snap_default(names(Y),
	                                                                    paste0("t", seq_len(Tt)))))
	for (t in seq_len(Tt)) {
		out <- matrix(NA_real_, n, n, dimnames = list(actors, actors))
		out[rnames[[t]], cnames[[t]]] <- Ylist[[t]]
		Ybal[[t]] <- out
		actor_presence[rnames[[t]], t] <- TRUE
	}
	if (is.null(Xdyad)) {
		return(list(Y = Ybal, Xdyad = NULL, actor_presence = actor_presence,
		            balanced_actor_set = TRUE))
	}
	if (!(is.list(Xdyad) && length(Xdyad) == Tt)) {
		cli::cli_abort(c(
			"Changing actor composition with dyadic covariates requires {.arg Xdyad} as a per-period list.",
			"i" = "Each entry should match that period's named {.arg Y} slice."))
	}
	X3 <- lapply(Xdyad, .lame_snap_to_3d)
	p <- dim(X3[[1L]])[3L]
	xn <- dimnames(X3[[1L]])[[3L]]
	if (is.null(xn)) xn <- paste0("dyad", seq_len(p))
	Xbal <- vector("list", Tt)
	names(Xbal) <- names(Xdyad)
	for (t in seq_len(Tt)) {
		dd <- dim(X3[[t]])
		if (dd[1L] != dims[1L, t] || dd[2L] != dims[2L, t]) {
			cli::cli_abort(c(
				"{.arg Xdyad} dimensions must match {.arg Y} before actor-set padding.",
				"i" = "Y slice {t} is {dims[1L, t]}x{dims[2L, t]}; Xdyad slice {t} is {dd[1L]}x{dd[2L]}."))
		}
		if (dd[3L] != p) {
			cli::cli_abort("All {.arg Xdyad} slices must hold the same number of covariates.")
		}
		x_rn <- dimnames(X3[[t]])[[1L]]
		x_cn <- dimnames(X3[[t]])[[2L]]
		Xt <- X3[[t]]
		if (!is.null(x_rn) || !is.null(x_cn)) {
			if (is.null(x_rn) || is.null(x_cn)) {
				cli::cli_abort(c(
					"{.arg Xdyad} slice {t} has incomplete dyad names.",
					"i" = "Provide both row and column names, or provide neither and use the same order as {.arg Y}."))
			}
			if (!setequal(x_rn, rnames[[t]]) || !setequal(x_cn, cnames[[t]])) {
				cli::cli_abort(c(
					"{.arg Xdyad} names must match {.arg Y} names within each period.",
					"i" = "Slice {t} has the right dimensions but different dyad names."))
			}
			Xt <- Xt[rnames[[t]], cnames[[t]], , drop = FALSE]
		}
		out <- array(0, dim = c(n, n, p), dimnames = list(actors, actors, xn))
		out[rnames[[t]], cnames[[t]], ] <- Xt
		Xbal[[t]] <- out
	}
	list(Y = Ybal, Xdyad = Xbal, actor_presence = actor_presence,
	     balanced_actor_set = TRUE)
}

.lame_snap_procrustes <- function(A, ref) {
	if (is.null(A) || is.null(ref) || ncol(A) != ncol(ref)) return(A)
	sv <- tryCatch(svd(crossprod(A, ref)), error = function(e) NULL)
	if (is.null(sv)) return(A)
	A %*% (sv$u %*% t(sv$v))
}

.lame_snap_factorize <- function(M, R, symmetric = FALSE) {
	n <- nrow(M)
	if (symmetric) {
		eg <- eigen((M + t(M)) / 2, symmetric = TRUE)
		idx <- order(eg$values, decreasing = TRUE)[seq_len(min(R, length(eg$values)))]
		vals <- pmax(eg$values[idx], 0)
		U <- eg$vectors[, idx, drop = FALSE] %*% diag(sqrt(vals), nrow = length(vals))
		if (ncol(U) < R) U <- cbind(U, matrix(0, n, R - ncol(U)))
		return(list(U = U[, seq_len(R), drop = FALSE],
		            V = U[, seq_len(R), drop = FALSE]))
	}
	sv <- svd(M)
	rr <- min(R, length(sv$d))
	dd <- sqrt(pmax(sv$d[seq_len(rr)], 0))
	U <- sweep(sv$u[, seq_len(rr), drop = FALSE], 2, dd, "*")
	V <- sweep(sv$v[, seq_len(rr), drop = FALSE], 2, dd, "*")
	if (rr < R) {
		U <- cbind(U, matrix(0, nrow(M), R - rr))
		V <- cbind(V, matrix(0, ncol(M), R - rr))
	}
	list(U = U[, seq_len(R), drop = FALSE], V = V[, seq_len(R), drop = FALSE])
}

.lame_snap_initial_factors <- function(prep, static_fit, R, symmetric, align, seed) {
	set.seed(seed)
	n <- prep$nr
	Tt <- prep$Tt
	U <- array(0, dim = c(n, R, Tt),
	           dimnames = list(prep$row_names, paste0("dim", seq_len(R)),
	                           prep$time_names))
	V <- array(0, dim = c(n, R, Tt),
	           dimnames = list(prep$col_names, paste0("dim", seq_len(R)),
	                           prep$time_names))
	beta <- static_fit$beta[seq_len(prep$p_dyad)]
	if (prep$p_dyad == 0L) beta <- numeric(0)
	Xb <- .lame_snap_xbeta(beta, prep$X)
	ab <- outer(static_fit$a, static_fit$b, "+")
	refU <- static_fit$U
	refV <- static_fit$V
	for (t in seq_len(Tt)) {
		M <- prep$Yarr[, , t] - static_fit$mu - Xb[, , t] - ab
		M[!is.finite(M)] <- 0
		if (symmetric) M <- (M + t(M)) / 2
		fac <- .lame_snap_factorize(M, R, symmetric)
		Ut <- fac$U
		Vt <- if (symmetric) Ut else fac$V
		if (identical(align, "global")) {
			Ut <- .lame_snap_procrustes(Ut, refU)
			Vt <- if (symmetric) Ut else .lame_snap_procrustes(Vt, refV)
		} else if (identical(align, "sequential") && t > 1L) {
			Uref <- U[, , t - 1L, drop = FALSE]
			Vref <- V[, , t - 1L, drop = FALSE]
			dim(Uref) <- c(n, R)
			dim(Vref) <- c(n, R)
			Ut <- .lame_snap_procrustes(Ut, Uref)
			Vt <- if (symmetric) Ut else .lame_snap_procrustes(Vt, Vref)
		}
		if (!all(is.finite(Ut)) || max(abs(Ut), na.rm = TRUE) == 0) {
			Ut <- refU + matrix(stats::rnorm(n * R, 0, 1e-3), n, R)
		}
		if (!all(is.finite(Vt)) || max(abs(Vt), na.rm = TRUE) == 0) {
			Vt <- if (symmetric) Ut else refV + matrix(stats::rnorm(n * R, 0, 1e-3), n, R)
		}
		U[, , t] <- Ut
		V[, , t] <- Vt
	}
	list(U = U, V = V)
}

.lame_snap_active <- function(mask, symmetric) {
	n <- dim(mask)[1L]
	Tt <- dim(mask)[3L]
	row_active <- matrix(FALSE, n, Tt)
	col_active <- matrix(FALSE, n, Tt)
	for (t in seq_len(Tt)) {
		row_active[, t] <- rowSums(mask[, , t], na.rm = TRUE) > 0
		col_active[, t] <- colSums(mask[, , t], na.rm = TRUE) > 0
		if (symmetric) {
			any_active <- row_active[, t] | col_active[, t]
			row_active[, t] <- any_active
			col_active[, t] <- any_active
		}
	}
	list(row = row_active, col = col_active)
}

.lame_snap_static_update <- function(Yarr, X, U, V, mu, beta, a, b,
                                     mask, symmetric, linear_solver = "eigen") {
	d <- dim(Yarr)
	n <- d[1L]
	Tt <- d[3L]
	p <- dim(X)[3L]
	Warr <- array(0, dim = d)
	Warr[mask] <- 1
	Yf <- Yarr
	Yf[!mask] <- 0
	O <- .lame_snap_uv_cube(U, V, symmetric)
	lin <- .ae_linsolver(linear_solver, X, mask, Warr, d[1L], d[2L], Tt, p)
	ab <- outer(a, b, "+")
	coef_mb <- lin$solve(Yf - c(ab) - O)
	mu <- coef_mb[1L]
	beta <- coef_mb[-1L]
	Xb <- .lame_snap_xbeta(beta, X)

	wsum <- apply(Warr, c(1L, 2L), sum)
	base <- array(0, dim = d)
	for (t in seq_len(Tt)) {
		base[, , t] <- Warr[, , t] * (Yf[, , t] - mu - Xb[, , t] - O[, , t])
	}
	base_rs <- apply(base, 1L, sum)
	base_cs <- apply(base, 2L, sum)
	cnt_row <- rowSums(wsum)
	cnt_col <- colSums(wsum)
	cnt_row_pos <- pmax(cnt_row, 1)
	cnt_col_pos <- pmax(cnt_col, 1)
	for (inner in seq_len(200L)) {
		a_prev <- a
		b_prev <- b
		if (symmetric) {
			a <- (base_rs - as.vector(wsum %*% a)) / cnt_row_pos
			a[cnt_row == 0] <- 0
			b <- a
		} else {
			a <- (base_rs - as.vector(wsum %*% b)) / cnt_row_pos
			a[cnt_row == 0] <- 0
			b <- (base_cs - as.vector(crossprod(wsum, a))) / cnt_col_pos
			b[cnt_col == 0] <- 0
		}
		if (max(abs(a - a_prev), na.rm = TRUE) < 1e-10 &&
		    max(abs(b - b_prev), na.rm = TRUE) < 1e-10) break
	}
	adir <- .ae_additive_solve(base_rs, base_cs, wsum, cnt_row, cnt_col, symmetric)
	a <- adir$a
	b <- adir$b
	mu <- mu + mean(a) + mean(b)
	a <- a - mean(a)
	b <- b - mean(b)
	list(mu = mu, beta = beta, a = a, b = b, Xb = .lame_snap_xbeta(beta, X))
}

.lame_snap_q_value <- function(q, i, t) {
	if (is.null(q) || is.na(q[i, t])) 0 else pmin(pmax(q[i, t], 0), 1)
}

.lame_snap_damp_q <- function(q_new, q_old, snap_damping) {
	if (is.null(q_old) || snap_damping >= 1) return(q_new)
	keep <- is.finite(q_new) & is.finite(q_old)
	q_new[keep] <- snap_damping * q_new[keep] + (1 - snap_damping) * q_old[keep]
	q_new
}

.lame_snap_weighted_quantile <- function(x, w, probs) {
	keep <- is.finite(x) & is.finite(w) & w > 0
	if (!any(keep)) return(rep(NA_real_, length(probs)))
	x <- x[keep]
	w <- w[keep]
	ord <- order(x)
	x <- x[ord]
	w <- w[ord]
	cw <- cumsum(w) / sum(w)
	vapply(probs, function(p) {
		x[which(cw >= pmin(pmax(p, 0), 1))[1L]]
	}, numeric(1L))
}

.lame_snap_resid_s2 <- function(Yarr, Xb, U, V, mu, a, b, mask, symmetric) {
	O <- .lame_snap_uv_cube(U, V, symmetric)
	ab <- outer(a, b, "+")
	sse <- 0
	n_obs <- 0
	for (t in seq_len(dim(Yarr)[3L])) {
		r <- Yarr[, , t] - mu - Xb[, , t] - ab - O[, , t]
		sse <- sse + sum(r[mask[, , t]]^2)
		n_obs <- n_obs + sum(mask[, , t])
	}
	if (n_obs <= 0L || !is.finite(sse)) return(1)
	max(sse / n_obs, 1e-10)
}

.lame_snap_temporal_terms <- function(A, cvec, Xcube, q, i, t,
                                      rho_uv, sigma_uv, snap_kappa,
                                      eligible) {
	R <- dim(Xcube)[2L]
	sigma2 <- sigma_uv * sigma_uv
	kappa2 <- snap_kappa * snap_kappa
	if (t > 1L && isTRUE(eligible[i, t])) {
		qt <- .lame_snap_q_value(q, i, t)
		A <- A + diag((1 - qt) / sigma2 + qt / kappa2, R)
		cvec <- cvec + (1 - qt) * rho_uv * Xcube[i, , t - 1L] / sigma2
	}
	if (t < dim(Xcube)[3L] && isTRUE(eligible[i, t + 1L])) {
		qf <- .lame_snap_q_value(q, i, t + 1L)
		A <- A + diag((1 - qf) * rho_uv * rho_uv / sigma2, R)
		cvec <- cvec + (1 - qf) * rho_uv * Xcube[i, , t + 1L] / sigma2
	}
	list(A = A, c = cvec)
}

.lame_snap_update_uv <- function(Yarr, Xb, U, V, mu, a, b, mask,
                                 q_u, q_v, rho_uv, sigma_uv, snap_kappa,
                                 eligible_u, eligible_v, symmetric,
                                 ridge = 1e-8, obs_var = 1) {
	n <- dim(Yarr)[1L]
	R <- dim(U)[2L]
	Tt <- dim(Yarr)[3L]
	w_obs <- 1 / max(obs_var, 1e-10)
	for (t in seq_len(Tt)) {
		for (i in seq_len(n)) {
			A <- diag(ridge, R)
			cvec <- numeric(R)
			obs <- which(mask[i, , t])
			if (length(obs) > 0L) {
				for (j in obs) {
					x <- V[j, , t]
					y <- Yarr[i, j, t] - mu - Xb[i, j, t] - a[i] - b[j]
					A <- A + w_obs * tcrossprod(x)
					cvec <- cvec + w_obs * y * x
				}
			}
			if (symmetric) {
				obs_col <- setdiff(which(mask[, i, t]), obs)
				if (length(obs_col) > 0L) {
					for (j in obs_col) {
						x <- U[j, , t]
						y <- Yarr[j, i, t] - mu - Xb[j, i, t] - a[j] - b[i]
						A <- A + w_obs * tcrossprod(x)
						cvec <- cvec + w_obs * y * x
					}
				}
			}
			tt <- .lame_snap_temporal_terms(A, cvec, U, q_u, i, t,
			                                rho_uv, sigma_uv, snap_kappa,
			                                eligible_u)
			U[i, , t] <- .ae_safe_solve(tt$A, tt$c)
		}
		if (symmetric) {
			V[, , t] <- U[, , t]
		} else {
			for (j in seq_len(n)) {
				A <- diag(ridge, R)
				cvec <- numeric(R)
				obs <- which(mask[, j, t])
				if (length(obs) > 0L) {
					for (i in obs) {
						x <- U[i, , t]
						y <- Yarr[i, j, t] - mu - Xb[i, j, t] - a[i] - b[j]
						A <- A + w_obs * tcrossprod(x)
						cvec <- cvec + w_obs * y * x
					}
				}
				tt <- .lame_snap_temporal_terms(A, cvec, V, q_v, j, t,
				                                rho_uv, sigma_uv, snap_kappa,
				                                eligible_v)
				V[j, , t] <- .ae_safe_solve(tt$A, tt$c)
			}
		}
	}
	list(U = U, V = V)
}

.lame_snap_local_data_terms <- function(Yarr, Xb, U, V, mu, a, b, mask,
                                        symmetric, side, i, t, obs_var,
                                        ridge = 1e-8) {
	R <- if (identical(side, "U")) dim(U)[2L] else dim(V)[2L]
	A <- diag(ridge, R)
	cvec <- numeric(R)
	w_obs <- 1 / max(obs_var, 1e-10)
	if (identical(side, "U")) {
		obs <- which(mask[i, , t])
		if (length(obs) > 0L) {
			for (j in obs) {
				x <- V[j, , t]
				y <- Yarr[i, j, t] - mu - Xb[i, j, t] - a[i] - b[j]
				A <- A + w_obs * tcrossprod(x)
				cvec <- cvec + w_obs * y * x
			}
		}
		if (symmetric) {
			obs_col <- setdiff(which(mask[, i, t]), obs)
			if (length(obs_col) > 0L) {
				for (j in obs_col) {
					x <- U[j, , t]
					y <- Yarr[j, i, t] - mu - Xb[j, i, t] - a[j] - b[i]
					A <- A + w_obs * tcrossprod(x)
					cvec <- cvec + w_obs * y * x
				}
			}
		}
	} else {
		obs <- which(mask[, i, t])
		if (length(obs) > 0L) {
			for (j in obs) {
				x <- U[j, , t]
				y <- Yarr[j, i, t] - mu - Xb[j, i, t] - a[j] - b[i]
				A <- A + w_obs * tcrossprod(x)
				cvec <- cvec + w_obs * y * x
			}
		}
	}
	list(A = A, c = cvec, const = 0)
}

.lame_snap_quad_min <- function(A, cvec, const = 0) {
	xhat <- .ae_safe_solve(A, cvec)
	as.numeric(const - sum(cvec * xhat))
}

.lame_snap_profile_score_one <- function(Yarr, Xb, U, V, mu, a, b, mask,
                                         Xcube, q, eligible, symmetric,
                                         side, rho_uv, sigma_uv, snap_kappa,
                                         pi_snap, snap_update, threshold,
	                                     iter, max_iter, obs_var, ridge) {
	n <- dim(Xcube)[1L]
	Tt <- dim(Xcube)[3L]
	R <- dim(Xcube)[2L]
	q_prev <- q
	q <- matrix(NA_real_, n, Tt)
	drift_d2 <- matrix(NA_real_, n, Tt)
	snap_d2 <- matrix(NA_real_, n, Tt)
	logit_score <- matrix(NA_real_, n, Tt)
	profile_drift_obj <- matrix(NA_real_, n, Tt)
	profile_snap_obj <- matrix(NA_real_, n, Tt)
	profile_delta_obj <- matrix(NA_real_, n, Tt)
	sigma2 <- sigma_uv * sigma_uv
	kappa2 <- snap_kappa * snap_kappa
	for (t in 2:Tt) {
		for (i in seq_len(n)) {
			if (!isTRUE(eligible[i, t])) next
			terms <- .lame_snap_local_data_terms(
				Yarr, Xb, U, V, mu, a, b, mask, symmetric,
				side, i, t, obs_var, ridge)
			A_base <- terms$A
			c_base <- terms$c
			const_base <- terms$const
			if (t < Tt && isTRUE(eligible[i, t + 1L])) {
				qf <- .lame_snap_q_value(q_prev, i, t + 1L)
				wf <- (1 - qf) / sigma2
				if (wf > 0) {
					x_next <- Xcube[i, , t + 1L]
					A_base <- A_base + diag(wf * rho_uv * rho_uv, R)
					c_base <- c_base + wf * rho_uv * x_next
					const_base <- const_base + wf * sum(x_next * x_next)
				}
			}
			x_prev <- Xcube[i, , t - 1L]
			A_drift <- A_base + diag(1 / sigma2, R)
			c_drift <- c_base + rho_uv * x_prev / sigma2
			const_drift <- const_base + rho_uv * rho_uv * sum(x_prev * x_prev) / sigma2
			A_snap <- A_base + diag(1 / kappa2, R)
			c_snap <- c_base
			const_snap <- const_base
			obj_drift <- .lame_snap_quad_min(A_drift, c_drift, const_drift)
			obj_snap <- .lame_snap_quad_min(A_snap, c_snap, const_snap)
			drift <- sum((Xcube[i, , t] - rho_uv * Xcube[i, , t - 1L])^2)
			snap <- sum(Xcube[i, , t]^2)
			score <- .lame_snap_logit(pi_snap) +
				R * log(pmax(sigma_uv, 1e-12) / snap_kappa) +
				0.5 * (obj_drift - obj_snap)
			score <- pmin(pmax(score, -35), 35)
			if (identical(snap_update, "hard")) {
				qq <- as.numeric(stats::plogis(score) > threshold)
			} else if (identical(snap_update, "annealed")) {
				temp <- max(1, 2 - iter / max(1, max_iter))
				qq <- stats::plogis(score / temp)
			} else {
				qq <- stats::plogis(score)
			}
			q[i, t] <- qq
			drift_d2[i, t] <- drift
			snap_d2[i, t] <- snap
			logit_score[i, t] <- score
			profile_drift_obj[i, t] <- obj_drift
			profile_snap_obj[i, t] <- obj_snap
			profile_delta_obj[i, t] <- obj_drift - obj_snap
		}
	}
	list(q = q, drift_d2 = drift_d2, snap_d2 = snap_d2,
	     logit = logit_score, profile_drift_obj = profile_drift_obj,
	     profile_snap_obj = profile_snap_obj,
	     profile_delta_obj = profile_delta_obj)
}

.lame_snap_collect_transitions <- function(U, V, q_u, q_v, eligible_u,
                                           eligible_v, symmetric,
                                           rho_ref, q_power = 2) {
	R <- dim(U)[2L]
	Tt <- dim(U)[3L]
	collect_side <- function(Xcube, q, eligible, side) {
		out <- vector("list", 0L)
		idx <- 0L
		for (t in 2:Tt) {
			for (i in seq_len(dim(Xcube)[1L])) {
				if (!isTRUE(eligible[i, t])) next
				qq <- .lame_snap_q_value(q, i, t)
				x <- Xcube[i, , t]
				xlag <- Xcube[i, , t - 1L]
				if (!all(is.finite(x)) || !all(is.finite(xlag))) next
				idx <- idx + 1L
				out[[idx]] <- list(
					side = side,
					i = i,
					t = t,
					q = qq,
					w = pmax(0, 1 - qq)^q_power,
					num = sum(x * xlag),
					den = sum(xlag * xlag),
					ref_d2 = sum((x - rho_ref * xlag)^2) / R,
					x = x,
					xlag = xlag)
			}
		}
		out
	}
	vals <- collect_side(U, q_u, eligible_u, "U")
	if (!symmetric) {
		vals <- c(vals, collect_side(V, q_v, eligible_v, "V"))
	}
	vals
}

.lame_snap_update_hyper <- function(U, V, q_u, q_v, eligible_u, eligible_v,
                                    symmetric, estimate_rho_uv,
                                    estimate_sigma_uv, rho_uv, sigma_uv,
                                    snap_pi_prior, min_sigma,
                                    hyper_update = "robust",
                                    drift_quantile = 0.05,
                                    drift_min_transitions = 50,
                                    rho_prior_mean = 0.98,
                                    rho_prior_weight = 25,
                                    q_power = 2) {
	R <- dim(U)[2L]
	Tt <- dim(U)[3L]
	num <- 0
	den <- 0
	ss <- 0
	wcount <- 0
	accum_side <- function(Xcube, q, eligible) {
		out <- c(num = 0, den = 0, ss = 0, wcount = 0, qsum = 0, m = 0)
		for (t in 2:Tt) {
			for (i in seq_len(dim(Xcube)[1L])) {
				if (!isTRUE(eligible[i, t])) next
				qq <- .lame_snap_q_value(q, i, t)
				w <- 1 - qq
				x <- Xcube[i, , t]
				xlag <- Xcube[i, , t - 1L]
				out["num"] <- out["num"] + w * sum(x * xlag)
				out["den"] <- out["den"] + w * sum(xlag * xlag)
				out["ss"] <- out["ss"] + w * sum((x - rho_uv * xlag)^2)
				out["wcount"] <- out["wcount"] + w
				out["qsum"] <- out["qsum"] + qq
				out["m"] <- out["m"] + 1
			}
		}
		out
	}
	u_stats <- accum_side(U, q_u, eligible_u)
	v_stats <- if (symmetric) rep(0, length(u_stats)) else accum_side(V, q_v, eligible_v)
	all_stats <- u_stats + v_stats
	used_robust <- FALSE
	if (identical(hyper_update, "robust")) {
		rho_ref <- if (is.finite(rho_prior_mean)) rho_prior_mean else rho_uv
		transitions <- .lame_snap_collect_transitions(
			U, V, q_u, q_v, eligible_u, eligible_v, symmetric,
			rho_ref = rho_ref, q_power = q_power)
		if (length(transitions) > 0L) {
			ref_d2 <- vapply(transitions, `[[`, numeric(1L), "ref_d2")
			w <- vapply(transitions, `[[`, numeric(1L), "w")
			n_eff <- sum(w[is.finite(w) & w > 0])
			if (n_eff >= drift_min_transitions) {
				q_eff <- max(drift_quantile,
				             min(1, drift_min_transitions / max(n_eff, 1)))
				cut <- .lame_snap_weighted_quantile(ref_d2, w, q_eff)
				keep <- is.finite(ref_d2) & is.finite(w) & w > 0 & ref_d2 <= cut[1L]
				if (!any(keep)) keep <- is.finite(ref_d2) & is.finite(w) & w > 0
				if (estimate_rho_uv && any(keep)) {
					num_k <- sum(vapply(transitions[keep], `[[`, numeric(1L), "num") *
					             w[keep])
					den_k <- sum(vapply(transitions[keep], `[[`, numeric(1L), "den") *
					             w[keep])
					if (den_k > 0) {
						rho_raw <- num_k / den_k
						n_eff_keep <- sum(w[keep])
						if (is.finite(rho_prior_mean) && rho_prior_weight > 0) {
							rho_raw <- (n_eff_keep * rho_raw +
							            rho_prior_weight * rho_prior_mean) /
								(n_eff_keep + rho_prior_weight)
						}
						rho_uv <- pmin(0.995, pmax(-0.995, rho_raw))
					}
				}
				if (estimate_sigma_uv && any(keep)) {
					innov <- vapply(transitions[keep], function(z) {
						sum((z$x - rho_uv * z$xlag)^2) / R
					}, numeric(1L))
					wk <- w[keep]
					sig2 <- stats::weighted.mean(innov, wk)
					if (is.finite(sig2) && sig2 > 0) sigma_uv <- sqrt(sig2)
				}
				used_robust <- TRUE
			}
		}
		if (used_robust) sigma_uv <- max(min_sigma, sigma_uv)
	}
	if (!used_robust) {
		if (estimate_rho_uv && all_stats["den"] > 0) {
			rho_uv <- pmin(0.99, pmax(-0.99, all_stats["num"] / all_stats["den"]))
		}
		if (estimate_sigma_uv && all_stats["wcount"] > 0) {
			ss2 <- 0
			wcount <- 0
			for (t in 2:Tt) {
				for (i in seq_len(dim(U)[1L])) {
					if (isTRUE(eligible_u[i, t])) {
						w <- 1 - .lame_snap_q_value(q_u, i, t)
						ss2 <- ss2 + w * sum((U[i, , t] - rho_uv * U[i, , t - 1L])^2)
						wcount <- wcount + w
					}
					if (!symmetric && isTRUE(eligible_v[i, t])) {
						w <- 1 - .lame_snap_q_value(q_v, i, t)
						ss2 <- ss2 + w * sum((V[i, , t] - rho_uv * V[i, , t - 1L])^2)
						wcount <- wcount + w
					}
				}
			}
			if (wcount > 0) sigma_uv <- sqrt(ss2 / max(1, wcount * R))
			sigma_uv <- max(min_sigma, sigma_uv)
		}
	}
	a_pi <- snap_pi_prior[1L]
	b_pi <- snap_pi_prior[2L]
	pi_snap <- (a_pi + all_stats["qsum"]) / (a_pi + b_pi + all_stats["m"])
	list(rho_uv = rho_uv, sigma_uv = sigma_uv, pi_snap = as.numeric(pi_snap))
}

.lame_snap_objective <- function(Yarr, Xb, U, V, mu, a, b, mask, q_u, q_v,
                                 rho_uv, sigma_uv, snap_kappa, pi_snap,
                                 eligible_u, eligible_v, symmetric,
                                 obs_var = 1) {
	O <- .lame_snap_uv_cube(U, V, symmetric)
	ab <- outer(a, b, "+")
	sse <- 0
	for (t in seq_len(dim(Yarr)[3L])) {
		r <- Yarr[, , t] - mu - Xb[, , t] - ab - O[, , t]
		sse <- sse + sum(r[mask[, , t]]^2)
	}
	sse <- sse / max(obs_var, 1e-10)
	temporal <- 0
	add_side <- function(Xcube, q, eligible) {
		val <- 0
		for (t in 2:dim(Xcube)[3L]) {
			for (i in seq_len(dim(Xcube)[1L])) {
				if (!isTRUE(eligible[i, t])) next
				qq <- .lame_snap_q_value(q, i, t)
				drift <- sum((Xcube[i, , t] - rho_uv * Xcube[i, , t - 1L])^2)
				snap <- sum(Xcube[i, , t]^2)
				val <- val + (1 - qq) * drift / (sigma_uv * sigma_uv) +
					qq * snap / (snap_kappa * snap_kappa)
				val <- val - 2 * ((1 - qq) * log1p(-pi_snap) + qq * log(pi_snap))
			}
		}
		val
	}
	temporal <- add_side(U, q_u, eligible_u)
	if (!symmetric) temporal <- temporal + add_side(V, q_v, eligible_v)
	sse + temporal
}

.lame_snap_fitted_lists <- function(Yarr, Xb, U, V, mu, a, b, mask, symmetric,
                                    row_names, col_names, time_names) {
	Tt <- dim(Yarr)[3L]
	ab <- outer(a, b, "+")
	O <- .lame_snap_uv_cube(U, V, symmetric)
	EZ <- fitted <- residuals <- vector("list", Tt)
	sse <- 0
	n_obs <- 0
	for (t in seq_len(Tt)) {
		ez <- mu + Xb[, , t] + ab + O[, , t]
		if (symmetric) ez <- (ez + t(ez)) / 2
		if (nrow(ez) == ncol(ez)) diag(ez) <- NA_real_
		r <- Yarr[, , t] - ez
		sse <- sse + sum(r[mask[, , t]]^2)
		n_obs <- n_obs + sum(mask[, , t])
		dimnames(ez) <- list(row_names, col_names)
		dimnames(r) <- list(row_names, col_names)
		EZ[[t]] <- ez
		fitted[[t]] <- ez
		residuals[[t]] <- r
	}
	names(EZ) <- names(fitted) <- names(residuals) <- time_names
	list(EZ = EZ, fitted = fitted, residuals = residuals,
	     s2 = if (n_obs > 0L) sse / n_obs else NA_real_,
	     deviance = sse)
}

.lame_snap_top_unstable <- function(delta_u, delta_v, q_u, q_v,
                                    prev_q_u = NULL, prev_q_v = NULL,
                                    symmetric, threshold, top_n = 10L) {
	make_side <- function(delta, q, side) {
		if (is.null(delta) || is.null(q)) return(NULL)
		keep <- is.finite(delta)
		if (!any(keep)) return(NULL)
		idx <- which(keep, arr.ind = TRUE)
		prev_q <- if (identical(side, "sender")) prev_q_u else prev_q_v
		prev_score <- if (!is.null(prev_q)) prev_q[keep] else NA_real_
		class_changed <- if (!is.null(prev_q)) {
			(q[keep] > threshold) != (prev_q[keep] > threshold)
		} else {
			NA
		}
		data.frame(
			side = side,
			actor = rownames(delta)[idx[, 1L]],
			time = colnames(delta)[idx[, 2L]],
			delta = delta[keep],
			prev_score = prev_score,
			snap_score = q[keep],
			class_changed = class_changed,
			stringsAsFactors = FALSE)
	}
	tab <- make_side(delta_u, q_u, "sender")
	if (!symmetric) tab <- rbind(tab, make_side(delta_v, q_v, "receiver"))
	if (is.null(tab) || nrow(tab) == 0L) {
		return(data.frame(side = character(0), actor = character(0),
		                  time = character(0), delta = numeric(0),
		                  prev_score = numeric(0), snap_score = numeric(0),
		                  class_changed = logical(0)))
	}
	tab <- tab[order(tab$delta, decreasing = TRUE), , drop = FALSE]
	tab[seq_len(min(top_n, nrow(tab))), , drop = FALSE]
}

.lame_snap_delta_summary <- function(delta_u, delta_v, q_u, q_v,
                                     prev_q_u, prev_q_v, symmetric,
                                     threshold, snap_stability_tol,
                                     snap_delta_quantile) {
	make_side <- function(delta, q, prev_q) {
		if (is.null(delta) || is.null(q) || is.null(prev_q)) {
			return(list(delta = numeric(0), changed = logical(0)))
		}
		keep_delta <- is.finite(delta)
		keep_class <- is.finite(q) & is.finite(prev_q)
		list(
			delta = delta[keep_delta],
			changed = (q[keep_class] > threshold) !=
				(prev_q[keep_class] > threshold))
	}
	u <- make_side(delta_u, q_u, prev_q_u)
	v <- if (symmetric) {
		list(delta = numeric(0), changed = logical(0))
	} else {
		make_side(delta_v, q_v, prev_q_v)
	}
	delta <- c(u$delta, v$delta)
	changed <- c(u$changed, v$changed)
	if (length(delta) == 0L) {
		return(list(max = NA_real_, mean = NA_real_, q95 = NA_real_,
		            q99 = NA_real_, convergence_delta = NA_real_,
		            n_unstable = NA_integer_, share_unstable = NA_real_,
		            n_class_changed = NA_integer_,
		            share_class_changed = NA_real_,
		            n_finite = 0L))
	}
	qs <- stats::quantile(
		delta,
		probs = unique(c(0.95, 0.99, snap_delta_quantile)),
		na.rm = TRUE, names = TRUE, type = 7)
	q_lookup <- function(p) {
		nm <- as.character(p)
		if (nm %in% names(qs)) as.numeric(qs[[nm]]) else
			as.numeric(stats::quantile(delta, p, na.rm = TRUE, type = 7))
	}
	list(
		max = max(delta, na.rm = TRUE),
		mean = mean(delta, na.rm = TRUE),
		q95 = q_lookup(0.95),
		q99 = q_lookup(0.99),
		convergence_delta = q_lookup(snap_delta_quantile),
		n_unstable = sum(delta > snap_stability_tol, na.rm = TRUE),
		share_unstable = mean(delta > snap_stability_tol, na.rm = TRUE),
		n_class_changed = sum(changed, na.rm = TRUE),
		share_class_changed = if (length(changed) > 0L) {
			mean(changed, na.rm = TRUE)
		} else {
			NA_real_
		},
		n_finite = length(delta))
}

#' Fast approximate dynamic snap-shift AME estimator
#'
#' @description
#' Fits a fast point-estimator approximation to the dynamic snap-shift
#' latent-factor model for longitudinal normal-valued networks.
#' It targets the same drift-versus-reset transition estimand as
#' \code{\link{lame}} with \code{dynamic_uv = TRUE} and
#' \code{dynamic_uv_kind = "snap"}, but returns ALS snap scores rather than
#' MCMC draws.
#'
#' @details
#' On longer panels, the default \code{align = "sequential"} keeps the
#' period-to-period movement in view. \code{align = "global"} can help on very
#' short panels, but on longer panels it can pull each period back toward the
#' pooled static fit and make early jumps look too strong. If you have a fixed
#' drift scale, pass \code{rho_uv} and \code{sigma_uv} with
#' \code{estimate_rho_uv = FALSE} and \code{estimate_sigma_uv = FALSE}.
#'
#' The convergence output separates the usual score summary from the worst
#' moving actor-period. By default, convergence uses the 95th percentile of
#' score changes and the share of class changes. \code{final_max_snap_delta}
#' and \code{unstable_transitions} still show the largest local moves. Read
#' \code{snap_prob} as an ALS snap score. It is useful for rankings and
#' heuristic classifications; it is not a Bayesian posterior probability.
#'
#' In bipartite mode, the fitted multiplicative term is
#' \eqn{U_t G V_t'} with one static interaction matrix \eqn{G}. The returned
#' \code{snap_prob} matrix contains row-actor snap scores and
#' \code{snap_prob_v} contains column-actor snap scores. Bipartite snap ALS
#' does not estimate node covariates, dynamic coefficients, or a dynamic
#' \eqn{G_t}; those combinations need a separate model.
#'
#' @param Y a list of relational matrices, or a 3D array. Unipartite fits use
#'   \code{[n, n, T]} panels. Bipartite fits use rectangular
#'   \code{[n_row, n_col, T]} panels. For named lists, slices may have changing
#'   actor composition; the estimator pads to the union actor set and treats
#'   absent actor-periods as unobserved and snap-ineligible.
#' @param Xdyad optional list of dyadic covariate matrices/arrays, or
#'   \code{NULL}.
#' @param Xrow,Xcol not supported by this fast snap-shift estimator; pass
#'   \code{NULL}. Use \code{\link{lame}} with \code{method = "mcmc"} for
#'   node-covariate dynamic snap models.
#' @param R positive integer latent rank. For bipartite fits, used as the
#'   default for \code{R_row} and \code{R_col}.
#' @param R_row,R_col positive integer latent ranks for bipartite row and
#'   column positions. Defaults to \code{R}.
#' @param family currently only \code{"normal"}.
#' @param mode \code{"unipartite"} or \code{"bipartite"}.
#' @param symmetric logical; if \code{TRUE}, fit a symmetric latent-factor
#'   approximation and report one snap-probability matrix.
#' @param max_iter maximum block-coordinate iterations.
#' @param tol convergence tolerance on the approximate penalized objective.
#' @param snap_kappa diffuse snap-prior standard deviation.
#' @param snap_pi_prior length-two vector \code{c(a, b)} for the beta prior
#'   used to regularize the snap rate.
#' @param snap_update one of \code{"soft"}, \code{"hard"}, or
#'   \code{"annealed"}.
#' @param snap_damping scalar in \code{(0, 1]}; damping applied to soft
#'   snap-score updates. Values below 1 mix the new profiled score with the
#'   previous score to reduce oscillation on large panels. Ignored for
#'   \code{snap_update = "hard"}.
#' @param estimate_rho_uv,estimate_sigma_uv logical flags for updating the
#'   AR drift persistence and innovation scale.
#' @param rho_uv,sigma_uv optional fixed/initial AR drift parameters.
#' @param hyper_update one of \code{"robust"} or \code{"em"}. The default
#'   \code{"robust"} estimates drift hyperparameters from the lower tail of
#'   transition innovations, which prevents broad ruptures from being absorbed
#'   into an overly diffuse drift process. \code{"em"} uses the untrimmed
#'   soft-classification moment update.
#' @param drift_quantile lower-tail transition quantile used by
#'   \code{hyper_update = "robust"} for estimating smooth-drift
#'   hyperparameters.
#' @param drift_min_transitions minimum effective number of transitions retained
#'   by the lower-tail update. This keeps the default from overfitting
#'   the smooth-drift scale on small panels.
#' @param rho_prior_mean,rho_prior_weight weak regularization for the
#'   persistence estimate. The defaults encode the snap-shift model's intended
#'   persistent, low-innovation drift baseline.
#' @param align initialization alignment mode; \code{"global"} aligns each
#'   period to the pooled static ALS fit, \code{"sequential"} aligns each
#'   period to the previous initialized period, and \code{"none"} leaves
#'   per-period factors unaligned.
#' @param threshold hard-classification threshold for \code{snap_class}.
#' @param min_sigma lower bound for the drift scale.
#' @param sigma_floor_fraction for robust hyperparameter updates, the iterative
#'   drift scale cannot fall below this fraction of the initial
#'   \code{sigma_uv}. This prevents broad snap assignments from collapsing the
#'   smooth-drift variance to \code{min_sigma}.
#' @param ridge small ridge added to latent-position normal equations.
#' @param snap_stability_tol tolerance for declaring soft snap scores stable.
#'   This is separate from \code{tol}, which tracks the penalized objective.
#' @param snap_convergence convergence criterion for snap scores.
#'   \code{"quantile"} uses \code{snap_delta_quantile} of the absolute
#'   score changes plus the class-change share; \code{"max"} uses the
#'   worst actor-period score change; \code{"classification"} uses only
#'   class-change stability. The worst-case max is always stored
#'   and warned about when high.
#' @param snap_delta_quantile quantile of actor-period score changes used by
#'   \code{snap_convergence = "quantile"}.
#' @param snap_class_change_tol maximum share of eligible actor-periods whose
#'   \code{snap_class} may change between iterations while still declaring
#'   snap convergence.
#' @param unstable_top_n number of largest actor-period snap-score changes to
#'   store in convergence diagnostics.
#' @param stability optional start-sensitivity preset. \code{"none"} runs one
#'   fit. \code{"quick"} and \code{"validation"} rerun the same estimator from
#'   additional seeds and attach \code{fit$stability} with snap-score,
#'   classification, top-period, and fitted-surface comparisons.
#' @param verbose logical; print progress.
#' @param seed integer seed for initialization perturbations.
#'
#' @return An object of class \code{"lame_snap_als"} with dynamic latent
#'   positions, snap scores, fitted values, residuals, and
#'   convergence diagnostics.
#'
#' @examples
#' set.seed(1)
#' Y_bip <- lapply(seq_len(3), function(t) {
#'   m <- matrix(rnorm(6 * 5), 6, 5)
#'   rownames(m) <- paste0("r", seq_len(6))
#'   colnames(m) <- paste0("c", seq_len(5))
#'   m
#' })
#' fit_bip <- lame_snap_als(Y_bip, R = 1, mode = "bipartite",
#'                          max_iter = 3, verbose = FALSE)
#' dim(fit_bip$snap_prob_v)
#'
#' @export
lame_snap_als <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
                          R = 2L, R_row = NULL, R_col = NULL,
                          family = "normal",
                          mode = c("unipartite", "bipartite"),
                          symmetric = FALSE,
                          max_iter = 200L, tol = 1e-6,
                          snap_kappa = 2,
                          snap_pi_prior = c(a = 1, b = 9),
                          snap_update = c("soft", "hard", "annealed"),
                          snap_damping = 0.7,
                          estimate_rho_uv = TRUE,
                          estimate_sigma_uv = TRUE,
                          rho_uv = NULL, sigma_uv = NULL,
                          hyper_update = c("robust", "em"),
                          drift_quantile = 0.05,
                          drift_min_transitions = 50,
                          rho_prior_mean = 0.98,
                          rho_prior_weight = 25,
                          align = c("sequential", "global", "none"),
                          threshold = 0.5,
                          min_sigma = 1e-4,
                          sigma_floor_fraction = 0.75,
                          ridge = 1e-8,
                          snap_stability_tol = 0.05,
                          snap_convergence = c("quantile", "max",
                                               "classification"),
                          snap_delta_quantile = 0.95,
                          snap_class_change_tol = 0.005,
                          unstable_top_n = 10L,
                          stability = c("none", "quick", "validation"),
                          verbose = TRUE,
                          seed = 6886) {
	mode <- match.arg(mode)
	align <- match.arg(align)
	snap_update <- match.arg(snap_update)
	hyper_update <- match.arg(hyper_update)
	snap_convergence <- match.arg(snap_convergence)
	stability <- match.arg(stability)
	if (identical(family, "nrm")) {
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val nrm} (amen-style) accepted as alias for {.val normal}.",
			"i" = "Update your script to {.code family = \"normal\"} when convenient."))
		family <- "normal"
	}
	if (length(family) != 1L || !identical(family, "normal")) {
		cli::cli_abort(c(
			"{.fn lame_snap_als} supports only {.code family = \"normal\"}.",
			"i" = "Binary/poisson snap ALS require an IRLS outer loop and are not supported by this estimator."))
	}
	if (!is.null(Xrow) || !is.null(Xcol)) {
		cli::cli_abort(c(
			"{.fn lame_snap_als} does not support node covariates.",
			"i" = "Pass dyadic covariates via {.arg Xdyad}; use {.fn lame} MCMC for node-covariate dynamic snap models."))
	}
	if (is.null(R) || length(R) != 1L || anyNA(R) ||
	    !isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be an integer scalar.")
	}
	if (identical(mode, "unipartite") && R < 1L) {
		cli::cli_abort("{.arg R} must be a positive integer.")
	}
	if (identical(mode, "bipartite") && R < 0L) {
		cli::cli_abort("{.arg R} must be non-negative for bipartite fits.")
	}
	R <- as.integer(R)
	if (length(max_iter) != 1L || max_iter < 1L) {
		cli::cli_abort("{.arg max_iter} must be a positive integer.")
	}
	if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
		cli::cli_abort("{.arg tol} must be a positive finite scalar.")
	}
	if (!is.numeric(snap_pi_prior) || length(snap_pi_prior) != 2L ||
	    any(!is.finite(snap_pi_prior)) || any(snap_pi_prior <= 0)) {
		cli::cli_abort("{.arg snap_pi_prior} must be a positive length-two numeric vector.")
	}
	if (!is.numeric(snap_kappa) || length(snap_kappa) != 1L ||
	    !is.finite(snap_kappa) || snap_kappa <= 0) {
		cli::cli_abort("{.arg snap_kappa} must be a positive finite scalar.")
	}
	if (!is.numeric(snap_damping) || length(snap_damping) != 1L ||
	    !is.finite(snap_damping) || snap_damping <= 0 || snap_damping > 1) {
		cli::cli_abort("{.arg snap_damping} must be a finite scalar in (0, 1].")
	}
	if (identical(snap_update, "hard")) snap_damping <- 1
	if (!is.numeric(drift_quantile) || length(drift_quantile) != 1L ||
	    !is.finite(drift_quantile) || drift_quantile <= 0 ||
	    drift_quantile >= 1) {
		cli::cli_abort("{.arg drift_quantile} must be a finite scalar in (0, 1).")
	}
	if (!is.numeric(drift_min_transitions) ||
	    length(drift_min_transitions) != 1L ||
	    !is.finite(drift_min_transitions) || drift_min_transitions < 1) {
		cli::cli_abort("{.arg drift_min_transitions} must be a finite scalar >= 1.")
	}
	if (!is.numeric(rho_prior_mean) || length(rho_prior_mean) != 1L ||
	    !is.finite(rho_prior_mean) || abs(rho_prior_mean) >= 1) {
		cli::cli_abort("{.arg rho_prior_mean} must be a finite scalar with absolute value below 1.")
	}
	if (!is.numeric(rho_prior_weight) || length(rho_prior_weight) != 1L ||
	    !is.finite(rho_prior_weight) || rho_prior_weight < 0) {
		cli::cli_abort("{.arg rho_prior_weight} must be a non-negative finite scalar.")
	}
	if (!is.numeric(snap_stability_tol) || length(snap_stability_tol) != 1L ||
	    !is.finite(snap_stability_tol) || snap_stability_tol <= 0 ||
	    snap_stability_tol > 1) {
		cli::cli_abort("{.arg snap_stability_tol} must be a finite scalar in (0, 1].")
	}
	if (!is.numeric(snap_delta_quantile) ||
	    length(snap_delta_quantile) != 1L ||
	    !is.finite(snap_delta_quantile) || snap_delta_quantile <= 0 ||
	    snap_delta_quantile > 1) {
		cli::cli_abort("{.arg snap_delta_quantile} must be a finite scalar in (0, 1].")
	}
	if (!is.numeric(snap_class_change_tol) ||
	    length(snap_class_change_tol) != 1L ||
	    !is.finite(snap_class_change_tol) || snap_class_change_tol < 0 ||
	    snap_class_change_tol > 1) {
		cli::cli_abort("{.arg snap_class_change_tol} must be a finite scalar in [0, 1].")
	}
	if (!is.numeric(sigma_floor_fraction) ||
	    length(sigma_floor_fraction) != 1L ||
	    !is.finite(sigma_floor_fraction) || sigma_floor_fraction < 0 ||
	    sigma_floor_fraction > 1) {
		cli::cli_abort("{.arg sigma_floor_fraction} must be a finite scalar in [0, 1].")
	}
	if (length(unstable_top_n) != 1L || anyNA(unstable_top_n) ||
	    unstable_top_n < 0 || !isTRUE(all.equal(unstable_top_n, round(unstable_top_n)))) {
		cli::cli_abort("{.arg unstable_top_n} must be a non-negative integer.")
	}
	unstable_top_n <- as.integer(unstable_top_n)
	stability_args <- list(
		Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
		R = R, R_row = R_row, R_col = R_col,
		family = family, mode = mode, symmetric = symmetric,
		max_iter = max_iter, tol = tol,
		snap_kappa = snap_kappa,
		snap_pi_prior = snap_pi_prior,
		snap_update = snap_update,
		snap_damping = snap_damping,
		estimate_rho_uv = estimate_rho_uv,
		estimate_sigma_uv = estimate_sigma_uv,
		rho_uv = rho_uv, sigma_uv = sigma_uv,
		hyper_update = hyper_update,
		drift_quantile = drift_quantile,
		drift_min_transitions = drift_min_transitions,
		rho_prior_mean = rho_prior_mean,
		rho_prior_weight = rho_prior_weight,
		align = align,
		threshold = threshold,
		min_sigma = min_sigma,
		sigma_floor_fraction = sigma_floor_fraction,
		ridge = ridge,
		snap_stability_tol = snap_stability_tol,
		snap_convergence = snap_convergence,
		snap_delta_quantile = snap_delta_quantile,
		snap_class_change_tol = snap_class_change_tol,
		unstable_top_n = unstable_top_n,
		stability = "none",
		verbose = FALSE)

	if (identical(mode, "bipartite")) {
		out <- .lame_snap_bipartite_als(
			Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
			R = R, R_row = R_row, R_col = R_col,
			family = family,
			max_iter = max_iter, tol = tol,
			snap_kappa = snap_kappa,
			snap_pi_prior = snap_pi_prior,
			snap_update = snap_update,
			snap_damping = snap_damping,
			estimate_rho_uv = estimate_rho_uv,
			estimate_sigma_uv = estimate_sigma_uv,
			rho_uv = rho_uv, sigma_uv = sigma_uv,
			align = align,
			threshold = threshold,
			min_sigma = min_sigma,
		ridge = ridge,
		snap_stability_tol = snap_stability_tol,
		snap_convergence = snap_convergence,
		snap_delta_quantile = snap_delta_quantile,
		snap_class_change_tol = snap_class_change_tol,
		hyper_update = hyper_update,
		drift_quantile = drift_quantile,
		drift_min_transitions = drift_min_transitions,
		sigma_floor_fraction = sigma_floor_fraction,
		unstable_top_n = unstable_top_n,
		rho_prior_mean = rho_prior_mean,
		rho_prior_weight = rho_prior_weight,
			verbose = verbose,
			seed = seed)
		out$call <- match.call()
		return(.lame_snap_attach_stability(out, stability, stability_args,
		                                   seed = seed, verbose = verbose))
	}

	panel <- .lame_snap_balance_panel(Y, Xdyad)
	Y_fit <- panel$Y
	Xdyad_fit <- panel$Xdyad
	prep <- .ae_prepare(Y_fit, Xdyad_fit, NULL, NULL, mode, longitudinal = TRUE)
	if (prep$nr != prep$nc) {
		cli::cli_abort("{.fn lame_snap_als} requires square unipartite network slices.")
	}
	if (prep$Tt < 2L) {
		cli::cli_abort("{.fn lame_snap_als} requires at least two time periods.")
	}
	if (R >= prep$nr) {
		cli::cli_abort("{.arg R} must be smaller than the number of actors.")
	}
	if (symmetric) {
		for (t in seq_len(prep$Tt)) {
			Yt <- prep$Yarr[, , t]
			d <- Yt - t(Yt)
			if (any(is.finite(d)) && max(abs(d), na.rm = TRUE) > 1e-8) {
				cli::cli_abort("{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.")
			}
			if (prep$p_dyad > 0L) for (k in seq_len(prep$p_dyad)) {
				prep$X[, , k, t] <- .ae_symmetrize(prep$X[, , k, t])
			}
		}
	}
	mask <- is.finite(prep$Yarr)
	active <- .lame_snap_active(mask, symmetric)
	eligible_u <- matrix(FALSE, prep$nr, prep$Tt)
	eligible_v <- matrix(FALSE, prep$nr, prep$Tt)
	for (t in 2:prep$Tt) {
		eligible_u[, t] <- active$row[, t] & active$row[, t - 1L]
		eligible_v[, t] <- active$col[, t] & active$col[, t - 1L]
	}
	if (symmetric) eligible_v <- eligible_u

	static_fit <- suppressWarnings(lame_als(
		Y_fit, Xdyad = Xdyad_fit, R = R, family = "normal",
		mode = "unipartite", symmetric = symmetric,
		max_iter = max(50L, min(100L, max_iter)),
		verbose = FALSE, seed = seed))
	init <- .lame_snap_initial_factors(prep, static_fit, R, symmetric, align, seed)
	U <- init$U
	V <- if (symmetric) init$U else init$V
	mu <- static_fit$mu
	beta <- static_fit$beta[seq_len(prep$p_dyad)]
	if (prep$p_dyad == 0L) beta <- numeric(0)
	a <- static_fit$a
	b <- static_fit$b
	Xb <- .lame_snap_xbeta(beta, prep$X)
	q_u <- matrix(NA_real_, prep$nr, prep$Tt,
	              dimnames = list(prep$row_names, prep$time_names))
	q_v <- matrix(NA_real_, prep$nr, prep$Tt,
	              dimnames = list(prep$col_names, prep$time_names))
	q_u[eligible_u] <- 0
	q_v[eligible_v] <- 0
	q_u[, 1L] <- NA_real_
	q_v[, 1L] <- NA_real_

	if (is.null(rho_uv)) {
		rho_uv <- .lame_snap_update_hyper(
			U, V, q_u, q_v, eligible_u, eligible_v, symmetric,
			estimate_rho_uv = TRUE, estimate_sigma_uv = FALSE,
			rho_uv = 0.8, sigma_uv = 1, snap_pi_prior = snap_pi_prior,
			min_sigma = min_sigma, hyper_update = hyper_update,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions,
			rho_prior_mean = rho_prior_mean,
			rho_prior_weight = rho_prior_weight)$rho_uv
	}
	if (is.null(sigma_uv)) {
		sigma_uv <- .lame_snap_update_hyper(
			U, V, q_u, q_v, eligible_u, eligible_v, symmetric,
			estimate_rho_uv = FALSE, estimate_sigma_uv = TRUE,
			rho_uv = rho_uv, sigma_uv = 1, snap_pi_prior = snap_pi_prior,
			min_sigma = min_sigma, hyper_update = hyper_update,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions,
			rho_prior_mean = rho_prior_mean,
			rho_prior_weight = rho_prior_weight)$sigma_uv
	}
	rho_uv <- as.numeric(rho_uv)
	sigma_uv <- max(min_sigma, as.numeric(sigma_uv))
	n_hyper_transitions <- sum(eligible_u, na.rm = TRUE)
	if (!symmetric) n_hyper_transitions <- n_hyper_transitions + sum(eligible_v, na.rm = TRUE)
	robust_hyper_active <- identical(hyper_update, "robust") &&
		n_hyper_transitions >= drift_min_transitions
	sigma_update_floor <- if (robust_hyper_active && estimate_sigma_uv) {
		max(min_sigma, sigma_floor_fraction * sigma_uv)
	} else {
		min_sigma
	}
	obs_var <- if (!is.null(static_fit$s2) && is.finite(static_fit$s2)) {
		max(static_fit$s2, 1e-10)
	} else {
		.lame_snap_resid_s2(prep$Yarr, Xb, U, V, mu, a, b, mask, symmetric)
	}
	pi_snap <- snap_pi_prior[1L] / sum(snap_pi_prior)
	objective_trace <- numeric(0)
	param_trace <- data.frame(iter = integer(0), objective = numeric(0),
	                          rho_uv = numeric(0), sigma_uv = numeric(0),
	                          pi_snap = numeric(0), obs_var = numeric(0),
	                          max_snap_delta = numeric(0),
	                          snap_delta_q95 = numeric(0),
	                          snap_delta_q99 = numeric(0),
	                          snap_delta_conv = numeric(0),
	                          mean_snap_delta = numeric(0),
	                          n_snap_delta_gt_tol = integer(0),
	                          share_snap_delta_gt_tol = numeric(0),
	                          n_snap_class_changed = integer(0),
	                          share_snap_class_changed = numeric(0))
	converged <- FALSE
	prev_obj <- Inf
	last_q_u <- q_u
	last_q_v <- q_v
	prev_delta_q_u <- q_u
	prev_delta_q_v <- q_v
	last_delta_u <- matrix(NA_real_, prep$nr, prep$Tt,
	                       dimnames = list(prep$row_names, prep$time_names))
	last_delta_v <- matrix(NA_real_, prep$nr, prep$Tt,
	                       dimnames = list(prep$col_names, prep$time_names))
	last_delta_summary <- NULL

	if (verbose) {
		cli::cli_h3("Dynamic snap ALS")
		cli::cli_text("Network: {prep$nr} actors, {prep$Tt} periods, R = {R}")
		cli::cli_text("Update: {.val {snap_update}} | hyper: {.val {hyper_update}} | align: {.val {align}} | symmetric: {.val {symmetric}}")
	}
	for (iter in seq_len(max_iter)) {
		stat <- .lame_snap_static_update(prep$Yarr, prep$X, U, V, mu, beta, a, b,
		                                 mask, symmetric)
		mu <- stat$mu
		beta <- stat$beta
		a <- stat$a
		b <- stat$b
		Xb <- stat$Xb
		uv <- .lame_snap_update_uv(prep$Yarr, Xb, U, V, mu, a, b, mask,
		                            q_u, q_v, rho_uv, sigma_uv, snap_kappa,
		                            eligible_u, eligible_v, symmetric, ridge,
		                            obs_var)
		U <- uv$U
		V <- if (symmetric) uv$U else uv$V
		obs_var <- .lame_snap_resid_s2(prep$Yarr, Xb, U, V, mu, a, b, mask,
		                               symmetric)
		qu <- .lame_snap_profile_score_one(
			prep$Yarr, Xb, U, V, mu, a, b, mask, U, q_u, eligible_u,
			symmetric, "U", rho_uv, sigma_uv, snap_kappa, pi_snap,
			snap_update, threshold, iter, max_iter, obs_var, ridge)
		qu$q <- .lame_snap_damp_q(qu$q, q_u, snap_damping)
		q_u <- qu$q
		qv <- if (symmetric) {
			list(q = q_u, drift_d2 = qu$drift_d2, snap_d2 = qu$snap_d2,
			     logit = qu$logit,
			     profile_drift_obj = qu$profile_drift_obj,
			     profile_snap_obj = qu$profile_snap_obj,
			     profile_delta_obj = qu$profile_delta_obj)
		} else {
			.lame_snap_profile_score_one(
				prep$Yarr, Xb, U, V, mu, a, b, mask, V, q_v, eligible_v,
				symmetric, "V", rho_uv, sigma_uv, snap_kappa, pi_snap,
				snap_update, threshold, iter, max_iter, obs_var, ridge)
		}
		if (!symmetric) qv$q <- .lame_snap_damp_q(qv$q, q_v, snap_damping)
		q_v <- qv$q
		dimnames(q_u) <- list(prep$row_names, prep$time_names)
		dimnames(q_v) <- list(prep$col_names, prep$time_names)
		hyp <- .lame_snap_update_hyper(U, V, q_u, q_v, eligible_u, eligible_v,
		                               symmetric, estimate_rho_uv,
		                               estimate_sigma_uv, rho_uv, sigma_uv,
		                               snap_pi_prior, min_sigma,
		                               hyper_update = hyper_update,
		                               drift_quantile = drift_quantile,
		                               drift_min_transitions = drift_min_transitions,
		                               rho_prior_mean = rho_prior_mean,
		                               rho_prior_weight = rho_prior_weight)
		rho_uv <- hyp$rho_uv
		sigma_uv <- max(sigma_update_floor, hyp$sigma_uv)
		pi_snap <- hyp$pi_snap
		obj <- .lame_snap_objective(prep$Yarr, Xb, U, V, mu, a, b, mask,
		                             q_u, q_v, rho_uv, sigma_uv, snap_kappa,
		                             pi_snap, eligible_u, eligible_v, symmetric,
		                             obs_var)
		objective_trace <- c(objective_trace, obj)
		prev_delta_q_u <- last_q_u
		prev_delta_q_v <- last_q_v
		last_delta_u <- abs(q_u - prev_delta_q_u)
		last_delta_v <- abs(q_v - prev_delta_q_v)
		last_delta_summary <- .lame_snap_delta_summary(
			last_delta_u, last_delta_v, q_u, q_v,
			prev_delta_q_u, prev_delta_q_v, symmetric, threshold,
			snap_stability_tol, snap_delta_quantile)
		q_delta <- last_delta_summary$max
		if (!is.finite(q_delta)) q_delta <- 0
		class_stable <- is.finite(last_delta_summary$share_class_changed) &&
			last_delta_summary$share_class_changed <= snap_class_change_tol
		snap_delta_for_convergence <- switch(
			snap_convergence,
			max = last_delta_summary$max,
			quantile = last_delta_summary$convergence_delta,
			classification = if (isTRUE(class_stable)) 0 else Inf)
		if (!is.finite(snap_delta_for_convergence)) {
			snap_delta_for_convergence <- Inf
		}
		snap_score_stable <- switch(
			snap_convergence,
			max = snap_delta_for_convergence <= snap_stability_tol &&
				class_stable,
			quantile = snap_delta_for_convergence <= snap_stability_tol &&
				class_stable,
			classification = class_stable)
		param_trace <- rbind(param_trace,
		                     data.frame(iter = iter, objective = obj,
		                                rho_uv = rho_uv, sigma_uv = sigma_uv,
		                                pi_snap = pi_snap, obs_var = obs_var,
		                                max_snap_delta = q_delta,
		                                snap_delta_q95 = last_delta_summary$q95,
		                                snap_delta_q99 = last_delta_summary$q99,
		                                snap_delta_conv = snap_delta_for_convergence,
		                                mean_snap_delta = last_delta_summary$mean,
		                                n_snap_delta_gt_tol = last_delta_summary$n_unstable,
		                                share_snap_delta_gt_tol = last_delta_summary$share_unstable,
		                                n_snap_class_changed = last_delta_summary$n_class_changed,
		                                share_snap_class_changed = last_delta_summary$share_class_changed))
		if (iter > 2L && is.finite(prev_obj) &&
		    abs(prev_obj - obj) / (abs(prev_obj) + 1) < tol &&
		    snap_score_stable) {
			converged <- TRUE
			break
		}
		prev_obj <- obj
		last_q_u <- q_u
		last_q_v <- q_v
	}
	fitbits <- .lame_snap_fitted_lists(prep$Yarr, Xb, U, V, mu, a, b, mask,
	                                   symmetric, prep$row_names,
	                                   prep$col_names, prep$time_names)
	coefficients <- c(intercept = mu, beta)
	if (length(beta) > 0L) {
		names(beta) <- prep$x_names_dyad
		names(coefficients)[-1L] <- prep$x_names_dyad
	}
	names(a) <- prep$row_names
	names(b) <- prep$col_names
	dimnames(U) <- list(prep$row_names, paste0("dim", seq_len(R)), prep$time_names)
	dimnames(V) <- list(prep$col_names, paste0("dim", seq_len(R)), prep$time_names)
	snap_class <- q_u > threshold
	snap_class[, 1L] <- NA
	if (!symmetric) {
		snap_class_v <- q_v > threshold
		snap_class_v[, 1L] <- NA
	} else {
		snap_class_v <- NULL
	}
	if (is.null(last_delta_summary)) {
		last_delta_summary <- .lame_snap_delta_summary(
			last_delta_u, last_delta_v, q_u, q_v,
			prev_delta_q_u, prev_delta_q_v, symmetric, threshold,
			snap_stability_tol, snap_delta_quantile)
	}
	unstable_transitions <- .lame_snap_top_unstable(
		last_delta_u, last_delta_v, q_u, q_v,
		prev_delta_q_u, prev_delta_q_v, symmetric, threshold,
		unstable_top_n)
	final_max_snap_delta <- if (nrow(param_trace) > 0L) {
		utils::tail(param_trace$max_snap_delta, 1L)
	} else {
		NA_real_
	}
	final_snap_delta_conv <- if (nrow(param_trace) > 0L) {
		utils::tail(param_trace$snap_delta_conv, 1L)
	} else {
		NA_real_
	}
	final_share_class_changed <- if (nrow(param_trace) > 0L) {
		utils::tail(param_trace$share_snap_class_changed, 1L)
	} else {
		NA_real_
	}
	max_snap_stable <- if (is.finite(final_max_snap_delta)) {
		final_max_snap_delta <= snap_stability_tol
	} else {
		NA
	}
	class_stable <- if (is.finite(final_share_class_changed)) {
		final_share_class_changed <= snap_class_change_tol
	} else {
		NA
	}
	snap_stable <- if (identical(snap_convergence, "classification")) {
		class_stable
	} else if (is.finite(final_snap_delta_conv) && !is.na(class_stable)) {
		final_snap_delta_conv <= snap_stability_tol && class_stable
	} else {
		NA
	}
	objective_diff <- diff(objective_trace)
	objective_increase_tol <- tol * (abs(utils::head(objective_trace, -1L)) + 1)
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
	transition_diagnostics <- list(
		actor_presence = panel$actor_presence,
		balanced_actor_set = panel$balanced_actor_set,
		eligible = eligible_u,
		eligible_v = if (symmetric) NULL else eligible_v,
		last_snap_delta = last_delta_u,
		last_snap_delta_v = if (symmetric) NULL else last_delta_v,
		unstable_transitions = unstable_transitions,
		drift_dist2 = qu$drift_d2,
		snap_dist2 = qu$snap_d2,
		logit_score = qu$logit,
		profile_drift_obj = qu$profile_drift_obj,
		profile_snap_obj = qu$profile_snap_obj,
		profile_delta_obj = qu$profile_delta_obj,
		drift_dist2_v = if (symmetric) NULL else qv$drift_d2,
		snap_dist2_v = if (symmetric) NULL else qv$snap_d2,
		logit_score_v = if (symmetric) NULL else qv$logit,
		profile_drift_obj_v = if (symmetric) NULL else qv$profile_drift_obj,
		profile_snap_obj_v = if (symmetric) NULL else qv$profile_snap_obj,
		profile_delta_obj_v = if (symmetric) NULL else qv$profile_delta_obj)
	out <- list(
		call = match.call(),
		method = "als_snap",
		family = "normal",
		mode = "unipartite",
		symmetric = symmetric,
		R = R,
		longitudinal = TRUE,
		mu = mu,
		beta = beta,
		a = a,
		b = b,
		U = U,
		V = V,
		L = NULL,
		G = NULL,
		coefficients = coefficients,
		VC = c(va = stats::var(a), cab = if (symmetric) NA_real_ else stats::cov(a, b),
		       vb = if (symmetric) NA_real_ else stats::var(b), rho = 0,
		       ve = fitbits$s2),
		s2 = fitbits$s2,
		EZ = fitbits$EZ,
		fitted = fitbits$fitted,
		YPM = fitbits$fitted,
		residuals = fitbits$residuals,
		deviance = fitbits$deviance,
		objective_trace = objective_trace,
		param_trace = param_trace,
		iterations = nrow(param_trace),
		converged = converged,
		convergence = list(converged = converged,
		                   iterations = nrow(param_trace),
		                   tolerance = tol,
		                   snap_stability_tol = snap_stability_tol,
		                   snap_convergence = snap_convergence,
		                   snap_delta_quantile = snap_delta_quantile,
		                   snap_class_change_tol = snap_class_change_tol,
		                   final_max_snap_delta = final_max_snap_delta,
		                   final_snap_delta_q95 = last_delta_summary$q95,
		                   final_snap_delta_q99 = last_delta_summary$q99,
		                   final_snap_delta_conv = final_snap_delta_conv,
		                   final_mean_snap_delta = last_delta_summary$mean,
		                   final_n_snap_delta_gt_tol = last_delta_summary$n_unstable,
		                   final_share_snap_delta_gt_tol = last_delta_summary$share_unstable,
		                   final_n_snap_class_changed = last_delta_summary$n_class_changed,
		                   final_share_snap_class_changed = final_share_class_changed,
		                   snap_stable = snap_stable,
		                   max_snap_stable = max_snap_stable,
		                   snap_class_stable = class_stable,
		                   unstable_transitions = unstable_transitions,
		                   objective_increases = objective_increases,
		                   max_objective_increase = max_objective_increase,
		                   final_objective = utils::tail(objective_trace, 1)),
		rho_uv = rho_uv,
		sigma_uv = sigma_uv,
		obs_var = fitbits$s2,
		pi_snap = pi_snap,
		snap_prob = q_u,
		snap_prob_v = if (symmetric) NULL else q_v,
		snap_class = snap_class,
		snap_class_v = snap_class_v,
		transition_diagnostics = transition_diagnostics,
		n_time = prep$Tt,
			dims = list(n_row = prep$nr, n_col = prep$nc, p = length(beta)),
			row_names = prep$row_names,
			col_names = prep$col_names,
			x_names = names(beta),
			time_names = prep$time_names,
		actor_presence = panel$actor_presence,
		balanced_actor_set = panel$balanced_actor_set,
		Y = prep$Yarr,
		X = prep$X,
		W_row = NULL,
		W_col = NULL,
		dynamic_uv = TRUE,
		dynamic_ab = FALSE,
		dynamic_beta = FALSE,
		dynamic_G = FALSE,
		dynamic = list(uv = TRUE, uv_kind = "snap", ab = FALSE,
		               beta = FALSE, beta_kind = NULL, G = FALSE,
		               snap_model_estimated = TRUE,
		               t_model_estimated = FALSE,
		               mode = "unipartite"),
		meta = list(
			sampler = "als_snap",
			algorithm = "profiled block coordinate descent / EM snap-shift ALS",
			align = align,
			hyper_update = hyper_update,
			snap_kappa = snap_kappa,
			snap_pi_prior = snap_pi_prior,
			rho_prior_mean = rho_prior_mean,
			rho_prior_weight = rho_prior_weight,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions,
			sigma_update_floor = sigma_update_floor,
			snap_convergence = snap_convergence,
			snap_delta_quantile = snap_delta_quantile,
			snap_class_change_tol = snap_class_change_tol,
			uncertainty_available = FALSE,
			approximation_note = paste0(
				"Fast dynamic snap-shift ALS fit. Snap values are ALS ",
				"scores, not posterior draws.")))
	class(out) <- c("lame_snap_als", "lame_als", "ame_als")
	if (verbose && !isTRUE(out$convergence$snap_stable)) {
		cli::cli_warn(c(
			"{.fn lame_snap_als} stopped while snap scores were still moving.",
			"i" = "Final {out$convergence$snap_convergence} delta: {round(out$convergence$final_snap_delta_conv, 4)}; tolerance: {out$convergence$snap_stability_tol}.",
			"i" = "{out$convergence$final_n_snap_class_changed} actor-period classification{?s} changed in the final step.",
			"i" = "See {.code fit$convergence$unstable_transitions} for the largest moves."))
	}
	if (verbose && isTRUE(out$convergence$snap_stable) &&
	    !isTRUE(out$convergence$max_snap_stable)) {
		cli::cli_warn(c(
			"The main snap summary is stable, but a few actor-period scores still moved.",
			"i" = "Final max snap-score delta was {round(out$convergence$final_max_snap_delta, 4)}.",
			"i" = "See {.code fit$convergence$unstable_transitions} for the largest moves."))
	}
	if (verbose && out$convergence$objective_increases > 0L) {
		cli::cli_warn(c(
			"{.fn lame_snap_als} objective was not monotone in {out$convergence$objective_increases} step{?s}.",
			"i" = "Largest objective increase was {round(out$convergence$max_objective_increase, 4)}; inspect {.code fit$param_trace}."))
	}
	.lame_snap_attach_stability(out, stability, stability_args,
	                            seed = seed, verbose = verbose)
}

.lame_snap_format_param <- function(x, digits = 4) {
	if (length(x) == 0L) return("")
	val <- format(round(x, digits), trim = TRUE)
	nm <- names(x)
	if (!is.null(nm) && any(nzchar(nm))) {
		paste0(nm, "=", val, collapse = ", ")
	} else {
		paste(val, collapse = ", ")
	}
}

.lame_snap_score_vector <- function(fit) {
	c(.lame_snap_matrix_vector(fit$snap_prob, "row"),
	  .lame_snap_matrix_vector(fit$snap_prob_v, "col"))
}

.lame_snap_class_vector <- function(fit) {
	c(.lame_snap_matrix_vector(fit$snap_class, "row"),
	  .lame_snap_matrix_vector(fit$snap_class_v, "col"))
}

.lame_snap_matrix_vector <- function(M, side) {
	if (is.null(M)) return(numeric(0))
	out <- as.vector(M)
	rn <- rownames(M)
	cn <- colnames(M)
	if (!is.null(rn) && !is.null(cn)) {
		names(out) <- as.vector(outer(rn, cn, function(a, tt) {
			paste(side, a, tt, sep = "::")
		}))
	}
	out
}

.lame_snap_fitted_vector <- function(fit) {
	ff <- fit$fitted
	if (is.null(ff)) ff <- fit$YPM
	if (is.null(ff)) return(numeric(0))
	if (is.list(ff)) {
		out <- unlist(lapply(seq_along(ff), function(t) {
			M <- ff[[t]]
			val <- as.vector(M)
			rn <- rownames(M)
			cn <- colnames(M)
			tn <- names(ff)[t]
			if (is.null(tn) || !nzchar(tn)) tn <- paste0("t", t)
			if (!is.null(rn) && !is.null(cn)) {
				names(val) <- as.vector(outer(rn, cn, function(a, b) {
					paste(a, b, tn, sep = "::")
				}))
			}
			val
		}), use.names = TRUE)
		return(out)
	}
	as.vector(ff)
}

.lame_snap_period_means <- function(fit) {
	time_names <- fit$time_names
	if (is.null(time_names)) {
		time_names <- paste0("t", seq_len(fit$n_time %||% ncol(fit$snap_prob)))
	}
	sums <- rep(0, length(time_names))
	counts <- rep(0, length(time_names))
	names(sums) <- names(counts) <- time_names
	add_scores <- function(M) {
		if (is.null(M)) return()
		for (t in seq_len(ncol(M))) {
			vals <- M[, t]
			keep <- is.finite(vals)
			if (any(keep)) {
				sums[t] <<- sums[t] + sum(vals[keep])
				counts[t] <<- counts[t] + sum(keep)
			}
		}
	}
	add_scores(fit$snap_prob)
	add_scores(fit$snap_prob_v)
	out <- sums / counts
	out[counts == 0] <- NA_real_
	out
}

.lame_snap_compare_fits <- function(base, alt, threshold = 0.5) {
	score_base <- .lame_snap_score_vector(base)
	score_alt <- .lame_snap_score_vector(alt)
	if (!is.null(names(score_base)) && !is.null(names(score_alt))) {
		common_score <- intersect(names(score_base), names(score_alt))
		score_base <- score_base[common_score]
		score_alt <- score_alt[common_score]
	} else {
		n <- min(length(score_base), length(score_alt))
		score_base <- score_base[seq_len(n)]
		score_alt <- score_alt[seq_len(n)]
	}
	keep <- is.finite(score_base) & is.finite(score_alt)
	diff <- abs(score_base[keep] - score_alt[keep])
	snap_cor <- if (sum(keep) > 2L &&
	    stats::sd(score_base[keep]) > 0 && stats::sd(score_alt[keep]) > 0) {
		stats::cor(score_base[keep], score_alt[keep])
	} else {
		NA_real_
	}
	class_base <- .lame_snap_class_vector(base)
	class_alt <- .lame_snap_class_vector(alt)
	if (!is.null(names(class_base)) && !is.null(names(class_alt))) {
		common_class <- intersect(names(class_base), names(class_alt))
		class_base <- class_base[common_class]
		class_alt <- class_alt[common_class]
	} else {
		nc <- min(length(class_base), length(class_alt))
		class_base <- class_base[seq_len(nc)]
		class_alt <- class_alt[seq_len(nc)]
	}
	class_keep <- !is.na(class_base) & !is.na(class_alt)
	class_change_share <- if (any(class_keep)) {
		mean(class_base[class_keep] != class_alt[class_keep])
	} else {
		NA_real_
	}
	f_base <- .lame_snap_fitted_vector(base)
	f_alt <- .lame_snap_fitted_vector(alt)
	fitted_rmse <- NA_real_
	if (!is.null(names(f_base)) && !is.null(names(f_alt))) {
		common_f <- intersect(names(f_base), names(f_alt))
		f_base <- f_base[common_f]
		f_alt <- f_alt[common_f]
	} else {
		nf <- min(length(f_base), length(f_alt))
		f_base <- f_base[seq_len(nf)]
		f_alt <- f_alt[seq_len(nf)]
	}
	if (length(f_base) > 0L) {
		f_keep <- is.finite(f_base) & is.finite(f_alt)
		if (any(f_keep)) {
			fitted_rmse <- sqrt(mean((f_base[f_keep] - f_alt[f_keep])^2))
		}
	}
	base_period <- .lame_snap_period_means(base)
	alt_period <- .lame_snap_period_means(alt)
	period_keep <- is.finite(base_period) & is.finite(alt_period)
	top_base <- top_alt <- NA_character_
	top_period_agrees <- NA
	if (any(period_keep)) {
		top_base <- names(base_period)[which.max(ifelse(period_keep, base_period, -Inf))]
		top_alt <- names(alt_period)[which.max(ifelse(period_keep, alt_period, -Inf))]
		top_period_agrees <- identical(top_base, top_alt)
	}
	data.frame(
		n_snap_compared = sum(keep),
		snap_cor = snap_cor,
		mean_abs_snap_diff = if (length(diff) > 0L) mean(diff) else NA_real_,
		q95_abs_snap_diff = if (length(diff) > 0L) {
			as.numeric(stats::quantile(diff, 0.95, na.rm = TRUE, names = FALSE))
		} else {
			NA_real_
		},
		max_abs_snap_diff = if (length(diff) > 0L) max(diff) else NA_real_,
		class_change_share = class_change_share,
		fitted_rmse = fitted_rmse,
		top_period_baseline = top_base,
		top_period_run = top_alt,
		top_period_agrees = top_period_agrees,
		stringsAsFactors = FALSE)
}

.lame_snap_start_stability <- function(fit, args, preset, seed) {
	n_runs <- switch(preset, quick = 2L, validation = 4L, 0L)
	base_seed <- if (length(seed) == 1L && is.finite(seed)) as.integer(seed) else 6886L
	seeds <- as.integer(base_seed + 7919L * seq_len(n_runs))
	runs <- vector("list", n_runs)
	for (i in seq_len(n_runs)) {
		run_args <- args
		run_args$seed <- seeds[i]
		run_args$stability <- "none"
		run_args$verbose <- FALSE
		alt <- tryCatch(do.call(lame_snap_als, run_args),
		                error = function(e) e)
		if (inherits(alt, "error")) {
			runs[[i]] <- data.frame(
				run = i, seed = seeds[i], ok = FALSE,
				error = conditionMessage(alt),
				n_snap_compared = 0L,
				snap_cor = NA_real_, mean_abs_snap_diff = NA_real_,
				q95_abs_snap_diff = NA_real_, max_abs_snap_diff = NA_real_,
				class_change_share = NA_real_, fitted_rmse = NA_real_,
				top_period_baseline = NA_character_,
				top_period_run = NA_character_,
				top_period_agrees = NA,
				converged = NA, iterations = NA_integer_,
				stringsAsFactors = FALSE)
		} else {
			cmp <- .lame_snap_compare_fits(fit, alt)
			runs[[i]] <- cbind(
				data.frame(run = i, seed = seeds[i], ok = TRUE,
				           error = NA_character_, stringsAsFactors = FALSE),
				cmp,
				data.frame(converged = isTRUE(alt$converged),
				           iterations = alt$iterations,
				           stringsAsFactors = FALSE))
		}
	}
	tab <- do.call(rbind, runs)
	ok <- tab$ok
	comparable_ok <- !is.na(tab$n_snap_compared) & tab$n_snap_compared > 0L
	score_ok <- is.na(tab$snap_cor) | tab$snap_cor >= 0.80
	class_ok <- !is.na(tab$class_change_share) &
		tab$class_change_share <= 0.10
	q95_ok <- !is.na(tab$q95_abs_snap_diff) &
		tab$q95_abs_snap_diff <= 0.25
	top_ok <- !is.na(tab$top_period_agrees) & tab$top_period_agrees
	stable <- any(ok) &&
		all(ok & comparable_ok & score_ok & class_ok & q95_ok & top_ok)
	list(
		preset = preset,
		baseline_seed = base_seed,
		seeds = seeds,
		runs = tab,
		stable = stable,
		criteria = list(min_snap_cor = 0.80,
		                max_class_change_share = 0.10,
		                max_q95_abs_snap_diff = 0.25,
		                require_top_period_agreement = TRUE))
}

.lame_snap_attach_stability <- function(fit, stability, args, seed, verbose) {
	fit$meta$stability <- stability
	if (identical(stability, "none")) return(fit)
	st <- .lame_snap_start_stability(fit, args, stability, seed)
	fit$stability <- st
	fit$convergence$start_stable <- st$stable
	if (verbose && !isTRUE(st$stable)) {
		cli::cli_warn(c(
			"{.fn lame_snap_als} start-sensitivity check found movement across starts.",
			"i" = "Inspect {.code fit$stability$runs}; use {.arg stability = \"validation\"} and a larger {.arg max_iter} for a slower check."))
	}
	fit
}

#' @export
print.lame_snap_als <- function(x, digits = 4, ...) {
	cli::cli_text("{.strong Dynamic snap-shift AME fit via ALS} (fast approximation)")
	rank_text <- if (identical(x$mode, "bipartite")) {
		paste0("R_row = ", x$R_row, ", R_col = ", x$R_col)
	} else {
		paste0("R = ", x$R)
	}
	cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}}{ifelse(x$symmetric, ' (symmetric)', '')} | {rank_text}")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} periods")
	cli::cli_text("Converged: {.val {x$converged}} in {x$iterations} iteration{?s}")
	rho_text <- .lame_snap_format_param(x$rho_uv, digits)
	sigma_text <- .lame_snap_format_param(x$sigma_uv, digits)
	pi_text <- .lame_snap_format_param(x$pi_snap, digits)
	cli::cli_text("rho_uv = {.val {rho_text}}, sigma_uv = {.val {sigma_text}}, pi_snap = {.val {pi_text}}")
	all_scores <- c(as.vector(x$snap_prob), as.vector(x$snap_prob_v))
	eligible <- !is.na(all_scores)
	hi <- sum(all_scores[eligible] > 0.5, na.rm = TRUE)
	cli::cli_text("Snap scores > 0.5: {hi} of {sum(eligible)} eligible actor-periods")
	if (!isTRUE(x$convergence$snap_stable)) {
		msg <- paste0("Snap scores are still moving: final ",
		              x$convergence$snap_convergence, " delta = ",
		              round(x$convergence$final_snap_delta_conv, digits),
		              " (tolerance ", x$convergence$snap_stability_tol, ").")
		cli::cli_text(cli::col_yellow(msg))
		cli::cli_text(cli::col_yellow("Final max delta = {round(x$convergence$final_max_snap_delta, digits)}; class changes = {x$convergence$final_n_snap_class_changed}."))
		cli::cli_text(cli::col_yellow("See {.code fit$convergence$unstable_transitions} for the largest moves."))
	} else if (!isTRUE(x$convergence$max_snap_stable)) {
		cli::cli_text(cli::col_yellow("The main snap summary is stable, but final max delta = {round(x$convergence$final_max_snap_delta, digits)}."))
		cli::cli_text(cli::col_yellow("See {.code fit$convergence$unstable_transitions} for the largest moves."))
	}
	if (!is.null(x$stability)) {
		cli::cli_text("Start sensitivity: {.val {x$stability$stable}} ({x$stability$preset}, {nrow(x$stability$runs)} rerun{?s})")
	}
	cli::cli_text(cli::col_grey("Snap scores are ALS scores, not posterior draws."))
	invisible(x)
}

#' @export
summary.lame_snap_als <- function(object, ...) {
	tab <- cbind(Estimate = object$coefficients)
	all_scores <- c(as.vector(object$snap_prob), as.vector(object$snap_prob_v))
	finite_scores <- all_scores[is.finite(all_scores)]
	out <- list(
		call = object$call,
		family = object$family,
		mode = object$mode,
		symmetric = object$symmetric,
		R = object$R,
		R_row = object$R_row,
		R_col = object$R_col,
		n_time = object$n_time,
		dims = object$dims,
		coefficients = tab,
		VC = object$VC,
		rho_uv = object$rho_uv,
		sigma_uv = object$sigma_uv,
		pi_snap = object$pi_snap,
		convergence = object$convergence,
		final_max_snap_delta = object$convergence$final_max_snap_delta,
		final_snap_delta_q95 = object$convergence$final_snap_delta_q95,
		final_snap_delta_q99 = object$convergence$final_snap_delta_q99,
		final_snap_delta_conv = object$convergence$final_snap_delta_conv,
		final_n_snap_delta_gt_tol = object$convergence$final_n_snap_delta_gt_tol,
		final_n_snap_class_changed = object$convergence$final_n_snap_class_changed,
		final_share_snap_class_changed = object$convergence$final_share_snap_class_changed,
		snap_stability_tol = object$convergence$snap_stability_tol,
		snap_convergence = object$convergence$snap_convergence,
		unstable_transitions = object$convergence$unstable_transitions,
		stability = object$stability,
		objective_increases = object$convergence$objective_increases,
		max_objective_increase = object$convergence$max_objective_increase,
		mean_snap = if (length(finite_scores) > 0L) mean(finite_scores) else NA_real_,
		max_snap = if (length(finite_scores) > 0L) max(finite_scores) else NA_real_,
		approximation_note = object$meta$approximation_note)
	class(out) <- "summary.lame_snap_als"
	out
}

#' @export
print.summary.lame_snap_als <- function(x, digits = 4, ...) {
	cli::cli_h2("Dynamic snap-shift ALS summary")
	rank_text <- if (identical(x$mode, "bipartite")) {
		paste0("R_row = ", x$R_row, ", R_col = ", x$R_col)
	} else {
		paste0("R = ", x$R)
	}
	cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}}{ifelse(x$symmetric, ' (symmetric)', '')} | {rank_text}")
	cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} periods")
	cli::cli_text("Converged: {.val {x$convergence$converged}} ({x$convergence$iterations} iteration{?s})")
	rho_text <- .lame_snap_format_param(x$rho_uv, digits)
	sigma_text <- .lame_snap_format_param(x$sigma_uv, digits)
	pi_text <- .lame_snap_format_param(x$pi_snap, digits)
	cli::cli_text("rho_uv = {.val {rho_text}}, sigma_uv = {.val {sigma_text}}, pi_snap = {.val {pi_text}}")
	cli::cli_text("Mean snap score = {round(x$mean_snap, digits)}, max = {round(x$max_snap, digits)}")
	if (!isTRUE(x$convergence$snap_stable)) {
		cli::cli_text(cli::col_yellow("Snap scores are still moving; final {x$snap_convergence} delta = {round(x$final_snap_delta_conv, digits)} (tolerance {x$snap_stability_tol})."))
		cli::cli_text(cli::col_yellow("Final max delta = {round(x$final_max_snap_delta, digits)}; class changes = {x$final_n_snap_class_changed}."))
		if (!is.null(x$unstable_transitions) && nrow(x$unstable_transitions) > 0L) {
			cli::cli_text(cli::col_yellow("Largest moving actor-periods:"))
			print(utils::head(x$unstable_transitions, 5L), row.names = FALSE)
		}
	} else if (!isTRUE(x$convergence$max_snap_stable)) {
		cli::cli_text(cli::col_yellow("The main snap summary is stable, but final max delta = {round(x$final_max_snap_delta, digits)}."))
	}
	if (!is.null(x$stability)) {
		cli::cli_text("Start sensitivity: {.val {x$stability$stable}} ({x$stability$preset}, {nrow(x$stability$runs)} rerun{?s})")
		if (!isTRUE(x$stability$stable)) {
			print(x$stability$runs, row.names = FALSE)
		}
	}
	if (x$objective_increases > 0L) {
		cli::cli_text(cli::col_yellow("Objective was not monotone in {x$objective_increases} step{?s}; max increase = {round(x$max_objective_increase, digits)}."))
	}
	cli::cli_rule()
	cli::cli_text("{.strong Regression coefficients}")
	print(round(x$coefficients, digits))
	cli::cli_text(cli::col_grey(x$approximation_note))
	invisible(x)
}
