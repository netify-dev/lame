.lame_bsnap_rect_identity <- function(nr, nc) {
	out <- matrix(0, nr, nc)
	for (k in seq_len(min(nr, nc))) out[k, k] <- 1
	out
}

.lame_bsnap_uvg_cube <- function(U, V, G) {
	du <- dim(U)
	dv <- dim(V)
	Tt <- du[3L]
	out <- array(0, dim = c(du[1L], dv[1L], Tt),
	             dimnames = list(dimnames(U)[[1L]],
	                             dimnames(V)[[1L]],
	                             dimnames(U)[[3L]]))
	for (t in seq_len(Tt)) {
		Ut <- matrix(U[, , t], du[1L], du[2L])
		Vt <- matrix(V[, , t], dv[1L], dv[2L])
		out[, , t] <- Ut %*% G %*% t(Vt)
	}
	out
}

.lame_bsnap_normalize_uvg <- function(U, V, G) {
	su <- sqrt(mean(U^2, na.rm = TRUE))
	sv <- sqrt(mean(V^2, na.rm = TRUE))
	if (is.finite(su) && su > 1e-8) {
		U <- U / su
		G <- su * G
	}
	if (is.finite(sv) && sv > 1e-8) {
		V <- V / sv
		G <- G * sv
	}
	list(U = U, V = V, G = G)
}

.lame_bsnap_initial_factors <- function(prep, mu, beta, a, b,
                                        R_row, R_col,
                                        align = "sequential",
                                        seed = 6886) {
	set.seed(seed)
	nr <- prep$nr
	nc <- prep$nc
	Tt <- prep$Tt
	U <- array(0, dim = c(nr, R_row, Tt),
	           dimnames = list(prep$row_names,
	                           paste0("row_dim", seq_len(R_row)),
	                           prep$time_names))
	V <- array(0, dim = c(nc, R_col, Tt),
	           dimnames = list(prep$col_names,
	                           paste0("col_dim", seq_len(R_col)),
	                           prep$time_names))
	G <- .lame_bsnap_rect_identity(R_row, R_col)
	dimnames(G) <- list(dimnames(U)[[2L]], dimnames(V)[[2L]])
	Xb <- .lame_snap_xbeta(beta, prep$X)
	ab <- outer(a, b, "+")
	Uref <- NULL
	Vref <- NULL
	for (t in seq_len(Tt)) {
		M <- prep$Yarr[, , t] - mu - Xb[, , t] - ab
		M[!is.finite(M)] <- 0
		sv <- tryCatch(svd(M), error = function(e) NULL)
		Ut <- matrix(stats::rnorm(nr * R_row, 0, 1e-3), nr, R_row)
		Vt <- matrix(stats::rnorm(nc * R_col, 0, 1e-3), nc, R_col)
		if (!is.null(sv) && length(sv$d) > 0L) {
			rr <- min(length(sv$d), R_row, R_col)
			if (rr > 0L) {
				dd <- sqrt(pmax(sv$d[seq_len(rr)], 0))
				Ut[, seq_len(rr)] <- sweep(sv$u[, seq_len(rr), drop = FALSE],
				                           2, dd, "*")
				Vt[, seq_len(rr)] <- sweep(sv$v[, seq_len(rr), drop = FALSE],
				                           2, dd, "*")
			}
		}
		if (identical(align, "sequential") && t > 1L) {
			Ut <- .lame_snap_procrustes(Ut, Uref)
			Vt <- .lame_snap_procrustes(Vt, Vref)
		} else if (identical(align, "global") && t > 1L) {
			Ut <- .lame_snap_procrustes(Ut, matrix(U[, , 1L], nr, R_row))
			Vt <- .lame_snap_procrustes(Vt, matrix(V[, , 1L], nc, R_col))
		}
		U[, , t] <- Ut
		V[, , t] <- Vt
		Uref <- Ut
		Vref <- Vt
	}
	n <- .lame_bsnap_normalize_uvg(U, V, G)
	list(U = n$U, V = n$V, G = n$G)
}

.lame_bsnap_active <- function(mask) {
	nr <- dim(mask)[1L]
	nc <- dim(mask)[2L]
	Tt <- dim(mask)[3L]
	row_active <- matrix(FALSE, nr, Tt)
	col_active <- matrix(FALSE, nc, Tt)
	row_exposure <- matrix(0, nr, Tt)
	col_exposure <- matrix(0, nc, Tt)
	for (t in seq_len(Tt)) {
		row_exposure[, t] <- rowSums(mask[, , t], na.rm = TRUE)
		col_exposure[, t] <- colSums(mask[, , t], na.rm = TRUE)
		row_active[, t] <- row_exposure[, t] > 0
		col_active[, t] <- col_exposure[, t] > 0
	}
	list(row = row_active, col = col_active,
	     row_exposure = row_exposure, col_exposure = col_exposure)
}

.lame_bsnap_static_update <- function(Yarr, X, U, V, G, mu, beta, a, b,
                                      mask, linear_solver = "eigen") {
	d <- dim(Yarr)
	nr <- d[1L]
	nc <- d[2L]
	Tt <- d[3L]
	p <- dim(X)[3L]
	Warr <- array(0, dim = d)
	Warr[mask] <- 1
	Yf <- Yarr
	Yf[!mask] <- 0
	O <- .lame_bsnap_uvg_cube(U, V, G)
	lin <- .ae_linsolver(linear_solver, X, mask, Warr, nr, nc, Tt, p)
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
	ab_fit <- .ae_additive_solve(base_rs, base_cs, wsum, cnt_row, cnt_col,
	                             symmetric = FALSE)
	a <- ab_fit$a
	b <- ab_fit$b
	mu <- mu + mean(a) + mean(b)
	a <- a - mean(a)
	b <- b - mean(b)
	list(mu = mu, beta = beta, a = a, b = b, Xb = .lame_snap_xbeta(beta, X))
}

.lame_bsnap_temporal_terms <- function(A, cvec, Xcube, q, i, t,
                                       rho, sigma, snap_kappa, eligible) {
	R <- dim(Xcube)[2L]
	sigma2 <- sigma * sigma
	kappa2 <- snap_kappa * snap_kappa
	if (t > 1L && isTRUE(eligible[i, t])) {
		qt <- .lame_snap_q_value(q, i, t)
		A <- A + diag((1 - qt) / sigma2 + qt / kappa2, R)
		cvec <- cvec + (1 - qt) * rho * Xcube[i, , t - 1L] / sigma2
	}
	if (t < dim(Xcube)[3L] && isTRUE(eligible[i, t + 1L])) {
		qf <- .lame_snap_q_value(q, i, t + 1L)
		A <- A + diag((1 - qf) * rho * rho / sigma2, R)
		cvec <- cvec + (1 - qf) * rho * Xcube[i, , t + 1L] / sigma2
	}
	list(A = A, c = cvec)
}

.lame_bsnap_update_uvg <- function(Yarr, Xb, U, V, G, mu, a, b, mask,
                                   q_u, q_v, rho_u, sigma_u, rho_v, sigma_v,
                                   snap_kappa, eligible_u, eligible_v,
                                   ridge = 1e-8, obs_var = 1) {
	nr <- dim(Yarr)[1L]
	nc <- dim(Yarr)[2L]
	RA <- dim(U)[2L]
	RB <- dim(V)[2L]
	Tt <- dim(Yarr)[3L]
	w_obs <- 1 / max(obs_var, 1e-10)
	for (t in seq_len(Tt)) {
		for (i in seq_len(nr)) {
			A <- diag(ridge, RA)
			cvec <- numeric(RA)
			obs <- which(mask[i, , t])
			if (length(obs) > 0L) {
				for (j in obs) {
					x <- as.numeric(G %*% V[j, , t])
					y <- Yarr[i, j, t] - mu - Xb[i, j, t] - a[i] - b[j]
					A <- A + w_obs * tcrossprod(x)
					cvec <- cvec + w_obs * y * x
				}
			}
			tt <- .lame_bsnap_temporal_terms(A, cvec, U, q_u, i, t,
			                                 rho_u, sigma_u, snap_kappa,
			                                 eligible_u)
			U[i, , t] <- .ae_safe_solve(tt$A, tt$c)
		}
		for (j in seq_len(nc)) {
			A <- diag(ridge, RB)
			cvec <- numeric(RB)
			obs <- which(mask[, j, t])
			if (length(obs) > 0L) {
				for (i in obs) {
					x <- as.numeric(t(G) %*% U[i, , t])
					y <- Yarr[i, j, t] - mu - Xb[i, j, t] - a[i] - b[j]
					A <- A + w_obs * tcrossprod(x)
					cvec <- cvec + w_obs * y * x
				}
			}
			tt <- .lame_bsnap_temporal_terms(A, cvec, V, q_v, j, t,
			                                 rho_v, sigma_v, snap_kappa,
			                                 eligible_v)
			V[j, , t] <- .ae_safe_solve(tt$A, tt$c)
		}
	}
	target <- array(0, dim = dim(Yarr))
	ab <- outer(a, b, "+")
	for (t in seq_len(Tt)) target[, , t] <- Yarr[, , t] - mu - Xb[, , t] - ab
	W <- array(0, dim = dim(Yarr))
	W[mask] <- w_obs
	G <- .lda_update_static_g_bip(target, mask, W, U, V, G, ridge = ridge)
	n <- .lame_bsnap_normalize_uvg(U, V, G)
	list(U = n$U, V = n$V, G = n$G)
}

.lame_bsnap_resid_s2 <- function(Yarr, Xb, U, V, G, mu, a, b, mask) {
	O <- .lame_bsnap_uvg_cube(U, V, G)
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

.lame_bsnap_score_side <- function(Xcube, q_old, eligible, rho, sigma,
                                   snap_kappa, pi_snap, snap_update,
                                   snap_damping, threshold, iter, max_iter) {
	n <- dim(Xcube)[1L]
	R <- dim(Xcube)[2L]
	Tt <- dim(Xcube)[3L]
	q <- matrix(NA_real_, n, Tt, dimnames = list(dimnames(Xcube)[[1L]],
	                                            dimnames(Xcube)[[3L]]))
	drift_d2 <- snap_d2 <- logit_score <- q
	sigma2 <- sigma * sigma
	kappa2 <- snap_kappa * snap_kappa
	for (t in 2:Tt) {
		for (i in seq_len(n)) {
			if (!isTRUE(eligible[i, t])) next
			drift <- sum((Xcube[i, , t] - rho * Xcube[i, , t - 1L])^2)
			snap <- sum(Xcube[i, , t]^2)
			score <- .lame_snap_logit(pi_snap) +
				R * log(pmax(sigma, 1e-12) / snap_kappa) +
				0.5 * drift / max(sigma2, 1e-12) -
				0.5 * snap / max(kappa2, 1e-12)
			score <- pmin(pmax(score, -40), 40)
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
		}
	}
	q <- .lame_snap_damp_q(q, q_old, snap_damping)
	list(q = q, drift_d2 = drift_d2, snap_d2 = snap_d2,
	     logit = logit_score)
}

.lame_bsnap_update_hyper_side <- function(Xcube, q, eligible, estimate_rho,
                                          estimate_sigma, rho, sigma,
                                          snap_pi_prior, min_sigma,
                                          rho_prior_mean = 0.98,
                                          rho_prior_weight = 25,
                                          hyper_update = "robust",
                                          drift_quantile = 0.05,
                                          drift_min_transitions = 50,
                                          q_power = 2) {
	R <- dim(Xcube)[2L]
	Tt <- dim(Xcube)[3L]
	num <- den <- ss <- wcount <- qsum <- m <- 0
	transitions <- vector("list", 0L)
	idx <- 0L
	rho_ref <- if (is.finite(rho_prior_mean)) rho_prior_mean else rho
	for (t in 2:Tt) {
		for (i in seq_len(dim(Xcube)[1L])) {
			if (!isTRUE(eligible[i, t])) next
			qq <- .lame_snap_q_value(q, i, t)
			w <- pmax(0, 1 - qq)
			x <- Xcube[i, , t]
			xlag <- Xcube[i, , t - 1L]
			if (!all(is.finite(x)) || !all(is.finite(xlag))) next
			num <- num + w * sum(x * xlag)
			den <- den + w * sum(xlag * xlag)
			ss <- ss + w * sum((x - rho * xlag)^2)
			wcount <- wcount + w
			qsum <- qsum + qq
			m <- m + 1L
			idx <- idx + 1L
			transitions[[idx]] <- list(
				q = qq,
				w = pmax(0, 1 - qq)^q_power,
				num = sum(x * xlag),
				den = sum(xlag * xlag),
				ref_d2 = sum((x - rho_ref * xlag)^2) / R,
				x = x,
				xlag = xlag)
		}
	}
	used_robust <- FALSE
	if (identical(hyper_update, "robust") && length(transitions) > 0L) {
		ref_d2 <- vapply(transitions, `[[`, numeric(1L), "ref_d2")
		w_rob <- vapply(transitions, `[[`, numeric(1L), "w")
		n_eff <- sum(w_rob[is.finite(w_rob) & w_rob > 0])
		if (n_eff >= drift_min_transitions) {
			q_eff <- max(drift_quantile,
			             min(1, drift_min_transitions / max(n_eff, 1)))
			cut <- .lame_snap_weighted_quantile(ref_d2, w_rob, q_eff)
			keep <- is.finite(ref_d2) & is.finite(w_rob) & w_rob > 0 &
				ref_d2 <= cut[1L]
			if (!any(keep)) keep <- is.finite(ref_d2) & is.finite(w_rob) & w_rob > 0
			if (estimate_rho && any(keep)) {
				num_k <- sum(vapply(transitions[keep], `[[`, numeric(1L), "num") *
				             w_rob[keep])
				den_k <- sum(vapply(transitions[keep], `[[`, numeric(1L), "den") *
				             w_rob[keep])
				if (den_k > 0) {
					rho_raw <- num_k / den_k
					n_eff_keep <- sum(w_rob[keep])
					if (is.finite(rho_prior_mean) && rho_prior_weight > 0) {
						rho_raw <- (n_eff_keep * rho_raw +
						            rho_prior_weight * rho_prior_mean) /
							(n_eff_keep + rho_prior_weight)
					}
					rho <- pmin(0.995, pmax(-0.995, rho_raw))
				}
			}
			if (estimate_sigma && any(keep)) {
				innov <- vapply(transitions[keep], function(z) {
					sum((z$x - rho * z$xlag)^2) / R
				}, numeric(1L))
				wk <- w_rob[keep]
				sig2 <- stats::weighted.mean(innov, wk)
				if (is.finite(sig2) && sig2 > 0) sigma <- sqrt(sig2)
			}
			used_robust <- TRUE
		}
	}
	if (!used_robust && estimate_rho && den > 0) {
		rho_raw <- num / den
		if (is.finite(rho_prior_mean) && rho_prior_weight > 0) {
			rho_raw <- (wcount * rho_raw + rho_prior_weight * rho_prior_mean) /
				(wcount + rho_prior_weight)
		}
		rho <- pmin(0.995, pmax(-0.995, rho_raw))
	}
	if (!used_robust && estimate_sigma && wcount > 0) {
		ss2 <- 0
		for (t in 2:Tt) {
			for (i in seq_len(dim(Xcube)[1L])) {
				if (!isTRUE(eligible[i, t])) next
				qq <- .lame_snap_q_value(q, i, t)
				w <- pmax(0, 1 - qq)
				ss2 <- ss2 + w * sum((Xcube[i, , t] - rho * Xcube[i, , t - 1L])^2)
			}
		}
		sigma <- sqrt(ss2 / max(1, wcount * R))
	}
	sigma <- max(min_sigma, sigma)
	a_pi <- snap_pi_prior[1L]
	b_pi <- snap_pi_prior[2L]
	pi_snap <- (a_pi + qsum) / (a_pi + b_pi + m)
	list(rho = as.numeric(rho), sigma = as.numeric(sigma),
	     pi_snap = as.numeric(pi_snap), n = m)
}

.lame_bsnap_surface <- function(U, V, G, mask, rho_u, rho_v, eligible_u,
                                eligible_v) {
	nr <- dim(U)[1L]
	nc <- dim(V)[1L]
	Tt <- dim(U)[3L]
	surf_u <- matrix(NA_real_, nr, Tt, dimnames = list(dimnames(U)[[1L]],
	                                                  dimnames(U)[[3L]]))
	surf_v <- matrix(NA_real_, nc, Tt, dimnames = list(dimnames(V)[[1L]],
	                                                  dimnames(V)[[3L]]))
	for (t in 2:Tt) {
		for (i in seq_len(nr)) {
			if (!isTRUE(eligible_u[i, t])) next
			obs <- which(mask[i, , t])
			if (length(obs) == 0L) next
			d <- U[i, , t] - rho_u * U[i, , t - 1L]
			jump <- vapply(obs, function(j) {
				sum(d * as.numeric(G %*% V[j, , t]))
			}, numeric(1L))
			surf_u[i, t] <- sqrt(mean(jump^2))
		}
		for (j in seq_len(nc)) {
			if (!isTRUE(eligible_v[j, t])) next
			obs <- which(mask[, j, t])
			if (length(obs) == 0L) next
			d <- V[j, , t] - rho_v * V[j, , t - 1L]
			jump <- vapply(obs, function(i) {
				sum(d * as.numeric(t(G) %*% U[i, , t]))
			}, numeric(1L))
			surf_v[j, t] <- sqrt(mean(jump^2))
		}
	}
	list(row = surf_u, col = surf_v)
}

.lame_bsnap_objective <- function(Yarr, Xb, U, V, G, mu, a, b, mask,
                                  q_u, q_v, rho_u, sigma_u, rho_v, sigma_v,
                                  snap_kappa, pi_u, pi_v,
                                  eligible_u, eligible_v, obs_var = 1) {
	O <- .lame_bsnap_uvg_cube(U, V, G)
	ab <- outer(a, b, "+")
	sse <- 0
	for (t in seq_len(dim(Yarr)[3L])) {
		r <- Yarr[, , t] - mu - Xb[, , t] - ab - O[, , t]
		sse <- sse + sum(r[mask[, , t]]^2)
	}
	add_side <- function(Xcube, q, eligible, rho, sigma, pi_snap) {
		val <- 0
		for (t in 2:dim(Xcube)[3L]) {
			for (i in seq_len(dim(Xcube)[1L])) {
				if (!isTRUE(eligible[i, t])) next
				qq <- .lame_snap_q_value(q, i, t)
				drift <- sum((Xcube[i, , t] - rho * Xcube[i, , t - 1L])^2)
				snap <- sum(Xcube[i, , t]^2)
				val <- val + (1 - qq) * drift / (sigma * sigma) +
					qq * snap / (snap_kappa * snap_kappa)
				val <- val - 2 * ((1 - qq) * log1p(-pi_snap) +
					qq * log(pi_snap))
			}
		}
		val
	}
	sse / max(obs_var, 1e-10) +
		add_side(U, q_u, eligible_u, rho_u, sigma_u, pi_u) +
		add_side(V, q_v, eligible_v, rho_v, sigma_v, pi_v)
}

.lame_bsnap_fitted_lists <- function(prep, Xb, U, V, G, mu, a, b, mask) {
	Tt <- prep$Tt
	ab <- outer(a, b, "+")
	O <- .lame_bsnap_uvg_cube(U, V, G)
	EZ <- fitted <- residuals <- vector("list", Tt)
	sse <- 0
	n_obs <- 0
	for (t in seq_len(Tt)) {
		ez <- mu + Xb[, , t] + ab + O[, , t]
		r <- prep$Yarr[, , t] - ez
		sse <- sse + sum(r[mask[, , t]]^2)
		n_obs <- n_obs + sum(mask[, , t])
		dimnames(ez) <- list(prep$row_names, prep$col_names)
		dimnames(r) <- list(prep$row_names, prep$col_names)
		EZ[[t]] <- ez
		fitted[[t]] <- ez
		residuals[[t]] <- r
	}
	names(EZ) <- names(fitted) <- names(residuals) <- prep$time_names
	list(EZ = EZ, fitted = fitted, residuals = residuals,
	     s2 = if (n_obs > 0L) sse / n_obs else NA_real_,
	     deviance = sse)
}

.lame_bsnap_g_diagnostics <- function(G) {
	sv <- tryCatch(svd(G, nu = 0, nv = 0)$d, error = function(e) numeric(0))
	smin <- if (length(sv) > 0L) min(sv) else NA_real_
	smax <- if (length(sv) > 0L) max(sv) else NA_real_
	weak_signal <- is.finite(smax) && smax <= sqrt(.Machine$double.eps)
	condition <- if (is.finite(smax) && is.finite(smin)) {
		if (isTRUE(weak_signal)) Inf else smax / max(smin, .Machine$double.eps * max(1, smax))
	} else {
		NA_real_
	}
	list(singular_values = sv, min_singular_value = smin,
	     condition_number = condition,
	     weak_signal = weak_signal,
	     near_singular = isTRUE(weak_signal) ||
		     (is.finite(condition) && condition > 1e6))
}

.lame_snap_bipartite_als <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
                                     R = 2L, R_row = NULL, R_col = NULL,
                                     family = "normal", max_iter = 200L,
                                     tol = 1e-6, snap_kappa = 2,
                                     snap_pi_prior = c(a = 1, b = 9),
                                     snap_update = c("soft", "hard", "annealed"),
                                     snap_damping = 0.7,
                                     estimate_rho_uv = TRUE,
                                     estimate_sigma_uv = TRUE,
                                     rho_uv = NULL, sigma_uv = NULL,
                                     align = c("sequential", "global", "none"),
                                     threshold = 0.5, min_sigma = 1e-4,
                                     ridge = 1e-8,
                                     snap_stability_tol = 0.05,
	                                     snap_convergence = c("quantile", "max",
	                                                          "classification"),
	                                     snap_delta_quantile = 0.95,
	                                     snap_class_change_tol = 0.005,
	                                     hyper_update = c("robust", "em"),
	                                     drift_quantile = 0.05,
	                                     drift_min_transitions = 50,
	                                     sigma_floor_fraction = 0.75,
	                                     unstable_top_n = 10L,
	                                     rho_prior_mean = 0.98,
	                                     rho_prior_weight = 25,
                                     verbose = TRUE,
                                     seed = 6886) {
	snap_update <- match.arg(snap_update)
	align <- match.arg(align)
	snap_convergence <- match.arg(snap_convergence)
	hyper_update <- match.arg(hyper_update)
	if (!identical(family, "normal")) {
		cli::cli_abort(c(
			"{.fn lame_snap_als} supports bipartite snap ALS only for {.code family = \"normal\"}.",
			"i" = "Binary/poisson bipartite snap ALS require an IRLS layer and are not supported by this estimator."))
	}
	if (!is.null(Xrow) || !is.null(Xcol)) {
		cli::cli_abort(c(
			"{.fn lame_snap_als} does not support node covariates.",
			"i" = "Pass dyadic covariates via {.arg Xdyad}; use {.fn lame} MCMC for node-covariate snap models."))
	}
	R_row <- as.integer(.lame_snap_default(R_row, R))
	R_col <- as.integer(.lame_snap_default(R_col, R))
	if (length(R_row) != 1L || length(R_col) != 1L ||
	    anyNA(c(R_row, R_col)) || R_row < 1L || R_col < 1L) {
		cli::cli_abort("{.arg R_row} and {.arg R_col} must be positive integers.")
	}
	if (identical(snap_update, "hard")) snap_damping <- 1
	prep <- .ae_prepare(Y, Xdyad, NULL, NULL, "bipartite", longitudinal = TRUE,
	                    node_time_policy = "quiet")
	if (prep$Tt < 2L) {
		cli::cli_abort("{.fn lame_snap_als} requires at least two time periods.")
	}
	if (R_row >= prep$nr) {
		cli::cli_abort("{.arg R_row} must be smaller than the number of row actors.")
	}
	if (R_col >= prep$nc) {
		cli::cli_abort("{.arg R_col} must be smaller than the number of column actors.")
	}

	static_fit <- suppressWarnings(lame_als(
		Y, Xdyad = Xdyad, R = 0, family = "normal",
		mode = "bipartite", symmetric = FALSE,
		max_iter = max(50L, min(100L, max_iter)),
		verbose = FALSE, seed = seed))
	mu <- static_fit$mu
	beta <- static_fit$beta[seq_len(prep$p_dyad)]
	if (prep$p_dyad == 0L) beta <- numeric(0)
	a <- static_fit$a
	b <- static_fit$b
	init <- .lame_bsnap_initial_factors(prep, mu, beta, a, b,
	                                    R_row, R_col, align, seed)
	U <- init$U
	V <- init$V
	G <- init$G
	Xb <- .lame_snap_xbeta(beta, prep$X)
	mask <- is.finite(prep$Yarr)
	active <- .lame_bsnap_active(mask)
	eligible_u <- matrix(FALSE, prep$nr, prep$Tt,
	                     dimnames = list(prep$row_names, prep$time_names))
	eligible_v <- matrix(FALSE, prep$nc, prep$Tt,
	                     dimnames = list(prep$col_names, prep$time_names))
	for (t in 2:prep$Tt) {
		eligible_u[, t] <- prep$row_presence[, t] &
			prep$row_presence[, t - 1L] &
			active$row[, t] & active$row[, t - 1L]
		eligible_v[, t] <- prep$col_presence[, t] &
			prep$col_presence[, t - 1L] &
			active$col[, t] & active$col[, t - 1L]
	}
	q_u <- matrix(NA_real_, prep$nr, prep$Tt,
	              dimnames = list(prep$row_names, prep$time_names))
	q_v <- matrix(NA_real_, prep$nc, prep$Tt,
	              dimnames = list(prep$col_names, prep$time_names))
	q_u[eligible_u] <- 0
	q_v[eligible_v] <- 0
	rho_u <- rho_v <- if (is.null(rho_uv)) 0.8 else as.numeric(rho_uv[1L])
	if (!is.null(rho_uv) && length(rho_uv) > 1L) rho_v <- as.numeric(rho_uv[2L])
	sigma_u <- sigma_v <- if (is.null(sigma_uv)) 1 else max(min_sigma, as.numeric(sigma_uv[1L]))
	if (!is.null(sigma_uv) && length(sigma_uv) > 1L) sigma_v <- max(min_sigma, as.numeric(sigma_uv[2L]))
	pi_u <- pi_v <- snap_pi_prior[1L] / sum(snap_pi_prior)
	obs_var <- .lame_bsnap_resid_s2(prep$Yarr, Xb, U, V, G, mu, a, b, mask)
	if (is.null(rho_uv) || is.null(sigma_uv)) {
		hu0 <- .lame_bsnap_update_hyper_side(
			U, q_u, eligible_u, estimate_rho_uv, estimate_sigma_uv,
			rho_u, sigma_u, snap_pi_prior, min_sigma,
			rho_prior_mean, rho_prior_weight,
			hyper_update = hyper_update,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions)
		hv0 <- .lame_bsnap_update_hyper_side(
			V, q_v, eligible_v, estimate_rho_uv, estimate_sigma_uv,
			rho_v, sigma_v, snap_pi_prior, min_sigma,
			rho_prior_mean, rho_prior_weight,
			hyper_update = hyper_update,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions)
		if (is.null(rho_uv)) {
			rho_u <- hu0$rho
			rho_v <- hv0$rho
		}
		if (is.null(sigma_uv)) {
			sigma_u <- hu0$sigma
			sigma_v <- hv0$sigma
		}
		pi_u <- hu0$pi_snap
		pi_v <- hv0$pi_snap
	}
	sigma_floor_u <- max(min_sigma, sigma_floor_fraction * sigma_u)
	sigma_floor_v <- max(min_sigma, sigma_floor_fraction * sigma_v)

	objective_trace <- numeric(0)
	param_trace <- data.frame()
	converged <- FALSE
	prev_obj <- Inf
	last_q_u <- q_u
	last_q_v <- q_v
	last_delta_u <- matrix(NA_real_, prep$nr, prep$Tt,
	                       dimnames = list(prep$row_names, prep$time_names))
	last_delta_v <- matrix(NA_real_, prep$nc, prep$Tt,
	                       dimnames = list(prep$col_names, prep$time_names))
	last_delta_summary <- NULL
	if (verbose) {
		cli::cli_h3("Bipartite dynamic snap ALS")
		cli::cli_text("Network: {prep$nr} x {prep$nc}, {prep$Tt} periods, R_row = {R_row}, R_col = {R_col}")
	}
	for (iter in seq_len(max_iter)) {
		stat <- .lame_bsnap_static_update(prep$Yarr, prep$X, U, V, G, mu,
		                                  beta, a, b, mask)
		mu <- stat$mu
		beta <- stat$beta
		a <- stat$a
		b <- stat$b
		Xb <- stat$Xb
		uv <- .lame_bsnap_update_uvg(prep$Yarr, Xb, U, V, G, mu, a, b,
		                             mask, q_u, q_v, rho_u, sigma_u,
		                             rho_v, sigma_v, snap_kappa,
		                             eligible_u, eligible_v, ridge,
		                             obs_var)
		U <- uv$U
		V <- uv$V
		G <- uv$G
		obs_var <- .lame_bsnap_resid_s2(prep$Yarr, Xb, U, V, G, mu, a, b, mask)
		qu <- .lame_bsnap_score_side(U, q_u, eligible_u, rho_u, sigma_u,
		                             snap_kappa, pi_u, snap_update,
		                             snap_damping, threshold, iter, max_iter)
		qv <- .lame_bsnap_score_side(V, q_v, eligible_v, rho_v, sigma_v,
		                             snap_kappa, pi_v, snap_update,
		                             snap_damping, threshold, iter, max_iter)
		q_u <- qu$q
		q_v <- qv$q
		hu <- .lame_bsnap_update_hyper_side(
			U, q_u, eligible_u, estimate_rho_uv, estimate_sigma_uv,
			rho_u, sigma_u, snap_pi_prior, min_sigma,
			rho_prior_mean, rho_prior_weight,
			hyper_update = hyper_update,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions)
		hv <- .lame_bsnap_update_hyper_side(
			V, q_v, eligible_v, estimate_rho_uv, estimate_sigma_uv,
			rho_v, sigma_v, snap_pi_prior, min_sigma,
			rho_prior_mean, rho_prior_weight,
			hyper_update = hyper_update,
			drift_quantile = drift_quantile,
			drift_min_transitions = drift_min_transitions)
		rho_u <- hu$rho
		sigma_u <- max(sigma_floor_u, hu$sigma)
		pi_u <- hu$pi_snap
		rho_v <- hv$rho
		sigma_v <- max(sigma_floor_v, hv$sigma)
		pi_v <- hv$pi_snap
		obj <- .lame_bsnap_objective(
			prep$Yarr, Xb, U, V, G, mu, a, b, mask, q_u, q_v,
			rho_u, sigma_u, rho_v, sigma_v, snap_kappa, pi_u, pi_v,
			eligible_u, eligible_v, obs_var)
		objective_trace <- c(objective_trace, obj)
		last_delta_u <- abs(q_u - last_q_u)
		last_delta_v <- abs(q_v - last_q_v)
		last_delta_summary <- .lame_snap_delta_summary(
			last_delta_u, last_delta_v, q_u, q_v, last_q_u, last_q_v,
			symmetric = FALSE, threshold, snap_stability_tol,
			snap_delta_quantile)
		class_stable <- is.finite(last_delta_summary$share_class_changed) &&
			last_delta_summary$share_class_changed <= snap_class_change_tol
		snap_delta_for_convergence <- switch(
			snap_convergence,
			max = last_delta_summary$max,
			quantile = last_delta_summary$convergence_delta,
			classification = if (isTRUE(class_stable)) 0 else Inf)
		if (!is.finite(snap_delta_for_convergence)) snap_delta_for_convergence <- Inf
		snap_score_stable <- if (identical(snap_convergence, "classification")) {
			class_stable
		} else {
			snap_delta_for_convergence <= snap_stability_tol && class_stable
		}
		param_trace <- rbind(param_trace,
		                     data.frame(iter = iter, objective = obj,
		                                rho_u = rho_u, sigma_u = sigma_u,
		                                pi_u = pi_u,
		                                rho_v = rho_v, sigma_v = sigma_v,
		                                pi_v = pi_v,
		                                obs_var = obs_var,
		                                max_snap_delta = last_delta_summary$max,
		                                snap_delta_q95 = last_delta_summary$q95,
		                                snap_delta_q99 = last_delta_summary$q99,
		                                snap_delta_conv = snap_delta_for_convergence,
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

	fitbits <- .lame_bsnap_fitted_lists(prep, Xb, U, V, G, mu, a, b, mask)
	mu_id <- .lame_snap_identified_mu(
		mu, .lame_bsnap_uvg_cube(U, V, G), prep$time_names,
		exclude_diagonal = FALSE)
	names(beta) <- prep$x_names_dyad
	names(a) <- prep$row_names
	names(b) <- prep$col_names
	coefficients <- c(intercept = mu, beta)
	snap_class <- q_u > threshold
	snap_class[, 1L] <- NA
	snap_class_v <- q_v > threshold
	snap_class_v[, 1L] <- NA
	unstable_transitions <- .lame_snap_top_unstable(
		last_delta_u, last_delta_v, q_u, q_v, last_q_u, last_q_v,
		symmetric = FALSE, threshold, unstable_top_n)
	if (is.null(last_delta_summary)) {
		last_delta_summary <- .lame_snap_delta_summary(
			last_delta_u, last_delta_v, q_u, q_v, last_q_u, last_q_v,
			FALSE, threshold, snap_stability_tol, snap_delta_quantile)
	}
	final_snap_delta_conv <- if (nrow(param_trace) > 0L) {
		utils::tail(param_trace$snap_delta_conv, 1L)
	} else NA_real_
	final_share_class_changed <- if (nrow(param_trace) > 0L) {
		utils::tail(param_trace$share_snap_class_changed, 1L)
	} else NA_real_
	class_stable <- is.finite(final_share_class_changed) &&
		final_share_class_changed <= snap_class_change_tol
	snap_stable <- if (identical(snap_convergence, "classification")) {
		class_stable
	} else {
		is.finite(final_snap_delta_conv) &&
			final_snap_delta_conv <= snap_stability_tol && class_stable
	}
	objective_diff <- diff(objective_trace)
	objective_increase_tol <- tol * (abs(utils::head(objective_trace, -1L)) + 1)
	objective_increases <- if (length(objective_diff) > 0L) {
		sum(objective_diff > objective_increase_tol, na.rm = TRUE)
	} else 0L
	surface <- .lame_bsnap_surface(U, V, G, mask, rho_u, rho_v,
	                               eligible_u, eligible_v)
	g_diag <- .lame_bsnap_g_diagnostics(G)
	if (verbose && isTRUE(g_diag$near_singular)) {
		cli::cli_warn(c(
			"Bipartite snap ALS fitted a near-singular interaction matrix.",
			"i" = "Coordinate snap scores in weakly identified dimensions can be unstable; inspect {.code fit$G_diagnostics}."))
	}
	out <- list(
		call = match.call(),
		method = "als_snap",
		family = "normal",
		mode = "bipartite",
		symmetric = FALSE,
		R = max(R_row, R_col),
		R_row = R_row,
		R_col = R_col,
		longitudinal = TRUE,
		mu = mu,
		mu_t = mu_id$mu_t,
		mu_identified = mu_id$mu_identified,
		beta = beta,
		a = a,
		b = b,
		U = U,
		V = V,
		L = NULL,
		G = G,
		G_diagnostics = g_diag,
		coefficients = coefficients,
		VC = c(va = stats::var(a), cab = NA_real_,
		       vb = stats::var(b), rho = 0, ve = fitbits$s2),
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
		convergence = list(
			converged = converged,
			iterations = nrow(param_trace),
			tolerance = tol,
			snap_stability_tol = snap_stability_tol,
			snap_convergence = snap_convergence,
			snap_delta_quantile = snap_delta_quantile,
			snap_class_change_tol = snap_class_change_tol,
			final_max_snap_delta = last_delta_summary$max,
			final_snap_delta_q95 = last_delta_summary$q95,
			final_snap_delta_q99 = last_delta_summary$q99,
			final_snap_delta_conv = final_snap_delta_conv,
			final_n_snap_delta_gt_tol = last_delta_summary$n_unstable,
			final_share_snap_delta_gt_tol = last_delta_summary$share_unstable,
			final_n_snap_class_changed = last_delta_summary$n_class_changed,
			final_share_snap_class_changed = final_share_class_changed,
			snap_stable = snap_stable,
			max_snap_stable = is.finite(last_delta_summary$max) &&
				last_delta_summary$max <= snap_stability_tol,
			snap_class_stable = class_stable,
			unstable_transitions = unstable_transitions,
			objective_increases = objective_increases,
			max_objective_increase = if (length(objective_diff) > 0L)
				max(c(0, objective_diff), na.rm = TRUE) else 0,
			final_objective = utils::tail(objective_trace, 1)),
		rho_uv = c(row = rho_u, col = rho_v),
		sigma_uv = c(row = sigma_u, col = sigma_v),
		obs_var = fitbits$s2,
		pi_snap = c(row = pi_u, col = pi_v),
		snap_prob = q_u,
		snap_prob_v = q_v,
		snap_class = snap_class,
		snap_class_v = snap_class_v,
		surface_jump = surface$row,
		surface_jump_v = surface$col,
		transition_diagnostics = list(
			row_presence = prep$row_presence,
			col_presence = prep$col_presence,
			eligible = eligible_u,
			eligible_v = eligible_v,
			row_exposure = active$row_exposure,
			col_exposure = active$col_exposure,
			surface_jump = surface$row,
			surface_jump_v = surface$col,
			drift_dist2 = qu$drift_d2,
			snap_dist2 = qu$snap_d2,
			logit_score = qu$logit,
			drift_dist2_v = qv$drift_d2,
			snap_dist2_v = qv$snap_d2,
			logit_score_v = qv$logit,
			last_snap_delta = last_delta_u,
			last_snap_delta_v = last_delta_v,
			unstable_transitions = unstable_transitions),
		n_time = prep$Tt,
			dims = list(n_row = prep$nr, n_col = prep$nc, p = length(beta)),
			row_names = prep$row_names,
			col_names = prep$col_names,
			x_names = names(beta),
			time_names = prep$time_names,
		row_presence = prep$row_presence,
		col_presence = prep$col_presence,
		actor_presence = NULL,
		changing_composition = isTRUE(prep$changing_composition),
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
		               mode = "bipartite"),
		meta = list(
			sampler = "als_snap",
				algorithm = "bipartite profiled block coordinate descent / EM snap-shift ALS",
				align = align,
				snap_kappa = snap_kappa,
				snap_pi_prior = snap_pi_prior,
				hyper_update = hyper_update,
				drift_quantile = drift_quantile,
				drift_min_transitions = drift_min_transitions,
				rho_prior_mean = rho_prior_mean,
				rho_prior_weight = rho_prior_weight,
				sigma_floor_fraction = sigma_floor_fraction,
				snap_convergence = snap_convergence,
			snap_delta_quantile = snap_delta_quantile,
			snap_class_change_tol = snap_class_change_tol,
			uncertainty_available = FALSE,
			approximation_note = paste0(
				"Fast bipartite dynamic snap-shift ALS fit with static G. ",
				"Snap values are ALS segmentation scores, not posterior draws.")))
	class(out) <- c("lame_snap_als", "lame_als", "ame_als")
	if (verbose && !isTRUE(out$convergence$snap_stable)) {
		cli::cli_warn(c(
			"{.fn lame_snap_als} stopped while bipartite snap scores were still moving.",
			"i" = "Final {out$convergence$snap_convergence} delta: {round(out$convergence$final_snap_delta_conv, 4)}; tolerance: {out$convergence$snap_stability_tol}.",
			"i" = "See {.code fit$convergence$unstable_transitions} for the largest moves."))
	}
	out
}
