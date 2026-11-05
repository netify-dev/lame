# stationarity check of the observed likelihood plus the path-solve penalties;
# additive scores are projected onto their per-period sum-to-zero constraints
.lda_stationarity <- function(prep, state, settings, tol = 1e-6) {
	s <- settings
	v <- state
	nr <- prep$nr
	nc <- prep$nc
	nt <- prep$Tt
	eta <- .lda_eta_array(.lda_eta(
		prep,
		v$coef_path,
		v$a_mat,
		v$b_mat,
		v$Oarr,
		symmetric = s$symmetric,
		beta_dyn = s$beta_dyn
	))
	y <- prep$Yarr
	observed <- is.finite(y)
	if (s$family == "normal") {
		score <- y - eta
		info <- array(1, dim(y))
	} else if (s$family == "poisson") {
		mu <- exp(eta)
		score <- y - mu
		info <- mu
	} else if (s$link == "logit") {
		p <- stats::plogis(eta)
		q <- stats::plogis(-eta)
		score <- y * q - (1 - y) * p
		info <- p * q
	} else {
		sign_y <- 2 * y - 1
		score <- sign_y * exp(dnorm(eta, log = TRUE) - pnorm(sign_y * eta, log.p = TRUE))
		info <- exp(
			2 * dnorm(eta, log = TRUE) - pnorm(eta, log.p = TRUE) -
				pnorm(eta, lower.tail = FALSE, log.p = TRUE)
		)
	}
	score[!observed] <- 0
	info[!observed] <- 0
	gradients <- curvatures <- list()
	add <- function(name, gradient, information) {
		gradients[[name]] <<- 2 * gradient
		curvatures[[name]] <<- 2 * information
	}
	p <- nrow(v$coef_path)
	bg <- bi <- matrix(0, p, nt)
	for (tt in seq_len(nt)) if (any(observed[, , tt])) {
		x <- .lda_design_obs(
			prep,
			tt,
			observed[, , tt],
			rownames(v$coef_path),
			symmetric = s$symmetric
		)
		bg[, tt] <- -crossprod(x, score[, , tt][observed[, , tt]])
		bi[, tt] <- colSums(x^2 * info[, , tt][observed[, , tt]])
	}
	for (k in seq_len(p)) {
		if (s$beta_dyn$mask[k]) {
			pen <- .lda_path_penalty(nt, s$lambda_beta * .lda_beta_ms(s$beta_dyn, k), s$dynamic_beta_kind, s$rho_beta)
			add(paste0("beta", k), bg[k, ] + pen %*% v$coef_path[k, ], bi[k, ] + diag(pen))
		} else add(paste0("beta", k), sum(bg[k, ]), sum(bi[k, ]))
	}
	ag <- -apply(score, c(1, 3), sum)
	ai <- apply(info, c(1, 3), sum)
	bg <- -apply(score, c(2, 3), sum)
	bi <- apply(info, c(2, 3), sum)
	if (s$symmetric) {
		ag <- ag + bg
		ai <- ai + bi
	}
	for (side in if (s$symmetric) "a" else c("a", "b")) {
		g <- if (side == "a") ag else bg
		h <- if (side == "a") ai else bi
		x <- if (side == "a") v$a_mat else v$b_mat
		if (s$dynamic_ab) {
			# actors never observed and without a level prior stay at zero and leave the projection
			act <- if (isTRUE(s$level_prec_ab > 0)) rep(TRUE, nrow(h)) else rowSums(h) > 0
			tw <- .lda_transition_from_presence(if (side == "a") prep$row_presence else prep$col_presence)
			for (i in seq_len(nrow(x))) {
				pen <- .lda_path_penalty_weighted(
					nt,
					s$lambda_ab,
					if (is.null(tw)) NULL else tw[i, ],
					s$rho_ab
				) + diag(s$level_prec_ab, nt)
				g[i, ] <- g[i, ] + pen %*% x[i, ]
				h[i, ] <- h[i, ] + diag(pen)
			}
			# project the score over the actors present in each period, matching the block solves
			pres <- if (side == "a") prep$row_presence else prep$col_presence
			if ((is.null(pres) || all(pres)) && all(act)) {
				g <- sweep(g, 2, colMeans(g), "-")
			} else {
				for (tt in seq_len(nt)) {
					ok <- act & (if (is.null(pres)) TRUE else pres[, tt])
					g[!ok, tt] <- 0
					if (any(ok)) g[ok, tt] <- g[ok, tt] - mean(g[ok, tt])
				}
			}
		} else {
			g <- rowSums(g)
			h <- rowSums(h)
			# actors with no observed cell have no score and stay out of the sum-to-zero projection
			act <- h > 0
			g[!act] <- 0
			g[act] <- g[act] - mean(g[act])
		}
		add(side, g, h)
	}
	if (!is.null(v$U) && length(v$U)) {
		U <- v$U
		V <- v$V
		dynamic <- length(dim(U)) == 3L
		ru <- dim(U)[2L]
		rv <- dim(V)[2L]
		if (s$symmetric && !dynamic && !is.null(v$L)) {
			U <- sweep(U, 2, sqrt(abs(v$L)), "*")
			V <- sweep(U, 2, sign(v$L), "*")
		}
		ug <- ui <- array(0, c(nr, ru, nt))
		vg <- vi <- array(0, c(nc, rv, nt))
		gg <- gi <- if (prep$bip) matrix(0, ru, rv) else NULL
		for (tt in seq_len(nt)) {
			ut <- if (dynamic) matrix(U[, , tt], nr, ru) else U
			vt <- if (dynamic) matrix(V[, , tt], nc, rv) else V
			ue <- if (prep$bip) ut %*% v$G_static else ut
			ve <- if (prep$bip) vt %*% t(v$G_static) else vt
			ug[, , tt] <- -score[, , tt] %*% ve
			ui[, , tt] <- info[, , tt] %*% (ve^2)
			vg[, , tt] <- -t(score[, , tt]) %*% ue
			vi[, , tt] <- t(info[, , tt]) %*% (ue^2)
			if (prep$bip) {
				gg <- gg - crossprod(ut, score[, , tt]) %*% vt
				gi <- gi + crossprod(ut^2, info[, , tt]) %*% (vt^2)
			}
		}
		if (s$symmetric) {
			sgn <- sign(if (dynamic) apply(U * V, 2, sum) else colSums(U * V))
			sgn[sgn == 0] <- 1
			ug <- ug + sweep(vg, 2, sgn, "*")
			ui <- ui + vi
		}
		for (side in if (s$symmetric) "u" else c("u", "v")) {
			cube <- if (side == "u") U else V
			g <- if (side == "u") ug else vg
			h <- if (side == "u") ui else vi
			if (dynamic) {
				tw <- .lda_transition_from_presence(if (side == "u") prep$row_presence else prep$col_presence)
				if (identical(s$dynamic_uv_kind, "t")) {
					tw <- .lda_combine_transition_weights(
						.lda_t_transition_weights(
							cube,
							s$rho_uv,
							s$lambda_uv
						),
						tw
					)
				}
				for (i in seq_len(dim(cube)[1L])) {
					pen <- .lda_path_penalty_weighted(
						nt,
						s$lambda_uv,
						if (is.null(tw)) NULL else tw[i, ],
						s$rho_uv
					) + diag(s$level_prec_uv, nt)
					if (s$symmetric) pen <- 2 * pen
					for (k in seq_len(dim(cube)[2L])) {
						g[i, k, ] <- g[i, k, ] + pen %*% cube[i, k, ]
						h[i, k, ] <- h[i, k, ] + diag(pen)
					}
				}
			} else {
				# the symmetric eigenmodel uses one factor on both sides, so its ridge is doubled
				m <- if (s$symmetric) 2 else 1
				g <- apply(g, c(1, 2), sum) + m * s$level_prec_uv * cube
				h <- apply(h, c(1, 2), sum) + m * s$level_prec_uv
			}
			add(side, g, h)
		}
		if (prep$bip) add("G", gg + s$level_prec_g * v$G_static, gi + s$level_prec_g)
	}
	# scores are divided by the response unit, so the threshold is relative for normal data
	unit <- if (is.null(s$unit)) 1 else s$unit
	scaled <- Map(function(g, h) abs(g) / sqrt(pmax(h, .Machine$double.xmin)) / unit, gradients, curvatures)
	# an empty block (a rank-zero side in a bipartite fit) has no score
	maxima <- vapply(scaled, function(z) if (length(z)) max(z) else 0, numeric(1))
	maximum <- max(maxima)
	list(
		pass = all(is.finite(maxima)) && maximum <= sqrt(tol),
		max_score = maximum,
		threshold = sqrt(tol),
		by_block = maxima,
		gradient = gradients,
		information = curvatures
	)
}
