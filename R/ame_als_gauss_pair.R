# dyadic correlation of the normal empirical-bayes fit of a directed network (dcor = TRUE)
# the pair errors are bivariate normal with variance s2 and correlation rho, so every pair quantity is closed form

# posterior mean and sd of rho from the pair statistics A = sum E[e_ij^2 + e_ji^2], B = sum E[e_ij e_ji]
# and their count P, on a grid of rho = sin(theta), refined on a rho grid when sd < 0.05
.ae_gauss_rho_post <- function(A, B, P, ve, power = .lame_rho_prior_power(), ngrid = 1201L) {
	# log density of rho up to a constant, and its grid moments
	lden <- function(r) {
		q <- 1 - r^2
		-0.5 * P * log(q) - (A - 2 * r * B) / (2 * ve * q) + .lame_rho_prior_log(r, power)
	}
	moments <- function(r, lw) {
		w <- exp(lw - max(lw))
		m <- sum(w * r) / sum(w)
		list(mean = m, sd = sqrt(sum(w * (r - m)^2) / sum(w)))
	}
	th <- (seq_len(ngrid) - 0.5) / ngrid * pi - pi / 2
	r <- sin(th)
	out <- moments(r, lden(r) + log(cos(th)))
	# a posterior narrower than the grid cells is integrated again over 12 sd about the mean
	if (out$sd < 0.05) {
		lo <- max(out$mean - 12 * out$sd, -0.99999); hi <- min(out$mean + 12 * out$sd, 0.99999)
		rg <- seq(lo, hi, length.out = ngrid)
		out <- moments(rg, lden(rg))
	}
	out
}

# start of rho: the residual correlation of the pairs, kept inside (-0.9, 0.9), or the caller's value
.ae_gauss_rho_start <- function(cx, S) {
	if (!is.null(cx$pair_rho0) && is.finite(cx$pair_rho0)) return(max(min(cx$pair_rho0, 0.9), -0.9))
	pc <- cx$pair
	if (pc$n_pairs < 3L) return(0)
	eta <- .ae_mf_predictor(cx, S)$eta
	r1 <- pc$y1 - eta[pc$one]; r2 <- pc$y2 - eta[pc$two]
	est <- suppressWarnings(stats::cor(r1, r2))
	if (!is.finite(est)) return(0)
	max(min(est, 0.9), -0.9)
}

# expected working problem of every cell when dyads are pairs: the .ae_mf_cells arrays (pair log likelihood
# and whitened second moment in the first cell), the cross weight cw, G = W y, and the statistics, score and curvature of rho
.ae_gauss_pair_cells <- function(cx, S, eta, V3, rho) {
	pc <- cx$pair
	ve <- S$ve
	yy <- cx$y
	yy[c(pc$one, pc$two)] <- NA_real_
	wk <- .ae_gauss_cells(yy, eta, V3, ve)
	m1 <- eta[pc$one]; m2 <- eta[pc$two]
	v1 <- V3[pc$one]; v2 <- V3[pc$two]
	cvv <- .ae_pair_cov(cx, S)
	r1 <- pc$y1 - m1; r2 <- pc$y2 - m2
	A <- sum(r1^2 + r2^2 + v1 + v2)
	B <- sum(r1 * r2 + cvv)
	P <- pc$n_pairs
	q <- 1 - rho^2
	k <- 1 / (ve * q)
	z <- wk$z; w <- wk$w; ell <- wk$ell; ee <- wk$ee
	cw <- G <- array(0, dim(cx$y))
	# a pair cell's response is its observation less rho times the other cell's residual
	z[pc$one] <- pc$y1 - rho * r2
	z[pc$two] <- pc$y2 - rho * r1
	w[pc$one] <- k; w[pc$two] <- k
	cw[pc$one] <- -rho * k; cw[pc$two] <- -rho * k
	G[pc$one] <- k * (pc$y1 - rho * pc$y2)
	G[pc$two] <- k * (pc$y2 - rho * pc$y1)
	rest <- which(is.finite(yy))
	G[rest] <- wk$w[rest] * wk$z[rest]
	# log likelihood of a pair and its whitened second moment, in the first cell
	ee_pair <- ((r1^2 + r2^2 + v1 + v2) - 2 * rho * (r1 * r2 + cvv)) / q
	ell[pc$one] <- -log(2 * pi * ve) - 0.5 * log(q) - 0.5 * k * ((r1^2 + r2^2 + v1 + v2) - 2 * rho * (r1 * r2 + cvv))
	ell[pc$two] <- 0
	ee[pc$one] <- ee_pair; ee[pc$two] <- 0
	# score and curvature of rho
	N <- B * (1 + rho^2) - rho * A
	rs <- P * rho / q + N / (ve * q^2)
	rc <- P * (1 + rho^2) / q^2 + ((2 * rho * B - A) * q + 4 * rho * N) / (ve * q^3)
	post <- if (P >= 3L) .ae_gauss_rho_post(A, B, P, ve) else list(mean = 0, sd = NA_real_)
	list(
		wk = list(z = z, w = w, ell = ell, ee = ee),
		pair = list(
			cw = cw,
			G = G,
			rs = rs,
			rc = rc,
			rho = rho,
			cv = cvv,
			A = A,
			B = B,
			rho_next = post$mean,
			rho_sd = post$sd
		)
	)
}

# cavity weights of pair cells: the information of y ~ N(m, s2 [[1, rho], [rho, 1]] + [[v1, cv], [cv, v2]])
# returns the diagonal weights (single-cell weight when the partner is missing) and the cross weights
.ae_gauss_pair_info_w <- function(cx, eta, V3, cv, rho, ve) {
	pc <- cx$pair
	yy <- cx$y
	yy[c(pc$one, pc$two)] <- NA_real_
	w <- .ae_gauss_info_w(yy, eta, V3, ve)
	cw <- array(0, dim(cx$y))
	c11 <- ve + V3[pc$one]; c22 <- ve + V3[pc$two]; c12 <- rho * ve + cv
	det <- c11 * c22 - c12^2
	w[pc$one] <- c22 / det
	w[pc$two] <- c11 / det
	cw[pc$one] <- -c12 / det; cw[pc$two] <- -c12 / det
	list(w = w, cw = cw)
}

# refuse a directed normal network whose reciprocal pairs are all exactly equal, as the sampler does
# (.ame_check_normal_reciprocity): the posterior of rho is improper once the pairs exceed the mean rank bound
.ae_gauss_check_reciprocity <- function(prep, R) {
	Y <- prep$Yarr
	n <- prep$nr
	ut <- upper.tri(matrix(0, n, n))
	pairs <- NULL
	all_equal <- TRUE
	for (tt in seq_len(prep$Tt)) {
		y <- Y[, , tt]
		pr <- which(ut & !is.na(y) & !is.na(t(y)), arr.ind = TRUE)
		if (!nrow(pr)) next
		all_equal <- all_equal && all(y[pr] == t(y)[pr])
		pairs <- rbind(pairs, cbind(pr, period = tt))
		if (!all_equal) return(invisible(NULL))
	}
	m <- if (is.null(pairs)) 0L else nrow(pairs)
	if (m == 0L || !all_equal) return(invisible(NULL))
	# rank absorbed by the mean (.lame_panel_reciprocity_bound); the static effects, coefficients and factors count once
	Xlist <- if (prep$p_dyad > 0L) lapply(
		seq_len(prep$Tt),
		function(tt) array(prep$X[, , , tt], c(n, n, prep$p_dyad))
	) else
		lapply(seq_len(prep$Tt), function(tt) array(0, c(n, n, 0L)))
	bound <- .lame_panel_reciprocity_bound(Y, Xlist, TRUE, TRUE, R)$rank
	excess <- .lame_rho_improper_dyads()
	if (m - bound >= excess) {
		cli::cli_abort(
			c(
				"This directed normal model has an improper posterior under the reciprocity prior.",
				"i" = "All {m} observed reciprocal pairs are exactly equal; the actor, regression and factor rank bound is {bound}, and the posterior is improper once the pairs exceed it by {excess}.",
				"i" = "Use {.code symmetric = TRUE} for an undirected network, or {.code dcor = FALSE} for fixed zero reciprocity."
			),
			class = "lame_improper_normal_reciprocity"
		)
	}
	invisible(NULL)
}
