# centre the reported effects of the dynamic EB fit over the actors, as the
# samplers do: only mu + mean(a) + mean(b) is identified, so the report moves
# a -> a - mean(a), b -> b - mean(b), mu -> mu + mean(a) + mean(b), and the
# coefficient covariance takes that of this functional of the joint posterior
.lde_centre_report <- function(ctx, ev, ab, lv, th, cpath, vc, vcc, st = NULL) {
	lay <- ctx$beta
	na <- ctx$na; nb <- ctx$nb; Tt <- ctx$Tt; q <- ctx$q
	keep <- list(
		ab = ab,
		lv = lv,
		th = th,
		cpath = cpath,
		vc = vc,
		vcc = vcc,
		centre = c(a = 0, b = 0, o = 0),
		o_shift = 0
	)
	fzm <- ctx$fzm
	on_a <- !isTRUE(fzm[1L]); on_b <- !ctx$sym && !isTRUE(fzm[2L])
	abar <- if (on_a) mean(ab$a[ctx$ia, , drop = FALSE]) else 0
	bbar <- if (ctx$sym) abar else if (on_b) mean(ab$b[ctx$ib, , drop = FALSE]) else 0
	c0 <- abar + bbar
	if (!is.finite(c0) || is.null(ev$fac)) return(keep)
	# the functional: the intercept plus the mean of the paths of the enabled margins
	A <- .lde_avg_map(ctx)
	ct <- A[1L, ]
	wa <- if (on_a) 1 / (na * Tt) else 0
	wb <- if (ctx$sym) wa else if (on_b) 1 / (nb * Tt) else 0
	cs <- numeric(q)
	cs[seq_len(na)] <- if (ctx$sym) 2 * wa else wa
	if (!ctx$sym) cs[na + seq_len(nb)] <- wb
	rs <- rep(list(cs), Tt)
	rz <- c(numeric(q), ct)
	w <- .lde_solve(ev$fac, rs, rz)
	var_f <- sum(vapply(seq_len(Tt), function(t) sum(cs * w$s[[t]]), 0)) + sum(rz * w$z)
	it <- q + seq_len(ctx$mth)
	Szt <- ev$fac$Szi[it, it, drop = FALSE]
	d_var <- var_f - as.numeric(ct %*% Szt %*% ct)
	d_cov <- as.vector(A %*% (w$z[it] - Szt %*% ct))
	shift_cov <- function(V) {
		V[1L, 1L] <- V[1L, 1L] + d_var
		V[1L, -1L] <- V[1L, -1L] + d_cov[-1L]
		V[-1L, 1L] <- V[-1L, 1L] + d_cov[-1L]
		V
	}
	if (!is.finite(d_var) || !all(is.finite(d_cov)) || shift_cov(vc)[1L, 1L] <= 0) return(keep)
	# the common shifts of the positions are integrated out as well
	extra <- NULL
	o_shift <- 0
	if (!is.null(st) && identical(ctx$fmode, "dynamic")) {
		p <- nrow(A)
		cs_l <- lapply(seq_len(p), function(j) if (j == 1L) cs else numeric(q))
		rz_l <- lapply(seq_len(p), function(j) c(numeric(q), A[j, ]))
		col <- .lde_pos_collapse(ctx, ev, st, cs_l, rz_l)
		if (!is.null(col) && all(is.finite(col$extra)) && all(is.finite(col$shift))) {
			extra <- col$extra
			# undirected networks: the mean shift moves into the intercept, since the
			# common shift always adds to the mean of the product
			if (isTRUE(ctx$sym)) o_shift <- col$shift[1L]
		}
	}
	shift_cov0 <- shift_cov
	if (!is.null(extra)) shift_cov <- function(V) shift_cov0(V) + extra
	# moving intercept: the same exchange in every period, each period's standard
	# error from its own combination with the mean of the effects
	d_var_path <- NULL
	if (isTRUE(lay$any) && !(1L %in% lay$static_idx)) {
		m <- ctx$mth
		M <- vapply(seq_len(m), function(k) { e <- numeric(m); e[k] <- 1; .lde_coef_path(ctx, e)[1L, ] }, numeric(Tt))
		d_var_path <- vapply(seq_len(Tt), function(t) {
			ctt <- M[t, ]
			wt <- .lde_solve(ev$fac, rs, c(numeric(q), ctt))
			vt <- sum(vapply(seq_len(Tt), function(u) sum(cs * wt$s[[u]]), 0)) + sum(c(numeric(q), ctt) * wt$z)
			vt - as.numeric(ctt %*% Szt %*% ctt)
		}, 0)
		if (!all(is.finite(d_var_path))) d_var_path <- NULL
	}
	ab$a[ctx$ia, ] <- ab$a[ctx$ia, , drop = FALSE] - abar
	ab$b[ctx$ib, ] <- ab$b[ctx$ib, , drop = FALSE] - bbar
	lv$a[ctx$ia] <- lv$a[ctx$ia] - abar
	lv$b[ctx$ib] <- lv$b[ctx$ib] - bbar
	# a bipartite interaction with static factors is held at mean zero, as in the
	# samplers, its mean moving into the intercept
	if (identical(ctx$fmode, "fac") && isTRUE(ctx$bip) && !is.null(st$O)) o_shift <- o_shift + mean(st$O)
	th[1L] <- th[1L] + c0 + o_shift
	cpath[1L, ] <- cpath[1L, ] + c0 + o_shift
	vc_new <- shift_cov(vc)
	vcc_new <- shift_cov(vcc)
	# static factors: the intercept takes its variance conditional on the factors
	# (the marginal one is too wide), with its covariances rescaled to keep the
	# correlations
	if (identical(ctx$fmode, "fac") && isTRUE(ctx$nrm) && vc_new[1L, 1L] > 0 && vcc_new[1L, 1L] > 0) {
		k <- sqrt(vcc_new[1L, 1L] / vc_new[1L, 1L])
		vc_new[1L, ] <- vc_new[1L, ] * k
		vc_new[, 1L] <- vc_new[, 1L] * k
		vc_new[1L, 1L] <- vcc_new[1L, 1L]
	}
	list(
		ab = ab,
		lv = lv,
		th = th,
		cpath = cpath,
		vc = vc_new,
		vcc = vcc_new,
		centre = c(a = abar, b = bbar, o = o_shift),
		d_var_path = d_var_path,
		o_shift = o_shift
	)
}

# Laplace collapse of the common shifts of the moving positions, and its effect
# on the reported coefficients. the mean-field posterior misses the direction in
# which all positions of a side shift together in a period; the shifts z are
# added to the unknowns, linearised at the mean-field means. with H the
# information of the additive unknowns and coefficients, B its coupling with z
# and D the information of z, S = D - B' H^-1 B, and the covariance of a
# functional f grows by (G + h' grad delta)' S^-1 (G + h' grad delta) with
# G = B' H^-1 f. for a Gaussian response the profile mode z* = -S^-1 g shifts
# the predictor constant by delta(z) (c d directed, signed c^2 undirected),
# which the intercept absorbs. returns NULL without AR(1) positions
.lde_pos_collapse <- function(ctx, ev, st, cs_list, rz_list, force = TRUE, details = FALSE) {
	pos <- st$pos
	if (is.null(pos) || !identical(ctx$fmode, "dynamic") || !identical(ctx$uv_kind, "ar1") || is.null(ev$fac) || is.null(ev$wk$w)) return(NULL)
	Tt <- ctx$Tt; q <- ctx$q; nu <- ctx$nu; sym <- isTRUE(ctx$sym)
	W <- ev$wk$w
	CW <- ev$wk$cw
	qu <- pos$qu; qv <- if (sym) pos$qu else pos$qv
	Ru <- dim(qu$mean)[2L]; Rv <- dim(qv$mean)[2L]
	rmin <- min(Ru, Rv)
	sg <- if (sym) ctx$sgn else rep(1, rmin)
	# the means and covariances of the positions of period t
	tt <- function(qs, t) min(t, dim(qs$mean)[3L])
	pm <- function(qs, t) matrix(qs$mean[, seq_len(rmin), tt(qs, t)], dim(qs$mean)[1L], rmin)
	psg <- function(qs, t) qs$sig[seq_len(rmin), seq_len(rmin), , tt(qs, t), drop = FALSE]
	sides <- if (sym) "u" else c("u", "v")
	coords <- do.call(
		rbind,
		lapply(sides, function(sd) expand.grid(side = sd, t = seq_len(Tt), r = seq_len(rmin), stringsAsFactors = FALSE))
	)
	Kp <- nrow(coords)
	zc <- coords
	M <- diag(Kp)
	K <- nrow(zc)
	jac <- function(k) {
		t <- coords$t[k]; r <- coords$r[k]
		if (sym) {
			u <- pm(qu, t)[, r]
			return(sg[r] * outer(u, u, "+"))
		}
		if (coords$side[k] == "u") matrix(rep(pm(qv, t)[, r], each = ctx$na), ctx$na, ctx$nb)
		else matrix(rep(pm(qu, t)[, r], times = ctx$nb), ctx$na, ctx$nb)
	}
	zero_rs <- rep(list(numeric(q)), Tt)
	Bp_rs <- vector("list", Kp); Bp_rz <- vector("list", Kp); Jk <- vector("list", Kp)
	for (k in seq_len(Kp)) {
		t <- coords$t[k]
		J <- jac(k)
		Jk[[k]] <- J
		WJ <- W[, , t] * J
		rs <- zero_rs
		rs[[t]] <- .lde_stack_sums(WJ, ctx)
		rz <- numeric(q + ctx$mth)
		rz[q + ctx$thidx[[t]]] <- nu * as.vector(crossprod(ctx$Xc[[t]], as.vector(WJ)))
		if (!is.null(CW)) {
			# the dyadic correlation: the cell and its partner share the information C
			JC <- J * CW[, , t]
			rs[[t]] <- rs[[t]] + c(colSums(JC), rowSums(JC))
			Xt <- ctx$Xc[[t]]
			rz[q + ctx$thidx[[t]]] <- rz[q + ctx$thidx[[t]]] +
				nu * vapply(seq_len(ncol(Xt)), function(l) sum(JC * t(matrix(Xt[, l], ctx$na, ctx$nb))), 0)
		}
		Bp_rs[[k]] <- rs; Bp_rz[[k]] <- rz
	}
	# the coupling of a shift (the tie of a coordinate to itself)
	B_rs <- vector("list", K); B_rz <- vector("list", K)
	for (k in seq_len(K)) {
		ix <- which(M[, k] == 1)
		B_rs[[k]] <- lapply(seq_len(Tt), function(t) Reduce(`+`, lapply(ix, function(i) Bp_rs[[i]][[t]])))
		B_rz[[k]] <- Reduce(`+`, lapply(ix, function(i) Bp_rz[[i]]))
	}
	# the coupling of a constant added to every cell of period t with the unknowns
	R_rs <- R_rz <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		rs <- zero_rs
		rs[[t]] <- .lde_stack_sums(W[, , t], ctx)
		rz <- numeric(q + ctx$mth)
		rz[q + ctx$thidx[[t]]] <- nu * as.vector(crossprod(ctx$Xc[[t]], as.vector(W[, , t])))
		if (!is.null(CW)) {
			rs[[t]] <- rs[[t]] + c(colSums(CW[, , t]), rowSums(CW[, , t]))
			Xt <- ctx$Xc[[t]]
			rz[q + ctx$thidx[[t]]] <- rz[q + ctx$thidx[[t]]] +
				nu * vapply(seq_len(ncol(Xt)), function(l) sum(CW[, , t] * t(matrix(Xt[, l], ctx$na, ctx$nb))), 0)
		}
		R_rs[[t]] <- rs; R_rz[[t]] <- rz
	}
	Dp <- matrix(0, Kp, Kp)
	for (k in seq_len(Kp)) for (l in k:Kp) if (coords$t[k] == coords$t[l]) {
		Dp[k, l] <- Dp[l, k] <- nu * (sum(W[, , coords$t[k]] * Jk[[k]] * Jk[[l]]) +
				(if (is.null(CW)) 0 else sum(CW[, , coords$t[k]] * Jk[[k]] * t(Jk[[l]]))))
	}
	D <- crossprod(M, Dp %*% M)
	# the prior information of the shifts: the n actors' autoregressive paths
	for (sd in sides) {
		ix <- which(zc$side == sd)
		P0 <- (if (sd == "u") qu else qv)$P0
		Rs <- (if (sd == "u") Ru else Rv)
		n <- (if (sd == "u") ctx$na else ctx$nb)
		src <- (zc$t[ix] - 1L) * Rs + zc$r[ix]
		D[ix, ix] <- D[ix, ix] + n * P0[src, src, drop = FALSE]
	}
	dot <- function(rs, rz, sol) sum(vapply(seq_len(Tt), function(t) sum(rs[[t]] * sol$s[[t]]), 0)) + sum(rz * sol$z)
	sols <- lapply(seq_len(K), function(k) .lde_solve(ev$fac, B_rs[[k]], B_rz[[k]]))
	BHB <- matrix(0, K, K)
	for (k in seq_len(K)) for (l in seq_len(K)) BHB[l, k] <- dot(B_rs[[l]], B_rz[[l]], sols[[k]])
	S <- (D - BHB + t(D - BHB)) / 2
	Sinv <- tryCatch(solve(S), error = function(e) NULL)
	if (is.null(Sinv) || any(!is.finite(Sinv))) return(NULL)
	# the force of the variance of the other side on the means (a Gaussian response), and the profile mode of the shifts
	gp <- numeric(Kp)
	if (isTRUE(force) && isTRUE(ctx$nrm)) {
		for (k in seq_len(Kp)) {
			t <- coords$t[k]; r <- coords$r[k]
			if (sym) {
				Mt <- pm(qu, t); Wm <- crossprod(W[, , t], Mt); sgr <- sg[seq_len(rmin)]; So <- psg(qu, t)
				val <- 0
				for (j in seq_len(ctx$na)) val <- val + sum(matrix(So[, , j, 1L], rmin, rmin)[r, ] * sgr * Wm[j, ])
				gp[k] <- -2 * nu * sg[r] * val
			} else if (coords$side[k] == "u") {
				Wm <- crossprod(W[, , t], pm(qu, t)); So <- psg(qv, t)
				val <- 0
				for (j in seq_len(ctx$nb)) val <- val + sum(matrix(So[, , j, 1L], rmin, rmin)[r, ] * Wm[j, ])
				gp[k] <- -val
			} else {
				Wm <- W[, , t] %*% pm(qv, t); So <- psg(qu, t)
				val <- 0
				for (i in seq_len(ctx$na)) val <- val + sum(matrix(So[, , i, 1L], rmin, rmin)[r, ] * Wm[i, ])
				gp[k] <- -val
			}
		}
	}
	g <- as.vector(crossprod(M, gp))
	zstar <- -as.vector(Sinv %*% g)
	zp <- as.vector(M %*% zstar)
	Sp <- M %*% Sinv %*% t(M)
	# the constant that the shifts add to the predictor of every cell of a period
	delta_at <- function(z) vapply(seq_len(Tt), function(t) {
		if (sym) sum(sg[seq_len(rmin)] * z[which(coords$t == t)]^2)
		else sum(z[which(coords$side == "u" & coords$t == t)] * z[which(coords$side == "v" & coords$t == t)])
	}, 0)
	e_delta <- delta_at(zp) + vapply(seq_len(Tt), function(t) {
		if (sym) sum(sg[seq_len(rmin)] * diag(Sp)[which(coords$t == t)])
		else sum(Sp[cbind(which(coords$side == "u" & coords$t == t), which(coords$side == "v" & coords$t == t))])
	}, 0)
	grad_delta <- function(t) {
		gd <- numeric(Kp)
		if (sym) { ix <- which(coords$t == t); gd[ix] <- 2 * sg[seq_len(rmin)] * zp[ix] }
		else {
			iu <- which(coords$side == "u" & coords$t == t); iv <- which(coords$side == "v" & coords$t == t)
			gd[iu] <- zp[iv]; gd[iv] <- zp[iu]
		}
		as.vector(crossprod(M, gd))
	}
	gdl <- lapply(seq_len(Tt), grad_delta)
	J <- length(cs_list)
	grad <- matrix(0, K, J)
	shift <- numeric(J)
	for (j in seq_len(J)) {
		sol <- .lde_solve(ev$fac, rep(list(cs_list[[j]]), Tt), rz_list[[j]])
		Gj <- vapply(seq_len(K), function(k) dot(B_rs[[k]], B_rz[[k]], sol), 0)
		hj <- vapply(seq_len(Tt), function(t) dot(R_rs[[t]], R_rz[[t]], sol), 0)
		grad[, j] <- -Gj - Reduce(`+`, Map(function(h, gd) h * gd, hj, gdl))
		shift[j] <- -sum(Gj * zstar) - sum(hj * e_delta)
	}
	out <- list(extra = crossprod(grad, Sinv %*% grad), shift = shift, zstar = zstar, e_delta = e_delta)
	if (isTRUE(details)) out <- c(out, list(S = S, Sinv = Sinv, coords = zc, g = g, D = D, BHB = BHB))
	out
}
