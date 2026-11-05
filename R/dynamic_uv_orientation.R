# orientation moves for the directed rank-two Gaussian AR(1) path: a common orthogonal map of U_t and V_t
# keeps their product and, given the neighboring periods, follows an O(2) matrix-Fisher law (two von Mises
# components weighted by I0(kappa)); von Mises rejection draw of Best and Fisher (1979), doi:10.2307/2346732

.lame_uv_von_mises <- function(kappa, mu = 0) {
	if (kappa < 1e-5) {
		repeat {
			theta <- stats::runif(1, -pi, pi)
			if (log(stats::runif(1)) <= kappa * (cos(theta) - 1)) return(theta + mu)
		}
	}
	a <- 1 + sqrt(1 + 4 * kappa * kappa)
	b <- 2 * kappa / (a + sqrt(2 * a))
	r <- (1 + b * b) / (2 * b)
	repeat {
		z <- cos(pi * stats::runif(1))
		f <- (1 + r * z) / (r + z)
		c <- kappa * (r - f)
		u <- stats::runif(1)
		if (u < c * (2 - c) || log(c / u) + 1 - c >= 0) {
			angle <- acos(pmin(1, pmax(-1, f)))
			return(mu + if (stats::runif(1) > 0.5) angle else -angle)
		}
	}
}

# log(exp(-x) I_0(x)); above 1e4 a four-term large-x expansion (DLMF 10.40.1), since besselI underflows there
.lame_uv_log_i0_scaled <- function(x) {
	if (any(!is.finite(x)) || any(x < 0)) stop("Invalid circular concentration.")
	out <- numeric(length(x))
	large <- x > 1e4
	out[!large] <- log(besselI(x[!large], 0, expon.scaled = TRUE))
	if (any(large)) {
		y <- x[large]
		term <- correction <- 1 / (8 * y)
		for (k in 2:4) {
			term <- term * (2 * k - 1)^2 / (8 * k * y)
			correction <- correction + term
		}
		out[large] <- -.5 * (log(2 * pi) + log(y)) + log1p(correction)
	}
	out
}

.lame_uv_o2 <- function(M) {
	# at extreme concentration keep the current orientation; the norm of M is constant on the rotation
	# orbit, so this choice preserves the conditional target
	if (!all(is.finite(M)) || sqrt(sum(M * M)) > 1e6) return(diag(2))
	a <- c(M[1, 1] + M[2, 2], M[1, 1] - M[2, 2])
	b <- c(M[2, 1] - M[1, 2], M[2, 1] + M[1, 2])
	kappa <- sqrt(a * a + b * b)
	log_mass <- .lame_uv_log_i0_scaled(kappa) + kappa
	positive <- stats::runif(1) < stats::plogis(log_mass[1] - log_mass[2])
	j <- if (positive) 1L else 2L
	theta <- .lame_uv_von_mises(kappa[j], atan2(b[j], a[j]))
	Q <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
	if (!positive) Q[, 2] <- -Q[, 2]
	Q
}

.lame_uv_common_orientation <- function(UV) {
	# rotate all periods by one common Procrustes map, which leaves the AR(1) prior unchanged
	U <- UV$U
	V <- UV$V
	n <- dim(U)[1L]
	Tt <- dim(U)[3L]
	avg_product <- matrix(0, n, n)
	avg_u <- matrix(0, n, 2)
	for (t in seq_len(Tt)) {
		ut <- matrix(U[, , t], n, 2)
		vt <- matrix(V[, , t], n, 2)
		avg_product <- avg_product + tcrossprod(ut, vt) / Tt
		avg_u <- avg_u + ut / Tt
	}
	rotation <- tryCatch({
		s <- svd(avg_product, nu = min(2L, n), nv = 0)
		k <- min(2L, length(s$d))
		ref <- matrix(0, n, 2)
		ref[, seq_len(k)] <- sweep(
			s$u[, seq_len(k), drop = FALSE],
			2,
			sqrt(s$d[seq_len(k)]),
			"*"
		)
		m <- svd(crossprod(ref, avg_u))
		m$v %*% t(m$u)
	}, error = function(e) NULL)
	if (is.null(rotation) || !all(is.finite(rotation))) return(UV)
	for (t in seq_len(Tt)) {
		U[, , t] <- matrix(U[, , t], n, 2) %*% rotation
		V[, , t] <- matrix(V[, , t], n, 2) %*% rotation
	}
	UV$U <- U
	UV$V <- V
	UV
}

# exact Gaussian conditionals along two shear directions: U_l + q U_k and V_k - q V_l keep every
# product with unit Jacobian, and the stationary AR(1) prior gives q a Gaussian conditional
.lame_uv_shear_ar1 <- function(UV, rho_uv, sigma_uv, max_abs = Inf) {
	U <- UV$U
	V <- UV$V
	Tt <- dim(U)[3L]
	innovation <- function(X) {
		z <- X
		z[, , 1L] <- sqrt(1 - rho_uv^2) * X[, , 1L]
		if (Tt > 1L) z[, , -1L] <- X[, , -1L] - rho_uv * X[, , -Tt]
		z
	}
	for (pair in list(c(1L, 2L), c(2L, 1L))) {
		k <- pair[1L]
		l <- pair[2L]
		du <- innovation(U)
		dv <- innovation(V)
		A <- sum(du[, k, ]^2) + sum(dv[, l, ]^2)
		B <- sum(du[, l, ] * du[, k, ]) - sum(dv[, k, ] * dv[, l, ])
		if (is.finite(A) && A > 0 && is.finite(B)) {
			q <- stats::rnorm(1, -B / A, sigma_uv / sqrt(A))
			stopifnot(is.finite(q))
			u_new <- U[, l, ] + q * U[, k, ]
			v_new <- V[, k, ] - q * V[, l, ]
			if (max(abs(u_new), abs(v_new)) <= max_abs) {
				U[, l, ] <- u_new
				V[, k, ] <- v_new
			}
		}
	}
	UV$U <- U
	UV$V <- V
	UV
}

.lame_uv_orientation_ar1 <- function(UV, rho_uv, sigma_uv, enabled = TRUE, align = TRUE, max_abs = Inf) {
	if (!isTRUE(enabled)) return(UV)
	U <- UV$U
	V <- UV$V
	stopifnot(
		length(dim(U)) == 3L,
		identical(dim(U), dim(V)),
		dim(U)[1L] > 0L,
		dim(U)[2L] == 2L,
		dim(U)[3L] > 0L,
		length(rho_uv) == 1L,
		is.finite(rho_uv),
		abs(rho_uv) < 1,
		length(sigma_uv) == 1L,
		is.finite(sigma_uv),
		sigma_uv > 0
	)
	n <- dim(U)[1L]
	Tt <- dim(U)[3L]
	for (t in seq_len(Tt)) {
		ut <- matrix(U[, , t], n, 2)
		vt <- matrix(V[, , t], n, 2)
		neighbor_u <- neighbor_v <- matrix(0, n, 2)
		if (t > 1L) {
			neighbor_u <- neighbor_u + U[, , t - 1L]
			neighbor_v <- neighbor_v + V[, , t - 1L]
		}
		if (t < Tt) {
			neighbor_u <- neighbor_u + U[, , t + 1L]
			neighbor_v <- neighbor_v + V[, , t + 1L]
		}
		M <- rho_uv / sigma_uv^2 *
			(crossprod(ut, neighbor_u) + crossprod(vt, neighbor_v))
		Q <- .lame_uv_o2(M)
		u_new <- ut %*% Q
		v_new <- vt %*% Q
		# reject a period's orbit draw that leaves the coordinate box
		if (max(abs(u_new), abs(v_new)) <= max_abs) {
			U[, , t] <- u_new
			V[, , t] <- v_new
		}
	}
	UV$U <- U
	UV$V <- V
	UV <- .lame_uv_shear_ar1(UV, rho_uv, sigma_uv, max_abs)
	if (isTRUE(align) && !is.finite(max_abs)) .lame_uv_common_orientation(UV) else UV
}

# joint draw of period signs, which preserves U_t G V_t' and the coordinate bounds; adjacent relative signs
# have masses exp(kappa_t * sign[t-1] * sign[t]), and the weights cover Gaussian, reset and Student innovations
.lame_uv_sign_ar1 <- function(UV, rho_uv, sigma_uv, weight_u = NULL, weight_v = NULL) {
	U <- UV$U
	V <- UV$V
	du <- dim(U)
	dv <- dim(V)
	stopifnot(
		length(du) == 3L,
		length(dv) == 3L,
		all(du > 0L),
		all(dv > 0L),
		du[3L] == dv[3L],
		length(rho_uv) == 1L,
		is.finite(rho_uv),
		abs(rho_uv) < 1,
		length(sigma_uv) == 1L,
		is.finite(sigma_uv),
		sigma_uv > 0,
		all(is.finite(U)),
		all(is.finite(V))
	)
	nt <- du[3L]
	if (nt < 2L) return(UV)
	um <- matrix(U, ncol = nt)
	vm <- matrix(V, ncol = nt)
	if (is.null(weight_u) && is.null(weight_v)) {
		kappa <- rho_uv / sigma_uv^2 *
			(colSums(um[, -nt, drop = FALSE] * um[, -1L, drop = FALSE]) +
					colSums(vm[, -nt, drop = FALSE] * vm[, -1L, drop = FALSE]))
	} else {
		stopifnot(
			identical(dim(weight_u), c(du[1L], nt)),
			identical(dim(weight_v), c(dv[1L], nt))
		)
		wu <- weight_u[rep(seq_len(du[1L]), du[2L]), -1L, drop = FALSE]
		wv <- weight_v[rep(seq_len(dv[1L]), dv[2L]), -1L, drop = FALSE]
		stopifnot(all(is.finite(wu) & wu >= 0), all(is.finite(wv) & wv >= 0))
		kappa <- rho_uv / sigma_uv^2 *
			(colSums(wu * um[, -nt, drop = FALSE] * um[, -1L, drop = FALSE]) +
					colSums(wv * vm[, -nt, drop = FALSE] * vm[, -1L, drop = FALSE]))
	}
	signs <- c(1, cumprod(ifelse(stats::runif(nt - 1L) < stats::plogis(2 * kappa), 1, -1)))
	UV$U <- array(sweep(um, 2, signs, "*"), du, dimnames = dimnames(U))
	UV$V <- array(sweep(vm, 2, signs, "*"), dv, dimnames = dimnames(V))
	UV
}

# signs along the singular directions of G preserve each bipartite product, and the isotropic AR(1)
# prior is unchanged by the common singular-vector bases
.lame_uv_svd_sign_ar1 <- function(UV, G, rho_uv, sigma_uv, max_abs = Inf) {
	U <- UV$U
	V <- UV$V
	du <- dim(U)
	dv <- dim(V)
	nt <- du[3L]
	stopifnot(
		length(du) == 3L,
		length(dv) == 3L,
		du[2L] == nrow(G),
		dv[2L] == ncol(G),
		nt == dv[3L]
	)
	if (nt < 2L || min(du[2L], dv[2L]) < 2L) return(UV)
	ss <- svd(G, nu = du[2L], nv = dv[2L])
	A <- U
	B <- V
	for (tt in seq_len(nt)) {
		A[, , tt] <- matrix(U[, , tt], du[1L], du[2L]) %*% ss$u
		B[, , tt] <- matrix(V[, , tt], dv[1L], dv[2L]) %*% ss$v
	}
	for (k in seq_len(min(du[2L], dv[2L]))) {
		a <- matrix(A[, k, ], du[1L], nt)
		b <- matrix(B[, k, ], dv[1L], nt)
		kappa <- rho_uv / sigma_uv^2 *
			(colSums(a[, -nt, drop = FALSE] * a[, -1L, drop = FALSE]) +
					colSums(b[, -nt, drop = FALSE] * b[, -1L, drop = FALSE]))
		signs <- c(1, cumprod(ifelse(stats::runif(nt - 1L) < stats::plogis(2 * kappa), 1, -1)))
		aa <- A
		bb <- B
		aa[, k, ] <- sweep(a, 2, signs, "*")
		bb[, k, ] <- sweep(b, 2, signs, "*")
		u_new <- U
		v_new <- V
		for (tt in seq_len(nt)) {
			u_new[, , tt] <- matrix(aa[, , tt], du[1L], du[2L]) %*% t(ss$u)
			v_new[, , tt] <- matrix(bb[, , tt], dv[1L], dv[2L]) %*% t(ss$v)
		}
		# reject this orbit draw if it leaves the coordinate box
		if (all(is.finite(u_new)) && all(is.finite(v_new)) && max(abs(u_new), abs(v_new)) <= max_abs) {
			A <- aa
			B <- bb
			U <- u_new
			V <- v_new
		}
	}
	UV$U <- U
	UV$V <- V
	UV
}
