# variational (mean-field) posterior of the latent factors of the EB fit, started from the posterior mode
# (R/ame_als_lowrank.R). q is Gaussian and independent across actors and sides, q(u_i) = N(uhat_i, Sig_i),
# q(v_j) = N(vhat_j, Om_j); the additive effects are Gaussian given the coefficients, and the coefficients and
# scales are point estimates. a cell's predictor is Gaussian with the exact moments under q,
#   m_c = x_c' theta + a_i + b_j + uhat_i' vhat_j,  V_c = var(a_i + b_j) + vprod_ij,
#   vprod_ij = uhat_i' Om_j uhat_i + vhat_j' Sig_i vhat_j + tr(Sig_i Om_j),
# and its expected score s_c and information w_c (by quadrature) give a Gaussian working model with response
# z_c = m_c + s_c / w_c, whose fixed point is a stationary point of the lower bound. factor rows:
# Sig_i = (Lam_uu + sum_j C_ij (vhat_j vhat_j' + Om_j))^-1, uhat_i = Sig_i (sum_j C_ij Rbar_ij vhat_j - Lam_uv vhat_i),
# receivers likewise; an undirected fit has one factor F with signs S, P = F S F'. the plug-in scale rule
# updates Suv = (Psi + SE) / (df + n) from mean-field moments with a gauge step (.ae_mf_gauge); the cavity rule
# uses cavity covariances (.ae_mf_cavity_cov) and a Newton step (.ae_mf_scale_step). ranks are compared by
# the Laplace evidence (R/ame_als_eb.R), not `fit$eb$variational$elbo`; R/ame_als_pair.R handles dyadic correlation


# finest pair rule (nodes a side) the iteration refines to when a stall shows the rule's resolution
# limiting; fits start at five
.ae_pair_quad_max <- 9L

.ae_mf_gh_cache <- new.env(parent = emptyenv())

# Gauss-Hermite rule for the standard normal: nodes z and weights summing to one
.ae_mf_gh <- function(K = 9L) {
	key <- as.character(K)
	hit <- .ae_mf_gh_cache[[key]]
	if (!is.null(hit)) return(hit)
	i <- seq_len(K - 1L)
	J <- matrix(0, K, K)
	J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- sqrt(i / 2)
	e <- eigen(J, symmetric = TRUE)
	out <- list(z = sqrt(2) * e$values, w = e$vectors[1L, ]^2)
	.ae_mf_gh_cache[[key]] <- out
	out
}

# quadrature engine: compiled kernels (src/ae_kernels.cpp) by default; options(lame.als_engine = "r") runs the
# _r reference implementations, and options(lame.als_threads = k) sets OpenMP threads without changing the result
.ae_engine <- function() {
	opt <- getOption("lame.als_engine", "cpp")
	if (!is.character(opt) || length(opt) != 1L || !opt %in% c("cpp", "r")) {
		cli::cli_abort("{.code options(lame.als_engine)} must be {.val cpp} or {.val r}.")
	}
	opt
}
.ae_threads <- function() {
	opt <- getOption("lame.als_threads", 1L)
	if (!is.numeric(opt) || length(opt) != 1L || !is.finite(opt) || opt < 1) {
		cli::cli_abort("{.code options(lame.als_threads)} must be a positive whole number.")
	}
	as.integer(opt)
}

# expected working problem of the observed cells of `y` with predictor N(m, v): response z = m + E psi / w,
# weight w = -E psi' and expected log likelihood ell, zero elsewhere; the increment is capped at 30 as in IRLS
.ae_mf_cells_r <- function(y, m, v, link, K = 9L) {
	obs <- is.finite(y)
	ix <- which(obs)
	sgn <- 2 * y[ix] - 1
	mm <- m[ix]
	nodes <- .ae_mf_gh(K)
	sdv <- sqrt(pmax(v[ix], 0))
	probit <- identical(link, "probit")
	es <- ew <- el <- 0
	for (k in seq_along(nodes$z)) {
		e <- mm + sdv * nodes$z[k]
		se <- sgn * e
		if (probit) {
			lps <- stats::pnorm(se, log.p = TRUE)
			lam <- exp(-0.5 * e * e - 0.91893853320467267 - lps)
			psi <- sgn * lam
			w <- lam * (lam + se)
		} else {
			lps <- stats::plogis(se, log.p = TRUE)
			p1 <- exp(lps)
			psi <- sgn * (1 - p1)
			w <- p1 * (1 - p1)
		}
		es <- es + nodes$w[k] * psi
		ew <- ew + nodes$w[k] * w
		el <- el + nodes$w[k] * lps
	}
	delta <- es / pmax(ew, 1e-12)
	big <- abs(delta) > 30
	if (any(big)) {
		delta[big] <- sign(delta[big]) * 30
		ew[big] <- abs(es[big]) / 30
	}
	z <- w <- ell <- array(0, dim(y))
	z[ix] <- mm + delta
	w[ix] <- ew
	ell[ix] <- el
	list(z = z, w = w, ell = ell)
}

# the same by the compiled kernel (ae_cells_cpp)
.ae_mf_cells <- function(y, m, v, link, K = 9L, ve = NULL, kin = NULL) {
	# a count cell is Poisson-lognormal (R/ame_als_pln.R), with noise variance `ve` and an inner rule of `kin` nodes
	if (identical(link, "log")) return(.ae_pln_cells(y, m, v, ve, K, if (is.null(kin)) .ae_pln_kin(ve) else kin))
	# a normal cell is Gaussian, in closed form (R/ame_als_gauss.R)
	if (identical(link, "identity")) return(.ae_gauss_cells(y, m, v, ve))
	if (identical(.ae_engine(), "r")) return(.ae_mf_cells_r(y, m, v, link, K))
	out <- ae_cells_cpp(y, m, v, identical(link, "probit"), as.integer(K), .ae_threads())
	d <- dim(y)
	dim(out$z) <- d; dim(out$w) <- d; dim(out$ell) <- d
	out
}

# cavity information weight of a cell for the factor-row covariances, nu = -d^2 log Z / dm^2 with
# Z(m, v) = E p(y | eta), eta ~ N(m, v); the expected weight E w(eta) would count the other side's variance twice
.ae_mf_info_w_r <- function(y, m, v, link, K = 9L) {
	obs <- is.finite(y)
	ix <- which(obs)
	sgn <- 2 * y[ix] - 1
	mm <- m[ix]
	vv <- pmax(v[ix], 0)
	w <- if (identical(link, "probit")) {
		sq <- sqrt(1 + vv)
		t <- sgn * mm / sq
		lam <- exp(stats::dnorm(t, log = TRUE) - stats::pnorm(t, log.p = TRUE))
		out <- lam * (lam + t)
		far <- t < -30
		if (any(far)) out[far] <- 1 - 1 / t[far]^2 + 6 / t[far]^4
		out / (1 + vv)
	} else {
		nodes <- .ae_mf_gh(K)
		sdv <- sqrt(vv)
		z0 <- z1 <- z2 <- 0
		for (k in seq_along(nodes$z)) {
			p1 <- stats::plogis(sgn * (mm + sdv * nodes$z[k]))
			z0 <- z0 + nodes$w[k] * p1
			z1 <- z1 + nodes$w[k] * sgn * p1 * (1 - p1)
			z2 <- z2 + nodes$w[k] * p1 * (1 - p1) * (1 - 2 * p1)
		}
		(z1 / z0)^2 - z2 / z0
	}
	w[!is.finite(w) | w < 0] <- 0
	out <- array(0, dim(y))
	out[ix] <- w
	out
}

.ae_mf_info_w <- function(y, m, v, link, K = 9L, ve = NULL) {
	if (identical(link, "log")) return(.ae_pln_info_w(y, m, v, ve))
	if (identical(link, "identity")) return(.ae_gauss_info_w(y, m, v, ve))
	if (identical(.ae_engine(), "r")) return(.ae_mf_info_w_r(y, m, v, link, K))
	out <- ae_info_w_cpp(y, m, v, identical(link, "probit"), as.integer(K), .ae_threads())
	dim(out) <- dim(y)
	out
}

# the summed information weights of the cells of an evaluation (one value per dyad)
.ae_mf_info_sum <- function(cx, ev) {
	if (!is.null(ev$pair)) return(rowSums(.ae_mf_info_pair(cx, ev)$w, dims = 2L))
	rowSums(.ae_mf_info_w(cx$y, ev$eta, array(ev$V, dim(cx$y)), cx$link, ve = ev$ve), dims = 2L)
}

# the cavity weights of the cells of an evaluation of a fit with a dyadic correlation: the
# diagonal weights and the cross weights of the pairs (.ae_pair_info_w)
.ae_mf_info_pair <- function(cx, ev) {
	if (cx$gauss) return(.ae_gauss_pair_info_w(cx, ev$eta, array(ev$V, dim(cx$y)), ev$pair$cv, ev$pair$rho, ev$ve))
	.ae_pair_info_w(cx, ev$eta, array(ev$V, dim(cx$y)), ev$pair$cv, ev$pair$rho)
}


.ae_mf_array_to_rows <- function(A) {
	R <- dim(A)[1L]
	t(matrix(A, R * R, dim(A)[3L]))
}

# the row products of an n x R matrix, one row per actor, columns (k, l) with k fastest
.ae_mf_outer_rows <- function(M) {
	R <- ncol(M)
	M[, rep(seq_len(R), times = R), drop = FALSE] * M[, rep(seq_len(R), each = R), drop = FALSE]
}

# solve n small systems (G_i + L) x_i = rhs_i (G_i the rows of `Gflat`, L a common prior precision),
# returning solutions and inverses; ranks one and two are closed forms, larger ranks use ae_row_solve_cpp and
# fall back to .ae_safe_inverse for rows that are not positive definite
.ae_mf_row_solve <- function(Gflat, L, rhs) {
	n <- nrow(rhs); R <- ncol(rhs)
	if (R == 1L) {
		s <- 1 / (Gflat[, 1L] + L[1L, 1L])
		return(list(x = matrix(s * rhs[, 1L], n, 1L), Sig = array(s, c(1L, 1L, n))))
	}
	if (R == 2L) {
		g11 <- Gflat[, 1L] + L[1L, 1L]
		g22 <- Gflat[, 4L] + L[2L, 2L]
		g12 <- 0.5 * (Gflat[, 2L] + Gflat[, 3L]) + 0.5 * (L[1L, 2L] + L[2L, 1L])
		det <- g11 * g22 - g12 * g12
		s11 <- g22 / det; s22 <- g11 / det; s12 <- -g12 / det
		Sig <- array(0, c(2L, 2L, n))
		Sig[1L, 1L, ] <- s11; Sig[2L, 2L, ] <- s22
		Sig[1L, 2L, ] <- s12; Sig[2L, 1L, ] <- s12
		return(list(
			x = cbind(
				s11 * rhs[, 1L] + s12 * rhs[, 2L],
				s12 * rhs[, 1L] + s22 * rhs[, 2L]
			),
			Sig = Sig
		))
	}
	if (identical(.ae_engine(), "r")) return(.ae_mf_row_solve_r(Gflat, L, rhs))
	out <- ae_row_solve_cpp(Gflat, L, rhs, .ae_threads())
	for (i in out$bad) {
		g <- matrix(Gflat[i, ], R, R) + L
		g <- (g + t(g)) / 2
		Si <- .ae_safe_inverse(g)
		out$Sig[, , i] <- Si
		out$x[i, ] <- Si %*% rhs[i, ]
	}
	list(x = out$x, Sig = out$Sig)
}

# R reference for the general rank, one system at a time
.ae_mf_row_solve_r <- function(Gflat, L, rhs) {
	n <- nrow(rhs); R <- ncol(rhs)
	x <- matrix(0, n, R); Sig <- array(0, c(R, R, n))
	for (i in seq_len(n)) {
		g <- matrix(Gflat[i, ], R, R) + L
		g <- (g + t(g)) / 2
		Si <- .ae_safe_inverse(g)
		Sig[, , i] <- Si
		x[i, ] <- Si %*% rhs[i, ]
	}
	list(x = x, Sig = Sig)
}

# log determinants of the n small positive definite matrices of an (R, R, n) array
.ae_mf_logdet_rows <- function(A) {
	R <- dim(A)[1L]; n <- dim(A)[3L]
	if (R == 1L) return(log(A[1L, 1L, ]))
	if (R == 2L) return(log(A[1L, 1L, ] * A[2L, 2L, ] - A[1L, 2L, ] * A[2L, 1L, ]))
	vapply(seq_len(n), function(i) {
		as.numeric(determinant(matrix(A[, , i], R, R), logarithm = TRUE)$modulus)
	}, numeric(1))
}


# var(u_i' v_j) for every pair under independent Gaussian q, from the R x R row covariances Sig (senders)
# and Om (receivers)
.ae_mf_prod_var <- function(U, V, Sig, Om) {
	UU <- .ae_mf_outer_rows(U); VV <- .ae_mf_outer_rows(V)
	SG <- .ae_mf_array_to_rows(Sig); OM <- .ae_mf_array_to_rows(Om)
	tcrossprod(UU, OM) + tcrossprod(SG, VV) + tcrossprod(SG, OM)
}

# the same for an undirected network: P_ij = f_i' S f_j with the signs `sg`,
#   var = fhat_i' S Sig_j S fhat_i + fhat_j' S Sig_i S fhat_j + tr(Sig_i S Sig_j S)
.ae_mf_prod_var_sym <- function(Fm, sg, Sig) {
	R <- ncol(Fm)
	ss <- tcrossprod(sg)[cbind(rep(seq_len(R), times = R), rep(seq_len(R), each = R))]
	FF <- .ae_mf_outer_rows(Fm)
	SR <- .ae_mf_array_to_rows(Sig)
	SG <- sweep(SR, 2L, ss, "*")
	tcrossprod(FF, SG) + tcrossprod(SG, FF) + tcrossprod(SR, SG)
}


# the precision blocks (uu, vv, uv) of a prior
.ae_mf_prior_blocks <- function(prior, R) {
	if (identical(prior$kind, "joint")) {
		Lam <- prior$Lam
		iu <- seq_len(R); iv <- R + seq_len(R)
		return(list(
			uu = Lam[iu, iu, drop = FALSE],
			vv = Lam[iv, iv, drop = FALSE],
			uv = Lam[iu, iv, drop = FALSE]
		))
	}
	if (prior$kind %in% c("scalar", "columns")) {
		return(list(
			uu = diag(1 / rep_len(prior$tau, R), R),
			vv = diag(1 / rep_len(prior$tau, R), R),
			uv = matrix(0, R, R)
		))
	}
	list(uu = diag(1 / prior$su, R), vv = diag(1 / prior$sv, R), uv = matrix(0, R, R))
}

# the second moments of q the scales are updated from: the stacked second
# moment SE of (u_i, v_i) for the joint prior, the per-column sums otherwise
.ae_mf_second_moments <- function(S, R) {
	if (!is.null(S$F)) {
		return(list(
			total = sum(S$F^2) + sum(apply(S$Sig, 3L, function(s) sum(diag(s)))),
			n = nrow(S$F) * R,
			col = colSums(S$F^2) + vapply(seq_len(R), function(k) sum(S$Sig[k, k, ]), 0),
			rows = nrow(S$F)
		))
	}
	su <- colSums(S$U^2) + vapply(seq_len(R), function(k) sum(S$Sig[k, k, ]), 0)
	sv <- colSums(S$V^2) + vapply(seq_len(R), function(k) sum(S$Om[k, k, ]), 0)
	out <- list(su = su, sv = sv, nu = nrow(S$U), nv = nrow(S$V))
	if (identical(S$prior$kind, "joint")) {
		SE <- crossprod(cbind(S$U, S$V))
		iu <- seq_len(R); iv <- R + seq_len(R)
		SE[iu, iu] <- SE[iu, iu] + apply(S$Sig, c(1L, 2L), sum)
		SE[iv, iv] <- SE[iv, iv] + apply(S$Om, c(1L, 2L), sum)
		out$SE <- (SE + t(SE)) / 2
	}
	out
}

# the variational update of the scale(s) of the prior from the second moments
.ae_mf_prior_update <- function(S, R) {
	prior <- S$prior
	sm <- .ae_mf_second_moments(S, R)
	if (identical(prior$kind, "scalar")) {
		prior$tau <- max((prior$suv0 + sm$total) / (sm$n + prior$kappa0), 1e-10)
	} else if (identical(prior$kind, "joint")) {
		prior$Suv <- (prior$psi + sm$SE) / (sm$nu + prior$kappa0)
		prior$Suv <- (prior$Suv + t(prior$Suv)) / 2
		prior$Lam <- .ae_safe_inverse(prior$Suv)
	} else {
		prior$su <- pmax((prior$suv0 + sm$su) / (sm$nu + prior$kappa0), 1e-10)
		prior$sv <- pmax((prior$suv0 + sm$sv) / (sm$nv + prior$kappa0), 1e-10)
	}
	prior
}


# the conditional covariances of the factor rows at fixed means (the start)
.ae_mf_cov0 <- function(U, V, C, Luu, Lvv) {
	R <- ncol(U)
	Gu <- C %*% .ae_mf_outer_rows(V)
	Gv <- crossprod(C, .ae_mf_outer_rows(U))
	list(
		Sig = .ae_mf_row_solve(Gu, Luu, matrix(0, nrow(U), R))$Sig,
		Om = .ae_mf_row_solve(Gv, Lvv, matrix(0, nrow(V), R))$Sig
	)
}

# sweeps for a directed or bipartite network at the weights C and the targets
# Rbar, with the scales fixed: all senders given the receivers' means and
# covariances, then all receivers given the senders'. the means and covariances
# enter the other block's system through the second moment vhat vhat' + Om
.ae_mf_sweeps <- function(Rbar, C, U, V, Sig, Om, prior, inner = 6L, tol = 1e-9) {
	R <- ncol(U)
	pb <- .ae_mf_prior_blocks(prior, R)
	joint <- identical(prior$kind, "joint")
	RC <- Rbar * C
	RC[C <= 0] <- 0
	for (it in seq_len(inner)) {
		U0 <- U; V0 <- V
		EVV <- .ae_mf_outer_rows(V) + .ae_mf_array_to_rows(Om)
		rhs <- RC %*% V
		if (joint) rhs <- rhs - V %*% t(pb$uv)
		su <- .ae_mf_row_solve(C %*% EVV, pb$uu, rhs)
		U <- su$x; Sig <- su$Sig
		EUU <- .ae_mf_outer_rows(U) + .ae_mf_array_to_rows(Sig)
		rhs <- crossprod(RC, U)
		if (joint) rhs <- rhs - U %*% pb$uv
		sv <- .ae_mf_row_solve(crossprod(C, EUU), pb$vv, rhs)
		V <- sv$x; Om <- sv$Sig
		if (max(abs(U - U0), abs(V - V0)) < tol) break
	}
	list(U = U, V = V, Sig = Sig, Om = Om, sweeps = it)
}

# the same for an undirected network: P_ij = f_i' S f_j, the rows updated in
# turn (each row is the exact conditional maximiser given the others), one
# scalar scale. the weights C are those of the unordered dyads
.ae_mf_sweeps_sym <- function(Rbar, C, Fm, sg, Sig, prior, inner = 6L, tol = 1e-9) {
	n <- nrow(Fm); R <- ncol(Fm)
	Cn <- C; diag(Cn) <- 0
	RC <- Rbar * Cn
	RC[Cn <= 0] <- 0
	ridge <- diag(1 / prior$tau, R)
	ss <- tcrossprod(sg)[cbind(rep(seq_len(R), times = R), rep(seq_len(R), each = R))]
	for (it in seq_len(inner)) {
		F0 <- Fm
		Fs <- sweep(Fm, 2L, sg, "*")
		EFF <- .ae_mf_outer_rows(Fs) + sweep(.ae_mf_array_to_rows(Sig), 2L, ss, "*")
		Gf <- Cn %*% EFF
		rhs <- RC %*% Fs
		for (i in seq_len(n)) {
			Si <- .ae_safe_inverse(matrix(Gf[i, ], R, R) + ridge)
			fi <- as.vector(Si %*% rhs[i, ])
			Fm[i, ] <- fi
			Sig[, , i] <- Si
			fs_new <- sg * fi
			eff_new <- fs_new[rep(seq_len(R), times = R)] * fs_new[rep(seq_len(R), each = R)] +
				ss * as.vector(Si)
			Gf <- Gf + tcrossprod(Cn[, i], eff_new - EFF[i, ])
			rhs <- rhs + tcrossprod(RC[, i], fs_new - Fs[i, ])
			EFF[i, ] <- eff_new
			Fs[i, ] <- fs_new
		}
		if (max(abs(Fm - F0)) < tol) break
	}
	list(F = Fm, Sig = Sig, sweeps = it)
}


.ae_mf_logdet <- function(A) {
	ch <- tryCatch(chol(A), error = function(e) NULL)
	if (!is.null(ch)) return(2 * sum(log(diag(ch))))
	as.numeric(determinant((A + t(A)) / 2, logarithm = TRUE)$modulus)
}

# the variance of the additive part of the predictor of every cell from the
# covariance `Sc` of the stacked actor effects
.ae_mf_vact <- function(Sc, ab) {
	nr <- ab$nr
	dg <- diag(Sc)
	if (identical(ab$kind, "symmetric")) return(outer(dg, dg, "+") + 2 * Sc)
	ia <- seq_len(nr); ib <- nr + seq_len(ab$nc)
	outer(dg[ia], dg[ib], "+") + 2 * Sc[ia, ib, drop = FALSE]
}

# factor part of the objective with the scales integrated out: the entropy of q (up to a constant) plus
# the expected log prior of the rows under the inverse-Wishart (inverse-gamma) hyperprior, -(df + n) / 2 log |Psi + SE|
.ae_mf_factor_term <- function(S, R) {
	sm <- .ae_mf_second_moments(S, R)
	pr <- S$prior
	if (!is.null(S$F)) {
		return(
			-0.5 * (pr$kappa0 + sm$n) * log(pr$suv0 + sm$total) +
				0.5 * sum(.ae_mf_logdet_rows(S$Sig))
		)
	}
	ent <- 0.5 * (sum(.ae_mf_logdet_rows(S$Sig)) + sum(.ae_mf_logdet_rows(S$Om)))
	if (identical(pr$kind, "joint")) {
		return(-0.5 * (pr$kappa0 + sm$nu) * .ae_mf_logdet(pr$psi + sm$SE) + ent)
	}
	-0.5 * (pr$kappa0 + sm$nu) * sum(log(pr$suv0 + sm$su)) -
		0.5 * (pr$kappa0 + sm$nv) * sum(log(pr$suv0 + sm$sv)) + ent
}

# KL divergences of q from the priors at the plug-in scales, for the actor effects (covariance `Sc`,
# means `m_s`, prior mean `mu0`, precision `Lam`) and the factor rows
.ae_mf_kl_actors <- function(Sc, m_s, mu0, Lam) {
	d <- m_s - mu0
	0.5 * (sum(Lam * Sc) + sum(d * as.vector(Lam %*% d)) - length(d) -
			.ae_mf_logdet(Lam) - .ae_mf_logdet(Sc))
}
.ae_mf_kl_factors <- function(S, R) {
	pr <- S$prior
	sm <- .ae_mf_second_moments(S, R)
	if (!is.null(S$F)) {
		n <- nrow(S$F)
		if (identical(pr$kind, "columns")) {
			return(
				0.5 * (sum(sm$col / pr$tau) - n * R + n * sum(log(pr$tau)) -
						sum(.ae_mf_logdet_rows(S$Sig)))
			)
		}
		return(
			0.5 * (sm$total / pr$tau - n * R + n * R * log(pr$tau) -
					sum(.ae_mf_logdet_rows(S$Sig)))
		)
	}
	nu <- sm$nu; nv <- sm$nv
	ld <- sum(.ae_mf_logdet_rows(S$Sig)) + sum(.ae_mf_logdet_rows(S$Om))
	if (identical(pr$kind, "joint")) {
		Lam <- pr$Lam
		W <- cbind(S$U, S$V)
		quad <- sum((W %*% Lam) * W) +
			sum(Lam[seq_len(R), seq_len(R)] * apply(S$Sig, c(1L, 2L), sum)) +
			sum(Lam[R + seq_len(R), R + seq_len(R)] * apply(S$Om, c(1L, 2L), sum))
		return(0.5 * (quad - 2 * R * nu + nu * .ae_mf_logdet(pr$Suv) - ld))
	}
	0.5 * (sum(sm$su / pr$su) + sum(sm$sv / pr$sv) - R * (nu + nv) +
			nu * sum(log(pr$su)) + nv * sum(log(pr$sv)) - ld)
}


# penalised solve of (mu, beta, actor effects, gamma) at the expected weights and working response, with
# the conditional covariance `Sc` of the effects; no Laplace correction, since q's variance is in the weights
.ae_mf_actor_step <- function(ctx, ab, Omat, st, sab, cross = NULL) {
	nr <- ctx$nr; nc <- ctx$nc; pg <- ctx$pg
	nu <- ab$nu
	sym <- identical(ab$kind, "symmetric")
	Lam <- .ae_eb_precision(sab, ab)
	H <- .ae_eb_penalise(ctx, Lam)
	g0 <- st$xtz; rs <- st$szr; cs <- st$szc
	if (!is.null(Omat)) {
		g0 <- g0 - ae_bcd_gtc_cpp(numeric(nr), numeric(nc), Omat, st$gm)
		wo <- ctx$wsum * Omat
		rs <- rs - rowSums(wo); cs <- cs - colSums(wo)
		# the information of a pair couples a cell to the product of the other cell of its dyad
		if (!is.null(cross)) {
			g0 <- g0 - ae_bcd_gtc_cpp(numeric(nr), numeric(nc), t(Omat), cross$gmc)
			wo2 <- cross$Csum * t(Omat)
			rs <- rs - rowSums(wo2); cs <- cs - colSums(wo2)
		}
	}
	rhs <- c(nu * g0, if (sym) nu * (rs + cs) else c(rs, cs), numeric(pg))
	s <- ctx$i_s
	inv <- .ae_eb_inverse_pair(H, s)
	S <- inv$S
	x <- as.vector(S %*% rhs)
	gam <- x[ctx$i_g]
	alpha <- x[s] - if (pg > 0L) as.vector(ctx$M %*% gam) else 0
	list(
		theta = x[ctx$i_th],
		a = x[s][seq_len(nr)],
		b = if (sym) x[s][seq_len(nr)] else x[s][nr + seq_len(nc)],
		gamma = gam,
		alpha = alpha,
		Sc = inv$Sc,
		Lam = Lam,
		sab_new = .ae_eb_update_sab(ab, alpha, .ae_eb_css(ctx, S))
	)
}


# the constants of the model the iteration shares: the data, the design, the
# prior state of the effects, the rank and the quadrature. `hyper` overrides
# the factor prior (development only): kind ("joint", "diag"), suv0 and kappa0
.ae_mf_context <- function(prep, ab, R, link, quad = 9L, hyper = NULL, pair = NULL) {
	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt; p <- prep$p_dyad
	y <- prep$Yarr; X <- prep[["X"]]
	suv0 <- .ae_uv_prior_scale(y, if (identical(link, "identity")) "normal" else "binary")
	kappa0 <- 2 * R + 2
	# the stacked rows of a directed network take the inverse-Wishart prior of ame() and lame() (R/uv_prior.R);
	# the per-column scales of the other kinds keep the constants of the mode
	uvp <- if (R > 0L && identical(ab$kind, "directed")) {
		.lame_uv_prior(
			y,
			if (identical(link, "identity")) "normal" else if (identical(link, "log")) "poisson" else "binary",
			R
		)
	}
	suv_joint <- if (!is.null(uvp)) uvp$scale else if (identical(link, "identity")) .lame_uv_prior_scale(y, "normal") else suv0
	kappa0_joint <- if (!is.null(uvp)) uvp$df else .lame_uv_prior_df(R)
	psi_joint <- if (!is.null(uvp)) uvp$psi else kappa0_joint * suv_joint
	if (!is.null(hyper$suv0)) suv0 <- suv_joint <- hyper$suv0
	if (!is.null(hyper$kappa0)) kappa0 <- kappa0_joint <- hyper$kappa0
	if (!is.null(hyper$suv0) || !is.null(hyper$kappa0)) psi_joint <- kappa0_joint * suv_joint
	cx <- list(
		nr = nr,
		nc = nc,
		Tt = Tt,
		p = p,
		y = y,
		X = X,
		obs = is.finite(y),
		obs_dyad = rowSums(is.finite(y), dims = 2L) > 0,
		xmat = if (p > 0L) matrix(aperm(X, c(1L, 2L, 4L, 3L)), nr * nc * Tt, p),
		ab = ab,
		nu = ab$nu,
		sym = identical(ab$kind, "symmetric"),
		M = .ae_eb_mean_design(ab),
		R = R,
		link = link,
		quad = quad,
		pln = identical(link, "log"),
		gauss = identical(link, "identity"),
		has_ve = link %in% c("log", "identity"),
		n_obs = sum(is.finite(y)),
		# response scale the iteration residual is measured in
		ysd = if (identical(link, "identity")) max(stats::sd(y[is.finite(y)]), .Machine$double.eps) else 1,
		suv_joint = suv_joint,
		suv0 = suv0,
		kappa0 = kappa0,
		kappa0_joint = kappa0_joint,
		psi_joint = psi_joint,
		joint = identical(ab$kind, "directed") && !identical(hyper$kind, "diag")
	)
	# dyadic correlation of a directed probit fit: pairs of cells (R/ame_als_pair.R)
	if (!is.null(pair) && identical(ab$kind, "directed") && identical(link, "probit")) {
		cx$pair <- .ae_pair_context(y)
		cx$pair_quad <- if (is.null(pair$quad)) 5L else pair$quad
		cx$pair_rho0 <- pair$rho
		cx$pair_update <- !isFALSE(pair$update)
	}
	# pairs of a directed normal network are bivariate normal in closed form (R/ame_als_gauss_pair.R)
	if (!is.null(pair) && identical(ab$kind, "directed") && identical(link, "identity")) {
		cx$pair <- .ae_pair_context(y)
		cx$pair_quad <- 0L
		cx$pair_rho0 <- pair$rho
		cx$pair_update <- !isFALSE(pair$update)
	}
	cx
}

# restart a factor column the mode set to zero (a zero column stays zero under the product updates) from a
# ridge-penalised rank-one ALS fit to the cell scores at the mode; the variational iteration keeps it only
# where it raises the lower bound
.ae_mf_reseed <- function(cx, S, fit) {
	R <- cx$R
	nrm <- if (cx$sym) sqrt(colSums(S$F^2)) else sqrt(colSums(S$U^2)) * sqrt(colSums(S[["V"]]^2))
	dead <- which(!(nrm > 1e-6 * max(nrm)))
	if (!length(dead) || !is.finite(max(nrm)) || max(nrm) <= 0) return(S)
	pr <- .ae_mf_predictor(cx, S)
	tm <- .ae_pln_terms(cx$y, pr$eta, S$ve)
	sc <- tm$s; sc[!is.finite(sc)] <- 0
	wc <- tm$w; wc[!is.finite(wc)] <- 0
	Sc <- rowSums(sc, dims = 2L); Wc <- rowSums(wc, dims = 2L)
	s0 <- if (cx$sym) mean(S$prior$tau) else (sum(fit$U^2) + sum(fit[["V"]]^2)) / (2 * nrow(fit$U) * R)
	lam <- 1 / max(s0, 1e-3)
	if (cx$sym) {
		Sc <- (Sc + t(Sc)) / 2; Wc <- (Wc + t(Wc)) / 2; diag(Sc) <- 0; diag(Wc) <- 0
		for (k in dead) {
			eg <- eigen(Sc, symmetric = TRUE)
			top <- which.max(abs(eg$values))
			sg <- if (eg$values[top] < 0) -1 else 1
			f <- eg$vectors[, top] * sqrt(max(s0, 1e-3) * nrow(Sc))
			for (it in 1:10) f <- sg * as.vector(Sc %*% f) / (lam + as.vector(Wc %*% f^2))
			S$F[, k] <- f; S$sg[k] <- sg
		}
		Fs <- sweep(S$F, 2L, S$sg, "*")
		S$Sig <- .ae_mf_row_solve(Wc %*% .ae_mf_outer_rows(Fs), diag(1 / S$prior$tau, R), matrix(0, nrow(Fs), R))$Sig
		return(S)
	}
	for (k in dead) {
		sv <- svd(Sc, 1, 1)
		v <- sv$v[, 1L] * sqrt(max(s0, 1e-3) * length(sv$v[, 1L]))
		u <- numeric(nrow(Sc))
		for (it in 1:10) {
			u <- as.vector(Sc %*% v) / (lam + as.vector(Wc %*% v^2))
			v <- as.vector(crossprod(Sc, u)) / (lam + as.vector(crossprod(Wc, u^2)))
		}
		S$U[, k] <- u; S[["V"]][, k] <- v
	}
	pb <- .ae_mf_prior_blocks(S$prior, R)
	cv <- .ae_mf_cov0(S$U, S[["V"]], Wc, pb$uu, pb$vv)
	S$Sig <- cv$Sig; S$Om <- cv$Om
	S
}

# start at the mode, with the covariances it implies at its Fisher weights
.ae_mf_init <- function(cx, fit, W0) {
	ab <- cx$ab; R <- cx$R
	W0[!is.finite(W0)] <- 0
	C0 <- rowSums(W0, dims = 2L)
	S <- list(
		mu = fit$mu,
		beta = as.numeric(fit$beta),
		a = as.numeric(fit$a),
		b = as.numeric(fit$b),
		gamma = as.numeric(fit$gamma),
		sab = fit$sab
	)
	# noise variance of a count or normal cell, a hyperparameter of the state
	if (cx$has_ve) S$ve <- fit$ve
	S$Sc <- .ae_eb_inverse(.ae_eb_info_s(C0, ab) + .ae_eb_precision(S$sab, ab))
	S$vact <- .ae_mf_vact(S$Sc, ab)
	if (R == 0L) return(S)
	if (cx$sym) {
		S$sg <- ifelse(fit$L < 0, -1, 1)
		S$F <- sweep(fit$U, 2L, sqrt(abs(fit$L)), "*")
		S$prior <- list(
			kind = "scalar",
			tau = fit$sig2_u[1L],
			suv0 = cx$suv0,
			kappa0 = cx$kappa0
		)
		Cn <- C0; diag(Cn) <- 0
		Fs <- sweep(S$F, 2L, S$sg, "*")
		S$Sig <- .ae_mf_row_solve(
			Cn %*% .ae_mf_outer_rows(Fs),
			diag(1 / S$prior$tau, R),
			matrix(0, nrow(Fs), R)
		)$Sig
		return(if (cx$pln) .ae_mf_reseed(cx, S, fit) else S)
	}
	S$U <- fit$U; S$V <- fit[["V"]]
	S$prior <- if (cx$joint) {
		sc <- c(fit$sig2_u, fit$sig2_v)
		list(
			kind = "joint",
			Suv = diag(sc, 2L * R),
			Lam = diag(1 / sc, 2L * R),
			psi = cx$psi_joint * diag(2L * R),
			suv0 = cx$suv0,
			kappa0 = cx$kappa0_joint
		)
	} else {
		list(
			kind = "diag",
			su = fit$sig2_u,
			sv = fit$sig2_v,
			suv0 = cx$suv0,
			kappa0 = cx$kappa0
		)
	}
	pb <- .ae_mf_prior_blocks(S$prior, R)
	cv <- .ae_mf_cov0(S$U, S$V, C0, pb$uu, pb$vv)
	S$Sig <- cv$Sig; S$Om <- cv$Om
	if (cx$pln) .ae_mf_reseed(cx, S, fit) else S
}

# mean of the multiplicative term (NULL at rank 0) and the variance of every pair's product
.ae_mf_product <- function(cx, S) {
	if (cx$R == 0L) return(NULL)
	if (cx$sym) tcrossprod(sweep(S$F, 2L, S$sg, "*"), S$F) else tcrossprod(S$U, S$V)
}
.ae_mf_vprod <- function(cx, S) {
	if (cx$R == 0L) return(0)
	if (cx$sym) .ae_mf_prod_var_sym(S$F, S$sg, S$Sig) else .ae_mf_prod_var(S$U, S$V, S$Sig, S$Om)
}

# the mean and the variance of the predictor of every cell at a state, and the mean of the product
.ae_mf_predictor <- function(cx, S) {
	Pm <- .ae_mf_product(cx, S)
	lin <- S$mu + outer(S$a, S$b, "+")
	if (cx$R > 0L) lin <- lin + Pm
	eta <- array(lin, c(cx$nr, cx$nc, cx$Tt))
	if (!is.null(cx$xmat)) eta <- eta + as.vector(cx$xmat %*% S$beta)
	Vt <- S$vact
	if (cx$R > 0L) Vt <- Vt + .ae_mf_vprod(cx, S)
	list(Pm = Pm, eta = eta, Vt = Vt)
}

# quadrature at a state: the expected working problem of every cell, or of every dyad's pair of cells
# under a dyadic correlation (R/ame_als_pair.R)
.ae_mf_evaluate <- function(cx, S, pr = NULL) {
	# the predictor depends only on the means and covariances, so a normal pair fit reuses it at another correlation
	if (is.null(pr)) pr <- .ae_mf_predictor(cx, S)
	V3 <- array(pr$Vt, c(cx$nr, cx$nc, cx$Tt))
	if (is.null(cx$pair)) {
		wk <- .ae_mf_cells(cx$y, pr$eta, V3, cx$link, cx$quad, ve = S$ve, kin = cx$kin)
		return(list(wk = wk, Pm = pr$Pm, ell = sum(wk$ell), eta = pr$eta, V = pr$Vt, ve = S$ve))
	}
	pe <- if (cx$gauss) .ae_gauss_pair_cells(cx, S, pr$eta, V3, S$rho) else .ae_pair_cells(cx, S, pr$eta, V3, S$rho)
	list(
		wk = pe$wk,
		Pm = pr$Pm,
		ell = sum(pe$wk$ell),
		eta = pr$eta,
		V = pr$Vt,
		pair = pe$pair,
		ve = S$ve,
		pr = if (cx$gauss) pr
	)
}

# objective the iteration climbs: the expected log likelihood, minus the KL of the effects at the reference
# covariance `sab_ref`, plus the factor part with the scales integrated out
.ae_mf_merit <- function(cx, S, ev, sab_ref) {
	Lam_s <- .ae_eb_precision(sab_ref, cx$ab)
	d <- c(S$a, if (!cx$sym) S$b)
	if (!is.null(cx$M)) d <- d - as.vector(cx$M %*% S$gamma)
	kl_s <- 0.5 * (sum(Lam_s * S$Sc) + sum(d * as.vector(Lam_s %*% d)) -
			.ae_mf_logdet(S$Sc))
	cx$nu * ev$ell - kl_s + (if (cx$R > 0L) .ae_mf_factor_term(S, cx$R) else 0) +
		(if (cx$has_ve) .ae_pln_ve_logprior(S$ve, cx$ab$ve_prior) else 0)
}

# evidence lower bound at the plug-in scales (hyperpriors left out): expected log likelihood less the KL
# of the effects and factors from their priors
.ae_mf_elbo <- function(cx, S, ev) {
	Lam_s <- .ae_eb_precision(S$sab, cx$ab)
	mu0 <- if (!is.null(cx$M)) as.vector(cx$M %*% S$gamma) else 0
	kl_a <- .ae_mf_kl_actors(S$Sc, c(S$a, if (!cx$sym) S$b), mu0, Lam_s)
	kl_f <- if (cx$R > 0L) .ae_mf_kl_factors(S, cx$R) else 0
	ell <- cx$nu * ev$ell
	list(elbo = ell - kl_a - kl_f, ell = ell, kl_actors = kl_a, kl_factors = kl_f)
}


# the conditional optimum of every block at the weights and response of the
# state's quadrature: the actor block (coefficients, effects, covariance of the
# effects, its own covariance update), then the factor sweeps at fixed scales
.ae_mf_target <- function(cx, S, ev, inner) {
	ab <- cx$ab; R <- cx$R
	Warr <- ev$wk$w
	wsum <- rowSums(Warr, dims = 2L)
	# under a dyadic correlation the additive block solves the pair normal equations (vector G and the cross
	# weights); the factor rows use the cell response, whose score is the exact pair score
	pair <- ev$pair
	Zst <- if (is.null(pair)) ev$wk$z else ifelse(Warr > 0, pair$G / Warr, 0)
	st <- ae_bcd_stats_cpp(Warr, Zst, cx$X, cx$nr, cx$nc, cx$Tt, cx$p, .ae_threads())
	ctx <- .ae_eb_context(ab, st, wsum, cx$p)
	cross <- NULL
	if (!is.null(pair)) {
		cross <- .ae_pair_cross(cx, pair$cw)
		ctx <- .ae_pair_add_cross(ctx, cross)
	}
	stp <- .ae_mf_actor_step(ctx, ab, ev$Pm, st, S$sab, cross)
	tg <- list(
		theta = stp$theta,
		gamma = stp$gamma,
		a = stp$a,
		b = stp$b,
		Sc = stp$Sc,
		sab = stp$sab_new
	)
	# EM image of the noise variance of a count or normal cell, from the expected second moment of the noise
	if (cx$has_ve) tg$ve <- if (isTRUE(cx$ve_fixed)) S$ve else
		.ae_pln_ve_image(sum(ev$wk$ee), cx$n_obs, cx$nu, cx$ab$ve_prior, .ae_ve_floor(cx$ab))
	if (R > 0L) {
		phi <- matrix(st$gm %*% stp$theta, cx$nr, cx$nc)
		sz_f <- if (is.null(pair)) st$sz else rowSums(Warr * ev$wk$z, dims = 2L)
		Rbar <- ae_bcd_rbar_cpp(stp$a, stp$b, phi, 0, sz_f, wsum)
		if (cx$sym) {
			sw <- .ae_mf_sweeps_sym(Rbar, wsum, S$F, S$sg, S$Sig, S$prior, inner)
			tg$F <- sw$F; tg$Sig <- sw$Sig
			tg$Pm <- tcrossprod(sweep(tg$F, 2L, S$sg, "*"), tg$F)
		} else {
			sw <- .ae_mf_sweeps(Rbar, wsum, S$U, S$V, S$Sig, S$Om, S$prior, inner)
			tg$U <- sw$U; tg$V <- sw$V; tg$Sig <- sw$Sig; tg$Om <- sw$Om
			tg$Pm <- tcrossprod(tg$U, tg$V)
		}
	}
	tg
}

# the state a fraction `rho` of the way from S to the target
.ae_mf_blend <- function(cx, S, tg, rho) {
	mix <- function(x, y) x + rho * (y - x)
	B <- S
	B$mu <- mix(S$mu, tg$theta[1L]); B$beta <- mix(S$beta, tg$theta[-1L])
	B$a <- mix(S$a, tg$a); B$b <- mix(S$b, tg$b); B$gamma <- mix(S$gamma, tg$gamma)
	B$Sc <- mix(S$Sc, tg$Sc)
	if (cx$has_ve) B$ve <- mix(S$ve, tg$ve)
	B$vact <- .ae_mf_vact(B$Sc, cx$ab)
	if (cx$R > 0L) {
		if (cx$sym) {
			B$F <- mix(S$F, tg$F)
		} else {
			B$U <- mix(S$U, tg$U); B$V <- mix(S$V, tg$V)
			B$Om <- mix(S$Om, tg$Om)
		}
		B$Sig <- mix(S$Sig, tg$Sig)
	}
	B
}

# EM update of the factor prior scales and the effect covariance at the state reached; returns the state
# and the largest relative change of a scale
.ae_mf_hyper <- function(cx, S, tg, rho) {
	sab_old <- S$sab
	S$sab <- S$sab + rho * (tg$sab - S$sab)
	change <- max(abs(log(diag(S$sab) / diag(sab_old))))
	if (cx$R > 0L) {
		old <- S$prior
		S$prior <- .ae_mf_prior_update(S, cx$R)
		scales <- function(pr) switch(
			pr$kind,
			joint = diag(pr$Suv),
			scalar = pr$tau,
			diag = c(pr$su, pr$sv)
		)
		change <- max(change, max(abs(log(scales(S$prior) / scales(old)))))
	}
	list(S = S, change = change)
}


# gauge of the factors: U -> U A, V -> V A^-T leaves the likelihood unchanged but not the prior and entropy
# terms. with the scales integrated out the objective along A is, for per-column scales,
#   f(A) = -(df + n_u) / 2 sum_k log(psi + [A' E_u A]_kk)
#          -(df + n_v) / 2 sum_k log(psi + [A^-1 E_v A^-T]_kk) + (n_u - n_v) log |det A|,
# and -(df + n) / 2 log |Psi + D SE D'| for the joint prior (D = diag(A', A^-1)). EM approaches its maximum
# very slowly, so it is taken directly; this does not change the fixed point
.ae_mf_gauge <- function(cx, S) {
	R <- cx$R
	if (R == 0L || cx$sym) return(S)
	pr <- S$prior
	Eu <- crossprod(S$U) + apply(S$Sig, c(1L, 2L), sum)
	Ev <- crossprod(S$V) + apply(S$Om, c(1L, 2L), sum)
	nu <- nrow(S$U); nv <- nrow(S$V)
	joint <- identical(pr$kind, "joint")
	if (joint) {
		SE <- crossprod(cbind(S$U, S$V))
		iu <- seq_len(R); iv <- R + seq_len(R)
		SE[iu, iu] <- Eu; SE[iv, iv] <- Ev
		SE <- (SE + t(SE)) / 2
	}
	fn <- function(A) {
		B <- tryCatch(solve(A), error = function(e) NULL)
		if (is.null(B)) return(-Inf)
		if (joint) {
			D <- matrix(0, 2L * R, 2L * R)
			D[iu, iu] <- t(A); D[iv, iv] <- B
			return(-0.5 * (pr$kappa0 + nu) * .ae_mf_logdet(pr$psi + D %*% SE %*% t(D)))
		}
		du <- pr$suv0 + colSums(A * (Eu %*% A))
		dv <- pr$suv0 + rowSums((B %*% Ev) * B)
		-0.5 * (pr$kappa0 + nu) * sum(log(du)) - 0.5 * (pr$kappa0 + nv) * sum(log(dv)) +
			(nu - nv) * log(abs(det(A)))
	}
	gr <- function(A) {
		B <- solve(A)
		if (joint) {
			D <- matrix(0, 2L * R, 2L * R)
			D[iu, iu] <- t(A); D[iv, iv] <- B
			K <- SE %*% t(D) %*% .ae_safe_inverse(pr$psi + D %*% SE %*% t(D))
			return(-(pr$kappa0 + nu) * (K[iu, iu] - t(B %*% K[iv, iv] %*% B)))
		}
		du <- pr$suv0 + colSums(A * (Eu %*% A))
		dv <- pr$suv0 + rowSums((B %*% Ev) * B)
		-(pr$kappa0 + nu) * sweep(Eu %*% A, 2L, du, "/") +
			(pr$kappa0 + nv) * t(B) %*% sweep(B %*% Ev, 1L, dv, "/") %*% t(B) +
			(nu - nv) * t(B)
	}
	I <- diag(R)
	f0 <- fn(I)
	opt <- tryCatch(
		stats::optim(
			as.vector(I),
			function(x) -fn(matrix(x, R, R)),
			function(x) -as.vector(gr(matrix(x, R, R))),
			method = "BFGS",
			control = list(maxit = 200L, reltol = 1e-14)
		),
		error = function(e) NULL
	)
	if (is.null(opt) || !all(is.finite(opt$par))) return(S)
	A <- matrix(opt$par, R, R)
	if (!(-opt$value > f0 + 1e-12) || abs(det(A)) < 1e-8 || abs(det(A)) > 1e8) return(S)
	B <- solve(A)
	S$U <- S$U %*% A; S$V <- S$V %*% t(B)
	for (i in seq_len(nu)) S$Sig[, , i] <- t(A) %*% S$Sig[, , i] %*% A
	for (j in seq_len(nv)) S$Om[, , j] <- B %*% S$Om[, , j] %*% t(B)
	S
}


# covariance of a factor row from the information its own cells carry through the other side's means,
# (Lam + sum_j C_ij vhat_j vhat_j')^-1 with cavity weights C (.ae_mf_info_w); the mean-field covariance counts
# the predictor variance twice and understates the rows' second moments
.ae_mf_cavity_cov <- function(cx, S, C) {
	R <- cx$R
	pb <- .ae_mf_prior_blocks(S$prior, R)
	if (cx$sym) {
		Cn <- C; diag(Cn) <- 0
		Fs <- sweep(S$F, 2L, S$sg, "*")
		return(list(
			Sig = .ae_mf_row_solve(
				Cn %*% .ae_mf_outer_rows(Fs),
				pb$uu,
				matrix(0, nrow(Fs), R)
			)$Sig
		))
	}
	list(
		Sig = .ae_mf_row_solve(
			C %*% .ae_mf_outer_rows(S$V),
			pb$uu,
			matrix(0, nrow(S$U), R)
		)$Sig,
		Om = .ae_mf_row_solve(
			crossprod(C, .ae_mf_outer_rows(S$U)),
			pb$vv,
			matrix(0, nrow(S$V), R)
		)$Sig
	)
}

# the isotropic prior of the stacked rows of a directed network: one scale s for every coordinate and the
# correlation rho of a sender and a receiver coordinate of one actor
.ae_mf_prior_iso <- function(prior, R, s, rho) {
	Suv <- s * rbind(cbind(diag(R), rho * diag(R)), cbind(rho * diag(R), diag(R)))
	prior$Suv <- Suv
	prior$Lam <- .ae_safe_inverse(Suv)
	prior
}

# undirected sampler prior: interaction U L U' with L_r ~ N(0, 1) and column precisions psi_r ~ Gamma(1, rate 1/2).
# the likelihood sees f_r = sqrt(|L_r|) u_r with rows N(0, tau_r), tau_r = |L_r| / psi_r, whose prior is
# p(tau) = tau^-2 I(1 / (2 tau)), I(c) = (2 pi)^-1/2 - c exp(c^2 / 2) Phi(-c) (heavy tailed); I(c) uses its
# asymptotic series where the difference cancels
.ae_mf_sym_I <- function(cc) {
	out <- numeric(length(cc))
	big <- cc > 12
	if (any(big)) {
		c2 <- cc[big]^-2
		out[big] <- (2 * pi)^-0.5 * c2 * (1 - 3 * c2 + 15 * c2^2 - 105 * c2^3)
	}
	if (any(!big)) out[!big] <- (2 * pi)^-0.5 - cc[!big] * exp(cc[!big]^2 / 2 + stats::pnorm(-cc[!big], log.p = TRUE))
	pmax(out, 0)
}
.ae_mf_sym_logprior <- function(tau) -2 * log(tau) + log(.ae_mf_sym_I(1 / (2 * tau)))

# plug-in scale of an undirected column from the sum of squares m of its n rows: the sampler's conditional
# of tau_r on a log grid (1e-6 to 1e6 in units of `unit`), summarised by 1 / E[1 / tau]
.ae_mf_sym_tau_plug <- function(m, n, stat = .ae_mf_scale_stat(), unit = 1) {
	lt <- seq(log(1e-6 * unit), log(1e6 * unit), length.out = 2001L)
	tau <- exp(lt)
	# the weights are N(0, unit^2) for a normal fit, so the prior is the density above of tau / unit (R/normal_scale.R)
	lp0 <- .ae_mf_sym_logprior(tau / unit)
	vapply(m, function(mm) {
		lp <- lp0 - 0.5 * n * log(tau) - mm / (2 * tau) + lt
		w <- exp(lp - max(lp))
		if (identical(stat, "mean")) sum(w * tau) / sum(w) else sum(w) / sum(w / tau)
	}, 0)
}

# bipartite sampler prior: P = U G V' with standard normal U, V and G, so the sender side has covariance
# G'G ~ Wishart(R, I) and the sender scales are its eigenvalues, with joint density
#   prod_k lambda_k^-1/2 exp(-lambda_k / 2) prod_{k < l} |lambda_k - lambda_l|.
# given the sum of squares m_k of column k and the other scales `s`,
#   p(lambda | m_k) ~ lambda^(-(n + 1) / 2) exp(-lambda / 2 - m_k / (2 lambda)) prod_{l != k} |lambda - s_l|,
# summarised by E[1 / lambda]^-1 on a log grid
.ae_mf_bip_scale_plug <- function(m, n, s, stat = .ae_mf_scale_stat(), unit = 1) {
	# `unit` is the prior variance of an entry of G: the response variance for a normal fit, else 1
	lt <- seq(log(1e-6 * unit), log(1e6 * unit), length.out = 2001L)
	tau <- exp(lt)
	R <- length(m)
	vapply(seq_len(R), function(k) {
		lp <- 0.5 * (1 - n) * lt - tau / (2 * unit) - m[k] / (2 * tau)
		for (l in setdiff(seq_len(R), k)) lp <- lp + log(abs(tau - s[l]))
		w <- exp(lp - max(lp))
		if (identical(stat, "mean")) sum(w * tau) / sum(w) else sum(w) / sum(w / tau)
	}, 0)
}

# posterior statistic of a scale that sets the row prior: "plugin" (default) the inverse of the mean
# precision, as in the mean-field update; "mean" the posterior mean covariance, dividing by df + n - 2 R - 1.
# options(lame.eb_scale_stat = "mean") selects it
.ae_mf_scale_stat <- function() {
	stat <- getOption("lame.eb_scale_stat")
	if (is.null(stat)) return("plugin")
	if (!is.character(stat) || length(stat) != 1L || !stat %in% c("mean", "plugin")) {
		cli::cli_abort("{.code options(lame.eb_scale_stat)} must be {.val mean} or {.val plugin}.")
	}
	stat
}

# covariance of the stacked rows from their cavity second moments SE: the variational posterior
# IW(Psi + SE, df + n) has mean (Psi + SE) / (df + n - 2 R - 1) and mean-field update (Psi + SE) / (df + n)
# (.ae_mf_scale_stat), reduced to an isotropic scale and one sender-receiver correlation
.ae_mf_scale_target <- function(cx, S, C) {
	R <- cx$R
	pr <- S$prior
	rest <- .ae_mf_cavity_cov(cx, S, C)
	if (cx$sym) {
		S2 <- S
		S2$Sig <- rest$Sig
		sm <- .ae_mf_second_moments(S2, R)
		return(
			list(tau = .ae_mf_sym_tau_plug(sm$col, sm$rows, unit = if (isTRUE(cx$gauss)) cx$ysd else 1), col = sm$col)
		)
	}
	if (identical(pr$kind, "diag")) {
		S2 <- S
		S2$Sig <- rest$Sig; S2$Om <- rest$Om
		sm <- .ae_mf_second_moments(S2, R)
		return(
			list(
				su = .ae_mf_bip_scale_plug(sm$su, sm$nu, pr$su, unit = if (isTRUE(cx$gauss)) cx$ysd^2 else 1),
				col = sm$su
			)
		)
	}
	S2 <- S
	S2$Sig <- rest$Sig; S2$Om <- rest$Om
	sm <- .ae_mf_second_moments(S2, R)
	# keep at least one pseudo-observation beyond the rows in the posterior-mean divisor, so a supplied
	# kappa0 <= 2 R + 1 cannot make the update expand
	stat <- if (!is.null(cx$scale_stat)) cx$scale_stat else .ae_mf_scale_stat()
	den <- if (identical(stat, "mean")) sm$nu + max(pr$kappa0 - (2 * R + 1), 1) else sm$nu + pr$kappa0
	M <- (pr$psi + sm$SE) / max(den, 1)
	iu <- seq_len(R); iv <- R + seq_len(R)
	s <- mean(diag(M))
	rho <- mean(diag(M[iu, iv, drop = FALSE])) / s
	list(s = s, rho = max(min(rho, 0.9), -0.9), M = M)
}

# the conditional objective at the prior the state carries: the same as .ae_mf_merit without integrating
# the scale out, so two states are compared at one prior
.ae_mf_factor_term_fixed <- function(S, R) {
	pr <- S$prior
	sm <- .ae_mf_second_moments(S, R)
	if (!is.null(S$F)) {
		return(
			-0.5 * (sum(sm$col / pr$tau) + sm$rows * sum(log(pr$tau))) +
				0.5 * sum(.ae_mf_logdet_rows(S$Sig))
		)
	}
	ent <- 0.5 * (sum(.ae_mf_logdet_rows(S$Sig)) + sum(.ae_mf_logdet_rows(S$Om)))
	if (!identical(pr$kind, "joint")) {
		return(
			-0.5 * (sum(sm$su / pr$su) + sm$nu * sum(log(pr$su)) + sum(sm$sv / pr$sv) + sm$nv * sum(log(pr$sv))) + ent
		)
	}
	-0.5 * (sum(pr$Lam * sm$SE) + sm$nu * .ae_mf_logdet(pr$Suv)) + ent
}
.ae_mf_merit_fixed <- function(cx, S, ev, sab_ref) {
	Lam_s <- .ae_eb_precision(sab_ref, cx$ab)
	d <- c(S$a, if (!cx$sym) S$b)
	if (!is.null(cx$M)) d <- d - as.vector(cx$M %*% S$gamma)
	kl_s <- 0.5 * (sum(Lam_s * S$Sc) + sum(d * as.vector(Lam_s %*% d)) -
			.ae_mf_logdet(S$Sc))
	cx$nu * ev$ell - kl_s + .ae_mf_factor_term_fixed(S, cx$R) +
		(if (cx$has_ve) .ae_pln_ve_logprior(S$ve, cx$ab$ve_prior) else 0)
}

# one bounded step of a log scale x towards its plug-in value (residual res = log(target) - x): a Newton
# step with the slope from the previous iterate when that is a contraction's, else the relaxed step. the
# plug-in map converges slowly where the data are weak; the fixed point is unchanged
.ae_mf_scale_step <- function(x, res, acc, key, omega, smax) {
	prev <- acc[[key]]
	step <- omega * res
	if (!is.null(prev) && abs(x - prev$x) > 1e-10) {
		slope <- (res - prev$res) / (x - prev$x)
		if (is.finite(slope) && slope < -1e-3 && slope > -2) step <- -res / slope
	}
	acc[[key]] <- list(x = x, res = res)
	max(min(step, 10 * smax), -10 * smax)
}

# cavity-rule hyperparameter update: an EM step for the effect covariance, and bounded relaxed steps in
# the log for the scale and sender-receiver correlation towards their plug-in values (.ae_mf_scale_target)
.ae_mf_hyper_cavity <- function(cx, S, tg, rho, C, update = TRUE, acc = new.env(), omega = 0.5, smax = 0.05) {
	sab_old <- S$sab
	S$sab <- S$sab + rho * (tg$sab - S$sab)
	change <- max(abs(log(diag(S$sab) / diag(sab_old))))
	if (update && !cx$sym && identical(S$prior$kind, "diag")) {
		tgt <- .ae_mf_scale_target(cx, S, C)
		ld <- vapply(
			seq_len(cx$R),
			function(k) .ae_mf_scale_step(
				log(S$prior$su[k]),
				log(tgt$su[k]) - log(S$prior$su[k]),
				acc,
				paste0("su", k),
				omega,
				smax
			),
			0
		)
		S$prior$su <- S$prior$su * exp(ld)
		change <- max(change, abs(ld))
	} else if (update && cx$sym) {
		tgt <- .ae_mf_scale_target(cx, S, C)
		ld <- vapply(
			seq_len(cx$R),
			function(k) .ae_mf_scale_step(
				log(S$prior$tau[k]),
				log(tgt$tau[k]) - log(S$prior$tau[k]),
				acc,
				paste0("tau", k),
				omega,
				smax
			),
			0
		)
		S$prior$tau <- S$prior$tau * exp(ld)
		change <- max(change, abs(ld))
	} else if (update) {
		tgt <- .ae_mf_scale_target(cx, S, C)
		s_old <- S$prior$Suv[1L, 1L]
		rho_old <- S$prior$Suv[1L, cx$R + 1L] / s_old
		ds <- .ae_mf_scale_step(log(s_old), log(tgt$s) - log(s_old), acc, "s", omega, smax)
		s_new <- s_old * exp(ds)
		r_new <- rho_old + max(min(omega * (tgt$rho - rho_old), smax), -smax)
		S$prior <- .ae_mf_prior_iso(S$prior, cx$R, s_new, r_new)
		change <- max(change, abs(ds), abs(r_new - rho_old))
	}
	list(S = S, change = change)
}

# all the scales of a prior multiplied by `f` (the covariance of the stacked rows, the scale of the
# undirected columns, the sender scales of a bipartite factor: its receiver scale is fixed)
.ae_mf_prior_shift <- function(prior, f) {
	if (identical(prior$kind, "joint")) {
		prior$Suv <- prior$Suv * f
		prior$Lam <- .ae_safe_inverse(prior$Suv)
	} else if (identical(prior$kind, "columns")) {
		prior$tau <- prior$tau * f
	} else if (identical(prior$kind, "diag")) {
		prior$su <- prior$su * f
	}
	prior
}

# the derivative of the log marginal likelihood in the log of the common scale: half the
# difference between the second moments of the rows (restored covariances) and the number of
# coordinates the prior accounts for, d log p(W | f S) / d log f = (tr(S^-1 SE) - n p) / 2 at f = 1
.ae_mf_scale_score <- function(cx, S, C) {
	R <- cx$R
	rest <- .ae_mf_cavity_cov(cx, S, C)
	S2 <- S
	S2$Sig <- rest$Sig; if (!cx$sym) S2$Om <- rest$Om
	sm <- .ae_mf_second_moments(S2, R)
	pr <- S$prior
	if (identical(pr$kind, "joint")) return(0.5 * (sum(pr$Lam * sm$SE) - sm$nu * 2L * R))
	if (identical(pr$kind, "columns")) return(0.5 * (sum(sm$col / pr$tau) - sm$rows * R))
	0.5 * (sum(sm$su / pr$su) - sm$nu * R)
}

# prior curvature in the log of a common scale multiplier: p / (2 s) along the isotropic ray of the
# directed inverse-Wishart, sum(s) / 2 for the bipartite Wishart eigenvalues, by differences for the undirected prior
.ae_mf_scale_prior_curvature <- function(cx, S) {
	pr <- S$prior
	R <- cx$R
	if (identical(pr$kind, "joint")) return(2L * R / (2 * mean(diag(pr$Suv))))
	if (identical(pr$kind, "diag")) return(sum(pr$su) / 2)
	h <- 0.05
	lp <- function(x) sum(.ae_mf_sym_logprior(pr$tau * exp(x)) + log(pr$tau * exp(x)))
	-(lp(h) - 2 * lp(0) + lp(-h)) / h^2
}


# refine a converged EB mode (`fit` from .ae_irls_run, `W0` its Fisher weights) to the variational
# posterior of its factors; every parameter moves. each iteration steps a fraction of the way to the
# conditional optimum of every block (.ae_mf_target), halving it until the objective does not fall; the fixed
# point does not depend on the fraction. `pair` configures the dyadic correlation of a directed probit fit:
# `rho` its start, `quad` nodes a side (5), `update` whether it is estimated, `coarse` whether to start with
# three nodes. the prior of rho is the package prior (R/rho_prior.R)
.ae_mf_fit <- function(
	prep,
	fit,
	ab,
	R,
	symmetric,
	link,
	W0,
	max_iter = 400L,
	tol = 1e-6,
	inner = 3L,
	quad = 9L,
	verbose = FALSE,
	hyper = NULL,
	step_max = 1,
	start = NULL,
	scale = "plugin",
	pair = NULL,
	ve_fixed = FALSE,
	scale_stat = NULL
) {
	cx <- .ae_mf_context(prep, ab, R, link, quad, hyper, pair)
	# posterior statistic of the scale (.ae_mf_scale_stat) when the caller sets one
	cx$scale_stat <- scale_stat
	# hold a count network's noise variance at its start (the Laplace maximiser, .ae_pln_ve_search) rather
	# than taking its EM step
	cx$ve_fixed <- isTRUE(ve_fixed)
	S <- if (is.null(start)) .ae_mf_init(cx, fit, W0) else start
	# the inner rule of a count cell is fixed for the iteration, from the variance of the noise with a margin
	if (cx$pln) cx$kin <- .ae_pln_kin(1.5 * S$ve)
	# the dyadic correlation starts at the estimate at the marginal predictors of the start
	if (!is.null(cx$pair) && is.null(S$rho)) S$rho <- if (cx$gauss) .ae_gauss_rho_start(cx, S) else .ae_pair_rho_start(cx, S)
	# scale rule: "plugin" the EM update from the mean-field second moments, "cavity" the same from the cavity
	# moments (isotropic with one sender-receiver correlation for a directed prior), "fixed" keeps the state's prior
	cavity <- R > 0L && scale %in% c("cavity", "fixed") &&
		((!cx$sym && S$prior$kind %in% c("joint", "diag")) || (cx$sym && S$prior$kind %in% c("scalar", "columns")))
	if (cavity && identical(scale, "cavity") && is.null(start)) {
		if (cx$sym) {
			S$prior <- list(kind = "columns", tau = pmax(colSums(S$F^2) / nrow(S$F), 1e-3 * cx$ysd))
		} else if (identical(S$prior$kind, "diag")) {
			# the receiver factor has the fixed unit scale: the mode's factors are rescaled to it
			# (U -> U D, V -> V D^-1 leaves the product alone) and the sender scales start at the second moments
			dv <- sqrt(pmax(colSums(S[["V"]]^2) / nrow(S[["V"]]), 1e-8))
			S$U <- sweep(S$U, 2L, dv, "*"); S[["V"]] <- sweep(S[["V"]], 2L, dv, "/")
			for (k in seq_len(R)) { S$Sig[k, , ] <- S$Sig[k, , ] * dv[k]; S$Sig[, k, ] <- S$Sig[, k, ] * dv[k]
				S$Om[k, , ] <- S$Om[k, , ] / dv[k]; S$Om[, k, ] <- S$Om[, k, ] / dv[k] }
			S$prior <- list(
				kind = "diag",
				su = pmax(colSums(S$U^2) / nrow(S$U), 1e-3 * cx$ysd^2),
				sv = rep(1, R),
				suv0 = S$prior$suv0,
				kappa0 = S$prior$kappa0
			)
		} else {
			s0 <- (sum(fit$U^2) + sum(fit[["V"]]^2)) / (2 * nrow(fit$U) * R)
			S$prior <- .ae_mf_prior_iso(S$prior, R, max(s0, 1e-3 * cx$suv_joint), 0)
		}
	}
	if (R > 0L && !cx$sym && !cavity) {
		S <- .ae_mf_gauge(cx, S)
		S$prior <- .ae_mf_prior_update(S, R)
	}
	merit_of <- if (cavity) .ae_mf_merit_fixed else .ae_mf_merit
	acc <- new.env(parent = emptyenv())
	ve_hist <- numeric(0)
	# the pair quadrature dominates the cost, so the iteration converges with the three-node rule first and
	# then with the fit's rule, which has the same fixed point
	quad_final <- cx$pair_quad
	coarse <- !is.null(cx$pair) && is.null(start) && quad_final > 3L && !isFALSE(pair$coarse)
	if (coarse) cx$pair_quad <- 3L
	ev <- .ae_mf_evaluate(cx, S)
	ev_rho <- S$rho
	merit <- merit_of(cx, S, ev, S$sab)
	rho <- step_max
	rho_min <- 2^-7
	converged <- FALSE
	termination <- "iteration limit"
	trace <- numeric(0)
	# the pair rules the iteration was refined to after a stall (.ae_pair_quad_max)
	refined <- integer(0)
	res <- Inf
	gain <- Inf
	it <- 0L
	for (it in seq_len(max_iter)) {
		# a scale that moves between iterations changes the objective, so the state and its candidates
		# are compared at the state's current prior
		if (cavity || !is.null(cx$pair)) merit <- merit_of(cx, S, ev, S$sab)
		# when the correlation moved slightly, correct the merit by the quadratic model of its score and curvature
		if (!is.null(cx$pair) && !identical(ev_rho, S$rho)) {
			dl <- S$rho - ev_rho
			merit <- merit + cx$nu * (ev$pair$rs * dl + 0.5 * ev$pair$rc * dl^2)
		}
		tg <- .ae_mf_target(cx, S, ev, inner)
		res_main <- max(abs(tg$theta - c(S$mu, S$beta)), abs(tg$a - S$a), abs(tg$b - S$b))
		if (R > 0L) res_main <- max(res_main, abs(tg$Pm - ev$Pm)[cx$obs_dyad])
		res_main <- res_main / cx$ysd
		# the noise variance of a count cell moves with the other parameters; this is the log change of a full EM step
		dve <- if (cx$has_ve && !cx$ve_fixed) abs(log(tg$ve / S$ve)) else 0
		res <- res_main
		accepted <- FALSE
		rho_try <- min(step_max, 2 * rho)
		while (rho_try >= rho_min) {
			cand <- .ae_mf_blend(cx, S, tg, rho_try)
			evc <- .ae_mf_evaluate(cx, cand)
			mc <- merit_of(cx, cand, evc, S$sab)
			if (is.finite(mc) && mc >= merit - 1e-9 * (1 + abs(merit))) {
				accepted <- TRUE
				break
			}
			rho_try <- rho_try / 2
		}
		if (!accepted) {
			if (coarse) {
				# the coarse rule cannot improve its own objective further: continue with the fit's rule
				coarse <- FALSE
				cx$pair_quad <- quad_final
				ev <- .ae_mf_evaluate(cx, S)
				ev_rho <- S$rho
				merit <- merit_of(cx, S, ev, S$sab)
				next
			}
			# a pair rule resolves the objective only to the accuracy of its nodes, which can stall a weakly informed
			# fit near the fixed point; refine to seven, then nine nodes before stopping (a normal pair has no quadrature)
			if (!is.null(cx$pair) && cx$pair_quad > 0L && cx$pair_quad < .ae_pair_quad_max) {
				cx$pair_quad <- if (cx$pair_quad < 7L) 7L else .ae_pair_quad_max
				quad_final <- cx$pair_quad
				refined <- c(refined, cx$pair_quad)
				ev <- .ae_mf_evaluate(cx, S)
				ev_rho <- S$rho
				merit <- merit_of(cx, S, ev, S$sab)
				rho <- step_max
				next
			}
			termination <- "no improving step"
			break
		}
		rho <- rho_try
		gain <- mc - merit
		if (cavity) {
			hy <- .ae_mf_hyper_cavity(
				cx,
				cand,
				tg,
				rho,
				.ae_mf_info_sum(cx, evc),
				update = identical(scale, "cavity"),
				acc = acc
			)
		} else {
			cand <- .ae_mf_gauge(cx, cand)
			hy <- .ae_mf_hyper(cx, cand, tg, rho)
		}
		S <- hy$S; ev <- evc
		hy$change <- max(hy$change, dve)
		ev_rho <- S$rho
		# the noise variance of a count network converges slowly under EM, so jump along the geometric tail of
		# the accepted iterates (.ae_pln_aitken)
		if (cx$has_ve && !cx$ve_fixed) {
			ve_hist <- c(ve_hist, log(S$ve))
			jump <- .ae_pln_aitken(ve_hist)
			if (is.finite(jump)) {
				S$ve <- exp(log(S$ve) + jump)
				ev <- .ae_mf_evaluate(cx, S)
				ve_hist <- numeric(0)
			}
		}
		# scoring step for the dyadic correlation at the accepted state; re-evaluate there unless the step is tiny
		if (!is.null(cx$pair) && cx$pair_update) {
			rho_new <- if (cx$gauss) max(min(evc$pair$rho_next, 0.98), -0.98) else .ae_pair_rho_next(evc$pair, S$rho)
			drho <- abs(rho_new - S$rho)
			if (drho > 0) {
				S$rho <- rho_new
				# a normal pair is cheap, so it is always re-evaluated, keeping the noise variance and correlation in step
				if (cx$gauss || abs(S$rho - ev_rho) > 5e-4) {
					ev <- .ae_mf_evaluate(cx, S, pr = if (cx$gauss) evc$pr)
					ev_rho <- S$rho
				}
			}
			hy$change <- max(hy$change, drho)
		}
		merit <- merit_of(cx, S, ev, S$sab)
		trace <- c(trace, merit)
		if (verbose && it %% 10L == 0L) {
			cli::cli_text("variational iteration {it}: residual {signif(res_main, 3)}, step {rho}, objective {round(merit, 4)}")
		}
		# converged: a full step moves nothing by more than tol or, along a flat direction, by more than 100 tol
		# while raising the objective by under 1e-8 of its size
		flat <- res_main < 100 * tol && gain < 1e-8 * max(1, abs(merit))
		# the cavity scale is held to a tenth of the other hyperparameters' change tolerance
		conv_now <- (res_main < tol || flat) && hy$change < (if (cavity) 10 else 100) * tol
		if (coarse && (conv_now || (res_main < 3e-3 && hy$change < 3e-3))) {
			# the coarse rule is as close to its fixed point as to the fit's: continue with the fit's rule
			coarse <- FALSE
			cx$pair_quad <- quad_final
			ev <- .ae_mf_evaluate(cx, S)
			ev_rho <- S$rho
			merit <- merit_of(cx, S, ev, S$sab)
			next
		}
		if (conv_now) {
			converged <- TRUE
			termination <- "converged"
			break
		}
	}
	# the evaluation the fit reports is at the rule and the correlation of the state
	if (!is.null(cx$pair) && (coarse || !identical(ev_rho, S$rho))) {
		cx$pair_quad <- quad_final
		ev <- .ae_mf_evaluate(cx, S)
	}
	out <- list(
		S = S,
		ev = ev,
		cx = cx,
		converged = converged,
		termination = termination,
		iterations = it,
		residual = res,
		gain = gain,
		tol = tol,
		merit = merit,
		trace = trace,
		step = rho,
		elbo = .ae_mf_elbo(cx, S, ev),
		V = ev$V
	)
	if (!is.null(cx$pair)) {
		out$pair <- list(
			rho = S$rho,
			score = ev$pair$rs,
			curvature = -ev$pair$rc,
			sd = if (cx$gauss) ev$pair$rho_sd else .ae_pair_rho_sd(ev$pair, S$rho),
			n_pairs = cx$pair$n_pairs,
			quad = cx$pair_quad,
			refined = refined
		)
	}
	if (cavity) {
		C <- .ae_mf_info_sum(cx, ev)
		out$scale <- if (identical(S$prior$kind, "diag")) {
			list(
				rule = scale,
				s = mean(S$prior$su),
				rho = 0,
				su = S$prior$su,
				target = .ae_mf_scale_target(cx, S, C),
				C = C
			)
		} else if (cx$sym) {
			list(
				rule = scale,
				s = mean(S$prior$tau),
				rho = 0,
				tau = S$prior$tau,
				target = .ae_mf_scale_target(cx, S, C),
				C = C
			)
		} else {
			list(
				rule = scale,
				s = mean(diag(S$prior$Suv)),
				rho = mean(diag(S$prior$Suv[seq_len(R), R + seq_len(R), drop = FALSE])) / mean(diag(S$prior$Suv)),
				stat = if (!is.null(cx$scale_stat)) cx$scale_stat else .ae_mf_scale_stat(),
				target = .ae_mf_scale_target(cx, S, C),
				C = C
			)
		}
	}
	out
}
