# observed-data HMC with static centered actor effects and static factors;
# coefficient paths keep their AR(1)/RW1 prior and the inverse-Wishart factor
# covariance is integrated out
.lame_binary_static_hmc_context <- function(Y, Xlist, static_idx, dynamic_idx, rank, uv_prior = NULL) {
	n <- dim(Y)[1]
	nt <- dim(Y)[3]
	H <- qr.Q(qr(stats::contr.helmert(n)))
	k <- n - 1L
	p <- dim(Xlist[[1]])[3]
	obs <- is.finite(Y)
	for (t in 1:nt) diag(obs[, , t]) <- FALSE
	y <- Y
	y[!obs] <- 0
	sgn <- 2 * y - 1
	X <- lapply(Xlist, function(x) matrix(x, n * n, p))
	rate <- mean(Y[obs])
	rate <- max(0.001, min(0.999, rate))
	weight <- stats::dnorm(stats::qnorm(rate))^2 / (rate * (1 - rate))
	dims <- list(
		static = length(static_idx),
		path = c(nt, length(dynamic_idx)),
		level = c(
			k,
			2L
		),
		U = c(n, rank),
		V = c(n, rank)
	)
	sizes <- vapply(dims, prod, numeric(1))
	ends <- cumsum(sizes)
	starts <- ends - sizes + 1
	indices <- lapply(seq_along(sizes), function(i) seq.int(starts[i], ends[i]))
	names(indices) <- names(dims)
	# inverse-Wishart factor prior (R/uv_prior.R): the supplied `uv_prior`, else the full-information default
	uv_df <- if (is.null(uv_prior)) .lame_uv_prior_df(rank) else uv_prior$df
	uv_psi <- if (is.null(uv_prior)) .lame_uv_prior_df(rank) else uv_prior$psi
	list(
		n = n,
		nt = nt,
		k = k,
		p = p,
		rank = rank,
		H = H,
		obs = obs,
		sgn = sgn,
		X = X,
		static_idx = static_idx,
		dynamic_idx = dynamic_idx,
		dims = dims,
		indices = indices,
		weight = weight,
		uv_df = uv_df,
		uv_psi = uv_psi
	)
}

.lame_binary_static_hmc_unpack <- function(q, cx) {
	out <- lapply(seq_along(cx$dims), function(j) {
		d <- cx$dims[[j]]
		z <- q[cx$indices[[j]]]
		if (length(d) > 1L)
			dim(z) <- d
		z
	})
	names(out) <- names(cx$dims)
	out
}

.lame_binary_static_hmc_pack <- function(z) unlist(z, use.names = FALSE)

.lame_binary_static_hmc_ar <- function(path, rho, var, initial_precision) {
	nt <- nrow(path)
	p <- ncol(path)
	gr <- matrix(0, nt, p)
	gr[1, ] <- -drop(initial_precision %*% path[1, ])
	q <- -sum(gr[1, ] * path[1, ])
	if (nt > 1L) {
		e <- path[-1, , drop = FALSE] - sweep(path[-nt, , drop = FALSE], 2, rho, "*")
		v <- sweep(e, 2, var, "/")
		q <- q + sum(e * v)
		gr[-1, ] <- gr[-1, , drop = FALSE] - v
		gr[-nt, ] <- gr[-nt, , drop = FALSE] + sweep(v, 2, rho, "*")
	}
	list(logp = -q / 2, gradient = gr)
}

.lame_binary_static_hmc_evaluate <- function(q, cx, h) {
	z <- .lame_binary_static_hmc_unpack(q, cx)
	n <- cx$n
	nt <- cx$nt
	H <- cx$H
	lev <- H %*% z$level
	uv <- z$U %*% t(z$V)
	eta <- array(0, c(n, n, nt))
	gb <- matrix(0, nt, cx$p)
	for (t in seq_len(nt)) {
		b <- numeric(cx$p)
		b[cx$static_idx] <- z$static
		b[cx$dynamic_idx] <- z$path[t, ]
		eta[, , t] <- matrix(cx$X[[t]] %*% b, n, n) + outer(lev[, 1], lev[, 2], "+") + uv
	}
	v <- cx$sgn[cx$obs] * eta[cx$obs]
	lp <- stats::pnorm(v, log.p = TRUE)
	lg <- sum(lp)
	ge <- array(0, dim(eta))
	ge[cx$obs] <- cx$sgn[cx$obs] * exp(stats::dnorm(v, log = TRUE) - lp)
	ga <- apply(ge, c(1, 3), sum)
	gbactor <- apply(ge, c(2, 3), sum)
	for (t in seq_len(nt)) gb[t, ] <- drop(crossprod(cx$X[[t]], c(ge[, , t])))
	gs <- colSums(gb[, cx$static_idx, drop = FALSE]) - z$static / h$g
	lg <- lg - 0.5 * sum(z$static^2 / h$g)
	ar <- .lame_binary_static_hmc_ar(z$path, h$rho_beta, h$var_beta, h$iC0)
	gp <- gb[, cx$dynamic_idx, drop = FALSE] + ar$gradient
	lg <- lg + ar$logp
	gl <- crossprod(H, cbind(rowSums(ga), rowSums(gbactor))) - z$level %*% h$iSab
	lg <- lg - 0.5 * sum(z$level * (z$level %*% h$iSab))
	# the SVD gives W (I + W'W)^(-1) stably when a factor draw has a large singular value
	W <- cbind(z$U, z$V)
	sv <- svd(W)
	fac <- cx$n + cx$uv_df
	gW <- -fac * sweep(sv$u, 2, sv$d / (cx$uv_psi + sv$d^2), "*") %*% t(sv$v)
	lg <- lg - fac / 2 * (2 * cx$rank * log(cx$uv_psi) + sum(log1p(sv$d^2 / cx$uv_psi)))
	e <- apply(ge, c(1, 2), sum)
	gu <- e %*% z$V + gW[, seq_len(cx$rank), drop = FALSE]
	gv <- t(e) %*% z$U + gW[, cx$rank + seq_len(cx$rank), drop = FALSE]
	list(logp = lg, gradient = .lame_binary_static_hmc_pack(list(gs, gp, gl, gu, gv)), eta = eta)
}

# diagonal mass: prior precision plus a data-based Fisher scale, held fixed within a trajectory
.lame_binary_static_hmc_mass <- function(cx, h) {
	ardiag <- function(nt, rho, var, initial) {
		q <- matrix(0, nt, length(rho))
		q[1, ] <- initial
		if (nt > 1) {
			q[-1, ] <- sweep(q[-1, , drop = FALSE], 2, 1 / var, "+")
			q[-nt, ] <- sweep(q[-nt, , drop = FALSE], 2, rho^2 / var, "+")
		}
		q
	}
	xx <- t(vapply(cx$X, function(x) colSums(x^2), numeric(cx$p))) * cx$weight
	ms <- 1 / h$g + colSums(xx[, cx$static_idx, drop = FALSE])
	mp <- ardiag(cx$nt, h$rho_beta, h$var_beta, diag(h$iC0)) + xx[, cx$dynamic_idx, drop = FALSE]
	ml <- matrix(rep(diag(h$iSab) + cx$nt * (cx$n - 1) * cx$weight, each = cx$k), cx$k, 2)
	mu <- matrix((cx$n + cx$uv_df) / cx$uv_psi + cx$nt * cx$n * cx$weight, cx$n, cx$rank)
	pmax(.lame_binary_static_hmc_pack(list(ms, mp, ml, mu, mu)), 1e-08)
}

# leapfrog integrator; a trajectory that turns non-finite is rejected whole
.lame_binary_static_hmc_integrate <- function(q, p, cx, h, mass, step, L) {
	z <- .lame_binary_static_hmc_evaluate(q, cx, h)
	# bounded force keeps the proposal volume-preserving and reversible; the
	# Metropolis ratio still uses the exact target
	force <- function(g) 20 * sqrt(mass) * tanh(g / (20 * sqrt(mass)))
	p <- p + step / 2 * force(z$gradient)
	for (k in seq_len(L)) {
		q <- q + step * p / mass
		z <- tryCatch(.lame_binary_static_hmc_evaluate(q, cx, h), error = function(e) NULL)
		if (is.null(z) || !is.finite(z$logp) || any(!is.finite(z$gradient)))
			return(NULL)
		p <- p + step * (if (k == L)
				0.5
			else 1) * force(z$gradient)
	}
	list(q = q, p = p, evaluation = z)
}

.lame_binary_static_hmc_step <- function(q, cx, h, step, L = sample(12:20, 1)) {
	mass <- .lame_binary_static_hmc_mass(cx, h)
	p <- stats::rnorm(length(q)) * sqrt(mass)
	old <- .lame_binary_static_hmc_evaluate(q, cx, h)
	new <- .lame_binary_static_hmc_integrate(q, p, cx, h, mass, step, L)
	if (is.null(new))
		return(list(q = q, accepted = FALSE, prob = 0, nonfinite = TRUE))
	lr <- .lame_hmc_log_acceptance(new$evaluation$logp, old$logp, new$p, p, mass)
	if (!is.finite(lr))
		return(list(q = q, accepted = FALSE, prob = 0, nonfinite = TRUE))
	prob <- min(1, exp(lr))
	accept <- is.finite(prob) && stats::runif(1) < prob
	list(
		q = if (accept) new$q else q,
		accepted = accept,
		prob = if (is.finite(prob)) prob else 0,
		nonfinite = FALSE
	)
}
