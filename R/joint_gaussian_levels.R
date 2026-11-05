.lame_joint_levels_context <- function(X, Q) {
	n <- dim(X[[1L]])[1L]; p <- nrow(Q); nt <- length(X)
	H <- qr.Q(qr(stats::contr.helmert(n)))
	ij <- which(row(matrix(0, n, n)) != col(matrix(0, n, n)), arr.ind = TRUE)
	idx <- ij[, 1L] + (ij[, 2L] - 1L) * n
	mirror <- match(ij[, 2L] + (ij[, 1L] - 1L) * n, idx)
	D <- do.call(
		rbind,
		lapply(
			X,
			function(x)
				cbind(
					matrix(x, n * n, p)[idx, , drop = FALSE],
					H[ij[, 1L], , drop = FALSE],
					H[ij[, 2L], , drop = FALSE]
				)
		)
	)
	partner <- unlist(lapply(seq_len(nt), function(t) mirror + (t - 1L) * length(idx)))
	if (any(!is.finite(D)) || any(!is.finite(Q))) return(NULL)
	# symmetrize the coefficient precision used for the root and the score
	Q <- Q / 2 + t(Q) / 2
	native <- joint_levels_context_cpp(
		D,
		H,
		Q,
		as.integer(rep(ij[, 1L] - 1L, nt)),
		as.integer(rep(ij[, 2L] - 1L, nt)),
		as.integer(partner - 1L)
	)
	list(n = n, p = p, nt = nt, H = H, index = idx, native = native)
}

.lame_joint_levels_parameters <- function(Z, a_dev, b_dev, Sab, rho, cx) {
	if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1 || !all(is.finite(Sab)))
		stop("Invalid joint coefficient/level covariance.", call. = FALSE)
	residual <- unlist(lapply(seq_len(cx$nt), function(t) {
		z <- Z[, , t] - outer(a_dev[, t], b_dev[, t], "+")
		z[cx$index]
	}))
	if (any(!is.finite(residual))) stop("Nonfinite joint coefficient/level residual.", call. = FALSE)
	joint_levels_parameters_cpp(cx$native, residual, Sab, rho)
}

.lame_joint_levels_draw <- function(Z, a_dev, b_dev, Sab, rho, cx) {
	moments <- .lame_joint_levels_parameters(Z, a_dev, b_dev, Sab, rho, cx)
	q <- drop(moments$mean + moments$noise_root %*% rnorm(length(moments$mean)))
	p <- cx$p; k <- cx$n - 1L
	list(
		beta = q[seq_len(p)],
		a = drop(cx$H %*% q[p + seq_len(k)]),
		b = drop(cx$H %*% q[p + k + seq_len(k)]),
		q = q
	)
}
