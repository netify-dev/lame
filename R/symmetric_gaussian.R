.lame_symmetric_gaussian_context <- function(X, Q, actor = FALSE, center = TRUE) {
	n <- dim(X[[1L]])[1L]; p <- nrow(Q)
	ij <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
	idx <- ij[, 1L] + (ij[, 2L] - 1L) * n
	mirror <- ij[, 2L] + (ij[, 1L] - 1L) * n
	H <- if (!isTRUE(actor)) matrix(0, n, 0L) else if (isTRUE(center)) {
		if (n > 1L) qr.Q(qr(stats::contr.helmert(n))) else matrix(0, n, 0L)
	} else diag(n)
	D <- do.call(rbind, lapply(X, function(x) {
		x <- matrix(x, n * n, p)
		cbind(
			(x[idx, , drop = FALSE] + x[mirror, , drop = FALSE]) / 2,
			H[ij[, 1L], , drop = FALSE] + H[ij[, 2L], , drop = FALSE]
		)
	}))
	if (any(!is.finite(D)) || any(!is.finite(Q)))
		stop("Nonfinite symmetric Gaussian design or prior.", call. = FALSE)
	list(
		n = n,
		p = p,
		k = ncol(H),
		nt = length(X),
		H = H,
		index = idx,
		mirror = mirror,
		D = D,
		precision = crossprod(D),
		Q = Q
	)
}

.lame_symmetric_gaussian_parameters <- function(Z, variance, cx) {
	if (cx$k && (!is.finite(variance) || variance <= 0))
		stop("Invalid symmetric actor variance.", call. = FALSE)
	y <- unlist(lapply(seq_len(cx$nt), function(t) {
		z <- Z[, , t]; (z[cx$index] + z[cx$mirror]) / 2
	}))
	if (any(!is.finite(y))) stop("Nonfinite symmetric latent response.", call. = FALSE)
	q <- cx$precision
	if (cx$p) q[seq_len(cx$p), seq_len(cx$p)] <- q[seq_len(cx$p), seq_len(cx$p)] + cx$Q
	if (cx$k) {
		j <- cx$p + seq_len(cx$k)
		q[j, j] <- q[j, j] + diag(1 / variance, cx$k)
	}
	if (!nrow(q)) return(list(mean = numeric(), root = q))
	root <- chol((q + t(q)) / 2)
	rhs <- drop(crossprod(cx$D, y))
	list(mean = drop(backsolve(root, forwardsolve(t(root), rhs))), root = root)
}

.lame_symmetric_gaussian_draw <- function(Z, variance, cx) {
	z <- .lame_symmetric_gaussian_parameters(Z, variance, cx)
	q <- if (length(z$mean)) z$mean + drop(backsolve(z$root, rnorm(length(z$mean)))) else numeric()
	a <- if (cx$k) drop(cx$H %*% q[cx$p + seq_len(cx$k)]) else numeric(cx$n)
	list(beta = q[seq_len(cx$p)], a = a, b = a, q = q)
}
