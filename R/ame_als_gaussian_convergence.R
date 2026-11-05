# centered coordinates remove the additive and factor-location ambiguities
# without changing any fitted value; a balanced factorization removes scale gauges
.ae_gaussian_context <- function(Z, X, R, symmetric = FALSE) {
	nr <- dim(Z)[1]; nc <- dim(Z)[2]; nt <- dim(Z)[3]; p <- dim(X)[3]
	contrast <- function(n) {
		if (n == 1L) return(matrix(0, 1L, 0L))
		q <- stats::contr.helmert(n)
		sweep(q, 2L, sqrt(colSums(q^2)), "/")
	}
	ok <- is.finite(Z)
	Z[!ok] <- 0
	X[!is.finite(X)] <- 0
	sizes <- if (symmetric) {
		c(1L, p, nr - 1L, 0L, (nr - 1L) * R, 0L)
	} else {
		c(1L, p, nr - 1L, nc - 1L, (nr - 1L) * R, (nc - 1L) * R)
	}
	ends <- cumsum(sizes)
	ix <- lapply(seq_along(sizes), function(k) {
		if (sizes[k]) seq.int(ends[k] - sizes[k] + 1L, ends[k]) else integer()
	})
	# response unit (sd of the observed working response), so the score does not
	# depend on the units of y
	zo <- Z[ok]
	unit <- if (length(zo) > 1L) stats::sd(zo) else NA_real_
	if (!is.finite(unit) || unit <= 0) unit <- if (length(zo)) sqrt(mean(zo^2)) else 1
	if (!is.finite(unit) || unit <= 0) unit <- 1
	list(
		nr = nr,
		nc = nc,
		nt = nt,
		p = p,
		R = R,
		symmetric = symmetric,
		sign = rep(1, R),
		qr = contrast(nr),
		qc = contrast(nc),
		z = Z,
		x = X,
		ok = ok,
		ws = apply(1 * ok, c(1L, 2L), sum),
		ix = ix,
		unit = unit
	)
}

.ae_gaussian_encode <- function(fit, context) {
	cx <- context
	if (cx$symmetric) {
		factor <- if (cx$R) sweep(fit$U, 2L, sqrt(abs(fit$L)), "*") else
			matrix(0, cx$nr, 0L)
		level <- colMeans(factor)
		center <- sweep(factor, 2L, level)
		a <- fit$a - mean(fit$a) + as.vector(center %*% (cx$sign * level))
		return(c(
			fit$mu + 2 * mean(fit$a) + sum(cx$sign * level^2),
			fit$beta,
			c(crossprod(cx$qr, a)),
			c(crossprod(cx$qr, center))
		))
	}
	U <- if (cx$R) fit$U else matrix(0, cx$nr, 0L)
	V <- if (cx$R) fit[["V"]] else matrix(0, cx$nc, 0L)
	mean_u <- colMeans(U); mean_v <- colMeans(V)
	U <- sweep(U, 2L, mean_u); V <- sweep(V, 2L, mean_v)
	a <- fit$a - mean(fit$a) + as.vector(U %*% mean_v)
	b <- fit$b - mean(fit$b) + as.vector(V %*% mean_u)
	c(
		fit$mu + mean(fit$a) + mean(fit$b) + sum(mean_u * mean_v),
		fit$beta,
		c(crossprod(cx$qr, a)),
		c(crossprod(cx$qc, b)),
		c(crossprod(cx$qr, U)),
		c(crossprod(cx$qc, V))
	)
}

.ae_gaussian_decode <- function(theta, cx) {
	a <- as.vector(cx$qr %*% theta[cx$ix[[3L]]])
	U <- cx$qr %*% matrix(theta[cx$ix[[5L]]], cx$nr - 1L, cx$R)
	b <- if (cx$symmetric) a else as.vector(cx$qc %*% theta[cx$ix[[4L]]])
	V <- if (cx$symmetric) U else
		cx$qc %*% matrix(theta[cx$ix[[6L]]], cx$nc - 1L, cx$R)
	list(
		mu = theta[cx$ix[[1L]]],
		beta = theta[cx$ix[[2L]]],
		a = a,
		b = b,
		U = U,
		V = V
	)
}

# gradient and diagonal Gauss-Newton curvature over the observed cells of the
# sse, counting both copies of a symmetric response
.ae_gaussian_eval <- function(theta, cx) {
	fit <- .ae_gaussian_decode(theta, cx)
	O <- if (cx$symmetric) {
		tcrossprod(sweep(fit$U, 2L, cx$sign, "*"), fit$U)
	} else tcrossprod(fit$U, fit[["V"]])
	eta <- array(rep(fit$mu + outer(fit$a, fit$b, "+") + O, cx$nt), dim(cx$z))
	if (cx$p) for (k in seq_len(cx$p)) {
		eta <- eta + fit$beta[k] * array(cx$x[, , k, ], dim(eta))
	}
	error <- eta - cx$z
	error[!cx$ok] <- 0
	es <- apply(error, c(1L, 2L), sum); ws <- cx$ws
	gb <- hb <- numeric(cx$p)
	if (cx$p) for (k in seq_len(cx$p)) {
		x <- array(cx$x[, , k, ], dim(error))
		gb[k] <- 2 * sum(error * x)
		hb[k] <- 2 * sum(cx$ok * x^2)
	}
	if (cx$symmetric) {
		ga <- 2 * crossprod(cx$qr, rowSums(es) + colSums(es))
		gf <- 2 * crossprod(cx$qr, (es + t(es)) %*% sweep(fit$U, 2L, cx$sign, "*"))
		ha <- 2 * colSums(cx$qr * ((diag(rowSums(ws) + colSums(ws)) + ws + t(ws)) %*% cx$qr))
		hf <- matrix(0, cx$nr - 1L, cx$R)
		if (cx$R) for (k in seq_len(cx$R)) {
			vk <- fit$U[, k]
			K <- diag(as.vector((ws + t(ws)) %*% (vk^2))) + 2 * ws * outer(vk, vk)
			hf[, k] <- 2 * colSums(cx$qr * (K %*% cx$qr))
		}
		gradient <- c(2 * sum(error), gb, c(ga), c(gf))
		information <- c(2 * sum(cx$ok), hb, c(ha), c(hf))
	} else {
		gradient <- c(
			2 * sum(error),
			gb,
			c(2 * crossprod(cx$qr, rowSums(es))),
			c(2 * crossprod(cx$qc, colSums(es))),
			c(2 * crossprod(cx$qr, es %*% fit[["V"]])),
			c(2 * crossprod(cx$qc, t(es) %*% fit$U))
		)
		information <- c(
			2 * sum(cx$ok),
			hb,
			c(2 * crossprod(cx$qr^2, rowSums(ws))),
			c(2 * crossprod(cx$qc^2, colSums(ws))),
			c(2 * crossprod(cx$qr^2, ws %*% (fit[["V"]]^2))),
			c(2 * crossprod(cx$qc^2, t(ws) %*% (fit$U^2)))
		)
	}
	fit$Omat <- O
	fit$L <- if (cx$symmetric) cx$sign else NULL
	value <- sum(error^2)
	# score: largest gradient over the root of its curvature, divided by the
	# response sd so the threshold is unit free (the residual sd would vanish for a
	# saturated fit)
	sigma <- cx$unit
	list(
		value = value,
		gradient = gradient,
		information = information,
		score = max(abs(gradient) / sqrt(pmax(information, 1e-30))) / sigma,
		sigma = sigma,
		fit = fit
	)
}

.ae_gaussian_canonical <- function(fit, cx) {
	if (!cx$R) return(fit)
	O <- fit$Omat
	j <- seq_len(cx$R)
	if (cx$symmetric) {
		ee <- eigen((O + t(O)) / 2, symmetric = TRUE)
		j <- order(abs(ee$values), decreasing = TRUE)[j]
		fit$U <- ee$vectors[, j, drop = FALSE]
		fit[["V"]] <- fit$U
		fit$L <- ee$values[j]
	} else {
		ss <- svd(O)
		root <- sqrt(pmax(ss$d[j], 0))
		fit$U <- sweep(ss$u[, j, drop = FALSE], 2L, root, "*")
		fit[["V"]] <- sweep(ss$v[, j, drop = FALSE], 2L, root, "*")
		fit$L <- NULL
	}
	fit
}

# number of L-BFGS refinement rounds (max_iter iterations each) for a Gaussian
# fit that plateaus above its score threshold; flat objectives need several
.ae_gaussian_rounds <- 5L

.ae_gaussian_certify <- function(fit, Z, X, R, symmetric, max_iter, tol) {
	cx <- .ae_gaussian_context(Z, X, R, symmetric)
	canonical <- .ae_gaussian_canonical(fit, cx)
	if (symmetric) cx$sign <- ifelse(canonical$L < 0, -1, 1)
	check <- .ae_gaussian_eval(.ae_gaussian_encode(canonical, cx), cx)
	initial_score <- check$score
	initial_sse <- fit$sse
	threshold <- sqrt(tol)
	rounds <- list()
	# refine only fits that plateau above the score threshold; certified fits are
	# kept exactly and a single requested sweep stays an unfinished refit
	if (fit$iterations > 1L && is.finite(check$score) && check$score > threshold) {
		for (round in seq_len(.ae_gaussian_rounds)) {
			start <- check
			theta <- .ae_gaussian_encode(canonical, cx)
			scale <- sqrt(pmax(start$information, 1e-6))
			cache <- new.env(parent = emptyenv())
			cache$q <- NULL
			evaluate <- function(q) {
				if (!identical(q, cache$q)) {
					cache$value <- .ae_gaussian_eval(q / scale, cx)
					cache$q <- q
				}
				cache$value
			}
			opt <- tryCatch(
				stats::optim(
					theta * scale,
					function(q) evaluate(q)$value,
					function(q) evaluate(q)$gradient / scale,
					method = "L-BFGS-B",
					control = list(
						maxit = max_iter,
						factr = 1,
						pgtol = threshold * 0.1 * start$sigma
					)
				),
				error = identity
			)
			if (inherits(opt, "error")) {
				rounds[[round]] <- list(accepted = FALSE, message = conditionMessage(opt))
				break
			}
			end <- .ae_gaussian_eval(opt$par / scale, cx)
			accepted <- is.finite(end$value) && end$value <= start$value + 1e-10 * max(start$value, cx$unit^2)
			rounds[[round]] <- list(
				initial = start$value,
				final = end$value,
				code = opt$convergence,
				counts = opt$counts,
				accepted = accepted
			)
			if (!accepted) break
			# indefinite symmetric factor rotations can grow without changing their
			# product, so refactor the product before checking convergence
			canonical <- .ae_gaussian_canonical(end$fit, cx)
			if (symmetric) cx$sign <- ifelse(canonical$L < 0, -1, 1)
			check <- .ae_gaussian_eval(.ae_gaussian_encode(canonical, cx), cx)
			if (!is.finite(check$value) ||
					abs(check$value - end$value) > 1e-9 * max(cx$unit^2, end$value)) {
				cli::cli_abort("Gaussian factor refactoring changed the fitted objective.")
			}
			rounds[[round]]$score <- check$score
			for (name in c("mu", "beta", "a", "b", "U", "V", "L", "Omat")) {
				fit[name] <- canonical[name]
			}
			# rank-zero fits report absent factors as NULL
			if (!R) fit[c("U", "V", "L")] <- list(NULL, NULL, NULL)
			fit$sse <- check$value
			fit$dev_history[length(fit$dev_history)] <- fit$sse
			if (is.finite(check$score) && check$score <= threshold) break
		}
	}
	pass <- is.finite(check$score) && check$score <= threshold
	refined <- any(vapply(rounds, function(z) isTRUE(z$accepted), logical(1)))
	fit$converged <- (isTRUE(fit$converged) || refined) && pass
	fit$gaussian_diagnostics <- list(
		pass = pass,
		max_score = check$score,
		threshold = threshold,
		initial_score = initial_score,
		initial_sse = initial_sse,
		final_sse = fit$sse,
		refinement = rounds,
		criterion = "maximum absolute centered SSE gradient divided by square root of diagonal Gauss-Newton curvature"
	)
	fit
}
