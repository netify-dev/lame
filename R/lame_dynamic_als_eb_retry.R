# retry of a dynamic empirical-Bayes fit with static factors that ran out of
# passes on a slow tail

# a count fit whose equations are certified but whose noise variance or factor
# scales still move along a slow geometric tail (.lde_tail()) when the passes run
# out: those scalars are moved to the limits of their tails and the iteration
# resumes, up to `rounds` times. the fixed point is the same, only reached sooner
.lde_fac_tail_retry <- function(ctx, fit, max_iter, tol, rounds = 3L) {
	k <- 0L
	budget <- max(100L, as.integer(max_iter) %/% 2L)
	while (!isTRUE(fit$converged) && identical(fit$termination, "iteration limit") &&
			isTRUE(ctx$vehyp) && identical(ctx$fmode, "fac") && k < rounds) {
		H <- fit$tails
		if (is.null(H) || nrow(H) < 5L) break
		cur <- .lde_tail_scalars(ctx, fit$st)
		if (length(cur) != ncol(H)) break
		lim <- vapply(seq_len(ncol(H)), function(j) {
			tl <- .lde_tail(H[, j])
			if (!is.finite(tl$dist)) return(NA_real_)
			cur[j] + if (isTRUE(tl$mono)) sign(tl$move) * tl$dist else 0
		}, 0)
		if (anyNA(lim)) break
		st1 <- .lde_tail_put(ctx, fit$st, lim)
		if (is.null(st1)) break
		k <- k + 1L
		nxt <- tryCatch(.lde_fit_core_dyn(ctx, st1, max_iter = budget, tol = tol), error = function(e) NULL)
		if (is.null(nxt) || !is.finite(nxt$ev$cert$merit)) break
		tr <- nxt$trace
		tr$iter <- tr$iter + fit$iterations
		nxt$trace <- rbind(fit$trace, tr)
		nxt$iterations <- fit$iterations + nxt$iterations
		fit <- nxt
	}
	fit$tail_retries <- k
	fit
}

# a state with the watched scalars of .lde_tail_scalars() set to `v` (log ve, then
# the log scales of the factor prior); NULL when the prior would be invalid
.lde_tail_put <- function(ctx, st, v) {
	st$hyp$ve <- exp(v[1L])
	if (identical(ctx$fmode, "fac") && !is.null(st$hyp$fac) && length(v) > 1L) {
		pr <- st$hyp$fac
		sc <- exp(v[-1L])
		if (identical(pr$kind, "joint")) {
			d <- sqrt(sc / diag(pr$Suv))
			pr$Suv <- pr$Suv * outer(d, d)
		} else if (pr$kind %in% c("scalar", "columns")) {
			pr$tau <- sc
		} else {
			nu <- length(pr$su)
			pr$su <- sc[seq_len(nu)]
			pr$sv <- sc[-seq_len(nu)]
		}
		pr <- tryCatch(.lde_fac_fix(pr), error = function(e) NULL)
		if (is.null(pr)) return(NULL)
		st$hyp$fac <- pr
	}
	st
}
