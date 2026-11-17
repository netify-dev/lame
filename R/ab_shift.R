# exact draw of the shift a_i + c, b_j - c of the additive effects of a fit without an intercept.
# the likelihood sees only a_i + b_j, so c is normal under the prior of the effects:
# independent margins give precision P = nA / va + nB / vb and mean (sum(b) / vb - sum(a) / va) / P;
# paired effects (a_i, b_i) ~ N(0, Sab) give P = n d'Qd and mean -d'Q (sum(a), sum(b)) / P,
# with Q = Sab^-1 and d = (1, -1)

# `paired` uses the full Sab (unipartite), otherwise its diagonal (bipartite); 0 when undefined
.lame_ab_shift_draw <- function(a, b, Sab, paired = TRUE) {
	na <- length(a); nb <- length(b)
	if (na < 1L || nb < 1L || !all(is.finite(a)) || !all(is.finite(b))) return(0)
	if (!is.matrix(Sab) || !all(dim(Sab) == 2L) || !all(is.finite(Sab))) return(0)
	if (paired && na == nb) {
		# no shift, given Sab that this draw does not move, leaves the posterior invariant
		Q <- tryCatch(solve(Sab), error = function(e) NULL)
		if (is.null(Q)) {
			.lame_fallback_note('ab_shift')
			return(0)
		}
		qd <- drop(Q %*% c(1, -1))
		P <- na * (qd[1L] - qd[2L])
		m <- -(qd[1L] * sum(a) + qd[2L] * sum(b)) / P
	} else {
		va <- Sab[1L, 1L]; vb <- Sab[2L, 2L]
		if (va <= 0 || vb <= 0) return(0)
		P <- na / va + nb / vb
		m <- (sum(b) / vb - sum(a) / va) / P
	}
	if (!is.finite(P) || P <= 0 || !is.finite(m)) {
		.lame_fallback_note('ab_shift')
		return(0)
	}
	rnorm(1L, m, 1 / sqrt(P))
}
