# joint scale move of the predictor of a binary bipartite fit
# rare events identify the predictor scale weakly, so the intercept, effects and G trade along a ridge
# the move multiplies the whole predictor by c (factors held), a Metropolis step on log c
# with Z integrated out; a Gaussian block of dimension d and quadratic form q adds d log c - (c^2 - 1) q / 2
# Z is redrawn afterwards, which makes this a partially collapsed sweep

# observed cells of the panel, with the sign of each (+1 event, -1 non-event)
.lame_bpl_context <- function(Y) {
	i1 <- which(Y == 1)
	i0 <- which(Y == 0)
	list(
		idx = c(i1, i0),
		sY = c(rep(1, length(i1)), rep(-1, length(i0))),
		n1 = length(i1),
		n0 = length(i0)
	)
}

# one visit: `tries` Metropolis steps on the factor; `room` is the largest factor the G bound allows
# `sd` is the step of log c and `att` the proposal count read by the burn-in Robbins-Monro adaptation
.lame_bpl_scale <- function(
	ez_o,
	sY,
	beta,
	iV0,
	q_lev,
	n_lev,
	q_dev = 0,
	n_dev = 0L,
	q_g = 0,
	n_g = 0L,
	room = Inf,
	tries = 8L,
	sd = 0.02,
	adapt = FALSE,
	att = 0L
) {
	ctot <- 1
	ll_cur <- sum(pnorm(sY * ez_o, log.p = TRUE))
	acc <- 0L
	for (k in seq_len(tries)) {
		cc <- ctot * exp(rnorm(1, 0, sd))
		att <- att + 1L
		ok <- FALSE
		if (is.finite(cc) && cc > 0 && cc < room) {
			ll_new <- sum(pnorm(sY * (cc * ez_o), log.p = TRUE))
			lr <- ll_new - ll_cur + .lame_scale_prior_logratio(
				cc / ctot,
				ctot * beta,
				iV0,
				ctot * ctot * (q_lev + q_g),
				n_lev + n_g,
				ctot * ctot * q_dev,
				n_dev
			)
			ok <- is.finite(lr) && log(runif(1)) < lr
			if (ok) {
				ctot <- cc
				ll_cur <- ll_new
				acc <- acc + 1L
			}
		}
		# step adaptation in burn-in only
		if (adapt) sd <- min(max(exp(log(sd) + (as.numeric(ok) - 0.35) / sqrt(att)), 1e-5), 1)
	}
	list(c = ctot, sd = sd, att = att, acc = acc)
}
