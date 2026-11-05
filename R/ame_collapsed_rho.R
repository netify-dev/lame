.ame_binary_collapsed_rho_mh <- function(rho, eta, cx, tries = 4L) {
	# independence proposals from the prior of rho, so prior and proposal cancel
	old <- .lame_binary_dynamic_collapsed_rho_loglik(rho, eta, cx)
	if (!is.finite(old) || !is.finite(rho) || abs(rho) >= 1)
		stop("Invalid current AME reciprocity state.", call. = FALSE)
	accepted <- 0L
	for (k in seq_len(tries)) {
		proposal <- .lame_rho_prior_draw()
		if (abs(proposal) >= 1) next
		new <- .lame_binary_dynamic_collapsed_rho_loglik(proposal, eta, cx)
		if (!is.finite(new)) stop("Nonfinite observed AME reciprocity likelihood.", call. = FALSE)
		if (log(runif(1L)) < new - old) {
			rho <- proposal; old <- new; accepted <- accepted + 1L
		}
	}
	list(rho = rho, accepts = accepted, attempts = tries)
}
