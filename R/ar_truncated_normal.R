# scalar AR proposal from a normal truncated to a finite interval: inverse-cdf
# in the centre, log-scale tails. nothing is clamped, since a clamped value would
# put a point mass in the truncated-normal conditional
.lame_ar_truncnorm_draw <- function(mean, sd, lower, upper) {
	if (!all(is.finite(c(mean, sd, lower, upper))) || sd <= 0 || lower >= upper)
		stop("Invalid truncated-normal AR proposal parameters.", call. = FALSE)
	flo <- stats::pnorm(lower, mean, sd)
	fhi <- stats::pnorm(upper, mean, sd)
	for (attempt in seq_len(1000L)) {
		u <- stats::runif(1L)
		p <- flo + u * (fhi - flo)
		if (fhi > flo && p >= 1e-12 && p <= 1 - 1e-12) {
			value <- stats::qnorm(p, mean, sd)
		} else if (upper <= mean) {
			la <- stats::pnorm(lower, mean, sd, log.p = TRUE)
			lb <- stats::pnorm(upper, mean, sd, log.p = TRUE)
			logp <- lb + log(u + (1 - u) * exp(la - lb))
			value <- stats::qnorm(logp, mean, sd, log.p = TRUE)
		} else {
			la <- stats::pnorm(lower, mean, sd, lower.tail = FALSE, log.p = TRUE)
			lb <- stats::pnorm(upper, mean, sd, lower.tail = FALSE, log.p = TRUE)
			logp <- la + log((1 - u) + u * exp(lb - la))
			value <- stats::qnorm(logp, mean, sd, lower.tail = FALSE, log.p = TRUE)
		}
		# a rounded endpoint is redrawn, never projected into the support
		if (is.finite(value) && value > lower && value < upper) return(value)
	}
	stop("Truncated-normal AR interval cannot be resolved numerically.", call. = FALSE)
}
