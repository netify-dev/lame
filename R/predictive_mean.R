# replicates averaged at each stored draw where the replicate's mean has no closed form
.predictive_replicates <- 10L

# term added to the running sum YPS at a stored draw: the mean of the replicate
# given that draw where the family has one in closed form (binary, normal,
# explicit-cutpoint ordinal), otherwise the average of the replicate and m - 1
# further ones from `draw` (cbin, frn, ordinal with data-induced cutpoints), or
# the replicate itself when no `draw` is given. s2 is the residual variance the
# binary replicate was drawn with. the further replicates come from a stream
# seeded by `seed` and the stored-draw `index`, and the sampler's stream is
# restored afterwards. cells where the replicate is missing or the predictor is
# not finite keep the replicate's value
.predictive_mean_term <- function(
	Ys,
	EZ,
	family,
	s2 = 1,
	alpha = NULL,
	levels = NULL,
	symmetric = FALSE,
	draw = NULL,
	seed = NULL,
	index = NULL
) {
	explicit_ordinal <- identical(family, "ordinal") && !is.null(alpha)
	closed_form <- family %in% c("binary", "normal") || explicit_ordinal
	replicated <- !closed_form && family %in% c("cbin", "frn", "ordinal") && is.function(draw)
	if (!closed_form && !replicated) return(Ys)
	if (isTRUE(symmetric)) {
		if (length(dim(EZ)) == 3L) {
			for (k in seq_len(dim(EZ)[3L])) EZ[, , k] <- (EZ[, , k] + t(EZ[, , k])) / 2
		} else {
			EZ <- (EZ + t(EZ)) / 2
		}
	}
	keep <- is.finite(EZ) & !is.na(Ys)
	if (replicated) {
		mu <- .predictive_replicate_mean(Ys, draw, seed, index)
	} else {
		eta <- EZ
		eta[!keep] <- NA
		mu <- switch(
			family,
			binary = stats::pnorm(eta / sqrt(s2)),
			normal = eta,
			ordinal = .ordinal_response_mean(eta, alpha, levels)
		)
	}
	out <- Ys
	out[keep] <- mu[keep]
	out
}

# average of Ys and m - 1 further replicates drawn on their own stream; a cell
# an extra replicate leaves missing takes the value of Ys
.predictive_replicate_mean <- function(Ys, draw, seed, index, m = .predictive_replicates) {
	had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
	on.exit({
		if (had_seed) {
			assign(".Random.seed", old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	})
	set.seed(.predictive_stream_seed(seed, index))
	total <- Ys
	for (r in seq_len(m - 1L)) {
		y_r <- draw()
		gap <- is.na(y_r)
		if (any(gap)) y_r[gap] <- Ys[gap]
		total <- total + y_r
	}
	total / m
}

# seed of the extra stream at a stored draw: a mix of the chain's seed and the
# stored-draw index, so it does not depend on the sampler's state. the seeds of
# up to 64 chains (1000 apart) give distinct streams for indices under two million
.predictive_stream_seed <- function(seed, index) {
	p <- 2147483647
	s <- if (is.numeric(seed) && length(seed) == 1L && is.finite(seed)) seed else 6886
	i <- if (is.numeric(index) && length(index) == 1L && is.finite(index)) index else 0
	as.integer(((round(s) %% p) * 48271 %% p + (round(i) %% p) * 16807) %% p)
}
