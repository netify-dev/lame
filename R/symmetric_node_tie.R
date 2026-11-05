# ties the row and column coefficients of a nodal covariate in a symmetric binary or ordinal fit.
# the symmetric likelihood depends on the pair only through their mean, so the difference is
# unidentified and would make the carried predictor asymmetric by (b_row - b_col) (x_i - x_j);
# after every coefficient draw both are set to their mean, a Gibbs sampler for the model with one
# coefficient per nodal covariate (prior$symmetric_node_tie = FALSE turns it off)

# `beta` in the design layout; `rb` and `cb` the positions of the matching row and column slices
# of the nodal covariates
.lame_sym_tie <- function(beta, rb, cb) {
	if (length(rb) == 0L || length(rb) != length(cb) || length(beta) < max(rb, cb)) return(beta)
	m <- (beta[rb] + beta[cb]) / 2
	beta[rb] <- m
	beta[cb] <- m
	beta
}

# the same for a dynamic_beta state: `path` the N x length(dyn_idx) dynamic paths, `stat` the static
# coefficients, `dyn_idx` and `stat_idx` their positions; a nodal covariate's two slices are both
# dynamic or both static
.lame_sym_tie_state <- function(path, stat, dyn_idx, stat_idx, rb, cb) {
	for (k in seq_along(rb)) {
		jr <- match(rb[k], dyn_idx)
		jc <- match(cb[k], dyn_idx)
		if (!is.na(jr) && !is.na(jc)) {
			m <- (path[, jr] + path[, jc]) / 2
			path[, jr] <- m
			path[, jc] <- m
			next
		}
		sr <- match(rb[k], stat_idx)
		sc <- match(cb[k], stat_idx)
		if (!is.na(sr) && !is.na(sc)) {
			m <- (stat[sr] + stat[sc]) / 2
			stat[sr] <- m
			stat[sc] <- m
		}
	}
	list(path = path, stat = stat)
}

# number of parameter dimensions the tie removes from the scale moves (R/binary_scale_orbit.R,
# R/ordinal_scale_move.R), which count coefficients in their Jacobian c^(d - 1): a tied pair is one
# parameter (a static pair one dimension, a pair of paths one per period)
.lame_sym_tie_dims <- function(rb, cb, dyn_idx, stat_idx, nt, any_dyn) {
	if (!any_dyn) return(length(rb))
	sum(vapply(seq_along(rb), function(k) {
		if (!is.na(match(rb[k], dyn_idx)) && !is.na(match(cb[k], dyn_idx))) nt
		else if (!is.na(match(rb[k], stat_idx)) && !is.na(match(cb[k], stat_idx))) 1L
		else 0L
	}, 0))
}
