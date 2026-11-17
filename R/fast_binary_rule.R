# the fast route's density warning for binary networks: its 95% interval coverage (intercept / dyadic
# slope) on rank-2 directed binary networks of 30 to 100 actors was 0.22-0.33 / 0.74-0.89 at density
# 0.10, 0.73-0.85 / 0.93-0.95 at 0.20 and 0.91-0.94 / 0.93-0.97 at 0.30, so a user-facing fast call
# (ame_als(), lame_als(), ame() and lame() with method = "als") warns once when min(d, 1 - d) is below
# the threshold. the fit is the fast route's either way; internal fits (the sampler's warm start,
# bootstrap replicates, refits, restarts) do not warn

# min(d, 1 - d) below which a fast binary fit warns
.lame_fast_binary_min_density <- 0.20

# observed density of a binary network: the share of observed cells that are ties, the diagonal of a
# one-mode network left out, pooled over the periods of a panel (a list of matrices or a 3-d array)
.lame_fast_binary_density <- function(Y, mode = "unipartite") {
	slices <- if (is.list(Y) && !is.data.frame(Y)) {
		Y
	} else if (is.array(Y) && length(dim(Y)) == 3L) {
		lapply(seq_len(dim(Y)[3L]), function(t) Y[, , t])
	} else {
		list(Y)
	}
	ties <- 0
	cells <- 0
	for (y in slices) {
		y <- as.matrix(y)
		if (!is.numeric(y) && !is.logical(y)) return(NA_real_)
		ok <- is.finite(y)
		if (!identical(mode, "bipartite") && nrow(y) == ncol(y)) diag(ok) <- FALSE
		ties <- ties + sum(y[ok] > 0)
		cells <- cells + sum(ok)
	}
	if (cells == 0) NA_real_ else ties / cells
}

# warn when a user-facing fast fit of a binary network has a density outside the band; silent when the
# density cannot be computed. `sampler` names the sampler for the same model
.lame_fast_binary_warn <- function(Y, family, mode, sampler = "ame") {
	if (!identical(family, "binary")) return(invisible(NULL))
	d <- tryCatch(.lame_fast_binary_density(Y, as.character(mode)[1L]), error = function(e) NA_real_)
	lo <- .lame_fast_binary_min_density
	if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || min(d, 1 - d) >= lo) return(invisible(NULL))
	d_txt <- format(round(d, 3), nsmall = 3)
	cli::cli_warn(c(
		"!" = "Binary network with density {d_txt}: outside {lo} to {1 - lo} the fast route tends to pull slopes toward zero and bias the intercept, and its intervals cover less often than stated.",
		"i" = "The sampler ({.fn {sampler}}) is more accurate there."
	), class = "lame_fast_binary_density")
	invisible(d)
}
