# explicit ordinal predictions use the cutpoints of the same posterior draw.
.ordinal_validate_cutpoints <- function(alpha, levels) {
	if (!is.numeric(levels) || length(levels) < 2L || any(!is.finite(levels)) ||
			any(diff(levels) <= 0) || !is.numeric(alpha) ||
			length(alpha) != length(levels) - 1L || any(!is.finite(alpha)) ||
			alpha[1L] != 0 || any(diff(alpha) <= 0))
		cli::cli_abort("Explicit ordinal prediction requires ordered category levels and matching finite cutpoints anchored at zero.")
	invisible(NULL)
}

.ordinal_saved_cutpoints <- function(fit, s) {
	levels <- fit$ordinal_levels
	if (length(s) != 1L || !is.finite(s) || s != floor(s) || s < 1L ||
			s > dim(fit$BETA)[1L])
		cli::cli_abort("Ordinal prediction requires a valid matched posterior draw index.")
	alpha <- 0
	if (length(levels) > 2L) {
		if (!is.matrix(fit$ALPHA) || nrow(fit$ALPHA) != dim(fit$BETA)[1L] ||
				ncol(fit$ALPHA) != length(levels) - 2L)
			cli::cli_abort("Explicit ordinal prediction requires aligned saved ALPHA cutpoint draws.")
		alpha <- c(0, fit$ALPHA[s, ])
	}
	.ordinal_validate_cutpoints(alpha, levels)
	alpha
}

.ordinal_categories <- function(Z, alpha, levels) {
	.ordinal_validate_cutpoints(alpha, levels)
	if (any(!is.finite(Z) & !is.na(Z)))
		cli::cli_abort("Ordinal prediction requires finite latent draws or missing cells.")
	out <- Z
	keep <- !is.na(Z)
	out[keep] <- levels[findInterval(Z[keep], alpha, left.open = TRUE) + 1L]
	out
}

.ordinal_response_mean <- function(eta, alpha, levels) {
	.ordinal_validate_cutpoints(alpha, levels)
	if (any(!is.finite(eta) & !is.na(eta)))
		cli::cli_abort("Ordinal prediction requires finite predictors or missing cells.")
	out <- eta * 0 + levels[1L]
	for (k in seq_along(alpha))
		out <- out + (levels[k + 1L] - levels[k]) * stats::pnorm(eta - alpha[k])
	out
}

.simY_ord_explicit <- function(EZ, rho, alpha, levels, symmetric = FALSE) {
	.ordinal_validate_cutpoints(alpha, levels)
	if (!is.matrix(EZ) || nrow(EZ) != ncol(EZ) || any(!is.finite(EZ) & !is.na(EZ)) ||
			length(rho) != 1L || !is.finite(rho) || abs(rho) > 1)
		cli::cli_abort("Explicit ordinal simulation requires a square predictor with finite or missing cells and valid reciprocity.")
	if (isTRUE(symmetric)) EZ <- EZ / 2 + t(EZ) / 2
	Z <- simZ(EZ, if (isTRUE(symmetric)) 0 else rho)
	diag(Z) <- NA_real_
	if (isTRUE(symmetric)) Z[lower.tri(Z)] <- t(Z)[lower.tri(Z)]
	.ordinal_categories(Z, alpha, levels)
}

.ordinal_joint_response_mean <- function(fit, design, t = 1L, rows = NULL, cols = NULL) {
	state <- .joint_predictive_state(fit)
	base <- if (inherits(fit, "lame")) fit$Xlist[[t]] else fit[["X"]]
	if (is.null(rows)) rows <- seq_len(dim(design)[1L])
	if (is.null(cols)) cols <- seq_len(dim(design)[2L])
	if (!is.null(base) && !identical(dim(base)[1:2], dim(design)[1:2])) {
		base[rows, cols, ] <- design
		design <- base
	}
	out <- matrix(0, length(rows), length(cols))
	for (s in state$index) {
		eta <- .joint_predictor(fit, state, s, t, design)
		if (isTRUE(fit$symmetric)) eta <- eta / 2 + t(eta) / 2
		out <- out + .ordinal_response_mean(
			eta[rows, cols, drop = FALSE],
			.ordinal_saved_cutpoints(fit, s),
			fit$ordinal_levels
		)
	}
	out / length(state$index)
}
