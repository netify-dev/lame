# the diagonal innovation geometry of the regular-grid coefficient sampler, recovered from
# the saved designs (expanded symmetric layout)
.lame_beta_transition_scale <- function(object) {
	B <- object$BETA
	names_beta <- dimnames(B)[[2L]]
	mask <- object$beta_dynamic_mask
	if (length(dim(B)) != 3L || is.null(names_beta) || anyDuplicated(names_beta) ||
			!is.logical(mask) || length(mask) != length(names_beta) || anyNA(mask)) {
		cli::cli_abort(
			"Coefficient transitions require named paths and their dynamic mask.",
			class = "lame_beta_transition"
		)
	}
	groups <- object$beta_dynamic_groups
	R <- object$RHO_BETA; V <- object$SIGMA_BETA
	if (length(groups) != length(names_beta) || anyNA(groups[mask]) ||
			!is.matrix(R) || !is.matrix(V) || nrow(R) != dim(B)[1L] || nrow(V) != dim(B)[1L] ||
			!all(groups[mask] %in% colnames(R)) || !all(groups[mask] %in% colnames(V)) ||
			any(!is.finite(R)) || any(!is.finite(V)) || any(abs(R) > 1) || any(V < 0)) {
		cli::cli_abort(
			"Coefficient paths and transition draws must be finite, nonnegative-scaled and aligned.",
			class = "lame_beta_transition"
		)
	}
	kind <- object$dynamic_beta_kind %||% "ar1"
	if (!kind %in% c("ar1", "rw1") || !is.null(object$dynamic_beta_per_actor)) {
		cli::cli_abort(
			"This diagnostic or forecast is not implemented for this coefficient transition model.",
			class = "lame_beta_transition"
		)
	}
	ti <- object$time_index
	if (!is.null(ti) && (!is.numeric(ti) || length(ti) != dim(B)[3L] ||
			any(!is.finite(ti)) || any(diff(ti) <= 0) ||
			length(unique(diff(ti))) > 1L)) {
		cli::cli_abort(
			"This diagnostic or forecast requires equally spaced coefficient periods.",
			class = "lame_beta_transition"
		)
	}
	X <- object$Xlist
	if (!is.list(X) || length(X) != dim(B)[3L] || !length(X)) {
		cli::cli_abort(
			"Saved training designs are required to recover coefficient innovation scales.",
			class = "lame_beta_transition"
		)
	}
	xn <- dimnames(X[[1L]])[[3L]]
	if (is.null(xn) || anyDuplicated(xn) || any(
		!vapply(
			X,
			function(x)
				is.numeric(x) && length(dim(x)) == 3L && identical(dimnames(x)[[3L]], xn),
			logical(1)
		)
	)) {
		cli::cli_abort(
			"Saved training designs must have consistent named coefficient columns.",
			class = "lame_beta_transition"
		)
	}
	stored <- match(xn, names_beta)
	if (isTRUE(object$symmetric)) {
		collapsed <- gsub("row|col", "node", xn)
		stored[is.na(stored)] <- match(collapsed[is.na(stored)], names_beta)
	}
	if (anyNA(stored) || !all(seq_along(names_beta) %in% stored)) {
		cli::cli_abort(
			"The saved design cannot be matched to the stored coefficient paths.",
			class = "lame_beta_transition"
		)
	}
	expanded <- which(mask[stored])
	counts <- tabulate(stored, length(names_beta))
	if (any(counts > if (isTRUE(object$symmetric)) 2L else 1L)) {
		cli::cli_abort("The saved coefficient layout is ambiguous.", class = "lame_beta_transition")
	}
	XtX <- lapply(X, function(x) crossprod(matrix(x, prod(dim(x)[1:2]), dim(x)[3L])))
	ncells <- lapply(X, function(x) prod(dim(x)[1:2]))
	lambda <- build_beta_state_scales(
		list(dynamic_idx = expanded),
		XtX,
		ncells
	)[[1L]]
	variance <- numeric(length(names_beta))
	for (j in seq_along(expanded)) {
		k <- stored[expanded[j]]
		variance[k] <- variance[k] + lambda[j, j] / counts[k]^2
	}
	if (any(!is.finite(variance)) || any(variance[mask] <= 0)) {
		cli::cli_abort("The fitted coefficient innovation scales are invalid.", class = "lame_beta_transition")
	}
	list(sd = stats::setNames(sqrt(variance), names_beta))
}
