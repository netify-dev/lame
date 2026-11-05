.lame_explicit_ordinal_support <- function(family, cutpoints, bip, symmetric, dcor, rho = NULL) {
	if (!identical(family, "ordinal") || !identical(cutpoints, "explicit")) return(invisible(NULL))
	if (isTRUE(bip)) {
		cli::cli_abort(
			c(
				"Explicit ordinal cutpoints are not implemented for bipartite fits.",
				"i" = "Use ordinal_cutpoints = 'data_induced' for the existing bipartite rank model."
			),
			class = "lame_explicit_ordinal_mode"
		)
	}
	if (!isTRUE(symmetric) && !isTRUE(dcor) &&
			(!is.null(rho) && any(!is.finite(rho) | rho != 0))) {
		cli::cli_abort(
			c(
				"Explicit ordinal dcor = FALSE requires zero directed-dyad correlation.",
				"i" = "Use dcor = TRUE for the correlated explicit-cutpoint model, or initialize rho at zero."
			),
			class = "lame_explicit_ordinal_correlation"
		)
	}
	invisible(NULL)
}
