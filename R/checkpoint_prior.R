# resolved prior from a checkpoint; the caller-supplied settings are kept so
# training-window refits do not inherit defaults estimated on the full panel
.lame_checkpoint_prior <- function(checkpoint, envir) {
	if ("resolved_prior" %in% names(checkpoint)) {
		prior <- checkpoint[["resolved_prior"]]
		input <- checkpoint[["requested_prior"]]
		if (!is.list(prior) || (!is.null(input) && !is.list(input))) {
			cli::cli_abort(
				"The checkpoint carries an invalid prior specification.",
				class = "lame_checkpoint_prior"
			)
		}
		attr(prior, "lame_checkpoint_input_prior") <- input %||% list()
		return(prior)
	}
	# checkpoints without a stored prior: Gaussian bipartite AR(1) fits take
	# these rho_uv defaults unless the saved call supplied them
	layout <- checkpoint$start_vals[["checkpoint_state"]][["layout"]]
	if (!is.list(layout) || !identical(layout$family, "normal") ||
			!identical(layout$mode, "bipartite") || !isTRUE(layout$dynamic_uv) ||
			!identical(layout$dynamic_uv_kind, "ar1")) return(NULL)
	input <- checkpoint$call_args$prior
	input <- if (is.null(input)) list() else if (is.list(input)) input else
		eval(input, envir = envir)
	if (!is.list(input)) {
		cli::cli_abort(
			"The saved checkpoint prior must evaluate to a list.",
			class = "lame_checkpoint_prior"
		)
	}
	prior <- input
	if (is.null(prior$rho_uv_mean)) prior$rho_uv_mean <- .9
	if (is.null(prior$rho_uv_sd)) prior$rho_uv_sd <- .1
	attr(prior, "lame_checkpoint_input_prior") <- input
	prior
}
