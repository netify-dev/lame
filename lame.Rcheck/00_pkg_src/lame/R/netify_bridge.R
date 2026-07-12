.lame_is_netify <- function(x) {
	inherits(x, "netify")
}

.lame_prepare_netify_input <- function(Y, Xdyad = NULL, Xrow = NULL, Xcol = NULL,
                                       mode, family,
                                       mode_missing, family_missing,
                                       longitudinal,
                                       caller = "lame",
                                       fit_method = c("gibbs", "als"),
                                       symmetric = FALSE,
                                       symmetric_missing = TRUE) {
	`%||%` <- function(x, y) if (is.null(x)) y else x
	fit_method <- match.arg(fit_method)
	if (!.lame_is_netify(Y)) {
		return(list(converted = FALSE))
	}
	if (!is.null(Xdyad) || !is.null(Xrow) || !is.null(Xcol)) {
		cli::cli_abort(c(
			"When {.arg Y} is a {.cls netify} object, pass covariates through {.pkg netify}.",
			"i" = "Use {.fn netify::netify} with {.arg dyad_vars} and {.arg nodal_vars}, or convert first with {.fn netify::to_lame} and pass the returned pieces explicitly."))
	}
	if (!requireNamespace("netify", quietly = TRUE)) {
		cli::cli_abort(c(
			"Package {.pkg netify} is required for {.cls netify} inputs.",
			"i" = "Install {.pkg netify} or pass ordinary {.pkg lame} matrices, arrays, and lists."))
	}

	netify_type <- attr(Y, "netify_type")
	is_longit <- identical(netify_type, "longit_array") ||
		identical(netify_type, "longit_list")
	if (isTRUE(longitudinal) && !is_longit) {
		cli::cli_abort(c(
			"{.fn {caller}} expects longitudinal data, but this {.cls netify} object is {.val {netify_type %||% 'cross_sec'}}.",
			"i" = "Use {.fn ame} or {.fn ame_als} for a cross-sectional netify object."))
	}
	if (!isTRUE(longitudinal) && is_longit) {
		cli::cli_abort(c(
			"{.fn {caller}} expects one network, but this {.cls netify} object is longitudinal.",
			"i" = "Use {.fn lame} or {.fn lame_als} for longitudinal netify objects."))
	}

	mode_arg <- if (isTRUE(mode_missing)) NULL else {
		match.arg(mode, c("unipartite", "bipartite"))
	}
	family_arg <- if (isTRUE(family_missing)) NULL else family
	# netify stores whether the network is undirected on the object; honor
	# it so a symmetric netify object is not silently fit as directed
	netify_symmetric <- isTRUE(attr(Y, "symmetric"))
	if (!isTRUE(symmetric_missing) &&
	    !identical(isTRUE(symmetric), netify_symmetric)) {
		cli::cli_abort(c(
			"{.arg symmetric} = {.val {isTRUE(symmetric)}} conflicts with the {.cls netify} object (symmetric = {.val {netify_symmetric}}).",
			"i" = "Drop {.arg symmetric} to use the value stored on the netify object, or rebuild the netify object with the intended symmetry."))
	}
	nl <- netify::to_lame(
		Y,
		lame = isTRUE(longitudinal),
		family = family_arg,
		pad = FALSE,
		fit_method = fit_method
	)
	# multilayer netify objects convert to one result per layer; there is
	# no single Y to fit, so stop with a diagnosis that names the layers
	if (is.null(nl$Y) && is.list(nl) && length(nl) > 0L &&
	    all(vapply(nl, function(x) is.list(x) && !is.null(x$Y), logical(1)))) {
		cli::cli_abort(c(
			"Multilayer {.cls netify} objects cannot be fit directly ({length(nl)} layers: {.val {names(nl)}}).",
			"i" = "Extract one layer with {.fn netify::subset_netify} and fit each layer separately."))
	}
	netify_mode <- nl$mode %||% if (identical(attr(Y, "mode"), "bipartite")) {
		"bipartite"
	} else {
		"unipartite"
	}
	netify_family <- nl$family %||% family

	if (!isTRUE(mode_missing) && !identical(mode_arg, netify_mode)) {
		cli::cli_abort(c(
			"{.arg mode} = {.val {mode_arg}} conflicts with the {.cls netify} object mode {.val {netify_mode}}.",
			"i" = "Drop {.arg mode} or rebuild the netify object with the intended mode."))
	}

	Y_out <- nl$Y
	if (isTRUE(longitudinal)) {
		if (!is.list(Y_out) || is.array(Y_out)) {
			cli::cli_abort(c(
				"{.fn netify::to_lame} did not return the longitudinal list format expected by {.fn {caller}}.",
				"i" = "This usually means the netify object was not built as a longitudinal panel."))
		}
		if (identical(netify_mode, "bipartite")) {
			has_names <- vapply(Y_out, function(yt) {
				!is.null(rownames(yt)) && !is.null(colnames(yt))
			}, logical(1))
			if (!all(has_names)) {
				cli::cli_abort(c(
					"Longitudinal bipartite netify inputs must carry row and column names on every slice.",
					"i" = "These names let {.pkg lame} align changing row and column actor sets separately."))
			}
		}
	} else if (!is.matrix(Y_out)) {
		cli::cli_abort(c(
			"{.fn netify::to_lame} did not return the matrix format expected by {.fn {caller}}.",
			"i" = "Use {.fn lame} or {.fn lame_als} if this is a longitudinal netify object."))
	}

	family_resolved <- if (isTRUE(family_missing)) netify_family else family

	# when the family was inferred (not user-supplied) and came back "normal",
	# but the weights look like counts, point the user at poisson/ordinal
	if (isTRUE(family_missing) && identical(family_resolved, "normal")) {
		v <- unlist(lapply(if (is.list(Y_out)) Y_out else list(Y_out),
		                   function(m) m[is.finite(m)]))
		if (length(v) > 0L && all(v >= 0) && all(v %% 1 == 0) &&
		    length(unique(v)) > 2L) {
			cli::cli_inform(c(
				"i" = "{.pkg netify} inferred {.arg family} = {.val normal}, but {.arg Y} looks like counts (non-negative integers with more than 2 distinct values).",
				"i" = "Pass {.code family = \"poisson\"} (or {.code family = \"ordinal\"} for bounded ratings) if that is what the weights are."))
		}
	}

	list(
		converted = TRUE,
		Y = Y_out,
		Xdyad = nl$Xdyad,
		Xrow = nl$Xrow,
		Xcol = nl$Xcol,
		mode = if (isTRUE(mode_missing)) netify_mode else mode_arg,
		family = family_resolved,
		symmetric = netify_symmetric
	)
}
