# dynamic regression coefficients for lame(): parse `dynamic_beta`, label the coefficient blocks
# (intercept, dyad, row, col), set per-block AR(1) starting values and assemble per-period beta vectors.
# the FFBS sampler is in src/dynamic_beta.cpp; none of this runs when dynamic_beta = FALSE

# block label of each coefficient (intercept, dyad, row, col) from the names lame() gives beta: a `_row` /
# `_col` suffix (or `.row` / `.col`) or a matching Xrow / Xcol column name, else dyad. shared by the parser and
# the assembly so both see the same labels
.lame_coef_blocks <- function(
	coef_names,
	Xrow,
	Xcol,
	Xdyad,
	intercept,
	symmetric,
	bipartite
) {
	if (is.null(coef_names) || length(coef_names) == 0L) return(character(0))
	p <- length(coef_names)
	out <- rep("dyad", p)
	if (isTRUE(intercept)) {
		out[1] <- "intercept"
	}
	# row and column covariate names from the second dimnames of Xrow / Xcol
	row_nms <- if (!is.null(Xrow)) {
		dn <- dimnames(Xrow)
		if (!is.null(dn) && length(dn) >= 2L) dn[[2]] else NULL
	} else NULL
	col_nms <- if (!is.null(Xcol)) {
		dn <- dimnames(Xcol)
		if (!is.null(dn) && length(dn) >= 2L) dn[[2]] else NULL
	} else NULL
	# suffix match, or a bare name found in Xrow / Xcol
	is_row <- grepl("[\\._]row$", coef_names) |
		(!is.null(row_nms) & coef_names %in% row_nms)
	is_col <- grepl("[\\._]col$", coef_names) |
		(!is.null(col_nms) & coef_names %in% col_nms)
	out[is_row] <- "row"
	out[is_col] <- "col"
	# a symmetric fit names its nodal covariates `node`; they count as row coefficients, and H_b is tied to H_a
	if (symmetric) {
		is_node <- grepl("[\\._]node$|^node[\\._]", coef_names)
		out[is_node] <- "row"
	}
	# keep the intercept label even if its name matches a row or column covariate
	if (isTRUE(intercept)) out[1] <- "intercept"
	out
}

# parse the `dynamic_beta` argument: FALSE / NULL (nothing dynamic), TRUE (everything), a character vector of
# coefficient names or block shortcuts (intercept, dyad, row, col), 1-based indices, or a logical mask of length p.
# each block with a dynamic coefficient gets its own rho_beta and sigma_beta. returns any, mask, static_idx,
# dynamic_idx, groups ("" for static columns), group_names and n_groups
parse_dynamic_beta <- function(
	dynamic_beta,
	coef_names,
	coef_block,
	intercept,
	Tn,
	family,
	mode
) {
	p <- length(coef_names)
	stopifnot(length(coef_block) == p)
	# FALSE or NULL: nothing is dynamic
	if (is.null(dynamic_beta) ||
			(is.logical(dynamic_beta) && length(dynamic_beta) == 1L &&
					isFALSE(dynamic_beta))) {
		return(list(
			any = FALSE,
			mask = rep(FALSE, p),
			static_idx = seq_len(p),
			dynamic_idx = integer(0),
			groups = rep("", p),
			group_names = character(0),
			n_groups = 0L
		))
	}
	# the AR(1) prior needs at least two periods
	if (!is.null(Tn) && Tn < 2L) {
		cli::cli_abort(c(
			"{.arg dynamic_beta} requires at least 2 time periods.",
			"i" = "Got {Tn}; use {.fn ame} for a single cross-section or pass more periods."
		))
	}
	mask <- rep(FALSE, p)
	if (isTRUE(dynamic_beta)) {
		mask[] <- TRUE
	} else if (is.logical(dynamic_beta)) {
		if (length(dynamic_beta) != p) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} as a logical vector must have length {p} (one entry per coefficient).",
				"i" = "Got length {length(dynamic_beta)}.",
				"i" = "Coefficient layout: {.val {coef_names}}."
			))
		}
		mask <- dynamic_beta
	} else if (is.numeric(dynamic_beta)) {
		# indices must be integer-valued; fractional ones are refused, not truncated
		if (any(dynamic_beta != as.integer(dynamic_beta))) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} numeric indices must be integer-valued.",
				"i" = "Got: {.val {dynamic_beta}}.",
				"i" = "Use {.code as.integer(...)} explicitly if you mean to truncate."
			))
		}
		idx <- as.integer(dynamic_beta)
		if (any(idx < 1L) || any(idx > p)) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} integer indices must lie in 1..{p}.",
				"i" = "Got: {.val {dynamic_beta}}."
			))
		}
		mask[idx] <- TRUE
	} else if (is.character(dynamic_beta)) {
		# block shortcuts and coefficient names
		shortcuts <- c("intercept", "dyad", "row", "col")
		req <- unique(dynamic_beta)
		for (r in req) {
			if (r %in% shortcuts) {
				hit <- coef_block == r
				if (!any(hit)) {
					avail <- unique(coef_block)
					cli::cli_warn(c(
						"{.arg dynamic_beta} requested block {.val {r}} but no coefficient has that block label.",
						"i" = if (length(avail)) "Available blocks: {.val {avail}}." else "This fit has no regression coefficients."
					))
				}
				mask[hit] <- TRUE
			} else if (r %in% coef_names) {
				mask[match(r, coef_names)] <- TRUE
			} else {
				cli::cli_abort(c(
					"{.arg dynamic_beta} name {.val {r}} is not a coefficient or block shortcut.",
					"i" = "Recognised coefficient names: {.val {coef_names}}.",
					"i" = "Recognised block shortcuts: {.val {shortcuts}}."
				))
			}
		}
	} else {
		cli::cli_abort(c(
			"{.arg dynamic_beta} must be FALSE, TRUE, a character vector, an integer index, or a logical mask.",
			"i" = "Got class {.cls {class(dynamic_beta)[1]}}."
		))
	}
	# a request that selected no coefficient is an error
	if (!any(mask)) {
		# TRUE on a fit without coefficients is not an error
		non_trivial_request <-
			(is.character(dynamic_beta) && length(dynamic_beta) > 0L) ||
			(is.numeric(dynamic_beta) && length(dynamic_beta) > 0L) ||
			(is.logical(dynamic_beta) && length(dynamic_beta) == p && any(dynamic_beta))
		if (non_trivial_request) {
			# ordinal fits have no intercept to make dynamic
			if (family %in% c("ordinal") &&
					is.character(dynamic_beta) && "intercept" %in% dynamic_beta) {
				cli::cli_abort(c(
					"{.arg dynamic_beta} requested block {.val intercept} but family {.val {family}} has no identified intercept.",
					"i" = "Choose a different block (e.g. {.val dyad}) or pick a different family."
				))
			}
			cli::cli_abort(c(
				"{.arg dynamic_beta} = {.val {dynamic_beta}} selected zero coefficients.",
				"i" = "Available coefficient names: {.val {coef_names}}.",
				"i" = "Available block labels: {.val {unique(coef_block)}}.",
				"i" = "Pass {.arg dynamic_beta} = {.val FALSE} (the default) for a static fit."
			))
		}
		return(list(
			any = FALSE,
			mask = rep(FALSE, p),
			static_idx = seq_len(p),
			dynamic_idx = integer(0),
			groups = rep("", p),
			group_names = character(0),
			n_groups = 0L
		))
	}
	# ordinal has no identified intercept, so it cannot be dynamic (see ?lame)
	if (family %in% c("ordinal")) {
		if (any(mask & coef_block == "intercept")) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} cannot include the intercept for family {.val {family}}.",
				"i" = "An intercept is not identified for {.val {family}} (see ?lame).",
				"i" = "Drop {.val intercept} from {.arg dynamic_beta}."
			))
		}
	}
	groups <- ifelse(mask, coef_block, "")
	group_names <- sort(unique(groups[groups != ""]))
	list(
		any = TRUE,
		mask = mask,
		static_idx = which(!mask),
		dynamic_idx = which(mask),
		groups = groups,
		group_names = group_names,
		n_groups = length(group_names)
	)
}

# starting values of the per-block AR(1) parameters: rho_beta at 0.8 and sigma_beta at 0.25, set by
# prior$rho_beta_mean and prior$sigma_beta_init; they only affect early mixing
init_rho_beta <- function(group_names, prior = list()) {
	pm <- if (!is.null(prior$rho_beta_mean)) prior$rho_beta_mean else 0.8
	if (length(group_names) == 0L) return(numeric(0))
	out <- rep(pm, length(group_names))
	names(out) <- group_names
	out
}
init_sigma_beta <- function(group_names, prior = list()) {
	si <- if (!is.null(prior$sigma_beta_init)) prior$sigma_beta_init else 0.25
	if (length(group_names) == 0L) return(numeric(0))
	out <- rep(si, length(group_names))
	names(out) <- group_names
	out
}

# per-period innovation scales Lambda of the coefficient AR(1), with Q_t = sigma_beta^2 * Lambda: a diagonal
# g-prior analogue, the mean number of dyads per period over the period-averaged diagonal of X'X
build_beta_state_scales <- function(
	beta_dyn,
	XtX_per_t,
	n_dyads_per_t
) {
	# returns one p_dyn x p_dyn matrix per period
	dyn_idx <- beta_dyn$dynamic_idx
	p_dyn <- length(dyn_idx)
	Tn <- length(XtX_per_t)
	out <- vector("list", Tn)
	# average the dynamic block of X'X over the periods with a finite design
	XtX_dyn_avg <- matrix(0, p_dyn, p_dyn)
	have <- 0L
	for (t in seq_len(Tn)) {
		Xt <- XtX_per_t[[t]]
		if (is.null(Xt) || any(!is.finite(Xt))) next
		XtX_dyn_avg <- XtX_dyn_avg + Xt[dyn_idx, dyn_idx, drop = FALSE]
		have <- have + 1L
	}
	if (have > 0L) XtX_dyn_avg <- XtX_dyn_avg / have
	# a diagonal Lambda keeps innovations independent across coefficients; an all-zero column gets diagonal 1
	d_avg <- diag(XtX_dyn_avg)
	d_avg[!is.finite(d_avg) | d_avg <= 0] <- 1
	n_d_mean <- mean(unlist(n_dyads_per_t), na.rm = TRUE)
	if (!is.finite(n_d_mean) || n_d_mean <= 0) n_d_mean <- 1
	scale_vec <- n_d_mean / d_avg
	# no upper cap: the scale is in the units of each covariate, so an absolute cap would bind for small units only
	Lambda <- diag(scale_vec, p_dyn)
	# the same Lambda in every period; sigma_beta is shared within a block
	for (t in seq_len(Tn)) out[[t]] <- Lambda
	out
}

# inverse of the diagonal innovation scale, with no condition-number cutoff across coefficient units
.lame_beta_scale_precision <- function(Lambda) {
	if (!is.matrix(Lambda) || !is.numeric(Lambda) || nrow(Lambda) != ncol(Lambda) ||
			any(!is.finite(Lambda)) || any(Lambda[row(Lambda) != col(Lambda)] != 0) ||
			any(diag(Lambda) <= 0)) {
		stop("Coefficient innovation scale must be a finite positive diagonal matrix.")
	}
	precision <- 1 / diag(Lambda)
	if (any(!is.finite(precision))) {
		stop("Coefficient innovation precision is outside numerical range; rescale the covariates.")
	}
	diag(precision, nrow(Lambda))
}

# the p-length beta of one period from its dynamic and static parts, used in the MCMC loop and at fit assembly
assemble_beta_path <- function(BETA_dyn_t, beta_static, dyn_idx, static_idx, p) {
	out <- numeric(p)
	if (length(dyn_idx) > 0L) out[dyn_idx] <- as.numeric(BETA_dyn_t)
	if (length(static_idx) > 0L) out[static_idx] <- as.numeric(beta_static)
	out
}

# sum-to-zero bases for the additive effects when the intercept is dynamic: a_i + intercept_t is unidentified
# per period, so a = H_a alpha with H_a an orthonormal basis of the complement of 1_n, and b likewise unless
# symmetric. returns H_a, H_b (NULL when not needed) and `active`
build_a_b_constraint_bases <- function(
	beta_dyn,
	coef_block,
	n_a,
	n_b,
	bipartite,
	symmetric,
	intercept = FALSE
) {
	need_constrain_a <- FALSE
	need_constrain_b <- FALSE
	if (beta_dyn$any && intercept) {
		need_constrain_a <- TRUE
		if (!symmetric) need_constrain_b <- TRUE
	}
	# a nodal coefficient without an intercept is identified by its proper prior and the actor prior, so it
	# needs no constraint

	make_contrast <- function(n) {
		# columns of H span the orthogonal complement of 1_n
		if (n < 2L) return(matrix(0, n, 0))
		I_n <- diag(n)
		ones <- matrix(1, n, n) / n
		P <- I_n - ones
		eg <- eigen(P, symmetric = TRUE)
		# eigenvalue 1 with multiplicity n - 1, eigenvalue 0 once
		keep <- order(eg$values, decreasing = TRUE)[seq_len(n - 1L)]
		eg$vectors[, keep, drop = FALSE]
	}
	H_a <- if (need_constrain_a) make_contrast(n_a) else NULL
	H_b <- if (need_constrain_b) make_contrast(n_b) else NULL
	# a symmetric fit has a == b, so H_b is H_a
	if (symmetric && !is.null(H_a)) H_b <- H_a
	list(
		H_a = H_a,
		H_b = H_b,
		active = !is.null(H_a) || !is.null(H_b)
	)
}
