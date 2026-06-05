# dynamic regression coefficients for lame.
#
# % %||% (r 4.4+ has this as a base operator; keep a local fallback so the
# package keeps installing on r 4.3 which still ships in stable distros).
`%||%` <- function(x, y) if (is.null(x)) y else x
#
# this file ships the un-exported r-side scaffolding that lame() calls when
# the user passes `dynamic_beta`. nothing here is user-facing: the parser
# decides which coefficient columns become ar(1)-evolving over time, builds
# block-level grouping (intercept/dyad/row/col), initialises rho_beta and
# sigma_beta hyperparameters per block, and assembles per-time-period beta
# vectors from the dynamic-path posterior. the actual ffbs sampler lives in
# src/dynamic_beta.cpp.
#
# the byte-identical-default contract: when `dynamic_beta = false` (the
# default), none of these helpers is ever called and lame() takes its
# existing rbeta_ab_rep_fc / rbeta_ab_bip_gibbs_cpp path. that's the
# regression-test guarantee that lives in tests/testthat/test_dynamic_beta_backcompat.r.

####
# infer the per-coefficient block label ("intercept", "dyad", "row", "col")
# from the colnames(beta) layout that lame() emits. lame()'s design
# construction labels them via the convention:
#   - "intercept" -> the leading column when intercept = true
#   - "<name>_row" / "<name>" matching a colnames(xrow[,,1]) -> "row"
#   - "<name>_col" / "<name>" matching a colnames(xcol[,,1]) -> "col"
#   - otherwise -> "dyad" (xdyad slice names, suffix "_dyad")
# the dot-form suffixes (".row", ".col", ".dyad") are recognised too.
# this helper is a single source of truth so the parser and the assembly
# downstream see the same block labels.
.lame_coef_blocks <- function(coef_names, Xrow, Xcol, Xdyad, intercept,
                              symmetric, bipartite) {
	if (is.null(coef_names) || length(coef_names) == 0L) return(character(0))
	p <- length(coef_names)
	out <- rep("dyad", p)  # default
	if (isTRUE(intercept)) {
		# first column carries the intercept
		out[1] <- "intercept"
	}
	# row covariate names: xrow is a 3-d array [n, p_row, t] for unipartite,
	# or a list for the bipartite-input path that's already been converted.
	row_nms <- if (!is.null(Xrow)) {
		dn <- dimnames(Xrow)
		if (!is.null(dn) && length(dn) >= 2L) dn[[2]] else NULL
	} else NULL
	col_nms <- if (!is.null(Xcol)) {
		dn <- dimnames(Xcol)
		if (!is.null(dn) && length(dn) >= 2L) dn[[2]] else NULL
	} else NULL
	# lame() appends "_row" / "_col" suffixes (".row" / ".col" also
	# accepted). match either underscore or dot form, then fall back to
	# bare-name matches.
	is_row <- grepl("[\\._]row$", coef_names) |
	          (!is.null(row_nms) & coef_names %in% row_nms)
	is_col <- grepl("[\\._]col$", coef_names) |
	          (!is.null(col_nms) & coef_names %in% col_nms)
	out[is_row] <- "row"
	out[is_col] <- "col"
	# symmetric case: the colnames are merged ("node" instead of "row"/"col");
	# treat them as "row" coefficients (which feeds into the constraint logic
	# correctly because symmetric implies a == b and the constraint helper
	# ties h_b to h_a).
	if (symmetric) {
		is_node <- grepl("[\\._]node$|^node[\\._]", coef_names)
		out[is_node] <- "row"
	}
	# preserve "intercept" label even if "intercept" appears among row/col
	# matches (defensive: in practice colnames(beta)[1] is "intercept" and
	# never matches a xrow/xcol column name)
	if (isTRUE(intercept)) out[1] <- "intercept"
	out
}

####
# parse the `dynamic_beta` argument and the coefficient layout.
#
# `dynamic_beta` may be:
#   false / null          -> nothing is dynamic (the default)
#   true                  -> every coefficient is dynamic
#   character vector      -> the named coefficients become dynamic; allowed
#                              shortcuts: "intercept", "dyad", "row", "col",
#                              or any actual column name from `coef_names`
#   integer vector        -> 1-based column indices into the beta layout
#   logical vector len p  -> per-coefficient mask
#
# `coef_block` is the per-column block label ("intercept", "dyad", "row",
# "col") that the lame() design step already computes when it assembles x;
# we use it to (a) translate shortcuts and (b) build the `groups` vector
# that determines which rho_beta / sigma_beta each coefficient is governed
# by. one rho_beta and one sigma_beta is sampled per distinct block that
# has at least one dynamic coefficient -- all dyadic dynamic coefficients
# share a rho, etc.
#
# returns a list with: any (logical scalar), mask (length-p logical),
# static_idx, dynamic_idx (integer vectors, 1-based), groups (length-p
# character; "" for static columns), group_names (the distinct dynamic
# block labels), n_groups (length(group_names)).
parse_dynamic_beta <- function(dynamic_beta, coef_names, coef_block,
                               intercept, T, family, mode) {
	p <- length(coef_names)
	stopifnot(length(coef_block) == p)
	# the empty / false / null case wires through as "nothing is dynamic"
	if (is.null(dynamic_beta) ||
	    (is.logical(dynamic_beta) && length(dynamic_beta) == 1L &&
	     isFALSE(dynamic_beta))) {
		return(list(any = FALSE,
		            mask = rep(FALSE, p),
		            static_idx = seq_len(p),
		            dynamic_idx = integer(0),
		            groups = rep("", p),
		            group_names = character(0),
		            n_groups = 0L))
	}
	# ar(1) prior is unidentified with a single time slice
	if (!is.null(T) && T < 2L) {
		cli::cli_abort(c(
			"{.arg dynamic_beta} requires at least 2 time periods.",
			"i" = "Got T = {T}; use {.fn ame} for a single cross-section or pass more periods."))
	}
	# resolve the mask
	mask <- rep(FALSE, p)
	if (isTRUE(dynamic_beta)) {
		mask[] <- TRUE
	} else if (is.logical(dynamic_beta)) {
		if (length(dynamic_beta) != p) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} as a logical vector must have length {p} (one entry per coefficient).",
				"i" = "Got length {length(dynamic_beta)}.",
				"i" = "Coefficient layout: {.val {coef_names}}."))
		}
		mask <- dynamic_beta
	} else if (is.numeric(dynamic_beta)) {
		# require integer-valued numeric; abort rather than silently
		# truncate fractional indices like 1.5 to 1.
		if (any(dynamic_beta != as.integer(dynamic_beta))) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} numeric indices must be integer-valued.",
				"i" = "Got: {.val {dynamic_beta}}.",
				"i" = "Use {.code as.integer(...)} explicitly if you mean to truncate."))
		}
		idx <- as.integer(dynamic_beta)
		if (any(idx < 1L) || any(idx > p)) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} integer indices must lie in 1..{p}.",
				"i" = "Got: {.val {dynamic_beta}}."))
		}
		mask[idx] <- TRUE
	} else if (is.character(dynamic_beta)) {
		# expand shortcuts and per-name lookups
		shortcuts <- c("intercept", "dyad", "row", "col")
		req <- unique(dynamic_beta)
		for (r in req) {
			if (r %in% shortcuts) {
				hit <- coef_block == r
				if (!any(hit)) {
					cli::cli_warn(c(
						"{.arg dynamic_beta} requested block {.val {r}} but no coefficient has that block label.",
						"i" = "Available blocks: {.val {unique(coef_block)}}."))
				}
				mask[hit] <- TRUE
			} else if (r %in% coef_names) {
				mask[match(r, coef_names)] <- TRUE
			} else {
				cli::cli_abort(c(
					"{.arg dynamic_beta} name {.val {r}} is not a coefficient or block shortcut.",
					"i" = "Recognised coefficient names: {.val {coef_names}}.",
					"i" = "Recognised block shortcuts: {.val {shortcuts}}."))
			}
		}
	} else {
		cli::cli_abort(c(
			"{.arg dynamic_beta} must be FALSE, TRUE, a character vector, an integer index, or a logical mask.",
			"i" = "Got class {.cls {class(dynamic_beta)[1]}}."))
	}
	# nothing was actually selected: the user explicitly asked for something
	# non-false but no coefficient ended up flagged. abort with a clear
	# message naming the offending request rather than silently falling
	# back to a static fit.
	if (!any(mask)) {
		# distinguish a character/integer/logical request from a true / false / null
		non_trivial_request <-
			(is.character(dynamic_beta) && length(dynamic_beta) > 0L) ||
			(is.numeric(dynamic_beta)   && length(dynamic_beta) > 0L) ||
			(is.logical(dynamic_beta)   && length(dynamic_beta) == p && any(dynamic_beta))
		if (non_trivial_request) {
			# special case: family-level prohibition is what zeroed out the mask
			if (family %in% c("ordinal", "rrl") &&
			    is.character(dynamic_beta) && "intercept" %in% dynamic_beta) {
				cli::cli_abort(c(
					"{.arg dynamic_beta} requested block {.val intercept} but family {.val {family}} has no identified intercept.",
					"i" = "Choose a different block (e.g. {.val dyad}) or pick a different family."))
			}
			cli::cli_abort(c(
				"{.arg dynamic_beta} = {.val {dynamic_beta}} selected zero coefficients.",
				"i" = "Available coefficient names: {.val {coef_names}}.",
				"i" = "Available block labels: {.val {unique(coef_block)}}.",
				"i" = "Pass {.arg dynamic_beta} = {.val FALSE} (the default) for a static fit."))
		}
		return(list(any = FALSE,
		            mask = rep(FALSE, p),
		            static_idx = seq_len(p),
		            dynamic_idx = integer(0),
		            groups = rep("", p),
		            group_names = character(0),
		            n_groups = 0L))
	}
	# family-level prohibitions matching ?lame: ordinal and rrl have no
	# identified intercept, so a time-varying intercept is meaningless. the
	# downstream nodal-coef-dynamic-on-rrl ban is the same.
	if (family %in% c("ordinal", "rrl")) {
		if (any(mask & coef_block == "intercept")) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} cannot include the intercept for family {.val {family}}.",
				"i" = "An intercept is not identified for {.val {family}} (see ?lame).",
				"i" = "Drop {.val intercept} from {.arg dynamic_beta}."))
		}
		if (family == "rrl" && any(mask & coef_block == "row")) {
			cli::cli_abort(c(
				"{.arg dynamic_beta} cannot include row-nodal coefficients for {.val rrl}.",
				"i" = "Row-nodal coefficients are not identified for {.val rrl} (see ?lame)."))
		}
	}
	groups <- ifelse(mask, coef_block, "")
	group_names <- sort(unique(groups[groups != ""]))
	list(any = TRUE,
	     mask = mask,
	     static_idx = which(!mask),
	     dynamic_idx = which(mask),
	     groups = groups,
	     group_names = group_names,
	     n_groups = length(group_names))
}

####
# default initial values for the per-block ar(1) parameters.
#
# rho_beta is initialised at the prior mean (0.8); sigma_beta is started at a
# slightly larger innovation (0.25) than dynamic_ab's 0.1 because coefficients
# typically live on a wider scale than additive effects. both can be overridden
# via `prior$rho_beta_mean` / `prior$sigma_beta_init`. these are prior
# defaults -- after burn-in the sampler is driving them, so the choice only
# matters for transient mixing.
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

####
# build the per-coefficient state-prior scale lambda_b that defines the
# innovation covariance for the ar(1) on beta_b,t. the default is a
# diagonal lambda equal to the inverse of the coefficient-wise design
# precision (the g-prior analogue), keeping the scale of ar(1) innovations
# commensurate with the data-level information about each coefficient. for
# an intercept-only block (1 col) lambda_b is a 1x1 matrix equal to the
# median per-period inverse-xtx entry.
#
# `xtx_per_t` is a length-t list of (p x p) matrices, one per period; the
# parser supplies the dynamic indices so the lambda blocks are built only
# for the dynamic block columns. when a coefficient never varies across
# periods (e.g. a constant covariate after centering) the scale falls back
# to 1 to avoid a divide-by-zero blow-up.
build_beta_state_scales <- function(beta_dyn, XtX_per_t, n_dyads_per_t,
                                    g = NULL, beta0_var = 10) {
	# returns a length-t list of (p_dyn x p_dyn) lambda matrices; the ffbs
	# uses q_t = sigma^2 * lambda. a global `g` scales the lambda matrices;
	# otherwise the n-dyads-based default (the static path's g-prior choice).
	# the per-coefficient scale upper bound is tied to beta0_var (specifically
	# 100 * beta0_var) so a more diffuse / less diffuse prior moves the cap
	# via prior$beta0_var.
	dyn_idx <- beta_dyn$dynamic_idx
	p_dyn   <- length(dyn_idx)
	T       <- length(XtX_per_t)
	out     <- vector("list", T)
	# stack and average the dynamic-block sub-xtx across periods so a
	# period with all-zero design (rare but possible) still gets a sane
	# scale matrix from its neighbours
	XtX_dyn_avg <- matrix(0, p_dyn, p_dyn)
	have <- 0L
	for (t in seq_len(T)) {
		Xt <- XtX_per_t[[t]]
		if (is.null(Xt) || any(!is.finite(Xt))) next
		XtX_dyn_avg <- XtX_dyn_avg + Xt[dyn_idx, dyn_idx, drop = FALSE]
		have <- have + 1L
	}
	if (have > 0L) XtX_dyn_avg <- XtX_dyn_avg / have
	# diagonalise: a diagonal lambda avoids cross-coefficient leakage of
	# the ar(1) innovation and matches the ruv_dynamic_fc_cpp prior pattern.
	# jitter prevents infinite scale on an all-zero column.
	d_avg <- diag(XtX_dyn_avg)
	d_avg[!is.finite(d_avg) | d_avg <= 0] <- 1
	# g-prior style scaling: scale = mean(n_dyads_per_t)/d_avg, fallback
	# to a fixed 1 if d_avg is tiny
	n_d_mean <- mean(unlist(n_dyads_per_t), na.rm = TRUE)
	if (!is.finite(n_d_mean) || n_d_mean <= 0) n_d_mean <- 1
	scale_vec <- n_d_mean / d_avg
	# enforce a sane upper bound: an extreme outlier scale destabilises ffbs.
	# tie the cap to beta0_var (100 * beta0_var) so the prior's diffuseness
	# scales with the user's choice. defaults: beta0_var=10 -> cap=1000 (was
	# 1e4 fixed). users with very-low-precision designs can raise beta0_var.
	scale_cap <- max(100 * beta0_var, 100)
	scale_vec <- pmin(scale_vec, scale_cap)
	Lambda <- diag(scale_vec, p_dyn)
	# use the same lambda for every period (the t-varying part of q_t comes
	# from sigma_beta which is time-shared across periods within a block)
	for (t in seq_len(T)) out[[t]] <- Lambda
	out
}

####
# assemble the per-period beta vector from the dynamic path + static block.
# `beta_path_t` is the dynamic-only piece of length p_dyn; `beta_static` is
# the static-only piece of length p_static; `dyn_idx` / `static_idx` tell us
# where each goes in the final p-length beta. used both inside the mcmc
# loop (to feed z|beta updates) and at fit-assembly time (to expand a
# posterior-mean dynamic path to a full [p x t] table).
assemble_beta_path <- function(BETA_dyn_t, beta_static, dyn_idx, static_idx, p) {
	out <- numeric(p)
	if (length(dyn_idx) > 0L)   out[dyn_idx]   <- as.numeric(BETA_dyn_t)
	if (length(static_idx) > 0L) out[static_idx] <- as.numeric(beta_static)
	out
}

####
# constraint basis for the additive-effects (a, b) sampler when the
# intercept or a nodal coefficient is dynamic.
#
# the identifiability problem: a_i + intercept_t is unidentified for any
# constant shift c since (a_i + c) + (intercept_t - c) gives the same
# fitted value. with a static intercept the lame sampler already takes
# care of this by summing-to-zero a; with a dynamic intercept the shift is
# per-period and the standard sum-to-zero on a doesn't pin it.
#
# the cleanest fix is the contrast basis: we sample a in the n-1 dimensional
# orthogonal complement of 1_n. that is, instead of sampling a_i directly,
# we sample alpha_i in r^{n-1} and set a = h * alpha where h is an n x (n-1)
# matrix whose columns are orthogonal to 1_n. then sum(a) = 1' h alpha = 0
# automatically. the same trick handles dynamic row-nodal coefficients (the
# total effect "intercept + row_coef * row_x" can absorb a constant per
# period unless we constrain a). bipartite: only the row block needs
# constraining if intercept-or-row-coef is dynamic, only the column block if
# intercept-or-col-coef is dynamic.
#
# this builder returns a list with $h_a, $h_b (or null when no constraint is
# needed). the c++ ffbs doesn't see these directly -- they live on the r
# side because the a/b sampler is the unipartite rbeta_ab_rep_fc_cpp pathway
# and the bipartite rbeta_ab_bip_gibbs_cpp pathway, both of which are easier
# to wrap on the r side. when both h_a and h_b are null nothing changes.
build_a_b_constraint_bases <- function(beta_dyn, coef_block, n_a, n_b,
                                       bipartite, symmetric,
                                       intercept = FALSE) {
	# identify which blocks need an a or b constraint.
	# the universal rule under dynamic_beta: when an intercept is in the
	# model, (intercept, a, b) is rank-deficient without a sum-to-zero
	# constraint on a (and b for asymmetric models). the static lame path
	# implicitly handles this via its joint beta-a-b marginalised sampler;
	# the dynamic-beta path splits beta from (a, b) and so must impose the
	# constraint explicitly.
	need_constrain_a <- FALSE
	need_constrain_b <- FALSE
	if (beta_dyn$any && intercept) {
		# intercept present -> always sum-to-zero a (and b unless symmetric)
		need_constrain_a <- TRUE
		if (!symmetric) need_constrain_b <- TRUE
	}
	if (beta_dyn$any) {
		dyn_blocks <- unique(coef_block[beta_dyn$mask])
		# intercept-dynamic -> still constrain (same rule, redundant)
		if ("intercept" %in% dyn_blocks) {
			need_constrain_a <- TRUE
			if (!symmetric) need_constrain_b <- TRUE
		}
		# row-coef-dynamic -> constrain a (row effects)
		if ("row" %in% dyn_blocks) need_constrain_a <- TRUE
		# col-coef-dynamic -> constrain b (col effects)
		if ("col" %in% dyn_blocks) need_constrain_b <- TRUE
			# dyadic-only dynamic without intercept needs no a/b constraint
	}
	make_contrast <- function(n) {
		# columns of h span the orthogonal complement of 1_n (helmert-style)
		if (n < 2L) return(matrix(0, n, 0))
		I_n <- diag(n)
		ones <- matrix(1, n, n) / n
		P <- I_n - ones
		# columns of the (n-1) leading eigenvectors of p
		eg <- eigen(P, symmetric = TRUE)
		# eigenvalues 1 with multiplicity n-1, eigenvalue 0 once
		keep <- order(eg$values, decreasing = TRUE)[seq_len(n - 1L)]
		eg$vectors[, keep, drop = FALSE]
	}
	H_a <- if (need_constrain_a) make_contrast(n_a) else NULL
	H_b <- if (need_constrain_b) make_contrast(n_b) else NULL
	# symmetric unipartite ties h_b to h_a since a == b
	if (symmetric && !is.null(H_a)) H_b <- H_a
	list(H_a = H_a, H_b = H_b,
	     active = !is.null(H_a) || !is.null(H_b))
}
