# shared-grid quadrature for the expected pair problem of the dynamic binary
# empirical-Bayes estimator when the predictor noise is wide (src/dynamic_eb_pair.cpp).
# the 19 functions of the Fisher pair problem (.ae_pair_quad_fisher()) are tabulated once per correlation
# on the grid x_k = cc sinh(k du) (lde_pair_table_cpp()), and a pair's expectation
# contracts its two Gaussian weight vectors with the table (lde_pair_apply_cpp()).
# pairs the grid does not resolve go to the direct kernel (lde_pair_direct_cpp()), and
# a call without wide noise uses the shared Gauss-Hermite rule throughout. the inner
# fits of the scale search use the coarse tier (.lde_pair_tier1), reported fits the fine one

# a table depends on the correlation alone, so a family of tables at the
# Chebyshev-Lobatto nodes of a window [rc - h, rc + h] serves every correlation inside
# it by barycentric interpolation (lde_pair_table_blend_cpp()). a fit starts and ends
# with an empty cache, so its results do not depend on earlier fits

# the cache of the pair quadrature, in an environment of the namespace: the family of tables (`fam`), the grid of the last plan
# (`hint`) and the number of plans (`calls`)
.lde_pair_env <- new.env(parent = emptyenv())

# an empty cache (the start and the end of a fit)
.lde_pair_reset <- function() {
	.lde_pair_env$fam <- NULL
	.lde_pair_env$hint <- NULL
	.lde_pair_env$calls <- NULL
	invisible(NULL)
}

# tables in the correlation: options(lame.pair_family = FALSE) builds one exact table a correlation
.lde_pair_family_on <- function() !identical(getOption("lame.pair_family", TRUE), FALSE)

# options read by .lde_pair_quad() and the planner: lame.pair_grid = FALSE uses the
# shared Gauss-Hermite rule for every pair, "force" plans from the whole menu at every
# call, and lame.pair_du_scale multiplies the largest grid spacing

# the coarse tier for the inner fits of the scale search, whose objective must be right
# to 1e-3 deviance units: twice the spacing, table pairs from 0.7 of the spacing, and
# windows of 6.5 standard deviations
.lde_pair_tier1 <- list(du = 2, ratio = 0.7, cutoff = 6.5)

# the grids of the planner (the map constant cc and the spacing du in u), from coarse to fine
.lde_pair_menu_cc <- c(3, 3, 3, 2, 2, 1.5, 1.5)
.lde_pair_menu_du <- c(0.1, 0.07, 0.05, 0.05, 0.035, 0.035, 0.025)

# the largest spacing du that resolves the pair functions at the largest noise sd of the call
.lde_pair_du_cap <- function(smax) {
	if (smax <= 2.5) 0.07 else if (smax <= 8) 0.05 else if (smax <= 12) 0.035 else if (smax <= 20) 0.025 else 0.0125
}

# Gauss-Hermite nodes for one predictor in the direct kernel (one node without noise):
# the smallest rung that integrates the pair functions to about 1e-9
.lde_gh_rung <- function(sd) {
	rungs <- c(5L, 6L, 7L, 8L, 9L, 11L, 13L, 21L, 31L, 41L, 61L)
	k <- rungs[pmin(
		findInterval(sd, c(0.12, 0.16, 0.2, 0.25, 0.3, 0.4, 0.55, 0.8, 1.05, 1.3), left.open = TRUE) + 1L,
		length(rungs)
	)]
	k[!(sd > 1e-6)] <- 1L
	k
}

# a predictor wider than this, and resolved by the grid, is integrated on the grid nodes
# in the direct kernel; a narrower one by its Gauss-Hermite rule
.lde_axis_sd <- 1

# the plan of a call: NULL when no predictor has an sd of `sd_min` or more (the shared
# rule throughout); otherwise pairs whose two predictors the grid (cc, du) resolves use
# its table and the rest the direct kernel. the grid is the menu entry with the smallest
# predicted cost on a sample of at most 1500 pairs. `variant` is the observed outcome
# (0 .. 3), `ncol` the number of functions (6, 12 or 19), `tier` 2 for the fixed point
# and 1 for the scale search
.lde_pair_plan <- function(
	m1,
	m2,
	sd1,
	sd2,
	variant,
	rho,
	ncol,
	tier = 2L,
	force = FALSE,
	cutoff = 7.5,
	sd_min = 0.08,
	pad = 4L
) {
	n <- length(m1)
	if (tier <= 1L) cutoff <- .lde_pair_tier1$cutoff
	fin <- is.finite(m1) & is.finite(m2) & is.finite(sd1) & is.finite(sd2)
	if (!any(fin & (sd1 >= sd_min | sd2 >= sd_min))) return(NULL)
	ratio <- if (tier <= 1L) .lde_pair_tier1$ratio else 1
	smax <- max(pmax(sd1, sd2)[fin])
	cap <- .lde_pair_du_cap(smax) * (if (tier <= 1L) .lde_pair_tier1$du else 1) * getOption("lame.pair_du_scale", 1)
	cc_m <- .lde_pair_menu_cc[.lde_pair_menu_du <= cap + 1e-12]
	du_m <- .lde_pair_menu_du[.lde_pair_menu_du <= cap + 1e-12]
	if (!length(du_m) || !any(abs(du_m - cap) < 1e-9)) {
		cc_m <- c(cc_m, 3, 2)
		du_m <- c(du_m, cap, cap)
	}
	# the outcomes of the pairs of the call (the data's: the same in every call of a fit)
	vall <- 0L
	for (v in unique(variant)) vall <- bitwOr(vall, bitwShiftL(1L, v))
	# the family of tables of the cache serves the call when its window holds the correlation, its grid is fine enough for the
	# call and it has the outcomes
	cached <- .lde_pair_env$fam
	use_cached <- !is.null(cached) && abs(rho - cached$rc) <= cached$h + 1e-15 && cached$du <= cap + 1e-12 &&
		bitwAnd(vall, cached$vmask) == vall && cached$tier >= tier
	fam_on <- .lde_pair_family_on()
	nfam <- if (!fam_on) 1L else if (tier <= 1L) 4L else 6L
	# the grid of the last plan is reused for seven calls in eight (pairs change slowly between passes)
	hint <- .lde_pair_env$hint
	.lde_pair_env$calls <- (if (is.null(.lde_pair_env$calls)) 0L else .lde_pair_env$calls) + 1L
	menu_call <- .lde_pair_env$calls %% 8L == 0L
	if (!is.null(hint) && !menu_call && identical(hint$tier, tier) && abs(hint$cap - cap) < 1e-12 && !force) {
		cc_m <- hint$cc
		du_m <- hint$du
	}
	pad <- pad + (if (fam_on) 2L else 0L)
	idx <- if (n > 1500L) unique(round(seq(1, n, length.out = 1500L))) else seq_len(n)
	scale <- n / length(idx)
	# cost model in seconds: a table node, a contraction multiply-add, a weight exponential,
	# a direct-kernel node, a family blend; a new family is amortised over ten evaluations
	c_tab <- 1.5e-6; c_mac <- 0.5e-9; c_exp <- 25e-9; c_node <- 2e-6; c_blend <- 0.3e-9; life <- 10
	table_cost <- function(N, kd) {
		if (kd == "cached") {
			return(if (identical(cached$cur$rho, rho)) 0 else cached$N^2 * .lde_pair_tab_cols * cached$n * c_blend)
		}
		if (nfam > 1L) N^2 * (c_tab * nfam / life + .lde_pair_tab_cols * nfam * c_blend) else N^2 * c_tab
	}
	cost_of <- function(cc, du, kd, ix) {
		a1 <- m1[ix]; a2 <- m2[ix]; b1 <- sd1[ix]; b2 <- sd2[ix]
		lo1 <- a1 - cutoff * b1; hi1 <- a1 + cutoff * b1; lo2 <- a2 - cutoff * b2; hi2 <- a2 + cutoff * b2
		sp1 <- du * sqrt(cc^2 + (abs(a1) + 2 * b1)^2); sp2 <- du * sqrt(cc^2 + (abs(a2) + 2 * b2)^2)
		e1 <- b1 >= ratio * sp1; e2 <- b2 >= ratio * sp2
		tab <- if (kd == "direct") rep(FALSE, length(ix)) else fin[ix] & b1 >= sd_min & b2 >= sd_min & e1 & e2
		outside <- 0
		if (kd == "cached") {
			inside <- lo1 >= cached$xlo & hi1 <= cached$xhi & lo2 >= cached$xlo & hi2 <= cached$xhi
			outside <- sum(tab & !inside) / length(ix)
			tab <- tab & inside
		}
		cost <- 0
		if (any(tab)) {
			if (kd == "cached") {
				N <- cached$N
			} else {
				klo <- floor(asinh(min(c(lo1[tab], lo2[tab])) / cc) / du) - pad
				khi <- ceiling(asinh(max(c(hi1[tab], hi2[tab])) / cc) / du) + pad
				N <- khi - klo + 1
			}
			n1 <- (asinh(hi1[tab] / cc) - asinh(lo1[tab] / cc)) / du + 1
			n2 <- (asinh(hi2[tab] / cc) - asinh(lo2[tab] / cc)) / du + 1
			cost <- table_cost(N, kd) + sum(ncol * n1 * n2 * c_mac + (n1 + n2) * c_exp) * scale
		}
		dr <- !tab
		if (any(dr)) {
			ax1 <- e1[dr] & b1[dr] > .lde_axis_sd; ax2 <- e2[dr] & b2[dr] > .lde_axis_sd
			nn1 <- ifelse(ax1, (asinh(hi1[dr] / cc) - asinh(lo1[dr] / cc)) / du + 1, .lde_gh_rung(b1[dr]))
			nn2 <- ifelse(ax2, (asinh(hi2[dr] / cc) - asinh(lo2[dr] / cc)) / du + 1, .lde_gh_rung(b2[dr]))
			cost <- cost + sum(nn1 * nn2) * c_node * scale
		}
		attr(cost, "outside") <- outside
		cost
	}
	# a family that serves the call is the only candidate between menu searches, unless more
	# than one pair in fifty falls outside its grid
	only_cached <- use_cached && !menu_call && !force
	repeat {
		cc_e <- c(if (use_cached) cached$cc, if (!only_cached) cc_m, if (!only_cached) 3)
		du_e <- c(if (use_cached) cached$du, if (!only_cached) du_m, if (!only_cached) cap)
		kind <- c(if (use_cached) "cached", if (!only_cached) rep("fresh", length(cc_m)), if (!only_cached) "direct")
		cl <- lapply(seq_along(cc_e), function(g) cost_of(cc_e[g], du_e[g], kind[g], idx))
		costs <- vapply(cl, function(z) as.numeric(z), numeric(1))
		if (only_cached && attr(cl[[1L]], "outside") > 0.02) {
			only_cached <- FALSE
			next
		}
		break
	}
	# a family in use stays unless another grid is clearly cheaper (it is a table of the call already paid for)
	adj <- costs
	if (use_cached) adj[1L] <- adj[1L] * 0.7
	g <- which.min(adj)
	cc <- cc_e[g]; du <- du_e[g]; kd <- kind[g]
	# the classification of every pair at the chosen grid
	lo1 <- m1 - cutoff * sd1; hi1 <- m1 + cutoff * sd1; lo2 <- m2 - cutoff * sd2; hi2 <- m2 + cutoff * sd2
	sp1 <- du * sqrt(cc^2 + (abs(m1) + 2 * sd1)^2); sp2 <- du * sqrt(cc^2 + (abs(m2) + 2 * sd2)^2)
	e1 <- sd1 >= ratio * sp1; e2 <- sd2 >= ratio * sp2
	tab <- if (kd == "direct") rep(FALSE, n) else fin & sd1 >= sd_min & sd2 >= sd_min & e1 & e2
	if (kd == "cached") {
		tab <- tab & lo1 >= cached$xlo & hi1 <= cached$xhi & lo2 >= cached$xlo & hi2 <= cached$xhi
	}
	plan <- list(
		cc = cc,
		du = du,
		tab = tab,
		cached = identical(kd, "cached"),
		cost = costs[g],
		vmask = vall,
		cutoff = cutoff
	)
	if (identical(kd, "fresh")) .lde_pair_env$hint <- list(cc = cc, du = du, tier = tier, cap = cap)
	if (any(tab)) {
		if (plan$cached) {
			plan$klo <- cached$klo; plan$khi <- cached$khi
		} else {
			plan$klo <- floor(asinh(min(c(lo1[tab], lo2[tab])) / cc) / du) - pad
			plan$khi <- ceiling(asinh(max(c(hi1[tab], hi2[tab])) / cc) / du) + pad
		}
		plan$N <- plan$khi - plan$klo + 1
	} else {
		plan$cached <- FALSE
	}
	dr <- which(!tab)
	if (length(dr)) {
		ax1 <- e1[dr] & sd1[dr] > .lde_axis_sd & is.finite(m1[dr]); ax2 <- e2[dr] & sd2[dr] > .lde_axis_sd & is.finite(m2[dr])
		plan$direct <- dr
		plan$type1 <- as.integer(ax1); plan$type2 <- as.integer(ax2)
		plan$K1 <- ifelse(ax1, 1L, .lde_gh_rung(sd1[dr]))
		plan$K2 <- ifelse(ax2, 1L, .lde_gh_rung(sd2[dr]))
	}
	plan
}

# the columns of a table (the information functions and the functions of the four outcomes)
.lde_pair_tab_cols <- 40L

# the family of tables of a plan's grid at the Chebyshev-Lobatto nodes of
# [rho - h, rho + h]; one exact table when the family is off, for |rho| > 0.9, or beyond
# options(lame.pair_family_bytes)
.lde_pair_family <- function(plan, rho, vmask, threads, tier) {
	N <- plan$N
	# refuse a table of more than 2 GB (the positions reproduce the panel exactly)
	if (as.double(N)^2 * .lde_pair_tab_cols * 8 > 2e9) {
		cli::cli_abort(c(
			"The pair quadrature needs a table of {N} x {N} nodes: the noise of the predictors is too wide.",
			"i" = "The positions or factors reproduce the panel exactly (their scale is not identified); lower {.arg R} or fit with {.code dcor = FALSE}."
		))
	}
	n <- if (!.lde_pair_family_on()) 1L else if (tier <= 1L) 4L else 6L
	if (n > 1L && (abs(rho) > 0.9 || n * as.double(N)^2 * .lde_pair_tab_cols * 8 > getOption("lame.pair_family_bytes", 3e8))) n <- 1L
	xmax <- plan$cc * sinh(max(abs(plan$klo), abs(plan$khi)) * plan$du)
	h <- if (n > 1L) min(if (tier <= 1L) 0.012 else 0.015, 0.05 * (1 - abs(rho)), 2 / xmax) else 0
	if (n > 1L) {
		j <- 0:(n - 1L)
		nodes <- rho - h * cos(pi * j / (n - 1L))
		bw <- (-1)^j
		bw[c(1L, n)] <- bw[c(1L, n)] / 2
	} else {
		nodes <- rho
		bw <- 1
	}
	tabs <- lapply(nodes, function(r) {
		lde_pair_table_cpp(
			r,
			plan$cc,
			plan$du,
			as.integer(plan$klo),
			as.integer(plan$khi),
			as.integer(vmask),
			as.integer(threads)
		)$table
	})
	x <- plan$cc * sinh(c(plan$klo, plan$khi) * plan$du)
	list(
		rc = rho,
		h = h,
		n = n,
		nodes = nodes,
		bw = bw,
		tabs = tabs,
		cc = plan$cc,
		du = plan$du,
		klo = plan$klo,
		khi = plan$khi,
		N = N,
		vmask = vmask,
		xlo = x[1L],
		xhi = x[2L],
		tier = tier,
		cur = NULL
	)
}

# the table of the plan for the correlation `rho`, as the list that lde_pair_apply_cpp() reads: a member of the family where the
# correlation is a node, else the combination of the members (the barycentric formula of the interpolation in the correlation)
.lde_pair_blend <- function(rho, threads = 1L) {
	fam <- .lde_pair_env$fam
	cur <- fam$cur
	if (!is.null(cur) && identical(cur$rho, rho)) return(cur$tab)
	d <- rho - fam$nodes
	hit <- which(abs(d) < 1e-13)
	vec <- if (length(hit)) {
		fam$tabs[[hit[1L]]]
	} else {
		lam <- fam$bw / d
		lde_pair_table_blend_cpp(fam$tabs, lam / sum(lam), as.integer(threads))
	}
	tab <- list(
		table = vec,
		N = fam$N,
		cc = fam$cc,
		du = fam$du,
		klo = fam$klo,
		khi = fam$khi,
		rho = rho,
		vmask = fam$vmask
	)
	.lde_pair_env$fam$cur <- list(rho = rho, tab = tab)
	tab
}

# the table of the plan for the correlation `rho` (a family built for the window when the plan has no family to serve
# it, else the combination of its members at `rho`), as the list that lde_pair_apply_cpp() reads
.lde_pair_table <- function(plan, rho, vmask, threads = 1L, tier = 2L) {
	if (!isTRUE(plan$cached)) {
		.lde_pair_env$fam <- .lde_pair_family(plan, rho, vmask, threads, tier)
	}
	.lde_pair_blend(rho, threads)
}
