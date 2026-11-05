# fast point estimation for the AME model by iterative block coordinate descent
# (BCD), adapted from the social influence regression estimator of Minhas and
# Hoff (2025, political analysis) as implemented in sir::sir_alsfit()

# internal helpers

.ae_default <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# pairwise symmetrization that keeps a dyad observed in only one direction;
# only dyads missing in both directions stay NA
.ae_symmetrize <- function(M) {
	Mt <- t(M)
	out <- M
	both <- is.finite(M) & is.finite(Mt)
	onlyt <- !is.finite(M) & is.finite(Mt)
	out[both] <- (M[both] + Mt[both]) / 2
	out[onlyt] <- Mt[onlyt]
	out
}

# minimum-norm solution of a symmetric positive semidefinite system A x = b by
# eigendecomposition, dropping negligible eigenvalues; MASS::ginv() is the
# fallback when eigen() fails
.ae_safe_solve <- function(A, bvec) {
	A <- (A + t(A)) / 2
	eg <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
	if (is.null(eg)) {
		out <- tryCatch(
			as.numeric(MASS::ginv(A) %*% bvec),
			error = function(e) rep(0, length(bvec))
		)
		return(as.numeric(out))
	}
	lam <- eg$values
	thr <- max(lam, 0) * sqrt(.Machine$double.eps)
	keep <- is.finite(lam) & lam > thr
	if (!any(keep)) return(rep(0, length(bvec)))
	Qk <- eg$vectors[, keep, drop = FALSE]
	as.numeric(Qk %*% (crossprod(Qk, bvec) / lam[keep]))
}

# .ae_safe_solve with the eigendecomposition done once; returns a function of
# the right-hand side, for a matrix that is fixed across BCD sweeps
.ae_safe_solver <- function(A) {
	A <- (A + t(A)) / 2
	eg <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
	if (is.null(eg)) return(function(bvec) .ae_safe_solve(A, bvec))
	lam <- eg$values
	thr <- max(lam, 0) * sqrt(.Machine$double.eps)
	keep <- is.finite(lam) & lam > thr
	if (!any(keep)) return(function(bvec) rep(0, length(bvec)))
	Qk <- eg$vectors[, keep, drop = FALSE]
	lk <- lam[keep]
	function(bvec) as.numeric(Qk %*% (crossprod(Qk, bvec) / lk))
}

# families the BCD estimator supports
.ae_supported_families <- c("normal", "binary", "poisson")

# fewest resampling blocks for a block bootstrap; with fewer, the parametric
# scheme is used with a warning
.ae_min_boot_blocks <- 10L

# identify node-covariate effects by orthogonality: regress the additive effect
# e on [1, W] over the active actors (those with observed dyads) and return the
# slopes and the residual heterogeneity, orthogonal to 1 and W
.ae_decompose <- function(e, W, active = NULL) {
	if (is.null(W) || ncol(W) == 0L) {
		return(list(intercept = 0, beta = numeric(0), resid = e))
	}
	if (is.null(active)) active <- rep(1, length(e))
	D <- cbind(1, W)
	Dw <- D * active
	g <- .ae_safe_solve(crossprod(Dw, D), crossprod(Dw, e))
	bcoef <- g[-1]
	names(bcoef) <- colnames(W)
	list(intercept = g[1], beta = bcoef, resid = as.numeric(e - D %*% g))
}

# connected components of the observed-dyad graph (union-find), as labels for
# row and column nodes; a bipartite graph has n_row + n_col nodes
.ae_components <- function(obs_count, bip) {
	nr <- nrow(obs_count); ncc <- ncol(obs_count)
	if (bip) {
		n <- nr + ncc
		ij <- which(obs_count > 0, arr.ind = TRUE)
		edges <- if (nrow(ij) == 0L) matrix(0L, 0, 2) else
			cbind(ij[, 1], nr + ij[, 2])
	} else {
		n <- nr
		adj <- (obs_count > 0) | t(obs_count > 0)
		adj[lower.tri(adj, diag = TRUE)] <- FALSE
		ij <- which(adj, arr.ind = TRUE)
		edges <- if (nrow(ij) == 0L) matrix(0L, 0, 2) else ij
	}
	parent <- seq_len(n)
	find <- function(x) {
		root <- x
		while (parent[root] != root) root <- parent[root]
		while (parent[x] != root) { nx <- parent[x]; parent[x] <<- root; x <- nx }
		root
	}
	for (e in seq_len(nrow(edges))) {
		ri <- find(edges[e, 1]); rj <- find(edges[e, 2])
		if (ri != rj) parent[ri] <- rj
	}
	root <- vapply(seq_len(n), find, integer(1))
	lab <- match(root, unique(root))
	if (bip) {
		list(
			row = lab[seq_len(nr)],
			col = lab[nr + seq_len(ncc)],
			n = length(unique(lab))
		)
	} else {
		list(row = lab, col = lab, n = length(unique(lab)))
	}
}

# number of connected components
.ae_n_components <- function(obs_count, bip) .ae_components(obs_count, bip)$n

# per-component additive gauge for directed networks: on a disconnected graph
# each component's free shift a += d, b -= d is set to its precision-weighted
# minimum-norm value, so a, b, va, vb and node coefficients are reproducible
.ae_component_gauge <- function(a, b, comp, w_row, w_col) {
	if (comp$n <= 1L) return(list(a = a, b = b))
	for (k in seq_len(comp$n)) {
		ri <- which(comp$row == k); ci <- which(comp$col == k)
		den <- sum(w_row[ri]) + sum(w_col[ci])
		if (den <= 0) next
		d <- (sum(w_col[ci] * b[ci]) - sum(w_row[ri] * a[ri])) / den
		a[ri] <- a[ri] + d
		b[ci] <- b[ci] - d
	}
	list(a = a, b = b)
}

# Gaussian working response on the latent scale: exact for normal, the sampler
# start-value transforms (see get_start_vals.R) otherwise
.ae_working_response <- function(Yarr, family) {
	d <- dim(Yarr)
	if (family == "normal") {
		return(Yarr)
	}
	if (family == "poisson") {
		return(log(Yarr + 1))
	}
	# binary: rank-based normal scores within each slice
	Z <- array(NA_real_, dim = d)
	for (t in seq_len(d[3])) {
		Yt <- Yarr[, , t]
		zt <- matrix(zscores(c(Yt)), d[1], d[2])
		if (family == "binary") {
			has0 <- any(Yt == 0, na.rm = TRUE)
			has1 <- any(Yt == 1, na.rm = TRUE)
			if (has0 && has1) {
				z01 <- 0.5 * (max(zt[Yt == 0], na.rm = TRUE) +
						min(zt[Yt == 1], na.rm = TRUE))
				zt <- zt - z01
			}
		}
		Z[, , t] <- zt
	}
	Z
}

# name a 2-d Xdyad slice after the argument expression, as ame() does; only a
# plain variable name is used, otherwise the positional "dyad<k>" default
.ae_dyad_expr_name <- function(expr) {
	if (!is.symbol(expr)) return(NULL)
	nm <- as.character(expr)
	if (length(nm) != 1L || is.na(nm) || !nzchar(nm) || nchar(nm) >= 40L) {
		return(NULL)
	}
	nm
}

# wrap a 2-d xdyad as a 1-slice 3-d array carrying `nm` as the slice name
.ae_name_dyad_slice <- function(x, nm) {
	if (is.null(x) || is.null(nm) || !is.matrix(x)) return(x)
	dn <- dimnames(x)
	out <- array(x, dim = c(nrow(x), ncol(x), 1L))
	dimnames(out) <- list(
		if (!is.null(dn)) dn[[1L]] else NULL,
		if (!is.null(dn)) dn[[2L]] else NULL,
		nm
	)
	out
}

# coerce one dyadic covariate slice to an n_row x n_col x pd array
.ae_dyad_to_3d <- function(x, nr, nc) {
	if (is.null(x)) return(NULL)
	if (length(dim(x)) == 3L) return(x)
	if (is.matrix(x)) {
		dn <- dimnames(x)
		out <- array(x, dim = c(nrow(x), ncol(x), 1L))
		dimnames(out) <- list(
			if (!is.null(dn)) dn[[1L]] else NULL,
			if (!is.null(dn)) dn[[2L]] else NULL,
			NULL
		)
		return(out)
	}
	cli::cli_abort("Each {.arg Xdyad} element must be a matrix or a 3D array.")
}

.ae_check_names <- function(nm, arg, slice, margin) {
	if (is.null(nm)) {
		cli::cli_abort(c(
			"Changing actor composition requires named network slices.",
			"i" = "{.arg {arg}} slice {slice} is missing {margin} names."
		))
	}
	if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
		cli::cli_abort(c(
			"{.arg {arg}} slice {slice} has invalid {margin} names.",
			"i" = "Actor names must be unique, non-missing, and non-empty."
		))
	}
	nm
}

.ae_union_names <- function(nms) {
	unique(unlist(nms, use.names = FALSE))
}

.ae_covariate_names <- function(x, p, prefix) {
	nm <- colnames(as.matrix(x))
	if (is.null(nm)) paste0(prefix, seq_len(p)) else nm
}

.ae_xdyad_list_for_balance <- function(Xdyad, Tt) {
	if (is.null(Xdyad)) return(NULL)
	if (is.array(Xdyad) && length(dim(Xdyad)) == 4L) {
		cli::cli_abort(c(
			"Changing actor composition with dyadic covariates requires {.arg Xdyad} as a per-period list.",
			"i" = "A 4D array has one fixed actor layout, so it cannot be aligned to changing actor sets."
		))
	}
	if (is.array(Xdyad) && length(dim(Xdyad)) == 3L && !is.list(Xdyad)) {
		cli::cli_abort(c(
			"{.arg Xdyad} for changing-composition longitudinal ALS is ambiguous as a 3D array.",
			"i" = "Supply a list of {Tt} matrices/arrays, one per time slice."
		))
	}
	if (!is.list(Xdyad)) Xdyad <- list(Xdyad)
	if (length(Xdyad) != Tt) {
		cli::cli_abort("{.arg Xdyad} must have one entry per time slice ({Tt}).")
	}
	lapply(Xdyad, .ae_dyad_to_3d, nr = NA_integer_, nc = NA_integer_)
}

.ae_node_list_for_balance <- function(Xn, Tt, arg) {
	if (is.null(Xn)) return(NULL)
	if (!is.list(Xn)) Xn <- list(Xn)
	if (length(Xn) != Tt) {
		cli::cli_abort("{.arg {arg}} must have one entry per time slice ({Tt}).")
	}
	lapply(Xn, as.matrix)
}

.ae_match_cov_names <- function(got, expected, arg, slice, margin) {
	if (is.null(got)) return(expected)
	if (!setequal(got, expected)) {
		cli::cli_abort(c(
			"{.arg {arg}} slice {slice} has {margin} names that do not match {.arg Y}.",
			"i" = "Covariates must be named with the same actors as the corresponding network slice."
		))
	}
	got
}

.ae_balance_longitudinal_inputs <- function(Ylist, Xdyad, Xrow, Xcol, bip) {
	Tt <- length(Ylist)
	dims <- vapply(Ylist, function(y) paste(dim(y), collapse = "x"), character(1))
	same_dims <- all(dims == dims[1L])
	rnames <- lapply(Ylist, rownames)
	cnames <- lapply(Ylist, colnames)
	all_named <- all(vapply(seq_len(Tt), function(t) {
		!is.null(rnames[[t]]) && !is.null(cnames[[t]])
	}, logical(1)))
	any_named <- any(vapply(seq_len(Tt), function(t) {
		!is.null(rnames[[t]]) || !is.null(cnames[[t]])
	}, logical(1)))
	names_same <- all_named && all(vapply(seq_len(Tt), function(t) {
		identical(rnames[[t]], rnames[[1L]]) &&
			identical(cnames[[t]], cnames[[1L]])
	}, logical(1)))
	unipartite_margin_mismatch <- !bip && all_named &&
		any(vapply(seq_len(Tt), function(t) {
			!identical(rnames[[t]], cnames[[t]])
		}, logical(1)))
	row_sets_same <- all_named && all(vapply(seq_len(Tt), function(t) {
		setequal(rnames[[t]], rnames[[1L]])
	}, logical(1)))
	col_sets_same <- all_named && all(vapply(seq_len(Tt), function(t) {
		setequal(cnames[[t]], cnames[[1L]])
	}, logical(1)))
	actor_sets_change <- !same_dims || (all_named && !(row_sets_same && col_sets_same))
	needs_balance <- !same_dims || (all_named && !names_same) ||
		unipartite_margin_mismatch

	if (!needs_balance) {
		if (any_named && !all_named) {
			cli::cli_warn(c(
				"Some network slices have actor names and others do not.",
				"i" = "ALS will use positional matching because the dimensions are fixed; give every slice names to align by actor."
			))
		}
		nr <- nrow(Ylist[[1L]])
		nc <- ncol(Ylist[[1L]])
		time_names <- .ae_default(names(Ylist), paste0("t", seq_len(Tt)))
		row_presence <- matrix(
			TRUE,
			nr,
			Tt,
			dimnames = list(rownames(Ylist[[1L]]), time_names)
		)
		col_presence <- matrix(
			TRUE,
			nc,
			Tt,
			dimnames = list(colnames(Ylist[[1L]]), time_names)
		)
		return(list(
			Y = Ylist,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			changing_composition = FALSE,
			row_presence = row_presence,
			col_presence = col_presence,
			actor_presence = if (!bip && nr == nc) row_presence else NULL
		))
	}

	if (!all_named) {
		cli::cli_abort(c(
			"Changing actor composition requires row and column names on every network slice.",
			"i" = "Slice dimensions or actor names differ across time, so actors must be matched by name."
		))
	}
	for (t in seq_len(Tt)) {
		rnames[[t]] <- .ae_check_names(rnames[[t]], "Y", t, "row")
		cnames[[t]] <- .ae_check_names(cnames[[t]], "Y", t, "column")
		if (!bip && !setequal(rnames[[t]], cnames[[t]])) {
			cli::cli_abort(c(
				"Unipartite changing-composition panels require the same row and column actor set within each slice.",
				"i" = "Slice {t} has row and column names that do not match."
			))
		}
	}

	row_names <- .ae_union_names(rnames)
	col_names <- if (bip) .ae_union_names(cnames) else row_names
	time_names <- .ae_default(names(Ylist), paste0("t", seq_len(Tt)))
	nr <- length(row_names)
	nc <- length(col_names)
	row_presence <- matrix(FALSE, nr, Tt, dimnames = list(row_names, time_names))
	col_presence <- matrix(FALSE, nc, Tt, dimnames = list(col_names, time_names))

	Ybal <- vector("list", Tt)
	names(Ybal) <- names(Ylist)
	for (t in seq_len(Tt)) {
		yt <- as.matrix(Ylist[[t]])
		storage.mode(yt) <- "double"
		out <- matrix(NA_real_, nr, nc, dimnames = list(row_names, col_names))
		ri <- match(rnames[[t]], row_names)
		ci <- match(cnames[[t]], col_names)
		out[ri, ci] <- yt
		Ybal[[t]] <- out
		row_presence[ri, t] <- TRUE
		col_presence[ci, t] <- TRUE
	}

	Xdyad_bal <- NULL
	if (!is.null(Xdyad)) {
		Xd <- .ae_xdyad_list_for_balance(Xdyad, Tt)
		pd1 <- dim(Xd[[1L]])[3L]
		dn1 <- dimnames(Xd[[1L]])
		xn <- if (length(dn1) >= 3L && !is.null(dn1[[3L]])) {
			dn1[[3L]]
		} else paste0("dyad", seq_len(pd1))
		Xdyad_bal <- vector("list", Tt)
		names(Xdyad_bal) <- names(Ylist)
		for (t in seq_len(Tt)) {
			xt <- Xd[[t]]
			dd <- dim(xt)
			if (dd[1L] != nrow(Ylist[[t]]) || dd[2L] != ncol(Ylist[[t]])) {
				cli::cli_abort(c(
					"{.arg Xdyad} dimensions must match {.arg Y} before actor-set balancing.",
					"i" = "Y slice {t} is {nrow(Ylist[[t]])}x{ncol(Ylist[[t]])}; Xdyad slice {t} is {dd[1L]}x{dd[2L]}."
				))
			}
			if (dd[3L] != pd1) {
				cli::cli_abort(c(
					"All {.arg Xdyad} slices must hold the same number of covariates.",
					"i" = "Slice 1 has {pd1}; slice {t} has {dd[3L]}."
				))
			}
			dn <- dimnames(xt)
			xdn <- if (length(dn) >= 3L) dn[[3L]] else NULL
			if (!is.null(xdn) && !identical(xdn, xn)) {
				cli::cli_abort(c(
					"All {.arg Xdyad} slices must have the same covariate names in the same order.",
					"i" = "Slice {t} has different covariate names than slice 1."
				))
			}
			xr <- if (length(dn) >= 1L) dn[[1L]] else NULL
			xc <- if (length(dn) >= 2L) dn[[2L]] else NULL
			xr <- .ae_match_cov_names(xr, rnames[[t]], "Xdyad", t, "row")
			xc <- .ae_match_cov_names(xc, cnames[[t]], "Xdyad", t, "column")
			out <- array(
				0,
				dim = c(nr, nc, pd1),
				dimnames = list(row_names, col_names, xn)
			)
			out[match(xr, row_names), match(xc, col_names), ] <- xt
			Xdyad_bal[[t]] <- out
		}
	}

	pad_node <- function(Xn, arg, actor_names_t, actor_union, prefix) {
		if (is.null(Xn)) return(NULL)
		Xn <- .ae_node_list_for_balance(Xn, Tt, arg)
		p1 <- ncol(Xn[[1L]])
		xn <- .ae_covariate_names(Xn[[1L]], p1, prefix)
		out_list <- vector("list", Tt)
		names(out_list) <- names(Ylist)
		for (t in seq_len(Tt)) {
			xt <- as.matrix(Xn[[t]])
			if (nrow(xt) != length(actor_names_t[[t]]) || ncol(xt) != p1) {
				cli::cli_abort(c(
					"{.arg {arg}} dimensions must match the corresponding {.arg Y} actor margin before balancing.",
					"i" = "Slice {t} should be {length(actor_names_t[[t]])} x {p1}; got {nrow(xt)} x {ncol(xt)}."
				))
			}
			got_cols <- colnames(xt)
			if (!is.null(got_cols) && !identical(got_cols, xn)) {
				cli::cli_abort(c(
					"All {.arg {arg}} slices must have the same covariate columns.",
					"i" = "Slice {t} has different column names than slice 1."
				))
			}
			xr <- rownames(xt)
			xr <- .ae_match_cov_names(xr, actor_names_t[[t]], arg, t, "actor")
			out <- matrix(
				NA_real_,
				length(actor_union),
				p1,
				dimnames = list(actor_union, xn)
			)
			out[match(xr, actor_union), ] <- xt
			out_list[[t]] <- out
		}
		out_list
	}
	Xrow_bal <- pad_node(Xrow, "Xrow", rnames, row_names, "row")
	Xcol_bal <- pad_node(Xcol, "Xcol", cnames, col_names, "col")

	list(
		Y = Ybal,
		Xdyad = Xdyad_bal,
		Xrow = Xrow_bal,
		Xcol = Xcol_bal,
		changing_composition = isTRUE(actor_sets_change),
		row_presence = row_presence,
		col_presence = col_presence,
		actor_presence = if (!bip) row_presence else NULL
	)
}

# refuse a covariate that is missing everywhere, with the messages of ame()
# (cross-section) and lame() (panel)
.ae_check_all_na_cov <- function(Xdyad, Xrow, Xcol, longitudinal) {
	check <- function(Xl, nm, dyadic) {
		if (is.null(Xl)) return(invisible(NULL))
		p <- if (dyadic) dim(Xl[[1L]])[3L] else ncol(Xl[[1L]])
		all_na <- vapply(seq_len(p), function(k) {
			all(vapply(Xl, function(x) {
				if (dyadic) all(is.na(x[, , k])) else all(is.na(x[, k]))
			}, logical(1L)))
		}, logical(1L))
		bad <- which(all_na)
		if (!length(bad)) return(invisible(NULL))
		if (longitudinal) {
			j <- bad[1L]
			cli::cli_abort(c(
				"{.arg {nm}} covariate {.val {j}} is entirely NA in every period.",
				"i" = "An all-NA covariate cannot be identified; drop it or supply values."
			))
		}
		what <- if (length(bad) == 1L) "covariate" else "covariates"
		verb <- if (length(bad) == 1L) "is" else "are"
		cli::cli_abort(c(
			"{.arg {nm}} {what} {.val {bad}} {verb} entirely NA.",
			"i" = "An all-NA covariate cannot be identified; drop it or supply values."
		))
	}
	check(Xdyad, "Xdyad", TRUE)
	check(Xrow, "Xrow", FALSE)
	check(Xcol, "Xcol", FALSE)
	invisible(NULL)
}

# cells whose dyadic, sender or receiver covariate is missing in their period,
# [n_row, n_col, t]; as in ame() and lame() these cells are unobserved. a
# symmetric fit drops a dyad when either of its cells is missing
.ae_cov_na_cells <- function(Xdyad, Xrow, Xcol, nr, nc, Tt, symmetric = FALSE) {
	miss <- array(FALSE, dim = c(nr, nc, Tt))
	for (tp in seq_len(Tt)) {
		m <- matrix(FALSE, nr, nc)
		if (!is.null(Xdyad)) m <- m | (rowSums(is.na(Xdyad[[tp]]), dims = 2L) > 0)
		if (!is.null(Xrow)) m <- m | (rowSums(is.na(Xrow[[tp]])) > 0)
		if (!is.null(Xcol)) {
			m <- m | matrix(rowSums(is.na(Xcol[[tp]])) > 0, nr, nc, byrow = TRUE)
		}
		if (symmetric && nr == nc) m <- m | t(m)
		miss[, , tp] <- m
	}
	miss
}

# canonical covariate representation: dyadic covariates form the design array
# X (n_row x n_col x pd x t); node covariates stay as per-actor matrices W_row
# and W_col, since they are collinear with the additive effects and are
# estimated afterwards in .ae_assemble() by between-actor regression
.ae_build_design <- function(
	Xdyad,
	Xrow,
	Xcol,
	nr,
	nc,
	Tt,
	node_time_policy = c("mean", "dynamic", "quiet"),
	row_presence = NULL,
	col_presence = NULL
) {
	node_time_policy <- match.arg(node_time_policy)
	# dyadic covariates to the design array
	dcols <- list(); dnms <- character(0)
	if (!is.null(Xdyad)) {
		pd <- dim(Xdyad[[1]])[3]
		dn0 <- dimnames(Xdyad[[1]])
		dnm <- if (length(dn0) >= 3L) dn0[[3]] else NULL
		raw_dnm <- .ae_default(dnm, paste0("dyad", seq_len(pd)))
		for (t in seq_len(Tt)) {
			dnt <- dimnames(Xdyad[[t]])
			dnmt <- if (length(dnt) >= 3L) dnt[[3]] else NULL
			if (!is.null(dnmt) && !identical(dnmt, raw_dnm)) {
				cli::cli_abort(c(
					"All {.arg Xdyad} slices must have the same covariate names in the same order.",
					"i" = "Slice {t} has different covariate names than slice 1."
				))
			}
		}
		for (k in seq_len(pd)) {
			arr <- array(0, dim = c(nr, nc, Tt))
			for (t in seq_len(Tt)) arr[, , t] <- Xdyad[[t]][, , k]
			dcols[[length(dcols) + 1L]] <- arr
			# dyadic covariate names follow the MCMC suffix convention
			base_nm <- .ae_default(dnm[k], paste0("dyad", k))
			dnms <- c(dnms, .lame_apply_suffix(base_nm, "dyad"))
		}
	}
	pd <- length(dcols)
	X <- array(0, dim = c(nr, nc, pd, Tt))
	for (k in seq_len(pd)) X[, , k, ] <- dcols[[k]]
	# missing covariate cells are unobserved (.ae_cov_na_cells), so this 0 is never used
	X[is.na(X)] <- 0
	if (pd > 0) dimnames(X) <- list(NULL, NULL, dnms, NULL)

	# node covariates to per-actor means over observed periods, plus period arrays;
	# an actor missing in every period gets 0
	node_means <- function(Xn, n, suffix, lbl, presence = NULL) {
		if (is.null(Xn)) return(NULL)
		p <- ncol(as.matrix(Xn[[1]]))
		nm0 <- colnames(as.matrix(Xn[[1]]))
		W <- matrix(0, n, p)
		W_arr <- array(0, dim = c(n, p, length(Xn)))
		varies <- FALSE
		for (k in seq_len(p)) {
			vals <- matrix(
				vapply(Xn, function(m) as.matrix(m)[, k], numeric(n)),
				n,
				length(Xn)
			)
			# missing values of present actors (absent actors are padded, not missing)
			miss <- !is.finite(vals)
			if (!is.null(presence) && all(dim(presence) == dim(vals))) {
				miss <- miss & presence
			}
			W[, k] <- rowMeans(vals, na.rm = TRUE)
			W_arr[, k, ] <- vals
			if (any(miss)) {
				W[!is.finite(W[, k]), k] <- 0
				filled <- vals
				filled[miss] <- matrix(W[, k], n, ncol(vals))[miss]
				W_arr[, k, ] <- filled
			}
			# flag within-actor variation over time that the per-actor mean discards
			if (length(Xn) > 1L && !varies) {
				within <- mean(
					apply(vals, 1, stats::sd, na.rm = TRUE),
					na.rm = TRUE
				)
				if (is.finite(within) && within > 1e-8) varies <- TRUE
			}
		}
		if (varies && !identical(node_time_policy, "quiet")) {
			msg <- if (identical(node_time_policy, "dynamic")) {
				c(
					"A {lbl} covariate varies within actor across time slices.",
					"i" = "Dynamic ALS uses period-specific node values for node coefficients selected by {.arg dynamic_beta}.",
					"i" = "Static node coefficients still use the per-actor mean."
				)
			} else {
				c(
					"A {lbl} covariate varies within actor across time slices.",
					"i" = "ALS summarises each node covariate by its per-actor mean, so within-actor covariate changes are not used.",
					"i" = "Use {.fn lame} with {.arg method = \"als\"} and a node block in {.arg dynamic_beta} for period-specific node values on the dynamic ALS path."
				)
			}
			cli::cli_warn(msg)
		}
		W[!is.finite(W)] <- 0
		W_arr[!is.finite(W_arr)] <- 0
		# add row/col suffixes without doubling existing ones
		base_nm <- .ae_default(nm0, paste0(lbl, seq_len(p)))
		colnames(W) <- .lame_apply_suffix(base_nm, sub("^_", "", suffix))
		dimnames(W_arr) <- list(NULL, colnames(W), NULL)
		list(mean = W, arr = W_arr)
	}
	row_node <- node_means(Xrow, nr, "_row", "row", row_presence)
	col_node <- node_means(Xcol, nc, "_col", "col", col_presence)
	W_row <- if (is.null(row_node)) NULL else row_node$mean
	W_col <- if (is.null(col_node)) NULL else col_node$mean
	W_row_arr <- if (is.null(row_node)) NULL else row_node$arr
	W_col_arr <- if (is.null(col_node)) NULL else col_node$arr

	list(
		X = X,
		p_dyad = pd,
		x_names_dyad = dnms,
		W_row = W_row,
		W_col = W_col,
		W_row_arr = W_row_arr,
		W_col_arr = W_col_arr
	)
}

# normalize raw inputs to the canonical representation:
#   Yarr : n_row x n_col x t (unipartite diagonal set to NA)
#   X    : n_row x n_col x p x t
.ae_prepare <- function(
	Y,
	Xdyad,
	Xrow,
	Xcol,
	mode,
	longitudinal,
	node_time_policy = c("mean", "dynamic", "quiet"),
	symmetric = FALSE
) {
	node_time_policy <- match.arg(node_time_policy)

	# Y to a list of matrices
	if (longitudinal) {
		if (is.array(Y) && length(dim(Y)) == 3L) {
			Ylist <- lapply(seq_len(dim(Y)[3]), function(t) Y[, , t])
			names(Ylist) <- dimnames(Y)[[3]]
		} else if (is.list(Y)) {
			Ylist <- Y
		} else {
			cli::cli_abort(c(
				"{.arg Y} for {.fn lame_als} must be a list of matrices or a 3D array.",
				"i" = "Received an object of class {.cls {class(Y)[1]}}."
			))
		}
	} else {
		if (!is.matrix(Y)) {
			cli::cli_abort(c(
				"{.arg Y} for {.fn ame_als} must be a matrix.",
				"i" = "For longitudinal data use {.fn lame_als}."
			))
		}
		Ylist <- list(Y)
	}
	Tt <- length(Ylist)
	if (Tt < 1L) cli::cli_abort("{.arg Y} contains no network slices.")

	bip <- identical(mode, "bipartite")
	panel <- .ae_balance_longitudinal_inputs(Ylist, Xdyad, Xrow, Xcol, bip)
	Ylist <- panel$Y
	Xdyad <- panel$Xdyad
	Xrow <- panel$Xrow
	Xcol <- panel$Xcol

	nr <- nrow(Ylist[[1]]); nc <- ncol(Ylist[[1]])
	rn1 <- rownames(Ylist[[1]]); cn1 <- colnames(Ylist[[1]])
	for (t in seq_len(Tt)) {
		rnt_check <- rownames(Ylist[[t]])
		cnt_check <- colnames(Ylist[[t]])
		if (!is.null(rnt_check)) {
			.ae_check_names(rnt_check, "Y", t, "row")
		}
		if (!is.null(cnt_check)) {
			.ae_check_names(cnt_check, "Y", t, "column")
		}
		if (nrow(Ylist[[t]]) != nr || ncol(Ylist[[t]]) != nc) {
			cli::cli_abort(c(
				"All network slices must share the same dimensions.",
				"i" = "Slice 1 is {nr}x{nc} but slice {t} is {nrow(Ylist[[t]])}x{ncol(Ylist[[t]])}."
			))
		}
		if (longitudinal && t > 1L) {
			rnt <- rownames(Ylist[[t]]); cnt <- colnames(Ylist[[t]])
			if ((!is.null(rn1) && !is.null(rnt) && !identical(rnt, rn1)) ||
					(!is.null(cn1) && !is.null(cnt) && !identical(cnt, cn1))) {
				cli::cli_abort(c(
					"Network slices were not aligned to a common actor order.",
					"i" = "This is an internal input-preparation error; please report it."
				))
			}
		}
	}

	row_names <- .ae_default(rownames(Ylist[[1]]), paste0("row", seq_len(nr)))
	col_names <- .ae_default(colnames(Ylist[[1]]), paste0("col", seq_len(nc)))
	time_names <- .ae_default(names(Ylist), paste0("t", seq_len(Tt)))

	# a square unipartite network has one node set; rectangular inputs are left
	# for the later mode and dimension checks
	if (!bip && nr == nc) col_names <- row_names
	if (!is.null(panel$row_presence)) {
		dimnames(panel$row_presence) <- list(row_names, time_names)
	}
	if (!is.null(panel$col_presence)) {
		dimnames(panel$col_presence) <- list(col_names, time_names)
	}
	if (!bip && nr == nc && !is.null(panel$actor_presence)) {
		dimnames(panel$actor_presence) <- list(row_names, time_names)
	}

	# Y to the canonical array with the unipartite diagonal blanked
	Yarr <- array(NA_real_, dim = c(nr, nc, Tt))
	for (t in seq_len(Tt)) {
		yt <- as.matrix(Ylist[[t]])
		storage.mode(yt) <- "double"
		if (!bip) diag(yt) <- NA_real_
		Yarr[, , t] <- yt
	}
	dimnames(Yarr) <- list(row_names, col_names, time_names)

	# warn on a period with no observed off-diagonal entries
	if (Tt > 1L) {
		empty <- which(vapply(
			seq_len(Tt),
			function(t) !any(is.finite(Yarr[, , t])),
			logical(1)
		))
		if (length(empty) > 0L) {
			cli::cli_warn(c(
				"Time slice{?s} {.val {empty}} {?has/have} no observed entries.",
				"i" = "Empty slices contribute nothing to the fit; check for an all-missing network."
			))
		}
	}

	# covariates to per-period lists
	norm_dyad <- function(Xd) {
		if (is.null(Xd)) return(NULL)
		if (!longitudinal) {
			Xd <- list(Xd)
		} else if (is.array(Xd) && length(dim(Xd)) == 4L) {
			Xd <- lapply(seq_len(dim(Xd)[4]), function(t) Xd[, , , t])
		} else if (is.array(Xd) && length(dim(Xd)) == 3L && !is.list(Xd)) {
			# a 3-d array is ambiguous for a panel (covariates or time?)
			cli::cli_abort(c(
				"{.arg Xdyad} for {.fn lame_als} is a 3D array, which is ambiguous.",
				"i" = "Supply a list of {Tt} matrices/arrays (one per time slice),",
				"i" = "or a 4D array {.code [n_row, n_col, p, T]}."
			))
		} else if (!is.list(Xd)) {
			Xd <- list(Xd)
		}
		if (length(Xd) != Tt) {
			cli::cli_abort("{.arg Xdyad} must have one entry per time slice ({Tt}).")
		}
		lapply(Xd, .ae_dyad_to_3d, nr = nr, nc = nc)
	}
	norm_node <- function(Xn, arg) {
		if (is.null(Xn)) return(NULL)
		if (!longitudinal || !is.list(Xn)) Xn <- list(Xn)
		if (length(Xn) != Tt) {
			cli::cli_abort("{.arg {arg}} must have one entry per time slice ({Tt}).")
		}
		lapply(Xn, as.matrix)
	}

	Xdyad <- norm_dyad(Xdyad)
	Xrow <- norm_node(Xrow, "Xrow")
	Xcol <- norm_node(Xcol, "Xcol")

	y_has_row_names <- !is.null(rownames(Ylist[[1L]]))
	y_has_col_names <- !is.null(colnames(Ylist[[1L]]))
	align_dyad_slice <- function(xt, t) {
		dn <- dimnames(xt)
		xr <- if (length(dn) >= 1L) dn[[1L]] else NULL
		xc <- if (length(dn) >= 2L) dn[[2L]] else NULL
		if (is.null(xr) && is.null(xc)) return(xt)
		if (!y_has_row_names || !y_has_col_names) return(xt)
		if (is.null(xr) || is.null(xc)) {
			cli::cli_abort(c(
				"{.arg Xdyad} slice {t} has incomplete dyad names.",
				"i" = "Provide both row and column names, or provide neither and use the same order as {.arg Y}."
			))
		}
		.ae_match_cov_names(xr, row_names, "Xdyad", t, "row")
		.ae_match_cov_names(xc, col_names, "Xdyad", t, "column")
		xt[row_names, col_names, , drop = FALSE]
	}
	align_node_slice <- function(xn, expected, arg, t, margin, y_named) {
		got <- rownames(xn)
		if (is.null(got)) return(xn)
		if (!y_named) return(xn)
		.ae_match_cov_names(got, expected, arg, t, margin)
		xn[expected, , drop = FALSE]
	}
	if (!is.null(Xdyad)) {
		Xdyad <- lapply(seq_len(Tt), function(t) align_dyad_slice(Xdyad[[t]], t))
	}
	if (!is.null(Xrow)) {
		Xrow <- lapply(seq_len(Tt), function(t) {
			align_node_slice(
				Xrow[[t]],
				row_names,
				"Xrow",
				t,
				"row",
				y_has_row_names
			)
		})
	}
	if (!is.null(Xcol)) {
		Xcol <- lapply(seq_len(Tt), function(t) {
			align_node_slice(
				Xcol[[t]],
				col_names,
				"Xcol",
				t,
				"column",
				y_has_col_names
			)
		})
	}

	if (!is.null(Xdyad)) {
		pd1 <- dim(Xdyad[[1]])[3]
		for (t in seq_len(Tt)) {
			dd <- dim(Xdyad[[t]])
			if (dd[1] != nr || dd[2] != nc) {
				cli::cli_abort(c(
					"{.arg Xdyad} dimensions must match {.arg Y}.",
					"i" = "Y is {nr}x{nc}; Xdyad slice {t} is {dd[1]}x{dd[2]}."
				))
			}
			if (dd[3] != pd1) {
				cli::cli_abort(c(
					"All {.arg Xdyad} slices must hold the same number of covariates.",
					"i" = "Slice 1 has {pd1}; slice {t} has {dd[3]}."
				))
			}
		}
	}
	if (!is.null(Xrow)) {
		pr1 <- ncol(Xrow[[1]])
		for (t in seq_len(Tt)) {
			if (nrow(Xrow[[t]]) != nr || ncol(Xrow[[t]]) != pr1) {
				cli::cli_abort(c(
					"All {.arg Xrow} slices must be {nr} x {pr1} (rows = row/sender nodes).",
					"i" = "Slice {t} is {nrow(Xrow[[t]])} x {ncol(Xrow[[t]])}."
				))
			}
		}
	}
	if (!is.null(Xcol)) {
		pc1 <- ncol(Xcol[[1]])
		for (t in seq_len(Tt)) {
			if (nrow(Xcol[[t]]) != nc || ncol(Xcol[[t]]) != pc1) {
				cli::cli_abort(c(
					"All {.arg Xcol} slices must be {nc} x {pc1} (rows = column/receiver nodes).",
					"i" = "Slice {t} is {nrow(Xcol[[t]])} x {ncol(Xcol[[t]])}."
				))
			}
		}
	}

	# covariate missingness: refuse covariates missing everywhere, otherwise mark
	# cells with a missing covariate unobserved and record actors left with no cell
	.ae_check_all_na_cov(Xdyad, Xrow, Xcol, longitudinal)
	# a symmetric fit needs symmetric dyadic covariates: fill a hole in one triangle
	# from the other and refuse a slice that stays asymmetric
	if (isTRUE(symmetric) && !bip && !is.null(Xdyad)) {
		slice_nms <- dimnames(Xdyad[[1L]])
		xs <- array(
			unlist(Xdyad, use.names = FALSE),
			c(dim(Xdyad[[1L]]), length(Xdyad)),
			dimnames = list(NULL, NULL, if (length(slice_nms) >= 3L) slice_nms[[3L]], time_names)
		)
		if (any(vapply(Xdyad, anyNA, logical(1L)))) {
			xs <- .lame_fill_symmetric_dyadic(xs)
			for (t in seq_along(Xdyad)) Xdyad[[t]][] <- xs[, , , t]
		}
		.lame_check_symmetric_xdyad(
			if (length(Xdyad) > 1L) xs else
				array(xs, dim(xs)[1:3], dimnames = dimnames(xs)[1:3])
		)
	}
	cov_na <- list(n = 0L, row = integer(0), col = integer(0))
	if (any(vapply(c(Xdyad, Xrow, Xcol), anyNA, logical(1L)))) {
		miss <- .ae_cov_na_cells(
			Xdyad,
			Xrow,
			Xcol,
			nr,
			nc,
			Tt,
			symmetric = isTRUE(symmetric) && !bip
		)
		lost <- miss & !is.na(Yarr)
		if (any(lost)) {
			obs0 <- !is.na(Yarr)
			Yarr[lost] <- NA_real_
			obs1 <- !is.na(Yarr)
			cov_na$n <- sum(lost)
			cov_na$row <- unname(which(rowSums(obs0) > 0 & rowSums(obs1) == 0))
			cov_na$col <- unname(which(
				rowSums(colSums(obs0)) > 0 &
					rowSums(colSums(obs1)) == 0
			))
		}
	}

	des <- .ae_build_design(
		Xdyad,
		Xrow,
		Xcol,
		nr,
		nc,
		Tt,
		node_time_policy = node_time_policy,
		row_presence = panel$row_presence,
		col_presence = panel$col_presence
	)
	if (!is.null(des$X) && length(dim(des$X)) == 4L) {
		dimnames(des$X) <- list(
			row_names,
			col_names,
			dimnames(des$X)[[3L]],
			time_names
		)
	}
	if (!is.null(des$W_row)) rownames(des$W_row) <- row_names
	if (!is.null(des$W_col)) rownames(des$W_col) <- col_names
	if (!is.null(des$W_row_arr)) {
		dimnames(des$W_row_arr) <- list(
			row_names,
			dimnames(des$W_row_arr)[[2L]],
			time_names
		)
	}
	if (!is.null(des$W_col_arr)) {
		dimnames(des$W_col_arr) <- list(
			col_names,
			dimnames(des$W_col_arr)[[2L]],
			time_names
		)
	}

	list(
		Yarr = Yarr,
		X = des$X,
		p_dyad = des$p_dyad,
		x_names_dyad = des$x_names_dyad,
		W_row = des$W_row,
		W_col = des$W_col,
		W_row_arr = des$W_row_arr,
		W_col_arr = des$W_col_arr,
		cov_na = cov_na,
		nr = nr,
		nc = nc,
		Tt = Tt,
		bip = bip,
		row_names = row_names,
		col_names = col_names,
		time_names = time_names,
		changing_composition = panel$changing_composition,
		row_presence = panel$row_presence,
		col_presence = panel$col_presence,
		actor_presence = panel$actor_presence
	)
}

# core estimator: iterative block coordinate descent

# exact solve of the two-way weighted additive normal equations
#   [ diag(cnt_row)  wsum         ] [a]   [base_rs]
#   [ wsum'          diag(cnt_col)] [b] = [base_cs]
# (symmetric: (diag(cnt_row) + wsum) a = base_rs, b = a), the fallback when
# Gauss-Seidel crawls under very unbalanced weights. the system has a one-
# dimensional gauge; the eigen pseudo-inverse gives the minimum-norm solution
.ae_additive_solve <- function(
	base_rs,
	base_cs,
	wsum,
	cnt_row,
	cnt_col,
	symmetric
) {
	psolve <- function(K, rhs) {
		eg <- eigen((K + t(K)) / 2, symmetric = TRUE)
		thr <- max(abs(eg$values)) * sqrt(.Machine$double.eps)
		keep <- abs(eg$values) > thr
		Qk <- eg$vectors[, keep, drop = FALSE]
		as.vector(Qk %*% (crossprod(Qk, rhs) / eg$values[keep]))
	}
	if (symmetric) {
		K <- diag(x = cnt_row, nrow = length(cnt_row)) + wsum
		a <- psolve(K, base_rs)
		a[cnt_row == 0] <- 0
		return(list(a = a, b = a))
	}
	np <- length(cnt_row); nq <- length(cnt_col)
	K <- matrix(0, np + nq, np + nq)
	K[seq_len(np), seq_len(np)] <- diag(x = cnt_row, nrow = np)
	K[np + seq_len(nq), np + seq_len(nq)] <- diag(x = cnt_col, nrow = nq)
	K[seq_len(np), np + seq_len(nq)] <- wsum
	K[np + seq_len(nq), seq_len(np)] <- t(wsum)
	sol <- psolve(K, c(base_rs, base_cs))
	a <- sol[seq_len(np)]; b <- sol[np + seq_len(nq)]
	a[cnt_row == 0] <- 0; b[cnt_col == 0] <- 0
	list(a = a, b = b)
}

# Cholesky factor (one column effect pinned) for a direct start of the additive
# block, since the weights are fixed within a BCD fit. the solution is moved
# along the gauge (a_i + t, b_j - t) to keep the count-weighted column sum of the
# Gauss-Seidel recursion. returns NULL for a disconnected graph or a large,
# singular or ill-conditioned system, and the caller uses the plain recursion
.ae_additive_factor <- function(
	wsum,
	cnt_row,
	cnt_col,
	symmetric,
	max_dim = 2500L,
	min_rcond = 2e-4
) {
	nr <- length(cnt_row); nc <- length(cnt_col)
	rows <- which(cnt_row > 0)
	factor_ok <- function(K) {
		ch <- tryCatch(chol(K), error = function(e) NULL)
		if (is.null(ch)) return(NULL)
		rc <- tryCatch(rcond(ch, triangular = TRUE), error = function(e) 0)
		if (!is.finite(rc) || rc < min_rcond) return(NULL)
		ch
	}
	if (symmetric) {
		if (!length(rows) || length(rows) > max_dim) return(NULL)
		K <- wsum[rows, rows, drop = FALSE]
		diag(K) <- diag(K) + cnt_row[rows]
		ch <- factor_ok(K)
		if (is.null(ch)) return(NULL)
		return(list(ch = ch, rows = rows))
	}
	cols <- which(cnt_col > 0)
	if (!length(rows) || !length(cols) ||
			length(rows) + length(cols) - 1L > max_dim) return(NULL)
	pin <- cols[which.max(cnt_col[cols])]
	free <- cols[cols != pin]
	nrow_act <- length(rows); ncol_free <- length(free)
	K <- matrix(0, nrow_act + ncol_free, nrow_act + ncol_free)
	r_idx <- seq_len(nrow_act); c_idx <- nrow_act + seq_len(ncol_free)
	K[cbind(r_idx, r_idx)] <- cnt_row[rows]
	K[cbind(c_idx, c_idx)] <- cnt_col[free]
	if (ncol_free) {
		wf <- wsum[rows, free, drop = FALSE]
		K[r_idx, c_idx] <- wf
		K[c_idx, r_idx] <- t(wf)
	}
	ch <- factor_ok(K)
	if (is.null(ch)) return(NULL)
	list(ch = ch, rows = rows, free = free, cols = cols, cnt_act = cnt_col[cols])
}

# .ae_additive_solve with its eigendecomposition computed once per BCD fit
.ae_additive_exact_solver <- function(wsum, cnt_row, cnt_col, symmetric) {
	np <- length(cnt_row); nq <- length(cnt_col)
	if (symmetric) {
		K <- diag(x = cnt_row, nrow = np) + wsum
	} else {
		K <- matrix(0, np + nq, np + nq)
		K[seq_len(np), seq_len(np)] <- diag(x = cnt_row, nrow = np)
		K[np + seq_len(nq), np + seq_len(nq)] <- diag(x = cnt_col, nrow = nq)
		K[seq_len(np), np + seq_len(nq)] <- wsum
		K[np + seq_len(nq), seq_len(np)] <- t(wsum)
	}
	eg <- eigen((K + t(K)) / 2, symmetric = TRUE)
	thr <- max(abs(eg$values)) * sqrt(.Machine$double.eps)
	keep <- abs(eg$values) > thr
	Qk <- eg$vectors[, keep, drop = FALSE]
	lam <- eg$values[keep]
	function(base_rs, base_cs) {
		if (symmetric) {
			a <- as.vector(Qk %*% (crossprod(Qk, base_rs) / lam))
			a[cnt_row == 0] <- 0
			return(list(a = a, b = a))
		}
		sol <- as.vector(Qk %*% (crossprod(Qk, c(base_rs, base_cs)) / lam))
		a <- sol[seq_len(np)]; b <- sol[np + seq_len(nq)]
		a[cnt_row == 0] <- 0; b[cnt_col == 0] <- 0
		list(a = a, b = b)
	}
}

# Z        : n_row x n_col x t Gaussian working response (NA = unobserved)
# X        : n_row x n_col x p x t design array
# uv_prior : NULL (unpenalized factors), or list(Suv0, kappa0, sig2_u, sig2_v)
#            for the penalized MAP factor block with EM-updated scales (IRLS)
# ab_prior : NULL (unpenalized additive effects), or the EB actor-effect prior
#            (bias_correct = "eb"): regression and additive blocks are solved
#            jointly and Sab is updated by EM across calls
# returns  : list(mu, beta, a, b, u, v, l, omat, sig2_u, sig2_v, sse,
#                 iterations, converged, dev_history)
# the sweeps use only sums over time and cells, formed once by ae_bcd_stats_cpp
.ae_bcd_fit <- function(
	Z,
	X,
	R,
	symmetric,
	max_iter,
	tol,
	init = NULL,
	lowrank_method = "mm",
	weights = NULL,
	linear_solver = "eigen",
	uv_prior = NULL,
	certify = FALSE,
	ab_prior = NULL,
	unit_rules = FALSE
) {

	# validate tol and max_iter before they reach the low-level code
	if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1) {
		cli::cli_abort("{.arg max_iter} must be a single positive integer.")
	}
	if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
		cli::cli_abort("{.arg tol} must be a single positive number.")
	}

	d <- dim(Z)
	nr <- d[1]; nc <- d[2]; Tt <- d[3]
	p <- dim(X)[3]

	mask <- is.finite(Z)
	if (!any(mask)) {
		cli::cli_abort(c(
			"No observed cells: the working response is entirely missing.",
			"i" = "Check that {.arg Y} has observed (non-NA, off-diagonal) entries."
		))
	}
	Zf <- Z; Zf[!mask] <- 0

	# observation weights: unit by default, family working weights under IRLS
	Warr <- array(0, dim = d)
	if (is.null(weights)) {
		Warr[mask] <- 1
	} else {
		ww <- weights
		ww[!mask | !is.finite(ww)] <- 0
		ww[ww < 0] <- 0
		Warr <- ww
	}
	# weighted observed-cell total per dyad
	wsum <- rowSums(Warr, dims = 2L)
	cnt_row <- rowSums(wsum)
	cnt_col <- colSums(wsum)

	# initialization, warm started when supplied
	mu <- .ae_default(init$mu, 0)
	beta <- .ae_default(init$beta, rep(0, p))
	if (length(beta) != p) beta <- rep(0, p)
	a <- .ae_default(init$a, rep(0, nr))
	b <- .ae_default(init$b, rep(0, nc))
	if (length(a) != nr) a <- rep(0, nr)
	if (length(b) != nc) b <- rep(0, nc)
	U <- init$U; V <- init[["V"]]; L <- init$L
	Omat <- matrix(0, nr, nc)
	if (R > 0 && !is.null(U)) {
		if (symmetric && !is.null(L)) {
			Omat <- U %*% (L * t(U))
		} else if (!is.null(V)) {
			Omat <- tcrossprod(U, V)
		}
	}
	# EM-updated factor prior scales, carried across IRLS calls
	sig2_u <- uv_prior$sig2_u
	sig2_v <- uv_prior$sig2_v

	pa <- p + 1L
	# sums of the fixed data over time and cells:
	# sz  : sum_t w z                      (n_row x n_col)
	# gm  : sum_t w [1, x_1, ..., x_p]     (n_row * n_col x (p + 1)); first column wsum
	# xtz : X'W z over all cells, with X = [1, x_1, ..., x_p]
	# zz  : sum w z^2 ;  xtx : X'W X
	# gr, gc, szr, szc : row and column totals of gm and sz (all a rank-0 sweep needs)
	# unobserved cells have zero weight
	st <- ae_bcd_stats_cpp(Warr, Zf, X, nr, nc, Tt, p, .ae_threads())
	sz <- st$sz; gm <- st$gm; xtz <- st$xtz; zz <- st$zz; xtx <- st$xtx
	gr <- st$gr; gc <- st$gc; szr <- st$szr; szc <- st$szc

	# regression-block solver for [intercept | dyadic covariates], factored once
	# because the design and mask are fixed across sweeps
	lin <- .ae_linsolver(linear_solver, X, mask, Warr, nr, nc, Tt, p, XtX = xtx)
	# refuse an exactly rank-deficient pooled design, naming the aliased columns
	if (lin$rank_deficient) {
		cn <- c(
			"intercept",
			if (length(dimnames(X)) >= 3L && length(dimnames(X)[[3L]]) == p) dimnames(X)[[3L]] else if (p > 0L) paste0("x", seq_len(p))
		)
		# the shared refusal of ame() and lame() (R/design_rank.R), intercept first
		xlist <- lapply(seq_len(Tt), function(t) {
			xt <- array(1, c(nr, nc, p + 1L), dimnames = list(NULL, NULL, cn))
			if (p > 0L) xt[, , seq_len(p) + 1L] <- X[, , , t]
			xt
		})
		yobs <- Z; yobs[!mask] <- NA_real_
		.lame_check_design_rank(xlist, yobs)
		# fallback refusal when that check passes at its own tolerance
		qx <- qr(xtx / max(diag(xtx), 1e-300), tol = 1e-10)
		aliased <- cn[qx$pivot[-seq_len(qx$rank)]]
		cli::cli_abort(
			c(
				"The regression design (intercept + dyadic covariates) is exactly rank-deficient.",
				"x" = "Aliased column{?s}: {.val {aliased}}.",
				"i" = "Drop a constant covariate (the intercept already covers it) or one covariate of each linearly dependent set; near-collinearity is not refused."
			),
			class = "lame_design_rank"
		)
	}

	# EB actor-effect prior: normal equations of the stacked linear parameters,
	# the sums the Laplace correction needs, and the EM state of Sab
	eb_ctx <- NULL
	sab <- gamma <- eb_kd <- hpg <- NULL
	sab_change <- 0
	if (!is.null(ab_prior)) {
		eb_ctx <- .ae_eb_context(ab_prior, st, wsum, p)
		sab <- ab_prior$sab
		gamma <- ab_prior$gamma
		if (!is.null(ab_prior$hprime)) {
			hpg <- .ae_eb_hprime_sums(ab_prior$hprime, X, nr, nc, Tt, p)
		}
	}

	# weighted residual sum of squares from the cells, used for the final value and
	# when the sums-based form would lose digits
	sse_arrays <- function(mu, beta, cmat) {
		ae_bcd_sse_exact_cpp(Warr, Zf, X, nr, nc, Tt, p, mu, beta, cmat)
	}

	obs_dyad <- wsum > 0
	off_dyad <- which(!obs_dyad)
	# floor of the relative change in sse, so a perfect fit does not divide by zero;
	# on the response's own scale for a normal least-squares fit, 0.1 otherwise
	sse_floor <- if (unit_rules) 0.1 * zz / max(sum(wsum), 1) else 0.1
	sse_old <- Inf
	coef_mb_old <- rep(Inf, pa)
	dev_history <- numeric(0)
	converged <- FALSE
	add_converged <- TRUE # whether the final additive solve converged
	iter <- 0L
	# direct additive-block solver, factored once the fit has run long enough to repay it
	add_direct <- NULL
	add_direct_tried <- FALSE
	add_exact <- NULL # exact re-solve of the additive block, built if needed

	for (iter in seq_len(max_iter)) {

		if (!is.null(eb_ctx)) {
			# EB path (bias_correct = "eb"): intercept, coefficients and penalized actor
			# effects are solved jointly, then Sab is updated by EM; the penalty fixes the
			# level of the effects, so there is no centering
			stp <- .ae_eb_step(eb_ctx, ab_prior, if (R > 0L) Omat, st, sab, hpg)
			coef_mb <- stp$theta
			coef_delta <- max(abs(coef_mb - coef_mb_old))
			coef_mb_old <- coef_mb
			mu <- coef_mb[1]
			beta <- coef_mb[-1]
			a <- stp$a; b <- stp$b; gamma <- stp$gamma
			eb_kd <- stp$kd
			sab_change <- max(
				abs(log(diag(stp$sab_new) / diag(sab))),
				if (identical(ab_prior$kind, "directed"))
					abs(stp$sab_new[1, 2] - sab[1, 2]) / sqrt(prod(diag(sab))) else 0
			)
			sab <- stp$sab_new
			add_converged <- TRUE
			coef_fit <- coef_mb
			dmu <- 0
			if (R > 0L) rs <- list(phi = matrix(gm %*% coef_mb, nr, nc))
		} else {
			# intercept and regression coefficients by joint weighted least squares, so a
			# constant covariate shows as collinear with the intercept; the solver maps the
			# dyad-level sums to (mu, beta)
			coef_mb <- if (is.null(lin$solve_xty)) {
				lin$solve(Zf - c(outer(a, b, "+") + Omat)) # broadcast over t
			} else if (R == 0L) {
				lin$solve_xty(xtz - as.vector(crossprod(gr, a) + crossprod(gc, b)))
			} else {
				lin$solve_xty(xtz - ae_bcd_gtc_cpp(a, b, Omat, gm))
			}
			# change in the gauge-invariant (mu, beta), a second convergence signal
			coef_delta <- max(abs(coef_mb - coef_mb_old))
			coef_mb_old <- coef_mb
			mu <- coef_mb[1]
			beta <- coef_mb[-1]

			# additive sender and receiver effects: row and column totals of the weighted
			# residual without a and b, sz - gm [mu, beta] - wsum * o
			if (R == 0L) {
				base_rs <- szr - as.vector(gr %*% coef_mb)
				base_cs <- szc - as.vector(gc %*% coef_mb)
			} else {
				rs <- ae_bcd_rs_cpp(coef_mb, gm, sz, wsum, Omat)
				base_rs <- rs$base_rs
				base_cs <- rs$base_cs
			}
			# Gauss-Seidel sweeps over observed cells to the exact additive minimizer
			# (single-vector recursion when symmetric), started after a few sweeps from the
			# direct solution of the factored system
			if (!add_direct_tried && iter > 8L) {
				add_direct_tried <- TRUE
				add_direct <- .ae_additive_factor(wsum, cnt_row, cnt_col, symmetric)
			}
			add <- ae_additive_cpp(
				wsum,
				base_rs,
				base_cs,
				a,
				b,
				cnt_row,
				cnt_col,
				symmetric,
				add_direct
			)
			a <- add$a; b <- add$b
			add_scale <- add$scale
			add_resid <- add$resid
			# fall back to an exact direct solve when Gauss-Seidel crawls under very
			# unbalanced weights
			if (add_resid > 1e-6 * add_scale) {
				if (is.null(add_exact)) {
					add_exact <- .ae_additive_exact_solver(
						wsum,
						cnt_row,
						cnt_col,
						symmetric
					)
				}
				adir <- add_exact(base_rs, base_cs)
				a <- adir$a; b <- adir$b
				add_resid <- ae_additive_resid_cpp(
					wsum,
					base_rs,
					base_cs,
					a,
					b,
					cnt_row,
					cnt_col
				)
			}
			# stop when the normal-equation residual stays above tolerance after the exact solve
			if (!isTRUE(add_resid <= 1e-6 * add_scale)) {
				cli::cli_abort(
					c(
						"The additive sender and receiver effects could not be solved exactly.",
						"x" = "The residual of their normal equations after the exact solve was {signif(add_resid, 3)} against a tolerance of {signif(1e-6 * add_scale, 3)}.",
						"i" = "The observed-dyad weights make the system numerically singular beyond what its pseudo-inverse resolves; check for actors whose cells carry (nearly) no weight."
					),
					class = "lame_als_numerical"
				)
			}
			add_converged <- TRUE
			# sum-to-zero effects, with their means folded into mu
			mean_a <- mean(a); mean_b <- mean(b)
			mu <- mu + mean_a + mean_b
			a <- a - mean_a
			b <- b - mean_b
			coef_fit <- c(mu, beta)
			dmu <- mu - coef_mb[1] # the grand mean folded into the intercept
		}

		# multiplicative factors: minimize the weighted low-rank objective
		#   sum_ij c_ij (rbar_ij - o_ij)^2 ,  rank(o) <= r,
		# with c_ij = wsum[i,j] and rbar the weighted mean residual per dyad. without
		# `uv_prior` the solver is `lowrank_method` ("mm" or "als"/"hybrid"); with it
		# (IRLS) the block solves the penalized MAP problem, since the unpenalized
		# rank-r MLE diverges under quasi-separation (see ame_als_lowrank.R)
		if (R > 0) {
			Rbar <- ae_bcd_rbar_cpp(a, b, rs$phi, dmu, sz, wsum)
			# EB: remove the Laplace cell correction from the factor target
			if (!is.null(eb_kd)) Rbar[obs_dyad] <- Rbar[obs_dyad] - eb_kd[obs_dyad] / wsum[obs_dyad]
			if (is.null(uv_prior)) {
				# normal least squares measures low-rank tolerances in the response's units
				lr <- .ae_lowrank_update(
					Omat,
					Rbar,
					wsum,
					obs_dyad,
					R,
					symmetric,
					lowrank_method,
					rel = unit_rules
				)
			} else {
				lr <- .ae_lowrank_update_pen(
					Omat,
					Rbar,
					wsum,
					obs_dyad,
					R,
					symmetric,
					U,
					V,
					sig2_u,
					sig2_v,
					uv_prior$suv0,
					uv_prior$kappa0,
					off = off_dyad,
					second_moment = !is.null(ab_prior),
					cfloor = if (!is.null(ab_prior) && identical(ab_prior$family, "normal")) 0 else 1
				)
				sig2_u <- lr$sig2_u; sig2_v <- lr$sig2_v
			}
			U <- lr$U; V <- lr[["V"]]; L <- lr$L; Omat <- lr$Omat
		}

		# convergence on the relative change in weighted sse, expanded as
		#   zz - 2 theta'xtz + theta'X'WX theta - 2 c'(sz - gm theta) + c'(wsum * c)
		# with c = a_i + b_j + o_ij; the arrays are used when that loses digits
		cross <- if (R == 0L) {
			# c = a_i + b_j: both sums reduce to row and column totals
			c(
				sum(a * (szr - as.vector(gr %*% coef_fit))) +
					sum(b * (szc - as.vector(gc %*% coef_fit))),
				sum(cnt_row * a * a) + 2 * sum(a * as.vector(wsum %*% b)) +
					sum(cnt_col * b * b)
			)
		} else {
			ae_bcd_sse_cpp(a, b, Omat, rs$phi, dmu, sz, wsum)
		}
		sse <- zz - 2 * sum(coef_fit * xtz) +
			sum(coef_fit * (xtx %*% coef_fit)) - 2 * cross[1] + cross[2]
		if (!is.finite(sse) || sse < 1e-4 * zz) {
			sse <- sse_arrays(mu, beta, outer(a, b, "+") + Omat)
		}
		dev_history <- c(dev_history, sse)
		# from the second sweep, converge on an sse plateau with stable (mu, beta), so
		# a flat-sse saddle with drifting coefficients does not stop the fit
		coef_stable <- coef_delta <= tol * 100 * (max(abs(coef_mb)) + 1e-8)
		if (iter > 1L && is.finite(sse) && coef_stable &&
				abs(sse_old - sse) / (abs(sse_old) + sse_floor) < tol) {
			converged <- TRUE
			sse_old <- sse
			break
		}
		sse_old <- sse
	}

	# reported objective: the final sse from the arrays
	sse_old <- sse_arrays(mu, beta, outer(a, b, "+") + Omat)
	dev_history[length(dev_history)] <- sse_old

	out <- list(
		mu = mu,
		beta = beta,
		a = a,
		b = b,
		U = U,
		V = V,
		L = L,
		Omat = Omat,
		sig2_u = sig2_u,
		sig2_v = sig2_v,
		sse = sse_old,
		iterations = iter,
		converged = converged,
		add_converged = add_converged,
		dev_history = dev_history
	)
	if (!is.null(ab_prior)) {
		out$sab <- sab; out$gamma <- gamma; out$sab_change <- sab_change
	}
	# certificate for the unweighted, unpenalized objective only
	if (isTRUE(certify)) {
		if (!is.null(weights) || !is.null(uv_prior)) {
			cli::cli_abort("The Gaussian convergence certificate requires an unweighted, unpenalized fit.")
		}
		out <- .ae_gaussian_certify(out, Z, X, R, symmetric, max_iter, tol)
	}
	out
}

# mean-model degrees of freedom: regression block, gauged additive effects and
# the rank-r factor manifold; used for the df-corrected residual variance and
# by the bootstrap
.ae_df_mean <- function(p_dyad, n_row_active, n_col_active, R, symmetric) {
	df_reg <- p_dyad + 1L
	df_add <- if (symmetric) {
		max(n_row_active - 1L, 0L)
	} else {
		max(n_row_active - 1L, 0L) + max(n_col_active - 1L, 0L)
	}
	df_mult <- if (R > 0L) {
		if (symmetric) {
			R * n_row_active - R * (R + 1L) %/% 2L
		} else {
			R * (n_row_active + n_col_active - R)
		}
	} else {
		0L
	}
	df_reg + df_add + df_mult
}

# reported variance components: noise-corrected additive variances and the
# latent-scale error variance of a Poisson IRLS fit

# centered trace of a sampling covariance, trace(P S P) / (n - 1) with
# P = I - 11'/n: the excess that var() of centered estimates carries
.ae_centered_trace <- function(S) {
	n <- nrow(S)
	if (n < 2L) return(NA_real_)
	(sum(diag(S)) - sum(S) / n) / (n - 1L)
}

# inverse of an additive-block matrix that is singular in the gauge directions;
# `null_vec` spans the a_i + d / b_j - d shift, whose projector stabilizes the
# solve without changing the centered trace. a Cholesky factor with pivots
# spread beyond the eigenvalue cutoff (a disconnected graph) is refused and the
# minimum-norm pseudo-inverse used instead
.ae_gauge_inverse <- function(K, null_vec = NULL) {
	n <- nrow(K)
	scale <- mean(diag(K))
	if (!is.finite(scale) || scale <= 0) return(NULL)
	Kg <- K
	if (!is.null(null_vec)) {
		Kg <- Kg + (scale / n) * tcrossprod(null_vec)
	}
	ch <- tryCatch(chol(Kg), error = function(e) NULL)
	if (!is.null(ch)) {
		pivot <- diag(ch)^2
		if (all(is.finite(pivot)) && min(pivot) > 1e-10 * max(pivot)) {
			return(chol2inv(ch))
		}
	}
	# minimum-norm pseudo-inverse, the gauge of the reported a, b
	eg <- tryCatch(eigen(Kg, symmetric = TRUE), error = function(e) NULL)
	if (is.null(eg)) return(NULL)
	keep <- eg$values > 1e-10 * max(eg$values)
	if (!any(keep)) return(NULL)
	eg$vectors[, keep, drop = FALSE] %*%
		(t(eg$vectors[, keep, drop = FALSE]) / eg$values[keep])
}

# first-order sampling variance of the estimated additive effects, removed from
# var(a) and var(b) to report va and vb. the additive normal equations have
#   K = [[diag(rowsums(w)), w], [w', diag(colsums(w))]]   (a_i + b_j)
#   K = diag(rowsums(w)) + w                              (symmetric a_i + a_j)
# with w the pooled observed-cell weights, so the excess is the centered trace
# of ve * K^+, conditional on the other blocks. above `max_dim` the diagonal
# approximation ve * mean(1 / rowsums(w)) is used
.ae_effect_noise_var <- function(
	mask,
	Wobs,
	ve,
	symmetric,
	active_row,
	active_col,
	max_dim = 2000L
) {
	out <- list(va = NA_real_, vb = NA_real_)
	if (!is.finite(ve) || ve <= 0) return(out)
	d <- dim(mask)
	if (is.null(d) || length(d) != 3L) return(out)
	nr <- d[1L]; nc <- d[2L]; Tt <- d[3L]
	W <- matrix(0, nr, nc)
	for (t in seq_len(Tt)) {
		wt <- if (is.null(Wobs)) matrix(1, nr, nc) else matrix(Wobs[, , t], nr, nc)
		wt[!mask[, , t]] <- 0
		wt[!is.finite(wt)] <- 0
		W <- W + wt
	}
	ia <- which(active_row == 1)
	ib <- which(active_col == 1)
	diag_only <- function(w_tot, n) {
		w_tot <- w_tot[w_tot > 0]
		if (length(w_tot) < 2L) return(NA_real_)
		ve * mean(1 / w_tot)
	}
	if (symmetric) {
		if (length(ia) < 2L) return(out)
		Ws <- W[ia, ia, drop = FALSE]
		rs <- rowSums(Ws)
		if (length(ia) > max_dim) {
			out$va <- diag_only(rs, length(ia))
			return(out)
		}
		K <- diag(rs, nrow = length(ia)) + Ws
		S <- .ae_gauge_inverse(K)
		if (is.null(S)) return(out)
		out$va <- ve * .ae_centered_trace(S)
		return(out)
	}
	if (length(ia) < 2L && length(ib) < 2L) return(out)
	Wab <- W[ia, ib, drop = FALSE]
	na <- length(ia); nb <- length(ib)
	rs <- rowSums(Wab); cs <- colSums(Wab)
	if (na + nb > max_dim) {
		out$va <- if (na > 1L) diag_only(rs, na) else NA_real_
		out$vb <- if (nb > 1L) diag_only(cs, nb) else NA_real_
		return(out)
	}
	K <- rbind(
		cbind(diag(rs, nrow = na), Wab),
		cbind(t(Wab), diag(cs, nrow = nb))
	)
	S <- .ae_gauge_inverse(K, null_vec = c(rep(1, na), rep(-1, nb)))
	if (is.null(S)) return(out)
	if (na > 1L) out$va <- ve * .ae_centered_trace(S[seq_len(na), seq_len(na), drop = FALSE])
	if (nb > 1L) out$vb <- ve * .ae_centered_trace(S[na + seq_len(nb), na + seq_len(nb), drop = FALSE])
	out
}

# lag-k cross-actor covariance profile of a per-period path matrix (actors by
# periods), over actors present in both periods
.ae_path_lag_profile <- function(M, presence = NULL) {
	Tt <- ncol(M)
	if (is.null(presence) || !all(dim(presence) == dim(M))) {
		presence <- matrix(TRUE, nrow(M), Tt)
	}
	S <- matrix(NA_real_, Tt, Tt)
	for (s in seq_len(Tt)) for (t in s:Tt) {
		ok <- presence[, s] & presence[, t] & is.finite(M[, s]) & is.finite(M[, t])
		if (sum(ok) > 2L) {
			v <- stats::cov(M[ok, s], M[ok, t])
			S[s, t] <- v
			S[t, s] <- v
		}
	}
	g <- rep(NA_real_, Tt)
	cnt <- rep(0L, Tt)
	for (k in seq_len(Tt) - 1L) {
		idx <- cbind(seq_len(Tt - k), seq_len(Tt - k) + k)
		v <- S[idx]
		v <- v[is.finite(v)]
		if (length(v)) {
			g[k + 1L] <- mean(v)
			cnt[k + 1L] <- length(v)
		}
	}
	list(g = g, cnt = cnt)
}

# AR(1) persistence of within-actor variation from the decay of successive lag
# differences, g_k - g_(k+1) = tau2 rho^k (1 - rho), which do not involve the
# level; clipped away from the boundary, where the level is unstable
.ae_path_decay_rho <- function(profiles, lower = 0.05, upper = 0.85) {
	num <- 0; den <- 0
	for (pf in profiles) {
		g <- pf$g
		g <- g[-1L] # lag 0 carries the white noise
		g <- g[is.finite(g)]
		if (length(g) < 3L) next
		dd <- -diff(g)
		num <- num + sum(dd[-1L])
		den <- den + sum(dd[-length(dd)])
	}
	if (!is.finite(num) || !is.finite(den) || den <= 0) return(NA_real_)
	min(max(num / den, lower), upper)
}

# level variance of a fitted additive path: the intercept of a count-weighted
# regression of lag-k covariances (k >= 1) on rho^k under
#   cov(a_.s, a_.t) = va + tau2 rho^|s-t|   for s != t,
# which nets out estimation noise; with two periods the lag-1 covariance is returned
.ae_path_level_var <- function(profile, rho) {
	g <- profile$g; cnt <- profile$cnt
	Tt <- length(g)
	ks <- seq_len(Tt - 1L)
	ok <- is.finite(g[ks + 1L]) & cnt[ks + 1L] > 0L
	ks <- ks[ok]
	if (!length(ks)) return(list(va = NA_real_, tau2 = NA_real_, rho = rho))
	gk <- g[ks + 1L]; w <- cnt[ks + 1L]
	if (length(ks) < 2L || !is.finite(rho)) {
		return(list(va = max(sum(w * gk) / sum(w), 0), tau2 = NA_real_, rho = rho))
	}
	D <- cbind(1, rho^ks)
	WD <- D * w
	sol <- tryCatch(
		solve(crossprod(D, WD), crossprod(WD, gk)),
		error = function(e) NULL
	)
	if (is.null(sol) || anyNA(sol)) {
		return(list(va = max(sum(w * gk) / sum(w), 0), tau2 = NA_real_, rho = rho))
	}
	if (sol[2L] < 0) { # tau2 >= 0: a flat profile is all level
		return(list(va = max(sum(w * gk) / sum(w), 0), tau2 = 0, rho = rho))
	}
	list(va = max(sol[1L], 0), tau2 = sol[2L], rho = rho)
}

# latent-scale error variance of a Poisson IRLS fit: the lognormal-Poisson
# model has var[y] = mu + (exp(ve) - 1) mu^2, so ve is the root in
# phi = exp(ve) - 1 of the pearson identity
#   sum (y - mu)^2 / (mu + phi mu^2) = n_obs - df
.ae_pois_latent_var <- function(Yarr, EZ, mask, df_mean) {
	y <- as.vector(Yarr)
	eta <- unlist(EZ, use.names = FALSE)
	ok <- as.vector(mask) & is.finite(y) & is.finite(eta)
	y <- y[ok]
	eta <- eta[ok]
	n_obs <- length(y)
	if (n_obs < 3L) return(NA_real_)
	target <- log(max(n_obs - df_mean, 1))
	# work on the log scale so large means neither cap the rate nor overflow
	log_residual <- eta
	positive <- y > 0
	log_y <- log(y[positive])
	log_residual[positive] <- pmax(log_y, eta[positive]) +
		log(-expm1(-abs(log_y - eta[positive])))
	pearson <- function(variance) {
		log_phi <- if (variance == 0) -Inf else if (variance > 50) variance else log(expm1(variance))
		other <- log_phi + 2 * eta
		largest <- pmax(eta, other)
		log_denominator <- largest + log1p(exp(-abs(eta - other)))
		terms <- 2 * log_residual - log_denominator
		peak <- max(terms)
		if (peak == -Inf) return(-Inf)
		peak + log(sum(exp(terms - peak))) - target
	}
	if (pearson(0) <= 0) return(0)
	hi <- 1
	while (pearson(hi) > 0) {
		hi <- 2 * hi
		if (!is.finite(hi)) return(NA_real_)
	}
	root <- tryCatch(
		stats::uniroot(pearson, c(0, hi), tol = 1e-10)$root,
		error = function(e) NA_real_
	)
	if (is.finite(root)) max(root, 0) else NA_real_
}

# assemble the fitted-model object from a BCD solution
# Wobs: final IRLS weights (NULL = unit), for the working-scale dispersion
# estimate_rho: FALSE skips the reciprocity estimate (multistart summaries)
.ae_assemble <- function(
	fit,
	prep,
	family,
	mode,
	symmetric,
	R,
	longitudinal,
	call,
	Zwork = NULL,
	link = NULL,
	Wobs = NULL,
	estimate_rho = TRUE
) {

	nr <- prep$nr; nc <- prep$nc; Tt <- prep$Tt
	# an EB fit keeps its penalized optimum: no re-gauging, node coefficients from the joint solve
	eb <- !is.null(fit$sab)
	# `Zwork`, the final IRLS working response, overrides the one-shot transform
	Z <- if (is.null(Zwork)) .ae_working_response(prep$Yarr, family) else Zwork
	# symmetrize the working response as the symmetric fit does
	if (symmetric) for (t in seq_len(Tt)) {
		Z[, , t] <- .ae_symmetrize(Z[, , t])
	}
	mask <- is.finite(Z)

	mu <- fit$mu
	beta_dyad <- fit$beta # dyadic-covariate coefficients
	a_raw <- fit$a; b_raw <- fit$b # raw additive effects (absorb node structure)
	if (length(beta_dyad) > 0) names(beta_dyad) <- prep$x_names_dyad

	U <- fit$U; V <- fit[["V"]]; L <- fit$L
	Omat <- fit$Omat

	# gauge: double-center the product o and fold its row means into a, column
	# means into b and grand mean into mu, which leaves EZ unchanged and makes
	# a, b, u, v unique whichever BCD fixed point was reached
	if (R > 0 && !is.null(Omat) && !eb) {
		rmn <- rowMeans(Omat); cmn <- colMeans(Omat); gmn <- mean(Omat)
		Omat <- Omat - outer(rmn, rep(1, nc)) - outer(rep(1, nr), cmn) + gmn
		a_raw <- a_raw + rmn
		b_raw <- b_raw + cmn
		mu <- mu - gmn
		# re-impose sum-to-zero, folding the grand mean into mu
		mu <- mu + mean(a_raw) + mean(b_raw)
		a_raw <- a_raw - mean(a_raw)
		b_raw <- b_raw - mean(b_raw)
		# refactor the double-centered o, which keeps rank <= r (o = u v', or u l u' when symmetric)
		if (symmetric) {
			eg <- eigen((Omat + t(Omat)) / 2, symmetric = TRUE)
			idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(R)]
			U <- eg$vectors[, idx, drop = FALSE]
			L <- eg$values[idx]
			V <- U
		} else {
			sv <- svd(Omat)
			dd <- sqrt(pmax(sv$d[seq_len(R)], 0))
			U <- sweep(sv$u[, seq_len(R), drop = FALSE], 2, dd, "*")
			V <- sweep(sv$v[, seq_len(R), drop = FALSE], 2, dd, "*")
			L <- NULL
		}
		# keep the product equal to its factors after truncation
		O_R <- if (symmetric) U %*% (L * t(U)) else tcrossprod(U, V)
		if (max(abs(Omat - O_R)) > 1e-8 * max(1, max(abs(Omat)))) Omat <- O_R
	}
	if (!is.null(U)) dimnames(U) <- list(prep$row_names, paste0("dim", seq_len(R)))
	if (!is.null(V)) dimnames(V) <- list(prep$col_names, paste0("dim", seq_len(R)))

	# active actors (with observed dyads, symmetrized when undirected) and the
	# per-component additive gauge
	obs_count <- apply(!is.na(prep$Yarr), c(1, 2), sum)
	if (symmetric) obs_count <- obs_count + t(obs_count)
	cnt_row <- rowSums(obs_count); cnt_col <- colSums(obs_count)
	active_row <- as.numeric(cnt_row > 0)
	active_col <- as.numeric(cnt_col > 0)
	# pin each component's free shift of a directed fit (.ae_component_gauge)
	if (!symmetric && !eb) {
		comp <- .ae_components(obs_count, prep$bip)
		if (comp$n > 1L) {
			cg <- .ae_component_gauge(a_raw, b_raw, comp, cnt_row, cnt_col)
			a_raw <- cg$a; b_raw <- cg$b
		}
	}

	# linear predictor on the link scale from the raw additive effects
	EZ <- vector("list", Tt)
	ab_raw <- outer(a_raw, b_raw, "+")
	resid_list <- vector("list", Tt)
	ve_num <- 0; ve_den <- 0
	vew_num <- 0 # weighted (IRLS dispersion) numerator
	for (t in seq_len(Tt)) {
		xb <- matrix(0, nr, nc)
		if (prep$p_dyad > 0) for (k in seq_len(prep$p_dyad)) {
			xb <- xb + beta_dyad[k] * prep$X[, , k, t]
		}
		ez <- mu + xb + ab_raw + Omat
		EZ[[t]] <- ez
		r <- Z[, , t] - ez
		resid_list[[t]] <- r
		mt <- mask[, , t]
		ve_num <- ve_num + sum(r[mt]^2)
		ve_den <- ve_den + sum(mt)
		if (!is.null(Wobs)) {
			wt <- Wobs[, , t]
			ok <- mt & is.finite(wt)
			vew_num <- vew_num + sum(wt[ok] * r[ok]^2)
		}
	}
	# df-corrected residual variance sse / (n_obs - df_mean), the counterpart of the
	# MCMC posterior mean s2; symmetric fits count each dyad twice, so df is doubled
	df_mean <- .ae_df_mean(
		prep$p_dyad,
		sum(active_row),
		sum(active_col),
		R,
		symmetric
	)
	if (symmetric) df_mean <- 2L * df_mean
	ve <- if (ve_den > 0) ve_num / max(ve_den - df_mean, 1) else NA_real_
	# dispersion of the weighted working problem, the scale of the additive-block
	# sampling covariance; under IRLS the weighted mean square (quasi-likelihood dispersion)
	ve_noise <- if (!is.null(Wobs) && ve_den > 0) {
		vew_num / max(ve_den - df_mean, 1)
	} else {
		ve
	}
	# binary IRLS: the link fixes the latent error variance at 1; the working-scale
	# dispersion is kept as the diagnostic ve_working
	ve_working <- NULL
	if (identical(family, "binary") && !is.null(link)) {
		ve_working <- if (!is.null(Wobs) && ve_den > 0) {
			vew_num / ve_den
		} else if (ve_den > 0) {
			ve_num / ve_den
		} else {
			NA_real_
		}
		ve <- 1
	}
	# a normal EB fit estimates the noise variance as a hyperparameter (ame_als_gauss.R)
	if (eb && identical(family, "normal") && is.finite(fit$ve)) {
		ve_working <- ve
		ve <- fit$ve
		ve_noise <- ve
	}
	# Poisson IRLS fits the marginal log-mean, the latent predictor plus ve / 2: ve
	# comes from the count overdispersion and the level is moved to the latent scale
	# by subtracting ve / 2, so response-scale predictions are unchanged
	pois_latent_var <- 0
	if (identical(family, "poisson") && identical(link, "log")) {
		ve_working <- if (ve_den > 0) ve_num / max(ve_den - df_mean, 1) else NA_real_
		if (eb && is.finite(fit$ve)) {
			# an EB fit estimates the lognormal variance in the Poisson-lognormal model
			# (ame_als_pln.R), with its predictor already on the latent scale
			ve <- fit$ve
			pois_latent_var <- ve
		} else {
			ve_lat <- .ae_pois_latent_var(prep$Yarr, EZ, mask, df_mean)
			if (is.finite(ve_lat)) {
				ve <- ve_lat
				pois_latent_var <- ve_lat
				mu <- mu - ve_lat / 2
				for (t in seq_len(Tt)) EZ[[t]] <- EZ[[t]] - ve_lat / 2
			}
		}
	}

	# node-covariate coefficients by the orthogonality constraint: between-actor
	# regression slopes, with a and b the residual heterogeneity
	if (eb) {
		# EB: node coefficients are fixed effects of the joint solve; effects keep
		# their prior mean and are not re-centered
		ab_p <- fit$ab_prior
		pr <- if (is.null(ab_p$Wr)) 0L else ncol(ab_p$Wr)
		pc <- if (is.null(ab_p$Wc)) 0L else ncol(ab_p$Wc)
		g_r <- fit$gamma[seq_len(pr)]; g_c <- fit$gamma[pr + seq_len(pc)]
		if (pr) names(g_r) <- colnames(prep$W_row)
		if (pc) names(g_c) <- colnames(prep$W_col)
		if (symmetric) {
			# symmetric: one effect per actor, the covariate enters through both members of a dyad
			Mn <- .ae_eb_mean_design(ab_p)
			a <- a_raw - if (is.null(Mn)) 0 else as.vector(Mn %*% fit$gamma)
			# actors with no observed cell keep the prior mean
			a[!ab_p$act_a] <- 0
			b <- a
			mu <- mu - 2 * (sum(ab_p$wbar_r * g_r) + sum(ab_p$wbar_c * g_c))
		} else {
			a <- a_raw - if (pr) as.vector(ab_p$Wr %*% g_r) else 0
			b <- b_raw - if (pc) as.vector(ab_p$Wc %*% g_c) else 0
			a[!ab_p$act_a] <- 0
			b[!ab_p$act_b] <- 0
			mu <- mu - sum(ab_p$wbar_r * g_r) - sum(ab_p$wbar_c * g_c)
		}
		beta <- c(beta_dyad, g_r, g_c)
	} else if (symmetric) {
		W_node <- prep$W_row
		if (!is.null(prep$W_col)) {
			W_node <- if (is.null(W_node)) prep$W_col else cbind(W_node, prep$W_col)
		}
		dec <- .ae_decompose(a_raw, W_node, active_row)
		a <- dec$resid; b <- dec$resid
		mu <- mu + 2 * dec$intercept # the intercept enters g_i and g_j
		beta <- c(beta_dyad, dec$beta)
	} else {
		dec_r <- .ae_decompose(a_raw, prep$W_row, active_row)
		dec_c <- .ae_decompose(b_raw, prep$W_col, active_col)
		a <- dec_r$resid; b <- dec_c$resid
		mu <- mu + dec_r$intercept + dec_c$intercept
		beta <- c(beta_dyad, dec_r$beta, dec_c$beta)
	}
	names(a) <- prep$row_names
	names(b) <- prep$col_names

	# variance components, with VC columns as in ame() and lame(); va and vb are
	# variances of the covariate-orthogonal effects over non-isolated actors
	va_raw <- if (sum(active_row) > 1) stats::var(a[active_row == 1]) else NA_real_
	vb_raw <- if (symmetric) {
		NA_real_
	} else if (sum(active_col) > 1) {
		stats::var(b[active_col == 1])
	} else {
		NA_real_
	}
	# subtract the first-order sampling variance so va and vb estimate the effect
	# variances; the empirical values are kept as vc_raw
	vc_noise <- .ae_effect_noise_var(
		mask,
		Wobs,
		ve_noise,
		symmetric,
		active_row,
		active_col
	)
	va <- if (is.finite(va_raw) && is.finite(vc_noise$va)) {
		max(va_raw - vc_noise$va, 0)
	} else va_raw
	vb <- if (is.finite(vb_raw) && is.finite(vc_noise$vb)) {
		max(vb_raw - vc_noise$vb, 0)
	} else vb_raw
	# no sender-receiver covariance for symmetric or bipartite networks
	cab <- if (symmetric || prep$bip) {
		NA_real_
	} else {
		both <- active_row == 1 & active_col == 1
		if (sum(both) > 2) stats::cov(a[both], b[both]) else NA_real_
	}
	if (eb) {
		# EB: the actor-effect covariance is the variance component itself
		va <- fit$sab[1L, 1L]
		vb <- if (symmetric) NA_real_ else fit$sab[2L, 2L]
		cab <- if (symmetric || prep$bip) NA_real_ else fit$sab[1L, 2L]
	}
	rho <- 0
	if (!prep$bip && !symmetric) {
		# pool (i,j) and (j,i) residual pairs across periods
		e_ij <- numeric(0); e_ji <- numeric(0)
		for (t in seq_len(Tt)) {
			rt <- resid_list[[t]]
			ut <- upper.tri(rt)
			e_ij <- c(e_ij, rt[ut])
			e_ji <- c(e_ji, t(rt)[ut])
		}
		ok <- is.finite(e_ij) & is.finite(e_ji)
		if (sum(ok) > 2) {
			rho <- suppressWarnings(stats::cor(e_ij[ok], e_ji[ok]))
			if (!is.finite(rho)) rho <- 0
		}
	}
	VC <- c(va = va, cab = cab, vb = vb, rho = rho, ve = ve)

	# response-scale fitted values with the unipartite diagonal blanked; IRLS fits
	# use their link, the transform path its own inverse maps
	fitted_list <- vector("list", Tt)
	for (t in seq_len(Tt)) {
		fitted_list[[t]] <- if (!is.null(link)) {
			switch(
				paste(family, link),
				# EZ is on the latent scale, so the count mean adds ve / 2
				"poisson log" = exp(EZ[[t]] + pois_latent_var / 2),
				"binary logit" = stats::plogis(EZ[[t]]),
				"binary probit" = stats::pnorm(EZ[[t]]),
				EZ[[t]]
			)
		} else switch(
			family,
			normal = EZ[[t]],
			poisson = .poisson_inverse_log1p(EZ[[t]]),
			binary = stats::pnorm(EZ[[t]])
		)
		if (!prep$bip) { # self-ties are undefined
			diag(EZ[[t]]) <- NA_real_
			diag(fitted_list[[t]]) <- NA_real_
		}
		dimnames(fitted_list[[t]]) <- list(prep$row_names, prep$col_names)
		dimnames(EZ[[t]]) <- list(prep$row_names, prep$col_names)
	}
	names(EZ) <- prep$time_names
	names(fitted_list) <- prep$time_names
	# deviance of the reported fit
	dev_reported <- .ae_fit_deviance(
		prep$Yarr,
		fitted_list,
		resid_list,
		family,
		link,
		eta = EZ,
		latent_var = pois_latent_var
	)

	coefficients <- c(intercept = mu, beta)

	approx_note <- if (family == "normal" && eb) {
		paste0(
			"Empirical-Bayes estimate (normal family): the posterior means of the actor effects",
			if (R > 0) " and the latent factors (a Gaussian variational posterior)" else "",
			", the coefficients, the variance of the noise and the variance components by empirical Bayes",
			if (!is.null(fit$eb_mf$record$pair)) "; the two cells of a dyad are a bivariate normal pair with the dyadic correlation rho" else "",
			". For posterior summaries use ame()/lame() (MCMC)."
		)
	} else if (family == "normal") {
		if (R > 0) {
			paste0(
				"Least-squares (independent-error Gaussian ML) AME point ",
				"estimate; with R > 0 the low-rank term is non-convex, so this ",
				"is a local optimum."
			)
		} else {
			paste0(
				"Least-squares (independent-error Gaussian ML) AME point ",
				"estimate; the global solution for R = 0."
			)
		}
	} else if (eb && !is.null(fit$eb_mf)) {
		paste0(
			"Empirical-Bayes variational estimate (",
			family,
			", ",
			link,
			" link): ",
			"the posterior means of the actor effects and the latent factors under ",
			"a Gaussian variational posterior, the coefficients and the variance ",
			"components by empirical Bayes",
			if (family == "poisson") "; a count is Poisson with a lognormal noise of variance ve" else "",
			if (!is.null(fit$eb_mf$record$pair)) {
				paste0(
					"; the two cells of a dyad are a bivariate ",
					if (family == "poisson") "Poisson-lognormal" else "probit",
					" pair with the dyadic correlation rho"
				)
			} else "",
			". For posterior summaries use ",
			"ame()/lame() (MCMC)."
		)
	} else if (eb) {
		paste0(
			"Empirical-Bayes point estimate: penalised iteratively reweighted ",
			"least squares (",
			family,
			", ",
			link,
			" link) with a Gaussian prior ",
			"on the actor effects, its covariance estimated by EM, and ",
			"Laplace-marginal coefficients. For posterior summaries use ",
			"ame()/lame() (MCMC)."
		)
	} else if (!is.null(link)) {
		paste0(
			"Fixed-effects point estimate: iteratively reweighted least squares ",
			"(",
			family,
			", ",
			link,
			" link) fit of the GLM AME model. ",
			"For calibrated family-specific inference use ame()/lame() (MCMC)."
		)
	} else {
		paste0(
			"Working-response point estimate: iterative block coordinate descent ",
			"run on a Gaussian working response (latent ",
			family,
			" scores). For calibrated family-specific inference use ame()/lame() (MCMC)."
		)
	}

	out <- list(
		call = call,
		family = family,
		mode = mode,
		symmetric = symmetric,
		R = R,
		longitudinal = longitudinal,
		mu = mu,
		beta = beta,
		a = a,
		b = b,
		U = U,
		V = V,
		L = L,
		G = if (prep$bip && R > 0) diag(R) else NULL,
		# point estimate of the double-centered multiplicative surface (U L U' when
		# symmetric, U V' otherwise), named as in the MCMC fits
		UVPM = if (R > 0 && !is.null(fit$Omat)) {
			m <- Omat
			if (!is.null(prep$row_names) && nrow(m) == length(prep$row_names))
				rownames(m) <- prep$row_names
			if (!is.null(prep$col_names) && ncol(m) == length(prep$col_names))
				colnames(m) <- prep$col_names
			m
		} else NULL,
		ULUPM = if (R > 0 && isTRUE(symmetric) && !is.null(fit$Omat)) {
			m <- Omat
			if (!is.null(prep$row_names) && nrow(m) == length(prep$row_names))
				dimnames(m) <- list(prep$row_names, prep$row_names)
			m
		} else NULL,
		coefficients = coefficients,
		VC = VC,
		s2 = ve,
		ve_working = ve_working,
		# empirical variances before the noise correction, and the correction
		vc_raw = c(va = va_raw, vb = vb_raw),
		vc_noise = c(va = vc_noise$va, vb = vc_noise$vb),
		pois_latent_var = if (pois_latent_var > 0) pois_latent_var else NULL,
		# mean-model df of the corrected variances and the Poisson overdispersion test
		df_mean = df_mean,
		EZ = EZ,
		fitted = fitted_list,
		YPM = fitted_list,
		residuals = resid_list,
		iterations = fit$iterations,
		converged = fit$converged,
		deviance = dev_reported,
		objective = fit$sse,
		dev_history = fit$dev_history,
		n_time = Tt,
		dims = list(n_row = nr, n_col = nc, p = length(beta)),
		row_names = prep$row_names,
		col_names = prep$col_names,
		x_names = names(beta),
		time_names = prep$time_names,
		Y = prep$Yarr,
		X = prep$X,
		W_row = prep$W_row,
		W_col = prep$W_col,
		changing_composition = isTRUE(prep$changing_composition),
		row_presence = prep$row_presence,
		col_presence = prep$col_presence,
		actor_presence = prep$actor_presence,
		dynamic_uv = FALSE,
		dynamic_ab = FALSE,
		dynamic_beta = FALSE,
		dynamic_G = FALSE,
		dynamic = list(
			uv = FALSE,
			uv_kind = NULL,
			ab = FALSE,
			beta = FALSE,
			beta_kind = NULL,
			G = FALSE,
			snap_model_estimated = FALSE,
			t_model_estimated = FALSE,
			mode = mode
		),
		meta = list(
			sampler = "bcd",
			algorithm = "iterative block coordinate descent (Hoff & Minhas SIR estimator)",
			uncertainty_available = FALSE,
			approximation_note = approx_note
		)
	)
	if (eb) {
		ev <- fit$eb_vcov
		gam <- fit$gamma
		if (length(gam) == length(c(colnames(prep$W_row), colnames(prep$W_col)))) {
			names(gam) <- c(colnames(prep$W_row), colnames(prep$W_col))
		}
		out$eb <- list(
			sab = fit$sab,
			Sab0 = fit$ab_prior$Sab0,
			eta0 = fit$ab_prior$eta0,
			gamma = gam,
			sab_change = fit$sab_change,
			sig2_u = fit$sig2_u,
			sig2_v = fit$sig2_v,
			mu_int = fit$mu,
			a_raw = fit$a,
			b_raw = fit$b,
			log_marginal = if (!is.null(ev$evidence)) ev$evidence$log_marginal else NA_real_,
			vcov = ev$vcov,
			vcov_cond = ev$vcov_cond,
			se = if (!is.null(ev$vcov)) sqrt(diag(ev$vcov)) else NULL,
			se_method = if (is.null(fit$eb_mf)) {
				"laplace (penalised fisher information, factors included)"
			} else {
				"laplace (penalised expected information at the variational posterior, factor second moments included)"
			},
			vcov_error = if (is.null(ev$vcov)) attr(ev, "message") else NULL,
			estimator = if (is.null(fit$eb_mf)) "laplace" else "variational"
		)
		# columns retired for lack of support (R/ame_als_retire.R)
		if (!is.null(fit$retirement)) out$eb$retirement <- fit$retirement
		if (!is.null(fit$eb_mf)) {
			out$eb$variational <- fit$eb_mf$record
			out$eb$laplace <- ev$evidence
			if (R > 0L) {
				out$eb$factor_posterior <- fit$eb_mf$posterior
				out$eb$suv <- fit$eb_mf$posterior$prior$Suv
			}
			# dyadic correlation of a directed network
			pr <- fit$eb_mf$record$pair
			if (!is.null(pr)) {
				out$eb$pair <- c(pr, list(cross = fit$eb_mf$cross))
			}
		}
		out$meta$algorithm <- "iterative block coordinate descent with an empirical-bayes actor-effect prior"
		# count network: lognormal variance and the Laplace marginal likelihood profile it was chosen on
		if (!is.null(fit$ve)) out$eb$ve <- fit$ve
		if (!is.null(fit$ve_search)) out$eb$ve_search <- fit$ve_search
		if (!is.null(fit$ab_prior$ve_prior)) out$eb$ve_prior <- fit$ab_prior$ve_prior
	}
	# longitudinal fits are "lame_als" and both inherit from "ame_als"
	class(out) <- if (isTRUE(longitudinal)) c("lame_als", "ame_als") else "ame_als"
	if (isTRUE(estimate_rho)) .ae_update_response_rho(out, link) else out
}

# attach a bootstrap result to a fresh fit, for the `bootstrap = n` argument
.ae_attach_bootstrap <- function(
	fit,
	R_boot,
	type,
	block_length,
	seed,
	verbose
) {
	# refuse too few replicates before running
	.ae_check_bootstrap_n(R_boot)
	bt <- ame_als_bootstrap(
		fit,
		R = R_boot,
		type = type,
		block_length = block_length,
		seed = seed,
		verbose = verbose
	)
	fit$bootstrap <- bt
	# record that uncertainty is available for confint, summary and print
	if (is.null(fit$meta)) fit$meta <- list()
	fit$meta$uncertainty_available <- TRUE
	fit$meta$uncertainty_source <- "bootstrap"
	fit$meta$bootstrap_R <- R_boot
	# the scheme that ran, which differs from the request on a panel too short for blocks
	fit$meta$bootstrap_type <- bt$type
	fit
}

# validate the dyadic correlation switch of ame_als() and lame_als()
.ae_check_dcor <- function(dcor) {
	if (!is.logical(dcor) || length(dcor) != 1L || is.na(dcor)) {
		cli::cli_abort("{.arg dcor} must be {.val TRUE} or {.val FALSE}.")
	}
	dcor
}

# validate the prior list of the static estimators; only s20, the inverse-gamma
# prior scale of the noise variance, is used, and other keys are reported
.ae_check_prior <- function(prior, family) {
	if (is.null(prior)) return(NULL)
	if (!is.list(prior) || (length(prior) > 0L && is.null(names(prior)))) {
		cli::cli_abort("{.arg prior} must be NULL or a named list.")
	}
	if (!is.null(prior$s20)) {
		prior$s20 <- .lame_validate_s20(prior$s20)
		if (!family %in% c("poisson", "pois", "normal", "nrm")) {
			cli::cli_warn("{.code prior$s20} is read by the Poisson and normal estimators only; {.arg family} is {.val {family}}, so it is ignored.")
			prior$s20 <- NULL
		}
	}
	extra <- setdiff(names(prior), "s20")
	if (length(extra)) {
		cli::cli_warn(c(
			"{.arg prior}: {.val {extra}} {?has/have} no counterpart in the ALS estimators and {?is/are} ignored.",
			"i" = "{.fn ame_als} and {.fn lame_als} read {.code prior$s20} (Poisson and normal); use {.fn ame} / {.fn lame} (MCMC) for the other keys."
		))
	}
	prior
}

# shared driver for ame_als() and lame_als()
.ae_run <- function(
	Y,
	Xdyad,
	Xrow,
	Xcol,
	R,
	family,
	mode,
	symmetric,
	max_iter,
	tol,
	verbose,
	seed,
	longitudinal,
	call,
	lowrank_method = "mm",
	non_normal_method = "irls",
	link = "probit",
	linear_solver = "eigen",
	multistart = "none",
	node_time_policy = c("mean", "dynamic", "quiet"),
	bias_correct = "none",
	certify_gaussian = TRUE,
	report_cov_na = TRUE,
	eb_posterior = c("variational", "mode"),
	dcor = FALSE,
	dcor_given = FALSE,
	prior = NULL
) {
	node_time_policy <- match.arg(node_time_policy)
	eb_posterior <- match.arg(eb_posterior)
	# NULL takes the family default once the estimation path is known
	bias_correct <- .ae_check_bias_correct(bias_correct)

	# amen-era family aliases, as in ame() and lame()
	amen_aliases <- c(
		nrm = "normal", bin = "binary", ord = "ordinal",
		pois = "poisson"
	)
	if (length(family) == 1L && is.character(family) &&
			family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
			"i" = "Update your script to {.code family = \"{new_family}\"} when convenient."
		))
		family <- new_family
	}
	# exact family matching with the error message of ame() and lame()
	valid_families <- c("normal", "binary", "ordinal", "cbin", "frn", "poisson")
	if (length(family) != 1L || !is.character(family) ||
			!family %in% valid_families) {
		cli::cli_abort(c(
			"{.arg family} = {.val {family}} is not a recognised family.",
			"i" = "Choose one of: {.val {valid_families}}."
		))
	}
	if (!family %in% .ae_supported_families) {
		supported <- .ae_supported_families
		cli::cli_abort(c(
			"family = {.val {family}} is not supported by the ALS estimator.",
			"i" = "Supported families: {.val {supported}}.",
			"i" = "{.val {family}} is a rank/censoring likelihood with no well-defined Gaussian working response; use {.fn ame} or {.fn lame} (MCMC) instead."
		))
	}
	mode <- match.arg(mode, c("unipartite", "bipartite"))
	# BCD is deterministic; `seed` only drives the random starts of `multistart`

	# Y must be numeric with no Inf or NaN
	y_vals <- if (is.list(Y)) unlist(Y, use.names = FALSE) else Y
	if (!is.numeric(y_vals) && !is.logical(y_vals)) {
		cli::cli_abort("{.arg Y} must be numeric.")
	}
	if (any(is.infinite(y_vals) | is.nan(y_vals))) {
		cli::cli_abort(c(
			"{.arg Y} contains infinite or NaN values.",
			"i" = "Only finite values and {.val NA} (missing) are allowed."
		))
	}
	# Poisson responses must be non-negative integer counts
	if (identical(family, "poisson")) .lame_check_poisson_counts(Y)
	# refuse Inf or NaN in covariates
	for (cov_nm in c("Xdyad", "Xrow", "Xcol")) {
		cov_val <- get(cov_nm)
		if (!is.null(cov_val)) {
			vals <- if (is.list(cov_val)) {
				unlist(lapply(cov_val, function(z) as.numeric(z)), use.names = FALSE)
			} else as.numeric(cov_val)
			if (any(is.infinite(vals) | is.nan(vals))) {
				cli::cli_abort(c(
					"{.arg {cov_nm}} contains infinite or NaN values.",
					"i" = "Only finite values and {.val NA} (missing) are allowed."
				))
			}
		}
	}
	.lame_als_check_covariate_scale(Xdyad, Xrow, Xcol)

	prep <- .ae_prepare(
		Y,
		Xdyad,
		Xrow,
		Xcol,
		mode,
		longitudinal,
		node_time_policy = node_time_policy,
		symmetric = isTRUE(symmetric)
	)
	if (isTRUE(report_cov_na)) .lame_inform_cov_na(prep$cov_na$n)

	if (mode == "unipartite" && prep$nr != prep$nc) {
		cli::cli_abort(c(
			"Network is rectangular ({prep$nr}x{prep$nc}) but mode = {.val unipartite}.",
			"i" = "Use mode = {.val bipartite} for rectangular networks."
		))
	}
	if (symmetric && prep$bip) {
		cli::cli_warn("{.arg symmetric} is ignored for bipartite networks.")
		symmetric <- FALSE
	}
	if (symmetric && prep$nr != prep$nc) {
		cli::cli_abort("Symmetric models require a square network.")
	}
	# a symmetric fit needs a symmetric Y
	if (symmetric && prep$nr == prep$nc) {
		for (t in seq_len(prep$Tt)) {
			Yt <- prep$Yarr[, , t]
			off <- Yt - t(Yt)
			if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
				cli::cli_abort(c(
					"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
					"i" = "Symmetrize the input or set {.arg symmetric} = FALSE."
				))
			}
		}
	}

	if (length(R) != 1L || anyNA(R) || R < 0 ||
			!isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	R <- as.integer(R)
	rmax <- min(prep$nr, prep$nc)
	if (R >= rmax) {
		cli::cli_abort(c(
			"{.arg R} = {R} must be smaller than min(n_row, n_col) = {rmax}."
		))
	}

	# symmetrize dyadic covariates so beta' x_ij is symmetric
	if (symmetric && prep$p_dyad > 0) {
		for (k in seq_len(prep$p_dyad)) for (t in seq_len(prep$Tt)) {
			prep$X[, , k, t] <- .ae_symmetrize(prep$X[, , k, t])
		}
	}

	# EB actor-effect prior (ame_als_eb.R): the default for binary, Poisson and
	# normal fits; a normal fit uses the Gaussian cell of ame_als_gauss.R
	eb_eligible <- (family %in% c("binary", "poisson") && identical(non_normal_method, "irls")) ||
		identical(family, "normal")
	if (identical(bias_correct, "eb") && !eb_eligible) {
		why <- if (!family %in% c("binary", "poisson", "normal")) {
			"The empirical-Bayes actor-effect prior is defined for the binary, Poisson and normal likelihoods; {.arg family} is {.val {family}}."
		} else {
			"{.code non_normal_method = \"transform\"} fits a fixed transformed response and has no likelihood to penalise; use {.code non_normal_method = \"irls\"}."
		}
		cli::cli_abort(
			c(
				"{.code bias_correct = \"eb\"} is available for binary, Poisson and normal fits only.",
				"i" = why
			),
			class = "lame_als_eb"
		)
	}
	use_eb <- eb_eligible && (is.null(bias_correct) || identical(bias_correct, "eb"))

	# warn on a disconnected observed-dyad graph, where each component has its own additive level
	obs_count <- apply(!is.na(prep$Yarr), c(1, 2), sum)
	# actors left without a cell by missing covariates are not components of their own
	if (length(prep$cov_na$row) || length(prep$cov_na$col)) {
		gone_r <- prep$cov_na$row; gone_c <- prep$cov_na$col
		if (!prep$bip) gone_r <- gone_c <- intersect(gone_r, gone_c)
		obs_count <- obs_count[!seq_len(prep$nr) %in% gone_r,
			!seq_len(prep$nc) %in% gone_c, drop = FALSE]
	}
	ncomp <- .ae_n_components(obs_count, prep$bip)
	if (ncomp > 1L && !use_eb) {
		cli::cli_warn(c(
			"The observed-dyad graph has {ncomp} disconnected components.",
			"i" = "Additive effects are identified only up to a per-component level shift; compare fitted values / predictions rather than raw a, b across components."
		))
	}

	# family and data checks
	obs_y <- prep$Yarr[is.finite(prep$Yarr)]
	.ae_check_als_response(prep$Yarr, family)
	# an all-zero count network has no finite log-rate, since the intercept is unpenalized
	if (use_eb && identical(family, "poisson") && length(obs_y) && !any(obs_y > 0)) {
		cli::cli_abort(
			c(
				"Poisson ALS has no finite estimate: every observed count is zero.",
				"i" = "The log-rate of an all-zero network has no maximum (the likelihood rises as the intercept falls). Use {.fn ame} / {.fn lame} (MCMC), whose priors keep it finite."
			),
			class = "lame_als_separation"
		)
	}
	# a one-signed dyadic covariate nonzero only on zero counts has no finite coefficient
	if (use_eb && identical(family, "poisson") && prep$p_dyad > 0L && length(obs_y)) {
		xmat <- matrix(aperm(prep$X, c(1L, 2L, 4L, 3L)), length(prep$Yarr), prep$p_dyad)
		yv <- as.vector(prep$Yarr)
		for (k in seq_len(prep$p_dyad)) {
			on <- is.finite(yv) & is.finite(xmat[, k]) & xmat[, k] != 0
			if (any(on) && !any(yv[on] > 0) && (all(xmat[on, k] > 0) || all(xmat[on, k] < 0))) {
				nm <- prep$x_names_dyad[k]
				cli::cli_abort(
					c(
						"Poisson ALS has no finite estimate for the coefficient of {.field {nm}}: every observed count is zero where it is non-zero, and it has one sign.",
						"i" = "The likelihood rises as the coefficient runs to {if (all(xmat[on, k] > 0)) 'minus' else 'plus'} infinity. Use {.fn ame} / {.fn lame} (MCMC), whose coefficient prior keeps it finite, or drop the covariate."
					),
					class = "lame_als_separation"
				)
			}
		}
	}

	# a constant normal response has an unbounded noise-variance likelihood
	if (use_eb && identical(family, "normal") && length(obs_y) && !(diff(range(obs_y)) > 0)) {
		cli::cli_abort(
			c(
				"Normal ALS has no estimate: every observed response is {.val {obs_y[1L]}}.",
				"i" = "A constant response has no variation for the noise, the actor effects or the latent factors to explain. Check {.arg Y} (a constant network, or a slice of one)."
			),
			class = "lame_als_data"
		)
	}

	# refuse a normal response whose spread is outside [1e-6, 1e8], beyond the
	# conditioning of the compiled factor updates
	if (use_eb && identical(family, "normal") && length(obs_y) > 1L) {
		ysd <- stats::sd(obs_y)
		if (is.finite(ysd) && (ysd < 1e-6 || ysd > 1e8)) {
			cli::cli_abort(
				c(
					"The response has standard deviation {.val {signif(ysd, 3)}}, outside the range 1e-6 to 1e8 the normal ALS estimator is conditioned for.",
					"i" = "Rescale {.arg Y} (for example {.code Y / sd(Y, na.rm = TRUE)}); the coefficients and variances rescale with it."
				),
				class = "lame_als_scale"
			)
		}
	}

	# warn on a binary period with a single observed value
	if (family == "binary") {
		for (t in seq_len(prep$Tt)) {
			yt <- prep$Yarr[, , t]
			if (length(unique(yt[is.finite(yt)])) == 1L) {
				cli::cli_warn(c(
					"Time slice {t} has a single distinct observed {family} value.",
					"i" = "A degenerate slice carries no within-slice rank information; its working-response scores are constant."
				))
			}
		}
	}

	# binary IRLS: actors with a single outcome value on a side are left out on that
	# side and their effect is reported as not estimated
	db <- NULL
	if (family == "binary" && identical(non_normal_method, "irls") && !use_eb &&
			length(obs_y) && all(obs_y %in% c(0, 1))) {
		db <- .ae_drop_boundary(prep$Yarr, symmetric)
		if (!any(!is.na(db$Y))) {
			cli::cli_abort(
				c(
					"Binary ALS cannot estimate this network: every actor has a single outcome value on a side once the others are left out.",
					"i" = "A fixed-effects fit has no finite estimate of such an actor's effect, and leaving their cells out leaves none. Use {.fn ame} / {.fn lame} (MCMC)."
				),
				class = "lame_als_size"
			)
		}
		if (length(db$row) || length(db$col)) prep$Yarr <- db$Y else db <- NULL
	}
	# fixed-effects fits report effects of actors left without cells as not
	# estimated; report_cov_na is FALSE for the static pass of a dynamic fit
	if (!use_eb && isTRUE(report_cov_na)) {
		db <- .ae_add_unobserved(db, prep$cov_na, symmetric)
	}

	Z <- .ae_working_response(prep$Yarr, family)
	if (symmetric) for (t in seq_len(prep$Tt)) {
		Z[, , t] <- .ae_symmetrize(Z[, , t])
	}

	p_total <- prep$p_dyad +
		(if (is.null(prep$W_row)) 0L else ncol(prep$W_row)) +
		(if (is.null(prep$W_col)) 0L else ncol(prep$W_col))
	if (verbose) {
		cli::cli_h3("AME exploration (iterative block coordinate descent)")
		cli::cli_text("Network: {prep$nr} x {prep$nc} {ifelse(prep$bip, 'bipartite', 'unipartite')}, {prep$Tt} time slice{?s}, R = {R}")
		cli::cli_text("Family: {.field {family}}{ifelse(symmetric, ' (symmetric)', '')} | Covariates: {p_total}")
	}

	# als and hybrid low-rank solvers are for directed or bipartite networks only
	if (symmetric && lowrank_method != "mm" && R > 0) {
		cli::cli_warn(c(
			"{.arg lowrank_method} = {.val {lowrank_method}} is not available for symmetric models.",
			"i" = "The symmetric eigenmodel uses the majorise-minimise solver."
		))
		lowrank_method <- "mm"
	}

	# IRLS is the default for binary and Poisson
	use_irls <- (identical(non_normal_method, "irls") &&
			family %in% c("binary", "poisson")) || (use_eb && identical(family, "normal"))
	eff_link <- if (family == "poisson") "log" else if (use_eb && identical(family, "normal")) "identity" else link
	# a penalized factor block has its own inner solvers (see ame_als_lowrank.R);
	# record the one used
	if (use_irls && R > 0L) {
		lowrank_method <- if (symmetric) "mm" else "als"
	}

	# an EB fit is refined from its posterior mode to the variational posterior
	# (ame_als_meanfield.R); the certificate is the mode's, the reported convergence
	# the refinement's
	vb_follows <- use_irls && use_eb && identical(eb_posterior, "variational")
	# dyadic correlation of a directed binary probit fit: the variational stage
	# treats the two cells of a dyad as one observation (ame_als_pair.R)
	pair_cfg <- NULL
	if (isTRUE(dcor) && family == "binary" && !symmetric && identical(mode, "unipartite")) {
		if (identical(eff_link, "logit")) {
			if (dcor_given) {
				cli::cli_abort(c(
					"{.code dcor = TRUE} needs the probit link.",
					"i" = "The dyadic correlation is the correlation of the latent normal errors of the two cells of a dyad, as in {.fn ame} and {.fn lame}; use {.code link = \"probit\"} or {.code dcor = FALSE}."
				))
			}
		} else if (vb_follows) {
			pair_cfg <- list()
		}
	}
	# directed count network: dyadic correlation at the fitted state (ame_als_pln_pair.R)
	if (isTRUE(dcor) && family == "poisson" && !symmetric && identical(mode, "unipartite") && vb_follows) {
		pair_cfg <- list(count = TRUE)
	}
	# directed normal network: bivariate normal dyads, correlation estimated with
	# the other hyperparameters (ame_als_gauss_pair.R)
	if (isTRUE(dcor) && family == "normal" && !symmetric && identical(mode, "unipartite") && vb_follows) {
		.ae_gauss_check_reciprocity(prep, R)
		pair_cfg <- list(gauss = TRUE)
	}
	# requested rank; unsupported columns may be retired (R/ame_als_retire.R)
	R_requested <- R
	retirement <- NULL
	if (use_irls) {
		# the EB prior keeps every actor, with finite effects for actors without ties
		ab_prior <- if (use_eb) .ae_eb_prior(prep, symmetric, family, s20 = prior$s20) else NULL
		# count network: noise variance maximizes the Laplace marginal likelihood
		# (.ae_pln_ve_search), found by continuation
		pln_ve <- use_eb && identical(family, "poisson")
		# the iteration at rank Rr: posterior mode and, by default, its variational
		# refinement; rank-0 iterations are cached per start and shared
		rank0 <- list()
		irls_once_rank <- function(init, Rr) {
			if (Rr == 0L) {
				for (m in rank0) if (identical(m$init, init)) return(m$ir)
			}
			if (pln_ve) {
				run <- function(ve, from) {
					.ae_irls_run(
						prep,
						family,
						eff_link,
						Rr,
						symmetric,
						max_iter,
						tol,
						lowrank_method = lowrank_method,
						init = if (is.null(from)) init else .ae_pln_mode_init(from$fit, ve),
						linear_solver = linear_solver,
						warm = !is.null(from),
						ab_prior = ab_prior,
						ve_fixed = ve,
						warm_retry = FALSE
					)
				}
				sr <- .ae_pln_ve_search(run, .ae_pln_ve_lo, ve_prior = ab_prior$ve_prior)
				ir <- sr$ir
				ir$fit$ve_search <- sr[c("ve", "objective", "curve", "evaluations")]
			} else {
				# normal: the first mode starts the variational stage, and the mode at its
				# settled noise variance is the one certified (.ae_gauss_evidence)
				gauss_start <- vb_follows && identical(family, "normal")
				ir <- .ae_irls_run(
					prep,
					family,
					eff_link,
					Rr,
					symmetric,
					max_iter,
					tol,
					lowrank_method = lowrank_method,
					init = init,
					linear_solver = linear_solver,
					ab_prior = ab_prior,
					irls_maxit = if (gauss_start) min(max_iter, .ae_gauss_start_iters),
					vcov_mode = if (gauss_start) "none" else "full"
				)
				if (gauss_start) {
					ir$stationarity$start_pass <- ir$stationarity$pass
					ir$stationarity$pass <- TRUE
					ir$fit$converged <- TRUE
				}
			}
			# skip the refinement for an unconverged rank-1 mode, whose column is retired anyway
			if (vb_follows && !(Rr == 1L && !isTRUE(ir$stationarity$pass))) {
				ir <- .ae_mf_stage(
					prep,
					ir,
					ab_prior,
					Rr,
					symmetric,
					eff_link,
					max_iter,
					tol,
					pair = pair_cfg,
					ve_fixed = pln_ve,
					defer_uncertainty = TRUE
				)
				# normal: evidence of the mode at the settled noise variance
				if (identical(family, "normal")) {
					ir <- .ae_gauss_evidence(
						prep,
						ir,
						ab_prior,
						Rr,
						symmetric,
						eff_link,
						max_iter,
						tol,
						lowrank_method,
						linear_solver
					)
				}
			}
			if (Rr == 0L) rank0[[length(rank0) + 1L]] <<- list(init = init, ir = ir)
			ir
		}
		# add the factor-scale uncertainty to the reported iteration only
		irls_finish <- function(ir) {
			if (!vb_follows || !isTRUE(ir$stationarity$pass) || !isTRUE(ir$fit$eb_mf$record$converged)) {
				ir$mf_pending <- NULL
				return(ir)
			}
			.ae_mf_finish(ir)
		}
		# retire the weakest column and refit until the iteration settles (R/ame_als_retire.R)
		irls_once <- function(init) .ae_retire_columns(irls_once_rank, init, R, symmetric, vb_follows, max_iter)
		if (multistart != "none" && R > 0) {
			# multistart on the IRLS path: keep the converged start with the lowest
			# objective at its own scales (.ae_irls_objective), or the highest Laplace
			# marginal likelihood for EB; raw deviance would favour the least-shrunk optimum
			n_starts <- if (multistart == "full") 8L else 4L
			starts <- .ae_multistart_starts(Z, R, symmetric, n_starts, seed)
			irs <- lapply(starts, irls_once)
			sel <- .ae_irls_select(irs, ab_prior)
			obj <- sel$value
			converged <- sel$converged
			best <- sel$best
			ir <- irls_finish(irs[[best]])
			fit <- ir$fit
			fit$multistart_sse <- obj
			fit$multistart_deviance <- vapply(irs, function(z) z$data_deviance, numeric(1))
			fit$multistart_distance <- sel$distance
			summ <- .ae_multistart_summary(function(j) {
				z <- irs[[j]]
				.ae_assemble(
					z$fit,
					prep,
					family,
					mode,
					symmetric,
					.ae_fit_rank(z$fit),
					longitudinal,
					call,
					Zwork = z$Z,
					link = eff_link,
					Wobs = z$W,
					estimate_rho = FALSE
				)
			}, n_starts, best)
			fit$multistart_coefficients <- summ$coefficients
			fit$multistart_product_distance <- summ$product_distance
			fit$multistart_selected <- best
			fit$multistart_converged <- converged
			fit$multistart_criterion <- sel$criterion
			if (verbose) {
				cli::cli_alert_info("Multi-start IRLS: {n_starts} starts; kept start {best}.")
			}
			# a converged start is a different optimum when its fitted predictor differs
			# by more than one deviance unit from the kept start's
			if (sel$material) {
				crit <- sel$criterion
				n_far <- sum(sel$different)
				n_cmp <- sum(sel$converged) - 1L
				gap <- signif(max(sel$excess[sel$different]), 2)
				cli::cli_warn(c(
					"Multi-start found materially different local optima.",
					"i" = "{n_far} of the {n_cmp} other converged start{?s} ended at a fitted predictor more than one deviance unit from the kept start's; their {crit} is higher by up to {gap}.",
					"i" = "The converged start with the lowest {crit} is kept; inspect {.code fit$multistart_sse}, {.code fit$multistart_deviance} and {.code fit$multistart_distance}."
				))
			}
		} else {
			ir <- irls_finish(irls_once(NULL))
			fit <- ir$fit
		}
		if (verbose) {
			cli::cli_alert_info("IRLS ({family}, {eff_link} link): {ir$irls_iters} reweighting iteration{?s}.")
		}
		# rank of the fit, after any retirement
		retirement <- ir$retirement
		if (!is.null(retirement)) {
			R <- retirement$effective_rank
			fit$retirement <- retirement
		}
	} else if (multistart != "none" && R > 0) {
		n_starts <- if (multistart == "full") 8L else 4L
		fit <- .ae_multistart_fit(
			Z,
			prep$X,
			R,
			symmetric,
			max_iter,
			tol,
			lowrank_method,
			linear_solver,
			n_starts,
			seed,
			certify = certify_gaussian,
			unit_rules = identical(family, "normal"),
			assemble = function(f) {
				.ae_assemble(
					f,
					prep,
					family,
					mode,
					symmetric,
					R,
					longitudinal,
					call,
					estimate_rho = FALSE
				)
			}
		)
		if (verbose) {
			cli::cli_alert_info("Multi-start: {n_starts} starts; kept start {fit$multistart_selected}.")
		}
		ss <- fit$multistart_sse
		spread <- (max(ss) - min(ss)) / max(1, abs(min(ss)))
		if (spread > 1e-4) {
			cli::cli_warn(c(
				"Multi-start found materially different local optima (relative SSE spread {signif(spread, 2)}).",
				"i" = "The lowest-SSE converged fit is preferred; inspect the recorded start diagnostics."
			))
		}
	} else {
		# suggest multistart when R > 0 (nonconvex objective)
		if (R > 0L && multistart == "none" && verbose) {
			cli::cli_inform(c(
				"i" = "Single-start ALS at {.code R = {R}}. The latent objective is nonconvex.",
				"i" = "{.code multistart = \"cheap\"} (4 starts) or {.code \"full\"} (8 starts) keeps the converged start with the lowest criterion."
			))
		}
		fit <- .ae_bcd_fit(
			Z,
			prep$X,
			R,
			symmetric,
			max_iter,
			tol,
			lowrank_method = lowrank_method,
			linear_solver = linear_solver,
			certify = certify_gaussian,
			unit_rules = identical(family, "normal")
		)
	}

	# convergence of the posterior mode, or of the whole fit when no variational stage follows
	mode_ok <- if (vb_follows) isTRUE(ir$stationarity$pass) else isTRUE(fit$converged)
	if (verbose) {
		if (mode_ok) {
			n_it <- if (use_eb) ir$irls_iters else fit$iterations
			cli::cli_alert_success("Converged in {n_it} iteration{?s} (objective = {round(fit$sse, 4)}).")
		} else {
			cli::cli_alert_warning("Did not converge in {max_iter} iterations.")
		}
	}
	if (!mode_ok) {
		if (vb_follows) {
			# the rank evidence is the Laplace marginal likelihood at the mode, so an
			# unsettled mode is refused
			stz <- ir$stationarity
			# name the requested rank and the last one tried when columns were retired
			rank_note <- if (is.null(retirement)) sprintf("rank %d", R) else {
				k <- R_requested - R
				sprintf(
					"rank %d requested; %d column%s retired, the last attempt was at rank %d",
					R_requested,
					k,
					if (k == 1L) "" else "s",
					R
				)
			}
			cli::cli_abort(
				c(
					"The posterior mode did not reach its likelihood and variance convergence criteria within {.arg max_iter} = {max_iter} reweighting steps ({rank_note}).",
					"x" = "It stopped at {.val {if (is.null(stz$termination)) \"the iteration limit\" else stz$termination}} after {ir$irls_iters} steps; the Laplace marginal likelihood at the mode, which the rank evidence reads, is not defined.",
					"i" = "Increase {.arg max_iter}, fit a smaller rank, or use {.fn ame} / {.fn lame} (MCMC) at this rank."
				),
				class = "lame_als_convergence"
			)
		} else if (use_irls && use_eb) {
			cli::cli_warn(
				c(
					"IRLS did not reach its likelihood and variance convergence criteria.",
					"i" = "Inspect {.code fit$irls_diagnostics} and increase {.arg max_iter}."
				),
				class = "lame_als_convergence"
			)
		} else if (use_irls) {
			# refuse rather than return the last iterate
			.ae_abort_unconverged("irls", ir$stationarity, max_iter)
		} else if (isTRUE(certify_gaussian)) {
			.ae_abort_unconverged(
				"bcd",
				c(
					fit$gaussian_diagnostics,
					list(iterations = fit$iterations)
				),
				max_iter
			)
		} else {
			# the static pass that seeds a dynamic fit has its own certificate
			cli::cli_warn(
				c(
					"Block coordinate descent did not converge within {.arg max_iter} = {max_iter} iterations.",
					"i" = "Inspect {.code fit$gaussian_diagnostics} and increase {.arg max_iter}."
				),
				class = "lame_als_convergence"
			)
		}
	}

	if (vb_follows) {
		rec <- fit$eb_mf$record
		if (verbose) {
			cli::cli_alert_info("Variational posterior: {rec$iterations} iteration{?s}, {rec$termination}; lower bound {round(rec$elbo, 3)}.")
		}
		if (!isTRUE(rec$converged)) .ae_mf_refuse(rec, R, R_requested)
	}

	out <- if (use_irls) {
		.ae_assemble(
			fit,
			prep,
			family,
			mode,
			symmetric,
			R,
			longitudinal,
			call,
			Zwork = ir$Z,
			link = eff_link,
			Wobs = ir$W
		)
	} else {
		.ae_assemble(fit, prep, family, mode, symmetric, R, longitudinal, call)
	}
	# final objective: for IRLS the penalized deviance at the stopping iterate
	if (use_irls) {
		out$objective <- ir$best_obj
		out$irls_iterations <- ir$irls_iters
		out$irls_diagnostics <- ir$stationarity
		# EB: iterations count the reweighting steps
		if (use_eb) out$iterations <- ir$irls_iters
	}
	out$gaussian_diagnostics <- fit$gaussian_diagnostics
	out$lowrank_method <- lowrank_method
	out$meta$lowrank_method <- lowrank_method
	out$non_normal_method <- if (use_irls) "irls" else "transform"
	out$link <- if (use_irls) eff_link else NULL
	out$linear_solver <- linear_solver
	out$multistart <- multistart
	# requested rank, when fewer columns were supported
	if (!is.null(retirement)) out$R_requested <- R_requested
	# final IRLS weights (NULL = unit), used by the weighted sandwich
	out$obs_weights <- if (use_irls) ir$W else NULL
	if (!is.null(fit$multistart_sse)) {
		out$multistart_sse <- fit$multistart_sse
		if (!is.null(fit$multistart_deviance)) {
			out$multistart_deviance <- fit$multistart_deviance
			out$multistart_distance <- fit$multistart_distance
		}
		out$multistart_coefficients <- fit$multistart_coefficients
		out$multistart_product_distance <- fit$multistart_product_distance
		out$multistart_selected <- fit$multistart_selected
		out$multistart_converged <- fit$multistart_converged
		out$multistart_criterion <- fit$multistart_criterion
		out$multistart_gaussian_diagnostics <- fit$multistart_gaussian_diagnostics
	}

	if (use_eb) {
		# EB is itself the estimator; nothing is refit or subtracted
		out$bias_correct <- list(
			method = "eb",
			applied = FALSE,
			requested = if (identical(bias_correct, "eb")) "eb",
			default = is.null(bias_correct)
		)
	} else {
		# split-sample jackknife: two half-sample refits combined as 2 * full - mean(halves)
		correction <- .ae_resolve_bias_correct(
			bias_correct,
			family,
			out$non_normal_method,
			prep$Tt
		)
		if (identical(correction, "bootstrap")) {
			if (!identical(family, "binary") || !identical(out$non_normal_method, "irls")) {
				cli::cli_abort(
					c(
						"The parametric-bootstrap bias correction is available for binary IRLS fits only.",
						"i" = "Use {.code bias_correct = \"none\"} or the family default."
					),
					class = "lame_jk_split"
				)
			}
			out <- .ae_boot_bias(out, max_iter, tol)
			if (verbose) {
				n_reps <- .ae_bb_reps
				cli::cli_alert_info("Parametric-bootstrap bias correction: two rounds of {n_reps} simulated refits; uncorrected values in {.code fit$raw}.")
			}
		} else if (identical(correction, "jackknife")) {
			out <- .ae_jackknife(
				out,
				max_iter,
				tol,
				explicit = !is.null(bias_correct)
			)
			if (identical(out$bias_correct$method, "jackknife")) {
				# an unconverged half refit is an error of class "lame_jk_nonconvergence"
				if (verbose) {
					cli::cli_alert_info("Split-sample bias correction: two half-sample refits combined as 2 x full - mean(halves); uncorrected values in {.code fit$raw}.")
				}
			} else if (identical(out$bias_correct$method, "bootstrap") && verbose) {
				n_reps <- .ae_bb_reps
				cli::cli_alert_info("Parametric-bootstrap bias correction (the split-sample correction does not reach this network): two rounds of {n_reps} simulated refits; uncorrected values in {.code fit$raw}.")
			}
		} else if (identical(bias_correct, "none")) {
			# the caller opted out
			out$bias_correct <- list(
				method = "none", applied = FALSE,
				requested = "none"
			)
		} else {
			out$bias_correct <- list(method = "none", applied = FALSE)
		}
	}
	out <- .ae_mark_not_estimated(out, db)
	.ae_inform_not_estimated(out$not_estimated)
	.ae_inform_retirement(out$eb$retirement)
	out
}

# public estimators

#' Fast (MCMC-free) AME estimation for a cross-sectional network
#'
#' @description
#' Fits an additive and multiplicative effects (AME) model to a single
#' cross-sectional network by \strong{iterative block coordinate descent},
#' producing a fast point estimate with no MCMC and no credible intervals. For
#' normal, binary and Poisson data the default is an empirical-Bayes / variational
#' estimator of the model that \code{\link{ame}} samples (see
#' \code{bias_correct}); its \code{vcov()} and \code{confint()} are Laplace
#' approximations.
#'
#' The estimation algorithm adapts the \emph{iterative block coordinate
#' descent} estimator of the Social Influence Regression (SIR) model of
#' Hoff & Minhas (the algorithm implemented in \code{sir::sir_alsfit()} and
#' nicknamed "ALS" in that package) to the AME model. It is a port and
#' adaptation, not original \pkg{lame} methodology.
#'
#' @details
#' The AME decomposition
#' \deqn{z_{ij} = \mu + \beta' x_{ij} + a_i + b_j + u_i' v_j + \epsilon_{ij}}
#' is conditionally linear: it is a linear regression in \eqn{(\mu, \beta, a, b)}
#' for fixed multiplicative factors \eqn{(U, V)}, and the optimal rank-\code{R}
#' \eqn{(U, V)} for fixed \eqn{(\mu, \beta, a, b)} is the truncated SVD of the
#' residual matrix. The estimator therefore cycles, until the residual sum of
#' squares stabilises, through three blocks --- each a monotone
#' (objective-non-increasing) update: a joint least-squares solve for the
#' intercept \eqn{\mu} and regression coefficients \eqn{\beta} (an exactly
#' rank-deficient design is refused with an error of class
#' \code{"lame_design_rank"} that names the aliased columns);
#' an exact solve of the additive effects \eqn{(a, b)} (Gauss-Seidel sweeps
#' started from the direct solution of their normal equations, and an exact
#' eigen solve wherever the sweeps do not reach the normal equations, as on a
#' weakly connected observed-dyad graph; a residual above the tolerance after
#' the exact solve stops the fit with an error of class
#' \code{"lame_als_numerical"});
#' and a weighted low-rank update for the multiplicative factors \eqn{(U, V)}
#' via SVD (eigen-decomposition when \code{symmetric = TRUE}). The weighting
#' makes the multiplicative step correct for unbalanced longitudinal panels
#' (dyads observed at unequal numbers of time points); severely unbalanced
#' panels converge in more iterations. A fixed-effects fit that
#' has not reached its criterion at \code{max_iter} is refitted once with four
#' times that limit, and refused with an error of class
#' \code{"lame_als_convergence"} if it still has not: the message gives the
#' stop (iteration limit, or no step that lowers the objective), the score
#' against its threshold, and the remedy. The last iterate is not returned. An
#' empirical-Bayes fit that does not settle is refused with the same class at once, with the
#' state it stopped in and the remedy.
#'
#' For fixed-effects Gaussian fits (\code{bias_correct = "none"}) and fixed-transform fits, convergence additionally requires a
#' centered first-order SSE score no larger than \code{sqrt(tol)}. A premature
#' plateau or exhausted outer budget receives up to five joint refinements,
#' each bounded by \code{max_iter} optimization iterations. These
#' preserve the same objective and observation mask; a requested single sweep
#' remains unfinished. \code{gaussian_diagnostics} records the initial and final
#' scores, objective values and refinement counts. Multi-start fits retain these
#' checks for every start in \code{multistart_gaussian_diagnostics}; warm refits
#' and the Gaussian bootstrap use the same certificate. A first-order certificate
#' attests a stationary point, a local optimum of the non-convex factor model
#' when \code{R > 0}.
#'
#' With \code{bias_correct = "none"} a \code{family = "normal"} fit is a
#' least-squares (Gaussian maximum-likelihood) fit --- the exact global solution
#' when \code{R = 0}, and a local optimum of the non-convex low-rank objective
#' when \code{R > 0}; the default for a normal response is the empirical-Bayes
#' estimator of \strong{Normal data (the Gaussian empirical-Bayes estimator)}
#' below.
#' For \code{"binary"} and \code{"poisson"}, the estimator uses an IRLS
#' working-response path by default, so coefficients are on the requested
#' link scale (binary probit/logit, Poisson log). The fixed-transform path is
#' still available as a faster exploratory score, but its coefficients are not
#' calibrated effect sizes. For \code{R > 0} the multiplicative factor block
#' of the IRLS families solves a penalized (MAP) sub-problem: each factor row
#' carries the ridge implied by the AME prior \eqn{[u_i, v_i] \sim N(0,
#' \Sigma_{uv})}, with the per-column prior scale estimated by an EM step (the
#' analog of the MCMC \eqn{\Sigma_{uv}} draw). This is what makes the estimate
#' well-defined: the \emph{unpenalized} rank-\code{R} binary MLE does not exist
#' under quasi-separation, and fitting it directly inflates every
#' coefficient (a runaway that grows with \code{R}, with the latent variance
#' diverging). With the penalty the latent variance does not run away.
#'
#' \strong{Empirical-Bayes estimator (the default of binary, Poisson and normal
#' fits).} A binary or Poisson IRLS fit, or a normal fit,
#' with \code{bias_correct = NULL} (or \code{"eb"}) treats the sender and
#' receiver effects as the MCMC \code{\link{ame}} / \code{\link{lame}} do, as
#' draws from a Gaussian prior \eqn{(a_i, b_i) \sim N(0, \Sigma_{ab})}{(a_i, b_i) ~ N(0, Sigma_ab)}
#' (an undirected fit has one effect per actor and a scalar variance, a
#' bipartite fit independent sender and receiver variances), instead of as
#' free fixed effects, and estimates the model by empirical Bayes in two
#' stages. The first is the posterior mode with a Laplace correction, in three
#' parts. (1) The IRLS working problem carries the penalty
#' \eqn{\sum_i \alpha_i' \Sigma_{ab}^{-1} \alpha_i}{sum_i alpha_i' Sigma_ab^(-1) alpha_i}, where \eqn{\alpha_i}{alpha_i}
#' is the actor's effects less the part its node covariates explain (the
#' node-covariate coefficients are fixed effects of the prior mean, as in the
#' MCMC model); the factor block is the penalised ridge above, whose scales are
#' updated by EM from the posterior second moments of the factors, an undirected
#' fit's too. At every
#' reweighting step the intercept, the coefficients and all the effects are
#' solved jointly and exactly from dyad-level sums, in one dense system with no
#' centring step. (2) \eqn{\Sigma_{ab}}{Sigma_ab} is updated by an EM step from the
#' effects and their Laplace covariance, under the inverse-Wishart prior that
#' \code{lame()} puts on it for binary data, so that the penalty is the MCMC
#' prior; without that prior the update drives the variance to zero on sparse
#' networks and with a latent term. (3) The coefficients solve the
#' Laplace-marginal estimating equations instead of the penalised-fit ones: the
#' score of every cell loses half the derivative of the log determinant of the
#' effects' information (the integral over the effects). The penalised fit
#' alone had a slope too small by up to 4.3 percent with 24 percent ties and 6.8 percent with 7 percent
#' ties (30 networks a cell). The latent factors
#' get the same cell correction, so one objective, the penalised deviance plus
#' that log determinant, is descended. This stage is what
#' \code{fit$irls_diagnostics} certifies and what the marginal likelihood
#' \code{fit$eb$log_marginal} (see \strong{Choosing the rank}) is evaluated at.
#'
#' \strong{Variational posterior (the values an empirical-Bayes fit reports).}
#' The mode is the starting point of the second stage, not its result. The mode
#' of a latent factor that the data barely identify lies at the edge of a broad
#' posterior and is too large (it fits noise), which inflates the coefficient of
#' a covariate that follows the product and, one rank above the generating one,
#' the slope. The reported values (\code{a}, \code{b}, \code{U}, \code{V},
#' \code{UVPM}, the coefficients and the variance components) are those of a
#' Gaussian variational posterior of the sender and receiver effects and of the
#' factors, independent across actors and between the sender and the receiver
#' side of the factors; the coefficients and the variance components remain
#' empirical-Bayes point estimates. Under this posterior the predictor of a cell
#' has a mean, the fitted predictor, and a variance, that of the effects plus
#' that of the product,
#' \eqn{\hat u_i' \Omega_j \hat u_i + \hat v_j' \Sigma_i \hat v_j + \mathrm{tr}(\Sigma_i \Omega_j)}{u_i' Om_j u_i + v_j' Sig_i v_j + tr(Sig_i Om_j)}
#' (\eqn{\Sigma_i}{Sig_i} and \eqn{\Omega_j}{Om_j} are the posterior covariances
#' of a sender's and a receiver's factor row), and the score and the information
#' of the likelihood are expected over that Gaussian (nine-point Gauss-Hermite
#' quadrature) instead of taken at its mean. Every equation of the fit uses the
#' expected score and information, and the system of a factor row adds the
#' second moment of the other side (mean product plus covariance). The means
#' and these covariances are the stationary point of one objective, the expected
#' log likelihood less the Kullback-Leibler divergence to the prior, at the scale
#' of the prior of the factors; the step of the iteration is halved when the
#' objective would fall, which leaves the fixed point unchanged. A fit without a
#' latent rank is variational too (the effects only). \code{fit$eb$estimator}
#' says which estimator produced the values (\code{"variational"}, or
#' \code{"laplace"} for the mode alone), \code{fit$eb$variational} records the
#' iteration (convergence, stop reason, number of iterations, the largest step a
#' full update would still take, the lower bound and its parts, and the scale of
#' the prior of the factors) and \code{fit$eb$factor_posterior} the covariances.
#' A fit whose mode or whose iteration stops short of its fixed point is not
#' returned: the call stops with an error of class
#' \code{"lame_als_convergence"} that reports the state and the remedy (a
#' smaller rank, a larger \code{max_iter}, or \code{\link{ame}} /
#' \code{\link{lame}}); the values of an iteration that has not settled are
#' never reported as the estimate, nor replaced by the mode.
#'
#' \strong{The scale of the prior of the factors.} The rows of the factors have
#' a Gaussian prior whose covariance the data inform weakly where there are few
#' events per loading. The fit estimates it by the variational update of the
#' inverse-Wishart prior \code{ame()} and \code{lame()} put on it (\code{kappa0}
#' degrees of freedom and scale matrix \code{kappa0 * Suv0}, the data-scaled default of
#' \code{\link{ame}}), the inverse of the
#' posterior mean precision,
#' \eqn{(\Psi + SE) / (\mathrm{df} + n)}{(Psi + SE) / (df + n)}, with
#' \eqn{SE} the second moments of the stacked sender and receiver rows. What
#' second moments to use decides the result. The mean-field covariance of a row
#' counts the uncertainty of the other side twice (the second moment of the
#' other side in its system is the mean product plus the covariance) and is below
#' the covariance of the exact posterior; the update evaluated there drives the
#' scale to the floor of its hyperprior wherever the data are weak, and the
#' product and the coefficients that follow it with it. Against an exact Gibbs
#' sampler of the same model run at the same scale (fixed hyperparameters, 8
#' networks of 32 to 100 actors), the second moments of the rows were 4 to 19
#' percent below the sampler's at the mean-field covariances (on average 9.5).
#' The default evaluates the update at the cavity covariances,
#' \eqn{\Sigma_i = (\Lambda + \sum_j \nu_{ij} \hat v_j \hat v_j')^{-1}}{Sig_i = (Lam + sum_j nu_ij v_j v_j')^-1},
#' where \eqn{\nu_{ij}} is the cavity weight of the cell, the curvature of the log
#' likelihood integrated over its Gaussian predictor (the information the cell
#' adds about the row once the rest of the predictor is integrated out), and
#' these were from 9 percent below to 5 percent above the sampler's (on average
#' 3 below); they are the covariances
#' \code{fit$eb$factor_posterior} reports (the covariances of the fixed point of
#' the updates, from which the lower bound and a refit start, are kept beside
#' them as \code{Sig_mf} and \code{Om_mf}). The covariance of the stacked rows of a
#' directed network is isotropic, one scale and one correlation between the
#' sender and receiver position of an actor (the sampler's prior is a full
#' \eqn{2R \times 2R} matrix: the scale carries one dimension of its
#' uncertainty, see below); a bipartite fit has the priors of its sampler, the
#' receiver factor at unit scale and, for the sender columns, the scales of
#' \eqn{P = U G V'} with standard normal \eqn{G}: the sender rows seen by the
#' likelihood have covariance \eqn{G'G}, a Wishart matrix with \eqn{R} degrees of
#' freedom, and the scales of the columns are its ordered eigenvalues, which
#' repel each other (the density of the eigenvalues carries the product of their
#' differences), so a column the data say nothing about falls to the small
#' eigenvalues of that matrix rather than staying at the scale of the prior; and an undirected fit the
#' prior of its sampler on each column (\eqn{U L U'} with \eqn{L_r \sim N(0, 1)} and
#' column precision \eqn{\psi_r} (gamma, shape 1, rate 1/2), so that a column of
#' the factor seen by the likelihood has variance \eqn{|L_r| / \psi_r}). The scale
#' is recorded in \code{fit$eb$variational$scale} (the rule, the scale, the
#' correlation of a directed fit, and the standard deviation of its logarithm).
#' \code{options(lame.eb_scale = "plugin")} selects the update at the mean-field
#' covariances with a full \eqn{2R \times 2R} covariance and the closed-form
#' rescaling of the factors that it needs, the fit of the earlier release. The
#' option, and \code{options(lame.eb_scale_uncertainty)} below, also govern the
#' static factors of a normal, binary or Poisson dynamic fit
#' (\code{lame(method = "als")}), not the persistence of its moving positions.
#'
#' \strong{Dyadic correlation.} A directed binary probit fit with the
#' empirical-Bayes estimator and \code{dcor = TRUE} (the default for directed
#' networks) takes the two cells of a dyad as one observation: the latent errors
#' of \eqn{y_{ij}} and \eqn{y_{ji}} are bivariate normal with correlation
#' \eqn{\rho}, and the likelihood of the pair is the bivariate normal
#' probability \eqn{\Phi_2(s_1 \eta_{ij}, s_2 \eta_{ji}; s_1 s_2 \rho)},
#' \eqn{s = 2y - 1}. The variational posterior is taken in that model. The
#' expected log likelihood of a pair, its score and its information (a
#' \eqn{2 \times 2} matrix: the cells of a dyad are not independent given the
#' predictors) are computed under the Gaussian posterior of the two predictors,
#' which are correlated through the effects of the two actors, by a
#' Gauss-Hermite rule of five nodes a side. The equations of the coefficients and
#' of the actor effects carry the cross information of the pairs, the rows of the
#' factors see the exact score of the pair, the cavity weights that set the scale
#' of the factors are those of the pair probability integrated over the predictors
#' (a bivariate normal probability in closed form), and \eqn{\rho} is estimated
#' with the other hyperparameters: at every iteration it moves to the posterior
#' mean of the quadratic approximation of the expected log likelihood in
#' \eqn{\rho} under the prior of the samplers, \eqn{2\,\mathrm{Beta}(2, 2) - 1}
#' (density proportional to \eqn{1 - \rho^2}), the one \code{\link{ame}} and
#' \code{\link{lame}} share and the fast fit takes from the same source. The
#' value is \code{fit$VC["rho"]}, its posterior standard
#' deviation and the method are in \code{fit$rho_estimation}, and the
#' record of the pair model (the correlation, its score and curvature, the
#' number of pairs, the quadrature) is \code{fit$eb$pair}. Against the
#' sampler (80 simulated networks of 60 and 100 actors, one and four periods, ranks
#' 0 and 2, with correlations of 0, 0.3, 0.6 and 0.8) the correlation is within 0.02
#' of the sampler's posterior mean on average and within one posterior standard
#' deviation (0.04 to 0.1) in 79 of the 80 networks, where the post hoc correlation of
#' the independent-cell fit is 0.02 to 0.08 too high; the slope of a covariate is
#' within 0.2 sampler standard deviations of the sampler's on average, and its
#' standard error 1.01 to 1.02 times the sampler's standard deviation, with or
#' without the pair model. The factors of an independent-cell fit absorb reciprocity,
#' which the pair model attributes to \eqn{\rho}. The cost is the pair quadrature, in
#' proportion to the number of dyads and periods: 1 to 2.5 times the independent-cell
#' fit for one period (60 to 150 actors, ranks 0 and 2) and 2 to 4.5 times for four.
#'
#' \strong{Counts (Poisson).} A Poisson fit with the default \code{bias_correct}
#' (or \code{"eb"}) is the same estimator for the model that \code{\link{ame}} and
#' \code{\link{lame}} sample: a count is
#' \eqn{y_{ij} \sim \mathrm{Poisson}(e^{z_{ij}})}{y_ij ~ Poisson(exp(z_ij))} with
#' \eqn{z_{ij} = \eta_{ij} + e_{ij}}{z_ij = eta_ij + e_ij} and
#' \eqn{e_{ij} \sim N(0, v_e)}{e_ij ~ N(0, ve)}, where \eqn{\eta} is the predictor of the
#' sections above (coefficients, sender and receiver effects, latent product). The
#' fixed-effects Poisson fit (\code{bias_correct = "none"}, the default of the earlier
#' releases) has no finite estimate for the effect of an actor without a count on a side,
#' which runs to minus infinity (-17.5 after 71 reweighting steps on a network of 30 actors
#' with one empty row and one empty column), and it fits the marginal mean with a post hoc
#' overdispersion. The empirical-Bayes fit gives such an actor a finite effect, shrunk
#' towards the prior mean as in the sampler. The effects have the Gaussian prior of the
#' binary fit, with the scale and the degrees of freedom that \code{ame()} and
#' \code{lame()} give counts: the variance across actors of the centred log row and column
#' rates (a half-count pseudo-observation keeps an empty actor finite) and
#' \eqn{\mathrm{round}(4 + 3 n / 100)} degrees of freedom, \eqn{n} the number of actors;
#' the factors have the priors of the binary fit. The variance of the noise \eqn{v_e} is a
#' hyperparameter, estimated as the covariance of the effects is, by type-II maximum
#' likelihood: it is the value that maximises the Laplace marginal likelihood of the model
#' (\code{fit$eb$log_marginal}, the mode of the effects and the factors solved at that
#' value) plus the log density of the prior that \code{ame()} and \code{lame()} give it
#' (the inverse-gamma of one degree of freedom on the scale one). It is
#' \code{fit$VC["ve"]} and \code{fit$eb$ve}, and the objective it was found on is
#' \code{fit$eb$ve_search$curve}. The stationarity condition of the maximum is the EM step
#' that \code{ame()} and \code{lame()} draw \eqn{v_e} from,
#' \eqn{(1 + \sum E[e^2]) / (N - 1)}{(1 + sum E[e^2]) / (N - 1)} (the sum running over the
#' observed cells, the dyads once for an undirected network, of the second moment of the
#' noise given the count), with the second moment taken, to the order of the Laplace
#' approximation, under the posterior of the predictor of the cell given all the
#' parameters; the step at the mode alone, or under the mean-field posterior of the
#' predictor, which is shrunk where the data inform the latent product weakly, is too
#' large (a count network with fewer than three nonzero counts per loading: 1.3 to 4.6 posterior
#' standard deviations above the sampler's, where the maximum of the marginal likelihood
#' is within 1). A count network without noise beyond the Poisson has a positive estimate, as the sampler's posterior
#' does: the inverse-gamma prior holds \eqn{v_e}{ve} off zero (on four plain Poisson panels of 22 actors and five periods the
#' estimate was 0.034 to 0.046, within 0.0015 of the sampler's posterior mean; 0.040 on average over 60 panels). The maximum is found by continuation
#' on \eqn{\log v_e}: the mode is solved at 0.05, where every column of the latent product
#' is alive, and then at larger values until the objective falls, then at the vertex of the
#' parabola through the three points about its maximum, each warm from the nearest solved
#' mode (typically 8 to 10 modes, 5 to 30 reweighting steps each); a mode solved cold at a
#' large \eqn{v_e} loses a weakly supported column and keeps it lost at every smaller
#' value, on a lower branch of the marginal likelihood. The variational stage holds
#' \eqn{v_e} at the maximum, and a column the mode lost is restarted from the score of the
#' cells before it iterates. The cell has
#' everything the binary cell has. Its noise is integrated out, so that the log likelihood,
#' the score and the information of a count at the predictor \eqn{\eta} are those of the
#' Poisson-lognormal distribution (the information is positive: the log likelihood is
#' concave in \eqn{\eta}); they are the cumulants of the posterior of the noise given the
#' count, computed by an adaptive Gauss-Hermite rule centred at its mode (14 to 64 nodes,
#' by the size of \eqn{v_e}), with the derivative of the information that the Laplace
#' correction needs. The variational stage takes the expectation of the cell over the
#' Gaussian posterior of the predictor with the rule of the binary fit (nine nodes where
#' the predictor is uncertain, fewer where its variance is below 0.3), and the cavity
#' weight of a cell is the information of the Poisson-lognormal cell whose noise variance
#' is \eqn{v_e} plus the variance of the predictor. The linear predictor
#' \code{EZ} is on the latent log scale of the model (there is no shift by half the noise
#' variance), the mean of a count is \eqn{e^{\eta + v_e / 2}}{exp(eta + ve / 2)}, which
#' is what \code{fit$fitted}, \code{fit$YPM}, \code{fitted()} and \code{predict()} hold (the
#' plug-in mean, as the samplers report it; the posterior predictive mean of a count is not
#' finite in this model, and \code{plug_in = FALSE} returns the value the integration over the
#' predictor gives, \eqn{e^{m + (v + v_e) / 2}}{exp(m + (v + ve) / 2)}).
#' A network without a count above zero has no finite log rate, and a dyadic covariate
#' of one sign that is non-zero only where every count is zero has no finite coefficient:
#' both are refused (class \code{"lame_als_separation"}). Rank evidence
#' (\code{fit$eb$log_marginal}) is the Laplace value of the Poisson-lognormal model.
#'
#' \strong{Dyadic correlation of counts.} A directed count network with
#' \code{dcor = TRUE} (the default for directed networks, as in \code{\link{ame}} and
#' \code{\link{lame}}) takes the two counts of a dyad as one observation at the fitted
#' state: the noises \eqn{(e_{ij}, e_{ji})} are bivariate normal with variance \eqn{v_e} and
#' correlation \eqn{\rho}, and the pair has the bivariate Poisson-lognormal probability. The
#' effects, the factors and \eqn{v_e} are fitted with independent cells (they agree with the
#' sampler's, see below); the pair stage then gives (1) \eqn{\rho}, the posterior mean, under
#' the prior the samplers share, of the profile of the pair probability at the fitted predictors,
#' reported as \code{fit$VC["rho"]} with its standard deviation in \code{fit$rho_estimation}, and (2)
#' the information of the pairs: the pair probability integrated over the Gaussian posterior
#' of the two predictors is the pair probability with the noise covariance plus their
#' covariance, and its curvature in the two means (two diagonal weights and a cross weight,
#' negative for a positive correlation) enters the covariance of the coefficients, so that the
#' standard error of a covariate that is the same in the two cells of a dyad is not understated
#' (\code{fit$eb$pair} records the correlation, its precision, the number of pairs and the
#' cross weights). The pair stage costs 13 to 22 evaluations of the pair probability of every
#' dyad; \code{dcor = FALSE} leaves it out and keeps the post hoc moment estimate of the
#' correlation of the latent errors (\code{fit$rho_working} holds it in either case). Undirected
#' and bipartite networks have no pair of cells.
#'
#' \strong{Evidence for counts.} Against the sampler (29 simulated count networks of 30 to
#' 50 actors, one and four periods; directed, bipartite and undirected, ranks 0 to 3, mean counts
#' 0.3 to 31, noise variances 0.05 to 2, empty actors, node covariates and reciprocity 0 to 0.8;
#' \code{ame()} / \code{lame()}, 30,000 scans after the burn-in) the slopes are within
#' 0.35 posterior standard deviations of the posterior means (median 0.05; the fixed-effects fit
#' 2.8 and 0.6) and the intercepts within 0.38, the standard errors are 0.99 to 1.09 times the
#' posterior standard deviations of the slopes (mean 1.04) and 0.97 to 1.02 times those of the
#' intercept of a cross-section, \code{va} and \code{vb} are within 0.3 posterior standard
#' deviations and \code{ve} within 0.5 (median 0.07; the fixed-effects fit 6.9 and 1.7), the actor
#' effects differ from the posterior means by a median 0.029 of their standard deviation (the
#' fixed-effects fit 0.17, and 2.3 with empty actors), and the correlation of a directed network is
#' within 0.027 of the sampler's (median 0.010 and at most 1.2 of its posterior standard
#' deviations). Without the pair stage the
#' standard error of a slope was 0.77 to 0.87 of the sampler's for dense counts with a covariate
#' that is the same in both cells of a dyad and a reciprocity of 0.6 to 0.8. On 1,500 simulated
#' networks (25 designs of 60: 30 to 60 actors, mean counts 0.1 to 20, noise variances 0.1 to 1.5,
#' ranks 0 to 3, directed, bipartite and undirected, one and four periods, reciprocal designs
#' with a symmetric covariate, and ranks above the generating one) the fit converged in every
#' network (the fixed-effects fit in 94.9 percent, 67 percent at worst), the slope had a bias of
#' 0.6 percent on average (4.7 at most, with 30 actors, mean count 0.1 and noise variance 1.5; the
#' fixed-effects fit -1.8, and -12.1 at worst, at rank 3) and its 95 percent interval covered the
#' truth in 95.5 percent of the networks (88 to 100 by design); its sampling standard deviation was
#' 0.76 times the fixed-effects fit's (median) and the root mean squared error of the actor effects
#' 0.79 times (0.18 at the sparsest). The central 50 and 90 percent predictive intervals of 17,150
#' held-out cells (10 percent of the cells of 70 networks of seven designs) covered 51 to 54 and 90
#' to 91 percent of the observed counts, with no departure of the randomised probability integral
#' transform from uniformity except at rank 3 with 3.0 nonzero counts per loading (Kolmogorov-Smirnov
#' p = 0.007). The lowest rank within one unit of the largest Laplace value was the generating rank
#' in 142 of 144 networks of 50 actors, and the rank of the largest value in 143 (each rank is
#' compared at the noise variance that maximises its own marginal likelihood).
#'
#' \strong{Counts with little information about the latent product.} With fewer than about
#' three nonzero counts per factor loading (\code{fit$eb$variational$information["per_loading"]}) the
#' variational posterior shrinks the product. Below 2.1 units of the information-weighted signal per loading
#' (\code{fit$eb$variational$information["signal_per_loading"]}) the fit withholds the product: \code{fit$U},
#' \code{fit$V} and \code{fit$UVPM} are \code{NA}, and \code{print()} states the statistic and the threshold (see
#' \strong{Choosing the rank of a normal fit}). The noise variance does not
#' take up what the shrunk product leaves: on 22 networks of 40 and 50 actors at ranks 2 and 3
#' (1.1 to 6.6 nonzero counts per loading) \code{ve} is within 0.41 posterior standard deviations of
#' the sampler's (mean 0.18; the step of the variational iteration itself was above it by 1.3 to 4.6,
#' 2.6 on average, below three nonzero counts per loading), on 60 networks of 50 actors at rank 3 and 2.8 nonzero
#' counts per loading it averaged 0.412 for a truth of 0.40 (0.406 at 4.2), and at a rank above the
#' generating rank it is within 0.6 posterior standard deviations of the sampler's (32 networks; 4 percent below it on average). The variational
#' posterior shrinks a product the data determine a little more than the sampler does: its norm is 0.95 to 0.98
#' of the norm of the sampler's posterior mean product (\code{lame()}, whose prior of the factors is
#' the fit's) from 2.6 nonzero counts per loading and 0.80 to 0.90 up to 2.1 (against \code{ame()},
#' whose prior of the factors is tighter: 0.89 to 0.98 and 0.65 to 0.85); it is as close to the
#' generating product as the sampler's posterior mean is (root mean squared difference 0.99 times the
#' sampler's on average, 0.92 to 1.06 on the 22 networks), the mean of the predictor of a cell differs
#' from the sampler's by 0.15 of the sampler's posterior standard deviation (root mean square over the
#' cells, averaged over the networks; at most 0.23), and its variance is 0.92 of the sampler's (0.74
#' to 1.03). With about one nonzero count per loading or fewer the iteration still settles: on 54 networks of 40 actors at
#' ranks 1 to 3 with 0.3 to 3.6 nonzero counts per loading every fit converges, none collapses a column of its product, and one, at rank 3
#' with 0.6 nonzero counts per loading, retires a column (see \strong{A rank the data do not support}); \code{ve} averages
#' 0.39, 0.35 and 0.39 at ranks 1, 2 and 3 for a truth of 0.40. The product correlates with the truth at 0.83 on average at
#' rank 3 and 2.8 nonzero counts per loading (the fixed-effects fit 0.71).
#'
#' \strong{Normal data (the Gaussian empirical-Bayes estimator).} A normal fit with the default \code{bias_correct}
#' (or \code{"eb"}) is the same estimator for the model that \code{\link{ame}} and \code{\link{lame}} sample,
#' \eqn{y_{ij} = \mu + \beta' x_{ij} + a_i + b_j + u_i' v_j + e_{ij}}{y_ij = mu + beta' x_ij + a_i + b_j + u_i' v_j + e_ij} with
#' \eqn{e_{ij} \sim N(0, s^2)}{e_ij ~ N(0, s2)}, and it is the estimator of the binary and count fits with the one cell in
#' which everything is exact. The cell of a normal response has the working response \eqn{y} itself, the weight
#' \eqn{1 / s^2}{1 / s2} and the score \eqn{(y - \eta) / s^2}{(y - eta) / s2}; its information does not depend on the predictor, so
#' the Laplace correction of the coefficients vanishes and the Laplace approximation of the coefficients and the effects is
#' not an approximation: at \code{R = 0} the fit is the posterior mean (the best linear unbiased prediction) of the
#' linear Gaussian model at the fitted hyperparameters, and \code{vcov()} is its generalised least squares covariance
#' (checked against an independent computation to 1e-6 in the tests). The effects have the Gaussian prior
#' \eqn{(a_i, b_i) \sim N(0, \Sigma_{ab})}{(a_i, b_i) ~ N(0, Sigma_ab)} and the inverse-Wishart prior on \eqn{\Sigma_{ab}}{Sigma_ab}
#' that \code{ame()} and \code{lame()} give a normal response (the scale is the mean variance of the row and the column means
#' of the centred response, \eqn{\mathrm{round}(4 + 3 n / 100)} degrees of freedom), the coefficients are flat, and the node
#' covariates are fixed effects of the prior mean of the effects. The variance of the noise \eqn{s^2}{s2} is a hyperparameter with the
#' sampler's inverse-gamma prior of one degree of freedom on the scale \code{prior$s20}, by default the variance of the
#' observed responses. \eqn{s^2}{s2} and \eqn{\Sigma_{ab}}{Sigma_ab} are the fixed points of the EM images that the samplers draw from
#' (the posterior mean of the inverse-gamma and inverse-Wishart full conditionals), formed from the posterior second
#' moments, so that they are net of the degrees of freedom the effects and the factors use:
#' \eqn{s^2 = (s_{20} + \sum E[e^2]) / (N - 1)}{s2 = (s20 + sum E[e^2]) / (N - 1)} with
#' \eqn{E[e^2] = (y - m)^2 + v}{E[e^2] = (y - m)^2 + v} (\eqn{m} and \eqn{v} the mean and variance of the predictor of the cell
#' under the posterior). The factors of a latent rank are the variational posterior of the binary fit, with its cavity rule for the scale of their
#' prior and the priors of the samplers (a directed network the shared prior of \code{ame()} and \code{lame()}, whose scale is
#' that of the effects). The fitted values are the posterior means of the predictor; the response is Gaussian, so they are
#' also the predictive means.
#'
#' \strong{Dyadic correlation of normal data.} A directed normal network with \code{dcor = TRUE} (the default) takes the two
#' cells of a dyad as one bivariate normal observation with covariance
#' \eqn{s^2 [[1, \rho], [\rho, 1]]}{s2 [[1, rho], [rho, 1]]}. Everything is in closed form: the information of a pair is the constant matrix
#' \eqn{[[1, -\rho], [-\rho, 1]] / (s^2 (1 - \rho^2))}{[[1, -rho], [-rho, 1]] / (s2 (1 - rho^2))}, whose cross weight enters the equations of
#' the coefficients and of the effects (the estimates are generalised least squares), the rows of the factors see the exact score of
#' the pair, the cavity weights of the factor scale are those of the pair density integrated over the predictors, and
#' \eqn{\rho}{rho} takes its posterior mean under the prior of the samplers (density proportional to \eqn{1 - \rho^2}{1 - rho^2}), computed on a
#' fine grid from the sufficient statistics of the pairs and not by a quadratic approximation. The variance of the noise is estimated from
#' the whitened pairs. The value is \code{fit$rho} (\code{fit$VC["rho"]}), its posterior standard deviation and method are in
#' \code{fit$rho_estimation}, and the record of the pair model in \code{fit$eb$pair}. \code{dcor = FALSE} fits independent
#' cells and reports the correlation of the residuals of the pairs as before; undirected and bipartite networks have no
#' pairs.
#'
#' \strong{Standard errors of normal fits.} \code{vcov()} is the Laplace covariance of the coefficients jointly with the
#' effects and the factors, built from the information of the Gaussian cell: the weight \eqn{1 / s^2}{1 / s2} and the cross weight of the
#' pairs, not the cavity weight that sets the scale of the factors. The information of a Gaussian cell is the same whatever
#' the predictor, so integrating the variance of the predictor into the weight counts it twice: with the cavity weights the standard errors of the
#' dyadic coefficients were 1.05 to 1.13 times the sampler's posterior standard deviations at ranks 1 to 3, with the
#' information of the cell 0.97 to 1.02. Node-covariate and intercept coefficients are covered too, and the part of
#' the uncertainty that comes from the scale of the prior of the factors is added as for the other families.
#' Beyond the designs above the intervals keep their level: on 1,000 simulated normal panels (200 a cell: 40 actors with 8 and
#' 12 periods at ranks 0 and 2, and 120 actors with 4 periods at rank 2; a dyadic covariate, a row and a column covariate, dyadic
#' correlation 0.3) the 95 percent Laplace intervals covered the dyadic coefficient in 0.943 of the panels (cells 0.915 to
#' 0.965), the row and column covariates in 0.959 and 0.956 (0.935 to 0.975) and the intercept in 0.950 (0.935 to 0.955); the Monte Carlo
#' standard error of a pooled figure is 0.007 and of a cell 0.015 to 0.020.
#'
#' \strong{Choosing the rank of a normal fit.} \code{fit$eb$log_marginal} is the Laplace approximation to the log marginal
#' likelihood of the independent-cell model at the posterior mode, solved from the variational state with the variance of
#' the noise held at its fitted value and to a scaled score of 0.01 (the log marginal likelihood moves by at most
#' \eqn{m s^2 / 2} from the mode, 0.05 for the 1,000 effects and factor entries of 170 actors at rank 2); it is exact at
#' \code{R = 0}, and it is the value of \strong{Choosing the rank}. The mode is the one \code{fit$irls_diagnostics} certifies.
#' \code{fit$eb$variational$information} measures the signal of the fitted product in noise variances per loading (for a binary
#' or count fit, \code{signal_per_loading}: the product weighted by the information of each cell, which is the same quantity
#' with the fisher weight of the cell in place of \code{1 / ve}). Below 2.1, in every family,
#' the data do not determine a product of that rank: the sampler's posterior mean product is shrunk towards zero and the
#' variational mean is not, so its size is not the posterior's. Such a fit withholds the product: \code{fit$U}, \code{fit$V}
#' (and \code{fit$L}), \code{fit$UVPM} and the positions of \code{\link{latent_positions}} are \code{NA}, \code{\link{uv_plot}}
#' stops with an error of class \code{"lame_product_withheld"}, and \code{fit$eb$product_withheld} records the reason, the
#' statistic, the threshold and the log marginal likelihood of every rank from 0 to the fitted one with the rank it chooses,
#' which \code{print()} and \code{summary()} state. \code{fitted()}, \code{predict()} and \code{residuals()} are those of
#' the fit (the predictor of the variational posterior), and so are the coefficients, their standard errors and the variance
#' components. The parametric bootstrap draws its replicates from the fit and reports no interval for the positions.
#' The dynamic estimator of a normal panel (\code{lame(method = "als")} with \code{dynamic_ab}, \code{dynamic_beta} or
#' \code{dynamic_uv}) records the same statistic in \code{fit$eb$information} and withholds in the same way, with the same
#' threshold for static factors and 1.3 per loading of a period for moving positions (the products of every period:
#' \code{fit$UVPM}, \code{fit$Omat}, \code{fit$Oarr} and the positions with their standard deviations); it has no rank
#' evidence, so its record has none.
#'
#' \strong{Scale of a normal response.} The iteration works in the units of the response: its tolerance is relative to the
#' standard deviation of the response and the floors of the variances are relative to its variance. The priors of the factors
#' are the samplers' and, like them, written in the unit of the response (\code{R/normal_scale.R}): the rows of a directed
#' network have the typical size of the additive effects (the square root of the mean variance of the row and column means of
#' the centred response) as the scale of their prior, the interaction matrix of a bipartite network and the weights of an
#' undirected one have the variance of the response as the variance of theirs. So the fit of \code{c * y} is the fit of
#' \code{y} in other units: on simulated networks of every kind at ranks 0 and 2 the coefficients, the variance components,
#' the effects and the product of \code{c * y} are \code{c} (or \code{c^2}) times those of \code{y} to 1e-13 for \code{c} of
#' 1e-3 and 1e3, an undirected panel without factors to 5e-5 (its stop rule), and the standard errors, which carry two refits
#' of the scale of the factors, to 1e-3. A response
#' without variation, or with a standard deviation outside 1e-6 to 1e8, is refused with the remedy (classes \code{"lame_als_data"} and
#' \code{"lame_als_scale"}).
#'
#' \strong{The earlier normal fit.} \code{bias_correct = "none"} returns the fixed-effects fit of the earlier releases
#' unchanged: the least-squares fit by block coordinate descent (the global solution at \code{R = 0}, a local optimum of the
#' non-convex low-rank objective otherwise), with the conditional sandwich covariance of the intercept and the
#' dyadic coefficients, which holds the additive and multiplicative effects fixed.
#'
#' \strong{Validation of normal fits.} Against the samplers' posterior means (55 simulated networks of 30 to 50 actors a side,
#' directed, undirected and bipartite, ranks 0 to 3, 38 against \code{\link{ame}} and 17 against \code{\link{lame}} fits, 12 of one
#' period and 5 of four periods, of 20,000 scans after 3,000 of burn-in), the 245 coefficients lie within 0.28 posterior standard
#' deviations (mean 0.04), the Laplace standard errors are 0.91 to 1.05 times the posterior standard deviations (mean 1.00;
#' the posterior standard deviation of the intercept of a \code{lame()} fit is that of the level of the centred effects, 2.4 to
#' 5.7 times smaller than the standard error of the population intercept, and is left out of this comparison), the
#' variance of the noise is within 0.95 posterior standard deviations, the variances of the sender and receiver effects within 0.16 and
#' the dyadic correlation within 0.45. On 103 simulated networks of 40 actors a side (the 55 above with 48 more of weak and
#' moderate signal; 89 with a latent rank, 1 to 3, directed, undirected and bipartite; against \code{\link{ame}} and
#' \code{\link{lame}} fits of 20,000 scans) the norm of the double-centred product is 0.993 to 1.146 of the sampler's posterior
#' mean product where the fit has at least 2.1 noise variances of signal per loading (39 networks; correlation of the products
#' 0.984 to 1.000). Below (50 networks) the undirected fits are 1.37 to 6.9 times the sampler's (0.31 to 1.73 per loading),
#' the directed and bipartite ones 1.08 to 5.0 below 0.7 per loading and 1.01 to 1.25 from 0.7: the sampler's posterior mean is
#' shrunk towards zero and the mean field is not. Every network from 1.99 per loading has a norm within 0.9 to 1.15 of the
#' sampler's, and the largest that does not is at 1.73 (pooled with the binary and count networks below, the threshold is 2.1).
#' Evaluating the update of
#' the scale of the prior at the posterior mean of the inverse-Wishart (\code{options(lame.eb_scale_stat = "mean")}) or at the
#' mean-field second moments does not close the gap (median 2.25 and 1.74 times the sampler's below one per loading, against
#' 2.14). The coefficients of the fits below 2.1 lie within 0.30 posterior standard deviations of the sampler's, with standard
#' errors 0.93 to 1.05 times its posterior standard deviations. On 118 binary and count networks of 40 actors a side with a
#' latent product (directed, undirected and bipartite, ranks 1 and 2, against \code{\link{ame}} fits of 20,000 scans; 71
#' binary networks of 8 and 25 percent ties and 47 Poisson-lognormal ones) the training events per loading do not measure
#' the information about the product: a bipartite binary fit with 5.6 events per loading had 4.4 to 5.7 times the norm of the
#' sampler's product, and a count fit with 13 nonzero counts per loading 7.5 times it. The information-weighted signal does:
#' below 2.1 (79 networks) the binary products are up to 6.2 times the sampler's (median 1.47; undirected up to 6.2,
#' bipartite 5.7, directed 2.9) and the count products up to 7.5 (median 1.87), and these are withheld; from 2.1 the norm is
#' 0.83 to 1.00 times the sampler's for a binary fit (19 networks: the variational posterior shrinks the product somewhat
#' more than the sampler does) and 0.94 to 1.12 for counts (20 networks). The largest network outside the rule is at 2.03
#' per loading and the next at 2.17. In repeated sampling (21 designs of 100 to 200 networks,
#' 4,100 networks and 19,700 coefficient fits: every kind, ranks 0 to 3, 20 to 80 actors, one and four periods, 20 percent
#' missing cells, a rank above the generating one and a weak product) the 95 percent Laplace intervals cover 95.2 percent
#' (dyadic coefficients 95.1, node covariates 95.9, intercept 94.3; designs from 94.1 to 96.7). The rank with the largest
#' \code{fit$eb$log_marginal} is the generating rank in all of the networks (252 networks of 40 actors, generating ranks
#' 0 to 2) with at least 8 noise variances of signal per loading, in 95 percent from 4 to 8, 89 percent from 2 to 4, 50
#' percent from 1 to 2, and in none below 1, where it picks a smaller rank; it picked a larger rank than the generating one in
#' one network (a rank-0 network). Of 300 random designs (10 to 60 actors, one and three periods, up to 30 percent missing
#' cells, fitted ranks 0 to 3 against generating ranks 0 to 2), all 300 fit: the two that the iteration did not settle at rank 3
#' (undirected panels of 14 and 10 actors over generating ranks of 0 and 1, a limit cycle of the scales of the columns) retire two
#' columns (the rank-1 fit) and converge at rank 3 under the priors in the unit of the response.
#'
#' \strong{Speed.} The quadratures of the empirical-Bayes fit (the expected
#' working problem of a cell and of a pair of cells, the cavity weights, the
#' bivariate normal probability of the pair likelihood, the log-domain integral
#' for its small values) and the sums over the panel that its equations need run
#' in compiled code; the reference implementations are kept in R and agree with
#' the kernels to rounding (typically 1e-12, and 1e-9 at worst on realistic
#' predictors). The kernels split their cells and pairs over OpenMP threads where
#' the platform has them: \code{options(lame.als_threads = k)} sets the
#' number of threads (default 1), and the estimates do not depend on it. Do not
#' combine threads with forked workers such as \code{parallel::mclapply()}. The
#' dense solves of the effect step (a Cholesky inverse of the stacked
#' sender and receiver effects, 600 by 600 at 300 actors) are the cost of a large
#' panel and run in the BLAS of the R installation, whose speed sets that cost.
#'
#' \strong{Standard errors.} \code{vcov()} and \code{confint()} return the
#' Laplace approximation to the covariance of the coefficients: the inverse
#' penalised information of all the parameters (coefficients, effects and
#' factors) at the cavity weights (for a normal fit the information of the Gaussian cell,
#' \strong{Standard errors of normal fits}), with the second moments of the factors in the
#' curvature of the other side, plus the part of the uncertainty that comes from
#' the scale of the prior of the factors. The logarithm of the scale has, to the
#' second order, a normal posterior whose precision is the information of the
#' data about it (the derivative of the log marginal likelihood in the log scale,
#' \eqn{(\mathrm{tr}(\Sigma^{-1} SE) - n p) / 2} at the cavity second moments,
#' from two refits at the scale times \eqn{e^{\pm 0.15}}) plus that of the
#' hyperprior; the coefficients move with it by \eqn{J = d\theta / d \log s}, and
#' the covariance gains \eqn{J J' / c}. This is one dimension of the uncertainty of the
#' scale (the sampler's covariance of the rows is a full matrix, with 10 free
#' parameters at \eqn{R = 2}); \code{options(lame.eb_scale_uncertainty = FALSE)}
#' leaves it out.
#'
#' \strong{A rank the data do not support.} A sparse binary or count network at a rank above
#' what it supports (or a network with no signal at all) leaves the iteration at that rank
#' with no stationary point the data determine: a column the data barely inform keeps a
#' live scale, its rows are fitted at the scale of the prior and their second moments are the
#' prior's, so that the posterior mode drifts (the largest score of the stationarity
#' certificate rises with every accepted step) or the variational iteration stops improving
#' at a residual of a few thousandths, or enters a limit cycle in the scales of the columns.
#' Such a fit is not refused. The five-node rule of the
#' pair likelihood of a directed binary network integrates a dyad to less than the gain of a step
#' near the fixed point where the rows are weakly informed, which is a limit of the
#' quadrature as often as of the data: a stall of the pair rule is first refined, to seven
#' and then to nine nodes a side, from the state it reached (\code{fit$eb$pair$refined}). An
#' iteration that is still not settled retires its weakest column, the one with the smallest
#' norm of the product \eqn{\|u_k\| \|v_k\|} (the smallest \eqn{|l_k|} of an undirected fit),
#' and the fit of the next rank down starts from the other columns; this repeats while the
#' iteration does not settle. The fit is the one of the rank it reached: \code{fit$R} is that
#' rank, \code{fit$R_requested} the rank asked for, \code{fit$eb$retirement} the record, and
#' the evidence, the standard errors and the factors are those of the effective rank;
#' the fit says so with a message of class \code{"lame_als_retirement"} when it is made, and
#' \code{print()} and \code{summary()} say what was retired. A fit that still does not settle
#' with no column left is refused with an error of class \code{"lame_als_convergence"} that
#' names the rank asked for, and so is one that ran out of a \code{max_iter} below the default
#' of 200 (a bound the caller chose: the iteration may only need more steps, so it is not
#' reduced in rank). Three sparse binary networks of 50 actors at ranks 2 and 3,
#' \code{YX_bin} at rank 4 and a network without a product at rank 2 settle at the
#' rank asked for with the five-node rule, and a count network of 40 actors at
#' rank 3 settles by retiring a column (rank 2). The data determine none of the
#' six products (0.00 to 0.39 units of signal per loading), and every one is
#' withheld.
#'
#' \strong{Predictions.} The response of a normal fit is its linear predictor plus
#' Gaussian noise, so \code{fitted()} and \code{predict()} return the posterior mean of the
#' predictor (\code{fit$fitted}), which is the mean of the response; the rest of this paragraph is about
#' binary and count fits. \code{fitted()} and \code{predict()} of a binary
#' empirical-Bayes fit return predictive probabilities: the link integrated over
#' the posterior of the linear predictor rather than applied at its mean. The
#' predictor of a cell is uncertain, so the plug-in probability of a held-out
#' event can be orders of magnitude too small. A Poisson fit returns the plug-in
#' mean count \eqn{e^{\eta + v_e / 2}}{exp(eta + ve / 2)} by default, as the
#' samplers do: the posterior predictive mean of a count is not finite in this
#' model (the inverse-gamma tail of \eqn{v_e}{ve}, and the inverse-Wishart tail of
#' the scale of the factors for a cell the data do not inform), so the value of the
#' integrated mean \eqn{e^{m + (v + v_e) / 2}}{exp(m + (v + ve) / 2)} would be a
#' property of the quadrature and of the prior; it stays available as
#' \code{plug_in = FALSE}, and the predictive law of a held-out count (its
#' quantiles and event probabilities, which are finite) is what the held-out
#' scores below use.
#' On the 170-state panel of Chapter 11 (10 periods, 15 percent of the cells
#' held out, 39 and 45 events on the two masks) the log loss of the held-out
#' cells at \code{R = 1} to 3 was 0.00427, 0.00440 and 0.00450 (first mask)
#' and 0.00506 at each rank (second mask, where every fit at \code{R = 1} to 3
#' retires its columns and returns the rank-0 fit) under the plug-in values,
#' and 0.003930, 0.003934 and 0.003972, and 0.004568, under the predictive ones
#' (sampler, one chain, from its draws: 0.003902, 0.003886 and 0.003897, and
#' 0.004480, 0.004448 and 0.004449); the smallest probability of a held-out
#' event was 7e-7 to 3e-6 (first mask) and 7e-9 under the plug-in values and
#' 4.5e-5 to 5.9e-5 and 5e-6 under the predictive ones. At \code{R = 0} the
#' predictive log loss, 0.004054 and 0.004568 on the two masks, is the
#' sampler's (0.004049 and 0.004585). The sampler's own predictor, replaced by
#' a Gaussian with its mean and variance, loses 0.003963, 0.003942 and 0.003947
#' on the first mask against 0.003902, 0.003886 and 0.003897 from the draws, so
#' at ranks 1 and 2 the gap of the fast fit to the sampler is in the shape of
#' the sampler's posterior of the predictor, not in its moments (at rank 3 the
#' fast fit is 0.00003 above the Gaussian replacement).
#' The predictor is \eqn{N(m, v)}: \eqn{m} is the fitted predictor, the
#' posterior mean; \eqn{v} is the variance \eqn{g' S g} of the cell's additive
#' predictor (the coefficients and the two effects, with \eqn{g} its gradient
#' and \eqn{S} the inverse of their penalised information at the expected
#' weights) plus the variance of the product. The probability is
#' \eqn{\Phi(m / \sqrt{1 + v})} for a probit link and
#' \eqn{\mathrm{plogis}(m / \sqrt{1 + \pi v / 8})} for a logit link.
#' \code{plug_in = TRUE} returns the link applied to the fitted predictor (what
#' \code{fit$fitted} holds and \code{simulate()} uses). Forming the covariance
#' costs one dense inverse, of the size of the one behind \code{vcov()}, per
#' call.
#'
#' \strong{Choosing the rank.} \code{fit$eb$log_marginal} is the Laplace
#' approximation to the log marginal likelihood of the data with the effects
#' and factors integrated, at the posterior mode that the variational stage
#' starts from (for a normal fit the mode at the fitted variance of the noise,
#' \strong{Choosing the rank of a normal fit}), with that mode's coefficients and scales: the deviance and the
#' penalties of the effects and the factors, with the log determinant of the
#' information of all of them, on the scale of the likelihood (an undirected
#' network counts each dyad once; the hyperpriors of the scales are left out).
#' The deviance falls with the rank; this does not, so the difference between
#' two ranks fitted to the same data is an approximate log Bayes factor. The
#' sign and label multiplicity of the factors is not counted, which favours the
#' lower rank. On 1,040 simulated binary networks (directed, bipartite and
#' undirected; generating rank 0 to 3, fitted ranks 0 to 3; 60 and 100 actors,
#' one and four periods, 7 to 25 percent ties and two sparser designs) 17 of the
#' 4,160 fits
#' were refused as not converged (class \code{"lame_als_convergence"}), all in
#' the two sparsest designs without a latent term; of the 1,033 networks fitted
#' at all four ranks, the rank with the largest value was the generating rank in
#' all of the 309 with at least 10 training events per factor loading, in 96
#' percent of the 386 with 3 to 10, in 59 percent of the 234 with 1 to 3 (35
#' percent were too low) and in 38 percent of the 104 with fewer than 1 (27
#' percent were too high, 36 percent too low); the lowest rank within one unit
#' of the largest value was right in 100, 96, 62 and 43 percent. With little
#' information per loading the data do not determine the rank. Take the lowest
#' rank whose value is within one unit of the largest: it named the generating
#' rank more often than the largest value below 3 events per loading (62
#' against 59 and 43 against 38 percent) and as often from 3 to 10 (96
#' percent). The variational lower bound,
#' \code{fit$eb$variational$elbo}, is not the value to use for this. A Gaussian
#' posterior cannot hold the broad posterior of positions that the data barely
#' identify, and the bound's slack grows with the rank: on the 1,033 networks its
#' largest value was the generating rank in 35, 56, 96, 100 and 100 percent
#' (fewer than 1, 1 to 3, 3 to 10, 10 to 30 and more than 30 events per
#' loading), never too high but too low for 96 percent of the networks with a
#' latent term and fewer than 1 event per loading, where the Laplace value found
#' the rank in 38 percent (and was too high for 64 percent of the networks without
#' a latent term and fewer than 1 event per loading). On the 170-state panel of
#' Chapter 11 the Laplace value was -1137.8, -1119.3, -1117.2 and -1121.1 at
#' \code{R = 0} to 3 on the first mask, so ranks 1 to 3 are preferred to 0 (by
#' 17 to 21 units, most at rank 2), as the sampler's held-out log loss prefers
#' them, while the lower bound falls with the rank (-1144.7, -1152.9, -1170.6
#' and -1195.9). On the second mask every fit at \code{R = 1} to 3 retires its
#' columns and returns the rank-0 fit (Laplace value -1117.3, lower bound
#' -1125.1), so neither value can prefer a latent term there, although the
#' sampler's held-out log loss does. The value is the evidence of the independent-cell
#' model at the posterior mode: a fit with a dyadic correlation (\strong{Dyadic
#' correlation}) reports the same \code{log_marginal} as the fit with
#' \code{dcor = FALSE}. On 15 simulated networks of 60 actors with a latent
#' correlation of 0.6 (generating ranks 0 to 2) the rule above picked the
#' generating rank in 10 and a higher one in none.
#'
#' What picks the rank where the Laplace value cannot is the held-out predictive
#' log score: mask a share of the observed cells (15 percent), refit at each
#' rank, and score the masked cells with \code{fitted()} (the integrated
#' predictive; a variational fit's predictive is not shrunk to the plug-in
#' probability, and its surplus columns are not shrunk away, so the score falls
#' above the generating rank). On the 1,033 networks it found the generating rank
#' in 51, 72, 85, 91 and 92 percent (the Laplace value in 38, 59, 96, 100 and
#' 100), and where there are fewer than 3 events per loading it was right for 64
#' percent of the networks without a latent term (the Laplace value 72, with the
#' one-unit rule 87; the lower bound 100) and 66 percent of those with one
#' (Laplace 45, the lower bound 32): no other criterion that needs no refit was
#' right for both. The rule that uses the held-out score below 3 events per
#' loading (\code{fit$eb$variational$information["per_loading"]} of the fit at the
#' rank the one-unit rule names, rank 1 where it names 0) and the Laplace value
#' with the one-unit rule above found the generating rank in 85 percent of the
#' 1,033 networks (the Laplace value alone 83 and 84, the lower bound 82, the
#' held-out score alone 81). On the 170-state panel the held-out log loss at
#' \code{R = 0} to 3 was 0.004054, 0.003930, 0.003934 and 0.003972 on the first
#' mask and 0.004568 at every rank on the second, whose fits at \code{R = 1} to 3
#' retire to rank 0 (the sampler's, one chain: 0.004049, 0.003902, 0.003886 and
#' 0.003897, and 0.004585, 0.004480, 0.004448 and 0.004449; the sampler's loss is
#' lowest at rank 2 on both masks, with rank 3 within 0.00001 of it).
#'
#' The variational stage needs a median of 46 iterations (87 with fewer than 1
#' event per loading), and 9 in 10 fits finish within 130. The
#' 170-state panel of Chapter 11 takes 29 seconds at \code{R = 0} and 1.5 to 3
#' minutes at \code{R = 1} to 3 (one start, one thread; 8 minutes at \code{R = 2}
#' with \code{multistart = "full"}). The
#' joint solve is dense, so the cost grows with the
#' cube of the number of actors.
#' An actor without a tie (or without a
#' non-tie) on a side keeps a finite, shrunken effect and stays in the fit, so
#' there is no \code{fit$not_estimated}. The effects \code{a} and \code{b} are
#' not re-centred: the intercept is the population intercept (the level of a
#' dyad between two actors with zero effects, at the means of the node
#' covariates) and the effects have mean zero a priori; without node covariates
#' the posterior-mean effects of the variational fit sum to zero, because the
#' intercept is stationary. \code{va}, \code{cab}
#' and \code{vb} are the fitted covariance of the effects, the EM estimate of
#' \eqn{\Sigma_{ab}}{Sigma_ab}, not descriptive summaries of the estimated effects.
#'
#' \strong{Accuracy.} Against the truth, on the 1,040 simulated networks above
#' (dyadic slope 0.7; sender and receiver sd 0.6, correlation 0.3; variance 0.6 a
#' latent dimension), fitted at the generating rank (above zero), the relative
#' bias of the slope was, by training events per loading (fewer than 1, 1 to 3,
#' 3 to 10, 10 to 30, more than 30), -3.9, -2.4, -0.8, -0.1 and 0.4 percent, with
#' root mean square error 10.4, 7.1, 4.8, 2.8 and 1.8 percent (the posterior
#' mode: -0.7, 0.8, 1.1, 0.7 and 0.8; 9.8, 7.2, 4.9, 2.9 and 1.9). One rank above
#' the generating one the bias was -0.5, 1.4, 2.2, 1.3 and 1.0 percent (mode 0.5,
#' 5.1, 5.2, 2.6 and 1.8), and with no latent term in the data and a latent rank
#' fitted, with fewer than 3, 3 to 10, 10 to 30 and more than 30 events per
#' actor, it was 17.5, 5.7, 4.0 and 1.1 percent (mode 12.4, 8.5, 6.7 and 2.1). At
#' the generating rank the fit withheld the product in 85, 70 and 16 percent of
#' the fits with fewer than 1, 1 to 3 and 3 to 10 events per loading and in none
#' above; the products it kept correlated 0.57, 0.78, 0.87, 0.93 and 0.96 with
#' the generating one and had norm 43, 54, 57, 58 and 52 against the truth's 81,
#' 72, 67, 63 and 54 (the mode, every fit: correlation 0.48, 0.68, 0.85, 0.93 and
#' 0.96, norm 49, 60, 61, 62 and 54 against 69, 62, 63, 63 and 54), and the
#' held-out log loss was no higher than the mode's in any band. Against the
#' sampler's posterior mean, on 101 fits with fewer than 1 training event per
#' loading and 96 with 1 to 3 (directed, bipartite and undirected networks of 60
#' and 100 actors, ranks 1 to 3, one \code{lame()} run of two chains each), the
#' slope of the fit lay on average 0.09 (95 percent interval -0.001 to 0.18) and
#' 0.10 (0.05 to 0.16) sampler standard deviations from it (relative to the
#' slope 1.1 and 1.2 percent), and the standard error was 1.05 and 1.03 times the
#' sampler's standard deviation. The fit withheld the product in 90 and 88 of
#' these fits; of the 135 directed fits, 132 withheld it, two with fewer than 1
#' event per loading retired every column and one kept it. The 9 and 8 products
#' it kept (bipartite and undirected, and the one directed) had 0.74 to 1.08 and
#' 0.92 to 1.07 times the norm of the sampler's posterior mean product. For
#' undirected networks the slope lay 0.40 and 0.34 standard deviations above the
#' sampler's.
#' Where the sampler finds no structure the variational mean keeps a sign pattern
#' from the mode: in 4 of 24 networks of 60 actors its norm is 2.3 to 7.8 times
#' the sampler's posterior-mean norm, which is the small remainder (0.11 to 0.37)
#' of a cloud of mutually near-orthogonal, sign-symmetric draws of about the size
#' of the fitted product (mean cosine between draws 0.004 to 0.11); a unimodal
#' Gaussian posterior cannot report that mean, and the fit withholds these four
#' products (0.30 to 1.18 units of signal per loading). The slope on those four
#' (\eqn{|z| \le 0.83}{|z| <= 0.83}), its standard error (1.08 to 1.10 times the
#' sampler's standard deviation) and the held-out log loss (0.011 below to 0.002
#' above the sampler's) are within the usual bands. For bipartite
#' networks of 60 actors fitted at \code{R = 3} (5 networks, 0.65 to 1.22 events
#' per loading) the slope lies within 1.52 sampler standard deviations of the
#' sampler's (\eqn{|z| \le 0.41}{|z| <= 0.41} on four) and the fit withholds the
#' product, whose variational mean has 0.76 to 1.21 times the norm of the
#' sampler's posterior mean product. The
#' scales of the sender columns are the eigenvalues of a Wishart matrix, as in
#' the sampler, so that a column the data do not inform falls to the small
#' eigenvalues: on these 5 networks the variance of the predictor of the held-out
#' cells is 0.57 to 1.02 times the sampler's, and the held-out log loss 0.003
#' below to 0.014 above its. In
#' cells of 40 networks (directed 30 to 100 actors, bipartite and undirected 60
#' actors, one and four periods, ranks 0 to 2, 14 cells) the slope bias was -2.2
#' to 1.2 percent (mode -0.5 to 3.4), the 95 percent interval from
#' \code{vcov()} covered the slope in 80 to 100 percent of the networks (mean 95)
#' and the intercept in 88 to 98 (mean 92), and the standard error was 0.78 to
#' 1.43 times the standard deviation of the slope. On the 170-state panel the
#' coefficients of the full-system fit (\code{R = 2}) lay within 0.31 standard
#' deviations of the sampler's posterior mean, with standard errors 0.86 to 1.23
#' times the sampler's; the fit withholds the product (0.26 units of signal per
#' loading, below 2.1), whose variational mean has norm 14.4 against the 15.9 of
#' the sampler's posterior mean product (correlation 0.93). On the
#' collinear-covariate design (40 actors, 8 periods, \code{R = 2}, truth 0.5, 24
#' networks) the covariate that follows the latent product was 28.0 percent too
#' large (sampler 26.2, mode 13.5, no latent rank 93): the prior shrinks the
#' product and the covariate absorbs what it shrinks.
#' The independent-cell fit (\code{dcor = FALSE}) has no dyadic correlation. On a
#' strongly reciprocal panel (the school network of Chapter 9, 32 actors,
#' \code{R = 2}; \code{ame()} with the book's specification, two chains of 50,000
#' scans, puts \code{rho} at 0.67, posterior standard deviation 0.16) its
#' coefficients lie within 0.72 sampler standard deviations of the sampler's,
#' with standard errors 1.07 to 1.14 times the sampler's, and its post hoc
#' correlation is 0.98, 2.0 posterior standard deviations above the sampler's.
#' The pair model (\strong{Dyadic correlation}), the default for a directed
#' network, estimates \code{rho} at 0.74 (0.43 posterior standard deviations
#' from the sampler's), and its coefficients lie within 0.74 sampler standard
#' deviations, with standard errors 0.76 to 0.97 times the sampler's. Both fits
#' withhold the product of this panel. Choose the rank as described above.
#'
#' Without the empirical-Bayes prior (\code{bias_correct = "none"} and the two
#' explicit corrections below) a binary fit is the maximum-likelihood fixed-effects probit
#' (or logit)
#' fit, penalised at \code{R > 0}: each reweighting step caps the
#' working-response increment at 30 and gives a capped cell the weight at
#' which it still contributes its score, and the loop runs to its fixed
#' point (on 32 simulated networks of 60 and 120 actors the coefficient and
#' deviance matched \code{glm()} to 5e-6 and 3e-8). In a fixed-effects fit an
#' actor with a single
#' outcome value on a side (a sender without a tie, a receiver with only
#' ties, and so on) has its effect there at the boundary of the likelihood
#' and its cells on that side carry no information about the other
#' parameters: those cells are left out of the fit (and of \code{fit$Y}),
#' the effect is \code{NA} in \code{a} / \code{b} with its linear predictor
#' and fitted values, the actors are listed in \code{fit$not_estimated} with
#' a message, and \code{va} / \code{vb} describe the estimated effects. On
#' simulated 250-actor networks with 1 to 50 such senders the coefficient,
#' \code{va} and \code{vb} stayed where they were without them.
#'
#' \strong{Split-sample bias correction.} A binary IRLS fit is a
#' fixed-effects probit (or logit) fit: every additive effect and every
#' latent position is estimated from the cells of its own actor, and
#' maximising over all of them stretches the latent scale outward in finite
#' samples (the incidental-parameters problem), while the factor penalty
#' pulls the product term in. On simulated 30 x 20 bipartite panels with five
#' periods and \code{R = 2} (20 seeds), the uncorrected dyadic coefficient
#' was 3.3\% too large, the sender effects 5.8\% too spread, \code{va}
#' 15\% and \code{vb} 11\% too large, and the fitted product surface 13\%
#' too flat; at 60 x 40 the stretch was 2 to 5\%, and on directed 30-actor
#' panels 2 to 6\%. \code{bias_correct = "jackknife"} removes it by a
#' split-sample jackknife wherever the jackknife holds (see the size boundary
#' below; the empirical-Bayes default needs no such correction). The
#' observed cells are split into two interleaved halves, the model is refit
#' on each half, warm-started from the full fit with the same settings, and
#' every identified quantity is reported as \eqn{2 \times} full minus the
#' mean of the two halves: a half carries twice the stretch of the full
#' panel, so the combination cancels it to first order. A half that has not
#' converged within \code{max_iter} is continued from its last iterate with
#' five times that limit; if it still has not converged the call is an error of
#' class \code{"lame_jk_nonconvergence"} (the class of the dynamic count
#' jackknife), since the correction never combines a last iterate. The split is the
#' checkerboard of the cells with \eqn{i + j + t} even and odd (every dyad
#' keeps the periods of one parity in each half, every actor about half of
#' its cells); when it leaves some actor without a tie or without a non-tie
#' in a half, a random split of the dyad-periods (drawn from a fixed internal
#' seed; the caller's random number stream is left unchanged) and then the
#' split into odd and even periods are tried, the latter with weights for
#' its unequal halves. Corrected are the intercept and regression
#' coefficients, the additive effects \code{a} and \code{b}, the
#' double-centred product surface (\code{U}, \code{V}, \code{L} are the
#' rank-\code{R} factorisation of the corrected surface) and \code{va},
#' \code{vb}, \code{cab}; \code{EZ}, the fitted values and the residuals are
#' rebuilt from them; \code{ve} and the IRLS weights are those of the full
#' fit, and \code{rho} is re-estimated from the paired binary likelihood at
#' the corrected marginal probabilities. On the panels above the corrected coefficient ratio
#' was 0.990 (Monte Carlo standard error 0.010) at 30 x 20 and 1.002 at
#' 60 x 40, the effect spreads within 1.5\%, \code{va} and \code{vb} within
#' 2\% at 30 x 20 and 60 x 40 and within 3\% on directed panels, and the
#' product slope 0.97 to 1.04. The correction costs two half-data refits,
#' 1.6 to 3.4 times the time of the uncorrected fit, and widens the spread
#' across samples slightly (standard deviation of the coefficient 0.043 to
#' 0.047 at 30 x 20). \code{fit$raw} keeps the uncorrected coefficients,
#' \code{a}, \code{b}, \code{VC}, product surface \code{O} and deviance
#' (\code{fit$deviance} is that of the corrected \code{EZ});
#' \code{fit$bias_correct} records whether the correction was applied, the
#' split used and the splits tried, the half fits' convergence flags and
#' iteration counts, and the shrink factors of the effects and the product.
#' \code{bias_correct = "none"} returns the uncorrected fit and
#' \code{bias_correct = "jackknife"} requests the correction for another
#' family.
#'
#' A single cross-section at \code{R = 0} (\code{ame_als()}) is split by the
#' checkerboard of its cells (\eqn{i + j} even and odd, each actor keeping
#' half of its cells in each half). Each half's observed graph then falls
#' into two sender-receiver components whose relative level is not
#' identified, so every cell's half value is taken from the half that
#' observes that cell's block. On simulated directed networks (20 to 40
#' seeds) the corrected slope, \code{va} and \code{vb} were 0.990, 0.982 and
#' 0.954 of the truth at 75 actors with 36\% ties, 0.990, 0.989 and 0.965 at
#' 100 and 0.999, 0.979 and 0.996 at 150 actors with 24\% ties (uncorrected
#' 1.017 to 1.029 and up to 1.09). At \code{R > 0} a single cross-section is
#' not corrected: the checkerboard over-corrected wherever the uncorrected
#' estimate needs it (100 actors, \code{R = 2}: slope 1.024 to 0.971,
#' \code{vb} 1.090 to 0.950; 50 actors: \code{va} 1.209 to 0.791), and so did
#' a random split. Such a cross-section is corrected by the parametric
#' bootstrap instead (see \strong{Parametric-bootstrap bias correction}).
#'
#' \strong{Parametric-bootstrap bias correction.} \code{bias_correct =
#' "bootstrap"} simulates the bias instead of splitting the sample, which
#' reaches the panels the split-sample correction does not (a single
#' cross-section with \code{R > 0} or below 25 rarer-outcome cells per actor,
#' a panel of fewer than five periods at \code{R > 0}, an actor without both
#' outcomes in a half). Networks are drawn from the fitted model with fresh
#' additive effects from the fitted \code{va}, \code{vb}, \code{cab} and the
#' intercept, coefficients, node-covariate parts, product surface and
#' reciprocity held at their fitted values, and are fitted by the same
#' uncorrected estimator, warm-started from the fit. The bias of each
#' coefficient, of \code{va}, \code{vb}, \code{cab} and of each cell of the
#' product surface is the mean of its replicate values less the value it was
#' generated with (the variance and covariance of the drawn effects for
#' \code{va}, \code{vb}, \code{cab}); the additive effects are divided by the
#' mean slope of the replicate effects on the drawn ones. The bias grows with
#' the spread of the effects, so it is larger at the stretched estimate than
#' at the truth, and one subtraction over-corrects small networks (30 actors,
#' \code{R = 0}, 150 networks: slope -2.4\%, \code{vb} 0.91 of the truth;
#' subtracting the bias at the fitted effects themselves gave 0.85); a second
#' round of replicates at the once-corrected parameters, drawn with the same
#' random numbers, gives the slope \eqn{s} of the bias, and the estimate is
#' corrected by \eqn{b_0 / (1 + s)}, the secant step toward the parameters
#' whose simulated estimate equals the fitted one. Each round has 30
#' replicates from a fixed internal seed (the fit is reproducible and the
#' caller's random number stream is left unchanged), so the correction costs
#' about 60 refits. On simulated directed probit networks (mean -1, about 27\%
#' ties, dyadic slope 0.7, sender and receiver sd 0.6, factor variance 0.6 at
#' \code{R = 1}; 100 to 150 networks a cell) the relative bias of the slope
#' went from 12.2\% to -0.4\% at 30 actors and \code{R = 0}, 17.4\% to 0.0\%
#' at 30 and \code{R = 1}, 6.8\% to 0.6\% at 50 and \code{R = 0}, 8.4\% to
#' 0.2\% at 50 and \code{R = 1}, 4.0\% to 0.3\% at 100 and \code{R = 1},
#' 1.9\% to 0.1\% at 200 and \code{R = 1}, 3.4\% to 0.5\% on 30-actor
#' panels of four periods at \code{R = 1}, and on 40 x 30 bipartite
#' cross-sections (60 networks) 9.3\% to 0.1\% at \code{R = 0} and 15.8\% to
#' 2.4\% (Monte Carlo standard error 1.4\%) at \code{R = 1};
#' \code{va} and \code{vb} went from 1.10 to 1.51 of the truth to 0.94 to
#' 1.08 (Monte Carlo standard error 0.02 to 0.05), the intercept bias from
#' -0.03 to -0.11 to within 0.02, and the slope of the product surface on the
#' generating one from 0.71 to 0.79 (30 actors), 0.86 to 0.93 (50) and 0.94
#' to 0.98 (100); the rest of its flattening is the factor penalty. The 95\%
#' interval estimate plus or minus 1.96 replicate standard deviations of the
#' first round (\code{fit$bias_correct$replicate_sd}) covered the slope in
#' 93 to 99\% of the networks in these seven cells (88\% in the 69
#' bootstrap-corrected fits among 100-actor cross-sections at \code{R = 0}),
#' against 72 to 92\% with the same standard deviation around the
#' uncorrected estimate and 49 to 82\% with the conditional sandwich
#' (\code{vcov()}). The correction took 4 to 8 seconds at \code{R = 0}
#' and 40 seconds (30 actors) to 9 minutes (200 actors) at \code{R = 1}
#' with twelve such fits running at once, against 0.2 to 9 seconds
#' uncorrected. \code{fit$bias_correct} records \code{method = "bootstrap"},
#' the replicates of each round that refit and converged, the first-round
#' biases, the secant factors and the corrections applied; \code{fit$raw}
#' keeps the uncorrected values. A refit of the corrected fit
#' (\code{\link{ame_als_refit}}, and every replicate of
#' \code{\link{ame_als_bootstrap}}) removes these stored corrections rather
#' than simulating again. Fewer than 80\% of a round's replicates refitting
#' and converging is an error of class \code{"lame_als_size"} (1 of 100
#' networks at 30 actors and \code{R = 1}), and so is an undirected fit with
#' \code{R > 0}, whose fitted product is flattened beyond what the replicates
#' reproduce: on undirected 50-actor networks with a rank-1 term of variance
#' 1 (60 networks) the product slope was 0.68 and the correction moved the
#' dyadic slope from 1.0\% above to 5.8\% below the truth, and with a term of
#' variance 0.36 (50 networks) the factor penalty removed the product in
#' every fit and the slope went from 6.5\% to 12.6\% below the truth.
#' Undirected fits at \code{R = 0} are corrected (50 actors: slope 6.6\% to
#' -0.8\%, \code{va} 1.12 to 0.98 of the truth).
#'
#' \strong{Size boundary.} A fixed-effects binary fit is supported where the
#' bias of what it reports is within a tolerance, 2\% on the dyadic coefficients and 5\% on
#' \code{va} and \code{vb} (the internal constant \code{.ae_bias_tol}; the
#' floors below follow from it). The size that governs the bias is the
#' number of observed cells per actor (over all periods, on the side with
#' fewer: senders or receivers, rows or columns) times the share of the
#' rarer outcome (ties in a sparse network), called rarer-outcome cells per
#' actor below, and the bias also grows with the spread of the additive
#' effects, measured on the fit as the larger of the sd of the estimated
#' sender and receiver effects over the latent noise sd (1 for probit,
#' \eqn{\pi / \sqrt{3}}{pi / sqrt(3)} for logit). On simulated networks (single
#' cross-sections of 60 to 500 actors at 7 to 40\% ties, panels of two to
#' four periods, directed, bipartite and undirected, probit and logit,
#' additive-effect sd 0.5 to 2, 3 to 10 seeds a cell) the uncorrected bias
#' fell as \eqn{c /} size, with \eqn{c} growing with the rank and the
#' spread. The split-sample correction is used where it can be computed:
#' panels with two or more periods at \code{R = 0} and five or more at
#' \code{R > 0}, and single cross-sections at \code{R = 0}, in each case
#' with every actor keeping a tie and a non-tie in both halves and, at
#' \code{R = 0}, from 25 rarer-outcome cells per actor (below that it
#' over-corrected \code{va} and \code{vb} to 0.75 and 0.87 of the truth at
#' 14, and to 0.86 and 0.89 at 18). The empirical-Bayes default does not use
#' this rule, which governs the fixed-effects fits of the families and settings
#' the estimator does not cover: a fit the split-sample correction does not
#' reach is reported uncorrected when its size is at or above the floor for
#' its rank and spread, and refused below it (the parametric bootstrap does not
#' reach an undirected fit with a multiplicative term); a
#' fit whose spread is above 2.25, or whose rank is above 8, is refused,
#' corrected or not, since neither the uncorrected nor the corrected estimate
#' was measured there. The floors (rarer-outcome cells per actor; in
#' brackets those at 3\% / 8\%), where every measured cell was within the
#' tolerance up to one Monte Carlo standard error:
#' \tabular{llll}{
#'   rank \tab spread up to 0.7 \tab 0.7 to 1.25 \tab 1.25 to 2.25 \cr
#'   \code{R = 0} \tab 73 (46) \tab 73 (46) \tab 95 (60) \cr
#'   \code{R = 1, 2} \tab 78 (50) \tab 94 (63) \tab 155 (98) \cr
#'   \code{R = 3, 4} \tab 117 (74) \tab 124 (83) \tab refused (155) \cr
#'   \code{R = 5} to \code{8} \tab 117 (74) \tab 145 (97) \tab refused (refused)
#' }
#' "Refused" marks a band where no measured size (up to about 175 cells per
#' actor) met the tolerance; ranks above 8 are not measured and refused. For
#' a directed cross-section with 36\% ties and spread up to 0.7 the floors
#' are about 204, 218 and 326 actors at \code{R = 0}, \code{1-2} and
#' \code{3-8}, and half that with two periods. A fit above the floor
#' records in \code{fit$bias_correct} \code{method = "none"} with the
#' \code{size}, the \code{floor}, the \code{spread} and a \code{basis}
#' giving the measured bias band there, and \code{print()} /
#' \code{summary()} say that no correction is needed; there is no warning.
#' Above the measured spread or rank, binary ALS stops with an error of class
#' \code{"lame_als_size"} that names the spread or rank and the measured
#' band; the MCMC \code{ame()} / \code{lame()} integrate the effects and
#' positions out and are the calibrated route there. The parametric bootstrap
#' also stops with that class when fewer than 80\% of its simulated
#' replicates refit and converge (networks simulated from the fit too sparse
#' for the fixed-effects estimator). \code{bias_correct = "none"} returns the
#' uncorrected fit at any size, recorded as requested
#' (\code{fit$bias_correct$requested}); \code{bias_correct = "jackknife"}
#' requires the correction and is an error (class \code{"lame_jk_split"})
#' where it cannot be computed. A bootstrap of a binary fit that carries no
#' correction is an error (class \code{"lame_als_bootstrap_bias"}), whatever its
#' size: each replicate is refit at the biased estimate and reproduces the
#' bias, which falls more slowly than the standard error
#' (\code{\link{ame_als_bootstrap}}). The default empirical-Bayes fit and the fits
#' corrected by the jackknife or the parametric bootstrap have one.
#'
#' Rank and ordinal families
#' (\code{"ordinal"}, \code{"cbin"},
#' \code{"frn"}) are not supported by ALS and raise an informative error; use
#' \code{\link{ame}} or \code{\link{lame}} for those likelihoods.
#'
#' \strong{Row/column (node) covariates.} A node covariate broadcasts to a
#' per-actor constant, which is collinear with the additive sender/receiver
#' effect, so its coefficient is not identified by the objective alone. It is
#' identified here by an explicit constraint: the additive effects are taken
#' orthogonal to the node covariates, and \code{beta_row}/\code{beta_col} are
#' the corresponding between-actor regression coefficients (the additive
#' effects then carry only the residual heterogeneity). This is the standard
#' estimand under the assumption that the additive effects are uncorrelated
#' with the node covariates; if that assumption is doubtful the coefficient
#' absorbs the covariate-correlated part of the additive heterogeneity.
#' Time-varying node covariates are summarised by their per-actor mean; if a
#' node covariate varies within actor over time the discarded within-actor
#' variation triggers a warning.
#'
#' \strong{Missing covariates.} As in \code{\link{ame}} and \code{\link{lame}},
#' a cell whose dyadic covariate, sender covariate or receiver covariate is
#' missing in its period is unobserved: it is set to \code{NA} in \code{Y}
#' before the fit, a message gives the number of cells, and the coefficients
#' are estimated from the complete cells. A missing value is not modelled as
#' 0, which biases the coefficient of a covariate that is not centred. A
#' missing value on the diagonal of a square unipartite covariate is not
#' counted, since the diagonal is unobserved anyway. In an undirected fit
#' (\code{symmetric = TRUE}) a dyad has one covariate value: a hole in one
#' triangle of a dyadic covariate is filled from the other triangle, and a
#' dyad missing in both is unobserved. The per-actor mean of a node covariate
#' runs over the periods in which it is observed. An actor missing a node
#' covariate in every period has no
#' observed cell on that side: the empirical-Bayes fit leaves its effect there
#' to the prior (zero when the actor has no cell on either side), and a
#' fixed-effects fit reports it as not estimated (\code{NA} in \code{a} /
#' \code{b}, listed in \code{fit$not_estimated} as \code{"no observed
#' cell"}). A covariate that is missing everywhere stops the fit, as it does
#' in \code{\link{ame}} and \code{\link{lame}}.
#'
#' \strong{Identifiability of the additive and multiplicative terms.} For
#' \code{R > 0} the additive term \eqn{a_i + b_j} and the multiplicative term
#' \eqn{u_i' v_j} are not separately identified by the objective alone: a
#' broadcast (row- or column-constant) component can sit in either, since a
#' pure sender effect \eqn{a 1'} is itself rank one. The estimator imposes the
#' standard AME gauge --- the multiplicative term is double-centered (zero row
#' and column means), so all broadcast structure is carried by \eqn{a, b} ---
#' which makes the reported \eqn{a, b, U, V} unique given the fitted values and
#' independent of the optimisation path. For a unipartite network the
#' double-centering means are taken over the full matrix, which includes the
#' structurally unobserved self-tie diagonal that the model fills in by its
#' low-rank completion; the additive/multiplicative split therefore carries an
#' \eqn{O(1/n)} dependence on that completion. On a disconnected observed-dyad
#' graph the additive effects additionally have a per-component level shift,
#' which is pinned to a precision-weighted minimum-norm gauge so that
#' \code{a}, \code{b} and the variance components remain reproducible. A dyadic
#' covariate that is itself (near) low-rank can still be partially aliased with
#' the multiplicative term, so keep \code{R} modest relative to the covariate
#' structure. This residual aliasing is intrinsic to the AME model --- the MCMC
#' estimator resolves it only through its priors.
#'
#' \strong{Choosing R.} An empirical-Bayes fit (the default for normal, binary and
#' Poisson responses) has a rank criterion, \code{fit$eb$log_marginal} (see \strong{Choosing the
#' rank}). For the fixed-effects fits there is no
#' automatic order-selection criterion (the
#' working-response objective has no likelihood, so AIC/BIC do not apply). Fit a
#' few values --- e.g. \code{lapply(0:4, function(r) ame_als(Y, R = r,
#' ...))} --- and inspect \code{deviance} (the residual sum of squares for the
#' \code{normal} family and the transform route, the binomial or Poisson
#' deviance of the fitted values for an IRLS fit, in each case of the
#' reported \code{EZ}; the optimiser's final objective, which adds the
#' factor penalty of an IRLS fit, is \code{objective}): it
#' falls steeply while real multiplicative signal is being captured and then
#' flattens, so the \dQuote{elbow} of that curve is a reasonable choice. With
#' \code{R > 0} the objective is non-convex; use \code{multistart} to guard
#' against local optima.
#'
#' A completely separating unpenalized linear predictor has no finite estimate;
#' the IRLS route refuses it with class \code{lame_als_separation}. Use
#' \code{ame()} or \code{lame()} with proper priors for such data.
#'
#' \strong{Uncertainty.} An empirical-Bayes fit (the default for normal, binary
#' and Poisson responses) reports Laplace standard errors for every coefficient
#' (\code{vcov()}, \code{confint()}, \code{fit$eb$se}). For a fixed-effects fit
#' (\code{bias_correct = "none"}) uncertainty is obtained separately, by the
#' bootstrap; see \code{\link{ame_als_bootstrap}}. (The SIR paper's own primary standard
#' errors are Hessian-based, classical and sandwich/robust. For AME a
#' Hessian-based variance is awkward on two counts: without an explicit gauge
#' fix the rotational invariance of the multiplicative factors leaves the joint
#' Hessian rank-deficient, and even with a gauge fixed the Gaussian
#' working-response approximation used for the non-normal families leaves a
#' Hessian-based variance miscalibrated. The bootstrap side-steps both, so it is
#' preferred here.)
#'
#' @param Y a square (unipartite) or rectangular (bipartite) relational matrix.
#'   A cross-sectional \pkg{netify} object is also accepted and converted with
#'   \code{netify::to_lame()}; when \code{family} or \code{mode} is omitted,
#'   the value inferred by \pkg{netify} is used.
#' @param Xdyad an \code{n_row x n_col} matrix or \code{n_row x n_col x pd}
#'   array of dyadic covariates, or \code{NULL}. A 2-D matrix is wrapped as a
#'   single slice named after the argument expression, so
#'   \code{Xdyad = distance} yields the coefficient \code{distance_dyad}, the
#'   same name \code{\link{ame}} gives it.
#' @param Xrow an \code{n_row x pr} matrix of row/sender covariates, or \code{NULL}.
#' @param Xcol an \code{n_col x pc} matrix of column/receiver covariates, or \code{NULL}.
#' @param R integer dimension of the multiplicative effects (default \code{0}).
#'   \strong{The covariate coefficients are conditional on this choice.} The
#'   multiplicative term \eqn{u_i'v_j} is a flexible high-variance regressor
#'   that can correlate with the dyadic covariates, so the estimated
#'   \code{beta} can shift -- and occasionally change sign -- as \code{R}
#'   increases. Comparing against an \code{R = 0} fit is a useful check on
#'   whether the covariate story is being driven by the latent rank.
#' @param family one of \code{"normal"}, \code{"binary"}, or
#'   \code{"poisson"}. The rank and censoring families are MCMC-only.
#' @param mode \code{"unipartite"} (square) or \code{"bipartite"} (rectangular).
#' @param symmetric logical; fit a symmetric (undirected) model. Unipartite only.
#'   Every dyadic covariate slice must be symmetric: an asymmetric one is an
#'   error of class \code{lame_symmetric_xdyad} that names the slice (and, for a
#'   panel, its first period), as in \code{\link{ame}} and \code{\link{lame}}.
#' @param max_iter maximum number of block coordinate descent iterations
#'   (default 200), and maximum number of outer reweighting steps for IRLS. A
#'   fixed-effects fit that has not converged within it is
#'   refitted once with four times the limit, and is an error (class
#'   \code{"lame_als_convergence"}) if that run has not converged either, as is
#'   an empirical-Bayes fit that does not settle (without the second run); the
#'   limit also bounds each round of the joint refinement of a Gaussian fit. An
#'   empirical-Bayes reweighting step runs at most five sweeps (it refreshes
#'   the weights, the Laplace correction and the variance updates after them).
#'   The variational stage of an empirical-Bayes fit runs at most
#'   \code{2 * max_iter} iterations. An empirical-Bayes fit whose iteration does
#'   not settle retires its weakest latent column and fits the next rank down
#'   (see \strong{A rank the data do not support}); a \code{max_iter} below 200
#'   disables that for an iteration that merely ran out of it.
#' @param tol convergence tolerance on the relative change in residual sum of
#'   squares (default \code{1e-6}). IRLS additionally requires the largest
#'   penalized likelihood score, scaled by the square root of its Fisher
#'   information plus penalty curvature, and the largest absolute log variance
#'   update to be at most
#'   \code{sqrt(tol)}.
#' @param lowrank_method inner solver for the multiplicative (low-rank) block:
#'   \code{"mm"} (default) weighted majorise-minimise; \code{"als"} alternating
#'   least squares; \code{"hybrid"} runs both and keeps the lower-objective
#'   result. \code{"als"}/\code{"hybrid"} converge faster than \code{"mm"} on
#'   strongly unbalanced longitudinal panels and are available for directed and
#'   bipartite models (symmetric fits always use \code{"mm"}). All three
#'   minimise the same objective, so the point estimate is unchanged for
#'   balanced data.
#' @param non_normal_method for the non-normal ALS families (it does not apply to
#'   a normal fit), \code{"irls"}
#'   (default for \code{binary}/\code{poisson}) runs iteratively reweighted
#'   least squares, giving a fast GLM AME fit with coefficients on
#'   the requested link scale (Poisson log, binary logit/probit).
#'   \code{"transform"} fits one fixed Gaussian working response
#'   (\code{log(y+1)} for Poisson, rank-normal scores for binary); its
#'   coefficients are on an uncalibrated working scale and are mainly useful
#'   for direction/ranking checks. A directed \code{R > 0} IRLS fit uses the
#'   hybrid low-rank solver internally because the IRLS weights are unbalanced.
#'   Uncertainty for either path comes from the bootstrap or sandwich
#'   covariance.
#' @param link link for \code{non_normal_method = "irls"} with a \code{binary}
#'   family: \code{"probit"} (default; matches \code{\link{ame}}/\code{\link{lame}})
#'   or \code{"logit"}. \code{poisson} always uses the log link;
#'   ignored otherwise.
#' @param linear_solver solver for the regression block: \code{"eigen"}
#'   (default) eigendecomposes the normal equations; \code{"qr"} uses a QR
#'   factorisation of the observed design, which is more stable for
#'   ill-conditioned covariates; \code{"auto"} picks \code{"qr"} when the design
#'   is ill-conditioned but full rank. All give the same answer for
#'   well-conditioned designs.
#' @param multistart for \code{R > 0} (a non-convex objective), \code{"none"}
#'   (default) fits from a single deterministic start; \code{"cheap"} (4 starts)
#'   and \code{"full"} (8 starts) also try random low-rank starts and keep the
#'   converged start with the lowest value of the criterion, evaluated at that
#'   start's own scales: for an IRLS fit the penalised deviance, which carries
#'   the factor penalty with its scale parameters (\code{"penalised
#'   deviance"}), and for an empirical-Bayes fit with a latent rank the
#'   Laplace marginal likelihood of the start, minus twice its logarithm
#'   (\code{"Laplace marginal likelihood"}; see \strong{Choosing the rank}),
#'   which integrates the effects and the factors and does not move with the
#'   scale of a factor column the data do not use (with the default variational
#'   stage every start is refined, and starts of different modes end at the same
#'   variational fit, so they neither differ in \code{multistart_distance} nor
#'   warn). An empirical-Bayes start
#'   without it is rated by the penalised deviance plus the penalty of the
#'   actor effects, the log determinant of their information and the
#'   hyperprior of their covariance (\code{"empirical-Bayes objective"}). The
#'   raw deviance is not compared: it always favours the least-shrunk optimum,
#'   the one with the largest
#'   scales. Other fits keep the lowest working SSE. The criterion
#'   (\code{multistart_criterion}), its value for each start
#'   (\code{multistart_sse}; for an empirical-Bayes fit it includes the
#'   hyperprior, which \code{objective} leaves out, or is the marginal
#'   likelihood), each start's
#'   observed-data deviance (\code{multistart_deviance}, IRLS fits), its
#'   distance from the kept start (\code{multistart_distance}, IRLS fits) and
#'   its convergence are recorded. The distance is the squared difference of
#'   the two fitted linear predictors, weighted by the mean of the two IRLS
#'   weights and summed over the observed cells: an approximation of the
#'   change in the deviance between the two fits. An IRLS fit warns of
#'   materially different local optima when a converged start is more than one
#'   deviance unit away from the kept one. Starts that end at the same fitted
#'   values can still differ in the penalised objective, which depends on the
#'   scales of the factors, and those settle more slowly than the fitted values
#'   (the Laplace marginal likelihood differs far less); that difference does
#'   not warn. Every multistart fit also keeps the
#'   coefficients each start reports (\code{multistart_coefficients}) and the
#'   distance of each start's multiplicative surface from the kept start's
#'   (\code{multistart_product_distance}), so that the spread of the estimates
#'   across the starts can be read from the fit. Reproducible given
#'   \code{seed}; the global RNG stream is left unchanged.
#' @param bias_correct how an ALS fit deals with the finite-sample
#'   (incidental-parameter) bias of a fixed-effects fit.
#'   \code{NULL} (default) fits \code{normal} models, and \code{binary} and
#'   \code{poisson} models on the IRLS path, by the empirical-Bayes estimator
#'   (\code{"eb"}), for directed, bipartite and undirected networks, with or
#'   without a latent rank: a Gaussian prior on the additive effects with its
#'   covariance estimated from the data, coefficients that solve the
#'   Laplace-marginal equations, and Laplace standard errors from \code{vcov()}
#'   for every coefficient (see \strong{Empirical-Bayes estimator}, \strong{Normal data}
#'   and \strong{Counts} in Details). A normal response is fitted as the model that
#'   \code{\link{ame}} and \code{\link{lame}} sample, with an exact Gaussian cell, and a count
#'   network as the Poisson-lognormal model that \code{\link{ame}}
#'   and \code{\link{lame}} sample. Other families are not corrected.
#'   \code{"eb"} asks for the empirical-Bayes estimator explicitly; it is an
#'   error (class \code{"lame_als_eb"}) for another family and for a binary or count fit with
#'   \code{non_normal_method = "transform"}. (The dynamic routes of
#'   \code{lame(method = "als")} have their own empirical-Bayes estimator for
#'   normal, binary and Poisson panels; see \code{bias_correct} in \code{\link{lame}}.)
#'   \code{"none"} returns the fixed-effects fit at any size and records the
#'   request (for a normal response, the least-squares fit of the earlier releases
#'   with the conditional sandwich covariance of the dyadic coefficients; for counts, the fit of the earlier releases, in which the effect
#'   of an actor without a count on a side runs to minus infinity); \code{"jackknife"} requires the split-sample correction of that
#'   fit, for any family, and is an error where it cannot be computed;
#'   \code{"bootstrap"} requires the parametric-bootstrap correction of a
#'   binary IRLS fit at any size. See \strong{Split-sample bias correction},
#'   \strong{Parametric-bootstrap bias correction} and \strong{Size boundary}
#'   in Details. A corrected fit keeps the uncorrected values in
#'   \code{fit$raw} and records the correction in \code{fit$bias_correct}; an
#'   empirical-Bayes fit records \code{method = "eb"} (\code{default = TRUE}
#'   when it was the default) and refits nothing.
#' @param verbose logical; print progress (default \code{TRUE}).
#' @param seed random seed (default \code{6886}) for the random starts of
#'   \code{multistart = "cheap"} / \code{"full"}; the caller's random number
#'   state is restored afterwards. The block coordinate descent itself is
#'   deterministic, so with \code{multistart = "none"} the point estimate does
#'   not depend on it.
#' @param dcor logical; whether a directed binary probit, Poisson or normal fit models the
#'   dyadic correlation of the latent errors of the two cells of a dyad, as
#'   \code{\link{ame}} and \code{\link{lame}} do with their \code{dcor = TRUE}.
#'   The default is \code{!symmetric}, the default of the samplers. With the
#'   empirical-Bayes estimator (the default of a binary fit) the correlation is
#'   part of the model that is fitted: the two cells of a dyad are one bivariate
#'   probit observation, \code{fit$VC["rho"]} is the correlation of that model, and the
#'   factors, the standard errors and the predictions are those of the model
#'   (see \strong{Dyadic correlation} in Details). \code{dcor = FALSE} fits
#'   independent cells, the model of the earlier release, and reports the
#'   correlation of the pairs of cells at the fitted predictor as before. A
#'   directed count fit takes the pairs at its fitted state (see
#'   \strong{Dyadic correlation of counts} in Details), and a directed normal
#'   fit takes the two cells of a dyad as one bivariate normal observation, in
#'   closed form (see \strong{Dyadic correlation of normal data} in Details; it
#'   costs 1.3 to 1.6 times the independent-cell fit). A fit
#'   that has no pair likelihood (undirected and bipartite networks, the other
#'   estimators of a binary fit, and the fixed-effects normal fit) keeps its own estimate of
#'   \code{rho}; an explicit \code{dcor = TRUE} with the logit link is an error,
#'   because the latent correlation is that of normal errors. The pair
#'   likelihood of a binary or count fit costs 1 to 2.5 times the independent-cell fit for one period and
#'   2 to 4.5 times for four, in proportion to the number of dyads and periods;
#'   \code{dcor = FALSE} is the fast path for very large networks. A directed
#'   normal network whose reciprocal pairs are
#'   all exactly equal (an undirected network given as a directed one) has an
#'   improper posterior for the correlation under the prior of the samplers, and
#'   is refused as \code{\link{ame}} refuses it (class
#'   \code{"lame_improper_normal_reciprocity"}): fit it with \code{symmetric = TRUE},
#'   or with \code{dcor = FALSE}.
#' @param prior \code{NULL} or a named list of priors. The estimators read one key, \code{s20}, the scale of the
#'   inverse-gamma prior of the variance of the noise of a Poisson or a normal fit (one degree of freedom on \code{s20};
#'   the default is 1 for counts and the variance of the observed responses for a normal response, the priors \code{\link{ame}} and
#'   \code{\link{lame}} draw it under). The other keys of the
#'   samplers' prior list have no counterpart in the estimators; they are reported and ignored. The prior of the
#'   multiplicative effects of a directed network is the samplers' default (\code{\link{ame}}, \code{Suv0}
#'   and \code{kappa0}).
#'
#' @return An object of class \code{"ame_als"}: a list with the point
#'   estimates. When \code{R > 0} the fitted multiplicative surface is
#'   exposed as \code{UVPM} (\code{U L U'} for symmetric fits, which also
#'   carry it as \code{ULUPM}; \code{U V'} otherwise), matching the slot
#'   names on MCMC fits.
#'   estimates \code{mu}, \code{beta}, \code{a}, \code{b}, \code{U}, \code{V}
#'   (\code{L} for symmetric models), the linear predictor \code{EZ},
#'   response-scale \code{fitted} values, working-scale \code{residuals},
#'   convergence information, and the variance-component vector \code{VC}.
#'   \code{R} is the rank of the fit; where the data did not support the rank
#'   asked for and columns were retired, \code{R_requested} (present only
#'   then) is the rank that was asked for (see \strong{A rank the data do not
#'   support}).
#'   For IRLS, \code{irls_iterations} counts accepted outer reweighting steps;
#'   \code{iterations} counts the final inner solve (the reweighting steps, for
#'   an empirical-Bayes fit). \code{irls_diagnostics}
#'   records the score and variance checks, their threshold, the stopping
#'   reason, and each accepted step's improvement under a common factor
#'   penalty. These checks concern the uncorrected fit; a split-sample
#'   correction also records each half's checks in
#'   \code{bias_correct$irls_diagnostics}. A changing empirical factor
#'   variance means successive penalized objective values need not decrease;
#'   an empirical-Bayes step is also accepted when it lowers the residual of
#'   the estimating equations that the checks evaluate. A fixed-effects fit
#'   that does not pass the checks within \code{max_iter}, nor within four times
#'   it, is an error of class \code{"lame_als_convergence"}: no fit is returned.
#'   (An empirical-Bayes fit that does not pass them warns with that class and
#'   has \code{converged = FALSE}.)
#'   For an empirical-Bayes fit (\code{bias_correct = "eb"}) \code{eb} is a
#'   list: \code{sab}, the fitted covariance of the actor effects (a
#'   \eqn{2 \times 2} matrix; an undirected fit has one variance and repeats
#'   it on the diagonal), \code{Sab0} and \code{eta0}, the scale and degrees
#'   of freedom of its inverse-Wishart prior, \code{gamma}, the coefficients of
#'   the node covariates in the prior mean of the effects, \code{sab_change},
#'   the relative change of \code{sab} in the last update, \code{sig2_u} and
#'   \code{sig2_v}, the fitted scales of the factor columns (one value per
#'   column and side; the isotropic scale of a directed fit repeated, the unit
#'   receiver scale of a bipartite fit, one scale per column of an undirected
#'   fit),
#'   \code{mu_int}, \code{a_raw} and \code{b_raw}, the intercept and the
#'   effects as the posterior is built from them (the effects are not
#'   re-centred), \code{log_marginal}, the Laplace log marginal likelihood at
#'   the posterior mode (see \strong{Choosing the rank} in Details),
#'   \code{vcov} and \code{vcov_cond}, the Laplace covariance of the
#'   coefficients (intercept, dyadic and node-covariate) jointly with the
#'   effects and the latent factors and conditional on the factors, \code{se},
#'   the square roots of the diagonal of \code{vcov}, \code{se_method} and
#'   \code{vcov_error}, the reason when the covariance could not be computed
#'   (\code{vcov()} then falls back to the sandwich), \code{estimator},
#'   which estimator produced the reported values (\code{"variational"}, the
#'   default, or \code{"laplace"}, the posterior mode alone),
#'   \code{variational}, the record of the variational iteration
#'   (\code{converged}, \code{termination}, \code{iterations},
#'   \code{residual}, the largest step a full update would still take, and
#'   \code{gain}, the objective's gain at the last step; \code{step}, the
#'   step length; \code{tolerance}; \code{elbo}, the evidence lower bound at
#'   the plug-in scales, with \code{ell}, \code{kl_actors} and
#'   \code{kl_factors}; \code{information}, the training events, the loadings
#'   and the events per loading of a binary or count fit with a latent rank, with \code{signal}, the double-centred
#'   fitted product squared and weighted by the fisher weight of each observed cell, \code{signal_per_loading} and
#'   \code{threshold}, 2.1, below which the product is withheld (a normal fit: \code{signal}, the
#'   squared norm of the double-centred fitted product summed over the periods, \code{loadings}, \code{per_loading},
#'   the signal in noise variances per loading, and \code{threshold}); \code{scale}, the rule that set
#'   the scale of the prior of the factors (\code{"cavity"} or
#'   \code{"plugin"}), the scale \code{s} (the mean over the
#'   columns of a bipartite or undirected fit), the correlation \code{rho} of
#'   the sender and receiver position of an actor (directed) and \code{sd_log},
#'   the standard deviation of the log scale that \code{vcov()} carries;
#'   \code{start}, \code{mode_log_marginal}),
#'   \code{factor_posterior}, the covariances of the posterior of the factor
#'   rows (arrays \code{Sig} and \code{Om}, or \code{Sig} and the factor
#'   \code{F} with its signs \code{sg} for an undirected fit; the cavity
#'   covariances, with the covariances of the fixed point of the updates beside
#'   them as \code{Sig_mf} and \code{Om_mf}) and the prior of
#'   the fit, \code{suv}, the fitted \eqn{2R \times 2R} covariance of the factor
#'   rows of a directed network, \code{laplace}, the parts of
#'   \code{log_marginal}, and, for a directed binary probit fit with a dyadic
#'   correlation, \code{pair}: the correlation \code{rho}, its posterior
#'   standard deviation \code{sd}, the score and curvature of the expected log
#'   likelihood in it, the number of dyads with both cells observed
#'   \code{n_pairs}, the nodes a side of the quadrature \code{quad} (zero for a normal pair, which has none), the rules
#'   the quadrature was refined to after a stall \code{refined} (empty where it was not) and the
#'   expected cross weights of the pairs; \code{retirement}, present only when columns were retired: a list
#'   with the \code{requested_rank}, the \code{effective_rank}, whether the fit
#'   of the effective rank \code{settled}, and \code{steps}, one record per
#'   retired column (the \code{rank} it was retired at, its \code{column},
#'   its \code{strength} and \code{share} of the strongest column, and the
#'   \code{stage} and \code{termination} of the iteration that did not settle);
#'   for a normal fit, \code{ve}, the variance of the
#'   noise (also \code{VC["ve"]} and \code{fit$s2}) and \code{ve_prior}, the degrees of freedom and
#'   the scale of its prior; and, for a Poisson fit, \code{ve},
#'   the variance of the lognormal noise of the count (also \code{VC["ve"]},
#'   \code{fit$s2} and \code{fit$pois_latent_var}), and \code{ve_search}, the
#'   search that found it: \code{ve}, \code{objective} (the Laplace log marginal
#'   likelihood plus the log prior density of \code{ve} at the maximum),
#'   \code{curve} (a data frame of the values tried: \code{ve},
#'   \code{objective}, whether the mode there \code{solved} its equations and
#'   its \code{iterations}) and \code{evaluations}.
#'   A fit with \code{multistart = "cheap"} or
#'   \code{"full"} (and \code{R > 0}) also records, for each start,
#'   \code{multistart_converged}, \code{multistart_criterion} and
#'   \code{multistart_sse} (the value of that criterion), the kept start
#'   \code{multistart_selected} and, for IRLS fits, \code{multistart_deviance}
#'   and \code{multistart_distance} (the distance of the fitted predictor from
#'   the kept start's, in deviance units); and for every multistart fit
#'   \code{multistart_coefficients}, a starts by coefficients matrix of what
#'   each start reports (intercept, dyadic and node-covariate coefficients, with
#'   the names and on the scale of \code{coef(fit)}) whose row for the kept
#'   start is \code{coef(fit)} (the uncorrected \code{fit$raw$coefficients} of
#'   a bias-corrected fit, since the correction is made on the kept start
#'   only), and \code{multistart_product_distance}, the Frobenius norm of the
#'   difference between the \code{UVPM} each start reports (before any bias
#'   correction) and the kept start's, over the norm of the kept start's: zero
#'   for the kept start. A kept \code{UVPM} of no size, when the factors of the
#'   starts shrink to nothing, is measured against
#'   \code{sqrt(.Machine$double.eps)}, so that starts that agree there are at
#'   distance zero. The surfaces themselves are not stored.
#'   The variance-component vector has five descriptive entries:
#'   \describe{
#'     \item{\code{va}, \code{vb}}{variances of the sender and receiver
#'       additive effects. The empirical variance of the \emph{estimated}
#'       effects carries their estimation error as well, so the first-order
#'       sampling variance of the additive block --- \code{ve} times the
#'       centred trace of the inverse two-way weighted normal equations --- is
#'       subtracted, floored at zero, and \code{va}/\code{vb} estimate the
#'       variance of the effects themselves, the quantity \code{\link{ame}} /
#'       \code{\link{lame}} report. The uncorrected empirical variances are
#'       kept as \code{vc_raw} and the correction as \code{vc_noise}; it is
#'       computed conditional on the fitted multiplicative term. For an
#'       empirical-Bayes fit they are the fitted variances of the effects
#'       themselves (\code{fit$eb$sab}).}
#'     \item{\code{cab}}{covariance of the \emph{estimated} sender and
#'       receiver effects, uncorrected (\code{NA} for symmetric or bipartite
#'       models); for an empirical-Bayes fit, their fitted covariance.}
#'     \item{\code{rho}}{dyadic reciprocity, distinct from the sender/receiver
#'       correlation \code{cor(a, b)}. For a directed binary probit fit with the
#'       empirical-Bayes estimator and \code{dcor = TRUE} it is the correlation of
#'       the fitted model (\strong{Dyadic correlation}), and
#'       \code{rho_estimation} holds its posterior standard deviation and
#'       \code{method = "variational pair likelihood"}. For other directed binary
#'       IRLS fits, this is
#'       the Gaussian copula correlation estimated from the paired binary
#'       likelihood at the reported marginal probabilities. With a probit
#'       link it is the latent error correlation. The working-residual
#'       correlation is retained as \code{rho_working}; \code{rho_estimation}
#'       records the number of reciprocal pairs, likelihood, and boundary or
#'       unidentified state. Without informative reciprocal pairs \code{rho}
#'       is \code{NA} and conditional network simulation is refused. Poisson
#'       IRLS fits estimate latent correlation from the weighted count covariance
#'       identity \eqn{\mathrm{cov}(Y_{ij}/\mu_{ij}-1,Y_{ji}/\mu_{ji}-1) = e^{\rho v_e}-1}.
#'       Independent Poisson noise contributes to the marginal variances but
#'       not this covariance. (This is the estimate of a fit without the pair
#'       stage, \code{bias_correct = "none"} or \code{dcor = FALSE}; a
#'       directed empirical-Bayes count fit estimates the correlation from the
#'       pairs, see \strong{Dyadic correlation of counts}.) With zero latent
#'       variance, correlation is
#'       unidentified and independent count simulation remains valid. Normal
#'       and transformed-response fits use reciprocal residual correlation.}
#'     \item{\code{ve}}{residual variance, on the model's own scale and
#'       degrees-of-freedom-corrected (\eqn{\mathrm{SSE}/(n_{obs} -
#'       \mathrm{df})}) so it matches the MCMC posterior-mean \code{s2}. For a
#'       \code{binary} IRLS fit the probit/logit model fixes the latent error
#'       variance at 1 by identification, so \code{ve} is reported as
#'       \code{1}; the working-scale GLM dispersion (about 1 under correct
#'       specification) is kept separately as \code{ve_working}. For a
#'       \code{poisson} IRLS fit \code{ve} is the variance of the lognormal
#'       latent layer of the model \code{\link{ame}} / \code{\link{lame}} fit
#'       (\eqn{y \sim \mathrm{Pois}(e^{z})}, \eqn{z \sim N(\eta, v_e)}); an
#'       empirical-Bayes fit estimates it by maximising the Laplace marginal
#'       likelihood (see \strong{Counts}), and the fixed-effects fit
#'       (\code{bias_correct = "none"}) from the count overdispersion
#'       \eqn{\mathrm{var}[y] = \mathrm{E}[y] + (e^{v_e} - 1)\mathrm{E}[y]^2}
#'       and \code{0} for equidispersed counts; the working-response mean
#'       square is again \code{ve_working}.}
#'   }
#'   \strong{Poisson level scale.} The fixed-effects IRLS fit
#'   (\code{bias_correct = "none"}) solves the marginal log-mean
#'   \eqn{\mathrm{E}[y] = e^{\eta}}, which is the latent predictor plus
#'   \eqn{v_e/2}. The reported intercept and \code{EZ} subtract that Jensen
#'   term, so they are on the same latent scale as \code{\link{ame}} /
#'   \code{\link{lame}}; the shift is stored as \code{pois_latent_var} and
#'   added back by \code{fitted}, \code{predict(type = "response")} and
#'   \code{simulate}, which therefore do not change. The dyadic slopes are
#'   unaffected: the term is constant across dyads. The empirical-Bayes fit
#'   works on the latent scale throughout and needs no shift.
#'
#'   These are summaries of the point estimates, not draws from a posterior.
#'   See \code{\link{ame_als_bootstrap}} for
#'   uncertainty and \code{\link{vcov.ame_als}} for a fast analytic
#'   covariance of the regression coefficients.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The iterative block
#' coordinate descent estimator adapted here originates with that work
#' (implemented in \code{sir::sir_alsfit()}).
#'
#' @seealso \code{\link{lame_als}} for longitudinal networks,
#'   \code{\link{ame_als_bootstrap}} for bootstrap uncertainty,
#'   \code{\link{ame}} for the full MCMC estimator.
#'
#' @examples
#' Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
#' fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' coef(fit)
#'
#' @param bootstrap integer: if > 0, additionally run \code{bootstrap} replicates
#'   of the parametric or block bootstrap (via \code{\link{ame_als_bootstrap}})
#'   after the point fit and attach the result as \code{fit$bootstrap}. The
#'   downstream accessors (\code{\link{confint.ame_als}}, \code{summary},
#'   \code{print}) then surface bootstrap intervals instead of the
#'   anti-conservative sandwich Wald intervals. At least 30 replicates are
#'   needed: a smaller positive number is an error of class
#'   \code{"lame_als_bootstrap_replicates"} raised before anything is fitted.
#'   A binary network fitted with \code{bias_correct = "none"} has no
#'   bootstrap (an error of class \code{"lame_als_bootstrap_bias"}, also
#'   raised before anything is fitted): see \code{\link{ame_als_bootstrap}}.
#'   Default \code{0} (no
#'   bootstrap) -- bootstrap is expensive, so it is not imposed on a user
#'   who just wants a quick fit.
#' @param bootstrap_type character: \code{"parametric"} (default) or
#'   \code{"block"} -- the bootstrap scheme to use when \code{bootstrap > 0}.
#'   The block scheme returns no standard error or interval for the
#'   parameters that resampling time slices cannot vary (the node-covariate
#'   coefficients and the variance components of the actor effects, among
#'   others; \code{fit$bootstrap$withheld} lists them with the reason), and
#'   the parametric scheme gives them; see \code{\link{ame_als_bootstrap}}.
#' @param bootstrap_block_length integer: block length for the block
#'   bootstrap; only used when \code{bootstrap_type = "block"}.
#' @param bootstrap_seed optional integer seed for the bootstrap (the point
#'   fit uses \code{seed}).
#'
#' @section Coverage relative to \code{\link{ame}} / \code{\link{lame}}:
#' The ALS estimator is a fast, frequentist point estimator. It covers most
#' static AME workflows, and the top-level \code{lame(..., method = "als")}
#' dispatcher covers several dynamic workflows, but posterior-specific features
#' still require the MCMC estimator:
#' \itemize{
#'   \item \strong{Families.} ALS supports \code{normal}, \code{binary}, and
#'     \code{poisson}. For \code{ordinal}, \code{cbin},
#'     and \code{frn}, fall back to \code{\link{ame}} /
#'     \code{\link{lame}}.
#'   \item \strong{Dynamic effects.} \code{\link{lame_als}} itself fits a
#'     static model pooled across time slices. The top-level dispatcher
#'     \code{lame(..., method = "als")} routes supported dynamic requests to a
#'     dynamic estimator for normal, binary, and Poisson panels (an
#'     empirical-Bayes fit by default; see
#'     \code{bias_correct} in \code{\link{lame}}),
#'     including named panels where actors enter or exit: \code{dynamic_ab},
#'     selected intercept/dyadic/node
#'     \code{dynamic_beta}, AR(1) \code{dynamic_uv} for directed, symmetric,
#'     and bipartite panels, and Student-t \code{dynamic_uv} for directed and
#'     symmetric panels (bipartite t requests are refused).
#'     \code{dynamic_uv_kind = "snap"} has no ALS route and is refused.
#'     Node-covariate coefficients use the same orthogonal additive-effect
#'     decomposition as \code{\link{lame_als}}; dynamic node coefficients use
#'     period-specific node values when selected by \code{dynamic_beta}, while
#'     static node coefficients use per-actor means. A time-varying bipartite
#'     interaction matrix (\code{dynamic_G = TRUE}) is refused on every route:
#'     use \code{dynamic_uv = TRUE} and read \code{fit$UVPM}. A dynamic normal fit whose latent product the
#'     data do not determine withholds it (\code{fit$eb$product_withheld}; see \strong{Choosing the rank of a
#'     normal fit} in \code{\link{ame_als}}). Changing actor
#'     composition requires row and column names on every slice so actors can
#'     be aligned; smoothing penalties are broken across actor-entry gaps (the
#'     empirical-Bayes fit treats an absent actor-period as a free AR state).
#'     Rank/censored dynamic families remain on the MCMC path. Static
#'     \code{ame_als()} / \code{lame_als()} fits still use a
#'     single latent rank \code{R}; the dynamic bipartite ALS dispatcher
#'     honours separate \code{R_row} and \code{R_col} values.
#'   \item \strong{Priors.} \code{ame_als()}, \code{lame_als()} and
#'     \code{lame(..., method = "als")} read one entry of \code{prior},
#'     \code{s20}, the scale of the prior of the noise variance of a normal or
#'     Poisson fit, and the dynamic empirical-Bayes fits read the hyperprior
#'     entries of \code{\link{lame}} in their EM. Every other entry of the
#'     samplers' \code{prior} list and \code{g = ...} are MCMC-only; the
#'     dispatcher \code{ame(..., method = "als")} warns and ignores \code{prior}
#'     and \code{g} (pass \code{prior = list(s20 = )} to \code{ame_als()}). The
#'     empirical-Bayes fits apply the default additive-effect prior of
#'     \code{\link{ame}} / \code{\link{lame}} with its covariance estimated.
#'   \item \strong{Posterior quantities.} No \code{$BETA} / \code{$VC} posterior
#'     draws on the point fit. Use \code{bootstrap = N} for a sampling
#'     distribution analogue; the draws are bootstrap, not Bayesian.
#'   \item \strong{Multi-chain.} Not applicable: an ALS fit is one point fit, \code{chains} has no effect, and \code{n_chains} is ignored with a warning.
#'   \item \strong{Convergence diagnostics.} No Rhat / ESS / \code{\link{trace_plot}}.
#' }
#'
#' Features that work the same on ALS fits: \code{coef},
#' \code{vcov} (Laplace covariance of an empirical-Bayes fit, sandwich on the
#' regression block of a fixed-effects fit), \code{confint} (auto-routes
#' bootstrap intervals when present, Wald intervals otherwise), \code{predict},
#' \code{fitted}, \code{residuals}, \code{summary}, \code{print},
#' \code{nobs}, \code{\link{simulate.ame_als}},
#' \code{\link{gof_plot.ame_als}}, \code{\link{ab_plot.ame_als}},
#' \code{\link{uv_plot}}, \code{\link{latent_positions}}.
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
ame_als <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	R = 0,
	family = "normal",
	mode = c("unipartite", "bipartite"),
	symmetric = FALSE,
	max_iter = 200,
	tol = 1e-6,
	lowrank_method = c("mm", "als", "hybrid"),
	non_normal_method = c("irls", "transform"),
	link = c("probit", "logit"),
	linear_solver = c("eigen", "qr", "auto"),
	multistart = c("none", "cheap", "full"),
	bias_correct = NULL,
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	verbose = TRUE,
	seed = 6886,
	dcor = !symmetric,
	prior = NULL
) {
	family_missing <- missing(family)
	mode_missing <- missing(mode)
	symmetric_missing <- missing(symmetric)
	dcor_given <- !missing(dcor)
	xdyad_name <- .ae_dyad_expr_name(substitute(Xdyad))
	netify_input <- .lame_prepare_netify_input(
		Y = Y,
		Xdyad = Xdyad,
		Xrow = Xrow,
		Xcol = Xcol,
		mode = mode,
		family = family,
		mode_missing = mode_missing,
		family_missing = family_missing,
		symmetric = symmetric,
		symmetric_missing = symmetric_missing,
		longitudinal = FALSE,
		caller = "ame_als",
		fit_method = "als"
	)
	if (isTRUE(netify_input$converted)) {
		Y <- netify_input$Y
		Xdyad <- netify_input$Xdyad
		Xrow <- netify_input$Xrow
		Xcol <- netify_input$Xcol
		mode <- netify_input$mode
		family <- netify_input$family
		symmetric <- netify_input$symmetric
	}

	mode <- match.arg(mode)
	lowrank_method <- match.arg(lowrank_method)
	non_normal_method <- match.arg(non_normal_method)
	link <- match.arg(link)
	linear_solver <- match.arg(linear_solver)
	multistart <- match.arg(multistart)
	bias_correct <- .ae_check_bias_correct(bias_correct)
	bootstrap_type <- match.arg(bootstrap_type)
	if (!is.numeric(bootstrap) || length(bootstrap) != 1L || bootstrap < 0L ||
			!isTRUE(all.equal(bootstrap, round(bootstrap)))) {
		cli::cli_abort("{.arg bootstrap} must be a non-negative integer (0 = no bootstrap).")
	}
	bootstrap <- as.integer(bootstrap)
	.ae_check_bootstrap_n(bootstrap)
	.ae_check_bootstrap_bias(bootstrap, family, non_normal_method, bias_correct)
	# name a single 2-d dyadic covariate after the argument, as ame() does
	if (!is.null(Xdyad) && is.matrix(Xdyad) && !is.null(xdyad_name) &&
			!isTRUE(netify_input$converted)) {
		Xdyad <- .ae_name_dyad_slice(Xdyad, xdyad_name)
		if (verbose) {
			cli::cli_inform(c(
				"i" = "{.arg Xdyad} was a 2-D matrix; wrapping as {.code array(., dim = c(n, n, 1))} with slice name {.val {xdyad_name}}."
			))
		}
	}
	# remind interactive users that a binary rank > 0 fit is a penalized point
	# estimate, with intervals from MCMC or the bootstrap
	if (verbose && is.character(family) && length(family) == 1L &&
			family %in% c("binary", "bin") &&
			(!identical(bias_correct, "eb") && !is.null(bias_correct)) &&
			is.numeric(R) && length(R) == 1L && isTRUE(is.finite(R)) && R > 0) {
		cli::cli_inform(c(
			"i" = "Binary ALS with {.arg R} > 0 is a penalized (MAP-style) point estimate.",
			"i" = "For calibrated inference use {.code bootstrap = N} or {.code ame(..., method = \"mcmc\")}."
		))
	}
	dcor <- .ae_check_dcor(dcor)
	prior <- .ae_check_prior(prior, family)
	fit <- .ae_run_budget(
		Y,
		Xdyad,
		Xrow,
		Xcol,
		R,
		family,
		mode,
		symmetric,
		max_iter,
		tol,
		verbose,
		seed,
		longitudinal = FALSE,
		call = match.call(),
		lowrank_method = lowrank_method,
		non_normal_method = non_normal_method,
		link = link,
		linear_solver = linear_solver,
		multistart = multistart,
		bias_correct = bias_correct,
		dcor = dcor,
		dcor_given = dcor_given,
		prior = prior
	)
	if (bootstrap > 0L) {
		fit <- .ae_attach_bootstrap(
			fit,
			bootstrap,
			bootstrap_type,
			bootstrap_block_length,
			bootstrap_seed,
			verbose
		)
	}
	fit
}


#' Fast (MCMC-free) AME estimation for a longitudinal network
#'
#' @description
#' Fits an additive and multiplicative effects (AME) model to a longitudinal
#' (replicated) network by \strong{iterative block coordinate descent},
#' producing a fast point estimate with no MCMC and no credible intervals. For
#' normal, binary and Poisson data the default is an empirical-Bayes / variational
#' estimator of the model that \code{\link{lame}} samples (see
#' \code{bias_correct}); its \code{vcov()} and \code{confint()} are Laplace
#' approximations.
#' This is the longitudinal counterpart of \code{\link{ame_als}}; the
#' estimated effects (\code{mu}, \code{beta}, \code{a}, \code{b}, \code{U},
#' \code{V}) are \emph{static} (pooled across time), as in a non-dynamic
#' \code{\link{lame}} fit.
#'
#' The estimation algorithm adapts the iterative block coordinate descent
#' estimator of the Social Influence Regression model of Hoff & Minhas
#' (\code{sir::sir_alsfit()}) to the AME model; it is a port and adaptation,
#' not original \pkg{lame} methodology. See \code{\link{ame_als}} for the
#' algorithm details.
#'
#' @param Y a list of \code{t} relational matrices, or a 3d array
#'   \code{[n_row, n_col, t]}. A named list may have changing actor
#'   composition; every slice must carry row and column names so actors can be
#'   aligned to the union panel. Unnamed lists and arrays are treated
#'   positionally and must have one fixed layout. A longitudinal \pkg{netify}
#'   object is also accepted and converted with
#'   \code{netify::to_lame(lame = TRUE)}.
#' @param Xdyad a list of \code{T} dyadic covariate matrices/arrays, or \code{NULL}.
#' @param Xrow a list of \code{T} row/sender covariate matrices, or \code{NULL}.
#' @param Xcol a list of \code{T} column/receiver covariate matrices, or \code{NULL}.
#' @inheritParams ame_als
#'
#' @return An object of class \code{"ame_als"}; see \code{\link{ame_als}},
#'   including the per-start fields (\code{multistart_coefficients},
#'   \code{multistart_product_distance} and the others) that a fit with
#'   \code{multistart} records.
#'
#' @references
#' Minhas, S. and Hoff, P. D. (2025). Decomposing Network Dynamics: Social
#' Influence Regression. \emph{Political Analysis}. The iterative block
#' coordinate descent estimator adapted here originates with that work
#' (implemented in \code{sir::sir_alsfit()}).
#'
#' @seealso \code{\link{ame_als}}, \code{\link{ame_als_bootstrap}},
#'   \code{\link{lame}} for the full MCMC estimator,
#'   \code{\link{als_dynamic_beta}} for a regression-only smoother
#'   that estimates \emph{only} a time-varying \eqn{\beta_t} (no
#'   \eqn{a, b, U, V}; not a special case of this function).
#'
#' @examples
#' Y <- replicate(4, { m <- matrix(rnorm(400), 20, 20); diag(m) <- NA; m },
#'                simplify = FALSE)
#' fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' coef(fit)
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
lame_als <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	R = 0,
	family = "normal",
	mode = c("unipartite", "bipartite"),
	symmetric = FALSE,
	max_iter = 200,
	tol = 1e-6,
	lowrank_method = c("mm", "als", "hybrid"),
	non_normal_method = c("irls", "transform"),
	link = c("probit", "logit"),
	linear_solver = c("eigen", "qr", "auto"),
	multistart = c("none", "cheap", "full"),
	bias_correct = NULL,
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	verbose = TRUE,
	seed = 6886,
	dcor = !symmetric,
	prior = NULL
) {
	family_missing <- missing(family)
	mode_missing <- missing(mode)
	symmetric_missing <- missing(symmetric)
	dcor_given <- !missing(dcor)
	netify_input <- .lame_prepare_netify_input(
		Y = Y,
		Xdyad = Xdyad,
		Xrow = Xrow,
		Xcol = Xcol,
		mode = mode,
		family = family,
		mode_missing = mode_missing,
		family_missing = family_missing,
		symmetric = symmetric,
		symmetric_missing = symmetric_missing,
		longitudinal = TRUE,
		caller = "lame_als",
		fit_method = "als"
	)
	if (isTRUE(netify_input$converted)) {
		Y <- netify_input$Y
		Xdyad <- netify_input$Xdyad
		Xrow <- netify_input$Xrow
		Xcol <- netify_input$Xcol
		mode <- netify_input$mode
		family <- netify_input$family
		symmetric <- netify_input$symmetric
	}

	mode <- match.arg(mode)
	lowrank_method <- match.arg(lowrank_method)
	non_normal_method <- match.arg(non_normal_method)
	link <- match.arg(link)
	linear_solver <- match.arg(linear_solver)
	multistart <- match.arg(multistart)
	bias_correct <- .ae_check_bias_correct(bias_correct)
	bootstrap_type <- match.arg(bootstrap_type)
	if (!is.numeric(bootstrap) || length(bootstrap) != 1L || bootstrap < 0L ||
			!isTRUE(all.equal(bootstrap, round(bootstrap)))) {
		cli::cli_abort("{.arg bootstrap} must be a non-negative integer (0 = no bootstrap).")
	}
	bootstrap <- as.integer(bootstrap)
	.ae_check_bootstrap_n(bootstrap)
	.ae_check_bootstrap_bias(bootstrap, family, non_normal_method, bias_correct)
	dcor <- .ae_check_dcor(dcor)
	prior <- .ae_check_prior(prior, family)
	fit <- .ae_run_budget(
		Y,
		Xdyad,
		Xrow,
		Xcol,
		R,
		family,
		mode,
		symmetric,
		max_iter,
		tol,
		verbose,
		seed,
		longitudinal = TRUE,
		call = match.call(),
		lowrank_method = lowrank_method,
		non_normal_method = non_normal_method,
		link = link,
		linear_solver = linear_solver,
		multistart = multistart,
		bias_correct = bias_correct,
		dcor = dcor,
		dcor_given = dcor_given,
		prior = prior
	)
	if (bootstrap > 0L) {
		fit <- .ae_attach_bootstrap(
			fit,
			bootstrap,
			bootstrap_type,
			bootstrap_block_length,
			bootstrap_seed,
			verbose
		)
	}
	fit
}

# S3 methods for ame_als objects

#' Extract coefficients from a fast AME fit
#'
#' Returns the point-estimated regression coefficients (intercept first) of an
#' \code{ame_als} fit. There is no posterior here; for uncertainty use
#' \code{\link{ame_als_bootstrap}}.
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return A named numeric vector of coefficients.
#' @method coef ame_als
#' @export
coef.ame_als <- function(object, ...) {
	object$coefficients
}

#' Extract fitted values from a fast AME fit
#'
#' Returns response-scale fitted values: a single matrix for a cross-sectional
#' (\code{ame_als}) fit, or a list of matrices for a longitudinal
#' (\code{lame_als}) fit.
#'
#' @details
#' For a binary empirical-Bayes fit (the default) the fitted values are
#' predictive probabilities: the expectation of the response over the posterior
#' of the linear predictor, \eqn{\Phi(m / \sqrt{1 + v})} for a probit link, with
#' \eqn{m} the posterior mean and \eqn{v} the posterior variance of the predictor
#' (see \code{\link{predict.ame_als}}). For a Poisson fit they are the plug-in
#' mean counts, \eqn{e^{\eta + v_e / 2}}{exp(eta + ve / 2)} at the fitted
#' predictor \eqn{\eta} (the posterior mean) and the variance \eqn{v_e}{ve} of
#' the noise of the count, as the samplers report them and as the dynamic
#' empirical-Bayes fit does: the posterior predictive mean of a count is not
#' finite in this model (the inverse-gamma tail of \eqn{v_e}{ve}, and the
#' inverse-Wishart tail of the scale of the factors for a cell the data do not
#' inform), so its value would be a property of the numerical integration and of
#' the prior; the plug-in mean is finite everywhere. \code{plug_in = FALSE}
#' returns the predictive mean of a count that the integration over the
#' predictor gives, \eqn{e^{m + (v + v_e) / 2}}{exp(m + (v + ve) / 2)}, and
#' \code{plug_in = TRUE} the plug-in values of a binary fit (the response mean
#' at the fitted predictor); the plug-in values are what \code{fit$fitted} and
#' \code{fit$YPM} hold and \code{simulate()} uses. \code{residuals()} follows
#' \code{fitted()}. The mean of a normal response is its
#' predictor, so a normal fit's fitted values are the posterior means of the
#' predictor, with or without \code{plug_in}.
#'
#' A dynamic fit (\code{lame(method = "als")} with \code{dynamic_ab},
#' \code{dynamic_beta} or \code{dynamic_uv}) follows the same rules. The
#' empirical-Bayes binary fit returns the predictive probability, with \eqn{v}
#' the variance of the cell's predictor: the additive effects given the
#' coefficients (the prior variance for an actor with no observed cell), the
#' product of the positions, and the coefficients' part at the design asked
#' for; \code{plug_in = TRUE} returns \eqn{\Phi(m)}, which \code{fit$fitted}
#' keeps. A fit stored without the variance record is an error of class
#' \code{"lame_predictive_missing"}. The empirical-Bayes count fit returns the
#' plug-in mean \eqn{e^{\eta + v_e / 2}}{exp(eta + ve / 2)}, the normal fit the
#' predictor, and a penalised fit (\code{bias_correct = "none"}) the response
#' mean at its predictor (\eqn{\Phi(\eta)}, \eqn{e^{\eta}}). The dynamic route
#' has no exposure and no forecasts (\code{predict(h > 0)} is an error).
#'
#' @param object an \code{ame_als} object.
#' @param plug_in logical; for an empirical-Bayes binary or count fit, return the plug-in
#'   fitted values (\code{TRUE}) or the predictive means that integrate the link
#'   over the uncertainty of the predictor (\code{FALSE}). The default is
#'   \code{TRUE} for a Poisson fit and \code{FALSE} for a binary fit.
#' @param ... ignored.
#' @return A matrix or list of matrices of fitted values.
#' @method fitted ame_als
#' @export
fitted.ame_als <- function(object, plug_in = NULL, ...) {
	plug_in <- .ae_resolve_plug_in(object, plug_in)
	if (!isTRUE(plug_in) && .ae_is_eb_fit(object)) {
		pr <- .ae_eb_predict_response(object, X = NULL, type = "response")
		if (!is.null(pr)) return(pr)
	}
	if (isTRUE(object$longitudinal)) object$fitted else object$fitted[[1]]
}

# default plug_in for fitted(), residuals() and predict(): the plug-in mean for
# Poisson (its predictive mean is not finite), the predictive probability for binary
.ae_resolve_plug_in <- function(object, plug_in = NULL) {
	if (is.null(plug_in)) return(identical(object$family, "poisson"))
	if (!is.logical(plug_in) || length(plug_in) != 1L || is.na(plug_in)) {
		cli::cli_abort("{.arg plug_in} must be {.val TRUE}, {.val FALSE} or {.val NULL} (the default of the family).")
	}
	plug_in
}

# whether an EB fit carries what its posterior is built from
.ae_is_eb_fit <- function(object) {
	!inherits(object, "lame_dynamic_als") &&
		identical(object$bias_correct$method, "eb") && !is.null(object$eb$mu_int) &&
		!is.null(object$obs_weights) && object$family %in% c("binary", "poisson")
}

# predictive probabilities (type = "response") or posterior mean predictor
# (type = "link") of an EB fit, as a matrix or list of matrices
.ae_eb_predict_response <- function(object, X = NULL, type = "response") {
	pr <- tryCatch(.ae_eb_predictive(object, X), error = function(e) e)
	if (inherits(pr, "error")) {
		cli::cli_warn(c(
			"The predictive distribution of this empirical-Bayes fit could not be computed: {conditionMessage(pr)}",
			"i" = "Returning the plug-in values."
		))
		return(NULL)
	}
	out <- lapply(seq_len(object$n_time), function(t) {
		m <- pr$mean[, , t]
		val <- if (identical(type, "link")) m else
			.ae_eb_integrate(m, pr$var[, , t], .ae_default(object$link, "probit"), ve = .ae_default(object$eb$ve, 0))
		if (!isTRUE(object$mode == "bipartite") && nrow(val) == ncol(val)) diag(val) <- NA_real_
		dimnames(val) <- dimnames(object$EZ[[t]])
		val
	})
	names(out) <- names(object$EZ)
	if (isTRUE(object$longitudinal)) out else out[[1]]
}

# inverse link to the response scale; `latent_var` is the lognormal variance of
# a Poisson IRLS fit, so the count mean is exp(ez + latent_var / 2)
.ae_link_inverse <- function(ez, family, link, latent_var = 0) {
	if (!is.null(link)) {
		switch(
			paste(family, link),
			"poisson log" = exp(ez + latent_var / 2),
			"binary logit" = stats::plogis(ez),
			"binary probit" = stats::pnorm(ez),
			ez
		)
	} else switch(
		family,
		normal = ez,
		poisson = .poisson_inverse_log1p(ez),
		binary = stats::pnorm(ez),
		ez
	)
}

#' Residuals from a fast AME fit
#'
#' @description
#' Returns residuals from an \code{\link{ame_als}} or
#' \code{\link{lame_als}} fit.
#'
#' @details
#' \code{type = "response"} (default) returns observed \code{Y} minus the
#' response-scale fitted values, so \code{residuals()} reconciles with
#' \code{fitted()}: for a binary empirical-Bayes fit the predictive
#' probabilities, or with \code{plug_in = TRUE} the plug-in ones, and for a
#' Poisson fit the plug-in mean counts, or with \code{plug_in = FALSE} the
#' predictive means that integrate the predictor.
#' \code{type = "working"} returns the Gaussian working-scale
#' residuals the estimator actually minimised (identical to \code{"response"}
#' for the \code{normal} family).
#'
#' @param object an \code{ame_als} object.
#' @param type \code{"response"} (default) or \code{"working"}.
#' @param plug_in logical; for an empirical-Bayes fit, take the response
#'   residuals from the plug-in fitted values (\code{TRUE}) or from the
#'   predictive means (\code{FALSE}); the default is the one of
#'   \code{\link{fitted.ame_als}}: \code{TRUE} for a Poisson fit, \code{FALSE}
#'   for a binary fit.
#' @param ... ignored.
#' @return A matrix (cross-sectional) or list of matrices (longitudinal).
#' @method residuals ame_als
#' @export
residuals.ame_als <- function(
	object,
	type = c("response", "working"),
	plug_in = NULL,
	...
) {
	type <- match.arg(type)
	res <- if (type == "working") {
		object$residuals
	} else {
		fv <- fitted.ame_als(object, plug_in = plug_in)
		if (!isTRUE(object$longitudinal)) fv <- list(fv)
		lapply(seq_len(object$n_time), function(t) {
			r <- object$Y[, , t] - fv[[t]]
			dimnames(r) <- dimnames(fv[[t]])
			r
		})
	}
	if (isTRUE(object$longitudinal)) res else res[[1]]
}

#' Predictions from a fast AME fit
#'
#' @description
#' Returns predictions from an \code{\link{ame_als}} or
#' \code{\link{lame_als}} fit: the fitted linear predictor on the link
#' scale (\code{type = "link"}) or the response scale (\code{type = "response"},
#' the default).
#'
#' @details
#' With \code{newdata = NULL} the training-data fitted values are returned.
#' Supplying \code{newdata} -- a dyadic covariate array with the \emph{same}
#' actors and dimensions as the fitted design -- substitutes the dyadic
#' covariate contribution while holding the intercept, additive effects and
#' multiplicative term fixed; it predicts counterfactual dyadic covariates for
#' the same network, not out-of-sample actors.
#'
#' For a binary empirical-Bayes fit (the default) the response-scale
#' predictions are predictive probabilities: the response integrated over the
#' posterior of the linear predictor, \eqn{\Phi(m / \sqrt{1 + v})} for a probit
#' link, rather than applied at its mean (see \strong{Predictions} in
#' \code{\link{ame_als}}). For a Poisson fit they are the plug-in mean counts
#' \eqn{e^{m + v_e / 2}}{exp(m + ve / 2)}, as the samplers report them: the
#' posterior predictive mean of a count is not finite in this model (the
#' inverse-gamma tail of \eqn{v_e}{ve}, the inverse-Wishart tail of the scale of
#' the factors for a cell the data do not inform), and the plug-in mean is
#' finite everywhere. \code{plug_in = FALSE} returns, for a count, the mean that
#' integrates the predictor, \eqn{e^{m + (v + v_e) / 2}}{exp(m + (v + ve) / 2)},
#' and \code{plug_in = TRUE}, for a binary fit, the link applied to the fitted
#' predictor, which is what \code{fit$fitted} and \code{simulate()} hold.
#' \code{type = "link"} returns the posterior mean \eqn{m} of the predictor.
#' A normal fit's response is its predictor plus Gaussian noise: its
#' predictions are the posterior mean of the predictor, with or without
#' \code{plug_in}.
#' Predictions for a design given in \code{newdata} add the coefficients'
#' contribution to the same posterior.
#'
#' @param object an \code{ame_als} object.
#' @param newdata optional dyadic covariate array matching \code{object$X}
#'   (\code{[n_row, n_col, p, T]}); a matrix or 3D array is coerced.
#' @param type \code{"response"} (default) or \code{"link"}.
#' @param plug_in logical; for an empirical-Bayes fit, return the plug-in
#'   values (\code{TRUE}: the link at the fitted predictor, the posterior mean)
#'   or the predictive ones that integrate the link over the uncertainty of the
#'   predictor (\code{FALSE}). The default is \code{TRUE} for a Poisson fit and
#'   \code{FALSE} for a binary fit.
#' @param ... ignored.
#' @return A matrix (cross-sectional) or list of matrices (longitudinal).
#' @method predict ame_als
#' @export
predict.ame_als <- function(
	object,
	newdata = NULL,
	type = c("response", "link"),
	plug_in = NULL,
	...
) {
	type <- match.arg(type)
	plug_in <- .ae_resolve_plug_in(object, plug_in)
	# refuse MCMC-only predict arguments such as `h`, which an ALS fit cannot use
	dots <- list(...)
	dot_names <- names(dots)
	if ("h" %in% dot_names && isTRUE(any(unlist(dots[dot_names == "h"]) > 0))) {
		cli::cli_abort(c(
			"{.fn predict.ame_als} does not support h-step forecasting.",
			"i" = "ALS is a point estimator with no dynamic state to propagate.",
			"i" = "Refit with {.code method = \"mcmc\", dynamic_beta = TRUE} (cross-sectional) or {.fn lame} (longitudinal) for forecasts."
		))
	}
	if (any(c("interval", "by_draw", "newexposure") %in% dot_names)) {
		cli::cli_warn(c(
			"Argument{?s} {.arg {intersect(c('interval','by_draw','newexposure'), dot_names)}} only apply to MCMC fits and are ignored here.",
			"i" = "For uncertainty on an ALS fit, run {.fn ame_als_bootstrap} or pass {.arg bootstrap} to {.fn ame_als}."
		))
	}
	EZ <- object$EZ
	if (!is.null(newdata)) {
		X_old <- object[["X"]]
		p_dyad <- dim(X_old)[3]
		if (p_dyad == 0L) {
			cli::cli_abort("This fit has no dyadic covariates; {.arg newdata} cannot be applied.")
		}
		Xn <- newdata
		row_target <- dimnames(X_old)[[1L]]
		col_target <- dimnames(X_old)[[2L]]
		cov_target <- dimnames(X_old)[[3L]]
		align_covariate_names <- function(xt, got, t) {
			if (is.null(got) || is.null(cov_target)) return(xt)
			if (setequal(got, cov_target)) {
				return(xt[, , cov_target, drop = FALSE])
			}
			got_suffixed <- .lame_apply_suffix(got, "dyad")
			if (setequal(got_suffixed, cov_target)) {
				xt <- xt[, , match(cov_target, got_suffixed), drop = FALSE]
				dimnames(xt)[[3L]] <- cov_target
				return(xt)
			}
			cli::cli_abort(c(
				"{.arg newdata} slice {t} has dyadic-covariate names that do not match the fit.",
				"i" = "Use the same covariate names as the fitted {.arg Xdyad}."
			))
		}
		align_newdata_slice <- function(xt, t) {
			xt <- .ae_dyad_to_3d(xt, nr = dim(X_old)[1L], nc = dim(X_old)[2L])
			dd <- dim(xt)
			if (dd[1L] != dim(X_old)[1L] || dd[2L] != dim(X_old)[2L] ||
					dd[3L] != p_dyad) {
				cli::cli_abort(c(
					"{.arg newdata} slice {t} must match the fitted dyadic design.",
					"i" = "Expected [{paste(dim(X_old)[1:3], collapse = ', ')}], got [{paste(dd, collapse = ', ')}]."
				))
			}
			dn <- dimnames(xt)
			xr <- if (length(dn) >= 1L) dn[[1L]] else NULL
			xc <- if (length(dn) >= 2L) dn[[2L]] else NULL
			xk <- if (length(dn) >= 3L) dn[[3L]] else NULL
			if (!is.null(xr) && !is.null(row_target)) {
				.ae_match_cov_names(xr, row_target, "newdata", t, "row")
				xt <- xt[row_target, , , drop = FALSE]
			}
			if (!is.null(xc) && !is.null(col_target)) {
				.ae_match_cov_names(xc, col_target, "newdata", t, "column")
				xt <- xt[, col_target, , drop = FALSE]
			}
			xt <- align_covariate_names(xt, xk, t)
			dimnames(xt) <- dimnames(X_old)[1:3]
			xt
		}
		if (is.list(Xn) && !is.data.frame(Xn)) {
			if (length(Xn) != object$n_time) {
				cli::cli_abort("{.arg newdata} must have one dyadic covariate entry per time slice ({object$n_time}).")
			}
			Xlist <- lapply(seq_len(object$n_time), function(t) align_newdata_slice(Xn[[t]], t))
			Xarr <- array(NA_real_, dim = dim(X_old), dimnames = dimnames(X_old))
			for (t in seq_len(object$n_time)) Xarr[, , , t] <- Xlist[[t]]
			Xn <- Xarr
		} else {
			if (is.matrix(Xn)) Xn <- array(Xn, c(dim(Xn), 1L, 1L))
			if (length(dim(Xn)) == 3L) Xn <- array(Xn, c(dim(Xn), 1L))
		}
		if (!identical(dim(Xn), dim(X_old))) {
			cli::cli_abort(c(
				"{.arg newdata} must match the fitted dyadic design dimensions.",
				"i" = "Expected [{paste(dim(X_old), collapse = ', ')}], got [{paste(dim(Xn), collapse = ', ')}]."
			))
		}
		dn <- dimnames(Xn)
		if (!is.null(dn[[1L]]) && !is.null(row_target)) {
			.ae_match_cov_names(dn[[1L]], row_target, "newdata", "all", "row")
			Xn <- Xn[row_target, , , , drop = FALSE]
		}
		if (!is.null(dn[[2L]]) && !is.null(col_target)) {
			.ae_match_cov_names(dn[[2L]], col_target, "newdata", "all", "column")
			Xn <- Xn[, col_target, , , drop = FALSE]
		}
		if (!is.null(dn[[3L]]) && !is.null(cov_target)) {
			got <- dn[[3L]]
			if (setequal(got, cov_target)) {
				Xn <- Xn[, , cov_target, , drop = FALSE]
			} else {
				got_suffixed <- .lame_apply_suffix(got, "dyad")
				if (!setequal(got_suffixed, cov_target)) {
					cli::cli_abort(c(
						"{.arg newdata} dyadic-covariate names do not match the fit.",
						"i" = "Use the same covariate names as the fitted {.arg Xdyad}."
					))
				}
				Xn <- Xn[, , match(cov_target, got_suffixed), , drop = FALSE]
				dimnames(Xn)[[3L]] <- cov_target
			}
		}
		beta_dyad <- object$beta[seq_len(p_dyad)]
		EZ <- lapply(seq_len(object$n_time), function(t) {
			xb_old <- xb_new <- matrix(0, object$dims$n_row, object$dims$n_col)
			for (k in seq_len(p_dyad)) {
				xb_old <- xb_old + beta_dyad[k] * X_old[, , k, t]
				xb_new <- xb_new + beta_dyad[k] * Xn[, , k, t]
			}
			ez <- object$EZ[[t]] - xb_old + xb_new
			dimnames(ez) <- dimnames(object$EZ[[t]])
			ez
		})
		names(EZ) <- object$time_names %||% names(object$EZ)
	}
	pred <- if (type == "link") EZ else {
		lat <- .ae_default(object$pois_latent_var, 0)
		lapply(
			EZ,
			function(ez) .ae_link_inverse(
				ez,
				object$family,
				object$link,
				latent_var = lat
			)
		)
	}
	if (!isTRUE(plug_in) && .ae_is_eb_fit(object)) {
		eb_pred <- .ae_eb_predict_response(object, X = if (!is.null(newdata)) Xn, type = type)
		if (!is.null(eb_pred)) return(eb_pred)
	}
	if (isTRUE(object$longitudinal)) pred else pred[[1]]
}

#' Plot the convergence of a fast AME fit
#'
#' @description
#' Plots the block coordinate descent deviance/SSE history -- a convergence
#' diagnostic for an \code{\link{ame_als}} or \code{\link{lame_als}}
#' fit. For substantive plots see \code{\link{uv_plot}}, \code{\link{ab_plot}}
#' and \code{\link{gof_plot}}.
#'
#' @param x an \code{ame_als} object.
#' @param ... further arguments passed to the underlying plot.
#' @return \code{x}, invisibly.
#' @method plot ame_als
#' @export
plot.ame_als <- function(x, ...) {
	dh <- x$dev_history
	if (is.null(dh) || length(dh) < 2L) {
		cli::cli_abort(c(
			"No convergence history is stored for this fit.",
			"i" = "Use {.fn uv_plot}, {.fn ab_plot} or {.fn gof_plot} to visualise the fit."
		))
	}
	graphics::plot.default(
		seq_along(dh),
		dh,
		type = "b",
		pch = 19,
		xlab = "Iteration",
		ylab = "Deviance or residual sum of squares",
		main = "Convergence of the fast estimator",
		...
	)
	invisible(x)
}

#' Log-likelihood is not defined for a fast AME fit
#'
#' The fast estimator minimises a (working-response) least-squares objective,
#' not a family likelihood, so it has no log-likelihood and \code{AIC} / \code{BIC}
#' are undefined. Calling \code{logLik()} raises an informative error rather
#' than returning a misleading number.
#'
#' An empirical-Bayes fit (the default for normal, binary and Poisson
#' responses) carries a Laplace approximation, at its posterior
#' mode, to its log marginal likelihood, \code{fit$eb$log_marginal} (see
#' \strong{Choosing the rank} in \code{\link{ame_als}}), which compares fits of
#' different latent ranks.
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return Never returns; raises an error.
#' @method logLik ame_als
#' @export
logLik.ame_als <- function(object, ...) {
	bullets <- c(
		"A log-likelihood is not defined for the fast {.fn ame_als} estimator.",
		"i" = "It minimises a working-response least-squares objective, so {.fn AIC} / {.fn BIC} do not apply."
	)
	if (!is.null(object$eb$log_marginal) && is.finite(object$eb$log_marginal)) {
		bullets <- c(
			bullets,
			"i" = "To compare fits of different latent ranks use the Laplace log marginal likelihood (at the posterior mode) in {.code fit$eb$log_marginal}."
		)
	}
	cli::cli_abort(c(
		bullets,
		"i" = "For likelihood-based model comparison use the MCMC estimator {.fn ame} / {.fn lame}."
	))
}

#' Print an ame_als object
#'
#' @param x an \code{ame_als} object.
#' @param digits number of digits to display.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @method print ame_als
#' @export
print.ame_als <- function(x, digits = 4, ...) {
	.lame_cli_stdout({
		cli::cli_text("{.strong AME fit via iterative block coordinate descent} (fast, MCMC-free)")
		cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}}{ifelse(x$symmetric, ' (symmetric)', '')} | R = {x$R}")
		cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slice{?s}")
		cli::cli_text("Converged: {.val {x$converged}} in {x$iterations} iteration{?s}")
		.ae_print_retirement(x$eb$retirement)
		.ae_print_bias_correct(x$bias_correct, x$family, x$eb$estimator)
		.ae_print_information(x$eb$variational$information, x$eb$product_withheld)
		# note that lame_als() is static even on a panel
		if (!is.null(x$n_time) && x$n_time > 1L) {
			cli::cli_alert_warning(c(
				"This is a {.strong STATIC} fit pooled across {x$n_time} time periods: ",
				"{.code U}, {.code V}, {.code a}, {.code b} do not vary by time. ",
				"For a dynamic point fit, call {.fn lame} with {.code method = \"als\"} ",
				"and dynamic arguments such as {.code dynamic_uv = TRUE}, ",
				"{.code dynamic_ab = TRUE}, or {.code dynamic_beta = ...}. ",
				"Use {.code method = \"mcmc\"} when you need posterior draws or an unsupported family."
			))
		}
		cli::cli_text("")
		cli::cli_text("{.strong Coefficients:}")
		cli::cli_verbatim(utils::capture.output(print(round(x$coefficients, digits))))
		cli::cli_text("")
		cli::cli_text("{.strong Variance components:}")
		cli::cli_verbatim(utils::capture.output(print(round(x$VC, digits))))
		rho_note <- .ae_rho_not_estimated_text(x)
		if (!is.null(rho_note)) cli::cli_text("{rho_note}")
		cli::cli_text("")
		cli::cli_text(
			cli::col_grey("va/vb: variances of the sender/receiver effects | cab: their covariance | rho: dyadic residual reciprocity | ve: residual variance.")
		)
		if (identical(x$bias_correct$method, "eb")) {
			cli::cli_text(
				cli::col_grey("va, cab and vb are the fitted covariance of the actor effects (the model's random-effect variances). The effects a and b are not re-centred; the intercept is the population intercept.")
			)
		} else {
			cli::cli_text(cli::col_grey("(Descriptive summaries, not random-effect variance components.)"))
		}
		if (!is.null(x$eb$se)) {
			cli::cli_text("")
			cli::cli_text("{.strong Laplace standard errors} (penalised Fisher information, empirical-Bayes actor-effect prior):")
			cli::cli_verbatim(utils::capture.output(print(round(x$eb$se, digits))))
		}
		if (isTRUE(x$meta$uncertainty_available)) {
			nm <- x$meta$bootstrap_type %||% "parametric"
			Rb <- x$meta$bootstrap_R %||% NA_integer_
			cli::cli_text(cli::col_grey(
				"Bootstrap intervals available ({nm}, R = {Rb}); use {.fn confint} for CIs."
			))
		} else if (is.null(x$eb$se)) {
			cli::cli_text(
				cli::col_grey("Point estimate only; no uncertainty. Re-fit with {.code bootstrap = N} (e.g. {.code ame_als(Y, ..., bootstrap = 200)}) to attach bootstrap intervals in one call.")
			)
		}
	})
	invisible(x)
}

# print note on a withheld or weakly determined product
.ae_print_information <- function(info, withheld = NULL) {
	# a product the data do not determine is withheld (.ae_withhold_product)
	if (!is.null(withheld)) {
		cli::cli_text(.ae_product_withheld_text(withheld))
		return(invisible(NULL))
	}
	# normal: signal of the product in noise variances per loading (.ae_gauss_information)
	if (!is.null(info) && "signal" %in% names(info) && !"signal_per_loading" %in% names(info)) {
		if (!is.finite(info[["per_loading"]]) || info[["per_loading"]] >= info[["threshold"]]) return(invisible(NULL))
		cli::cli_text("Information about the latent factors: the fitted product has {signif(info[['per_loading']], 2)} noise variances of signal per loading ({info[['loadings']]} loadings); the data barely determine it at this rank, and a smaller rank may be preferred (compare {.code fit$eb$log_marginal} across ranks).")
		return(invisible(NULL))
	}
	# binary and count fits withhold a weak product, so there is nothing more to print
	invisible(NULL)
}

# one line on the bias correction of a static ALS fit
.ae_print_bias_correct <- function(bc, family, estimator = NULL) {
	if (is.null(bc)) return(invisible(NULL))
	if (identical(bc$target, "latent count variance and reciprocal covariance") &&
			isTRUE(bc$applied)) {
		cli::cli_text("Bias correction: split-dyad count variance and reciprocity (uncorrected summaries in {.code $raw})")
	} else if (identical(bc$method, "bootstrap") && !isFALSE(bc$applied)) {
		cli::cli_text("Bias correction: parametric bootstrap{if (isTRUE(bc$stored)) ' (the correction of the fit this one was refit from)' else ''} (uncorrected values in {.code $raw})")
	} else if (identical(bc$method, "eb")) {
		how <- if (identical(estimator, "variational")) {
			"variational posterior of the effects and factors, coefficients and scales by empirical Bayes"
		} else "Laplace-marginal coefficients"
		cells <- if (identical(family, "poisson")) "Poisson-lognormal counts with the variance of the noise by marginal likelihood, " else
		if (identical(family, "normal")) "exact Gaussian cell with the variance of the noise by EM, " else ""
		cli::cli_text("Estimator: empirical Bayes ({cells}Gaussian prior on the actor effects, covariance by EM, {how}); no bootstrap or jackknife correction is applied")
	} else if (identical(bc$method, "jackknife") && !isFALSE(bc$applied)) {
		split <- bc$split %||% "checkerboard"
		cli::cli_text("Bias correction: split-sample jackknife, {split} split (uncorrected values in {.code $raw})")
	} else if (!is.null(bc$basis)) {
		cli::cli_text("Bias correction: none needed at this size ({bc$basis})")
	} else if (identical(bc$requested, "none") && identical(family, "binary")) {
		cli::cli_text("Bias correction: none ({.code bias_correct = \"none\"} was requested)")
	} else if (identical(bc$requested, "jackknife")) {
		reason <- bc$reason %||% "the split was refused"
		cli::cli_alert_warning("Bias correction not applied: these {family} ALS estimates are uncorrected; {reason}.")
	}
	invisible(NULL)
}

#' Summarize an ame_als object
#'
#' @param object an \code{ame_als} object.
#' @param ... ignored.
#' @return A list of class \code{"summary.ame_als"}, printed as a table.
#' @method summary ame_als
#' @export
summary.ame_als <- function(object, ...) {
	# coefficient table, with bootstrap se and interval columns when attached
	bs <- object$bootstrap
	have_bs <- isTRUE(object$meta$uncertainty_available) &&
		!is.null(bs) && !is.null(bs$param_names) &&
		!is.null(bs$point_est)
	if (!have_bs && !is.null(object$eb$se)) {
		se <- object$eb$se[names(object$coefficients)]
		tab <- cbind(
			Estimate = object$coefficients,
			Std.Error = se,
			`CI 2.5%` = object$coefficients - 1.96 * se,
			`CI 97.5%` = object$coefficients + 1.96 * se
		)
		rownames(tab) <- names(object$coefficients)
	} else if (have_bs) {
		# align bootstrap columns with the coefficients by name
		coef_nms <- names(object$coefficients)
		match_idx <- match(coef_nms, bs$param_names)
		tab <- cbind(
			Estimate = object$coefficients,
			Std.Error = bs$se[match_idx],
			`CI 2.5%` = bs$ci_lo[match_idx],
			`CI 97.5%` = bs$ci_hi[match_idx]
		)
		rownames(tab) <- coef_nms
	} else {
		tab <- cbind(Estimate = object$coefficients)
	}
	out <- list(
		call = object$call,
		family = object$family,
		mode = object$mode,
		symmetric = object$symmetric,
		R = object$R,
		n_time = object$n_time,
		dims = object$dims,
		coefficients = tab,
		VC = object$VC,
		rho_note = .ae_rho_not_estimated_text(object),
		converged = object$converged,
		iterations = object$iterations,
		deviance = object$deviance,
		objective = object$objective,
		bias_correct = object$bias_correct,
		estimator = object$eb$estimator,
		information = object$eb$variational$information,
		product_withheld = object$eb$product_withheld,
		retirement = object$eb$retirement,
		approximation_note = object$meta$approximation_note,
		uncertainty_available = isTRUE(object$meta$uncertainty_available),
		uncertainty_source = object$meta$uncertainty_source,
		bootstrap_R = object$meta$bootstrap_R,
		bootstrap_type = object$meta$bootstrap_type,
		bootstrap_withheld = if (have_bs) bs$withheld
	)
	class(out) <- "summary.ame_als"
	out
}

#' Print a fast AME summary
#'
#' @param x a \code{summary.ame_als} object.
#' @param digits number of digits to display.
#' @param ... ignored.
#' @return \code{x}, invisibly.
#' @method print summary.ame_als
#' @export
print.summary.ame_als <- function(x, digits = 4, ...) {
	.lame_cli_stdout({
		cli::cli_h2("AME exploration summary")
		cli::cli_text("Algorithm: iterative block coordinate descent (Hoff & Minhas SIR estimator)")
		cli::cli_text("Family: {.field {x$family}} | Mode: {.field {x$mode}}{ifelse(x$symmetric, ' (symmetric)', '')} | R = {x$R}")
		cli::cli_text("Network: {x$dims$n_row} x {x$dims$n_col}, {x$n_time} time slice{?s}, {x$dims$p} covariate{?s}")
		cli::cli_text("Converged: {.val {x$converged}} ({x$iterations} iteration{?s}); deviance = {round(x$deviance, digits)}{if (!is.null(x$objective) && is.finite(x$objective) && abs(x$objective - x$deviance) > 1e-8 * max(1, abs(x$deviance))) paste0(' (optimiser objective ', round(x$objective, digits), ')') else ''}")
		.ae_print_retirement(x$retirement)
		.ae_print_bias_correct(x$bias_correct, x$family, x$estimator)
		.ae_print_information(x$information, x$product_withheld)
		if (!is.null(x$n_time) && x$n_time > 1L) {
			cli::cli_alert_warning(c(
				"This is a {.strong STATIC} fit pooled across {x$n_time} time periods: ",
				"{.code U}, {.code V}, {.code a}, {.code b} do not vary by time. ",
				"For a dynamic point fit, call {.fn lame} with {.code method = \"als\"} ",
				"and dynamic arguments such as {.code dynamic_uv = TRUE}, ",
				"{.code dynamic_ab = TRUE}, or {.code dynamic_beta = ...}. ",
				"Use {.code method = \"mcmc\"} when you need posterior draws or an unsupported family."
			))
		}
		cli::cli_rule()
		cli::cli_text("{.strong Regression coefficients}")
		cli::cli_verbatim(utils::capture.output(print(round(x$coefficients, digits))))
		cli::cli_text("")
		cli::cli_text("{.strong Variance components}")
		cli::cli_verbatim(utils::capture.output(print(round(x$VC, digits))))
		if (!is.null(x$rho_note)) cli::cli_text("{x$rho_note}")
		cli::cli_text(
			cli::col_grey("va/vb: variances of the sender/receiver effects | cab: their covariance | rho: dyadic residual reciprocity | ve: residual variance.")
		)
		if (identical(x$bias_correct$method, "eb")) {
			cli::cli_text(
				cli::col_grey("va, cab and vb are the fitted covariance of the actor effects (the model's random-effect variances). The effects a and b are not re-centred; the intercept is the population intercept.")
			)
		} else {
			cli::cli_text(cli::col_grey("(Descriptive summaries, not random-effect variance components.)"))
		}
		cli::cli_rule()
		cli::cli_text(cli::col_grey(x$approximation_note))
		if (!isTRUE(x$uncertainty_available) && "Std.Error" %in% colnames(x$coefficients)) {
			cli::cli_text(
				cli::col_grey("Standard errors: Laplace approximation (penalised Fisher information), conditional on the variance components; intervals are Wald +/- 1.96 se.")
			)
		} else if (isTRUE(x$uncertainty_available)) {
			nm <- x$uncertainty_source %||% "bootstrap"
			Rb <- x$bootstrap_R %||% NA_integer_
			tp <- x$bootstrap_type %||% "parametric"
			cli::cli_text(cli::col_grey(
				"Uncertainty available: {nm} ({tp}, R = {Rb}). Use {.fn confint} for intervals."
			))
			.ae_print_withheld(list(withheld = x$bootstrap_withheld, type = x$bootstrap_type))
		} else {
			cli::cli_text(
				cli::col_grey("No standard errors: this is a point estimate. Re-fit with {.code bootstrap = N} (e.g. {.code ame_als(Y, ..., bootstrap = 200)}) to attach bootstrap intervals.")
			)
		}
	})
	invisible(x)
}
