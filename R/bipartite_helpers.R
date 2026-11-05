#' Bipartite network helper functions
#'
#' @name bipartite_helpers
#' @keywords internal

#' @description Initialize bipartite starting values
#' @return List of starting values for bipartite MCMC
#' @keywords internal
init_bipartite_startvals <- function(
	Y,
	family,
	nA,
	nB,
	RA,
	RB,
	Tn,
	Xlist = NULL,
	odmax = NULL,
	unit_var = 1
) {

	# seed Z; the rank families (cbin, frn) need a seed that respects their bounds, since the row-cone
	# samplers condition on bounds built from the surrounding Z and an all-NA seed never moves
	Z <- array(dim = c(nA, nB, Tn))

	for (t in 1:Tn) {
		Yt <- Y[, , t]

		if (family == "normal") {
			Z[, , t] <- Yt
		} else if (family == "binary") {
			Z[, , t] <- matrix(qnorm((Yt + 0.5) / 2), nA, nB)

			if (sum(Yt == 0, na.rm = TRUE) > 0 && sum(Yt == 1, na.rm = TRUE) > 0) {
				z0_max <- max(Z[, , t][Yt == 0], na.rm = TRUE)
				z1_min <- min(Z[, , t][Yt == 1], na.rm = TRUE)
				z_split <- (z0_max + z1_min) / 2
				Z[, , t] <- Z[, , t] - z_split
			}
		} else if (family == "poisson") {
			# log(y + 1) on observed cells; missing cells are drawn around the median so the first sweep is finite
			Z[, , t] <- log(Yt + 1)
			miss <- is.na(Z[, , t])
			if (any(miss)) {
				slice <- Z[, , t]
				slice[miss] <- stats::rnorm(
					sum(miss),
					mean = stats::median(slice, na.rm = TRUE),
					sd = 1
				)
				Z[, , t] <- slice
			}
		} else if (family == "ordinal") {
			# rank-normal start on observed cells; missing cells get a standard normal draw
			Z[, , t] <- matrix(
				qnorm((Yt + 0.5) / (max(Yt, na.rm = TRUE) + 1)),
				nA,
				nB
			)
			miss <- is.na(Z[, , t])
			if (any(miss)) Z[, , t][miss] <- stats::rnorm(sum(miss), 0, 1)
		} else if (family == "cbin") {
			# y = 1 -> z > 0, y = 0 -> z < 0, NA unconstrained
			slice <- matrix(stats::rnorm(nA * nB), nA, nB)
			pos <- !is.na(Yt) & Yt == 1
			neg <- !is.na(Yt) & Yt == 0
			slice[pos] <- abs(slice[pos])
			slice[neg] <- -abs(slice[neg])
			Z[, , t] <- slice
		} else if (family == "frn") {
			# nominated cells start positive in rank order, observed non-nominations negative,
			# missing cells standard normal
			slice <- matrix(stats::rnorm(nA * nB), nA, nB)
			for (i in seq_len(nA)) {
				not_nominated <- which(Yt[i, ] == 0)
				slice[i, not_nominated] <- -abs(slice[i, not_nominated])
				ranked <- which(Yt[i, ] > 0)
				if (length(ranked)) {
					by_rank <- ranked[order(Yt[i, ranked], decreasing = TRUE)]
					slice[i, by_rank] <- sort(
						abs(stats::rnorm(length(ranked))),
						decreasing = TRUE
					)
				}
			}
			Z[, , t] <- slice
		}
	}

	# starting effects, positions and G
	a <- rep(0, nA)
	b <- rep(0, nB)

	if (RA > 0) {
		U <- array(rnorm(nA * RA * Tn, 0, 0.1), c(nA, RA, Tn))
	} else {
		U <- array(0, c(nA, 1, Tn))
	}

	if (RB > 0) {
		V <- array(rnorm(nB * RB * Tn, 0, 0.1), c(nB, RB, Tn))
	} else {
		V <- array(0, c(nB, 1, Tn))
	}

	if (RA > 0 && RB > 0) {
		min_rank <- min(RA, RB)
		G <- matrix(0, RA, RB)
		diag(G)[1:min_rank] <- sqrt(unit_var)
	} else {
		G <- matrix(0, max(1, RA), max(1, RB))
	}

	# coefficients and variances
	if (!is.null(Xlist)) {
		p <- dim(Xlist[[1]])[3]
		beta <- rep(0, p)
	} else {
		beta <- numeric(0)
	}

	# the starting variances are written in the unit of the response (1 where the link fixes the scale)
	sigma2_a <- unit_var
	sigma2_b <- unit_var
	s2 <- unit_var
	rho <- NULL

	list(
		Z = Z,
		a = a,
		b = b,
		U = U,
		V = V,
		G = G,
		beta = beta,
		sigma2_a = sigma2_a,
		sigma2_b = sigma2_b,
		s2 = s2,
		rho = rho
	)
}

# overlay user-supplied start values on the bipartite initialiser's
# output. a slot is taken only when its shape matches the panel; the
# initialiser's value is kept for the others and they are reported once
.lame_bip_overlay_start <- function(
	sv,
	start_vals,
	nA,
	nB,
	RA,
	RB,
	N,
	dynamic_ab = FALSE
) {
	if (!is.list(start_vals) || length(start_vals) == 0L) return(sv)
	p <- length(sv$beta)
	skipped <- character(0)
	dims_ok <- function(x, d) {
		is.numeric(x) && !is.null(dim(x)) &&
			identical(as.integer(dim(x)), as.integer(d))
	}
	vec_ok <- function(x, len) is.numeric(x) && is.null(dim(x)) && length(x) == len
	scalar_ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
	take <- function(nm, ok) {
		x <- start_vals[[nm]]
		if (is.null(x)) return(invisible(NULL))
		if (isTRUE(ok(x))) sv[[nm]] <<- x else skipped <<- c(skipped, nm)
		invisible(NULL)
	}
	take("Z", function(x) dims_ok(x, c(nA, nB, N)))
	take("beta", function(x) vec_ok(x, p))
	# additive effects: one value per actor, or one per actor and period
	# for a dynamic fit; a static fit takes the actor means of a per-period
	# start
	for (side in list(list(nm = "a", n = nA), list(nm = "b", n = nB))) {
		x <- start_vals[[side$nm]]
		if (is.null(x)) next
		if (dims_ok(x, c(side$n, N)) && !dynamic_ab) x <- rowMeans(x)
		if (vec_ok(x, side$n) || (dynamic_ab && dims_ok(x, c(side$n, N)))) {
			sv[[side$nm]] <- x
		} else {
			skipped <- c(skipped, side$nm)
		}
	}
	# latent positions: the initialiser carries one slice per period, so a
	# static matrix is replicated across periods
	for (side in list(
		list(nm = "U", n = nA, r = RA),
		list(nm = "V", n = nB, r = RB)
	)) {
		x <- start_vals[[side$nm]]
		if (is.null(x) || side$r == 0L) next
		if (dims_ok(x, c(side$n, side$r))) x <- array(x, c(side$n, side$r, N))
		if (dims_ok(x, c(side$n, side$r, N))) {
			sv[[side$nm]] <- x
		} else {
			skipped <- c(skipped, side$nm)
		}
	}
	g <- start_vals[["G"]]
	if (!is.null(g) && RA > 0L && RB > 0L) {
		# a per-period G enters as its period mean
		if (dims_ok(g, c(RA, RB, N))) g <- apply(g, c(1L, 2L), mean)
		if (dims_ok(g, c(RA, RB)) && all(is.finite(g))) {
			sv$G <- g
		} else {
			skipped <- c(skipped, "G")
		}
	}
	take("s2", function(x) scalar_ok(x) && x > 0)
	Sab <- start_vals[["Sab"]]
	if (!is.null(Sab)) {
		if (dims_ok(Sab, c(2L, 2L)) && all(is.finite(diag(Sab))) &&
				all(diag(Sab) > 0)) {
			sv$sigma2_a <- Sab[1L, 1L]
			sv$sigma2_b <- Sab[2L, 2L]
		} else {
			skipped <- c(skipped, "Sab")
		}
	}
	for (nm in c("rho_uv", "sigma_uv", "rho_ab", "sigma_ab", "level_move_sd")) take(nm, scalar_ok)
	if (length(skipped) > 0L) {
		cli::cli_warn(c(
			"{.arg start_vals} slot{?s} {.val {skipped}} do{?es/} not match this bipartite panel and {?is/are} ignored.",
			"i" = "Expected {.code Z} {nA} x {nB} x {N}, {.code beta} of length {p}, {.code a} of length {nA}, {.code b} of length {nB}, {.code U} {nA} x {RA}, {.code V} {nB} x {RB}, {.code G} {RA} x {RB}, {.code Sab} 2 x 2 and {.code s2} a positive scalar."
		))
	}
	sv
}

#' Convert bipartite list data to array format
#' @return List with components \code{Y}, \code{Xdyad}, \code{Xrow}, \code{Xcol} in array format
#' @keywords internal
list_to_array_bipartite <- function(
	rowActorSet,
	colActorSet,
	Y_list,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL
) {
	nA <- length(rowActorSet)
	nB <- length(colActorSet)
	Tn <- length(Y_list)

	# place Y on the sorted actor sets
	Y_array <- array(
		NA,
		c(nA, nB, Tn),
		dimnames = list(rowActorSet, colActorSet, names(Y_list))
	)

	for (t in 1:Tn) {
		Yt <- Y_list[[t]]
		row_idx <- match(rownames(Yt), rowActorSet)
		col_idx <- match(colnames(Yt), colActorSet)

		for (i in seq_along(row_idx)) {
			for (j in seq_along(col_idx)) {
				if (!is.na(row_idx[i]) && !is.na(col_idx[j])) {
					Y_array[row_idx[i], col_idx[j], t] <- Yt[i, j]
				}
			}
		}
	}

	# suffix rule shared with the other covariate-naming sites (R/globals.R)
	.bip_apply_suffix <- .lame_apply_suffix

	# every covariate below is placed by matching its dimnames against the
	# sorted actor sets, so a slice without names cannot be placed and stops
	# the fit here
	.bip_need_names <- function(nms, t, what, axis, expect) {
		if (!is.null(nms)) return(invisible(NULL))
		cli::cli_abort(c(
			"{.arg {what}} period {t} has no {axis} names.",
			"i" = "Bipartite covariates are aligned to the actor sets by name; a slice without names cannot be placed.",
			"i" = "Set {.code {axis}names()} to the {expect} actor labels used in {.arg Y}."
		))
	}

	# place Xdyad on the sorted actor sets
	if (!is.null(Xdyad)) {
		Xdyad <- .coerce_Xdyad_3d(Xdyad)
		Xdyad_out <- vector("list", Tn)
		for (t in 1:Tn) {
			Xt <- Xdyad[[t]]
			# the slot stays NULL, keeping one element per period
			if (is.null(Xt)) next
			pd_t <- dim(Xt)[3]
			# `_dyad` suffix on slice names, as in the unipartite design (R/design_array.R)
			slice_nms <- dimnames(Xt)[[3]]
			if (is.null(slice_nms)) slice_nms <- paste0("Xdyad", seq_len(pd_t))
			slice_nms <- .bip_apply_suffix(slice_nms, "dyad")
			Xt_new <- array(
				NA,
				dim = c(nA, nB, pd_t),
				dimnames = list(rowActorSet, colActorSet, slice_nms)
			)
			.bip_need_names(rownames(Xt), t, "Xdyad", "row", "row")
			.bip_need_names(colnames(Xt), t, "Xdyad", "col", "column")
			row_idx <- match(rownames(Xt), rowActorSet)
			col_idx <- match(colnames(Xt), colActorSet)
			for (i in seq_along(row_idx)) {
				for (j in seq_along(col_idx)) {
					if (!is.na(row_idx[i]) && !is.na(col_idx[j])) {
						Xt_new[row_idx[i], col_idx[j], ] <- Xt[i, j, ]
					}
				}
			}
			Xdyad_out[[t]] <- Xt_new
		}
		Xdyad <- Xdyad_out
	}

	# place Xrow on the sorted row actors as a 3-d array
	if (!is.null(Xrow)) {
		pr <- ncol(Xrow[[1]])
		row_nms <- .bip_apply_suffix(
			if (is.null(colnames(Xrow[[1]]))) paste0("Xrow", seq_len(pr))
			else colnames(Xrow[[1]]),
			"row"
		)
		Xrow_array <- array(
			NA,
			dim = c(nA, pr, Tn),
			dimnames = list(rowActorSet, row_nms, names(Y_list))
		)
		for (t in 1:Tn) {
			Xrt <- Xrow[[t]]
			if (is.null(Xrt)) next
			.bip_need_names(rownames(Xrt), t, "Xrow", "row", "row")
			row_idx <- match(rownames(Xrt), rowActorSet)
			for (i in seq_along(row_idx)) {
				if (!is.na(row_idx[i])) Xrow_array[row_idx[i], , t] <- Xrt[i, ]
			}
		}
		Xrow <- Xrow_array
	}

	# place Xcol on the sorted column actors as a 3-d array
	if (!is.null(Xcol)) {
		pc <- ncol(Xcol[[1]])
		col_nms <- .bip_apply_suffix(
			if (is.null(colnames(Xcol[[1]]))) paste0("Xcol", seq_len(pc))
			else colnames(Xcol[[1]]),
			"col"
		)
		Xcol_array <- array(
			NA,
			dim = c(nB, pc, Tn),
			dimnames = list(colActorSet, col_nms, names(Y_list))
		)
		for (t in 1:Tn) {
			Xct <- Xcol[[t]]
			if (is.null(Xct)) next
			.bip_need_names(rownames(Xct), t, "Xcol", "row", "column")
			col_idx <- match(rownames(Xct), colActorSet)
			for (i in seq_along(col_idx)) {
				if (!is.na(col_idx[i])) Xcol_array[col_idx[i], , t] <- Xct[i, ]
			}
		}
		Xcol <- Xcol_array
	}

	list(Y = Y_array, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol)
}
