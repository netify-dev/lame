#' Draw latent Z for censored-binary nomination data
#'
#' Gibbs refresh of the latent normal matrix \code{Z} for the
#' censored-binary (\code{cbin}) family, where each actor emits a capped
#' set of binary nominations. Within every sending row the update honours
#' three order constraints: a nominated tie has a positive latent value,
#' every nominated tie dominates every non-nominated tie, and a
#' non-nominated tie in a row that has not spent its full nomination
#' budget (\code{odobs < odmax}) is held below zero. Reciprocal cells are
#' coupled through the dyadic correlation \code{rho}; the self-loop
#' diagonal and any missing cells are redrawn from their unconstrained
#' conditionals.
#'
#' @param Z current latent matrix (square).
#' @param EZ conditional mean of \code{Z}.
#' @param rho dyadic (within-dyad) correlation.
#' @param Y observed nomination matrix, entries in \code{{0, 1}} with
#'   \code{NA} for unobserved cells.
#' @param odmax nomination cap, a scalar or length-\code{nrow(Z)} vector.
#' @param odobs per-actor count of nominations the actor has made: the
#'   observed nominations plus the hidden cells whose current latent value
#'   is positive (the cells of an actor absent from the slice never count;
#'   see \code{.lame_nom_room}). A row whose count falls short of \code{odmax} has unused
#'   nomination slots, so its non-nominated cells are held below zero.
#' @return the updated latent matrix, same shape and dimnames as \code{Z}.
#' @keywords internal
#' @author lame authors
rZ_cbin_fc <- function(Z, EZ, rho, Y, odmax, odobs) {
	Z <- .lame_nom_feasible(Z, Y, NULL, odmax, odobs, FALSE)
	# rows with hidden cells: label the hidden cells by the current latent values,
	# draw the rows given those labels, then redraw the hidden cells given the
	# observed cells alone (see R/rZ_nomination_hidden.R)
	hid <- .lame_nom_hidden_cells(Y, FALSE)
	if (any(hid)) {
		Ya <- .lame_cbin_implied_labels(Z, Y, odmax, hid)
		Z <- .rZ_cbin_sweep(Z, EZ, rho, Ya, odmax, rowSums(!is.na(Ya) & Ya > 0))
		return(.lame_cbin_hidden_moves(Z, EZ, rho, Y, odmax, hid, FALSE))
	}
	.rZ_cbin_sweep(Z, EZ, rho, Y, odmax, odobs)
}
# the complete-data draw of rZ_cbin_fc: every observed cell of a row from
# its truncated conditional given the row's nominations
.rZ_cbin_sweep <- function(Z, EZ, rho, Y, odmax, odobs) {
	n <- nrow(Z)
	cond_sd <- sqrt(1 - rho^2)

	# unobserved cells carry a sentinel class so class masks never index with NA
	Yc <- Y
	Yc[is.na(Yc)] <- -1L

	# row and column index of every cell, to split the two off-diagonal triangles
	# and look up per-row bounds
	ri <- row(Z)
	ci <- col(Z)
	tri_hi <- which(ri < ci) # (i, j), i < j
	tri_lo <- which(ri > ci) # (i, j), i > j

	# reciprocal linear position of a set of cells: (i, j) -> (j, i)
	mirror <- function(pos) ci[pos] + (ri[pos] - 1L) * n

	# rows with unused budget: their non-nominated ties are capped at zero rather
	# than at the weakest nominee
	has_room <- rep(FALSE, n)
	has_room[which(odobs < odmax)] <- TRUE

	# draw cells from N(ez, cond_sd^2) truncated to (lo, hi) on the log.p scale,
	# which stays finite deep in a tail; a non-finite draw keeps the current value
	refresh <- function(pos, lo, hi) {
		if (!length(pos)) return(invisible(NULL))
		mp <- mirror(pos)
		ez <- EZ[pos] + rho * (Z[mp] - EZ[mp])
		val <- cond_sd * rtnorm_interval_logp(
			ez / cond_sd,
			lo / cond_sd,
			hi / cond_sd
		)
		bad <- !is.finite(val)
		if (any(bad)) val[bad] <- Z[pos][bad]
		Z[pos] <<- val
	}

	# sweep the two nomination classes in random order; each class's bounds read
	# the current Z of the other
	for (cls in sample(c(0L, 1L))) {
		if (cls == 1L) {
			# nominee floor: positive and above the strongest non-nominee in the row
			zt <- Z
			zt[Yc != 0L] <- -Inf
			lo_row <- pmax(0, apply(zt, 1L, max))
			hi_row <- rep(Inf, n)
		} else {
			# non-nominee ceiling: below the weakest nominee in the row, and below zero
			# for rows with unused budget
			zt <- Z
			zt[Yc != 1L] <- Inf
			hi_row <- apply(zt, 1L, min)
			hi_row[has_room] <- 0
			lo_row <- rep(-Inf, n)
		}

		# upper triangle first, then lower, so the second block sees the updated
		# reciprocal values (a sequential Gibbs pass)
		hit_hi <- tri_hi[Yc[tri_hi] == cls]
		refresh(hit_hi, lo_row[ri[hit_hi]], hi_row[ri[hit_hi]])

		hit_lo <- tri_lo[Yc[tri_lo] == cls]
		refresh(hit_lo, lo_row[ri[hit_lo]], hi_row[ri[hit_lo]])
	}

	# self-loops carry the marginal dyad variance; missing cells are drawn from the
	# unit-variance prior (this second write sets any NA diagonal)
	diag(Z) <- stats::rnorm(n, diag(EZ), sqrt(1 + rho))
	gap <- is.na(Y)
	if (any(gap)) Z[gap] <- stats::rnorm(sum(gap), EZ[gap], 1)

	# free cells (absent actors and any cell no nomination can occupy) are
	# partners in their dyads: under a dyadic correlation redraw them given the
	# partner cell, upper triangle first, so each pair keeps its joint law
	if (rho != 0 && any(gap)) {
		off <- gap
		diag(off) <- FALSE
		hit_hi <- which(off & ri < ci)
		hit_lo <- which(off & ri > ci)
		refresh(hit_hi, rep(-Inf, length(hit_hi)), rep(Inf, length(hit_hi)))
		refresh(hit_lo, rep(-Inf, length(hit_lo)), rep(Inf, length(hit_lo)))
	}

	Z
}
