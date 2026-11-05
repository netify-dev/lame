#' Simulate a and Sab from full conditional distributions under frn likelihood
#'
#' Simulate a and Sab from full conditional distributions under frn likelihood
#'
#'
#' @usage raSab_frn_fc(Z, Y, YL, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL,
#' SS=round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
#' @param YL list of ranked individuals, from least to most preferred in each
#' row
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each individual
#' @param odobs per-actor count of nominations the actor has made: the
#'   observed nominations plus the hidden cells whose current latent value
#'   is positive (the cells of an actor absent from the slice never count;
#'   see \code{.lame_nom_room}). A row whose count falls short of \code{odmax} has unused
#'   nomination slots, so its non-nominated cells are held below zero.
#' @param Sab0 prior (inverse) scale matrix for the prior distribution
#' @param eta0 prior degrees of freedom for the prior distribution
#' @param SS number of iterations
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @keywords internal
#' @author lame authors
#' @export raSab_frn_fc
raSab_frn_fc <-
	function(Z, Y, YL, a, b, Sab, odmax, odobs, Sab0 = NULL, eta0 = NULL, SS = round(sqrt(nrow(Z)))) {

		n <- nrow(Z)
		if (is.null(Sab0)) { Sab0 <- diag(2) }
		if (is.null(eta0)) { eta0 <- 4 }

		# row-effect-free residual R_ij = Z_ij - a_i, held fixed while a is refreshed
		resid <- Z - outer(a, rep(1, n))
		rows <- seq_len(n)

		# lower bound for a_i: the least-preferred nominee YL[, 1] keeps Z_ij >= 0
		# rows without nominations are unconstrained
		lo <- -resid[cbind(rows, YL[, 1])]
		lo[is.na(lo)] <- -Inf
		# rows with hidden cells are bounded by the nominations their current latent
		# values imply (see R/rZ_nomination_hidden.R)
		hid <- .lame_nom_hidden_cells(Y, FALSE)
		if (any(hid)) {
			Ya <- .lame_cbin_implied_labels(Z, Y, odmax, hid)
			nom <- !is.na(Ya) & Ya > 0
			lo_m <- -resid; lo_m[!nom] <- -Inf
			lo <- apply(lo_m, 1, max)
			Y <- 1 * (Ya > 0)
			odobs <- rowSums(nom)
		}

		# upper bound for a_i: every observed non-nominee keeps Z_ij <= 0
		drop <- is.na(Y) | (Y != 0)
		cand <- resid
		cand[drop] <- -Inf
		hi <- -apply(cand, 1, max, na.rm = TRUE)
		## a saturated out-degree (odobs >= odmax) removes the ceiling
		hi[odobs >= odmax] <- Inf

		for (iter in seq_len(SS)) {

			# a | b ~ N((S_ab / S_bb) b, S_aa - S_ab^2 / S_bb) under the bivariate N(0, Sab)
			slope <- Sab[1, 2] / Sab[2, 2]
			mu_a <- slope * b
			sd_a <- sqrt(Sab[1, 1] - slope * Sab[1, 2])

			# draw each a_i from that normal truncated to [lo_i, hi_i], stable in the tails
			a_new <- .lame_checked_truncated_gaussian(mu_a, sd_a, lo, hi)
			a <- a_new

			# conjugate inverse-Wishart update of Sab, drawn as the inverse of a Wishart draw
			post_scale <- eta0 * Sab0 + crossprod(cbind(a, b))
			Sab <- solve(rwish(solve(post_scale), eta0 + n))
		}

		list(Z = resid + outer(a, rep(1, n)), a = a, Sab = Sab)
	}
