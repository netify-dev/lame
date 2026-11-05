#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood
#'
#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood
#'
#'
#' @usage raSab_cbin_fc(Z, Y, a, b, Sab, odmax, odobs, Sab0=NULL, eta0=NULL,SS =
#' round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
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
#' @export raSab_cbin_fc
raSab_cbin_fc <-
	function(Z, Y, a, b, Sab, odmax, odobs, Sab0 = NULL, eta0 = NULL, SS = round(sqrt(nrow(Z)))) {

		# default inverse-Wishart prior on Sab
		if (is.null(Sab0)) { Sab0 <- diag(2) }
		if (is.null(eta0)) { eta0 <- 4 }

		n <- nrow(Z)

		# rows with hidden cells move given the nominations their current latent values imply
		# (see rZ_nomination_hidden.R); the hidden-cell draw in rZ_cbin_fc lets those change
		hid <- .lame_nom_hidden_cells(Y, FALSE)
		if (any(hid)) {
			Y <- .lame_cbin_implied_labels(Z, Y, odmax, hid)
			odobs <- rowSums(!is.na(Y) & Y > 0)
		}

		# strip the current row effect: a proposed a_i gives latent value resid[i,j] + a_i
		resid <- Z - matrix(a, nrow = n, ncol = n)

		# each observed cell bounds a_i on one side: a_i > -resid[i,j] (Y = 1), a_i < -resid[i,j] (Y = 0)
		thr <- -resid

		# tightest interval [lo_i, hi_i] for each row effect; missing cells give no bound
		one_cell <- !is.na(Y) & (Y == 1)
		zero_cell <- !is.na(Y) & (Y == 0)

		lo_mat <- matrix(-Inf, nrow = n, ncol = n)
		hi_mat <- matrix(Inf, nrow = n, ncol = n)
		lo_mat[one_cell] <- thr[one_cell]
		hi_mat[zero_cell] <- thr[zero_cell]

		lo <- apply(lo_mat, 1, max)
		hi <- apply(hi_mat, 1, min)

		# censored binary: a row at its nomination cap (odobs >= odmax) may have suppressed
		# ties, so its zeros carry no information and the ceiling is released
		hi[odobs >= odmax] <- Inf

		# Gibbs sweeps alternating a | (b, Sab) and Sab | (a, b)
		for (iter in seq_len(SS)) {

			# conditional of a given b under (a, b) ~ N(0, Sab):
			#   E[a|b] = b * Sab_ab / Sab_bb,  Var[a|b] = Sab_aa - Sab_ab^2 / Sab_bb
			cond_mean <- b * (Sab[1, 2] / Sab[2, 2])
			cond_sd <- sqrt(Sab[1, 1] - Sab[1, 2]^2 / Sab[2, 2])

			# draw each a_i from that normal truncated to [lo_i, hi_i], stable in the tails
			a_new <- .lame_checked_truncated_gaussian(cond_mean, cond_sd, lo, hi)
			a <- a_new

			# inverse-Wishart update of Sab: scale eta0*Sab0 + [a b]'[a b], df eta0 + n
			ab_mat <- cbind(a, b)
			S_post <- eta0 * Sab0 + t(ab_mat) %*% ab_mat
			Sab <- solve(rwish(solve(S_post), eta0 + n))
		}

		# rebuild Z with the refreshed row effect
		list(Z = resid + matrix(a, nrow = n, ncol = n), a = a, Sab = Sab)
	}
