.lame_binary_row_truncnorm <- function(mean, sd, lower, upper) {
	n <- length(mean)
	if (!is.numeric(mean) || any(!is.finite(mean)) || length(sd) != 1L ||
			!is.finite(sd) || sd <= 0 || length(lower) != n || length(upper) != n ||
			anyNA(lower) || anyNA(upper) || any(lower >= upper)) {
		cli::cli_abort(
			"Binary row effects require a finite Gaussian mean, positive scale and nonempty truncation intervals.",
			class = "lame_binary_row_conditional_error"
		)
	}
	lo <- (lower - mean) / sd
	hi <- (upper - mean) / sd
	if (anyNA(lo) || anyNA(hi) || any(lo >= hi) ||
			any(is.finite(lower) & !is.finite(lo)) ||
			any(is.finite(upper) & !is.finite(hi))) {
		cli::cli_abort(
			"Binary row truncation cannot be represented at the supplied Gaussian scale.",
			class = "lame_binary_row_numerical_error"
		)
	}
	p_low <- pnorm(lo)
	p_high <- pnorm(hi)
	central <- lo < 5 & hi > -5 &
		p_high - p_low > 64 * .Machine$double.eps * pmax(p_low, p_high)
	if (all(central)) {
		# central intervals use the plain inverse-cdf draw
		out <- mean + sd * qnorm(runif(n, p_low, p_high))
	} else {
		# tails use survival and lower-tail log masses to avoid rounded cdf endpoints
		positive <- !central & lo >= 0
		negative <- !central & hi <= 0
		crossing <- !central & !positive & !negative
		log_low <- pnorm(lo, log.p = TRUE)
		log_high <- pnorm(hi, log.p = TRUE)
		log_surv_low <- pnorm(lo, lower.tail = FALSE, log.p = TRUE)
		log_surv_high <- pnorm(hi, lower.tail = FALSE, log.p = TRUE)
		if (any(positive & log_surv_low <= log_surv_high) ||
				any((negative | crossing) & log_low >= log_high)) {
			cli::cli_abort(
				"Binary row truncation has an unresolved log-probability interval.",
				class = "lame_binary_row_numerical_error"
			)
		}
		u <- runif(n)
		z <- numeric(n)
		z[central] <- qnorm(
			p_low[central] +
				(p_high[central] - p_low[central]) * u[central]
		)
		if (any(positive)) {
			mass <- -expm1(log_surv_high[positive] - log_surv_low[positive])
			lp <- log_surv_low[positive] + log1p(-u[positive] * mass)
			z[positive] <- qnorm(lp, lower.tail = FALSE, log.p = TRUE)
		}
		if (any(negative | crossing)) {
			idx <- negative | crossing
			mass <- -expm1(log_low[idx] - log_high[idx])
			lp <- log_high[idx] + log1p(-(1 - u[idx]) * mass)
			z[idx] <- qnorm(lp, log.p = TRUE)
		}
		out <- mean + sd * z
	}
	if (any(!is.finite(out)) || any(out <= lower | out >= upper)) {
		cli::cli_abort(
			"Binary row draw cannot be represented strictly inside its truncation interval.",
			class = "lame_binary_row_numerical_error"
		)
	}
	out
}

#' Sample additive row effects and their covariance for the binary family
#'
#' Joint Gibbs update of the social-relations additive row effects \code{a}
#' and the additive-effect covariance \code{Sab} under a probit (binary)
#' likelihood.  Given the latent normal scores \code{Z} and the observed 0/1
#' matrix \code{Y}, each row effect is redrawn from its Gaussian full
#' conditional (conditioned on the column effects \code{b} through \code{Sab})
#' truncated to the interval that keeps every latent score sign-compatible with
#' the observed responses.  The covariance \code{Sab} is then refreshed from its
#' inverse-Wishart full conditional.
#'
#' @param Z square matrix of current latent normal scores
#' @param Y square binary relational matrix (NA entries allowed)
#' @param a current value of the row effects
#' @param b current value of the column effects
#' @param Sab current 2x2 covariance of the additive effects
#' @param Sab0 prior inverse-scale matrix (defaults to the 2x2 identity)
#' @param eta0 prior degrees of freedom (defaults to 4)
#' @param SS number of inner Gibbs sweeps
#' @return list with the updated \code{Z}, \code{a} and \code{Sab}
#' @keywords internal
#' @author lame authors
#' @export raSab_bin_fc
raSab_bin_fc <-
	function(Z, Y, a, b, Sab, Sab0 = NULL, eta0 = NULL, SS = round(sqrt(nrow(Z)))) {
		if (is.null(Sab0)) { Sab0 <- diag(2) }
		if (is.null(eta0)) { eta0 <- 4 }

		n <- nrow(Z)

		# one inverse-Wishart draw with scale psi and df degrees of freedom, via a
		# Bartlett-decomposed Wishart on the inverse scale
		draw_inv_wishart <- function(psi, df) {
			p <- nrow(psi)
			root <- t(chol(solve(psi))) # lower factor: solve(psi) = root %*% t(root)
			bart <- matrix(0, p, p)
			diag(bart) <- sqrt(rchisq(p, df - seq_len(p) + 1))
			if (p > 1L) { bart[lower.tri(bart)] <- rnorm(p * (p - 1L) / 2L) }
			wfac <- root %*% bart # wfac %*% t(wfac) ~ Wishart(solve(psi), df)
			solve(tcrossprod(wfac)) # its inverse is the IW(psi, df) draw
		}

		# strip the row effect from each latent score; the remainder is held fixed
		# while a is resampled, so Z_ij = carry_ij + a_i throughout
		carry <- Z - matrix(a, n, n)

		# probit sign constraints on Z_ij = carry_ij + a_i give, per row i,
		#   Y_ij == 1  ->  a_i > -carry_ij   (a lower limit)
		#   Y_ij == 0  ->  a_i < -carry_ij   (an upper limit)
		# missing cells contribute no limit
		cut <- -carry
		low_src <- cut ; low_src[is.na(Y) | Y != 1] <- -Inf
		high_src <- cut ; high_src[is.na(Y) | Y != 0] <- Inf
		a_low <- apply(low_src, 1, max)
		a_high <- apply(high_src, 1, min)
		a_low[is.na(a_low)] <- -Inf
		a_high[is.na(a_high)] <- Inf

		for (ss in seq_len(SS)) {
			# a | b is Gaussian through the 2x2 covariance Sab
			slope <- Sab[1, 2] / Sab[2, 2]
			cond_mean <- slope * b
			cond_sd <- sqrt(Sab[1, 1] - Sab[1, 2]^2 / Sab[2, 2])

			# tail-stable truncated normal draw of a
			a <- .lame_binary_row_truncnorm(cond_mean, cond_sd, a_low, a_high)

			# Sab ~ inverse-Wishart(eta0 * Sab0 + [a b]'[a b], eta0 + n)
			ab <- cbind(a, b)
			Sab <- draw_inv_wishart(eta0 * Sab0 + crossprod(ab), eta0 + n)
		}

		list(Z = carry + matrix(a, n, n), a = a, Sab = Sab)
	}
