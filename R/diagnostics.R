# reference for coefficient-path roughness: the largest jump expected by chance (not a Bayes factor)

#' Detect potential change points in a dynamic_beta posterior path
#'
#' For each selected dynamic coefficient, compares the largest absolute first
#' difference divided by its fitted innovation standard deviation
#' \eqn{\sigma_\beta \sqrt{\Lambda_{jj}}} with a simulated reference.
#' Each reference path has the same number of periods as the fitted path,
#' starts from the first state of the matched posterior draw, and follows that
#' draw's persistence and innovation scale: it is the largest jump a path that
#' begins where the fitted path begins makes by chance. The design-based
#' innovation scale and the averaging of symmetric node coefficients are retained.
#' This is a posterior-conditioned transition reference, not an unconditional
#' prior predictive check or a marginal-likelihood Bayes factor.
#'
#' The \code{bf} column is the tail ratio
#' \eqn{\Pr(M^{post} > m^*) / 0.05}, where \eqn{m^*} is the reference's
#' 95th percentile. Its range is zero to twenty, in steps of 20 over the number
#' of stored draws. Values above \code{threshold_bf} flag unusually abrupt
#' posterior paths for inspection; they do not establish a change point.
#' The reference supports \code{"ar1"} and \code{"rw1"} on equally spaced
#' periods. RW2, Matern, irregular periods and actor-specific coefficient laws
#' are refused rather than assigned an AR(1) reference. Saved named training
#' designs, the dynamic mask, and aligned transition draws are required.
#'
#' @section Measured properties:
#' The table gives the percentage of simulated panels with
#' \code{bf > threshold_bf} at the default threshold, from 200 fits of each
#' path type and panel shape. Every panel has one standard-normal dyadic
#' covariate whose coefficient follows the named path, is fitted with
#' \code{dynamic_beta = "dyad"} (\code{dynamic_beta_kind = "rw1"} for
#' random-walk paths, \code{"ar1"} otherwise) and keeps 40 stored draws (12 for
#' the binary panels, which have 50 actors; the others have 15). The outcome is
#' normal, binary or Poisson. \emph{Smooth} paths are a constant, a linear
#' ramp and two stationary AR(1) paths (persistence 0.7, innovation standard
#' deviation 0.15 and 0.4); \emph{in-model} paths are random walks and AR(1)
#' paths drawn from the prior of the fit; a \emph{step} is a jump of 0.6 or 1.2
#' (from 0.2 to 0.8 or to 1.4) at the middle of the panel. Normal outcomes have
#' residual standard deviation 0.4; binary outcomes are probit with intercept
#' -0.3; Poisson outcomes have log-rate offset -0.5 and no AR(1) paths drawn
#' from the prior. \emph{Worst path type} is the highest rate among the smooth
#' and in-model path types of the row.
#' \tabular{rlrrrrr}{
#'   \strong{T} \tab \strong{outcome} \tab \strong{smooth} \tab \strong{in-model} \tab \strong{worst path type} \tab \strong{step 0.6} \tab \strong{step 1.2} \cr
#'   4 \tab binary \tab 1.2 \tab 2.5 \tab 5.0 \tab 2.5 \tab 20.5 \cr
#'   5 \tab normal \tab 1.2 \tab 3.2 \tab 4.0 \tab 4.0 \tab 20.0 \cr
#'   5 \tab Poisson \tab 0.8 \tab 0.0 \tab 3.0 \tab 2.5 \tab 20.5 \cr
#'   8 \tab normal \tab 1.1 \tab 1.8 \tab 3.0 \tab 56.0 \tab 85.0 \cr
#'   12 \tab normal \tab 1.9 \tab 3.0 \tab 5.5 \tab 100.0 \tab 100.0 \cr
#' }
#' Pooled over panel shapes, smooth paths were flagged in 1.2 percent of the
#' fits and in-model paths in 2.3 percent. At a threshold of 9 the pooled
#' in-model rate is 5.0 percent, with no margin, and at 10 it is 2.3 percent,
#' so 10 is the default. The tail ratio
#' judges a jump against the innovation scale the fit estimates from the same
#' path, so a panel of four or five periods rarely shows a single jump as
#' unusual (about a fifth of the panels with a step of 1.2 are flagged), while
#' panels of eight or more periods flag it in most cases.
#'
#' @param fit A fitted \code{lame} object with \code{dynamic_beta}.
#' @param coefs Optional character vector of coefficient names to check.
#'   Defaults to every dynamic coefficient on the fit.
#' @param n_prior_sims Number of simulations for the transition reference.
#'   Default \code{2000}.
#' @param threshold_bf Cutoff for the "warn" column. Default \code{10}, the
#'   the threshold at which the measured false-positive rate is below 5
#'   percent for smooth and for in-model paths (see Measured properties).
#' @param seed Optional seed for reproducibility.
#'
#' @return A data frame with one row per dynamic coefficient: \code{coef},
#'   \code{bf} (the heuristic tail-ratio score), \code{m_post_mean} (posterior
#'   mean of \eqn{M}), \code{m_prior_q95} (95\% quantile of the reference
#'   \eqn{M}), \code{t_hat} (period of the largest scaled jump), \code{warn}
#'   (logical, \code{bf > threshold_bf}).
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 60, burn = 15, odens = 5, verbose = FALSE, chains = 1)
#' detect_change_point(fit)
#' }
#'
#' @export
detect_change_point <- function(
	fit,
	coefs = NULL,
	n_prior_sims = 2000L,
	threshold_bf = 10,
	seed = NULL
) {
	if (!isTRUE(fit$dynamic_beta)) {
		cli::cli_abort(c(
			"{.fn detect_change_point} requires a fit with {.code dynamic_beta} active.",
			"i" = "Refit with {.code dynamic_beta = TRUE} or a per-coefficient subset."
		))
	}
	if (is.null(fit$BETA) || length(dim(fit$BETA)) != 3L) {
		cli::cli_abort("{.code fit$BETA} must be a 3-D array (iter x p x T).")
	}
	if (is.null(fit$RHO_BETA) || is.null(fit$SIGMA_BETA)) {
		cli::cli_abort("{.code fit$RHO_BETA} and {.code fit$SIGMA_BETA} required for the prior simulation.")
	}
	if (!is.null(seed)) set.seed(seed)

	if (!is.numeric(n_prior_sims) || length(n_prior_sims) != 1L ||
			!is.finite(n_prior_sims) || n_prior_sims < 1 || n_prior_sims != floor(n_prior_sims) ||
			!is.numeric(threshold_bf) || length(threshold_bf) != 1L ||
			!is.finite(threshold_bf) || threshold_bf < 0) {
		cli::cli_abort("Supply a positive integer n_prior_sims and a finite nonnegative threshold_bf.")
	}
	geometry <- .lame_beta_transition_scale(fit)
	B <- fit$BETA; S <- dim(B)[1L]; Tn <- dim(B)[3L]
	if (S < 1L || Tn < 2L || any(!is.finite(B))) {
		cli::cli_abort("Change-point diagnostics require finite paths with at least two periods.")
	}
	coef_nms <- dimnames(B)[[2L]]
	mask <- fit$beta_dynamic_mask
	if (is.null(coefs)) coefs <- coef_nms[mask]
	if (!is.character(coefs) || !length(coefs) || anyNA(coefs) || anyDuplicated(coefs) ||
			!all(coefs %in% coef_nms[mask])) {
		cli::cli_abort("coefs must select named dynamic coefficients from the fit.")
	}
	groups <- fit$beta_dynamic_groups
	R <- fit$RHO_BETA; V <- fit$SIGMA_BETA
	if (length(groups) != length(coef_nms) || anyNA(groups[mask]) ||
			!is.matrix(R) || !is.matrix(V) || nrow(R) != S || nrow(V) != S ||
			!all(groups[mask] %in% colnames(R)) || !all(groups[mask] %in% colnames(V)) ||
			any(!is.finite(R)) || any(!is.finite(V)) || any(V <= 0)) {
		cli::cli_abort(
			"Coefficient paths and transition draws must be finite, positive-scaled and aligned.",
			class = "lame_beta_transition"
		)
	}
	kind <- fit$dynamic_beta_kind %||% "ar1"
	if ((kind == "ar1" && any(abs(R) >= 1)) || (kind == "rw1" && any(R != 1))) {
		cli::cli_abort(
			"Persistence draws do not match the fitted coefficient transition model.",
			class = "lame_beta_transition"
		)
	}
	out_rows <- lapply(coefs, function(cf) {
		j <- match(cf, coef_nms); blk <- groups[j]
		beta_path <- matrix(B[, j, ], S, Tn)
		sigma <- V[, blk] * geometry$sd[j]
		diffs <- abs(beta_path[, -1L, drop = FALSE] - beta_path[, -Tn, drop = FALSE])
		scaled_diffs <- diffs / sigma
		M_post <- apply(scaled_diffs, 1L, max)
		t_jumps <- apply(scaled_diffs, 1L, which.max) + 1L
		t_hat <- as.integer(names(sort(table(t_jumps), decreasing = TRUE))[1L])
		# reference paths start at the draw's first state and follow its AR(1) transition
		idx <- sample.int(S, n_prior_sims, replace = TRUE)
		M_prior <- vapply(idx, function(k) {
			path <- numeric(Tn)
			path[1L] <- beta_path[k, 1L]
			for (tt in 2:Tn) path[tt] <- R[k, blk] * path[tt - 1L] + stats::rnorm(1L, 0, sigma[k])
			max(abs(diff(path)) / sigma[k])
		}, numeric(1))
		m_star <- unname(stats::quantile(M_prior, 0.95))
		bf <- mean(M_post > m_star) / 0.05
		data.frame(
			coef = cf,
			bf = bf,
			m_post_mean = mean(M_post),
			m_prior_q95 = m_star,
			t_hat = t_hat,
			warn = bf > threshold_bf,
			stringsAsFactors = FALSE
		)
	})
	do.call(rbind, out_rows)
}
