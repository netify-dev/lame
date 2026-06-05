# heuristic change-point diagnostic for dynamic_beta fits.
# compares the max absolute scaled first difference of the posterior
# beta path to a prior-null reference simulated from the fit's
# rho_beta / sigma_beta. large values flag abrupt breaks. this is a
# heuristic tail-ratio score, not a formal bayes factor.

#' Detect potential change points in a dynamic_beta posterior path
#'
#' For each dynamic coefficient in a \code{dynamic_beta} fit, compares the
#' posterior distribution of the maximum scaled first-difference
#' \eqn{M = \max_t |\beta_t - \beta_{t-1}| / \sigma_\beta} to a Monte-Carlo
#' approximation of the same statistic under the AR(1) prior. The returned
#' \code{bf} column is the tail-ratio score
#' \eqn{\Pr(M^{post} > m^*) / 0.05}, where \eqn{m^*} is the 95\%
#' quantile of \eqn{M^{prior}}. It is kept under the historical column name
#' for compatibility, but it is not a marginal-likelihood Bayes factor. Use it
#' to surface posterior temporal jumps that the AR(1) prior cannot comfortably
#' accommodate.
#'
#' @param fit A fitted \code{lame} object with \code{dynamic_beta}.
#' @param coefs Optional character vector of coefficient names to check.
#'   Defaults to every dynamic coefficient on the fit.
#' @param n_prior_sims Number of prior simulations for the null distribution.
#'   Default \code{2000}.
#' @param threshold_bf Cutoff for the "warn" column. Default \code{10};
#'   lower to 5 if the false-positive rate on AR(1) data is too high.
#' @param seed Optional seed for reproducibility.
#'
#' @return A data frame with one row per dynamic coefficient: \code{coef},
#'   \code{bf} (the heuristic tail-ratio score), \code{m_post_mean} (posterior
#'   mean of \eqn{M}), \code{m_prior_q95} (95\% quantile of the prior
#'   \eqn{M}), \code{t_hat} (period of the largest scaled jump), \code{warn}
#'   (logical, \code{bf > threshold_bf}).
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'             dynamic_beta = "dyad",
#'             nscan = 200, burn = 50, odens = 5, verbose = FALSE)
#' detect_change_point(fit)
#' }
#'
#' @export
detect_change_point <- function(
	fit,
	coefs        = NULL,
	n_prior_sims = 2000L,
	threshold_bf = 10,
	seed         = NULL
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

	# pick the dynamic coefficients (mask preferred over all-coefs)
	dn        <- dimnames(fit$BETA)
	coef_nms  <- if (!is.null(dn)) dn[[2L]] else paste0("v", seq_len(dim(fit$BETA)[2L]))
	mask      <- fit$beta_dynamic_mask
	default_c <- if (!is.null(mask)) coef_nms[mask] else coef_nms
	if (is.null(coefs)) coefs <- default_c

	# which block each requested coef belongs to (for rho/sigma sampling)
	groups <- fit$beta_dynamic_groups
	# fallback if missing: assume all in the first block
	if (is.null(groups)) groups <- rep(colnames(fit$RHO_BETA)[1], length(coef_nms))

	rho_draws_by_block   <- as.data.frame(fit$RHO_BETA, stringsAsFactors = FALSE)
	sigma_draws_by_block <- as.data.frame(fit$SIGMA_BETA, stringsAsFactors = FALSE)

	out_rows <- vector("list", length(coefs))
	for (i in seq_along(coefs)) {
		cf <- coefs[i]
		j  <- match(cf, coef_nms)
		if (is.na(j)) {
			cli::cli_warn("Coefficient {.val {cf}} not found in fit; skipping.")
			next
		}
		blk <- groups[j]
		if (!nzchar(blk) || !(blk %in% colnames(fit$RHO_BETA))) {
			# coef is static or block label missing; skip
			next
		}
		# posterior m
		beta_path <- fit$BETA[, j, , drop = FALSE]  # iter x 1 x t
		# squeeze the singleton coef dim out
		beta_path <- matrix(beta_path, nrow = dim(beta_path)[1L], ncol = dim(beta_path)[3L])
		# scale each iteration's diffs by that iteration's sigma_beta
		# (per-block sigma; align by row index)
		sb_block <- fit$SIGMA_BETA[, blk]
		# safety: align lengths
		n_iter <- min(nrow(beta_path), length(sb_block))
		beta_path <- beta_path[seq_len(n_iter), , drop = FALSE]
		sb_block  <- sb_block[seq_len(n_iter)]

		diffs <- abs(t(apply(beta_path, 1, diff)))   # iter x (t-1)
		scaled_diffs <- diffs / pmax(sb_block, 1e-8)
		M_post <- apply(scaled_diffs, 1, max)
		# which period was the largest jump? take posterior mode (= 1 + most common argmax)
		t_jumps <- apply(scaled_diffs, 1, which.max) + 1L
		t_hat <- as.integer(names(sort(table(t_jumps), decreasing = TRUE))[1L])

		# prior simulation: draw (rho, sigma) from the per-block posterior and
		# generate prior ar(1) paths of length t
		Tn <- ncol(beta_path) + 1L
		# resample (rho, sigma) from the posterior with replacement
		idx <- sample.int(n_iter, n_prior_sims, replace = TRUE)
		rho_s   <- fit$RHO_BETA[idx, blk]
		sigma_s <- fit$SIGMA_BETA[idx, blk]
		M_prior <- vapply(seq_len(n_prior_sims), function(k) {
			r <- rho_s[k]; s <- sigma_s[k]
			# stationary initial draw, then t-1 innovations
			path <- numeric(Tn)
			path[1] <- stats::rnorm(1, 0, s / sqrt(max(1 - r * r, 1e-8)))
			for (t in 2:Tn) {
				path[t] <- r * path[t - 1] + stats::rnorm(1, 0, s)
			}
			max(abs(diff(path)) / s)
		}, numeric(1))

		m_star <- stats::quantile(M_prior, 0.95)
		bf <- mean(M_post > m_star) / 0.05

		out_rows[[i]] <- data.frame(
			coef         = cf,
			bf           = bf,
			m_post_mean  = mean(M_post),
			m_prior_q95  = as.numeric(m_star),
			t_hat        = t_hat,
			warn         = bf > threshold_bf,
			stringsAsFactors = FALSE
		)
	}
	out_rows <- out_rows[!vapply(out_rows, is.null, logical(1))]
	do.call(rbind, out_rows)
}
