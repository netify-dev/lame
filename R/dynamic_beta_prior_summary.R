# prior-elicitation helper for dynamic_beta. draws sample paths from
# the prior and reports max first-difference, roughness, range, and
# trend correlation so the user can sanity-check before fitting.

#' Summarise the implied prior on a time-varying coefficient path
#'
#' Draws \code{ndraws} sample paths of length \code{T} from the
#' \code{dynamic_beta} prior (given a kind, AR(1) hyperparameters, and an
#' inverse-gamma on the innovation variance) and reports the implied
#' distribution of common summary statistics: maximum absolute first
#' difference, roughness (sum of squared first differences), path range, and
#' correlation with time. Use this before fitting to confirm that your prior
#' is loose / tight enough.
#'
#' For \code{kind = "ar1"}, \code{rho_mean}/\code{rho_sd} are translated to a
#' Beta prior on the standardised \eqn{(\rho - \rho_{lower}) /
#' (\rho_{upper} - \rho_{lower})}. For \code{kind = "rw1"},
#' \eqn{\rho = 1} is fixed and only \eqn{\sigma_\beta^2} is sampled. For
#' \code{kind = "rw2"}, the second-difference variance is sampled
#' and the first two values use a diffuse \eqn{N(0, 10)} initial prior.
#'
#' @param T Length of the path (default 10).
#' @param ndraws Number of paths to draw (default 5000).
#' @param kind One of \code{"ar1"}, \code{"rw1"}, \code{"rw2"},
#'   \code{"matern32"}. Default \code{"ar1"}. For \code{"matern32"} the
#'   path is drawn from a continuous-time Matérn-3/2 GP with length
#'   scale \code{matern32_length_scale} and marginal variance
#'   \eqn{\sigma_\beta^2}.
#' @param rho_mean,rho_sd Prior mean/SD on \eqn{\rho} (AR(1) only). Used to
#'   parameterise a Beta prior on the standardised \eqn{\rho}.
#' @param rho_lower,rho_upper Truncation bounds (AR(1) only). Default
#'   \code{[0, 0.999]}.
#' @param sigma_shape,sigma_scale Inverse-Gamma shape / scale on
#'   \eqn{\sigma_\beta^2}. Defaults \code{shape = 2, scale = 1}.
#' @param matern32_length_scale Length scale for the Matérn-3/2 covariance
#'   (\code{kind = "matern32"} only). Default \code{2}.
#' @param threshold Cutoff for the \code{prob_max_diff_gt_threshold}
#'   summary. Default \code{1}.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{summary}}{Data frame with quantile rows for max-first-diff,
#'         roughness, path range, trend correlation.}
#'   \item{\code{prob_max_diff_gt_threshold}}{Empirical probability that the
#'         maximum absolute first difference exceeds \code{threshold}.}
#'   \item{\code{rho}}{Length \code{ndraws} vector of sampled rho values
#'         (\code{NA} for kinds with fixed rho).}
#'   \item{\code{sigma}}{Length \code{ndraws} vector of sampled sigma values.}
#'   \item{\code{paths}}{\code{ndraws x T} matrix of sample paths.}
#' }
#'
#' @examples
#' \donttest{
#' # Default prior: AR(1) with rho_mean = 0.8, rho_sd = 0.15
#' s <- dynamic_beta_prior_summary(T = 10, ndraws = 2000, seed = 1)
#' s$summary
#' # what's the probability that consecutive beta_t differ by more than 1?
#' s$prob_max_diff_gt_threshold
#'
#' # Tighter prior on innovation variance
#' s_tight <- dynamic_beta_prior_summary(sigma_scale = 0.1, seed = 1)
#' s_tight$summary
#' }
#'
#' @export
dynamic_beta_prior_summary <- function(
	T            = 10L,
	ndraws       = 5000L,
	kind         = c("ar1", "rw1", "rw2", "matern32"),
	rho_mean     = 0.8,
	rho_sd       = 0.15,
	rho_lower    = 0,
	rho_upper    = 0.999,
	sigma_shape  = 2,
	sigma_scale  = 1,
	matern32_length_scale = 2,
	threshold    = 1,
	seed         = NULL
) {
	kind <- match.arg(kind)
	if (T < 2L) {
		stop("`T` must be at least 2.")
	}
	if (!is.null(seed)) set.seed(seed)

	# ar(1) only: translate (rho_mean, rho_sd) to beta(a, b) on the
	# standardised rho. uses the method-of-moments inversion:
	#   m_z = (rho_mean - rho_lower) / (rho_upper - rho_lower)
	#   s_z = rho_sd / (rho_upper - rho_lower)
	#   kappa = m_z(1 - m_z)/s_z^2 - 1
	#   a = m_z * kappa,  b = (1 - m_z) * kappa.
	ab <- NULL
	if (kind == "ar1") {
		m_z <- (rho_mean - rho_lower) / (rho_upper - rho_lower)
		s_z <- rho_sd / (rho_upper - rho_lower)
		if (m_z <= 0 || m_z >= 1) {
			stop("`rho_mean` must be strictly between `rho_lower` and `rho_upper`.")
		}
		kappa <- m_z * (1 - m_z) / (s_z^2) - 1
		if (kappa <= 0) {
			stop("`rho_mean`/`rho_sd` imply invalid Beta prior (variance too large).")
		}
		ab <- c(a = m_z * kappa, b = (1 - m_z) * kappa)
	}

	paths <- matrix(NA_real_, ndraws, T)
	rho   <- rep(NA_real_, ndraws)
	sigma <- rep(NA_real_, ndraws)

	for (m in seq_len(ndraws)) {
		# sample sigma from ig(shape, scale): var = scale / gamma(shape, 1)
		s2 <- sigma_scale / stats::rgamma(1, shape = sigma_shape, rate = 1)
		sigma[m] <- sqrt(s2)

		if (kind == "ar1") {
			z <- stats::rbeta(1, ab["a"], ab["b"])
			r <- rho_lower + (rho_upper - rho_lower) * z
			rho[m] <- r
			# stationary initial draw
			init_sd <- sigma[m] / sqrt(max(1 - r^2, 1e-8))
			paths[m, 1] <- stats::rnorm(1, 0, init_sd)
			for (t in 2:T) {
				paths[m, t] <- r * paths[m, t - 1] + stats::rnorm(1, 0, sigma[m])
			}
		} else if (kind == "rw1") {
			# diffuse initial state; rho = 1 fixed
			rho[m] <- 1
			paths[m, 1] <- stats::rnorm(1, 0, 10)
			for (t in 2:T) {
				paths[m, t] <- paths[m, t - 1] + stats::rnorm(1, 0, sigma[m])
			}
		} else if (kind == "rw2") {
			# rw2: second-difference innovation
			rho[m] <- NA_real_
			paths[m, 1] <- stats::rnorm(1, 0, 10)
			paths[m, 2] <- paths[m, 1] + stats::rnorm(1, 0, sigma[m])
			if (T >= 3) {
				for (t in 3:T) {
					paths[m, t] <- 2 * paths[m, t - 1] - paths[m, t - 2] +
						stats::rnorm(1, 0, sigma[m])
				}
			}
		} else {
			# matern32: draw from continuous-time matérn-3/2 gp with marginal
			# variance sigma^2 and length scale matern32_length_scale.
			# k(d) = sigma^2 * (1 + sqrt(3) d / l) * exp(-sqrt(3) d / l)
			rho[m] <- NA_real_
			ts <- seq_len(T)
			D <- abs(outer(ts, ts, "-"))
			l_scale <- max(matern32_length_scale, .Machine$double.eps)
			r3 <- sqrt(3) * D / l_scale
			K <- sigma[m]^2 * (1 + r3) * exp(-r3)
			# tiny jitter for numerical stability
			K <- K + diag(1e-8, T)
			L <- tryCatch(chol(K), error = function(e) {
				# fall back to a diagonal draw on degenerate cov
				diag(sqrt(diag(K)), T, T)
			})
			paths[m, ] <- as.numeric(crossprod(L, stats::rnorm(T)))
		}
	}

	# summary statistics on the sampled paths
	diffs       <- t(apply(paths, 1, diff))
	max_abs_d1  <- apply(abs(diffs), 1, max)
	roughness   <- rowSums(diffs^2)
	path_range  <- apply(paths, 1, function(x) diff(range(x)))
	trend_cor   <- apply(paths, 1, function(x) {
		# guard against constant paths (correlation undefined)
		if (stats::sd(x) < 1e-10) NA_real_ else stats::cor(x, seq_len(T))
	})

	qs <- function(v) {
		v <- v[is.finite(v)]
		if (length(v) == 0) return(c(q05 = NA_real_, q50 = NA_real_, q95 = NA_real_))
		stats::quantile(v, probs = c(0.05, 0.50, 0.95), names = FALSE)
	}

	summary_df <- data.frame(
		quantity = c("max_abs_first_diff", "roughness", "range", "trend_cor"),
		q05      = c(qs(max_abs_d1)[1], qs(roughness)[1], qs(path_range)[1], qs(trend_cor)[1]),
		q50      = c(qs(max_abs_d1)[2], qs(roughness)[2], qs(path_range)[2], qs(trend_cor)[2]),
		q95      = c(qs(max_abs_d1)[3], qs(roughness)[3], qs(path_range)[3], qs(trend_cor)[3]),
		row.names = NULL
	)

	list(
		summary                    = summary_df,
		prob_max_diff_gt_threshold = mean(max_abs_d1 > threshold, na.rm = TRUE),
		rho                        = rho,
		sigma                      = sigma,
		paths                      = paths,
		kind                       = kind
	)
}
