#' Convert an ALS fit to MCMC starting values
#'
#' Builds a \code{start_vals} list for \code{\link{ame}} or \code{\link{lame}}
#' from an ALS point estimate. This is useful when the ALS fit has already found
#' a good latent-space solution and the MCMC chain should start near that
#' solution rather than from diffuse random values.
#'
#' @param fit a fitted \code{\link{ame_als}}, \code{\link{lame_als}}, or
#'   dynamic ALS object returned by \code{\link{lame}} with
#'   \code{method = "als"}.
#' @param jitter non-negative standard deviation for independent Gaussian
#'   perturbations added to numeric starting values. Use a small positive value
#'   for multiple MCMC chains.
#' @param seed optional integer seed used for the jitter.
#'
#' @return A list that can be passed to \code{start_vals}. For a dynamic ALS
#'   fit it also carries \code{rho_uv}, \code{sigma_uv}, \code{rho_ab}, and
#'   \code{sigma_ab}, so the MCMC run starts its AR(1) hyperparameters from
#'   the ALS solution. The dyadic \code{rho} start is the single pooled
#'   correlation in \code{fit$VC["rho"]}, matching ALS simulation. Period-specific
#'   \code{rho_path} / \code{residual_rho_path} values are residual diagnostics,
#'   including in legacy objects labelled \code{dynamic_rho = TRUE}, and are
#'   not used as covariance-path starting values. A MCMC model requesting
#'   \code{dynamic_rho = TRUE} expands this scalar initial value; explicit
#'   MCMC path starts remain available through its \code{start_vals} argument.
#'   For a penalised dynamic fit (\code{bias_correct = "none"} or
#'   \code{"jackknife"}) read these as starting values
#'   rather than estimates:
#'   the persistences are the values the ALS smoothness penalty used (the
#'   prior means unless \code{prior} set them), and the innovation scales are
#'   implied by the fitted paths, which the penalty shrinks, so they run
#'   below the posterior scale. For a dynamic empirical-Bayes fit (the default
#'   for normal, binary and Poisson panels) they are estimates under the
#'   \code{lame()} hyperpriors. For a binary fit the intercept is shifted so
#'   that the implied tie rate, \code{mean(pnorm(EZ))}, equals the observed
#'   one, with everything else held at the ALS values; the shift is stored in
#'   \code{attr(, "level_shift")}. The ALS intercept is a working-response
#'   level: under the default IRLS it sits within a few hundredths of the
#'   level the sampler has to reach, and under \code{non_normal_method =
#'   "transform"} it can sit a probit unit or more above it, the sampler's
#'   slowest direction.
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' fit_als <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
#' sv <- als_start_vals(fit_als)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit_mcmc <- ame(Y, R = 1, start_vals = sv,
#'                 nscan = 20, burn = 5, odens = 1, verbose = FALSE, chains = 1)
#' @export
als_start_vals <- function(fit, jitter = 0, seed = NULL) {
	# restore a withheld product so starts use its variational means
	fit <- .ae_product_restore(fit)
	`%||%` <- function(x, y) if (is.null(x)) y else x

	if (!inherits(fit, c("ame_als", "lame_als"))) {
		cli::cli_abort("{.arg fit} must be an ALS fit from {.fn ame_als}, {.fn lame_als}, or {.code lame(method = \"als\")}.")
	}
	if (!is.numeric(jitter) || length(jitter) != 1L ||
			!is.finite(jitter) || jitter < 0) {
		cli::cli_abort("{.arg jitter} must be a non-negative finite number.")
	}
	if (!is.null(seed)) {
		old_seed_exists <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
		if (old_seed_exists) old_seed <- get(".Random.seed", envir = globalenv())
		on.exit({
			if (old_seed_exists) {
				assign(".Random.seed", old_seed, envir = globalenv())
			} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
				rm(".Random.seed", envir = globalenv())
			}
		}, add = TRUE)
		set.seed(seed)
	}

	add_jitter <- function(x, scale = jitter) {
		if (is.null(x) || scale <= 0) return(x)
		if (!is.numeric(x)) return(x)
		out <- x
		ok <- is.finite(out)
		out[ok] <- out[ok] + stats::rnorm(sum(ok), 0, scale)
		dim(out) <- dim(x)
		dimnames(out) <- dimnames(x)
		names(out) <- names(x)
		out
	}

	clean_vec <- function(x) {
		if (is.null(x)) return(NULL)
		nm <- names(x)
		x <- as.numeric(x)
		names(x) <- nm
		x[!is.finite(x)] <- 0
		x
	}

	beta <- fit$coefficients %||% fit$beta
	if (is.null(beta) && !is.null(fit$coef_path)) {
		beta <- rowMeans(fit$coef_path, na.rm = TRUE)
	}
	if (is.null(beta) && !is.null(fit$BETA)) {
		if (length(dim(fit$BETA)) == 3L) {
			beta <- apply(fit$BETA, 2L, mean, na.rm = TRUE)
		} else {
			beta <- colMeans(fit$BETA, na.rm = TRUE)
		}
	}
	beta <- clean_vec(beta)

	get_s2 <- function(fit) {
		# probit latent variance is 1; the ALS s2 is on the IRLS working scale
		if (identical(fit$family, "binary")) return(1)
		candidates <- c(fit$s2, fit$VC[["ve"]])
		candidates <- candidates[is.finite(candidates) & candidates > 0]
		if (length(candidates) == 0L) return(NULL)
		unname(candidates[[1L]])
	}

	get_sab <- function(fit) {
		vc <- fit$VC
		if (is.null(vc)) return(NULL)
		if (all(c("va", "cab", "vb") %in% names(vc))) {
			S <- matrix(
				c(
					vc[["va"]],
					vc[["cab"]],
					vc[["cab"]],
					vc[["vb"]]
				),
				2L,
				2L
			)
		} else if ("va" %in% names(vc)) {
			S <- matrix(vc[["va"]], 1L, 1L)
		} else {
			return(NULL)
		}
		# variance floor: 1e-6 of the response variance for normal, 1e-6 otherwise
		unit <- if (identical(fit$family, "normal") && !is.null(fit$Y)) .lame_unit_var(fit$Y, "normal") else 1
		.lame_start_sab_covariance(S, floor = 1e-6 * unit)
	}

	a_start <- if (isTRUE(fit$dynamic_ab)) {
		fit$a_dynamic %||% fit$APM %||% fit$a
	} else {
		fit$APM %||% fit$a
	}
	b_start <- if (isTRUE(fit$dynamic_ab)) {
		fit$b_dynamic %||% fit$BPM %||% fit$b
	} else {
		fit$BPM %||% fit$b
	}
	# effects a binary ALS fit leaves unestimated start at zero
	a_start <- .ae_na_zero(a_start)
	b_start <- .ae_na_zero(b_start)
	U_start <- fit$U
	# symmetric fits carry U and L; the sampler uses V = U L
	V_start <- fit[["V"]]
	if (is.null(V_start) && !is.null(U_start) && !is.null(fit$L)) {
		L_s <- if (is.matrix(fit$L)) fit$L else diag(fit$L, nrow = length(fit$L))
		V_start <- U_start %*% L_s
	}
	if (is.null(V_start)) V_start <- U_start

	out <- list(
		beta = add_jitter(beta, jitter / 5),
		a = add_jitter(a_start),
		b = add_jitter(b_start),
		U = add_jitter(U_start),
		V = add_jitter(V_start),
		# ALS fits carry a single pooled dyadic correlation
		rho = clean_vec(fit$VC[["rho"]] %||% 0),
		s2 = get_s2(fit),
		Sab = get_sab(fit)
	)

	# keep starting correlations inside (-1, 1)
	boundary_rho <- abs(out$rho) == 1
	out$rho[boundary_rho] <- sign(out$rho[boundary_rho]) * 0.99

	if (!is.null(fit$rho_uv)) out$rho_uv <- fit$rho_uv
	if (!is.null(fit$sigma_uv)) out$sigma_uv <- fit$sigma_uv
	if (!is.null(fit$rho_ab)) out$rho_ab <- fit$rho_ab
	if (!is.null(fit$sigma_ab)) out$sigma_ab <- fit$sigma_ab
	if (!is.null(fit[["G"]])) out[["G"]] <- add_jitter(fit[["G"]], jitter / 5)

	# shift a binary intercept so the implied tie rate matches the data; the
	# intercept is the sampler's slowest direction on sparse panels
	if (identical(fit$family, "binary") && !is.null(out$beta) && length(out$beta) >= 1L &&
			!is.null(fit$EZ) && !is.null(fit$Y)) {
		ez <- unlist(fit$EZ, use.names = FALSE)
		yy <- unlist(fit$Y, use.names = FALSE)
		ok <- is.finite(ez) & is.finite(yy)
		p_obs <- mean(yy[ok] > 0)
		if (sum(ok) > 0 && is.finite(p_obs) && p_obs > 0 && p_obs < 1) {
			shift <- tryCatch(
				stats::uniroot(
					function(d) mean(stats::pnorm(ez[ok] + d)) - p_obs,
					c(-20, 20)
				)$root,
				error = function(e) NA_real_
			)
			if (is.finite(shift)) {
				k <- which(names(out$beta) == "intercept")
				if (length(k) == 0L) k <- 1L
				out$beta[k[1L]] <- out$beta[k[1L]] + shift
				level_shift <- shift
			}
		}
	}

	out <- out[!vapply(out, is.null, logical(1))]
	attr(out, "source") <- class(fit)[[1L]]
	attr(out, "jitter") <- jitter
	if (exists("level_shift", inherits = FALSE)) attr(out, "level_shift") <- level_shift
	out
}

# make a starting Sab positive definite: floor the variances and drop a cross
# moment the ALS estimates leave incompatible (starts only)
.lame_start_sab_covariance <- function(S, active = rep(TRUE, nrow(S)), floor = 1e-6) {
	S[!is.finite(S)] <- 0
	D <- diag(S)
	D[active] <- pmax(D[active], floor)
	D[!active] <- 0
	diag(S) <- D
	S[!active, ] <- 0
	S[, !active] <- 0
	idx <- which(active)
	if (length(idx) > 1L && is.null(tryCatch(chol(S[idx, idx, drop = FALSE]), error = function(e) NULL))) {
		S[row(S) != col(S)] <- 0
	}
	S
}
