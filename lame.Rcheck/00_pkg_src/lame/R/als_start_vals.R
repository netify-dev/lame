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
#' @return A list that can be passed to \code{start_vals}.
#' @export
als_start_vals <- function(fit, jitter = 0, seed = NULL) {
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
		# for probit binary the latent error variance is fixed at 1; never
		# seed an MCMC chain with an ALS working-scale s2 (which lives on
		# the IRLS working response, not the probit latent scale)
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
			S <- matrix(c(vc[["va"]], vc[["cab"]],
			              vc[["cab"]], vc[["vb"]]), 2L, 2L)
		} else if ("va" %in% names(vc)) {
			S <- matrix(vc[["va"]], 1L, 1L)
		} else {
			return(NULL)
		}
		S[!is.finite(S)] <- 0
		diag(S) <- pmax(diag(S), 1e-6)
		S
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
	U_start <- fit$U
	V_start <- fit$V %||% fit$U

	out <- list(
		beta = add_jitter(beta, jitter / 5),
		a = add_jitter(a_start),
		b = add_jitter(b_start),
		U = add_jitter(U_start),
		V = add_jitter(V_start),
		rho = clean_vec(fit$rho_path %||% fit$VC[["rho"]] %||% 0),
		s2 = get_s2(fit),
		Sab = get_sab(fit)
	)

	if (!is.null(fit$rho_uv)) out$rho_uv <- fit$rho_uv
	if (!is.null(fit$sigma_uv)) out$sigma_uv <- fit$sigma_uv
	if (!is.null(fit$rho_ab)) out$rho_ab <- fit$rho_ab
	if (!is.null(fit$sigma_ab)) out$sigma_ab <- fit$sigma_ab
	if (!is.null(fit$G)) out$G <- add_jitter(fit$G, jitter / 5)
	if (!is.null(fit$G_cube)) out$G <- add_jitter(fit$G_cube, jitter / 5)

	out <- out[!vapply(out, is.null, logical(1))]
	attr(out, "source") <- class(fit)[[1L]]
	attr(out, "jitter") <- jitter
	out
}
