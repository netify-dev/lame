# loo / waic interop for ame / lame fits via the per-iteration log_lik matrix
# stored when the user passed save_log_lik = true.

#' Read the per-iteration log-lik matrix back from on-disk chunks
#'
#' For fits run with \code{save_log_lik = "chunked"}, the per-iteration
#' pointwise log-likelihood lives in per-column-chunk binary files instead
#' of \code{fit$log_lik}. This helper reconstitutes the
#' \code{[n_stored, n_obs]} matrix in memory by issuing one \code{readBin}
#' per chunk. \code{loo.lame()} / \code{waic.lame()} call this transparently.
#'
#' @param x A fitted \code{ame}/\code{lame} object with
#'   \code{$log_lik_chunks} (i.e. fit with \code{save_log_lik = "chunked"}).
#' @return The full \code{[n_stored, n_obs]} double matrix.
#' @export
read_log_lik <- function(x) {
	ll_mat <- x[["log_lik"]]
	if (!is.null(ll_mat)) return(ll_mat)
	meta <- x[["log_lik_chunks"]]
	if (is.null(meta)) {
		cli::cli_abort(c(
			"No log-likelihood stored on the fit object.",
			"i" = "Refit with {.code save_log_lik = TRUE} or {.code save_log_lik = \"chunked\"}."))
	}
	out <- matrix(NA_real_, nrow = meta$n_stored, ncol = meta$n_obs_total)
	for (cc in seq_len(meta$n_chunks)) {
		f <- meta$files[cc]
		if (!file.exists(f)) {
			cli::cli_abort("Log-lik chunk file {.file {f}} not found; was {.var tempdir()} cleared?")
		}
		w <- meta$chunk_ends[cc] - meta$chunk_starts[cc] + 1L
		vals <- readBin(f, what = numeric(), n = meta$n_stored * w, size = 8L)
		# values were written one iteration-row's chunk-slice at a time, so
		# the on-disk layout is naturally [n_stored, w] row-major.
		out[, meta$chunk_starts[cc]:meta$chunk_ends[cc]] <- matrix(
			vals, nrow = meta$n_stored, ncol = w, byrow = TRUE)
	}
	out
}

#' Approximate leave-one-out cross-validation for AME / LAME fits
#'
#' S3 method for \code{\link[loo]{loo}} that uses the per-iteration
#' pointwise log-likelihood stored on the fit object (\code{fit$log_lik}) when
#' the model was fit with \code{save_log_lik = TRUE}. Returns the standard
#' \code{loo} object with Pareto-k diagnostics.
#'
#' \strong{What `log_lik` measures.} For \code{family} in
#' \{normal, binary, cbin, tobit, poisson, ordinal\} the stored
#' pointwise log-likelihood is the exact family-specific Y density on
#' the response scale, so \code{elpd_loo} is directly comparable to a
#' \code{loo()} output from Stan / brms fit to the same family. For
#' the rank likelihoods \code{frn} and \code{rrl} the exact marginal
#' needs GHK Monte Carlo (Halton sequence); on the longitudinal
#' \code{lame()} path you can opt in with \code{log_lik_method =
#' "observed_ghk"}, on the cross-sectional \code{ame()} path the
#' fallback is the augmented-Z normal approximation (with a one-time
#' warning). Inspect \code{fit$log_lik_method} on any fit to see which
#' branch was used.
#'
#' \strong{Chunked log-lik portability.} When fit with
#' \code{save_log_lik = "chunked"}, the on-disk chunk files default to
#' \code{tempdir()}, which is cleared at the end of the R session.
#' If you intend to \code{saveRDS()} the fit and reload it in a fresh
#' session, supply an explicit persistent \code{log_lik_path} (e.g.
#' \code{"./loglik_chunks"}) so the chunks survive the round trip.
#'
#' @param x A fitted \code{ame} or \code{lame} object that has \code{$log_lik}.
#' @param ... Additional arguments forwarded to \code{loo::loo.matrix} (e.g.
#'   \code{cores}, \code{r_eff}).
#'
#' @return A \code{loo} object.
#'
#' @examples
#' \donttest{
#' data(YX_nrm)
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 0,
#'            nscan = 200, burn = 50, odens = 5,
#'            save_log_lik = TRUE, verbose = FALSE)
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   loo_res <- loo::loo(fit)
#'   print(loo_res)
#' }
#' }
#'
#' @method loo ame
#' @export
loo.ame <- function(x, ...) {
	# use [[ to avoid partial-name matching (e.g. $log_lik -> $log_lik_method
	# when both fields are present)
	ll_mat <- x[["log_lik"]]
	ll_chunks <- x[["log_lik_chunks"]]
	if (is.null(ll_mat) && is.null(ll_chunks)) {
		cli::cli_abort(c(
			"No log-likelihood stored on the fit object.",
			"i" = "Refit with {.code save_log_lik = TRUE} to enable {.fn loo::loo}.",
			"i" = "Example: {.code lame(..., save_log_lik = TRUE)}."))
	}
	if (!requireNamespace("loo", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg loo} is required. Install with {.code install.packages(\"loo\")}.")
	}
	ll <- if (!is.null(ll_mat)) ll_mat else read_log_lik(x)
	loo::loo(ll, ...)
}

#' @rdname loo.ame
#' @method loo lame
#' @export
loo.lame <- function(x, ...) {
	loo.ame(x, ...)
}

#' WAIC for AME / LAME fits
#'
#' S3 method for \code{\link[loo]{waic}} that uses the stored \code{fit$log_lik}.
#'
#' @param x A fitted \code{ame} or \code{lame} object with \code{$log_lik}.
#' @param ... Additional arguments forwarded to \code{loo::waic.matrix}.
#'
#' @return A \code{waic} object.
#'
#' @method waic ame
#' @export
waic.ame <- function(x, ...) {
	ll_mat <- x[["log_lik"]]
	ll_chunks <- x[["log_lik_chunks"]]
	if (is.null(ll_mat) && is.null(ll_chunks)) {
		cli::cli_abort(c(
			"No log-likelihood stored on the fit object.",
			"i" = "Refit with {.code save_log_lik = TRUE} to enable {.fn loo::waic}."))
	}
	if (!requireNamespace("loo", quietly = TRUE)) {
		cli::cli_abort("Package {.pkg loo} is required.")
	}
	ll <- if (!is.null(ll_mat)) ll_mat else read_log_lik(x)
	loo::waic(ll, ...)
}

#' @rdname waic.ame
#' @method waic lame
#' @export
waic.lame <- function(x, ...) {
	waic.ame(x, ...)
}

# loo / waic for ame_als abort with an explicit message: als is a point
# estimator with no per-iteration log-likelihood, so loo::loo() / waic()
# do not apply.

#' @rdname loo.ame
#' @method loo ame_als
#' @export
loo.ame_als <- function(x, ...) {
	cli::cli_abort(c(
		"{.fn loo::loo} is not available for {.cls ame_als} fits.",
		"i" = "ALS is a point estimator and has no per-iteration log-likelihood.",
		"i" = "For LOO / WAIC, refit with {.code method = \"mcmc\"} and {.code save_log_lik = TRUE}."))
}

#' @rdname waic.ame
#' @method waic ame_als
#' @export
waic.ame_als <- function(x, ...) {
	cli::cli_abort(c(
		"{.fn loo::waic} is not available for {.cls ame_als} fits.",
		"i" = "ALS is a point estimator and has no per-iteration log-likelihood.",
		"i" = "For LOO / WAIC, refit with {.code method = \"mcmc\"} and {.code save_log_lik = TRUE}."))
}

#' Generic dispatcher for loo / waic on ame / lame fits
#'
#' Lightweight S3 generics so calls of the form \code{loo(fit)} dispatch
#' through R's S3 system even when the \pkg{loo} package is not loaded.
#' When \pkg{loo} \emph{is} loaded, its generic resolves first; this
#' fallback only fires for bare-namespace use.
#'
#' @param x a fitted \code{ame} or \code{lame} object with \code{$log_lik}.
#' @param ... passed to the relevant method.
#' @return A \code{loo} / \code{waic} object.
#' @export
loo <- function(x, ...) UseMethod("loo")

#' @rdname loo
#' @export
waic <- function(x, ...) UseMethod("waic")

# .onload helper: register the methods against loo's generic when loo is
# available, so qualified loo::loo(fit) works. parallel pattern to
# .register_posterior_methods() in zzz.r.
.register_loo_methods <- function() {
	if (requireNamespace("loo", quietly = TRUE)) {
		registerS3method("loo",  "ame",     loo.ame,     envir = asNamespace("loo"))
		registerS3method("loo",  "lame",    loo.lame,    envir = asNamespace("loo"))
		registerS3method("loo",  "ame_als", loo.ame_als, envir = asNamespace("loo"))
		registerS3method("waic", "ame",     waic.ame,    envir = asNamespace("loo"))
		registerS3method("waic", "lame",    waic.lame,   envir = asNamespace("loo"))
		registerS3method("waic", "ame_als", waic.ame_als, envir = asNamespace("loo"))
	}
}
