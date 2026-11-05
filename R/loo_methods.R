# loo and waic methods for ame and lame fits, built on the per-iteration
# log_lik matrix stored with save_log_lik = TRUE

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
#' @examples
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, nscan = 20, burn = 5, odens = 1,
#'            save_log_lik = TRUE, verbose = FALSE, chains = 1)
#' ll <- read_log_lik(fit)
#' dim(ll)
#' @export
read_log_lik <- function(x) {
	ll_mat <- x[["log_lik"]]
	if (!is.null(ll_mat)) return(ll_mat)
	meta <- x[["log_lik_chunks"]]
	if (is.null(meta)) {
		cli::cli_abort(c(
			"No log-likelihood stored on the fit object.",
			"i" = "Refit with {.code save_log_lik = TRUE} or {.code save_log_lik = \"chunked\"}."
		))
	}
	out <- matrix(NA_real_, nrow = meta$n_stored, ncol = meta$n_obs_total)
	for (cc in seq_len(meta$n_chunks)) {
		f <- meta$files[cc]
		if (!file.exists(f)) {
			cli::cli_abort("Log-lik chunk file {.file {f}} not found; was {.var tempdir()} cleared?")
		}
		w <- meta$chunk_ends[cc] - meta$chunk_starts[cc] + 1L
		vals <- readBin(f, what = numeric(), n = meta$n_stored * w, size = 8L)
		# chunks were written one iteration row at a time, so the on-disk
		# layout is row-major [n_stored, w]
		out[, meta$chunk_starts[cc]:meta$chunk_ends[cc]] <- matrix(
			vals,
			nrow = meta$n_stored,
			ncol = w,
			byrow = TRUE
		)
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
#' \strong{Observed likelihood factors.} LOO and WAIC require a current
#' \code{fit$log_lik_contract} verifying observed-data conditional likelihood
#' factors. Latent augmented scores, approximate rank/cap diagnostics and
#' legacy fits without this contract are rejected. Directed cell marginals
#' do not factorize under nonzero or estimated residual dyad correlation;
#' joint dyad factors would be needed. Symmetric fits store each unordered
#' dyad once. Correct marginal predictive scores remain useful and readable
#' through \code{read_log_lik()} even when LOO or WAIC is unavailable.
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
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 0, dcor = FALSE,
#'            nscan = 60, burn = 15, odens = 5,
#'            save_log_lik = TRUE, verbose = FALSE, chains = 1)
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   loo_res <- loo::loo(fit)
#'   print(loo_res)
#' }
#' }
#'
#' @method loo ame
#' @export
loo.ame <- function(x, ...) {
	.lame_require_draws(x, "loo")
	# [[ avoids partial matching of $log_lik to $log_lik_method
	ll_mat <- x[["log_lik"]]
	ll_chunks <- x[["log_lik_chunks"]]
	if (is.null(ll_mat) && is.null(ll_chunks)) {
		cli::cli_abort(c(
			"No log-likelihood stored on the fit object.",
			"i" = "Refit with {.code save_log_lik = TRUE} to enable {.fn loo::loo}.",
			"i" = "Example: {.code lame(..., save_log_lik = TRUE)}."
		))
	}
	.lame_require_loo_contract(x)
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
	.lame_require_draws(x, "loo")
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
	.lame_require_draws(x, "waic")
	ll_mat <- x[["log_lik"]]
	ll_chunks <- x[["log_lik_chunks"]]
	if (is.null(ll_mat) && is.null(ll_chunks)) {
		cli::cli_abort(c(
			"No log-likelihood stored on the fit object.",
			"i" = "Refit with {.code save_log_lik = TRUE} to enable {.fn loo::waic}."
		))
	}
	.lame_require_loo_contract(x)
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
	.lame_require_draws(x, "waic")
	waic.ame(x, ...)
}

# ame_als fits are point estimates with no per-iteration log-likelihood, so
# loo and waic abort

#' @rdname loo.ame
#' @method loo ame_als
#' @export
loo.ame_als <- function(x, ...) {
	cli::cli_abort(c(
		"{.fn loo::loo} is not available for {.cls ame_als} fits.",
		"i" = "ALS is a point estimator and has no per-iteration log-likelihood.",
		"i" = "For LOO / WAIC, refit with {.code method = \"mcmc\"} and {.code save_log_lik = TRUE}."
	))
}

#' @rdname waic.ame
#' @method waic ame_als
#' @export
waic.ame_als <- function(x, ...) {
	cli::cli_abort(c(
		"{.fn loo::waic} is not available for {.cls ame_als} fits.",
		"i" = "ALS is a point estimator and has no per-iteration log-likelihood.",
		"i" = "For LOO / WAIC, refit with {.code method = \"mcmc\"} and {.code save_log_lik = TRUE}."
	))
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
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 12
#' Y <- matrix(rnorm(n * n), n, n); diag(Y) <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y, R = 1, dcor = FALSE, nscan = 100, burn = 5, odens = 1,
#'            save_log_lik = TRUE, verbose = FALSE, chains = 1)
#' if (requireNamespace("loo", quietly = TRUE)) {
#'   # a fit this small can trip loo's own low-sample Pareto-k diagnostic;
#'   # that diagnostic warning is expected here and not specific to lame
#'   print(suppressWarnings(loo(fit)))
#' }
#' }
#' @export
loo <- function(x, ...) UseMethod("loo")

#' @rdname loo
#' @export
waic <- function(x, ...) UseMethod("waic")

# register the methods on loo's generics at load time so loo::loo(fit) works
# (same pattern as .register_posterior_methods() in zzz.R)
.register_loo_methods <- function() {
	if (requireNamespace("loo", quietly = TRUE)) {
		registerS3method("loo", "ame", loo.ame, envir = asNamespace("loo"))
		registerS3method("loo", "lame", loo.lame, envir = asNamespace("loo"))
		registerS3method("loo", "ame_als", loo.ame_als, envir = asNamespace("loo"))
		registerS3method("waic", "ame", waic.ame, envir = asNamespace("loo"))
		registerS3method("waic", "lame", waic.lame, envir = asNamespace("loo"))
		registerS3method("waic", "ame_als", waic.ame_als, envir = asNamespace("loo"))
	}
}
