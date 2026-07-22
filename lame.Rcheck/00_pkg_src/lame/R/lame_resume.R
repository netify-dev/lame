# resume a lame() run from the checkpoint written by
# lame(..., checkpoint_path = ...). uses the saved beta / vc means as
# start values and re-invokes lame() with the original call args.

#' Resume a `lame()` MCMC run from a checkpoint
#'
#' For a fit started with \code{lame(..., checkpoint_path = "X.rds")}
#' that terminated early (either via \code{max_seconds} or an external
#' interruption), continues the chain from the most recent checkpoint.
#' The implementation is pragmatic: it loads the saved RNG state and
#' re-invokes \code{lame()} with the original call arguments. The
#' result is a fresh fit that picks up where the previous one left
#' off in the random-number stream; the underlying MCMC counter
#' restarts at 1.
#'
#' \strong{Equivalent consolidated entry point.}
#' \code{lame(resume_from = path, ...)} short-circuits to this
#' function with the user-supplied overrides forwarded; pass
#' \code{nscan = K} on the resume call to request \code{K}
#' additional stored draws. Both call shapes are supported. The
#' \code{lame_resume(path, ...)} form is safer for nested calls: the
#' consolidated form re-evaluates the
#' saved \code{lame()} call in \code{parent.frame()}, which is
#' the \code{lame()} frame whose required formal arguments
#' (\code{Y}, \code{Xdyad}, etc.) were not supplied on the
#' resume call. Prefer \code{lame_resume(path, ...)} when calling
#' from inside other functions or from non-global scopes.
#'
#' Use \code{nscan_more = K} to override the \code{nscan} value to
#' \code{K} for the continuation. Other arguments can be overridden by
#' passing them to \code{...}.
#'
#' @param path Checkpoint file path (the one passed as
#'   \code{checkpoint_path} to \code{lame()}).
#' @param nscan_more Optional integer, the new \code{nscan} for the
#'   continuation. If \code{NULL}, the original \code{nscan} is used.
#' @param ... Additional arguments forwarded to \code{lame()} (override
#'   the saved values).
#' @param .envir Environment in which to evaluate the saved call's data
#'   arguments. Defaults to the caller's frame; used internally when
#'   \code{lame()} forwards a resume.
#'
#' @return A fitted \code{lame} object.
#'
#' @examples
#' \donttest{
#' ck <- tempfile(fileext = ".rds")
#' data(YX_bin_list)
#' # short run with very-aggressive max_seconds to force early termination
#' fit1 <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
#'              nscan = 5000, burn = 50, odens = 5,
#'              checkpoint_path = ck, checkpoint_every = 50L,
#'              max_seconds = 0.5, verbose = FALSE)
#' if (isTRUE(fit1$terminated_early)) {
#'   fit2 <- lame_resume(ck, nscan_more = 200)
#'   dim(fit2$BETA)
#' }
#' }
#'
#' @export
lame_resume <- function(path, nscan_more = NULL, ..., .envir = NULL) {
	if (!file.exists(path)) {
		cli::cli_abort("Checkpoint file {.file {path}} not found.")
	}
	ck <- readRDS(path)
	if (!is.list(ck) || is.null(ck$call_args)) {
		cli::cli_abort("Checkpoint at {.file {path}} is malformed (no {.code call_args}).")
	}
	cl <- ck$call_args

	# user overrides
	usr <- list(...)
	if (!is.null(nscan_more)) usr$nscan <- as.integer(nscan_more)
	# patch the call: first scrub the checkpoint args from the saved call,
	# then layer user overrides on top, so user-supplied checkpoint_path /
	# max_seconds on the resume call take effect.
	cl2 <- cl
	cl2$checkpoint_path <- NULL
	cl2$max_seconds <- NULL
	for (nm in names(usr)) cl2[[nm]] <- usr[[nm]]
	# the original run paid the burn-in cost; the continuation re-uses the
	# saved rng stream and skips burn-in unless the user overrides it.
	if (!"burn" %in% names(usr)) cl2$burn <- 0L
	# hand the checkpoint's rng stream to the continuation, then put the
	# caller's stream back so resuming leaves their session untouched
	if (!is.null(ck$rng_state)) {
		.had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
		.old_seed <- if (.had_seed) get(".Random.seed", envir = globalenv()) else NULL
		on.exit({
			if (.had_seed) {
				assign(".Random.seed", .old_seed, envir = globalenv())
			} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
				rm(".Random.seed", envir = globalenv())
			}
		}, add = TRUE)
		assign(".Random.seed", ck$rng_state, envir = globalenv())
	}
	# evaluate in the user's frame. when called directly by a user from the
	# top level, parent.frame() is globalenv() and y/xdyad references in the
	# saved call resolve fine. when called via `do.call()` from the
	# consolidated `lame(resume_from = ...)` short-circuit, the lame() frame
	# has y as a missing promise; in that case the caller must pass `.envir`
	# explicitly (typically the user's frame captured at lame() entry).
	use_env <- if (!is.null(.envir)) .envir else parent.frame()
	eval(cl2, envir = use_env)
}
