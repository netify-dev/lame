# resume a lame() run from its checkpoint (sampler state, seed and original call), finishing
# any burn-in and advancing the seed by the completed iterations

#' Resume a `lame()` MCMC run from a checkpoint
#'
#' For a fit started with \code{lame(..., checkpoint_path = "X.rds")}
#' that terminated early (either via \code{max_seconds} or an external
#' interruption), continues the chain from the most recent checkpoint.
#' The checkpoint holds the full sampler state (latent matrix,
#' coefficients, additive and multiplicative effects, variance
#' components, dynamic coefficient paths, auxiliary states and adapted
#' proposal scales) together with the seed
#' and the original call. The continuation re-invokes \code{lame()} with
#' that state as \code{start_vals}, completing any unfinished burn-in,
#' and with \code{seed} advanced by the number of
#' completed iterations, so it does not replay the first run's random
#' stream. The result is a fresh fit whose stored draws start where the
#' checkpointed chain stopped; draws stored before the interruption are
#' not carried over. \code{nscan_more} requests additional production sweeps;
#' the saved thinning setting \code{odens} determines how many draws are retained.
#'
#' \strong{Equivalent consolidated entry point.}
#' \code{lame(resume_from = path, ...)} short-circuits to this
#' function with the user-supplied overrides forwarded; pass
#' \code{nscan = K} on the resume call to request \code{K}
#' additional post-burn sweeps. Both call shapes evaluate saved data arguments
#' in the caller's frame, so objects referenced by the original call must
#' remain available there.
#'
#' Use \code{nscan_more = K} to override the \code{nscan} value to
#' \code{K} for the continuation. A fit continued from the object runs with the
#' arguments of its call, \code{max_seconds} included: a fit that stopped at its
#' time limit stops again at the same limit unless \code{max_seconds = Inf} is passed.
#' The continuation holds \code{nscan_more / odens} draws per chain
#' (\code{chains * nscan_more / odens} pooled) when it runs to the end, plus one
#' block of each extension. Other arguments can be overridden by
#' passing them to \code{...}. Older checkpoints without the complete
#' transition state are treated as warm starts with a fresh burn-in and a
#' warning. A \code{start_vals}, \code{burn} or
#' \code{seed} passed there replaces the checkpoint's choice.
#'
#' @param path Checkpoint RDS file path, as written through
#' \code{checkpoint_path}, or a fitted \code{lame} / \code{ame} object
#' (for instance the \code{fit} carried by a \code{lame_chain_unsettled}
#' error). A fit is continued from the state it stopped in, with the call that
#' produced it, \code{burn = 0} and the seed advanced by the sweeps it ran.
#' For a pooled multi-chain fit every chain is continued from its own final state,
#' and the chains are checked and compared again (see \code{\link{ame}}); a single
#' chain of a refused call is continued with \code{lame_resume(e$fits[[k]])}.
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
#' fit1 <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X,
#'              family = "binary", R = 0, nscan = 5000, burn = 50, odens = 5,
#'              checkpoint_path = ck, checkpoint_every = 50L,
#'              max_seconds = 0.5, verbose = FALSE)
#' if (isTRUE(fit1$terminated_early)) {
#'   # each chain wrote its own checkpoint (ck with _chain1, _chain2 before
#'   # the extension): continue one chain from its file, or every chain
#'   # from the returned fit
#'   fit2 <- lame_resume(sub(".rds", "_chain1.rds", ck, fixed = TRUE), nscan_more = 200)
#'   dim(fit2$BETA)
#'   # the time limit of the call applies again unless it is lifted
#'   fit3 <- lame_resume(fit1, nscan_more = 200, max_seconds = Inf)
#'   nrow(fit3$BETA)  # 2 chains * 200 / 5 = 80
#' }
#' }
#'
#' @details
#' Tempered checkpoints restore every replica and the complete primary and
#' auxiliary RNG streams, completing any unfinished original burn-in. The
#' returned fit contains only newly retained draws, which exactly continue an
#' uninterrupted run. Data, priors, seed, thinning and sampler changes are
#' refused; start a new fit to change the model. Output settings, time limits
#' and the number of additional production sweeps may be changed.
#'
#' Checkpoints retain resolved priors independently of the original call.
#' Where the checkpoint type permits a \code{prior} override, it can be supplied
#' explicitly. Otherwise a continuation uses those
#' settings even if a caller variable or a package default has changed.
#' The original user settings are retained separately for training-window refits.
#' Complete older Gaussian bipartite AR(1) checkpoints without frozen priors
#' retain their former persistence defaults (mean 0.9, SD 0.1) when the saved
#' call did not supply them.
#'
#' @export
lame_resume <- function(path, nscan_more = NULL, ..., .envir = NULL) {
	if (inherits(path, "ame")) {
		return(.lame_resume_fit(path, nscan_more, list(...), if (is.null(.envir)) parent.frame() else .envir))
	}
	if (!file.exists(path)) {
		# a multi-chain run writes one checkpoint per chain (.lame_chain_path)
		per_chain <- vapply(seq_len(64L), function(k) .lame_chain_path(path, k), "")
		per_chain <- per_chain[file.exists(per_chain)]
		if (length(per_chain)) {
			cli::cli_abort(
				c(
					"Checkpoint file {.file {path}} not found; the run wrote one checkpoint per chain.",
					"i" = "Resume one chain from its file ({.file {per_chain}}), or every chain from the returned fit with {.code lame_resume(fit, nscan_more = )}."
				),
				class = "lame_resume_chain_files",
				files = per_chain
			)
		}
		cli::cli_abort("Checkpoint file {.file {path}} not found.")
	}
	ck <- readRDS(path)
	if (is.list(ck) && identical(ck$type, "lame_tempered"))
		return(.lame_tempered_resume(ck, nscan_more, list(...)))
	if (!is.list(ck) || is.null(ck$call_args)) {
		cli::cli_abort("Checkpoint at {.file {path}} is malformed (no {.code call_args}).")
	}
	cl <- ck$call_args

	usr <- list(...)
	if (!is.null(nscan_more)) usr$nscan <- as.integer(nscan_more)
	# drop the checkpoint arguments from the saved call, then apply user overrides
	cl2 <- cl
	cl2$checkpoint_path <- NULL
	cl2$max_seconds <- NULL
	cl2$resume_from <- NULL
	# a checkpoint is one chain
	cl2$chains <- 1L
	for (nm in names(usr)) cl2[[nm]] <- usr[[nm]]
	# an interrupted burn-in is completed before storing new draws
	if (!"burn" %in% names(usr)) {
		remaining <- ck[["burn_remaining"]]
		if (is.numeric(remaining) && length(remaining) == 1L &&
				is.finite(remaining) && remaining >= 0) {
			cl2$burn <- as.integer(remaining)
		} else {
			cl2$burn <- 0L
		}
	}
	if (is.null(ck$start_vals[["checkpoint_state"]]) && !"start_vals" %in% names(usr)) {
		if (!"burn" %in% names(usr)) {
			cl2$burn <- if (!is.null(cl$burn)) cl$burn else formals(lame)$burn
		}
		cli::cli_warn(
			c(
				"This older checkpoint lacks the complete transition state.",
				"i" = "It is used as a warm start with a fresh burn-in; supply {.arg burn} to choose its length."
			),
			class = "lame_checkpoint_legacy"
		)
	}
	# advance the saved seed by the completed iterations so the continuation does not replay the first run
	if (!"seed" %in% names(usr)) {
		base_seed <- if (is.numeric(ck$seed) && length(ck$seed) == 1L &&
				is.finite(ck$seed)) ck$seed else 6886
		iter_done <- if (is.numeric(ck$iter) && length(ck$iter) == 1L &&
				is.finite(ck$iter)) ck$iter else 0
		cl2$seed <- as.integer(base_seed + iter_done)
	}
	# evaluate in the user's frame; a call through lame(resume_from = ) must pass `.envir`, since Y is
	# a missing promise there
	use_env <- if (!is.null(.envir)) .envir else parent.frame()
	if (!"prior" %in% names(usr)) {
		resume_prior <- .lame_checkpoint_prior(ck, use_env)
		if (!is.null(resume_prior)) cl2$prior <- resume_prior
	}
	# the saved state enters start_vals through a binding in a child frame, so fit$call holds a name,
	# not the state
	eval_env <- new.env(parent = use_env)
	if (!"start_vals" %in% names(usr)) {
		if (is.list(ck$start_vals)) {
			assign(".lame_resume_state", ck$start_vals, envir = eval_env)
			cl2$start_vals <- as.name(".lame_resume_state")
		} else {
			cli::cli_warn(c(
				"Checkpoint at {.file {path}} carries no sampler state; the continuation starts from fresh start values.",
				"i" = "Pass {.arg start_vals} to seed it yourself."
			))
		}
	}
	eval(cl2, envir = eval_env)
}

# continue a fitted object from its final state; a fit refused by the chain policy carries its
# evaluated arguments in fit$auto_extend$call_args, any other fit re-evaluates fit$call in `env`
.lame_resume_fit <- function(fit, nscan_more, usr, env) {
	if (isTRUE(fit$is_combined) && !is.null(fit$chain_start_vals)) return(
		.lame_resume_chains(fit, nscan_more, usr, env)
	)
	if (is.null(fit$start_vals)) {
		cli::cli_abort("The fit carries no sampler state ({.code fit$start_vals}) to continue from.")
	}
	rec <- fit$auto_extend
	is_lame <- inherits(fit, "lame")
	args <- rec$call_args
	if (is.null(args)) {
		if (is.null(fit$call)) cli::cli_abort("The fit carries neither its call nor its arguments.")
		cl <- as.list(fit$call)[-1L]
		args <- lapply(cl, eval, envir = env)
	}
	args$chains <- 1L
	args$start_vals <- fit$start_vals
	args$burn <- 0L
	if (!is.null(nscan_more)) args$nscan <- as.integer(nscan_more)
	seed0 <- args$seed %||% 6886
	args$seed <- if (is.list(rec) && is.finite(rec$next_seed)) as.integer(rec$next_seed) else
		as.integer(seed0 + (fit$refit_options$burn %||% 0) + (fit$refit_options$nscan %||% 0))
	if (is_lame && is.null(usr$prior)) {
		pr <- fit$prior
		attr(pr, "lame_checkpoint_input_prior") <- fit$refit_options$prior %||% args$prior %||% list()
		args$prior <- pr
	}
	for (nm in names(usr)) args[nm] <- list(usr[[nm]])
	if ("max_seconds" %in% names(usr) && is.null(usr$max_seconds)) args$max_seconds <- Inf
	do.call(if (is_lame) lame else ame, args, envir = env)
}

# continue every chain of a pooled fit from its own final state, checked and extended together
.lame_resume_chains <- function(fit, nscan_more, usr, env) {
	rec <- fit$auto_extend
	is_lame <- inherits(fit, "lame")
	args <- rec$call_args
	if (is.null(args)) {
		if (is.null(fit$call)) cli::cli_abort("The fit carries neither its call nor its arguments.")
		cl <- as.list(fit$call)[-1L]
		args <- lapply(cl, eval, envir = env)
	}
	for (nm in c("chains", "n_chains", "cores", "start_vals", "burn", "seed")) args[[nm]] <- NULL
	if (!is.null(nscan_more)) args$nscan <- as.integer(nscan_more)
	for (nm in names(usr)) args[nm] <- list(usr[[nm]])
	if ("max_seconds" %in% names(usr) && is.null(usr$max_seconds)) args$max_seconds <- Inf
	chains <- fit$chains
	seeds <- if (!is.null(rec$next_seeds)) rec$next_seeds else fit$chain_seeds + 1000 * 0
	resolved <- if (is_lame) fit$prior else NULL
	input_prior <- if (is_lame) (fit$refit_options$prior %||% args$prior %||% list()) else NULL
	args$verbose <- isTRUE(args$verbose)
	.lame_chains_run(
		if (is_lame) "lame" else "ame",
		args,
		chains,
		as.numeric(seeds),
		fit$call,
		cores = NULL,
		caller = env,
		start_list = fit$chain_start_vals,
		resolved = resolved,
		input_prior = input_prior
	)
}
