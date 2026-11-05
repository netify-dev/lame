# original fit call with the arguments in `...` overwritten
#' @noRd
.update_ame_call <- function(object, ...) {
	if (is.null(object$call)) {
		cli::cli_abort(c(
			"No {.code call} stored on the fit object.",
			"i" = "{.fn update} needs the original call. Refit explicitly with the new arguments."
		))
	}
	cl <- object$call
	dots <- list(...)
	if (length(dots) > 0L) {
		nms <- names(dots)
		if (is.null(nms) || any(nms == "")) {
			cli::cli_abort("All arguments to {.fn update.ame} must be named.")
		}
		for (nm in nms) cl[[nm]] <- dots[[nm]]
	}
	cl
}

# call plus the environment to evaluate it in: one holding the data snapshot
# when the fit used `freeze_call = TRUE`, otherwise the caller's frame
.update_ame_eval_setup <- function(object, dots, caller_env) {
	cl <- do.call(.update_ame_call, c(list(object = object), dots))
	if (isTRUE(object$freeze_call) && !is.null(object$data_snapshot)) {
		# point the call at the frozen data, overriding any Y, Xdyad, Xrow or
		# Xcol passed through `...`
		snap_env <- new.env(parent = caller_env)
		snap_env$.lame_frozen_Y <- object$data_snapshot$Y
		snap_env$.lame_frozen_Xdyad <- object$data_snapshot$Xdyad
		snap_env$.lame_frozen_Xrow <- object$data_snapshot$Xrow
		snap_env$.lame_frozen_Xcol <- object$data_snapshot$Xcol
		cl$Y <- quote(.lame_frozen_Y)
		if (!is.null(object$data_snapshot$Xdyad)) cl$Xdyad <- quote(.lame_frozen_Xdyad)
		if (!is.null(object$data_snapshot$Xrow)) cl$Xrow <- quote(.lame_frozen_Xrow)
		if (!is.null(object$data_snapshot$Xcol)) cl$Xcol <- quote(.lame_frozen_Xcol)
		list(cl = cl, envir = snap_env)
	} else {
		list(cl = cl, envir = caller_env)
	}
}

#' Update an AME / LAME fit
#'
#' S3 method for \code{\link[stats]{update}} that re-fits the model with
#' modified arguments. Reuses the original call recorded in \code{fit$call}.
#' If the fit was produced with \code{freeze_call = TRUE}, the data snapshot
#' on \code{fit$data_snapshot} is used in place of looking up names in the
#' caller's environment.
#'
#' @param object A fitted \code{ame} or \code{lame} object.
#' @param ... Named arguments to overwrite in the original call (e.g.
#'   \code{nscan = 5000}, \code{dynamic_beta = TRUE}, \code{R = 2}).
#' @param evaluate Logical: if \code{TRUE} (default), evaluate the updated
#'   call and return the new fit; if \code{FALSE}, return the unevaluated
#'   call.
#'
#' @return A new fitted object, or (when \code{evaluate = FALSE}) the
#'   modified call.
#'
#' @examples
#' \donttest{
#' data(YX_bin_list)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- lame(lapply(YX_bin_list$Y, function(y) 1 * (y > 0)), YX_bin_list$X, family = "binary",
#'             burn = 5, nscan = 20, odens = 1, verbose = FALSE, chains = 1)
#' # toggle dynamic_beta on without rewriting the whole call
#' fit_dyn <- update(fit, dynamic_beta = "dyad")
#' dim(fit_dyn$BETA)  # 3-D now
#' }
#'
#' @method update ame
#' @export
update.ame <- function(object, ..., evaluate = TRUE) {
	setup <- .update_ame_eval_setup(object, list(...), parent.frame())
	if (!isTRUE(evaluate)) return(setup$cl)
	eval(setup$cl, envir = setup$envir)
}

#' @rdname update.ame
#' @method update lame
#' @export
update.lame <- function(object, ..., evaluate = TRUE) {
	# called inline so parent.frame() is the user's frame and names such as
	# `dat$Y` resolve
	setup <- .update_ame_eval_setup(object, list(...), parent.frame())
	if (!isTRUE(evaluate)) return(setup$cl)
	eval(setup$cl, envir = setup$envir)
}
