# internal helper: extract the original lame() call and overwrite the named
# arguments passed via `...`. not exported.
#' @noRd
.update_ame_call <- function(object, ...) {
	if (is.null(object$call)) {
		cli::cli_abort(c(
			"No {.code call} stored on the fit object.",
			"i" = "{.fn update} needs the original call. Refit explicitly with the new arguments."))
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

# build the call + an environment whose enclosing frame holds the
# data snapshot (if `freeze_call = true` was set at fit time). returns a list
# (cl, envir) ready for eval(). when no snapshot is present, `envir` is the
# caller's frame (parent.frame() of the s3 method dispatcher).
.update_ame_eval_setup <- function(object, dots, caller_env) {
	cl <- do.call(.update_ame_call, c(list(object = object), dots))
	if (isTRUE(object$freeze_call) && !is.null(object$data_snapshot)) {
		# substitute frozen data references in the call. the snapshot
		# overrides any y/xdyad/xrow/xcol the user passes through `...`
		# unless they explicitly thaw with `freeze_call = false`.
		snap_env <- new.env(parent = caller_env)
		snap_env$.lame_frozen_Y     <- object$data_snapshot$Y
		snap_env$.lame_frozen_Xdyad <- object$data_snapshot$Xdyad
		snap_env$.lame_frozen_Xrow  <- object$data_snapshot$Xrow
		snap_env$.lame_frozen_Xcol  <- object$data_snapshot$Xcol
		cl$Y     <- quote(.lame_frozen_Y)
		if (!is.null(object$data_snapshot$Xdyad)) cl$Xdyad <- quote(.lame_frozen_Xdyad)
		if (!is.null(object$data_snapshot$Xrow))  cl$Xrow  <- quote(.lame_frozen_Xrow)
		if (!is.null(object$data_snapshot$Xcol))  cl$Xcol  <- quote(.lame_frozen_Xcol)
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
#' fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary",
#'             burn = 5, nscan = 20, odens = 1, verbose = FALSE)
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
	# evaluate in the snapshot env if freeze_call was set; otherwise in the
	# caller's frame so user-supplied data objects (e.g. `dat$y`) are visible.
	eval(setup$cl, envir = setup$envir)
}

#' @rdname update.ame
#' @method update lame
#' @export
update.lame <- function(object, ..., evaluate = TRUE) {
	# inline (not a method call) so parent.frame() inside the eval()
	# resolves to the user's calling frame and references like `dat$y`
	# stay in scope.
	setup <- .update_ame_eval_setup(object, list(...), parent.frame())
	if (!isTRUE(evaluate)) return(setup$cl)
	eval(setup$cl, envir = setup$envir)
}
