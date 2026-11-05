# safe evaluation of user-supplied `custom_gof` functions; failures are counted and reported once
# after fitting

#' @noRd
.new_custom_gof_err_tracker <- function() {
	env <- new.env(parent = emptyenv())
	env$count <- 0L
	env$first_msg <- NULL
	env
}

#' @noRd
.record_custom_gof_err <- function(tracker, e) {
	tracker$count <- tracker$count + 1L
	if (is.null(tracker$first_msg)) {
		tracker$first_msg <- conditionMessage(e)
	}
}

# warn at most once per fit location so a broken custom_gof does not flood the console
#' @noRd
.maybe_warn_custom_gof <- function(tracker, where = "MCMC") {
	if (!is.environment(tracker) || isTRUE(tracker$count == 0L)) return(invisible(NULL))
	cli::cli_warn(c(
		"!" = "{.arg custom_gof} failed on {.val {tracker$count}} {where} call{?s}; affected GOF entries set to {.val NA}.",
		"i" = "First error: {tracker$first_msg}",
		"i" = "Check your {.arg custom_gof} function: it must accept a single matrix {.arg Y} and return a numeric scalar (or named numeric vector for multi-statistic mode)."
	))
	invisible(NULL)
}
