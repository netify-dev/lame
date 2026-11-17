# the factor by which the dynamic fit raises the sweep limit of a normal fit once (lame_dynamic_als())
.ae_budget_factor <- 4L

# numerical record of one ALS call: the parts of the fit that could not be
# computed (with the first error of each), and the counted exact fallbacks.
# the outermost entry point opens the record and attaches it to its fit
.ae_notes <- new.env(parent = emptyenv())

.ae_notes_open <- function() {
	if (isTRUE(.ae_notes$open)) return(FALSE)
	.ae_notes$open <- TRUE
	.ae_notes$failures <- list()
	.ae_notes$fallbacks <- numeric(0)
	TRUE
}

.ae_notes_close <- function() {
	.ae_notes$open <- FALSE
	invisible(NULL)
}

# a part of the fit left out after an error: `what` names it for the user
.ae_note_failure <- function(site, what, e) {
	if (!isTRUE(.ae_notes$open)) return(invisible(NULL))
	if (is.null(.ae_notes$failures[[site]])) {
		msg <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
		.ae_notes$failures[[site]] <- list(what = what, message = msg)
	}
	invisible(NULL)
}

# a step that took an exact alternative or a valid default
.ae_note_fallback <- function(site) {
	if (!isTRUE(.ae_notes$open)) return(invisible(NULL))
	counts <- .ae_notes$fallbacks
	counts[site] <- (if (site %in% names(counts)) counts[[site]] else 0) + 1
	.ae_notes$fallbacks <- counts
	invisible(NULL)
}

# attach the record to the fit of the outermost call and warn for each part left out
.ae_notes_attach <- function(fit, top) {
	if (!isTRUE(top)) return(fit)
	failures <- .ae_notes$failures %||% list()
	fallbacks <- .ae_notes$fallbacks %||% numeric(0)
	.ae_notes_close()
	if (!is.list(fit)) return(fit)
	if (length(fallbacks)) fit$exact_fallbacks <- fallbacks
	if (length(failures)) {
		fit$numerical_failures <- failures
		lines <- vapply(failures, function(f) {
			m <- gsub("}", "}}", gsub("{", "{{", f$message, fixed = TRUE), fixed = TRUE)
			paste0(f$what, ": ", m)
		}, character(1))
		cli::cli_warn(
			c(
				"Parts of the ALS fit could not be computed and are left out of it:",
				setNames(lines, rep("x", length(lines))),
				"i" = "They are listed in {.code fit$numerical_failures}."
			),
			class = "lame_als_numerical_fallback",
			failures = failures
		)
	}
	fit
}

# evaluate a refit (a bootstrap or bias-correction replicate, a stability restart)
# outside the record of the call around it: a refit keeps a record of its own
.ae_notes_isolated <- function(expr) {
	saved <- list(open = .ae_notes$open, failures = .ae_notes$failures, fallbacks = .ae_notes$fallbacks)
	.ae_notes$open <- FALSE
	on.exit({
		.ae_notes$open <- saved$open
		.ae_notes$failures <- saved$failures
		.ae_notes$fallbacks <- saved$fallbacks
	})
	expr
}
