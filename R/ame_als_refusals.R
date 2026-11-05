# refusal of an unconverged fixed-effects or normal-family ALS fit, and the
# single budget escalation before it. a fit stopped on its iteration limit is
# refitted once with a larger limit; a fit where no step lowers the objective
# is refused at once. the last iterate is never returned

# the factor by which the iteration limit is raised for the one escalated run
.ae_budget_factor <- 4L

# the refusal. `kind` is "irls" or "bcd"; `state` holds the termination,
# scores and iterations; `from` is the starting limit when this is the
# escalated run (NULL otherwise)
.ae_abort_unconverged <- function(kind, state, max_iter, from = NULL) {
	fmt <- function(x) format(signif(x, 3))
	retry <- FALSE
	if (identical(kind, "irls")) {
		term <- state$termination
		head <- "IRLS did not reach its likelihood and variance convergence criteria within {.arg max_iter} = {max_iter}."
		detail <- if (is.null(state$reason)) {
			sprintf(
				"It stopped on %s after %d accepted reweighting steps: the largest scaled likelihood score was %s against a threshold of %s, and the largest change of a log factor scale was %s.",
				term,
				as.integer(state$outer_iterations),
				fmt(state$max_score),
				fmt(state$threshold),
				fmt(state$scale_change)
			)
		} else {
			sprintf("It stopped on %s.", state$reason)
		}
		limit <- identical(term, "iteration limit")
		remedy <- if (limit) {
			"Increase {.arg max_iter}: the iteration limit, not a failed step, is what stopped the fit."
		} else {
			"The reweighted solve could not lower the penalised deviance, so a larger {.arg max_iter} does not help: lower {.arg R}, or use the default empirical-Bayes estimator ({.code bias_correct = NULL})."
		}
		retry <- limit
	} else {
		gd <- state
		head <- "Block coordinate descent did not converge within {.arg max_iter} = {max_iter} iterations."
		detail <- if (is.null(gd$max_score)) {
			"It stopped before the first-order criterion could be evaluated."
		} else {
			sprintf(
				"The largest scaled first-order score was %s against a threshold of %s after %d sweeps%s.",
				fmt(gd$max_score),
				fmt(gd$threshold),
				as.integer(gd$iterations),
				if (length(gd$refinement)) sprintf(
					" and %d joint refinement round%s",
					length(gd$refinement),
					if (length(gd$refinement) > 1L) "s" else ""
				) else ""
			)
		}
		remedy <- "Increase {.arg max_iter}; for {.code R > 0} the objective is non-convex and a flat one converges slowly, so a lower {.arg R} or {.code multistart = \"cheap\"} may also finish it."
		retry <- TRUE
	}
	msg <- c(head, "x" = detail)
	if (!is.null(from)) {
		msg <- c(
			msg,
			"i" = sprintf(
				"The limit was raised from %d to %d for one more run before the fit was refused.",
				as.integer(from),
				as.integer(max_iter)
			)
		)
	}
	msg <- c(msg, "i" = remedy)
	cli::cli_abort(
		msg,
		class = "lame_als_convergence",
		kind = kind,
		state = state,
		max_iter = max_iter,
		retry = retry && is.null(from)
	)
}

# shared driver of ame_als() and lame_als() with the budget escalation;
# warnings and messages the second run repeats are shown once
.ae_run_budget <- function(
	Y,
	Xdyad,
	Xrow,
	Xcol,
	R,
	family,
	mode,
	symmetric,
	max_iter,
	tol,
	verbose,
	seed,
	longitudinal,
	call,
	...
) {
	run <- function(mi, vb = verbose) {
		.ae_run(
			Y,
			Xdyad,
			Xrow,
			Xcol,
			R,
			family,
			mode,
			symmetric,
			mi,
			tol,
			vb,
			seed,
			longitudinal,
			call,
			...
		)
	}
	# log marginal likelihood of the same fit at a smaller rank r, reported by a
	# fit whose product is withheld (.ae_withhold_product); the prior-scale
	# uncertainty fits are skipped since they do not change the mode's evidence
	lower <- function(r) {
		op <- options(lame.eb_scale_uncertainty = FALSE)
		on.exit(options(op), add = TRUE)
		at <- function(post) {
			suppressWarnings(suppressMessages(.ae_run(
				Y,
				Xdyad,
				Xrow,
				Xcol,
				r,
				family,
				mode,
				symmetric,
				max_iter,
				tol,
				FALSE,
				seed,
				longitudinal,
				call,
				...,
				eb_posterior = post
			)))
		}
		# a converged binary or count mode gives the Laplace evidence directly; normal
		# fits and unconverged modes take it from the variational refinement
		fm <- tryCatch(at("mode"), error = function(e) NULL)
		lm <- .ae_default(fm$eb$log_marginal, NA_real_)
		if (!is.null(fm) && !identical(fm$family, "normal") && isTRUE(fm$converged) && is.finite(lm)) return(lm)
		f <- at("variational")
		.ae_default(f$eb$log_marginal, NA_real_)
	}
	seen <- character()
	first <- tryCatch(
		withCallingHandlers(
			run(max_iter),
			warning = function(w) seen <<- c(seen, conditionMessage(w)),
			message = function(m) seen <<- c(seen, conditionMessage(m))
		),
		error = function(e) e
	)
	if (!inherits(first, "error")) return(.ae_withhold_product(first, lower))
	if (!inherits(first, "lame_als_convergence") || !isTRUE(first$retry)) {
		stop(first)
	}
	budget <- as.integer(min(.Machine$integer.max, .ae_budget_factor * max_iter))
	if (isTRUE(verbose)) {
		cli::cli_alert_info("Not converged in {max_iter} iterations; refitting once with {budget}.")
	}
	out <- withCallingHandlers(
		tryCatch(run(budget, vb = FALSE), error = function(e) {
			if (inherits(e, "lame_als_convergence")) {
				.ae_abort_unconverged(e$kind, e$state, e$max_iter, from = max_iter)
			}
			stop(e)
		}),
		warning = function(w) {
			if (conditionMessage(w) %in% seen) invokeRestart("muffleWarning")
		},
		message = function(m) {
			if (conditionMessage(m) %in% seen) invokeRestart("muffleMessage")
		}
	)
	if (isTRUE(verbose)) {
		cli::cli_alert_success("Converged with the iteration limit raised to {budget}.")
	}
	.ae_withhold_product(out, lower)
}
