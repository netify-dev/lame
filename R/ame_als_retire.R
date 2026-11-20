# retirement of factor columns the data do not support (empirical-bayes route)
# a column the data barely inform leaves no settled iteration, so the weakest column is dropped
# and the fit refit at the next rank down, provided the lower rank's variational lower bound is
# at least as high; retired columns are in fit$eb$retirement. a rank whose bound is higher is
# kept: its last iterate is the fit, with a warning (fit$eb$unsettled), unless some reported
# quantity cannot be formed there, when the fit of the next rank down is returned

# the number of factor columns of a fit (0 when it has none)
.ae_fit_rank <- function(fit) {
	if (is.null(fit$U)) return(0L)
	ncol(as.matrix(fit$U))
}

# whether the empirical-bayes iteration `ir` meets its stationarity certificate and,
# when a variational refinement follows, the refinement has converged
.ae_ir_settled <- function(ir, vb_follows) {
	isTRUE(ir$stationarity$pass) && (!vb_follows || isTRUE(ir$fit$eb_mf$record$converged))
}

# strength of every factor column: ||u_k|| ||v_k||, or |l_k| for the undirected eigenmodel
.ae_column_strength <- function(fit, symmetric) {
	if (.ae_fit_rank(fit) == 0L) return(numeric(0))
	if (isTRUE(symmetric) && !is.null(fit$L)) return(abs(as.numeric(fit$L)))
	sqrt(colSums(as.matrix(fit$U)^2) * colSums(as.matrix(fit[["V"]])^2))
}

# start of a lower-rank run: the columns of a start that stay
.ae_init_keep <- function(init, keep) {
	if (is.null(init)) return(NULL)
	# no column left means the cold start, whose rank-0 fit has no factor matrices
	if (!length(keep)) return(NULL)
	for (nm in c("U", "V")) {
		if (!is.null(init[[nm]])) init[[nm]] <- as.matrix(init[[nm]])[, keep, drop = FALSE]
	}
	if (!is.null(init$L)) init$L <- init$L[keep]
	init
}

# what stopped the iteration, in words, for the record of a retirement
.ae_ir_stop <- function(ir, vb_follows) {
	if (!isTRUE(ir$stationarity$pass)) {
		return(list(stage = "posterior mode", termination = ir$stationarity$termination %||% "not converged"))
	}
	rec <- ir$fit$eb_mf$record
	list(stage = "variational refinement", termination = rec$termination %||% "not converged")
}

# default iteration budget of ame_als() and lame_als() (their `max_iter`)
.ae_retire_min_iter <- 200L

# variational lower bound of an iteration, NA when it has none (its variational stage did not
# produce one)
.ae_ir_elbo <- function(ir) {
	elbo <- ir$fit$eb_mf$record$elbo
	if (is.null(elbo)) NA_real_ else elbo
}

# run the empirical-bayes iteration at rank R and retire columns until it settles. an unsettled
# rank keeps its columns when its lower bound is finite and above that of the next rank down, or
# the next rank down did not settle and has no bound to compare with (a rank is never retired
# against a missing bound), and `formed` can form everything the fit reports at its last iterate: that iteration is
# returned with the record `unsettled`. when `formed` cannot (it returns NULL, or the reason as a
# string, such as "bootstrap" for a bootstrap the call asks for), the rank is
# retired as one that did not converge and the fit of the next rank down, already computed,
# continues the loop. an iteration that only ran out of a budget below the default is refused,
# not retired
.ae_retire_columns <- function(run, init, R, symmetric, vb_follows, max_iter = .ae_retire_min_iter, formed = identity) {
	ir <- run(init, R)
	if (!vb_follows || R == 0L) return(ir)
	budget_ok <- isTRUE(max_iter >= .ae_retire_min_iter)
	steps <- list()
	Rr <- R
	keep_init <- init
	while (!.ae_ir_settled(ir, vb_follows) && Rr > 0L) {
		strength <- .ae_column_strength(ir$fit, symmetric)
		if (length(strength) != Rr || !all(is.finite(strength))) break
		weakest <- which.min(strength)
		stop_now <- .ae_ir_stop(ir, vb_follows)
		if (!budget_ok && identical(stop_now$termination, "iteration limit")) break
		keep <- setdiff(seq_len(Rr), weakest)
		lower_init <- .ae_init_keep(keep_init, keep)
		lower <- run(lower_init, Rr - 1L)
		elbo_r <- .ae_ir_elbo(ir); elbo_lower <- .ae_ir_elbo(lower)
		reason <- "lower bound"
		# a lower rank without a bound is preferred only if it settled
		above <- if (is.finite(elbo_lower)) elbo_r > elbo_lower else !.ae_ir_settled(lower, vb_follows)
		if (is.finite(elbo_r) && above) {
			done <- formed(ir)
			if (is.list(done)) {
				ir <- done
				ir$unsettled <- c(
					list(rank = Rr, iterations = .ae_ir_stop_iterations(ir, stop_now$stage)),
					stop_now,
					list(elbo = elbo_r, elbo_lower = elbo_lower)
				)
				break
			}
			reason <- if (is.character(done)) done else "not formed"
		}
		steps[[length(steps) + 1L]] <- list(
			rank = Rr,
			column = weakest,
			strength = strength[weakest],
			share = strength[weakest] / max(strength),
			strengths = strength,
			stage = stop_now$stage,
			termination = stop_now$termination,
			iterations = .ae_ir_stop_iterations(ir, stop_now$stage),
			reason = reason
		)
		keep_init <- lower_init
		Rr <- Rr - 1L
		ir <- lower
	}
	if (length(steps)) {
		ir$retirement <- list(
			requested_rank = R,
			effective_rank = Rr,
			steps = steps,
			settled = .ae_ir_settled(ir, vb_follows)
		)
	}
	ir
}

# iterations of the stage that stopped (reweighting steps of the mode, or variational iterations)
.ae_ir_stop_iterations <- function(ir, stage) {
	n <- if (identical(stage, "posterior mode")) ir$irls_iters else ir$fit$eb_mf$record$iterations
	if (is.null(n)) NA_integer_ else as.integer(n)
}

# what stopped a stage, in words
.ae_stop_words <- function(termination, iterations) {
	it <- if (isTRUE(is.finite(iterations))) paste0(" after ", iterations, " iteration", if (iterations == 1L) "" else "s") else ""
	switch(
		termination,
		"iteration limit" = paste0("reached its iteration limit", it),
		"no improving step" = ,
		"no improving conditional step" = paste0("found no step that improves its objective", it),
		paste0("stopped (", termination, ")", it)
	)
}

# whether every quantity a fit reports is formed and finite: coefficients, standard errors,
# covariance, intervals, positions, fitted values, predictions and a simulated network
# (missing cells, such as a diagonal or an actor without ties, stay missing)
.ae_fit_formed <- function(out) {
	if (!is.list(out) || !inherits(out, "ame_als")) return(FALSE)
	num <- function(x) as.numeric(unlist(x, use.names = FALSE))
	all_finite <- function(x) { v <- num(x); length(v) > 0L && all(is.finite(v)) }
	no_bad <- function(x) { v <- num(x); !any(is.nan(v) | is.infinite(v)) }
	eb <- out$eb
	if (is.null(eb$vcov) || !is.null(eb$vcov_singular) || !all_finite(eb$vcov) || !all_finite(eb$se)) return(FALSE)
	ez <- num(out$EZ)
	if (!any(is.finite(ez)) || !no_bad(ez)) return(FALSE)
	# an observed cell with a finite predictor has a finite value (a list of periods is checked as a whole)
	y <- num(out$Y)
	cell <- if (length(y) == length(ez)) is.finite(ez) & !is.na(y) else is.finite(ez)
	cells_ok <- function(x) {
		v <- num(x)
		if (length(v) != length(ez)) return(any(is.finite(v)) && no_bad(v))
		no_bad(v) && all(is.finite(v[cell]))
	}
	for (nm in c("U", "V", "L", "UVPM", "ULUPM", "a", "b")) if (!no_bad(out[[nm]])) return(FALSE)
	# a method that refuses by design (class lame_als_prediction) refuses for a settled fit too
	can <- function(expr) tryCatch(isTRUE(expr), lame_als_prediction = function(e) TRUE, error = function(e) FALSE)
	ok <- can(all_finite(stats::coef(out))) && can(all_finite(stats::vcov(out))) &&
		can(all_finite(stats::confint(out))) && cells_ok(out$YPM) && can(cells_ok(stats::fitted(out))) &&
		can(cells_ok(stats::predict(out))) && can(cells_ok(stats::simulate(out, nsim = 1L, seed = 1L)$Y))
	ok
}

# the sentences that describe a retirement: the rule of the lower bound, or a rank whose
# estimates could not be formed after it did not converge
.ae_retirement_text <- function(ret) {
	req <- ret$requested_rank
	eff <- ret$effective_rank
	none <- if (eff == 0L) " (no latent factors)" else ""
	unformed <- Filter(function(st) isTRUE(st$reason %in% c("not formed", "bootstrap")), ret$steps)
	if (!length(unformed)) {
		return(list(
			head = paste0(
				"Rank ", req, " was requested; the data support ", eff, ". ",
				cli::pluralize("{req - eff} factor column{?s} not supported by the data {?was/were} retired"),
				", and the fit is the rank-", eff, " fit", none, "."
			),
			why = character(0),
			more = FALSE
		))
	}
	why <- vapply(unformed, function(st) {
		paste0(
			"Rank ", st$rank, " did not converge: its ", st$stage, " ",
			.ae_stop_words(st$termination, st$iterations),
			if (identical(st$reason, "bootstrap")) ", and the bootstrap asked for is not formed at an iterate that has not settled." else
				", and its estimates could not all be formed at the last iterate."
		)
	}, character(1))
	list(
		head = paste0("Rank ", req, " was requested; the fit is the rank-", eff, " fit", none, "."),
		why = why,
		more = any(vapply(unformed, function(st) identical(st$termination, "iteration limit"), logical(1)))
	)
}

# message for a fit with retired columns (print and summary repeat it)
.ae_inform_retirement <- function(ret) {
	if (is.null(ret)) return(invisible(NULL))
	tx <- .ae_retirement_text(ret)
	cli::cli_inform(
		c(
			"!" = tx$head,
			stats::setNames(tx$why, rep("x", length(tx$why))),
			if (tx$more) c("i" = "A larger {.arg max_iter} continues the fit at the rank asked for."),
			" " = "{.code fit$R} is the rank of the fit and {.code fit$R_requested} the rank asked for; the record is in {.code fit$eb$retirement}."
		),
		class = "lame_als_retirement"
	)
	invisible(NULL)
}

# the lines the print and summary methods show for a retirement
.ae_print_retirement <- function(ret) {
	if (is.null(ret)) return(invisible(NULL))
	tx <- .ae_retirement_text(ret)
	cli::cli_alert_warning(paste(c(tx$head, tx$why), collapse = " "))
	invisible(NULL)
}

# the sentence for a fit returned at a rank whose iteration did not settle (fit$eb$unsettled)
.ae_unsettled_text <- function(un) {
	paste0(
		"The rank-", un$rank, " fit did not converge: its ", un$stage, " ",
		.ae_stop_words(un$termination, un$iterations),
		". Its estimates are those at the last iterate."
	)
}

# warning for such a fit, given once when it is made
.ae_warn_unsettled <- function(un, log_marginal_missing) {
	cli::cli_warn(
		c(
			.ae_unsettled_text(un),
			"i" = paste0(
				if (isTRUE(is.finite(un$elbo_lower))) {
					paste0("Rank ", un$rank, " has a higher variational lower bound than rank ", un$rank - 1L)
				} else {
					paste0("Rank ", un$rank - 1L, " has no variational lower bound to compare with")
				},
				", so the rank-", un$rank, " fit is returned",
				if (log_marginal_missing) "; the Laplace marginal likelihood {.code fit$eb$log_marginal}, defined at a settled mode, is NA." else "."
			),
			"i" = paste0(
				if (identical(un$termination, "iteration limit")) "A larger {.arg max_iter} continues the fit; the record is in " else "The record is in ",
				"{.code fit$eb$unsettled}."
			)
		),
		class = "lame_als_unsettled"
	)
	invisible(NULL)
}

# the line the print and summary methods show for it
.ae_print_unsettled <- function(un) {
	if (is.null(un)) return(invisible(NULL))
	more <- if (identical(un$termination, "iteration limit")) "A larger max_iter continues the fit."
	cli::cli_alert_warning(paste(c(.ae_unsettled_text(un), more), collapse = " "))
	invisible(NULL)
}

# an unsettled fit stores as missing what is defined only at a settled mode: the Laplace
# evidence at the mode (the rank criterion), which a normal fit also reads at the settled
# noise variance of its variational stage
.ae_mark_unsettled <- function(out, un, family) {
	out$converged <- FALSE
	out$eb$unsettled <- un
	at_mode <- identical(un$stage, "posterior mode") || identical(family, "normal")
	if (at_mode) {
		out$eb$log_marginal <- NA_real_
		if (is.list(out$eb$laplace)) out$eb$laplace[] <- lapply(out$eb$laplace, function(x) if (is.numeric(x)) NA_real_ else x)
		if (!is.null(out$eb$variational$mode_log_marginal)) out$eb$variational$mode_log_marginal <- NA_real_
	}
	out
}
