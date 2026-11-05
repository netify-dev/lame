# retirement of factor columns the data do not support (empirical-bayes route)
# a column the data barely inform leaves no settled iteration, so the weakest column is dropped
# and the fit refit at the next rank down; retired columns are in fit$eb$retirement

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

# run the empirical-bayes iteration at rank R and retire columns until it settles
# an iteration that only ran out of a budget below the default is refused, not retired
.ae_retire_columns <- function(run, init, R, symmetric, vb_follows, max_iter = .ae_retire_min_iter) {
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
		steps[[length(steps) + 1L]] <- list(
			rank = Rr,
			column = weakest,
			strength = strength[weakest],
			share = strength[weakest] / max(strength),
			strengths = strength,
			stage = stop_now$stage,
			termination = stop_now$termination
		)
		keep <- setdiff(seq_len(Rr), weakest)
		keep_init <- .ae_init_keep(keep_init, keep)
		Rr <- Rr - 1L
		ir <- run(keep_init, Rr)
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

# message for a fit with retired columns (print and summary repeat it)
.ae_inform_retirement <- function(ret) {
	if (is.null(ret)) return(invisible(NULL))
	req <- ret$requested_rank
	eff <- ret$effective_rank
	cli::cli_inform(
		c(
			"!" = "Rank {req} was requested; the data support {eff}. {req - eff} factor column{?s} not supported by the data {?was/were} retired, and the fit is the rank-{eff} fit{if (eff == 0L) ' (no latent factors)' else ''}.",
			" " = "{.code fit$R} is the rank of the fit and {.code fit$R_requested} the rank asked for; the record is in {.code fit$eb$retirement}."
		),
		class = "lame_als_retirement"
	)
	invisible(NULL)
}

# the lines the print and summary methods show for a retirement
.ae_print_retirement <- function(ret, requested = NULL) {
	if (is.null(ret)) return(invisible(NULL))
	req <- ret$requested_rank
	eff <- ret$effective_rank
	cli::cli_alert_warning(c(
		"Rank {req} was requested; the data support {eff}. ",
		"{req - eff} factor column{?s} not supported by the data {?was/were} retired (the fit is the rank-{eff} fit)."
	))
	invisible(NULL)
}
