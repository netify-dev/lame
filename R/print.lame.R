#' Print method for LAME objects
#'
#' Provides a concise print output for fitted LAME models. When the fit has
#' two or more dynamic effects active simultaneously
#' (\code{dynamic_uv + dynamic_ab + dynamic_beta} in any
#' combination), a
#' compact joint table summarising the AR(1) hyperparameters of each
#' dynamic block is printed instead of a separate paragraph per component.
#' Disable the compact mode by passing \code{compact = FALSE}.
#' A fit that stores no posterior draws (\code{max_seconds} elapsed before the
#' first stored draw) prints its dimensions and says so.
#'
#' @param x an object of class "lame"
#' @param compact logical: when \code{TRUE} (default), use a compact joint
#'   table for fits with 2+ dynamic components; set to \code{FALSE} to
#'   force the per-component long-form display.
#' @param digits Number of digits to display in the compact table. Default 3.
#' @param ... additional arguments (not used)
#' @return the lame object invisibly
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @method print lame
#' @export
#' @import cli
print.lame <- function(x, compact = TRUE, digits = 3, ...) {
	.lame_cli_stdout(.lame_print_lame_body(x, compact = compact, digits = digits))
	invisible(x)
}

# the cli output of print.lame; the method writes it to standard output
.lame_print_lame_body <- function(x, compact = TRUE, digits = 3) {
	cli::cli_h1("Longitudinal additive and multiplicative effects (LAME) model")

	# number of periods, from n_time or else the stored responses
	n_periods <- x$n_time
	if (is.null(n_periods)) {
		n_periods <- if (length(x$Y_T) > 0) length(x$Y_T)
		else if (is.array(x$Y) && length(dim(x$Y)) == 3L) dim(x$Y)[3]
		else if (is.list(x$Y)) length(x$Y)
		else NA_integer_
	}
	# row and column counts, so a bipartite fit prints nA x nB
	n_row <- NULL; n_col <- NULL
	if (length(x$Y_T) > 0) {
		n_row <- nrow(x$Y_T[[1]]); n_col <- ncol(x$Y_T[[1]])
	} else if (is.array(x$Y) && length(dim(x$Y)) == 3L) {
		n_row <- dim(x$Y)[1]; n_col <- dim(x$Y)[2]
	} else if (is.list(x$Y) && length(x$Y) > 0) {
		n_row <- nrow(x$Y[[1]]); n_col <- ncol(x$Y[[1]])
	}
	if (!is.null(n_row)) {
		cli::cli_h3("Network structure")
		mode_str <- if (identical(x$mode, "bipartite")) " (bipartite)" else ""
		cli::cli_bullets(c(
			"*" = "Time periods: {.val {n_periods}}",
			"*" = "Network dimensions: {.val {n_row}} x {.val {n_col}} per period{mode_str}"
		))
	}

	# the first dimension of BETA indexes draws in both 2-d and 3-d storage
	n_stored <- dim(x$BETA)[1]
	cli::cli_bullets(c(
		"*" = "Stored posterior samples: {.val {n_stored}}"
	))
	if (.lame_has_no_draws(x)) {
		cli::cli_alert_warning("No posterior draw was stored: {.arg max_seconds} elapsed before the first stored draw, so the fit has no estimates.")
		cli::cli_alert_info("Continue the chains with {.code lame_resume(fit, nscan_more = N, max_seconds = Inf)}, or refit with a larger {.arg max_seconds} or a shorter {.arg burn}.")
		return(invisible(x))
	}

	if (!is.null(x$model.name)) {
		cli::cli_bullets(c(
			"*" = "Model type: {.field {x$model.name}}"
		))
	}

	# dynamic flags come from the fit, since a static bipartite fit can also store U as a 3-d array
	is_dynamic_uv <- isTRUE(x$dynamic_uv)
	is_dynamic_ab <- isTRUE(x$dynamic_ab) || !is.null(x$a_dynamic)
	is_dynamic_beta <- isTRUE(x$dynamic_beta) || length(dim(x$BETA)) == 3L

	cli::cli_h3("Model parameters")
	param_info <- character()
	n_coef <- if (length(dim(x$BETA)) == 3L) dim(x$BETA)[2] else ncol(x$BETA)
	param_info <- c(param_info, "*" = "Regression coefficients: {.val {n_coef}}")
	if (is_dynamic_beta) {
		n_dyn <- if (!is.null(x$beta_dynamic_mask)) sum(x$beta_dynamic_mask) else NA_integer_
		grp_str <- if (!is.null(x$beta_dynamic_groups))
			paste(unique(x$beta_dynamic_groups[x$beta_dynamic_groups != ""]), collapse = ", ")
		else ""
		param_info <- c(
			param_info,
			"*" = "Dynamic regression coefficients: {.val {n_dyn}}{if (nzchar(grp_str)) paste0(' (blocks: ', grp_str, ')') else ''}"
		)
	}

	if (!is.null(x$U)) {
		if (is_dynamic_uv) {
			R <- dim(x$U)[2]
			Tn <- dim(x$U)[3]
			param_info <- c(
				param_info,
				"*" = "Multiplicative effects: {.val {R}}-dimensional {.emph (dynamic over {Tn} periods)}"
			)
		} else {
			param_info <- c(
				param_info,
				"*" = "Multiplicative effects: {.val {ncol(x$U)}}-dimensional"
			)
		}
	}
	if (is_dynamic_ab) {
		param_info <- c(
			param_info,
			"*" = "Additive effects: {.emph dynamic} (per-actor level plus AR(1) deviation)"
		)
	}

	cli::cli_bullets(param_info)

	if (!is.null(.bound_line <- .lame_bound_share_line(x))) cli::cli_alert_info(.bound_line)
	if (!is.null(.chain_line <- .lame_chain_print_line(x$auto_extend, x))) cli::cli_alert_info(.chain_line)

	# what the fit does not report: the effect variances and intercept of a
	# sparse panel, the intercept of a one-period cap
	for (st_ in .lame_ab_statements(x)) cli::cli_alert_info("{st_}")

	# count of unique nodes when actors enter or exit
	if (!is.null(x$nodeID)) {
		total_nodes <- length(unique(unlist(x$nodeID)))
		cli::cli_bullets(c(
			"i" = "Total unique nodes across periods: {.val {total_nodes}}"
		))
	}

	# warn when a panel fit keeps U, a and b time-invariant
	if (!is.null(n_periods) && !is.na(n_periods) && n_periods > 1L &&
			!is_dynamic_uv && !is_dynamic_ab && !is_dynamic_beta) {
		cli::cli_alert_warning(c(
			"This is a {.strong STATIC} fit pooled across {n_periods} time periods: ",
			"{.code U}, {.code V}, {.code a}, {.code b} do not vary by time."
		))
		cli::cli_alert_info(c(
			"For time-varying latent positions / additive effects, refit with ",
			"{.code dynamic_uv = TRUE} and/or {.code dynamic_ab = TRUE}."
		))
	}

	# at rank two or more the snap and t samplers do not recover the position
	# hyperparameters (see .lame_uv_hyper_unidentified)
	uv_unident <- .lame_uv_hyper_unidentified(x)
	# the snap kind's reset rate where it is reported
	pi_snap_shown <- identical(x$dynamic_uv_kind, "snap") && is.numeric(x$pi_snap) &&
		length(x$pi_snap) > 0L && !.lame_snap_rate_unidentified(x)
	# compact table when two or more dynamic components are active
	if (isTRUE(compact)) {
		n_dyn <- sum(c(is_dynamic_uv, is_dynamic_ab, is_dynamic_beta))
		.qfmt <- function(v) {
			if (is.null(v) || length(v) == 0L) return("")
			q <- stats::quantile(as.numeric(v), c(0.025, 0.5, 0.975), na.rm = TRUE)
			sprintf("%.*f [%.*f, %.*f]", digits, q[2L], digits, q[1L], digits, q[3L])
		}
		if (n_dyn >= 2L) {
			cli::cli_h3("Dynamic components")
			rows <- list()
			if (is_dynamic_uv) {
				rows[[length(rows) + 1L]] <- data.frame(
					component = "U/V",
					dim = if (!is.null(x$U)) dim(x$U)[2] else NA_integer_,
					rho = if (uv_unident) "not reported" else .qfmt(x$rho_uv),
					sigma_innov = if (uv_unident) "not reported" else
					if (isTRUE(x$symmetric)) "sqrt(1 - rho^2) (gauge)" else .qfmt(x$sigma_uv),
					stringsAsFactors = FALSE
				)
			}
			if (is_dynamic_ab) {
				rows[[length(rows) + 1L]] <- data.frame(
					component = "a/b",
					dim = if (!is.null(x$APM)) length(x$APM) else NA_integer_,
					rho = .qfmt(x$rho_ab),
					sigma_innov = .qfmt(x$sigma_ab),
					stringsAsFactors = FALSE
				)
			}
			if (is_dynamic_beta) {
				n_dyn_beta <- if (!is.null(x$beta_dynamic_mask))
					sum(x$beta_dynamic_mask) else NA_integer_
				rb <- if (!is.null(x$RHO_BETA)) {
					# per-draw mean across blocks; summary() gives each block
					apply(x$RHO_BETA, 1, mean)
				} else x$rho_beta
				sb <- if (!is.null(x$SIGMA_BETA)) {
					apply(x$SIGMA_BETA, 1, mean)
				} else x$sigma_beta
				rows[[length(rows) + 1L]] <- data.frame(
					component = "beta",
					dim = n_dyn_beta,
					rho = .qfmt(rb),
					sigma_innov = .qfmt(sb),
					stringsAsFactors = FALSE
				)
			}
			if (length(rows) > 0L) {
				tbl <- do.call(rbind, rows)
				cli::cli_verbatim(utils::capture.output(print(tbl, row.names = FALSE)))
				cli::cli_alert_info("{.field sigma_innov} columns are AR(1) innovation SDs (not marginal/stationary SDs).")
				if (uv_unident) {
					cli::cli_alert_info("{.code rho_uv} / {.code sigma_uv} are not reported for {.code dynamic_uv_kind = {.val {x$dynamic_uv_kind %||% \"snap\"}}} at rank two or more (draws in {.code fit$rho_uv}, {.code fit$sigma_uv}).")
				}
				if (.lame_snap_rate_unidentified(x)) {
					cli::cli_alert_info("The reset rate {.code pi_snap} is not reported at rank two or more (draws in {.code fit$pi_snap}).")
				}
				if (pi_snap_shown) {
					pi_txt <- .qfmt(x$pi_snap)
					cli::cli_alert_info("Reset rate {.code pi_snap}: {pi_txt}.")
				}
			}
		} else if (n_dyn == 1L) {
			# one dynamic component: its AR(1) hyperparameters, one line each
			hp <- character()
			if (is_dynamic_uv) {
				hp <- c(
					hp,
					"*" = if (uv_unident) {
						paste0(
							"Latent positions AR(1): rho_uv and sigma_uv not reported for this transition kind at rank two or more",
							if (.lame_snap_rate_unidentified(x)) ", nor the reset rate pi_snap" else "",
							" (the chains do not recover them)"
						)
					} else paste0(
						"Latent positions AR(1): rho_uv ",
						.qfmt(x$rho_uv),
						if (isTRUE(x$symmetric)) ", innovation sd sqrt(1 - rho_uv^2) (gauge)" else paste0(
							", innovation sd ",
							.qfmt(x$sigma_uv)
						),
						if (pi_snap_shown) paste0(", reset rate pi_snap ", .qfmt(x$pi_snap)) else ""
					)
				)
			}
			if (is_dynamic_ab) {
				hp <- c(
					hp,
					"*" = paste0(
						"Additive effects AR(1): rho_ab ",
						.qfmt(x$rho_ab),
						", innovation sd ",
						.qfmt(x$sigma_ab)
					)
				)
			}
			if (is_dynamic_beta) {
				rb <- if (!is.null(x$RHO_BETA)) apply(x$RHO_BETA, 1, mean) else x$rho_beta
				sb <- if (!is.null(x$SIGMA_BETA)) apply(x$SIGMA_BETA, 1, mean) else x$sigma_beta
				hp <- c(hp, "*" = paste0("Coefficient AR(1): rho_beta ", .qfmt(rb), ", innovation sd ", .qfmt(sb)))
			}
			if (length(hp) > 0L) {
				cli::cli_h3("AR(1) hyperparameters")
				cli::cli_bullets(hp)
			}
		}
	}

	cli::cli_alert_info("Use {.code summary(object)} for detailed results")

	invisible(x)
}

# write the cli output of a display method to standard output, keeping its styling,
# so capture.output(), sink() and knitr see it
.lame_cli_stdout <- function(expr) {
	lines <- cli::cli_fmt(expr)
	if (length(lines)) cli::cat_line(lines)
	invisible(NULL)
}
