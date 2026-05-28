#' Print method for LAME objects
#'
#' Provides a concise print output for fitted LAME models. When the fit has
#' two or more dynamic effects active simultaneously
#' (\code{dynamic_uv + dynamic_ab + dynamic_beta} in any combination), a
#' compact joint table summarising the AR(1) hyperparameters of each
#' dynamic block is printed instead of a separate paragraph per component.
#' Disable the compact mode by passing \code{compact = FALSE}.
#'
#' @param x an object of class "lame"
#' @param compact logical: when \code{TRUE} (default), use a compact joint
#'   table for fits with 2+ dynamic components; set to \code{FALSE} to
#'   force the per-component long-form display.
#' @param digits Number of digits to display in the compact table. Default 3.
#' @param ... additional arguments (not used)
#' @return the lame object invisibly
#' @author Peter Hoff, Cassy Dorff, Shahryar Minhas
#' @method print lame
#' @export
#' @import cli
print.lame <- function(x, compact = TRUE, digits = 3, ...) {
	cli::cli_h1("Longitudinal Additive and Multiplicative Effects (LAME) Model")
	
	# network dimensions for display. fall back to Y_T when n_time and Y
	# are not populated.
	n_periods <- x$n_time
	if (is.null(n_periods)) {
		n_periods <- if (length(x$Y_T) > 0) length(x$Y_T)
		             else if (is.array(x$Y) && length(dim(x$Y)) == 3L) dim(x$Y)[3]
		             else if (is.list(x$Y)) length(x$Y)
		             else NA_integer_
	}
	# capture both row and column counts so a bipartite (rectangular) fit
	# prints "nA x nB".
	n_row <- NULL; n_col <- NULL
	if (length(x$Y_T) > 0) {
		n_row <- nrow(x$Y_T[[1]]); n_col <- ncol(x$Y_T[[1]])
	} else if (is.array(x$Y) && length(dim(x$Y)) == 3L) {
		n_row <- dim(x$Y)[1]; n_col <- dim(x$Y)[2]
	} else if (is.list(x$Y) && length(x$Y) > 0) {
		n_row <- nrow(x$Y[[1]]); n_col <- ncol(x$Y[[1]])
	}
	if (!is.null(n_row)) {
		cli::cli_h3("Network Structure")
		mode_str <- if (identical(x$mode, "bipartite")) " (bipartite)" else ""
		cli::cli_bullets(c(
			"*" = "Time periods: {.val {n_periods}}",
			"*" = "Network dimensions: {.val {n_row}} x {.val {n_col}} per period{mode_str}"
		))
	}
	
	# first dim of BETA is iter for both 2-D and 3-D (dynamic_beta) storage
	n_stored <- dim(x$BETA)[1]
	cli::cli_bullets(c(
		"*" = "Stored posterior samples: {.val {n_stored}}"
	))
	
	if (!is.null(x$model.name)) {
		cli::cli_bullets(c(
			"*" = "Model type: {.field {x$model.name}}"
		))
	}
	
	# detect dynamic UV from the FIT FLAG, not from U's array dimensions: the
	# bipartite static path stores U as [nA, R, T] too, so the shape-only check
	# misclassified static bipartite fits as dynamic and suppressed the STATIC
	# notice below. (Unipartite static stores U as 2-D, so the old check happened
	# to work there.) Prefer the explicit fit$dynamic_uv flag.
	is_dynamic_uv <- isTRUE(x$dynamic_uv)
	is_dynamic_ab <- isTRUE(x$dynamic_ab) || !is.null(x$a_dynamic)
	is_dynamic_beta <- isTRUE(x$dynamic_beta) || length(dim(x$BETA)) == 3L

	# summarize what the model is estimating
	cli::cli_h3("Model Parameters")
	param_info <- character()
	n_coef <- if (length(dim(x$BETA)) == 3L) dim(x$BETA)[2] else ncol(x$BETA)
	param_info <- c(param_info, "*" = "Regression coefficients: {.val {n_coef}}")
	if (is_dynamic_beta) {
		n_dyn <- if (!is.null(x$beta_dynamic_mask)) sum(x$beta_dynamic_mask) else NA_integer_
		grp_str <- if (!is.null(x$beta_dynamic_groups))
			paste(unique(x$beta_dynamic_groups[x$beta_dynamic_groups != ""]), collapse = ", ")
			else ""
		param_info <- c(param_info,
			"*" = "Dynamic regression coefficients: {.val {n_dyn}}{if (nzchar(grp_str)) paste0(' (blocks: ', grp_str, ')') else ''}")
	}

	if (!is.null(x$U)) {
		if (is_dynamic_uv) {
			R <- dim(x$U)[2]
			T <- dim(x$U)[3]
			param_info <- c(param_info,
				"*" = "Multiplicative effects: {.val {R}}-dimensional {.emph (dynamic over {T} periods)}")
		} else {
			param_info <- c(param_info,
				"*" = "Multiplicative effects: {.val {ncol(x$U)}}-dimensional")
		}
	}
	cli::cli_bullets(param_info)

	# helpful when nodes enter/exit across periods
	if (!is.null(x$nodeID)) {
		total_nodes <- length(unique(unlist(x$nodeID)))
		cli::cli_bullets(c(
			"i" = "Total unique nodes across periods: {.val {total_nodes}}"
		))
	}

	# loud notice for the most common novice surprise: a panel + default flags
	# fits a STATIC model. U / a / b are time-invariant; per-period EZ varies
	# only through per-period X covariates.
	if (!is.null(n_periods) && !is.na(n_periods) && n_periods > 1L &&
	    !is_dynamic_uv && !is_dynamic_ab && !is_dynamic_beta) {
		cli::cli_alert_warning(c(
			"This is a {.strong STATIC} fit pooled across {n_periods} time periods: ",
			"{.code U}, {.code V}, {.code a}, {.code b} do not vary by time."))
		cli::cli_alert_info(c(
			"For time-varying latent positions / additive effects, refit with ",
			"{.code dynamic_uv = TRUE} and/or {.code dynamic_ab = TRUE}."))
	}

	# compact joint table when 2+ dynamic components are active.
	# Replaces the per-component paragraphs that would otherwise become
	# repetitive in joint dynamic_uv + dynamic_ab + dynamic_beta fits.
	if (isTRUE(compact)) {
		n_dyn <- sum(c(is_dynamic_uv, is_dynamic_ab, is_dynamic_beta))
		if (n_dyn >= 2L) {
			cli::cli_h3("Dynamic components")
			.qfmt <- function(v) {
				if (is.null(v) || length(v) == 0L) return("")
				q <- stats::quantile(as.numeric(v), c(0.025, 0.5, 0.975), na.rm = TRUE)
				sprintf("%.*f [%.*f, %.*f]", digits, q[2L], digits, q[1L], digits, q[3L])
			}
			rows <- list()
			if (is_dynamic_uv) {
				rows[[length(rows) + 1L]] <- data.frame(
					component   = "U/V",
					dim         = if (!is.null(x$U)) dim(x$U)[2] else NA_integer_,
					rho         = .qfmt(x$rho_uv),
					sigma_innov = .qfmt(x$sigma_uv),
					stringsAsFactors = FALSE
				)
			}
			if (is_dynamic_ab) {
				rows[[length(rows) + 1L]] <- data.frame(
					component   = "a/b",
					dim         = if (!is.null(x$APM)) length(x$APM) else NA_integer_,
					rho         = .qfmt(x$rho_ab),
					sigma_innov = .qfmt(x$sigma_ab),
					stringsAsFactors = FALSE
				)
			}
			if (is_dynamic_beta) {
				n_dyn_beta <- if (!is.null(x$beta_dynamic_mask))
					sum(x$beta_dynamic_mask) else NA_integer_
				rb <- if (!is.null(x$RHO_BETA)) {
					# average of per-iteration medians across blocks for a single
					# compact row; user gets per-block detail via summary()
					apply(x$RHO_BETA, 1, mean)
				} else x$rho_beta
				sb <- if (!is.null(x$SIGMA_BETA)) {
					apply(x$SIGMA_BETA, 1, mean)
				} else x$sigma_beta
				rows[[length(rows) + 1L]] <- data.frame(
					component   = "beta",
					dim         = n_dyn_beta,
					rho         = .qfmt(rb),
					sigma_innov = .qfmt(sb),
					stringsAsFactors = FALSE
				)
			}
			if (length(rows) > 0L) {
				tbl <- do.call(rbind, rows)
				print(tbl, row.names = FALSE)
				cli::cli_alert_info("{.field sigma_innov} columns are AR(1) innovation SDs (not marginal/stationary SDs).")
			}
		}
	}

	cli::cli_alert_info("Use {.code summary(object)} for detailed results")

	invisible(x)
}