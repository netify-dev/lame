# family-specific predictive scores on held-out cells

#' Held-out predictive evaluation for an ame / lame fit
#'
#' Computes family-appropriate held-out predictive scores given a fit
#' and a logical mask of cells to score. The function does not split the
#' data for you: it expects you to have either (a) refit on a training
#' subset and now want to score the held-out cells of the same matrix,
#' or (b) have predicted probabilities you want scored. Both AUROC and
#' PR-AUC are reported when applicable, alongside Brier and mean log
#' density.
#'
#' \strong{Workflow.} The typical pattern is: mask a random sample of
#' dyads to \code{NA} in \code{Y}, refit (\code{lame()} handles \code{NA}
#' internally via data augmentation), call
#' \code{predict(fit, type = "response")}, then pass that prediction
#' alongside the original \code{Y} and the held-out mask to this
#' function. See the examples.
#'
#' \strong{Dependencies.} AUROC / PR-AUC use \pkg{precrec} when
#' available; if not installed, only Brier and mean log-density are
#' computed and a one-line note is emitted.
#'
#' @param y_obs Observed outcomes. For longitudinal fits, a list of
#'   per-period matrices; for cross-sectional, a single matrix. Cells
#'   not in \code{mask} are ignored.
#' @param y_pred Predicted probabilities / means on the response scale.
#'   Same shape as \code{y_obs}.
#' @param mask Logical mask of the same shape as \code{y_obs} marking
#'   cells to score (\code{TRUE} = include in evaluation).
#' @param family Family string; used to pick the scoring rule. Defaults
#'   to \code{"binary"}.
#' @param fit For \code{family = "normal"}: the fit that produced
#'   \code{y_pred} (an \code{ame}/\code{lame} MCMC fit or an
#'   \code{ame_als}/\code{lame_als} fit of the normal family). Its residual
#'   variance sets the predictive density of \code{mean_logdens}; without it
#'   \code{mean_logdens} is \code{NA}. Ignored for the other families.
#'
#' @return A one-row data frame with columns appropriate to the family:
#'   \code{n_eval}, plus \code{auroc} + \code{auprc} + \code{brier} +
#'   \code{logloss} (binary / cbin), or \code{rmse} + \code{mae} +
#'   \code{mean_logdens} (normal), or \code{mean_logdens} +
#'   \code{rmse} (poisson).
#'
#'   For the normal family \code{mean_logdens} is the mean over the scored
#'   cells of \eqn{\log \phi(y; \hat y, v_e)}, the normal log density at the
#'   prediction with the fit's residual variance \eqn{v_e}: the posterior mean
#'   of \code{VC[, "ve"]} for an MCMC fit, the fitted \code{VC["ve"]} for a
#'   fast fit. It is \code{NA} when \code{fit} is not supplied.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 25; Y <- matrix(rbinom(n*n, 1, 0.3), n, n); diag(Y) <- NA
#' rownames(Y) <- colnames(Y) <- paste0("a", sprintf("%02d", 1:n))
#' # mask 20% of dyads
#' mask <- matrix(FALSE, n, n)
#' obs_idx <- which(!is.na(Y))
#' set.seed(1)
#' mask[sample(obs_idx, floor(0.2 * length(obs_idx)))] <- TRUE
#' Y_train <- Y; Y_train[mask] <- NA
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(Y_train, R = 0, family = "binary",
#'            burn = 15, nscan = 60, odens = 5, verbose = FALSE, plot = FALSE, chains = 1)
#' y_pred <- predict(fit, type = "response")
#' evaluate_heldout(Y, y_pred, mask, family = "binary")
#' }
#' @export
evaluate_heldout <- function(y_obs, y_pred, mask, family = "binary", fit = NULL) {
	# flatten inputs to matching vectors
	to_vec <- function(z) {
		if (is.list(z)) unlist(lapply(z, as.vector)) else as.vector(z)
	}
	y <- to_vec(y_obs); p <- to_vec(y_pred); m <- to_vec(mask)
	if (length(y) != length(p) || length(y) != length(m)) {
		cli::cli_abort("{.arg y_obs}, {.arg y_pred} and {.arg mask} must have the same shape.")
	}
	keep <- isTRUE(any(m)) & is.finite(y) & is.finite(p) & as.logical(m)
	if (!any(keep)) {
		cli::cli_abort("{.arg mask} selects no observed, finite cells.")
	}
	y <- y[keep]; p <- p[keep]
	n_eval <- length(y)

	if (family %in% c("binary", "cbin")) {
		# keep p inside (0, 1) so the log loss is finite
		p_c <- pmin(pmax(p, 1e-12), 1 - 1e-12)
		brier <- mean((p_c - y)^2)
		log_l <- -mean(y * log(p_c) + (1 - y) * log(1 - p_c))
		out <- data.frame(
			n_eval = n_eval,
			auroc = NA_real_,
			auprc = NA_real_,
			brier = brier,
			logloss = log_l
		)
		if (requireNamespace("precrec", quietly = TRUE)) {
			mm <- precrec::evalmod(scores = p_c, labels = y)
			auc_df <- precrec::auc(mm)
			out$auroc <- as.numeric(auc_df$aucs[auc_df$curvetypes == "ROC"])
			out$auprc <- as.numeric(auc_df$aucs[auc_df$curvetypes == "PRC"])
		} else if (requireNamespace("pROC", quietly = TRUE)) {
			rr <- pROC::roc(response = y, predictor = p_c, quiet = TRUE)
			out$auroc <- as.numeric(pROC::auc(rr))
			# no pr-auc without precrec
		} else {
			cli::cli_inform(c("i" = "Install {.pkg precrec} or {.pkg pROC} for AUROC / PR-AUC."))
		}
		return(out)
	}

	if (family == "normal") {
		err <- y - p
		rmse <- sqrt(mean(err^2))
		mae <- mean(abs(err))
		# the predictive density uses the fit's own residual variance
		ve <- if (is.null(fit)) NA_real_ else .lame_heldout_ve(fit)
		mean_logdens <- if (is.na(ve)) NA_real_ else
			mean(stats::dnorm(y, mean = p, sd = sqrt(ve), log = TRUE))
		return(data.frame(
			n_eval = n_eval,
			rmse = rmse,
			mae = mae,
			mean_logdens = mean_logdens
		))
	}

	if (family == "poisson") {
		lam <- pmin(pmax(p, 1e-12), 1e8)
		err <- y - lam
		rmse <- sqrt(mean(err^2))
		mean_logdens <- mean(stats::dpois(round(y), lambda = lam, log = TRUE))
		return(data.frame(
			n_eval = n_eval,
			rmse = rmse,
			mean_logdens = mean_logdens
		))
	}

	cli::cli_warn("Family {.val {family}} is not handled by {.fn evaluate_heldout}; returning {.code n_eval} only.")
	data.frame(n_eval = n_eval)
}

# residual variance of a normal fit: posterior mean of ve (MCMC) or fitted ve (ALS)
.lame_heldout_ve <- function(fit) {
	if (!inherits(fit, c("ame", "lame", "ame_als"))) {
		cli::cli_abort(
			"{.arg fit} must be an {.cls ame}, {.cls lame} or {.cls ame_als} fit.", class = "lame_heldout_fit"
		)
	}
	if (!identical(fit$family %||% "normal", "normal")) {
		cli::cli_abort(
			"{.arg fit} is a {.val {fit$family}} fit; the normal held-out density needs a normal fit.", class = "lame_heldout_fit"
		)
	}
	VC <- fit$VC
	ve <- if (is.matrix(VC) && "ve" %in% colnames(VC)) mean(VC[, "ve"]) else if ("ve" %in% names(VC)) unname(VC[["ve"]]) else NA_real_
	if (!is.finite(ve) || ve <= 0) {
		cli::cli_abort(
			"{.arg fit} carries no finite positive residual variance {.code ve}.", class = "lame_heldout_fit"
		)
	}
	ve
}
