# held-out link-prediction evaluator
#
# given a fit, a held-out mask, and the original Y, compute family-aware
# predictive scores. AUROC + PR-AUC + Brier for binary; RMSE + MAE for
# normal; mean log-likelihood for any family. user is expected to have
# either refit on a masked Y or to pass in fitted predicted probabilities
# already.

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
#'
#' @return A one-row data frame with columns appropriate to the family:
#'   \code{n_eval}, plus \code{auroc} + \code{auprc} + \code{brier} +
#'   \code{logloss} (binary / cbin), or \code{rmse} + \code{mae} +
#'   \code{mean_logdens} (normal / tobit), or \code{mean_logdens} +
#'   \code{rmse} (poisson).
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
#' fit <- ame(Y_train, R = 0, family = "binary",
#'            burn = 50, nscan = 200, odens = 5, verbose = FALSE, plot = FALSE)
#' y_pred <- predict(fit, type = "response")
#' evaluate_heldout(Y, y_pred, mask, family = "binary")
#' }
#' @export
evaluate_heldout <- function(y_obs, y_pred, mask, family = "binary") {
	# normalise inputs to one big numeric vector pair (y, p)
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
		# clamp p strictly inside (0, 1) so log-loss is finite
		p_c <- pmin(pmax(p, 1e-12), 1 - 1e-12)
		brier  <- mean((p_c - y)^2)
		log_l  <- -mean(y * log(p_c) + (1 - y) * log(1 - p_c))
		out <- data.frame(n_eval = n_eval, auroc = NA_real_,
		                  auprc = NA_real_, brier = brier, logloss = log_l)
		if (requireNamespace("precrec", quietly = TRUE)) {
			mm <- precrec::evalmod(scores = p_c, labels = y)
			auc_df <- precrec::auc(mm)
			out$auroc <- as.numeric(auc_df$aucs[auc_df$curvetypes == "ROC"])
			out$auprc <- as.numeric(auc_df$aucs[auc_df$curvetypes == "PRC"])
		} else if (requireNamespace("pROC", quietly = TRUE)) {
			rr <- pROC::roc(response = y, predictor = p_c, quiet = TRUE)
			out$auroc <- as.numeric(pROC::auc(rr))
			# no PR-AUC without precrec
		} else {
			cli::cli_inform(c("i" = "Install {.pkg precrec} or {.pkg pROC} for AUROC / PR-AUC."))
		}
		return(out)
	}

	if (family %in% c("normal", "tobit")) {
		err <- y - p
		rmse <- sqrt(mean(err^2))
		mae <- mean(abs(err))
		# crude mean log-density using residual SD; the real fit-level s2 is
		# unknown to this helper. Reported with the caveat that it's a rough
		# approximation in family = "normal" / "tobit" mode.
		sd_resid <- max(stats::sd(err), 1e-8)
		mean_logdens <- mean(stats::dnorm(y, mean = p, sd = sd_resid, log = TRUE))
		return(data.frame(n_eval = n_eval, rmse = rmse, mae = mae,
		                  mean_logdens = mean_logdens))
	}

	if (family == "poisson") {
		lam <- pmin(pmax(p, 1e-12), 1e8)
		err <- y - lam
		rmse <- sqrt(mean(err^2))
		mean_logdens <- mean(stats::dpois(round(y), lambda = lam, log = TRUE))
		return(data.frame(n_eval = n_eval, rmse = rmse,
		                  mean_logdens = mean_logdens))
	}

	cli::cli_warn("Family {.val {family}} is not handled by {.fn evaluate_heldout}; returning {.code n_eval} only.")
	data.frame(n_eval = n_eval)
}
