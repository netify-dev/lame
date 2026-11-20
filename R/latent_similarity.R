#' Similarity of latent positions over every fitted dimension
#'
#' @description
#' The cosine similarity between the nodes of one side of the latent product:
#' between its rows (\code{nodes = "row"}) or its columns
#' (\code{nodes = "col"}). Two nodes with similarity near 1 connect to the
#' other side in the same pattern, near -1 in the opposite one. It uses every
#' fitted dimension, so it does not depend on the rotation of the positions.
#' \code{uv_plot(fit, layout = "heatmap", nodes = "row")} draws it.
#'
#' @param fit A fitted \code{ame}, \code{lame} or \code{ame_als} model with
#'   R > 0.
#' @param nodes \code{"row"} (row actors; senders in a one-mode network) or
#'   \code{"col"} (column actors; receivers).
#' @param time_point For a fit with dynamic positions, the period (an index or
#'   a period name; the default is the last). Fits with static positions take
#'   no \code{time_point}.
#'
#' @details
#' With \eqn{U}, \eqn{V} the reported positions and \eqn{M} the middle matrix
#' of the fit (the identity, \code{fit$G} for bipartite fits, \code{fit$L} for
#' symmetric fits), the product is \eqn{U M V'}. The similarity of its rows is
#' \eqn{U M V'V M' U'} scaled to unit diagonal, that of its columns
#' \eqn{V M' U'U M V'} scaled the same way. Both are computed from the factors,
#' without forming the product. A node whose position is not distinguishable
#' from zero (every coordinate within 1.96 posterior sds of zero) has no
#' defined direction: it is listed in the attribute \code{"weak"} and placed
#' last in \code{"order"}. At most 2,000 nodes are compared.
#'
#' The latent product is centred (zero row and column means), so positions
#' are relative to the average node of each side. A node with no structure of
#' its own does not sit at the origin but opposite the average pattern, and
#' several such nodes look similar to each other.
#'
#' @return A symmetric matrix of similarities with the node names as dimnames,
#'   with attributes \code{"order"} (the clustering order: average linkage on
#'   one minus the similarity, the weak nodes last) and \code{"weak"} (the
#'   names of the nodes not distinguishable from zero; \code{character(0)}
#'   when there are none or the fit has no posterior sds of its positions).
#'
#' @seealso \code{\link{uv_plot}}, \code{\link{latent_positions}}
#'
#' @examples
#' \donttest{
#' data(YX_nrm)
#' # chains = 1 keeps the example short; the default runs two chains and compares them
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 3,
#'            burn = 10, nscan = 100, odens = 1, verbose = FALSE, chains = 1)
#' sim <- latent_similarity(fit, "row")
#' round(sim[attr(sim, "order")[1:5], attr(sim, "order")[1:5]], 2)
#' }
#'
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @export
latent_similarity <- function(fit, nodes, time_point = NULL) {
	if (!inherits(fit, c("ame", "lame", "ame_als"))) {
		cli::cli_abort("{.arg fit} must be an object of class {.cls ame}, {.cls lame} or {.cls ame_als}.")
	}
	if (missing(nodes)) {
		cli::cli_abort("{.arg nodes} is required: {.val row} for the rows of the latent product or {.val col} for its columns.")
	}
	nodes <- match.arg(nodes, c("row", "col"))
	.lame_require_draws(fit, "latent_similarity")
	if (.ae_product_is_withheld(fit)) {
		pw <- fit$eb$product_withheld
		cli::cli_abort(
			c(
				"This fit has no latent positions: the latent product is not identified by the data at rank {pw$rank}.",
				"i" = "Compare ranks with {.code fit$eb$product_withheld$log_marginal}, or fit {.fn ame} / {.fn lame} for the posterior of the positions."
			),
			class = "lame_product_withheld"
		)
	}
	fit <- .lame_positions(fit)
	U <- fit[["U"]]
	if (is.null(U) || length(dim(U)) < 2L || dim(U)[2L] == 0L) {
		cli::cli_abort(c(
			"This fit has no latent positions to compare.",
			"i" = "It was fitted with {.code R = 0}, or with a one-sided bipartite rank, which drops the multiplicative term."
		))
	}
	parts <- .uv_similarity_parts(fit, time_point)
	A <- if (nodes == "row") parts$U else parts$V
	if (nrow(A) > 2000L) {
		other <- if (nodes == "row") "col" else "row"
		n_other <- if (nodes == "row") nrow(parts$V) else nrow(parts$U)
		cli::cli_abort(c(
			"Similarity is computed for at most 2,000 nodes; the {nodes} side has {format(nrow(A), big.mark = ',')}.",
			"i" = "Use {.code nodes = \"{other}\"} ({format(n_other, big.mark = ',')} nodes), or {.code uv_plot(fit, layout = \"heatmap\", nodes = \"both\")} for the latent product."
		))
	}
	.uv_side_similarity(parts, nodes)
}

# positions, middle matrix and posterior sds of one period, rows with finite
# positions only
.uv_similarity_parts <- function(fit, time_point = NULL) {
	sym <- isTRUE(fit$symmetric)
	U <- fit[["U"]]
	V <- if (sym) NULL else fit[["V"]]
	U_sd <- fit[["U_sd"]]
	V_sd <- if (sym) NULL else fit[["V_sd"]]
	period <- NULL
	if (length(dim(U)) == 3L) {
		n_t <- dim(U)[3L]
		t_labels <- dimnames(U)[[3L]]
		tp <- if (is.null(time_point)) n_t else time_point
		if (identical(tp, "average")) {
			cli::cli_abort(c(
				"Time-averaged positions are not the factors of any period's product, so the similarity needs one period.",
				"i" = "Pass {.arg time_point} as a period index or name."
			))
		}
		idx <- if (is.character(tp)) match(tp, t_labels) else tp
		if (length(idx) != 1L || is.na(idx) || !idx %in% seq_len(n_t)) {
			cli::cli_abort("{.arg time_point} must be one period: an index between 1 and {n_t}, or a period name.")
		}
		idx <- as.integer(idx)
		slice <- function(A) {
			if (is.null(A) || length(dim(A)) != 3L) return(NULL)
			matrix(A[, , idx], dim(A)[1L], dim(A)[2L], dimnames = list(dimnames(A)[[1L]], NULL))
		}
		period <- if (!is.null(t_labels)) t_labels[idx] else paste0("T", idx)
		U <- slice(U)
		V <- slice(V)
		U_sd <- slice(U_sd)
		V_sd <- slice(V_sd)
	} else if (!is.null(time_point)) {
		cli::cli_abort(c(
			"{.arg time_point} applies to fits with dynamic positions; this fit has one set of positions.",
			"i" = "Drop {.arg time_point}."
		))
	}
	if (sym) {
		L <- fit[["L"]]
		if (is.matrix(L)) L <- diag(L)
		M <- if (length(L) == ncol(U) && all(is.finite(L))) diag(L, nrow = length(L)) else diag(1, ncol(U))
		V <- U
		V_sd <- U_sd
	} else {
		G <- fit[["G"]]
		M <- if (is.matrix(G) && identical(dim(G), c(ncol(U), ncol(V))) && all(is.finite(G))) {
			G
		} else diag(1, ncol(U), ncol(V))
	}
	if (is.null(rownames(U))) rownames(U) <- paste0(if (identical(fit$mode, "bipartite")) "R" else "S", seq_len(nrow(U)))
	if (is.null(rownames(V))) rownames(V) <- if (sym) rownames(U) else paste0(if (identical(fit$mode, "bipartite")) "C" else "R", seq_len(nrow(V)))
	if (!identical(dim(U_sd), dim(U))) U_sd <- NULL
	if (!identical(dim(V_sd), dim(V))) V_sd <- NULL
	# actors absent in the period have no position
	ok_u <- rowSums(!is.finite(U)) == 0L
	ok_v <- rowSums(!is.finite(V)) == 0L
	list(
		U = U[ok_u, , drop = FALSE],
		V = V[ok_v, , drop = FALSE],
		M = M,
		U_sd = if (!is.null(U_sd)) U_sd[ok_u, , drop = FALSE],
		V_sd = if (!is.null(V_sd)) V_sd[ok_v, , drop = FALSE],
		k = if (sym) ncol(U) else min(ncol(U), ncol(V)),
		period = period
	)
}

# cosine similarity between the rows (or columns) of U M V' from the factors:
# with V'V = Q Q', row i of the product has the inner products of row i of U M Q
.uv_side_similarity <- function(parts, nodes) {
	if (nodes == "row") {
		A <- parts$U
		B <- parts$V
		Mx <- parts$M
		A_sd <- parts$U_sd
	} else {
		A <- parts$V
		B <- parts$U
		Mx <- t(parts$M)
		A_sd <- parts$V_sd
	}
	eb <- eigen(crossprod(B), symmetric = TRUE)
	Q <- eb$vectors %*% diag(sqrt(pmax(eb$values, 0)), nrow = length(eb$values))
	E <- A %*% Mx %*% Q
	nrm <- sqrt(rowSums(E^2))
	zero <- nrm <= sqrt(.Machine$double.eps) * max(nrm)
	E <- E / ifelse(zero, 1, nrm)
	S <- pmin(pmax(tcrossprod(E), -1), 1)
	S[zero, ] <- NA_real_
	S[, zero] <- NA_real_
	diag(S)[!zero] <- 1
	dimnames(S) <- list(rownames(A), rownames(A))
	# no defined direction: every coordinate within 1.96 posterior sds of zero
	weak <- zero
	if (!is.null(A_sd)) {
		inside <- abs(A) <= 1.96 * A_sd
		inside[is.na(inside)] <- FALSE
		weak <- weak | rowSums(!inside) == 0L
	}
	strong <- which(!weak)
	ord <- strong
	if (length(strong) > 2L) {
		hc <- stats::hclust(stats::as.dist(1 - S[strong, strong, drop = FALSE]), method = "average")
		ord <- strong[hc$order]
	}
	attr(S, "order") <- unname(c(ord, which(weak)))
	attr(S, "weak") <- if (is.null(A_sd) && !any(zero)) character(0) else rownames(A)[weak]
	S
}
