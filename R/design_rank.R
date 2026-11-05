# refuse a pooled design (observed cells, all periods) with a column that is an
# exact linear combination of earlier ones, to relative tolerance tol; nearly
# collinear designs are fitted. each period is reduced to its triangular factor
# so only one period's design is held in memory
.lame_check_design_rank <- function(xlist, y, tol = 1e-8) {
	if (is.array(xlist) && length(dim(xlist)) == 3L) xlist <- list(xlist)
	if (is.null(xlist) || !length(xlist)) return(invisible(NULL))
	if (length(dim(y)) == 2L) y <- array(y, c(dim(y), 1L))
	p <- dim(xlist[[1L]])[3L]
	if (is.null(p) || p == 0L) return(invisible(NULL))
	nm <- dimnames(xlist[[1L]])[[3L]]
	if (is.null(nm) || length(nm) != p) nm <- paste0("column ", seq_len(p))
	n_obs <- 0L
	blocks <- lapply(seq_along(xlist), function(t) {
		x <- matrix(xlist[[t]], ncol = p)
		keep <- !is.na(c(y[, , t])) & rowSums(!is.finite(x)) == 0L
		if (!any(keep)) return(NULL)
		n_obs <<- n_obs + sum(keep)
		q <- qr(x[keep, , drop = FALSE], LAPACK = TRUE)
		qr.R(q)[, order(q$pivot), drop = FALSE]
	})
	r <- do.call(rbind, blocks)
	if (is.null(r)) return(invisible(NULL))
	scale <- sqrt(colSums(r^2))
	zero <- scale == 0
	scale[zero] <- 1
	r <- sweep(r, 2L, scale, "/")
	q <- qr(r, tol = tol)
	if (q$rank == p) return(invisible(NULL))
	aliased <- sort(q$pivot[(q$rank + 1L):p])
	kept <- setdiff(seq_len(p), aliased)
	# each aliased column as the combination of the others it equals
	relation <- vapply(aliased, function(a) {
		if (zero[a]) return(paste0(nm[a], " is zero on every observed cell"))
		if (!length(kept)) return(paste0(nm[a], " has no other column to be told from"))
		b <- qr.coef(qr(r[, kept, drop = FALSE]), r[, a])
		b <- b * scale[a] / scale[kept]
		b[!is.finite(b) | abs(b) < 1e-8] <- 0
		use <- which(b != 0)
		if (!length(use)) return(paste0(nm[a], " is zero on every observed cell"))
		use <- use[order(-round(abs(b[use]), 6L), use)][seq_len(min(4L, length(use)))]
		paste0(nm[a], " = ", paste(paste(signif(b[use], 4L), nm[kept][use]), collapse = " + "))
	}, character(1L))
	relation <- gsub("+ -", "- ", relation, fixed = TRUE)
	one <- length(aliased) == 1L
	cols <- paste0("`", nm[aliased], "`", collapse = ", ")
	header <- paste0(
		"The regression design is exactly rank-deficient: ",
		cols,
		if (one) " is" else " are",
		" aliased with the other columns."
	)
	detail <- paste0(
		"On the ",
		n_obs,
		" observed cell",
		if (n_obs != 1L) "s",
		" of ",
		length(xlist),
		" period",
		if (length(xlist) != 1L) "s",
		" the design has rank ",
		q$rank,
		" for its ",
		p,
		if (p == 1L) " column" else " columns",
		", so the data cannot separate ",
		if (one) "that coefficient" else "those coefficients",
		" from the others and the posterior along ",
		if (one) "that direction" else "those directions",
		" would be the prior's alone."
	)
	# escape braces in covariate names so cli does not read them as markup
	esc <- function(x) gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
	cli::cli_abort(
		c(
			esc(header),
			stats::setNames(esc(relation), rep("x", length(relation))),
			"i" = esc(detail),
			"i" = "Drop the aliased column or the covariates it duplicates. A design that is only nearly collinear is fitted."
		),
		class = "lame_design_rank",
		aliased = nm[aliased],
		rank = q$rank,
		n_col = p
	)
}
