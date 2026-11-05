# trace of the joint coefficient / additive smoother on its centring
# constraints; the coefficient paths border the block-tridiagonal actor paths
.lda_mean_hat_trace <- function(
	K,
	P,
	n_beta,
	n_row,
	n_col,
	n_time,
	dynamic_ab,
	symmetric
) {
	q <- nrow(K)
	H <- Matrix::forceSymmetric(K + P + Matrix::Diagonal(q, 1e-8))
	inverse <- function(x) {
		x <- as.matrix(x)
		chol2inv(chol((x + t(x)) / 2))
	}
	block_trace <- function() {
		if (!dynamic_ab || q == n_beta) return(sum(inverse(H) * as.matrix(K)))
		n_a <- (n_row - 1L) * n_time
		blocks <- lapply(seq_len(n_time), function(tt) {
			c(
				n_beta + (tt - 1L) * (n_row - 1L) + seq_len(n_row - 1L),
				if (!symmetric) n_beta + n_a + (tt - 1L) * (n_col - 1L) +
					seq_len(n_col - 1L)
			)
		})
		border <- seq_len(n_beta)
		forward <- off <- rhs <- solution <- vector("list", n_time)
		data_diag <- data_border <- vector("list", n_time)
		for (tt in seq_len(n_time)) {
			index <- blocks[[tt]]
			A <- as.matrix(H[index, index, drop = FALSE])
			rhs[[tt]] <- as.matrix(H[index, border, drop = FALSE])
			data_diag[[tt]] <- as.matrix(K[index, index, drop = FALSE])
			data_border[[tt]] <- as.matrix(K[index, border, drop = FALSE])
			if (tt > 1L) {
				off[[tt - 1L]] <- as.matrix(H[blocks[[tt - 1L]], index, drop = FALSE])
				B <- off[[tt - 1L]]
				A <- A - crossprod(B, forward[[tt - 1L]] %*% B)
				rhs[[tt]] <- rhs[[tt]] - crossprod(
					B,
					forward[[tt - 1L]] %*% rhs[[tt - 1L]]
				)
			}
			forward[[tt]] <- inverse(A)
		}
		solution[[n_time]] <- forward[[n_time]] %*% rhs[[n_time]]
		D <- forward[[n_time]]
		total <- sum(D * data_diag[[n_time]])
		if (n_time > 1L) for (tt in (n_time - 1L):1L) {
			solution[[tt]] <- forward[[tt]] %*%
				(rhs[[tt]] - off[[tt]] %*% solution[[tt + 1L]])
			FB <- forward[[tt]] %*% off[[tt]]
			D <- forward[[tt]] + FB %*% D %*% t(FB)
			total <- total + sum(D * data_diag[[tt]])
		}
		S <- as.matrix(H[border, border, drop = FALSE])
		Q <- as.matrix(K[border, border, drop = FALSE])
		for (tt in seq_len(n_time)) {
			Z <- solution[[tt]]
			S <- S - crossprod(as.matrix(H[blocks[[tt]], border, drop = FALSE]), Z)
			cross <- crossprod(Z, data_border[[tt]])
			Q <- Q + crossprod(Z, data_diag[[tt]] %*% Z) - cross - t(cross)
		}
		total + sum(inverse(S) * Q)
	}
	# the full sparse solve also handles a nearly singular border without
	# subtracting two large information matrices in the Schur complement
	sparse_trace <- function() {
		factor <- Matrix::Cholesky(H, LDL = FALSE, perm = TRUE)
		total <- 0
		for (start in seq.int(1L, q, by = 256L)) {
			index <- seq.int(start, min(start + 255L, q))
			z <- as.matrix(Matrix::solve(
				factor,
				as.matrix(K[, index, drop = FALSE]),
				system = "A"
			))
			total <- total + sum(z[cbind(index, seq_along(index))])
		}
		total
	}
	value <- tryCatch(block_trace(), error = function(e) sparse_trace())
	if (!is.finite(value)) {
		cli::cli_abort(
			"The effective dimension of the fitted mean could not be evaluated.",
			class = "lame_als_numerical"
		)
	}
	value
}

.lda_joint_mean_df <- function(
	prep,
	coef_names,
	beta_dyn,
	weights,
	dynamic_ab,
	lambda_ab,
	rho_ab,
	lambda_beta,
	beta_kind,
	rho_beta,
	level_prec_ab = 0,
	symmetric = FALSE
) {
	nt <- prep$Tt
	nr <- prep$nr
	nc <- prep$nc
	p <- length(coef_names)
	index <- .lda_coef_index(beta_dyn$mask, nt)
	ta <- if (dynamic_ab) nt else 1L
	contrast <- function(n) {
		Matrix::sparseMatrix(
			i = c(seq_len(n - 1L), rep(n, n - 1L)),
			j = rep(seq_len(n - 1L), 2L),
			x = c(rep(1, n - 1L), rep(-1, n - 1L)),
			dims = c(n, n - 1L)
		)
	}
	qa <- contrast(nr)
	qb <- if (symmetric) NULL else contrast(nc)
	na <- (nr - 1L) * ta
	nb <- if (symmetric) 0L else (nc - 1L) * ta
	q <- index$n + na + nb
	k_i <- k_j <- k_x <- p_i <- p_j <- p_x <- list()
	add_block <- function(indices, value, penalty = FALSE) {
		if (inherits(value, "Matrix")) {
			entries <- Matrix::summary(Matrix::forceSymmetric(value))
			ii <- indices[entries$i]
			jj <- indices[entries$j]
			xx <- entries$x
		} else {
			entries <- which(upper.tri(value, diag = TRUE) & value != 0, arr.ind = TRUE)
			ii <- indices[entries[, 1L]]
			jj <- indices[entries[, 2L]]
			xx <- value[entries]
		}
		low <- pmin(ii, jj)
		high <- pmax(ii, jj)
		if (penalty) {
			p_i[[length(p_i) + 1L]] <<- low
			p_j[[length(p_j) + 1L]] <<- high
			p_x[[length(p_x) + 1L]] <<- xx
		} else {
			k_i[[length(k_i) + 1L]] <<- low
			k_j[[length(k_j) + 1L]] <<- high
			k_x[[length(k_x) + 1L]] <<- xx
		}
	}
	project <- function(M) {
		r <- nrow(M)
		c <- ncol(M)
		sweep(
			sweep(M[-r, -c, drop = FALSE], 1L, M[-r, c], "-"),
			2L,
			M[r, -c],
			"-"
		) + M[r, c]
	}
	for (tt in seq_len(nt)) {
		observed <- is.finite(prep$Yarr[, , tt]) &
			is.finite(weights[, , tt]) & weights[, , tt] > 0
		if (!any(observed)) next
		cells <- which(observed, arr.ind = TRUE)
		X <- .lda_design_obs(prep, tt, observed, coef_names, symmetric = symmetric)
		aj <- index$n + ((if (dynamic_ab) tt else 1L) - 1L) * (nr - 1L) +
			seq_len(nr - 1L)
		bj <- if (symmetric) integer(0) else index$n + na +
			((if (dynamic_ab) tt else 1L) - 1L) * (nc - 1L) + seq_len(nc - 1L)
		cj <- vapply(seq_len(p), function(k) {
			if (beta_dyn$mask[k]) index$dynamic_map[k, tt] else index$static_map[k]
		}, integer(1))
		W <- weights[, , tt]
		W[!observed] <- 0
		ca <- matrix(0, nr, p)
		cb <- matrix(0, nc, p)
		sa <- rowsum(X * W[observed], cells[, 1L], reorder = FALSE)
		sb <- rowsum(X * W[observed], cells[, 2L], reorder = FALSE)
		ca[as.integer(rownames(sa)), ] <- sa
		cb[as.integer(rownames(sb)), ] <- sb
		A <- if (symmetric) {
			project(diag(rowSums(W) + colSums(W)) + W + t(W))
		} else diag(rowSums(W)[-nr], nr - 1L) + rowSums(W)[nr]
		across <- if (symmetric) ca + cb else ca
		coef_a <- t(sweep(across[-nr, , drop = FALSE], 2L, across[nr, ], "-"))
		information <- rbind(
			cbind(crossprod(X * sqrt(W[observed])), coef_a),
			cbind(t(coef_a), A)
		)
		if (!symmetric) {
			B <- diag(colSums(W)[-nc], nc - 1L) + colSums(W)[nc]
			coef_b <- t(sweep(cb[-nc, , drop = FALSE], 2L, cb[nc, ], "-"))
			cross <- rbind(coef_b, project(W))
			information <- rbind(cbind(information, cross), cbind(t(cross), B))
		}
		add_block(c(cj, aj, bj), information)
	}
	if (length(index$dynamic_idx) && lambda_beta > 0) {
		for (k in index$dynamic_idx) {
			add_block(
				index$dynamic_map[k, ],
				.lda_path_penalty(nt, lambda_beta * .lda_beta_ms(beta_dyn, k), beta_kind, rho_beta),
				TRUE
			)
		}
	}
	actor_penalty <- function(Q, presence) {
		n <- nrow(Q)
		np <- ncol(Q)
		if (!dynamic_ab) return(Matrix::Matrix(0, np, np, sparse = TRUE))
		transition <- .lda_transition_from_presence(presence)
		D <- matrix(0, max(nt - 1L, 0L), nt)
		if (nt > 1L) for (tt in 2:nt) {
			D[tt - 1L, tt] <- 1
			D[tt - 1L, tt - 1L] <- -rho_ab
		}
		F <- kronecker(Matrix::Matrix(D, sparse = TRUE), Q)
		w <- if (is.null(transition)) rep(1, n * (nt - 1L)) else
			as.vector(transition[, -1L, drop = FALSE])
		F <- F * sqrt(w)
		lambda_ab * Matrix::crossprod(F) + level_prec_ab *
			kronecker(Matrix::Diagonal(nt), Matrix::crossprod(Q))
	}
	add_block(
		index$n + seq_len(na),
		(if (symmetric) 2 else 1) *
			actor_penalty(qa, prep$row_presence),
		TRUE
	)
	if (!symmetric) add_block(
		index$n + na + seq_len(nb),
		actor_penalty(qb, prep$col_presence),
		TRUE
	)
	assemble <- function(i, j, x) {
		Matrix::sparseMatrix(
			i = as.integer(unlist(i)),
			j = as.integer(unlist(j)),
			x = as.numeric(unlist(x)),
			dims = c(q, q),
			symmetric = TRUE
		)
	}
	K <- assemble(k_i, k_j, k_x)
	P <- assemble(p_i, p_j, p_x)
	(if (symmetric) 2 else 1) * .lda_mean_hat_trace(
		K,
		P,
		index$n,
		nr,
		nc,
		nt,
		dynamic_ab,
		symmetric
	)
}
