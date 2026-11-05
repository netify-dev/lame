# align latent positions over time to remove sign and rotation flips
# unipartite: U is [n, r, T]
align_over_time_unip <- function(U) {
	if (is.null(U) || length(dim(U)) != 3L) return(U)
	Tn <- dim(U)[3]
	if (Tn < 2L) return(U)
	U_al <- U
	for (t in 2:Tn) {
		M <- crossprod(U_al[, , t], U_al[, , t - 1])
		if (!all(is.finite(M))) {
			next
		}
		sv <- tryCatch(base::svd(M), error = function(e) NULL)
		if (!is.null(sv)) {
			Ropt <- sv$u %*% t(sv$v)
			U_al[, , t] <- U_al[, , t] %*% Ropt
		}
	}
	U_al
}

# symmetric: U [n, r, T] with U_t L U_t' fixed; only a sign change of a column
# keeps the product for any diagonal L, so each column takes the sign that best
# matches the period before it
align_over_time_sym <- function(U) {
	if (is.null(U) || length(dim(U)) != 3L) return(U)
	Tn <- dim(U)[3]
	if (Tn < 2L) return(U)
	U_al <- U
	for (t in 2:Tn) {
		U_t <- matrix(U_al[, , t], dim(U)[1])
		agree <- colSums(U_t * matrix(U_al[, , t - 1], dim(U)[1]), na.rm = TRUE)
		U_al[, , t] <- sweep(U_t, 2L, ifelse(agree < 0, -1, 1), "*")
	}
	U_al
}

# directed unipartite: U and V [n, r, T] share one rotation per period, found
# from the stacked positions, so U V' is unchanged
align_over_time_unip_pair <- function(U, V) {
	if (is.null(U) || is.null(V) || length(dim(U)) != 3L || length(dim(V)) != 3L) {
		return(list(U = U, V = V))
	}
	Tn <- dim(U)[3]
	if (Tn < 2L) return(list(U = U, V = V))
	U_al <- U; V_al <- V
	for (t in 2:Tn) {
		M <- crossprod(rbind(U_al[, , t], V_al[, , t]), rbind(U_al[, , t - 1], V_al[, , t - 1]))
		if (!all(is.finite(M))) next
		sv <- tryCatch(base::svd(M), error = function(e) NULL)
		if (is.null(sv)) next
		Q <- sv$u %*% t(sv$v)
		U_al[, , t] <- U_al[, , t] %*% Q
		V_al[, , t] <- V_al[, , t] %*% Q
	}
	list(U = U_al, V = V_al)
}

# bipartite: U [nA, rA, T], V [nB, rB, T], G [rA, rB]. one G serves every
# period, so U_t G V_t' is unchanged only under transformations that leave G
# fixed: in the singular-vector bases of G these are sign changes of a singular
# direction applied to U_t and V_t together. each period takes the sign
# changes that best match the period before it, and G is returned as given
align_over_time_bip <- function(U, V, G) {
	if (is.null(U) || is.null(V) || length(dim(U)) != 3L || length(dim(V)) != 3L) {
		return(list(U = U, V = V, G = G))
	}
	Tn <- dim(U)[3]
	if (Tn < 2L || is.null(G) || !all(is.finite(G))) return(list(U = U, V = V, G = G))
	G <- as.matrix(G)
	if (nrow(G) != dim(U)[2] || ncol(G) != dim(V)[2]) return(list(U = U, V = V, G = G))
	sg <- svd(G, nu = nrow(G), nv = ncol(G))
	k <- sum(sg$d > max(dim(G)) * max(sg$d, 0) * .Machine$double.eps)
	if (k == 0L) return(list(U = U, V = V, G = G))
	P <- sg$u[, seq_len(k), drop = FALSE]
	Q <- sg$v[, seq_len(k), drop = FALSE]
	nA <- dim(U)[1]; nB <- dim(V)[1]
	U_al <- U; V_al <- V
	for (t in 2:Tn) {
		U_t <- matrix(U_al[, , t], nA)
		V_t <- matrix(V_al[, , t], nB)
		# agreement of each singular direction with the previous period
		agree <- colSums((U_t %*% P) * (matrix(U_al[, , t - 1L], nA) %*% P)) +
			colSums((V_t %*% Q) * (matrix(V_al[, , t - 1L], nB) %*% Q))
		if (!all(is.finite(agree))) next
		flip <- which(agree < 0)
		if (length(flip) == 0L) next
		# reflections of the flipped directions; t(RU) %*% G %*% RV equals G
		RU <- diag(nrow(G)) - 2 * tcrossprod(P[, flip, drop = FALSE])
		RV <- diag(ncol(G)) - 2 * tcrossprod(Q[, flip, drop = FALSE])
		U_al[, , t] <- U_t %*% RU
		V_al[, , t] <- V_t %*% RV
	}
	list(U = U_al, V = V_al, G = G)
}

# AR(1)-type persistence of the aligned latent series
estimate_rho_from_U <- function(U) {
	if (is.null(U) || length(dim(U)) != 3L) return(NA_real_)
	Tn <- dim(U)[3]
	if (Tn < 2L) return(NA_real_)
	num <- 0; den <- 0
	for (t in 2:Tn) {
		Ut <- U[, , t]
		Ut1 <- U[, , t - 1]
		num <- num + sum(Ut * Ut1, na.rm = TRUE)
		den <- den + sum(Ut1 * Ut1, na.rm = TRUE)
	}
	if (den > 0) num / den else NA_real_
}
