# collapsed per-actor update of the binary sampler for bipartite networks, the
# rectangular counterpart of binary_actor_collapsed_move.R. with the latent array
# integrated out an observed cell contributes Phi(s eta), s = 2 y - 1, and each block
# takes an exact elliptical slice step (Murray, Adams and MacKay 2010), so actors
# without events are not held to steps of 1 / sqrt(cells) per sweep.
#
# rows then columns, in random order, with G and the hyperparameters held:
#   level translation a -> a + delta (e_i - 1/nA), beta_0 -> beta_0 + delta / nA,
#     which moves only row i's cells and keeps sum a (Gaussian prior in delta);
#   static position u_i or v_j under its N(0, I) prior;
#   under dynamic_ab, the half difference of a pair's deviation paths, N(0, Sigma / 2);
#   without centring (ame()), the ridge of the free level means and the intercept,
#     drawn from its Gaussian conditional since it leaves the likelihood unchanged.
# a step that does not terminate within 100 evaluations leaves its block in place.

# per-actor index sets: in the slice of the predictor that row i (column j)
# occupies, the positions of its observed cells and their signs (+1 event,
# -1 non-event); missing cells and absent actor-periods are left out
.lame_actor_move_context_bip <- function(Y) {
	d <- dim(Y)
	nA <- d[1L]
	nB <- d[2L]
	N <- if (length(d) >= 3L) d[3L] else 1L
	Y <- array(Y, dim = c(nA, nB, N))
	# the cells the latent draw constrains (other values are left unconstrained)
	obs <- is.finite(Y) & (Y == 0 | Y == 1)
	sgn <- array(0, dim = c(nA, nB, N))
	sgn[obs] <- 2 * Y[obs] - 1
	side <- function(i, row) {
		m <- if (row) nB else nA
		o <- if (row) matrix(obs[i, , ], m, N) else matrix(obs[, i, ], m, N)
		s <- if (row) matrix(sgn[i, , ], m, N) else matrix(sgn[, i, ], m, N)
		idx <- which(o)
		list(idx = idx, sg = as.integer(s[idx]))
	}
	list(
		nA = nA,
		nB = nB,
		N = N,
		rows = lapply(seq_len(nA), side, row = TRUE),
		cols = lapply(seq_len(nB), side, row = FALSE)
	)
}

# log prior of a level translation of actor i in its shift delta,
# -(a2 delta^2 + 2 a1 delta) / 2; with `compensate` the intercept absorbs delta / n,
# and `ic` is the intercept prior of a dynamic_beta fit (see .lame_actor_move_icpt)
.lame_actor_move_level_prior_bip <- function(i, x, v, iV0, beta, compensate, ic = NULL) {
	# the element is taken without its name (see .lame_actor_move_level_prior)
	xi <- unname(x[i])
	if (!compensate) return(c(a1 = xi / v, a2 = 1 / v))
	n <- length(x)
	lq <- .lame_actor_move_icpt(ic, iV0, beta)
	c(
		a1 = (xi - mean(x)) / v + lq[["l"]] / n,
		a2 = (n - 1) / (n * v) + lq[["q"]] / n^2
	)
}

# one visit of the collapsed move. `la`, `lb` are the additive levels, `a_mat` /
# `b_mat` the dynamic paths, `U`, `V`, `G` the static positions (NULL without a
# static multiplicative term), `dev` list(rho, sigma) of the deviation paths,
# `free_means` adds the ridge step and `ic` is the intercept of a dynamic_beta fit
.lame_actor_move_visit_bip <- function(
	EZ,
	beta,
	la,
	lb,
	a_mat,
	b_mat,
	U,
	V,
	G,
	va,
	vb,
	iV0,
	compensate,
	cx,
	move_a = TRUE,
	move_b = TRUE,
	dev = NULL,
	free_means = FALSE,
	ic = NULL
) {
	nA <- cx$nA
	nB <- cx$nB
	N <- cx$N
	EZ <- array(EZ, dim = c(nA, nB, N))
	uv <- !is.null(U) && !is.null(V) && !is.null(G) && ncol(U) > 0L && ncol(V) > 0L
	stats <- c(level_calls = 0, level_evals = 0, uv_calls = 0, uv_evals = 0, stuck = 0)
	if (!is.null(dev)) stats <- c(stats, dev_calls = 0, dev_evals = 0)
	if (free_means && compensate) stats <- c(stats, ridge_calls = 0)
	for (row in c(TRUE, FALSE)) {
		n_side <- if (row) nA else nB
		m_opp <- if (row) nB else nA
		v_side <- if (row) va else vb
		level <- (if (row) move_a else move_b) && n_side >= 2L && is.finite(v_side) && v_side > 0
		if (uv) {
			# the other side's positions mapped through G: row j holds the
			# direction along which this side's position enters the cell
			opp <- if (row) tcrossprod(V, G) else U %*% G
			r_side <- ncol(opp)
			root <- diag(r_side)
			mu0 <- numeric(r_side)
		}
		for (i in sample.int(n_side)) {
			sd_ <- if (row) cx$rows[[i]] else cx$cols[[i]]
			sg <- sd_$sg
			if (level) {
				x <- if (row) la else lb
				base <- (if (row) matrix(EZ[i, , ], m_opp, N) else matrix(EZ[, i, ], m_opp, N))[sd_$idx]
				ll_level <- function(d) sum(pnorm(sg * (base + d), log.p = TRUE))
				pr <- .lame_actor_move_level_prior_bip(i, x, v_side, iV0, beta, compensate, ic)
				ans <- .lame_ess_step(0, -pr[["a1"]] / pr[["a2"]], matrix(1 / sqrt(pr[["a2"]])), ll_level, ll_level(0))
				stats["level_calls"] <- stats["level_calls"] + 1
				stats["level_evals"] <- stats["level_evals"] + ans$evals
				stats["stuck"] <- stats["stuck"] + ans$stuck
				d <- ans$x
				if (d != 0) {
					w <- numeric(n_side)
					if (compensate) {
						w[] <- -1 / n_side
						w[i] <- 1 - 1 / n_side
						sh <- .lame_actor_move_icpt_shift(ic, beta, d / n_side)
						beta <- sh$beta; ic <- sh$ic
					} else w[i] <- 1
					if (row) {
						la <- la + w * d
						if (!is.null(a_mat)) a_mat <- a_mat + w * d
						EZ[i, , ] <- EZ[i, , ] + d
					} else {
						lb <- lb + w * d
						if (!is.null(b_mat)) b_mat <- b_mat + w * d
						EZ[, i, ] <- EZ[, i, ] + d
					}
				}
			}
			if (uv) {
				cur <- if (row) U[i, ] else V[i, ]
				partner <- ((sd_$idx - 1L) %% m_opp) + 1L
				base <- (if (row) matrix(EZ[i, , ], m_opp, N) else matrix(EZ[, i, ], m_opp, N))[sd_$idx] -
					drop(opp %*% cur)[partner]
				ll_uv <- function(x) sum(pnorm(sg * (base + drop(opp %*% x)[partner]), log.p = TRUE))
				ans <- .lame_ess_step(cur, mu0, root, ll_uv, ll_uv(cur))
				stats["uv_calls"] <- stats["uv_calls"] + 1
				stats["uv_evals"] <- stats["uv_evals"] + ans$evals
				stats["stuck"] <- stats["stuck"] + ans$stuck
				dx <- ans$x - cur
				if (any(dx != 0)) {
					shift <- drop(opp %*% dx)
					if (row) {
						U[i, ] <- ans$x
						EZ[i, , ] <- EZ[i, , ] + shift
					} else {
						V[i, ] <- ans$x
						EZ[, i, ] <- EZ[, i, ] + shift
					}
				}
			}
		}
	}
	if (free_means && compensate) {
		# ridge of the free (not centred) levels from its Gaussian conditional, no latent draw
		act <- c(move_a && is.finite(va) && va > 0, move_b && is.finite(vb) && vb > 0)
		if (any(act)) {
			m1 <- drop(iV0 %*% beta)[1L]
			v0 <- iV0[1L, 1L]
			P <- matrix(v0, 2L, 2L) + diag(c(nA / va, nB / vb), 2L)
			h <- c(m1 - sum(la) / va, m1 - sum(lb) / vb)
			P <- P[act, act, drop = FALSE]
			h <- h[act]
			ch <- tryCatch(chol(P), error = function(e) NULL)
			if (!is.null(ch) && all(is.finite(h))) {
				cc <- solve(P, h) + backsolve(ch, rnorm(length(h)))
				stats["ridge_calls"] <- stats["ridge_calls"] + 1
				if (act[1L]) { la <- la + cc[1L]; if (!is.null(a_mat)) a_mat <- a_mat + cc[1L] }
				if (act[2L]) { lb <- lb + cc[sum(act)]; if (!is.null(b_mat)) b_mat <- b_mat + cc[sum(act)] }
				beta[1L] <- beta[1L] - sum(cc)
			}
		}
	}
	if (!is.null(dev) && N >= 2L) {
		# deviation paths of a pair of actors with their sum held: half their difference
		# has prior N(0, Sigma / 2), Sigma the AR(1) covariance, and moves only their cells
		Sig <- dev$sigma^2 / (1 - dev$rho^2) * dev$rho^abs(outer(seq_len(N), seq_len(N), "-"))
		root <- tryCatch(t(chol(Sig / 2)), error = function(e) NULL)
		for (row in c(TRUE, FALSE)) {
			if (is.null(root)) break
			n_side <- if (row) nA else nB
			m_opp <- if (row) nB else nA
			if (!(if (row) move_a else move_b) || n_side < 2L || is.null(if (row) a_mat else b_mat)) next
			perm <- sample.int(n_side)
			for (q in seq_len(n_side %/% 2L)) {
				i <- perm[2L * q - 1L]
				k <- perm[2L * q]
				si <- if (row) cx$rows[[i]] else cx$cols[[i]]
				sk <- if (row) cx$rows[[k]] else cx$cols[[k]]
				bi <- (if (row) matrix(EZ[i, , ], m_opp, N) else matrix(EZ[, i, ], m_opp, N))[si$idx]
				bk <- (if (row) matrix(EZ[k, , ], m_opp, N) else matrix(EZ[, k, ], m_opp, N))[sk$idx]
				# the variable is the half difference of the two deviation paths
				mat <- if (row) a_mat else b_mat
				lev <- if (row) la else lb
				x0 <- ((mat[i, ] - lev[i]) - (mat[k, ] - lev[k])) / 2
				per_i <- ((si$idx - 1L) %/% m_opp) + 1L
				per_k <- ((sk$idx - 1L) %/% m_opp) + 1L
				ll_dev <- function(x) sum(pnorm(si$sg * (bi + (x - x0)[per_i]), log.p = TRUE)) +
					sum(pnorm(sk$sg * (bk - (x - x0)[per_k]), log.p = TRUE))
				ans <- .lame_ess_step(x0, numeric(N), root, ll_dev, ll_dev(x0))
				stats["dev_calls"] <- stats["dev_calls"] + 1
				stats["dev_evals"] <- stats["dev_evals"] + ans$evals
				stats["stuck"] <- stats["stuck"] + ans$stuck
				x <- ans$x - x0
				if (any(x != 0)) {
					shift <- rep(x, each = m_opp)
					if (row) {
						a_mat[i, ] <- a_mat[i, ] + x
						a_mat[k, ] <- a_mat[k, ] - x
						EZ[i, , ] <- EZ[i, , ] + shift
						EZ[k, , ] <- EZ[k, , ] - shift
					} else {
						b_mat[i, ] <- b_mat[i, ] + x
						b_mat[k, ] <- b_mat[k, ] - x
						EZ[, i, ] <- EZ[, i, ] + shift
						EZ[, k, ] <- EZ[, k, ] - shift
					}
				}
			}
		}
	}
	list(EZ = EZ, beta = beta, ic = ic, la = la, lb = lb, a_mat = a_mat, b_mat = b_mat, U = U, V = V, stats = stats)
}
