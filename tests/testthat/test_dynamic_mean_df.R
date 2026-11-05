skip_on_cran()

dynamic_df_dense_reference = function(prep, mask, w, dynamic_ab, lab, rab, lbeta, kind, rbeta, level, sym) {
	nr = prep$nr;nc = prep$nc;nt = prep$Tt;p = length(mask);ta = if (dynamic_ab) nt else 1L
	names = character();map = matrix(0, p, nt)
	for (k in seq_len(p)) if (mask[k]) {map[k, ] = length(names) + seq_len(nt);names = c(names, paste(k, seq_len(nt)))} else {map[k, ] = length(names) + 1L;names = c(names, paste(k, 'static'))}
	nbeta = length(names);na = nr * ta;nb = if (sym) 0L else nc * ta;q = nbeta + na + nb;rows = list();values = numeric();wt = numeric()
	for (t in seq_len(nt)) {
		ij = which(
			is.finite(prep$Yarr[, , t]) & is.finite(w[, , t]) & w[, , t] > 0,
			arr.ind = TRUE
		);obs = matrix(FALSE, nr, nc);obs[ij] = TRUE;nn = nrow(ij);A = matrix(0, nn, q)
		A[, map[, t]] = cbind(1, matrix(prep$X[, , , t], nr * nc, p - 1)[which(obs), , drop = FALSE])
		aidx = nbeta + ((if (dynamic_ab) t else 1L) - 1L) * nr + ij[, 1];A[cbind(seq_len(nn), aidx)] = A[cbind(seq_len(nn), aidx)] + 1
		bidx = if (sym) nbeta + ((if (dynamic_ab) t else 1L) - 1L) * nr + ij[, 2] else nbeta + na + ((if (dynamic_ab) t else 1L) - 1L) * nc + ij[, 2]
		A[cbind(seq_len(nn), bidx)] = A[cbind(seq_len(nn), bidx)] + 1
		rows[[t]] = A;wt = c(wt, w[, , t][obs])
	}
	A = do.call(rbind, rows);K = crossprod(A * sqrt(wt));P = matrix(0, q, q)
	temporal = function(kind, rho) {if (kind == 'rw2') return(
		diff(diag(nt), differences = 2)
	);D = matrix(0, nt - 1, nt);for (t in 2:nt) {D[t - 1, t] = 1;D[t - 1, t - 1] = if (kind == 'ar1') -rho else -1};D}
	for (k in which(mask)) {jj = map[k, ];P[jj, jj] = P[jj, jj] + lbeta * crossprod(temporal(kind, rbeta))}
	C = matrix(0, ta * (if (sym) 1 else 2), q)
	for (t in seq_len(ta)) {C[t, nbeta + (t - 1) * nr + seq_len(nr)] = 1;if (!sym) C[ta + t, nbeta + na + (t - 1) * nc + seq_len(nc)] = 1}
	if (dynamic_ab) {
		margins = if (sym) 'a' else c('a', 'b')
		for (margin in margins) {n = if (margin == 'a') nr else nc;offset = if (margin == 'a') nbeta else nbeta + na;pres = if (margin == 'a') prep$row_presence else prep$col_presence
			for (i in seq_len(n)) {jj = offset + (seq_len(nt) - 1L) * n + i;D = temporal('ar1', rab);tr = if (is.null(pres)) rep(1, nt - 1) else as.numeric(pres[i, -1] & pres[i, -nt]);P[jj, jj] = P[jj, jj] + (if (sym) 2 else 1) * (lab * crossprod(D * sqrt(tr)) + diag(level, nt))}
		}
	}
	H = rbind(
		cbind(K + P, t(C)),
		cbind(C, matrix(0, nrow(C), nrow(C)))
	);rhs = rbind(K, matrix(0, nrow(C), q));Z = solve(H, rhs);(if (sym) 2 else 1) * sum(
		diag(Z[seq_len(q), , drop = FALSE])
	)
}

test_that("changing coefficients retain the exact unpenalized residual dimension", {
	set.seed(2183001);n = 8L;nt = 8L;p = 4L;Y = X = vector('list', nt);a = rnorm(n);b = rnorm(n);beta = matrix(rnorm(p * nt, .5, .1), p, nt)
	for (t in seq_len(nt)) {
		x = array(rnorm(n * n * p), c(n, n, p));eta = outer(a, b, '+') + apply(
			sweep(x, 3, beta[, t], '*'),
			c(1, 2),
			sum
		)
		y = eta + matrix(rnorm(n * n), n);diag(y) = NA_real_;dimnames(y) = list(
			paste0('a', 1:n),
			paste0('a', 1:n)
		);dimnames(x) = c(
			dimnames(y),
			list(paste0('x', 1:p))
		);Y[[t]] = y;X[[t]] = x
	}
	design = do.call(rbind, lapply(seq_len(nt), function(t) {
		ij = which(is.finite(Y[[t]]), arr.ind = TRUE);period = matrix(0, nrow(ij), (p + 1) * nt)
		period[, ((t - 1) * (p + 1) + 1):(t * (p + 1))] = cbind(
			1,
			matrix(X[[t]], n * n, p)[which(is.finite(Y[[t]])), , drop = FALSE]
		)
		cbind(period, diag(n)[ij[, 1], -1, drop = FALSE], diag(n)[ij[, 2], -1, drop = FALSE])
	}))
	y = unlist(Y);keep = is.finite(y);y = y[keep];reference = qr(design);mean = qr.fitted(reference, y)
	fit = lame(
		Y,
		Xdyad = X,
		R = 0,
		family = 'normal',
		method = 'als',
		dynamic_beta = TRUE,
		bias_correct = 'none',
		prior = list(lambda_beta_als = 0),
		als_max_iter = 1000,
		als_tol = 1e-10,
		verbose = FALSE
	)
	expect_true(fit$converged)
	expect_equal(unname(unlist(fit$EZ)[keep]), unname(mean), tolerance = 1e-5)
	expect_equal(fit$df_mean, reference$rank, tolerance = 1e-6)
	expect_equal(unname(fit$VC['ve']), sum((y - mean)^2) / (length(y) - reference$rank), tolerance = 1e-6)
})

test_that("joint path information matches independent constrained coordinates", {
	configs = expand.grid(
		mode = c('directed', 'bipartite', 'symmetric'),
		dynamic_ab = c(FALSE, TRUE),
		variant = 1:2,
		stringsAsFactors = FALSE
	)
	for (k in seq_len(nrow(configs))) {
		g = configs[k, ];set.seed(2185000 + k);nr = 7L;nc = if (g$mode == 'bipartite') 5L else nr;nt = 4L;sym = g$mode == 'symmetric'
		prep = list(
			nr = nr,
			nc = nc,
			Tt = nt,
			p_dyad = 2L,
			x_names_dyad = c('x1', 'x2'),
			bip = g$mode == 'bipartite',
			Yarr = array(0, c(nr, nc, nt)),
			X = array(rnorm(nr * nc * 2 * nt), c(nr, nc, 2, nt)),
			row_presence = matrix(TRUE, nr, nt),
			col_presence = matrix(TRUE, nc, nt)
		)
		if (!prep$bip) for (t in seq_len(nt)) diag(prep$Yarr[, , t]) = NA_real_
		if (sym) for (t in seq_len(nt)) for (j in 1:2) prep$X[, , j, t] = (prep$X[, , j, t] + t(prep$X[, , j, t])) / 2
		if (g$variant == 2L && g$dynamic_ab) {
			prep$Yarr[nr, , 2] = NA_real_;prep$row_presence[nr, 2] = FALSE
			if (!prep$bip) {prep$Yarr[, nc, 2] = NA_real_;prep$col_presence[nc, 2] = FALSE}
		}
		weight = array(if (g$variant == 1L) 1 else exp(rnorm(nr * nc * nt, 0, .5)), c(nr, nc, nt))
		if (sym) for (t in seq_len(nt)) weight[, , t] = (weight[, , t] + t(weight[, , t])) / 2
		mask = if (g$variant == 1L) c(TRUE, TRUE, TRUE) else c(TRUE, FALSE, TRUE)
		bd = list(mask = mask, any = any(mask));lambda = if (g$variant == 1L) 0 else 10;level = if (g$variant == 2L && g$dynamic_ab) 1 else 0
		kind = if (k %% 3 == 0) 'rw2' else if (k %% 3 == 1) 'rw1' else 'ar1'
		expected = dynamic_df_dense_reference(
			prep,
			mask,
			weight,
			g$dynamic_ab,
			lambda,
			.7,
			lambda,
			kind,
			.8,
			level,
			sym
		)
		actual = lame:::.lda_joint_mean_df(
			prep,
			c('intercept', 'x1', 'x2'),
			bd,
			weight,
			g$dynamic_ab,
			lambda,
			.7,
			lambda,
			kind,
			.8,
			level,
			sym
		)
		expect_equal(actual, expected, tolerance = 2e-6, info = paste(g$mode, g$dynamic_ab, g$variant))
	}
})

test_that("dynamic count dispersion uses the remaining mean dimension", {
	set.seed(2185991);n = 24L;nt = 6L;Y = X = vector('list', nt)
	for (t in seq_len(nt)) {
		x = matrix(rnorm(n * n), n);eta = .5 + .4 * x + outer(rnorm(n, 0, .3), rnorm(n, 0, .3), '+')
		y = matrix(rpois(n * n, exp(eta + matrix(rnorm(n * n, 0, sqrt(.5)), n))), n);diag(y) = NA_real_
		dimnames(y) = list(paste0('a', 1:n), paste0('a', 1:n));Y[[t]] = y;X[[t]] = array(
			x,
			c(n, n, 1),
			dimnames = c(dimnames(y), list('x'))
		)
	}
	# the raw Pearson root, before the jackknife correction tested separately
	fit = lame(
		Y,
		Xdyad = X,
		R = 0,
		family = 'poisson',
		method = 'als',
		dynamic_ab = TRUE,
		bias_correct = 'none',
		als_max_iter = 1000,
		als_tol = 1e-6,
		verbose = FALSE
	)
	expect_true(fit$converged)
	expect_gt(fit$df_mean, 2 + 2 * (n - 1))
	variance = unname(fit$VC['ve']);expect_gt(variance, 0)
	y = unlist(Y);mu = unlist(fit$fitted);ok = is.finite(y)
	pearson = sum((y[ok] - mu[ok])^2 / (mu[ok] + expm1(variance) * mu[ok]^2))
	expect_equal(pearson, sum(ok) - fit$df_mean, tolerance = 1e-6)
})
