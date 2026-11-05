test_that("centered static actor solves agree with independent constrained equations", {
	for (sym in c(FALSE, TRUE)) for (replication in 1:3) {
		set.seed(2693000 + 100 * sym + replication);n = 7L;nc = if (sym) n else 5L;nt = 4L
		Z = array(
			rnorm(n * nc * nt),
			c(n, nc, nt)
		);W = array(exp(rnorm(length(Z))), dim(Z));Xb = array(rnorm(length(Z)), dim(Z));O = array(rnorm(length(Z)), dim(Z));Z[runif(length(Z)) < .1] = NA
		obs = which(is.finite(Z), arr.ind = TRUE);np = if (sym) n else n + nc;D = matrix(0, nrow(obs), np)
		for (j in seq_len(nrow(obs))) {D[j, obs[j, 1]] = D[j, obs[j, 1]] + 1;idx = if (sym) obs[j, 2] else n + obs[j, 2];D[j, idx] = D[j, idx] + 1}
		C = if (sym) matrix(1, 1, n) else rbind(c(rep(1, n), rep(0, nc)), c(rep(0, n), rep(1, nc)))
		H = crossprod(D * sqrt(W[obs]));K = rbind(cbind(H, t(C)), cbind(C, matrix(0, nrow(C), nrow(C))))
		ref = solve(K, c(crossprod(D, W[obs] * (Z - Xb - O)[obs]), rep(0, nrow(C))))[seq_len(np)]
		got = lame:::.lda_update_static_ab_centered(Z, W, Xb, O, sym);theta = if (sym) got$a[, 1] else c(got$a[, 1], got$b[, 1])
		expect_equal(theta, ref, tolerance = 1e-8)
		expect_lt(max(abs(colSums(got$a))), 1e-10)
		expect_lt(max(abs(colSums(got$b))), 1e-10)
	}
})

test_that("penalized changing intercepts and static actors reach the joint Gaussian optimum", {
	skip_on_cran()
	for (mode in c('directed', 'bipartite', 'symmetric')) {
		set.seed(2693951);nr = 8L;nc = if (mode == 'bipartite') 6L else nr;nt = 5L;sym = mode == 'symmetric';Y = X = vector('list', nt)
		a = rnorm(nr, 0, .3);b = if (sym) a else rnorm(nc, 0, .3)
		for (t in seq_len(nt)) {
			x = matrix(rnorm(nr * nc), nr);noise = matrix(rnorm(nr * nc), nr)
			if (sym) {x = (x + t(x)) / sqrt(2);noise[lower.tri(noise)] = t(noise)[lower.tri(noise)]}
			y = .8 + .1 * t + (.5 + .05 * t) * x + outer(a, b, '+') + noise;if (mode != 'bipartite') diag(y) = NA_real_
			dimnames(y) = list(
				paste0('r', 1:nr),
				if (mode == 'bipartite') paste0('c', 1:nc) else paste0('r', 1:nc)
			);Y[[t]] = y;X[[t]] = array(x, c(nr, nc, 1), dimnames = c(dimnames(y), list('x')))
		}
		fit = lame(
			Y,
			Xdyad = X,
			R = 0,
			family = 'normal',
			mode = if (mode == 'bipartite') 'bipartite' else 'unipartite',
			symmetric = sym,
			method = 'als',
			dynamic_beta = TRUE,
			als_max_iter = 200,
			als_tol = 1e-10,
			verbose = FALSE,
			bias_correct = 'none'
		)
		nbeta = 2 * nt;q = nbeta + nr + if (sym) 0L else nc;rows = list();values = numeric()
		for (t in seq_len(nt)) {
			ij = which(is.finite(Y[[t]]), arr.ind = TRUE);N = nrow(ij);D = matrix(0, N, q)
			D[, c(t, nt + t)] = cbind(1, matrix(X[[t]], nr, nc)[ij]);D[cbind(1:N, nbeta + ij[, 1])] = 1
			jj = nbeta + if (sym) ij[, 2] else nr + ij[, 2];D[cbind(1:N, jj)] = D[cbind(1:N, jj)] + 1
			rows[[t]] = D;values = c(values, Y[[t]][ij])
		}
		D = do.call(rbind, rows);P = matrix(0, q, q);transition = matrix(0, nt - 1, nt)
		for (t in 2:nt) {transition[t - 1, t] = 1;transition[t - 1, t - 1] = -fit$rho_beta}
		# smoothness weight of a moving coefficient: lambda_beta times the column's mean square (one for the intercept)
		ms = c(1, mean(unlist(lapply(seq_len(nt), function(t) {ok = is.finite(Y[[t]]);matrix(X[[t]], nr, nc)[ok]^2}))))
		for (k in 0:1) {jj = k * nt + seq_len(nt);P[jj, jj] = fit$lambda_beta_als * ms[k + 1] * crossprod(transition)}
		C = matrix(0, if (sym) 1L else 2L, q);C[1, nbeta + seq_len(nr)] = 1;if (!sym) C[2, nbeta + nr + seq_len(nc)] = 1
		K = rbind(cbind(crossprod(D) + P, t(C)), cbind(C, matrix(0, nrow(C), nrow(C))))
		ref = solve(K, c(crossprod(D, values), rep(0, nrow(C))))[seq_len(q)]
		expect_true(fit$converged, info = mode)
		expect_true(fit$diagnostics$stationarity$pass, info = mode)
		expect_lt(fit$iterations, 200)
		expect_equal(unname(unlist(fit$EZ)[is.finite(unlist(Y))]), as.vector(D %*% ref), tolerance = 2e-5, info = mode)
	}
})
