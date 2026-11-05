skip_on_cran()

test_that('marginal actor-factor blocks preserve the independently evaluated posterior', {
	skip_on_cran()
	bvn = .lame_binary_dynamic_bvn_log
	actor_factor_block_context = .lame_binary_actor_slice_context
	actor_factor_block_apply = .lame_binary_actor_slice_apply
	actor_factor_marginal_system = .lame_binary_actor_slice_system
	actor_factor_marginal_loglik = .lame_binary_actor_slice_loglik
	actor_factor_marginal_slice = .lame_binary_actor_slice_step
	set.seed(921742);n = 4L;nt = 2L;Sab = matrix(c(.3, .08, .08, .25), 2);records = list();stationarity = list()
	transpose = function(z) aperm(z, c(2, 1, 3))
	for (rank in 1:2) for (side in 1:2) for (cap in c(.65, Inf)) {
		rho = if (side == 1) .45 else -.3
		U = array(
			runif(n * rank * nt, -.2, .2),
			c(n, rank, nt)
		);V = array(runif(n * rank * nt, -.2, .2), c(n, rank, nt))
		levels = matrix(rnorm(n * 2, sd = .1), n, 2);levels = sweep(levels, 2, colMeans(levels), '-')
		A = B = matrix(
			rnorm(n * nt, sd = .15),
			n,
			nt
		);A = sweep(A, 2, colMeans(A), '-') + levels[, 1];B = sweep(B, 2, colMeans(B), '-') + levels[, 2]
		base = array(rnorm(n * n * nt, sd = .1), c(n, n, nt))
		predictor = function(U, V, A, B) {z = base;for (tt in 1:nt) z[, , tt] = z[, , tt] + outer(A[, tt], B[, tt], '+') + matrix(U[, , tt], n, rank) %*% t(matrix(V[, , tt], n, rank));z}
		eta = predictor(U, V, A, B);Y = array(NA_real_, c(n, n, nt));Y[2, 1, 1] = 1;Y[1, 2, 1] = 0;Y[4, 3, 2] = 1;Y[2, 4, 2] = 0
		pos = if (side == 1) U else V;opp = if (side == 1) V else U;mar = if (side == 1) A else B;lev = levels[, side];other = levels[, 3 - side]
		yy = if (side == 1) Y else transpose(Y);ez = if (side == 1) eta else transpose(eta)
		cx = actor_factor_block_context(rank, nt, Sab, side, .3, .25, .65, .22)
		ss = actor_factor_marginal_system(pos, opp, mar, lev, other, yy, ez, c(2L, 4L), rho, cx)
		full = function(ee) {ll = 0;for (tt in 1:nt) for (i in 1:(n - 1)) for (j in (i + 1):n) {
			y1 = Y[i, j, tt];y2 = Y[j, i, tt]
			if (is.finite(y1) && is.finite(y2)) ll = ll + log(
				pbivnorm::pbivnorm(
					(2 * y1 - 1) * ee[i, j, tt],
					(2 * y2 - 1) * ee[j, i, tt],
					rho * (2 * y1 - 1) * (2 * y2 - 1)
				)
			) else if (is.finite(y1)) ll = ll + pnorm(
				(2 * y1 - 1) * ee[i, j, tt],
				log.p = TRUE
			) else if (is.finite(y2)) ll = ll + pnorm((2 * y2 - 1) * ee[j, i, tt], log.p = TRUE)
		};ll}
		loglik = function(x) actor_factor_marginal_loglik(
			ss$offset1 + drop(ss$D1 %*% x),
			ss$offset2 + drop(ss$D2 %*% x),
			ss,
			bvn
		)
		initial = loglik(ss$current);full_initial = full(eta);errors = numeric(50);centering = numeric(50)
		for (k in seq_along(errors)) {
			xx = ss$prior_mean + drop(cx$root %*% rnorm(cx$p));zz = actor_factor_block_apply(
				xx,
				ss,
				pos,
				opp,
				mar,
				lev,
				ez
			)
			correct = if (side == 1) predictor(zz$position, V, zz$margin, B) else predictor(U, zz$position, A, zz$margin)
			errors[k] = max(
				abs((loglik(xx) - initial) - (full(correct) - full_initial)),
				abs(zz$eta - (if (side == 1) correct else transpose(correct)))
			)
			centering[k] = max(abs(colSums(zz$margin) - colSums(mar)), abs(sum(zz$level) - sum(lev)))
		}
		stopifnot(max(errors) < 1e-10, max(centering) < 1e-12, max(abs(opp)) <= cap)
		N = 6000L;delta = matrix(0, N, 7);weights = numeric(N);evals = bounds = numeric(N)
		summary = function(x) c(
			x[1],
			mean(x[cx$ix$dev]),
			mean(x[cx$ix$first]),
			mean(x[cx$ix$second]),
			mean(x[cx$ix$first]^2),
			mean(x[cx$ix$second]^2),
			x[1]^2
		)
		for (k in 1:N) {
			repeat{xx = ss$prior_mean + drop(cx$root %*% rnorm(cx$p));if (max(
				abs(xx[c(cx$ix$first, cx$ix$second)])
			) <= cap) break}
			ss$current = xx;weights[k] = exp(loglik(xx));z = actor_factor_marginal_slice(ss, max_abs = cap, bvn = bvn)
			delta[k, ] = (summary(z$x) - summary(xx)) * weights[k];evals[k] = z$evaluations;bounds[k] = z$bounds
		}
		z = colMeans(delta) / (apply(delta, 2, sd) / sqrt(N));ess = sum(weights)^2 / sum(weights^2)
		rec = list(
			rank = rank,
			side = side,
			cap = cap,
			likelihood_error = max(errors),
			centering_error = max(centering),
			stationarity_z = z,
			importance_ess = ess,
			max_evaluations = max(evals),
			bound_rejections = sum(bounds),
			pass = max(abs(z)) < 4.5 && ess > 1000 && all(is.finite(z))
		)
		records[[length(records) + 1L]] = rec;stopifnot(rec$pass)
	}
	expect_true(all(vapply(records, function(z) z$pass, logical(1))))
})

.actor_slice_fit_args = function() {
	set.seed(538135);n = 6L;Y = lapply(
		1:3,
		function(t) {y = matrix(
			rbinom(n * n, 1, .3),
			n,
			n,
			dimnames = list(as.character(1:n), as.character(1:n))
		);diag(y) = NA;y[1, 2] = 1;y[2, 1] = 0;y}
	)
	list(
		Y = Y,
		family = 'binary',
		R = 2L,
		dynamic_uv = TRUE,
		dynamic_ab = TRUE,
		burn = 40L,
		nscan = 80L,
		odens = 10L,
		seed = 228351L,
		verbose = FALSE,
		gof = FALSE,
		plot = FALSE,
		posterior_opts = list(save_UV = TRUE, thin_UV = 1)
	)
}
test_that('actor slice rounds refresh binary auxiliaries and preserve checkpoint state', {
	skip_on_cran()
	args = .actor_slice_fit_args();f = suppressMessages(do.call(lame, args));health = f$actor_factor_slice
	expect_equal(unname(health['rounds']), 5)
	expect_equal(unname(health['blocks']), 30)
	expect_equal(unname(health['numerical_rejections']), 0)
	expect_gt(unname(health['evaluations']), 0)
	state = f$start_vals$checkpoint_state$values
	expect_true(all(f$sampler_health$numerical_failures == 0))
	Y = f$Y;obs = is.finite(Y);for (t in 1:3) diag(obs[, , t]) = FALSE
	expect_true(all(state$Z[obs] * (2 * Y[obs] - 1) > 0))
	expect_true(all(is.finite(state$Z)))
	expect_lt(
		max(abs(colMeans(state$a_mat) - mean(state$a_lev)), abs(colMeans(state$b_mat) - mean(state$b_lev))),
		1e-10
	)
	ck = tempfile(fileext = '.rds');on.exit(unlink(ck))
	args$start_vals = f$start_vals;args$burn = 0L;args$nscan = 40L;args$checkpoint_path = ck;args$max_seconds = 1e-12
	suppressWarnings(suppressMessages(do.call(lame, args)));saved = readRDS(ck)$start_vals$checkpoint_state$values
	for (nm in c(
		'U_cube', 'V_cube', 'a_mat', 'b_mat', 'a_lev', 'b_lev', 'Z', 'rho_uv', 'sigma_uv'
	)) expect_identical(saved[[nm]], state[[nm]], info = nm)
	args$max_seconds = NULL;args$checkpoint_path = NULL
	resumed = suppressMessages(do.call(lame, args))
	expect_equal(unname(resumed$actor_factor_slice['rounds']), 2)
	expect_equal(unname(resumed$actor_factor_slice['numerical_rejections']), 0)
	args = .actor_slice_fit_args();args$prior = list(actor_factor_slice = FALSE)
	disabled = suppressMessages(do.call(lame, args));expect_null(disabled$actor_factor_slice)
	args$prior$actor_factor_slice = 1
	expect_error(do.call(lame, args), 'actor_factor_slice must be TRUE or FALSE', fixed = TRUE)
	pooled = .lame_pool_sampler_diagnostics(list(), list(f, f))
	expect_null(pooled$actor_factor_slice)
	expect_equal(pooled$sampler_diagnostics_by_chain[[2]]$actor_factor_slice, health)
})
test_that('actor slice numerical failures abort the fit', {
	local_mocked_bindings(
		.lame_binary_actor_slice_round = function(...) stop('injected numerical failure'),
		.package = 'lame'
	)
	e = expect_error(suppressMessages(do.call(lame, .actor_slice_fit_args())), class = 'lame_sampler_failure')
	expect_match(gsub('[[:space:]]+', ' ', conditionMessage(e)), 'failed numerically on 5 of 5 rounds')
	expect_no_error(
		.lame_report_actor_slice_health(
			list(actor_factor_slice = c(rounds = 5, numerical_rejections = 0, bound_rejections = 20))
		)
	)
})
