# independently seeded chains per fit (seed + 1000 * (k - 1)) and the between-chain gate:
# every chain is extended from its final state while one drifts or a monitored split
# R-hat exceeds its threshold, up to prior$chain_extensions extensions

# smallest split R-hat the gate flags; a series is flagged only above its noise bound as well
.lame_chain_rhat_max <- 1.05
# family-wise false-alarm level per fit of each check (drift, R-hat) on settled chains
.lame_chain_alpha <- 0.01
# chains are compared on a series only when its threshold flags chains this many posterior sds
# apart with at least this probability
.lame_chain_detect_sd <- 1
.lame_chain_detect_power <- 0.95
# an extension for series short of effective draws runs long enough to hold a quarter more than they
# need, up to this many times nscan
.lame_chain_ess_margin <- 1.25
.lame_chain_nscan_max <- 4
# stored draws per chain the gate needs (as for the drift check)
.lame_chain_gate_min_draws <- 60L

# number of chains when the user does not say: options(lame.chains), else 2
.lame_chains_default <- function() getOption("lame.chains", 2L)

#' @noRd
.lame_chains_check <- function(chains) {
	if (!is.numeric(chains) || length(chains) != 1L || !is.finite(chains) || chains < 1 || chains != floor(chains)) {
		cli::cli_abort("{.arg chains} must be a whole number of at least 1.", class = "lame_chains_arg")
	}
	as.integer(chains)
}

# the path of chain k when a single path was given: file_chain<k>.ext
.lame_chain_path <- function(path, k) {
	if (is.null(path)) return(NULL)
	if (grepl("\\.[^./\\\\]+$", path)) sub(
		"(\\.[^./\\\\]+)$",
		paste0("_chain", k, "\\1"),
		path
	) else paste0(path, "_chain", k)
}

# the operating system, as a function so tests can stand in for Windows
.lame_chains_os <- function() .Platform$OS.type

# the number of cores for the chains: the argument, else options(lame.cores), else all but one
.lame_chains_cores <- function(cores, chains) {
	if (is.null(cores)) cores <- getOption("lame.cores", NULL)
	if (is.null(cores)) cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
	if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1) {
		cli::cli_abort("{.arg cores} must be a whole number of at least 1.", class = "lame_chains_arg")
	}
	if (identical(.lame_chains_os(), "windows")) return(1L)
	as.integer(min(chains, cores))
}

# run `fn(j)` for j in seq_len(n) on `cores` forked workers (sequentially for one)
.lame_chains_map <- function(n, fn, cores) {
	if (cores > 1L && n > 1L) {
		out <- parallel::mclapply(seq_len(n), fn, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
		bad <- vapply(out, inherits, logical(1), "try-error")
		if (any(bad)) cli::cli_abort(
			c("A chain worker stopped.", "i" = "{as.character(out[[which(bad)[1L]]])}"),
			class = "lame_chains_worker"
		)
		out
	} else lapply(seq_len(n), fn)
}

# one chain call; its warnings and messages are kept for the parent to re-signal,
# and a chain that is still moving comes back as its condition
.lame_chain_worker <- function(fun, args) {
	conds <- list()
	res <- withCallingHandlers(
		tryCatch(do.call(fun, args), lame_chain_unsettled = function(e) e, error = function(e) e),
		warning = function(w) { conds[[length(conds) + 1L]] <<- w; invokeRestart("muffleWarning") },
		message = function(m) { conds[[length(conds) + 1L]] <<- m; invokeRestart("muffleMessage") }
	)
	list(res = res, conds = conds)
}

# re-signal the warnings and messages of the chains, each distinct one once
.lame_chains_resignal <- function(rs, seen) {
	for (r in rs) for (cnd in r$conds) {
		key <- paste(class(cnd)[1L], conditionMessage(cnd))
		if (key %in% seen$keys) next
		seen$keys <- c(seen$keys, key)
		if (inherits(cnd, "warning")) warning(cnd) else message(cnd)
	}
	invisible(seen)
}

# the monitored series of one chain, a matrix of stored draws by series
.lame_gate_series <- function(fit) {
	reason <- .lame_ab_var_unident_reason(fit)
	hidden <- character(0)
	if (!is.na(reason)) hidden <- .lame_ab_var_hidden_rows(reason)
	coef_reason <- .lame_ab_coef_reason(fit)
	cols <- list()
	B <- fit$BETA
	if (!is.null(B) && length(B) > 0L) {
		nm <- dimnames(B)[[2L]] %||% paste0("beta", seq_len(dim(B)[2L]))
		keep <- which(!nm %in% .lame_ab_coef_hidden_rows(coef_reason, nm))
		if (length(dim(B)) == 3L) {
			per <- dimnames(B)[[3L]] %||% paste0("p", seq_len(dim(B)[3L]))
			# a static coefficient of a fit with dynamic blocks is the same in every period: monitor it once
			for (j in keep) {
				if (length(fit$beta_dynamic_mask) == dim(B)[2L] && !isTRUE(fit$beta_dynamic_mask[j])) {
					cols[[nm[j]]] <- B[, j, 1L]
				} else {
					for (t in seq_len(dim(B)[3L])) cols[[paste0(nm[j], "[", per[t], "]")]] <- B[, j, t]
				}
			}
		} else for (j in keep) cols[[nm[j]]] <- B[, j]
	}
	if (!is.null(fit$VC) && !is.null(colnames(fit$VC))) {
		for (v in setdiff(
			intersect(colnames(fit$VC), c("va", "vb", "cab", "ve", "rho")),
			hidden
		)) cols[[v]] <- fit$VC[, v]
	}
	for (v in setdiff(c("rho_ab", "sigma_ab", "rho_uv", "sigma_uv", "pi_snap"), hidden)) {
		if (is.numeric(fit[[v]]) && length(fit[[v]]) > 1L) cols[[v]] <- as.numeric(fit[[v]])
	}
	for (v in c("RHO_BETA", "SIGMA_BETA", "MU_BETA", "RHO")) {
		# skip blocks whose persistence is 1 by construction (rw1, rw2, matern32)
		M <- if (identical(v, "RHO_BETA")) .lame_rho_beta_ar1(fit) else fit[[v]]
		if (is.numeric(M) && length(M) > 1L) {
			M <- as.matrix(M)
			nmv <- colnames(M) %||% if (ncol(M) > 1L) paste0("g", seq_len(ncol(M))) else ""
			for (j in seq_len(ncol(M))) cols[[if (nzchar(nmv[j])) paste0(tolower(v), "[", nmv[j], "]") else tolower(v)]] <- M[, j]
		}
	}
	S <- fit$chain_summary
	if (is.matrix(S)) for (j in seq_len(ncol(S))) cols[[colnames(S)[j]]] <- S[, j]
	if (length(cols) == 0L) return(NULL)
	n <- min(lengths(cols))
	out <- do.call(cbind, lapply(cols, function(x) as.numeric(x)[seq_len(n)]))
	out
}

# split rank-normalised R-hat of every monitored series across chains with the gate's verdict on it
# (see .lame_gate_verdict); empty when the runs are too short
.lame_gate_table <- function(fits, alpha = .lame_chain_alpha) {
	empty <- data.frame(
		parameter = character(0),
		rhat = numeric(0),
		ess = numeric(0),
		ess_chain = numeric(0),
		threshold = numeric(0),
		informative = logical(0),
		ess_needed = numeric(0),
		stringsAsFactors = FALSE
	)
	if (length(fits) < 2L) return(empty)
	ser <- lapply(fits, .lame_gate_series)
	if (any(vapply(ser, is.null, logical(1)))) return(empty)
	n <- min(vapply(ser, nrow, integer(1)))
	if (n < .lame_chain_gate_min_draws) return(empty)
	common <- Reduce(intersect, lapply(ser, colnames))
	if (length(common) == 0L) return(empty)
	parts <- vapply(
		common,
		function(v) .lame_rhat_parts(vapply(ser, function(m) m[seq_len(n), v], numeric(n))),
		numeric(5)
	)
	cbind(data.frame(parameter = common, stringsAsFactors = FALSE), .lame_gate_verdict(parts, alpha))
}

# the gate's verdict on the series of `parts` (columns from .lame_rhat_parts). The threshold is 1.05
# or the R-hat settled chains exceed with probability a, whichever is larger, with a the Sidak level
# of the bulk and tail R-hat of every compared series at the family-wise level alpha. The chains can
# be compared on a series when that threshold flags chains .lame_chain_detect_sd posterior sds apart
# with probability .lame_chain_detect_power, which takes ess_needed effective draws per split chain
# (pooled over the chains). A series passes when it can be compared and its R-hat is at most its
# threshold; chains constant at one value pass and chains constant at different values (R-hat Inf)
# fail
.lame_gate_verdict <- function(parts, alpha = .lame_chain_alpha) {
	rhat <- unname(parts["rhat", ])
	ok <- is.finite(rhat)
	a <- 1 - (1 - alpha)^(1 / max(1, 2 * sum(ok)))
	n <- unname(parts["n", ok])
	k <- unname(parts["chains", ok])
	e <- unname(parts["ess", ok])
	out <- data.frame(
		rhat = rhat,
		ess = unname(parts["ess", ]),
		ess_chain = unname(parts["ess_chain", ]),
		threshold = NA_real_,
		informative = ifelse(is.infinite(rhat), TRUE, NA),
		ess_needed = NA_real_
	)
	if (!any(ok)) return(out)
	out$threshold[ok] <- pmax(.lame_chain_rhat_max, .lame_rhat_bound(a, n, k, e))
	out$ess_needed[ok] <- .lame_rhat_ess_needed(a, n, k)
	out$informative[ok] <- e >= out$ess_needed[ok]
	out
}

# the R-hat every reported surface uses, so summary(), compute_mcmc_diagnostics() and
# rhat_dynamic_beta() match fit$rhat; NA under 8 draws per chain, for non-finite series and for
# chains constant at one value, Inf for chains constant at different values
.lame_rhat_chains <- function(x) {
	x <- as.matrix(x)
	if (nrow(x) < 8L || ncol(x) < 2L || !all(is.finite(x))) return(NA_real_)
	tol <- 1e-10 * max(1, max(abs(x)))
	if (all(apply(x, 2L, stats::sd) < tol)) return(if (diff(range(x)) > tol) Inf else NA_real_)
	.lame_rhat(x)
}
# split rank-normalised R-hat (Vehtari et al. 2021) of a draws-by-chains matrix: the larger
# of the bulk and folded (tail) values, as posterior::rhat() reports it
.lame_rhat <- function(x) {
	x <- as.matrix(x)
	half <- nrow(x) %/% 2L
	if (half < 4L) return(NA_real_)
	split <- cbind(x[seq_len(half), , drop = FALSE], x[nrow(x) - half + seq_len(half), , drop = FALSE])
	basic <- function(m) {
		n <- nrow(m)
		w <- mean(apply(m, 2L, stats::var))
		b <- n * stats::var(colMeans(m))
		if (!is.finite(w) || w <= 0) return(NA_real_)
		sqrt(((n - 1) / n * w + b / n) / w)
	}
	bulk <- basic(.lame_rank_z(split))
	tail <- basic(.lame_rank_z(abs(split - stats::median(split))))
	max(bulk, tail, na.rm = TRUE)
}

# normal scores of the pooled ranks of a matrix, keeping its shape
.lame_rank_z <- function(m) {
	r <- rank(m, ties.method = "average")
	matrix(stats::qnorm((r - 3 / 8) / (length(m) + 1 / 4)), nrow(m), ncol(m))
}

# R-hat of a draws-by-chains matrix with the effective draws per split chain (pooled over the chains,
# and of the chain with the fewest), the split length and the split count; the effective draws come
# from the autocorrelation of each whole chain, rank normalised, the smaller of the bulk and folded
.lame_rhat_parts <- function(x) {
	out <- c(rhat = NA_real_, ess = NA_real_, ess_chain = NA_real_, n = NA_real_, chains = NA_real_)
	x <- as.matrix(x)
	out["rhat"] <- .lame_rhat_chains(x)
	if (!is.finite(out["rhat"])) return(out)
	zb <- .lame_rank_z(x)
	zt <- .lame_rank_z(abs(x - stats::median(x)))
	tau <- function(cols) {
		t <- suppressWarnings(max(
			.lame_within_tau(zb[, cols, drop = FALSE]),
			.lame_within_tau(zt[, cols, drop = FALSE]),
			na.rm = TRUE
		))
		if (is.finite(t)) t else 1
	}
	half <- nrow(x) %/% 2L
	out["ess"] <- min(half, half / tau(seq_len(ncol(x))))
	out["ess_chain"] <- min(half, half / vapply(seq_len(ncol(x)), tau, numeric(1)))
	out["n"] <- half
	out["chains"] <- 2 * ncol(x)
	out
}

# integrated autocorrelation time of a draws-by-chains matrix from the autocovariances within
# chains (each about its own mean, averaged over chains), by Geyer's initial monotone sequence;
# at least 1
.lame_within_tau <- function(m) {
	n <- nrow(m)
	mc <- sweep(m, 2L, colMeans(m))
	acov <- function(t) {
		early <- mc[seq_len(n - t), , drop = FALSE]
		late <- mc[t + seq_len(n - t), , drop = FALSE]
		mean(colSums(early * late)) / n
	}
	v0 <- acov(0L)
	if (!is.finite(v0) || v0 <= 0) return(NA_real_)
	tau <- -1
	prev <- Inf
	for (k in seq(0L, n - 2L, by = 2L)) {
		p <- (acov(k) + acov(k + 1L)) / v0
		if (p <= 0) break
		p <- min(p, prev)
		tau <- tau + 2 * p
		prev <- p
	}
	max(tau, 1)
}

# split R-hat that settled chains exceed with probability a: for k split chains of n draws with e
# effective draws each, R-hat^2 = (n - 1) / n * (1 + F / (e - 1)) with F ~ F(k - 1, k (e - 1))
.lame_rhat_bound <- function(a, n, k, e) {
	e <- pmin(pmax(e, 2), n)
	sqrt((n - 1) / n * (1 + stats::qf(1 - a, k - 1, k * (e - 1)) / (e - 1)))
}

# probability that `threshold` flags chains delta posterior sds apart (two chains delta apart, any
# others between them): their squared split R-hat is about (n - 1) / n * (1 + F / (e - 1)), F
# noncentral F(k - 1, k (e - 1)) with noncentrality e * delta^2
.lame_rhat_power <- function(threshold, n, k, e, delta = .lame_chain_detect_sd) {
	e <- pmin(pmax(e, 2), n)
	f <- (e - 1) * (threshold^2 * n / (n - 1) - 1)
	stats::pf(f, k - 1, k * (e - 1), ncp = e * delta^2, lower.tail = FALSE)
}

# effective draws per split chain from which the threshold of a test at level a on k split chains
# of n draws flags chains .lame_chain_detect_sd posterior sds apart with probability
# .lame_chain_detect_power (the power rises with e)
.lame_rhat_ess_needed <- function(a, n, k) {
	gap <- function(e, n, k) {
		m <- max(n, e)
		thr <- max(.lame_chain_rhat_max, .lame_rhat_bound(a, m, k, e))
		.lame_rhat_power(thr, m, k, e) - .lame_chain_detect_power
	}
	key <- paste(n, k)
	u <- !duplicated(key)
	root <- mapply(function(n, k) stats::uniroot(gap, c(2, 1e5), n = n, k = k, tol = 1e-4)$root, n[u], k[u])
	unname(root[match(key, key[u])])
}

# the series the gate flags: R-hat above the threshold (1.05 for a table without one), chains
# constant at different values, and series with too few effective draws to compare the chains
.lame_gate_flagged <- function(tab, rhat_max = .lame_chain_rhat_max) {
	thr <- if (is.null(tab$threshold)) rep(rhat_max, nrow(tab)) else pmax(rhat_max, tab$threshold)
	thr[!is.finite(thr)] <- rhat_max
	few <- if (is.null(tab$informative)) logical(nrow(tab)) else tab$informative %in% FALSE
	tab[!is.na(tab$rhat) & (tab$rhat > thr | few), , drop = FALSE]
}

# flagged series (rows of .lame_gate_flagged) that are flagged only for too few effective draws
.lame_gate_short <- function(bad, rhat_max = .lame_chain_rhat_max) {
	thr <- if (is.null(bad$threshold)) rep(rhat_max, nrow(bad)) else pmax(rhat_max, bad$threshold)
	thr[!is.finite(thr)] <- rhat_max
	is.finite(bad$rhat) & bad$rhat <= thr
}

# why each flagged series (rows of .lame_gate_flagged) was flagged
.lame_gate_text <- function(bad, rhat_max = .lame_chain_rhat_max) {
	if (nrow(bad) == 0L) return(character(0))
	thr <- if (is.null(bad$threshold)) rep(rhat_max, nrow(bad)) else pmax(rhat_max, bad$threshold)
	thr[!is.finite(thr)] <- rhat_max
	out <- sprintf("%s (R-hat %.3f above %.3f)", bad$parameter, bad$rhat, thr)
	few <- .lame_gate_short(bad, rhat_max)
	if (any(few)) out[few] <- sprintf(
		"%s (%.0f effective draws per split chain, %.0f needed to compare the chains)",
		bad$parameter[few],
		bad$ess[few],
		ceiling(bad$ess_needed[few])
	)
	flat <- is.infinite(bad$rhat)
	out[flat] <- sprintf("%s (constant at different values in different chains)", bad$parameter[flat])
	out
}

# nscan at which every flagged series short of effective draws would hold a quarter more than it
# needs (at the same autocorrelation), rounded up to two significant digits and then to a multiple of
# odens; nscan when none is short
.lame_gate_nscan <- function(bad, nscan, odens = 1L) {
	few <- .lame_gate_short(bad) & bad$informative %in% FALSE & is.finite(bad$ess) & bad$ess > 0
	if (!any(few)) return(nscan)
	x <- nscan * max(1, .lame_chain_ess_margin * max(bad$ess_needed[few] / bad$ess[few]))
	step <- 10^max(0, floor(log10(x)) - 1)
	ceiling(ceiling(x / step) * step / odens) * odens
}

# the state a chain continues from: the exact checkpoint of a tempered chain, else its start values
.lame_chain_state <- function(fit) {
	if (!is.null(fit$tempered_checkpoint)) fit["tempered_checkpoint"] else fit$start_vals
}

# run `chains` chains of `fun` ("ame" or "lame") on the evaluated `args`; returns the pooled
# fit, or the chain list when `as_list`; `start_list` holds per-chain start values
.lame_chains_run <- function(
	fun,
	args,
	chains,
	seeds,
	call,
	cores = NULL,
	caller = parent.frame(),
	as_list = FALSE,
	start_list = NULL,
	verbose = FALSE,
	plot = FALSE,
	resolved = NULL,
	input_prior = NULL
) {
	fn <- get(fun, envir = asNamespace("lame"), mode = "function")
	prior <- args$prior %||% list()
	K <- as.integer(.lame_chain_limit(prior))
	nscan <- args$nscan %||% eval(formals(fn)$nscan)
	burn <- args$burn %||% eval(formals(fn)$burn)
	odens <- args$odens %||% eval(formals(fn)$odens)
	cores <- .lame_chains_cores(cores, chains)
	seen <- new.env(); seen$keys <- character(0)
	if (is.null(start_list)) start_list <- vector("list", chains)
	if (!is.null(start_list[[1L]])) burn <- 0L
	if (is.null(input_prior)) input_prior <- prior
	# chunked likelihood storage: one directory per chain, rewritten by each extension
	ll_base <- if (identical(fun, "lame") && identical(args$save_log_lik, "chunked")) {
		sub("[/\\\\]+$", "", args$log_lik_path %||% tempfile(pattern = "lame_loglik_"))
	} else NULL
	# the arguments of chain j in round r (0: the first run)
	chain_args <- function(j, r, prev) {
		a <- args
		a$chains <- 1L
		a$cores <- NULL
		a$seed <- seeds[j] + (if (r > 0L) burn + offset else 0)
		if (r > 0L) a$nscan <- nscan_r
		a$verbose <- FALSE
		a$plot <- FALSE
		for (nm in c("out_file", "checkpoint_path")) if (!is.null(a[[nm]])) a[[nm]] <- .lame_chain_path(a[[nm]], j)
		if (!is.null(ll_base)) a$log_lik_path <- paste0(ll_base, "_chain", j)
		cont <- r > 0L || !is.null(start_list[[j]])
		use_resolved <- cont && identical(fun, "lame") && !is.null(resolved)
		pr <- if (use_resolved) resolved else prior
		attr(pr, "lame_chain_block") <- list(mc = call, k = r, record = TRUE, nopolicy = TRUE)
		if (use_resolved) attr(pr, "lame_checkpoint_input_prior") <- input_prior
		a$prior <- pr
		if (r > 0L) { a$start_vals <- prev[[j]]; a$burn <- 0L }
		else if (!is.null(start_list[[j]])) { a$start_vals <- start_list[[j]]; a$burn <- 0L }
		a
	}
	run_round <- function(r, prev) {
		rs <- .lame_chains_map(chains, function(j) .lame_chain_worker(fn, chain_args(j, r, prev)), cores)
		.lame_chains_resignal(rs, seen)
		res <- lapply(rs, `[[`, "res")
		errs <- Filter(function(x) inherits(x, "error") && !inherits(x, "lame_chain_unsettled"), res)
		if (length(errs) > 0L) stop(errs[[1L]])
		res
	}
	say <- function(r, res) {
		if (!isTRUE(verbose)) return(invisible())
		for (j in seq_len(chains)) cli::cli_inform(
			c(
				"i" = "Chain {j} of {chains} (seed {seeds[j] + (if (r > 0L) burn + offset else 0)}): {if (r == 0L) 'run' else paste0('extension ', r)} complete."
			)
		)
	}
	triggers <- list()
	r <- 0L
	# sweeps of the discarded blocks, and the length of the current one
	offset <- 0
	nscan_r <- nscan
	prev <- NULL
	res <- run_round(0L, NULL)
	say(0L, res)
	repeat {
		own <- Filter(function(x) inherits(x, "lame_chain_unsettled"), res)
		fits <- lapply(res, function(x) if (inherits(x, "lame_chain_unsettled")) x$fit else x)
		if (identical(fun, "lame") && is.null(resolved)) {
			resolved <- fits[[1L]]$prior
			if (!is.null(fits[[1L]]$refit_options$prior)) input_prior <- fits[[1L]]$refit_options$prior
		}
		if (identical(fun, "lame") && is.null(resolved)) stop("internal: no resolved prior")
		early <- any(vapply(fits, function(f) isTRUE(f$terminated_early), logical(1)))
		moving <- lapply(fits, function(f) if (early) NULL else .lame_chain_moving(f, n_chains = length(fits)))
		gate <- if (early) .lame_gate_table(list()) else .lame_gate_table(fits)
		bad <- .lame_gate_flagged(gate)
		drift_txt <- unlist(
			lapply(
				seq_along(moving),
				function(j) if (!is.null(moving[[j]]) && nrow(moving[[j]]) > 0L) paste0(
					"chain ",
					j,
					": ",
					moving[[j]]$text
				)
			)
		)
		rhat_txt <- .lame_gate_text(bad)
		short <- .lame_gate_short(bad)
		if (length(drift_txt) == 0L && length(rhat_txt) == 0L) break
		if (isTRUE(verbose)) cli::cli_inform(
			c("i" = "The chains are not settled ({paste(c(drift_txt, rhat_txt), collapse = ', ')}); round {r} of {K}.")
		)
		# chains that do not disagree but are still too short to compare are returned, not refused
		if (r >= K && length(drift_txt) == 0L && all(short)) {
			cli::cli_warn(
				c(
					"The chains hold too few effective draws to be compared on {rhat_txt}, after {r} extension{?s}.",
					"i" = "Continue them with {.code lame_resume(fit, nscan = {(.lame_gate_nscan(bad, nscan_r, odens))})}."
				),
				class = "lame_chain_not_compared"
			)
			break
		}
		if (r >= K) {
			pooled <- .lame_chains_pool(
				fits,
				seeds,
				burn,
				nscan_r,
				r,
				triggers,
				call,
				args,
				gate,
				input_prior,
				as_list = FALSE,
				refused = TRUE,
				discarded = offset
			)
			cli::cli_abort(
				c(
					if (r == 0L && K == 0L) "The chains were not settled when they stopped and {.code prior$chain_extensions} is 0, so they were not extended: {c(drift_txt, rhat_txt)}."
					else if (r == 0L) "The chains were not settled when they stopped: {c(drift_txt, rhat_txt)}."
					else "The chains were still not settled after {r} extension{?s} ({format(offset, big.mark = ',')} sweeps discarded): {c(drift_txt, rhat_txt)}.",
					"i" = "The pooled fit of the last block is in the {.field fit} field of the condition and the per-chain fits in {.field fits}; the offending series are in {.field rhat}. Continue every chain with {.code lame_resume(e$fit, nscan = {(.lame_gate_nscan(bad, nscan_r, odens))})}, or one chain with {.code lame_resume(e$fits[[k]], nscan = {(.lame_gate_nscan(bad, nscan_r, odens))})}, where {.code e} is the caught condition."
				),
				class = "lame_chain_unsettled",
				fit = pooled,
				fits = fits,
				rhat = gate,
				moving = c(drift_txt, rhat_txt),
				extensions = as.integer(r),
				call = caller
			)
		}
		triggers[[length(triggers) + 1L]] <- c(drift_txt, rhat_txt)
		offset <- offset + nscan_r
		# series short of effective draws get a block long enough to hold them
		if (any(short)) nscan_r <- min(max(nscan_r, .lame_gate_nscan(bad[short, , drop = FALSE], nscan_r, odens)), .lame_chain_nscan_max * nscan)
		prev <- lapply(fits, .lame_chain_state)
		rm(fits, res, own)
		r <- r + 1L
		res <- run_round(r, prev)
		say(r, res)
	}
	pooled <- .lame_chains_pool(
		fits,
		seeds,
		burn,
		nscan_r,
		r,
		triggers,
		call,
		args,
		gate,
		input_prior,
		as_list = as_list,
		refused = FALSE,
		discarded = offset
	)
	pooled
}

# the pooled fit (or the list of chain fits) with the chain records
.lame_chains_pool <- function(
	fits,
	seeds,
	burn,
	nscan,
	extensions,
	triggers,
	call,
	args,
	gate,
	input_prior,
	as_list = FALSE,
	refused = FALSE,
	discarded = extensions * nscan
) {
	chains <- length(fits)
	K <- as.integer(.lame_chain_limit(args$prior %||% list()))
	next_seeds <- seeds + burn + discarded + nscan
	rec <- .lame_chain_record(K, extensions, nscan, triggers, NA_real_)
	rec$discarded_sweeps <- as.integer(discarded)
	rec$block_nscan <- as.integer(nscan)
	rec$next_seed <- NULL
	rec$rhat_max <- .lame_chain_rhat_max
	rec$detect_sd <- .lame_chain_detect_sd
	rec$detect_power <- .lame_chain_detect_power
	rec$rhat <- gate
	rec$chains <- chains
	rec$next_seeds <- next_seeds
	if (refused) { rec$refused <- TRUE; rec$call_args <- args }
	for (j in seq_len(chains)) fits[[j]]$chain_id <- j
	if (as_list) {
		for (j in seq_len(chains)) fits[[j]]$auto_extend <- rec
		return(fits)
	}
	start_list <- lapply(fits, .lame_chain_state)
	empty <- all(vapply(fits, .lame_has_no_draws, logical(1)))
	# a run that stored no draws has nothing to diagnose
	diag <- if (empty) NULL else tryCatch(
		suppressMessages(suppressWarnings(compute_mcmc_diagnostics(fits))),
		error = function(e) {
			cli::cli_warn(
				c(
					"The pooled convergence diagnostics could not be computed, so {.code fit$diagnostics} is empty.",
					"x" = "{conditionMessage(e)}"
				),
				class = "lame_diagnostics_unavailable"
			)
			NULL
		}
	)
	fits_slim <- fits
	for (j in seq_len(chains)) fits_slim[[j]]$chain_summary <- NULL
	pooled <- if (empty) {
		# no chain stored a draw: the fit of the first chain, with the chain records
		f <- fits_slim[[1L]]
		f$n_chains <- chains
		f$is_combined <- TRUE
		f$chain_indicator <- integer(0)
		f
	} else suppressWarnings(suppressMessages(combine_ame_chains(fits_slim, diagnostics = FALSE)))
	if (identical(pooled[["sampler"]], "tempered")) pooled <- .lame_tempered_pool(pooled, fits)
	pooled$chain_summary <- NULL
	pooled$diagnostics <- diag
	pooled$call <- call
	pooled$chains <- chains
	pooled$chain_seeds <- seeds
	pooled$rhat <- gate
	pooled$chain_start_vals <- start_list
	pooled$auto_extend <- rec
	pooled$terminated_early <- any(vapply(fits, function(f) isTRUE(f$terminated_early), logical(1)))
	if (empty) pooled <- .lame_empty_draws(pooled)
	pooled
}

# the chains of a user's ame() call; `env` is the frame of ame() after its early
# argument handling
.lame_chains_from_ame <- function(env, mc, chains) {
	fm <- setdiff(names(formals(ame)), c("...", "chains", "cores"))
	args <- mget(intersect(fm, ls(env, all.names = TRUE)), envir = env)
	cores <- get("cores", envir = env, inherits = FALSE)
	.lame_chains_entry("ame", args, chains, mc, cores, parent.frame(2L))
}

# the chains of a user's lame() call: `ctx` holds the evaluated arguments the user
# supplied (after any conversion of the input)
.lame_chains_from_lame <- function(env, ctx, chains) {
	args <- ctx$args
	for (nm in c("seed", "burn", "nscan", "odens", "prior", "verbose", "plot")) {
		args[nm] <- list(get(nm, envir = env, inherits = FALSE))
	}
	args$chains <- NULL
	cores <- get("cores", envir = env, inherits = FALSE)
	args$cores <- NULL
	.lame_chains_entry("lame", args, chains, ctx$mc, cores, ctx$caller)
}

.lame_chains_entry <- function(fun, args, chains, call, cores, caller) {
	prior <- args$prior %||% list()
	as_list <- isTRUE(attr(prior, "lame_chain_list", exact = TRUE))
	if (as_list) { attr(prior, "lame_chain_list") <- NULL; args$prior <- prior }
	if (!is.numeric(args$seed) || length(args$seed) != 1L || !is.finite(args$seed)) args$seed <- 6886
	seeds <- as.numeric(args$seed) + 1000 * (seq_len(chains) - 1L)
	verbose <- isTRUE(args$verbose)
	plot <- isTRUE(args$plot)
	out <- .lame_chains_run(
		fun,
		args,
		chains,
		seeds,
		call,
		cores = cores,
		caller = caller,
		as_list = as_list,
		verbose = verbose,
		plot = plot
	)
	if (plot && !as_list) try(plot(out), silent = TRUE)
	out
}
