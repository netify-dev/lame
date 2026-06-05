#!/usr/bin/env Rscript
# byte-identical-default verifier. each fit is run with default arguments
# so any change to the base mcmc kernel flips the hash. results are
# written to inst/verify_default_hash_results.csv; exit code signals pass
# (0) or fail (1). a one-line stderr summary is emitted for ci pipes.
#
# usage:
#   rscript inst/verify_default_hash.r baseline
#   rscript inst/verify_default_hash.r verify

suppressMessages({
	library(lame)
	if (!requireNamespace("digest", quietly = TRUE)) {
		stop("the 'digest' package is required for hash verification")
	}
})

make_uni_normal = function(seed) {
	set.seed(seed)
	n = 12L; T = 3L
	Y = array(stats::rnorm(n * n * T), c(n, n, T))
	X = array(stats::rnorm(n * n * T * 2L), c(n, n, T, 2L))
	for (t in seq_len(T)) diag(Y[, , t]) = NA
	list(Y = Y, X = X, family = "normal", mode = "unipartite")
}
make_uni_binary = function(seed) {
	set.seed(seed)
	n = 12L; T = 3L
	Y = array(rbinom(n * n * T, 1, 0.3), c(n, n, T))
	X = array(stats::rnorm(n * n * T * 2L), c(n, n, T, 2L))
	for (t in seq_len(T)) diag(Y[, , t]) = NA
	list(Y = Y, X = X, family = "binary", mode = "unipartite")
}
make_bip_normal = function(seed) {
	set.seed(seed)
	nA = 8L; nB = 6L; T = 3L
	Y = array(stats::rnorm(nA * nB * T), c(nA, nB, T))
	X = array(stats::rnorm(nA * nB * T * 2L), c(nA, nB, T, 2L))
	list(Y = Y, X = X, family = "normal", mode = "bipartite")
}

configs = list(
	uni_normal = make_uni_normal,
	uni_binary = make_uni_binary,
	bip_normal = make_bip_normal
)
seeds = c(1L, 2L, 3L, 4L, 5L)
hash_fields = c("BETA", "VC", "APM", "BPM", "U", "V", "UVPM", "EZ", "YPM")

compute_hash = function(cfg_name, seed) {
	dat = configs[[cfg_name]](seed)
	Y_list = lapply(seq_len(dim(dat$Y)[3L]), function(t) dat$Y[, , t])
	X_list = lapply(seq_len(dim(dat$X)[3L]), function(t) dat$X[, , t, ])
	args = list(
		Y = Y_list, Xdyad = X_list,
		family = dat$family,
		mode = dat$mode,
		R = 0L,
		nscan = 40L, burn = 10L, odens = 5L,
		seed = seed,
		verbose = FALSE, plot = FALSE
	)
	fit = suppressWarnings(suppressMessages(do.call(lame::lame, args)))
	picked = fit[hash_fields]
	digest::digest(picked, algo = "sha256")
}

build_matrix = function() {
	grid = expand.grid(cfg = names(configs), seed = seeds,
	                   stringsAsFactors = FALSE)
	grid$hash = NA_character_
	for (i in seq_len(nrow(grid))) {
		grid$hash[i] = tryCatch(
			compute_hash(grid$cfg[i], grid$seed[i]),
			error = function(e) paste0("ERR:", conditionMessage(e)))
	}
	grid
}

main = function(mode) {
	pkg_root = tryCatch(rprojroot::find_root(rprojroot::is_r_package),
	                    error = function(e) ".")
	baseline_path = file.path(pkg_root, "inst", "verify_default_hash_baseline.rds")
	results_path = file.path(pkg_root, "inst", "verify_default_hash_results.csv")
	current = build_matrix()
	utils::write.csv(current, results_path, row.names = FALSE)

	if (mode == "baseline") {
		saveRDS(current, baseline_path)
		writeLines("BASELINE_WRITTEN", con = stderr())
		return(invisible(current))
	}

	if (mode == "verify") {
		if (!file.exists(baseline_path)) {
			stop("baseline not found at ", baseline_path,
			     " - run with `baseline` first")
		}
		baseline = readRDS(baseline_path)
		merged = merge(baseline, current,
		               by = c("cfg", "seed"),
		               suffixes = c("_base", "_now"))
		mismatch = merged[merged$hash_base != merged$hash_now, ]
		mismatch_path = file.path(pkg_root, "inst", "verify_default_hash_mismatches.csv")
		utils::write.csv(mismatch, mismatch_path, row.names = FALSE)
		if (nrow(mismatch) == 0L) {
			writeLines(sprintf("ALL_HASHES_MATCH (%d configs)", nrow(merged)),
			           con = stderr())
			return(invisible(current))
		}
		writeLines(sprintf("HASH_MISMATCH (%d configs)", nrow(mismatch)),
		           con = stderr())
		quit(status = 1L)
	}
	stop("usage: Rscript inst/verify_default_hash.R [baseline|verify]")
}

args = commandArgs(trailingOnly = TRUE)
mode = if (length(args) >= 1L) args[[1L]] else "verify"
main(mode)
