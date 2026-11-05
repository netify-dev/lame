# rebuilds vignettes/articles/lame-overview-fits.rds, the stored converged fits
# that vignettes/articles/lame-overview.Rmd loads
#
# run from the package root with the package installed:
#   Rscript data-raw/lame-overview-fits.R [cores]
#
# every model runs two seeded chains and the script stops if a chain misses the
# R-hat and ESS thresholds below; a rerun on the same build and platform
# reproduces the file (compare the printed md5)
#
# chain length: at 60,000 post-burn-in iterations both seeds of the slowest model
# (R = 1, estimated reciprocity) meet the thresholds with room to spare; the R = 1
# model with independent residuals mixes more slowly, so both comparison models run 100,000
#
# stored, with the draw arrays at single precision:
#   fit          R = 1 fit, seed 1, without the pieces the vignette does not read
#   chain2       coefficient and variance draws of the same model, seed 2
#   fit_noUV     R = 0 fit, seed 1, slimmed the same way
#   loo          loo objects of the two independent-residual models (seed 1)
#   diagnostics  split R-hat and ESS of every chain and across the two chains

suppressPackageStartupMessages({
	library(lame)
	library(posterior)
	library(loo)
	library(parallel)
})
Sys.setenv(OMP_NUM_THREADS = "1")

args <- commandArgs(trailingOnly = TRUE)
cores <- if (length(args) >= 1L) as.integer(args[1L]) else 4L
out_file <- file.path("vignettes", "articles", "lame-overview-fits.rds")
stopifnot(dir.exists("vignettes"), cores >= 1L)

rhat_max <- 1.01
ess_min <- 400
seeds <- c(1L, 2L)
models <- c("fit", "fit_noUV", "loo_R0", "loo_R1")

# data, as built in the setup chunk of vignettes/articles/lame-overview.Rmd
data("dutchcollege", package = "lame")
n <- nrow(dutchcollege$Y)
T_all <- dim(dutchcollege$Y)[3]
actor_names <- sprintf("S%02d", seq_len(n))
Y <- lapply(2:T_all, function(t) {
	Yt <- (dutchcollege$Y[, , t] > 0) * 1
	diag(Yt) <- NA
	rownames(Yt) <- colnames(Yt) <- actor_names
	Yt
})
names(Y) <- paste0("t", 2:T_all)
X_node <- dutchcollege$X
rownames(X_node) <- actor_names
Xrow <- lapply(seq_along(Y), function(t) X_node)
Xcol <- Xrow
same_male <- outer(X_node[, "male"], X_node[, "male"], "==") * 1
same_smoker <- outer(X_node[, "smoker"], X_node[, "smoker"], "==") * 1
same_program <- outer(X_node[, "program"], X_node[, "program"], "==") * 1
Xdyad_one <- array(
	0,
	dim = c(n, n, 3),
	dimnames = list(
		actor_names,
		actor_names,
		c("same_male", "same_smoker", "same_program")
	)
)
Xdyad_one[, , 1] <- same_male
Xdyad_one[, , 2] <- same_smoker
Xdyad_one[, , 3] <- same_program
Xdyad <- lapply(seq_along(Y), function(t) Xdyad_one)

# the four models, with the arguments shown in the vignette
fit_model <- function(model, seed) {
	switch(
		model,
		fit = lame(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			family = "binary",
			rvar = TRUE,
			cvar = TRUE,
			dcor = TRUE,
			posterior_opts = posterior_options(save_ab = TRUE, save_UV = TRUE, thin_ab = 5),
			R = 1,
			symmetric = FALSE,
			burn = 2000,
			nscan = 60000,
			odens = 60,
			seed = seed,
			chains = 1,
			verbose = FALSE,
			plot = FALSE
		),
		fit_noUV = lame(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			family = "binary",
			rvar = TRUE,
			cvar = TRUE,
			dcor = TRUE,
			posterior_opts = posterior_options(save_ab = FALSE, save_UV = FALSE),
			R = 0,
			symmetric = FALSE,
			burn = 2000,
			nscan = 60000,
			odens = 60,
			seed = seed,
			chains = 1,
			gof = FALSE,
			verbose = FALSE,
			plot = FALSE
		),
		loo_R0 = ,
		loo_R1 = {
			loo_args <- list(
				Y = Y,
				Xdyad = Xdyad,
				Xrow = Xrow,
				Xcol = Xcol,
				family = "binary",
				rvar = TRUE,
				cvar = TRUE,
				dcor = FALSE,
				symmetric = FALSE,
				burn = 2000,
				nscan = 100000,
				odens = 100,
				seed = seed,
				chains = 1,
				save_log_lik = TRUE,
				verbose = FALSE,
				plot = FALSE
			)
			do.call(lame, c(loo_args, list(R = if (model == "loo_R0") 0 else 1)))
		}
	)
}

# large draw arrays at single precision (7 significant digits), far below the Monte Carlo error
single <- function(x) {
	y <- readBin(writeBin(as.numeric(x), raw(), size = 4L), "numeric", n = length(x), size = 4L)
	y[is.na(x)] <- x[is.na(x)]
	attributes(y) <- attributes(x)
	y
}

# drop what the vignette does not read; keep latent draws only for the iterations
# whose additive-effect draws are stored
slim_fit <- function(fit) {
	for (nm in c(
		"forecast_state", "start_vals", "refit_data", "refit_options", "log_lik", "log_lik_index"
	)) fit[[nm]] <- NULL
	if (!is.null(fit$U_samples)) {
		keep <- fit$ab_draw_index
		fit$U_samples <- fit$U_samples[, , keep, drop = FALSE]
		fit$V_samples <- fit$V_samples[, , keep, drop = FALSE]
		fit$uv_draw_index <- keep
	}
	for (nm in c("BETA", "VC", "a_samples", "b_samples", "U_samples", "V_samples")) {
		if (!is.null(fit[[nm]])) fit[[nm]] <- single(fit[[nm]])
	}
	fit$GOF <- lapply(fit$GOF, single)
	fit
}

# optional directory caching each raw chain, so the assembly can rerun without refitting
cache <- Sys.getenv("LAME_FITS_CACHE", unset = "")

run_job <- function(model, seed) {
	f <- if (nzchar(cache)) file.path(cache, sprintf("%s_s%d.rds", model, seed))
	if (!is.null(f) && file.exists(f)) return(readRDS(f))
	t0 <- Sys.time()
	fit <- fit_model(model, seed)
	draws <- cbind(fit$BETA, fit$VC)
	if (!is.null(fit$call$seed)) fit$call$seed <- seed
	raw <- if (seed == seeds[1L]) {
		if (model %in% c("fit", "fit_noUV")) fit else loo::loo(fit)
	}
	cat(sprintf("%-9s seed %d: %.0f s\n", model, seed, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
	r <- list(model = model, seed = seed, draws = draws, raw = raw)
	if (!is.null(f)) saveRDS(r, f)
	r
}

jobs <- expand.grid(model = models, seed = seeds, stringsAsFactors = FALSE)
res <- mclapply(
	seq_len(nrow(jobs)),
	function(i) run_job(jobs$model[i], jobs$seed[i]),
	mc.cores = min(cores, nrow(jobs)),
	mc.preschedule = FALSE
)
bad <- vapply(res, function(r) inherits(r, "try-error") || !is.list(r), logical(1))
if (any(bad)) stop("failed jobs: ", paste(jobs$model[bad], jobs$seed[bad], collapse = ", "))
pick <- function(model, seed) res[[which(jobs$model == model & jobs$seed == seed)]]
stored <- function(model) {
	raw <- pick(model, seeds[1L])$raw
	if (model %in% c("fit", "fit_noUV")) slim_fit(raw) else raw
}

# R-hat and ESS per chain and across the two chains, constant columns (ve) dropped
chain_table <- function(d) {
	s <- as.data.frame(posterior::summarise_draws(d, "rhat", "ess_bulk", "ess_tail"))
	s[!is.na(s$rhat), ]
}
as_chain <- function(m) posterior::as_draws_array(posterior::as_draws_matrix(m))
diagnostics <- do.call(rbind, lapply(models, function(m) {
	d1 <- as_chain(single(pick(m, seeds[1L])$draws))
	d2 <- as_chain(single(pick(m, seeds[2L])$draws))
	rbind(
		cbind(model = m, chain = "seed 1", chain_table(d1)),
		cbind(model = m, chain = "seed 2", chain_table(d2)),
		cbind(model = m, chain = "both", chain_table(posterior::bind_draws(d1, d2, along = "chain")))
	)
}))
rownames(diagnostics) <- NULL
passed <- with(
	diagnostics,
	tapply(
		rhat < rhat_max & ess_bulk >= ess_min & ess_tail >= ess_min,
		list(model, chain),
		all
	)
)
print(passed)
summary_by <- do.call(
	rbind,
	lapply(
		split(diagnostics, diagnostics[c("model", "chain")], drop = TRUE),
		function(d)
			data.frame(
				model = d$model[1L],
				chain = d$chain[1L],
				max_rhat = max(d$rhat),
				min_ess_bulk = min(d$ess_bulk),
				min_ess_tail = min(d$ess_tail)
			)
	)
)
print(summary_by, row.names = FALSE)
if (!all(passed)) stop("some chain misses the thresholds; lengthen the chains")

# the vignette's own table must agree with the one computed here
vt <- chain_table(posterior::as_draws(stored("fit")))
st <- subset(diagnostics, model == "fit" & chain == "seed 1")
stopifnot(isTRUE(all.equal(vt$rhat, st$rhat)), isTRUE(all.equal(vt$ess_bulk, st$ess_bulk)))

out <- list(
	meta = list(
		lame_version = as.character(utils::packageVersion("lame")),
		burn = 2000L,
		nscan = c(fit = 60000L, fit_noUV = 60000L, loo_R0 = 100000L, loo_R1 = 100000L),
		odens = c(fit = 60L, fit_noUV = 60L, loo_R0 = 100L, loo_R1 = 100L),
		stored = 1000L,
		thin_ab = 5L,
		seeds = seeds,
		rhat_max = rhat_max,
		ess_min = ess_min
	),
	fit = stored("fit"),
	chain2 = single(pick("fit", seeds[2L])$draws),
	fit_noUV = stored("fit_noUV"),
	loo = list(no_latent = stored("loo_R0"), latent_R1 = stored("loo_R1")),
	diagnostics = diagnostics
)
saveRDS(out, out_file, compress = "xz")
cat(sprintf("\n%s: %.0f KB, md5 %s\n", out_file, file.size(out_file) / 1024, unname(tools::md5sum(out_file))))
