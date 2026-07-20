## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>",
	fig.width = 6,
	fig.height = 4,
	dpi = 96,
	out.width = "100%",
	fig.align = "center",
	dev = "png"
)

## ----setup--------------------------------------------------------------------
library(lame)
library(ggplot2)
set.seed(6886)

data("dutchcollege")

n <- nrow(dutchcollege$Y)
T_all <- dim(dutchcollege$Y)[3]
actor_names <- sprintf("S%02d", seq_len(n))

# binarize and drop the cold-start wave (t = 1 has almost no ties)
Y <- lapply(2:T_all, function(t) {
	Yt <- (dutchcollege$Y[, , t] > 0) * 1
	diag(Yt) <- NA
	rownames(Yt) <- colnames(Yt) <- actor_names
	Yt
})
names(Y) <- paste0("t", 2:T_all)

# nodal covariates passed identically as sender and receiver
X_node <- dutchcollege$X
rownames(X_node) <- actor_names
Xrow <- lapply(seq_along(Y), function(t) X_node)
Xcol <- Xrow

# dyadic homophily indicators
same_male    <- outer(X_node[, "male"],    X_node[, "male"],    "==") * 1
same_smoker  <- outer(X_node[, "smoker"],  X_node[, "smoker"],  "==") * 1
same_program <- outer(X_node[, "program"], X_node[, "program"], "==") * 1
Xdyad_one <- array(0, dim = c(n, n, 3),
                   dimnames = list(actor_names, actor_names,
                                   c("same_male", "same_smoker", "same_program")))
Xdyad_one[, , 1] <- same_male
Xdyad_one[, , 2] <- same_smoker
Xdyad_one[, , 3] <- same_program
Xdyad <- lapply(seq_along(Y), function(t) Xdyad_one)

## -----------------------------------------------------------------------------
sapply(Y, function(y) round(mean(y, na.rm = TRUE), 2))

## ----message=FALSE, warning=FALSE---------------------------------------------
fit <- lame(
	Y = Y,
	Xdyad = Xdyad,           # dyadic homophily indicators
	Xrow = Xrow,             # sender covariates
	Xcol = Xcol,             # receiver covariates
	family = "binary",       # binary probit model
	rvar = TRUE,             # sender random effects
	cvar = TRUE,             # receiver random effects
	dcor = TRUE,             # dyadic correlation (reciprocity)
	R = 1,                   # 1-D multiplicative latent space
	symmetric = FALSE,       # friendships are directed
	burn = 30,               # short burn-in for this worked example
	nscan = 150,             # compact post-burn-in run for the vignette
	odens = 10,              # thinning -> 15 stored draws
	posterior_opts = list(save_UV = TRUE),  # keep U/V draws: gives
	                                        # latent_positions() real
	                                        # posterior SDs below
	save_log_lik = TRUE,     # pointwise log-lik: feeds loo() at the end
	verbose = FALSE,
	plot = FALSE
)

## -----------------------------------------------------------------------------
summary(fit)

## ----intercept-decomp---------------------------------------------------------
bhat        <- colMeans(fit$BETA)
slice_means <- apply(fit$X[[1]], 3, mean, na.rm = TRUE)   # mean of each design column
centroid_lp <- sum(bhat * slice_means[names(bhat)])       # predictor at the average profile
c(intercept     = unname(bhat["intercept"]),
  centroid_lp   = centroid_lp,
  qnorm_density = qnorm(mean(unlist(Y), na.rm = TRUE)))

## ----overview-trace, fig.width=8, fig.height=8, dpi=100, dev="png", fig.alt="MCMC trace plots (top) and marginal posterior density plots (bottom) for each regression coefficient of the Dutch college fit; well-mixed fuzzy-caterpillar traces around a stable mean and smooth unimodal densities are the signatures of a converged chain."----
trace_plot(fit, params = "beta")

## ----as-draws-demo, eval = requireNamespace("posterior", quietly = TRUE)------
library(posterior)
draws <- posterior::as_draws(fit)              # draws_array [iter, chain, var]
posterior::summarise_draws(draws)              # rhat, ess_bulk, ess_tail per param

## ----overview-gof, fig.width=8, fig.height=6, dpi=100, dev="png", fig.alt="Longitudinal goodness-of-fit panels for the Dutch college fit: each facet plots one network statistic (Sender Degree Heterogeneity, Receiver Degree Heterogeneity, Dyadic Dependence, Triadic Dependence, Transitivity) across the six time periods. The observed series is a solid orange line with filled orange points; the posterior-predictive median is a dashed dark line; the 95 percent credible interval is a grey ribbon. The three elements differ on both colour and linetype, so the comparison survives grayscale printing and colour-blind viewing. Observed series falling outside the band flag structural features the model is missing."----
gof_plot(fit)

## ----overview-gof-temporal----------------------------------------------------
# density rises across the panel on net (the network gets denser as
# students settle in, though the trajectory is not monotone), so it is the
# natural temporal-trend target.
gt_density <- gof_temporal(fit, stat = "density", n_rep = 100, seed = 6886)
gt_density   # prints stat, observed slope, n_rep, and the p_pp

## ----overview-uv, fig.width=7, fig.height=6.5, dpi=100, dev="png", fig.alt="Ranked dot plot of the one-dimensional latent positions for all 32 students: for each student, an orange point marks the posterior-mean sender position and a blue point the receiver position, each with a horizontal 95 percent interval; students are ordered by sender position, and a dashed vertical line marks zero."----
lp1 <- subset(latent_positions(fit), dimension == 1)
ordv <- with(subset(lp1, type == "U"), actor[order(value)])
lp1$actor <- factor(lp1$actor, levels = ordv)

ggplot(lp1, aes(x = value, y = actor, color = type)) +
	geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
	geom_errorbar(aes(xmin = value - 2 * posterior_sd,
	                  xmax = value + 2 * posterior_sd),
	              width = 0, alpha = 0.5, orientation = "y") +
	geom_point(size = 1.8) +
	scale_color_manual(values = c(U = "#D55E00", V = "#0072B2"),
	                   labels = c(U = "sender (u)", V = "receiver (v)")) +
	labs(x = "Latent position (dimension 1)", y = NULL, color = NULL) +
	theme_bw() +
	theme(panel.border = element_blank(), panel.grid.major.y = element_blank(),
	      axis.ticks = element_blank(), legend.position = "top")

## -----------------------------------------------------------------------------
lp <- latent_positions(fit)
head(lp)

## ----tidyverse-entry, eval = requireNamespace("netify", quietly = TRUE), message = FALSE----
# library(netify)
# 
# # a tidy edgelist is often how network data arrives. Here we melt the
# # Dutch-college wave list back into one (from, to, year, tie) row per
# # ordered dyad, carrying a dyadic covariate along for the ride.
# edges_long <- do.call(rbind, lapply(seq_along(Y), function(t) {
# 	Yt  <- Y[[t]]
# 	idx <- which(!is.na(Yt), arr.ind = TRUE)
# 	data.frame(
# 		from = rownames(Yt)[idx[, 1]],
# 		to   = colnames(Yt)[idx[, 2]],
# 		year = names(Y)[t],
# 		tie  = as.integer(Yt[idx]),
# 		same_program = Xdyad[[t]][, , "same_program"][idx]
# 	)
# }))
# head(edges_long)
# 
# netlet <- netify::netify(
# 	edges_long,
# 	actor1 = "from", actor2 = "to", time = "year",
# 	weight = "tie", dyad_vars = "same_program",
# 	symmetric = FALSE, mode = "unipartite",
# 	missing_to_zero = FALSE,
# 	output_format = "longit_array"
# )
# 
# # a netify object drops straight into lame(); reduced iterations here so the
# # vignette builds quickly (use burn >= 500, nscan >= 4000 for real analyses).
# fit_netlet <- lame(
# 	netlet,
# 	family = "binary", R = 1,
# 	burn = 50, nscan = 200, odens = 5,
# 	verbose = FALSE
# )
# 
# # to inspect what lame() will see, to_lame() returns the pieces directly
# pieces <- netify::to_lame(netlet, lame = TRUE)
# names(pieces)

## ----tidy-fit-----------------------------------------------------------------
# tidy(fit) -> term / estimate / std.error / statistic / p.value /
# conf.low / conf.high; one row per coefficient (or per
# coefficient x period for dynamic_beta fits). The generic ships with
# lame so this runs without broom installed; when broom is on the path
# `broom::tidy(fit)` dispatches to the same method via
# generics::tidy registration.
tidy(fit)

## ----modelsummary-fit, eval = requireNamespace("modelsummary", quietly = TRUE) && requireNamespace("broom", quietly = TRUE)----
# # glance() returns the one-row model summary modelsummary uses for its
# # lower panel: nobs, n_actors, n_periods, n_stored, family, mode, R,
# # dynamic_uv / dynamic_ab / dynamic_beta, elpd_loo (NA unless save_log_lik=TRUE).
# broom::glance(fit)
# 
# # pass the fit straight to modelsummary -- tidy() drives the upper
# # (coefficients) panel, glance() drives the lower (GOF) panel.
# modelsummary::modelsummary(
# 	list("Dutch college" = fit),
# 	statistic = "conf.int",
# 	gof_map = c("nobs", "n_actors", "n_periods", "n_stored",
# 	            "family", "R", "dynamic_uv", "dynamic_ab")
# )

## ----fit-nouv, message = FALSE, warning = FALSE-------------------------------
fit_noUV <- lame(
	Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
	family = "binary", rvar = TRUE, cvar = TRUE, dcor = TRUE,
	R = 0,                          # no latent space
	symmetric = FALSE,
	burn = 30, nscan = 150, odens = 10,
	save_log_lik = TRUE,            # reused by loo_compare() below
	verbose = FALSE, plot = FALSE
)

## ----modelsummary-multi, eval = requireNamespace("modelsummary", quietly = TRUE) && requireNamespace("broom", quietly = TRUE), message = FALSE, warning = FALSE----
# modelsummary::modelsummary(
# 	list("AME (R = 1)" = fit, "Additive only (R = 0)" = fit_noUV),
# 	statistic = "conf.int",
# 	gof_map = c("nobs", "n_actors", "n_periods", "n_stored",
# 	            "family", "R", "dynamic_uv", "dynamic_ab")
# )

## ----autoplot-fit, fig.width = 6, fig.height = 4, fig.alt="Horizontal coefficient plot returned by autoplot.lame: each row is a regression coefficient, with a point at the posterior mean and a horizontal segment marking the 95 percent credible interval; a dashed reference line at zero indicates where an effect is indistinguishable from null."----
# autoplot.lame returns a horizontal coefplot for static fits and a
# ribbon-per-period for dynamic_beta fits, so the same call works for
# both. ggplot layers compose on top.
autoplot(fit) +
	ggtitle("Dutch college friendship: posterior coefficients")

## ----autoplot-ab, fig.width = 6, fig.height = 4, fig.alt="Lollipop plot of sender random effects: one stem per student, length equal to the posterior-mean sender effect a_i; long positive stems mark students who nominate far more friends than their covariates predict."----
# which = "ab" delegates to ab_plot() (sender side by default; use
# ab_plot(fit, effect = "receiver") for the column side). Each lollipop
# is a student's sender random effect a_i: the longest positive stems
# are students who nominate far more friends than their covariates
# predict -- the heterogeneity that va summarises as one number.
autoplot(fit, which = "ab")

## ----pdl, message = FALSE-----------------------------------------------------
# prediction_draws_long() returns a long-format data frame with
# .chain / .iteration / .draw / period / period_label / i / j /
# actor_i / actor_j / .value -- ready for tidybayes / marginaleffects.
pdl <- prediction_draws_long(fit, type = "response", n_draws = 50)
# the diagonal (i == j) is a self-tie, undefined in a friendship network, but
# the linear predictor is still numerically defined there -- drop it so the
# frame covers only real (off-diagonal) pairs.
pdl <- pdl[pdl$i != pdl$j, ]
head(pdl, 3)
dim(pdl)

## ----loo-compare-demo, eval = requireNamespace("loo", quietly = TRUE), message = FALSE, warning = FALSE----
# both fits already carry save_log_lik = TRUE, so this reuses them -- no
# refitting. Pass a NAMED list so the rows are labelled by model rather
# than model1/model2.
cmp <- loo::loo_compare(list(no_latent = loo::loo(fit_noUV),
                             latent_R1 = loo::loo(fit)))
cmp

## ----loo-compare-read, echo = FALSE, eval = requireNamespace("loo", quietly = TRUE)----
best  <- rownames(cmp)[1]
worse <- rownames(cmp)[2]
dval  <- abs(round(cmp[2, "elpd_diff"], 1))
sval  <- round(cmp[2, "se_diff"], 1)

## ----repro-seed, eval = FALSE-------------------------------------------------
# fit <- lame(..., seed = 6886)          # the sampler seed lives in lame()
# saveRDS(list(fit = fit, session = sessionInfo()),
#         file = "lame_fit_replication.rds")

## ----repro-sessioninfo--------------------------------------------------------
sessionInfo()

