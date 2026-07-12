## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>",
	fig.width = 6,
	fig.height = 4,
	dpi = 96,
	out.width = "100%",
	fig.align = "center"
)

## ----netify-bridge, eval = FALSE----------------------------------------------
# library(netify)
# library(lame)
# 
# netlet <- netify::netify(
# 	edges_long,
# 	actor1 = "from", actor2 = "to", time = "year",
# 	weight = "tie",
# 	dyad_vars = "distance",
# 	output_format = "longit_list",
# 	missing_to_zero = FALSE
# )
# 
# fit <- lame(netlet, family = "binary", R = 2, verbose = FALSE)

## ----quick-start, message=FALSE-----------------------------------------------
library(lame)
set.seed(6886)

# simulate 3 time periods of a 25-node directed network with a
# dyadic similarity covariate driving the ties.
n <- 25; T_periods <- 3
true_intercept <- -1.0
true_beta      <- 0.7

# zero-padded names ("N01" ... "N25") sort the same way alphabetically as
# positionally, so `lame()`'s internal alphabetic actor sort leaves the
# row order unchanged. with non-padded names ("N1", "N10", "N2", ...) the
# stored output comes back in sorted order, not your input order -- the fit
# is still correct as long as Y and the X arrays share the same names.
# explicit dimnames on Y and every covariate array are the safe way to
# guarantee that alignment.
actor_names <- sprintf("N%02d", seq_len(n))

X_list <- lapply(seq_len(T_periods), function(t) {
	# 3-D array [actor x actor x covariate]; the third-dim name becomes
	# the coefficient label downstream
	x <- matrix(rnorm(n * n), n, n)
	array(x, dim = c(n, n, 1),
	      dimnames = list(actor_names, actor_names, "similarity"))
})
Y_list <- lapply(seq_len(T_periods), function(t) {
	# build the linear predictor eta = intercept + beta * X (the same
	# arithmetic a logistic / probit regression would do). then push it
	# through pnorm() -- the standard-normal CDF -- to turn it into a
	# tie probability in [0, 1]. that CDF link is what "probit" means,
	# and it is the link `family = "binary"` uses inside `lame()`. finally,
	# draw a 0/1 tie from a Bernoulli with that probability.
	#
	# the `[, , 1]` peels the first (and only) covariate slice off the
	# 3-D array, leaving an n x n matrix. third-dim index = covariate.
	eta <- true_intercept + true_beta * X_list[[t]][, , 1]
	Y   <- matrix(rbinom(n * n, 1, pnorm(eta)), n, n)
	diag(Y) <- NA   # self-ties are undefined in a unipartite network
	rownames(Y) <- colnames(Y) <- actor_names
	Y
})

# fit a longitudinal AME model with the dyadic similarity covariate.
fit <- lame(
	Y = Y_list,
	Xdyad = X_list,         # one similarity matrix per period
	R = 2,                  # 2D latent space
	family = "binary",      # probit for 0/1 networks
	burn = 100,             # burn-in (use 1000+ for real work)
	nscan = 1000,           # post-burn-in samples (use 5000+ for real work)
	odens = 10,             # thinning
	verbose = FALSE,        # suppress the progress bar / iteration log
	plot = FALSE            # don't pop up live MCMC diagnostic plots during sampling
	                        # `lame()` draws live diagnostics when plot = TRUE;
	                        #  `ame()` accepts plot = for signature parity but ignores it
	                        #  -- the single-period sampler has no live plotting
)

summary(fit)

## ----s3-methods---------------------------------------------------------------
# regression coefficients (posterior means)
coef(fit)

# 95% credible intervals
confint(fit)

# broom-style one-row-per-coefficient frame; ships with lame so it works
# without broom installed and dispatches through broom::tidy(fit) when
# broom is loaded. glance(fit) gives the one-row model summary that
# modelsummary uses for its lower panel. See the overview vignette for
# the full modelsummary / tidybayes / autoplot round-trip.
tidy(fit)

# n_row_actors / n_col_actors are bipartite-only and are NA on a unipartite
# fit like this one. elpd_loo is NA unless the model was fit with
# save_log_lik = TRUE and the loo result cached via fit$loo <- loo(fit).
glance(fit)

# predicted probabilities for every dyad at every time point.
# for a `lame()` fit (panel data), `predict()` returns a *list of length T*,
# where T is the number of time periods. Each element is an n x n matrix of
# posterior-mean predicted tie probabilities (between 0 and 1, since
# family = "binary"). For a single-period `ame()` fit, predict() returns
# the n x n matrix directly, not wrapped in a list.
Y_hat <- predict(fit, type = "response")
length(Y_hat)                      # 3, one matrix per period
dim(Y_hat[[1]])                    # 25 x 25
cat("Predicted probability range:",
		round(range(unlist(Y_hat), na.rm = TRUE), 3), "\n")

# residuals: same list-of-matrices shape, observed minus predicted
resid_list <- residuals(fit)
cat("Residual SD:", round(sd(unlist(resid_list), na.rm = TRUE), 3), "\n")

## ----trace-plot, fig.height=5, fig.alt="MCMC trace and density plots for the regression coefficients (intercept and similarity_dyad); well-mixed traces and unimodal densities indicate adequate convergence."----
trace_plot(fit, params = "beta")

## ----gof-plot, fig.width=8, fig.height=5, fig.alt="Longitudinal goodness-of-fit panels: each facet plots one network statistic (sender/receiver degree heterogeneity, dyadic dependence, triadic dependence, transitivity) across time. The observed series is a solid orange line with points (Okabe-Ito #D55E00); the posterior-predictive median is a dark dashed line; the 95 percent credible interval is a grey ribbon. The dual colour-plus-linetype encoding survives greyscale and colour-blind viewing."----
gof_plot(fit)

## ----uv-plot, fig.width=7, fig.height=6, fig.alt="Circular-layout latent-space plot: each actor's sender position is a triangle and receiver position a circle on concentric rings; actors placed near each other share similar tie patterns."----
uv_plot(fit)

## ----ab-plot, fig.height=4, fig.alt="Lollipop plot of sender effects, sorted by posterior mean; each actor is a stem from zero to its point estimate, with positive values marking unusually active senders. In this null example all stems are short and reflect sampling noise, not planted heterogeneity."----
ab_plot(fit, effect = "sender")

## ----cross-sectional, message=FALSE-------------------------------------------
fit_cs <- ame(
	Y = Y_list[[1]],         # just one time period
	Xdyad = X_list[[1]],     # 3-D array [n, n, 1] with "similarity" slice name
	R = 2,
	family = "binary",
	burn = 100,
	nscan = 1000,
	odens = 10,
	verbose = FALSE
)

coef(fit_cs)

## ----dynamic, message=FALSE---------------------------------------------------
fit_dyn <- lame(
	Y = Y_list,
	Xdyad = X_list,
	R = 2,
	dynamic_ab = TRUE,    # time-varying sociality/popularity
	dynamic_uv = TRUE,    # time-varying latent positions
	family = "binary",
	burn = 100,
	nscan = 1000,
	odens = 10,
	verbose = FALSE,
	plot = FALSE
)

summary(fit_dyn)

## ----bipartite, message=FALSE-------------------------------------------------
set.seed(42)          # seed the data simulation so the recovery is reproducible
nA <- 15; nB <- 10
row_names <- sprintf("R%02d", seq_len(nA))
col_names <- sprintf("C%02d", seq_len(nB))
X_bip <- lapply(1:3, function(t) {
	x <- matrix(rnorm(nA * nB), nA, nB)
	array(x, dim = c(nA, nB, 1),
	      dimnames = list(row_names, col_names, "similarity"))
})
Y_bip <- lapply(1:3, function(t) {
	eta <- -0.8 + 0.6 * X_bip[[t]][, , 1]
	Y   <- matrix(rbinom(nA * nB, 1, pnorm(eta)), nA, nB)
	rownames(Y) <- row_names; colnames(Y) <- col_names
	Y
})

fit_bip <- lame(
	Y = Y_bip,
	Xdyad = X_bip,
	mode = "bipartite",
	R = 2,
	family = "binary",
	burn = 100, nscan = 1000, odens = 10,
	verbose = FALSE, plot = FALSE
)

summary(fit_bip)

## ----checkpoint-demo, eval = FALSE--------------------------------------------
# ck <- tempfile(fileext = ".rds")
# fit1 <- lame(
# 	Y = Y_list, Xdyad = X_list, R = 2, family = "binary",
# 	nscan = 5000, burn = 200, odens = 25,
# 	max_seconds = 30,          # stop after 30 s wall clock
# 	checkpoint_path = ck,      # write state here
# 	verbose = FALSE
# )
# if (isTRUE(fit1$terminated_early)) {
# 	# `nscan` on the resume call = additional stored draws
# 	fit2 <- lame(resume_from = ck, nscan = 2000)
# }

