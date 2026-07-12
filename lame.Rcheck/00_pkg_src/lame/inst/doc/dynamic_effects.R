## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>",
	fig.width = 7,
	fig.height = 5
)

## ----dynamic-beta-example, message=FALSE--------------------------------------
library(lame)
set.seed(2026)

n_db <- 15
T_db <- 5
beta_t_true <- seq(-0.5, 1.0, length.out = T_db)   # truth: rising effect

X_db <- replicate(T_db, matrix(rnorm(n_db * n_db), n_db, n_db), simplify = FALSE)
Y_db <- vector("list", T_db)
for (t in seq_len(T_db)) {
	Yt <- beta_t_true[t] * X_db[[t]] + matrix(rnorm(n_db * n_db, 0, 0.4), n_db, n_db)
	diag(Yt) <- NA
	rownames(Yt) <- colnames(Yt) <- paste0("a", seq_len(n_db))
	Y_db[[t]] <- Yt
}
names(Y_db) <- paste0("t", seq_len(T_db))
X_db_arr <- lapply(X_db, function(x) array(x, c(n_db, n_db, 1)))

fit_db <- lame(Y_db, Xdyad = X_db_arr,
               family = "normal", R = 0,
               nscan = 200, burn = 50, odens = 5,
               dynamic_beta = "dyad", verbose = FALSE)

# fit$BETA is a 3-D array [n_stored x p x T]
dim(fit_db$BETA)

# coef() returns a p x T matrix of posterior means
coef(fit_db)

# confint() returns one row per coef[t] with 95% credible interval
head(confint(fit_db), 8)

## ----dynamic-beta-migration---------------------------------------------------
# check the shape before summarising coefficients
if (length(dim(fit_db$BETA)) == 3L) {
	# 3-D: take posterior mean over the iteration margin -> [p, T]
	pm <- apply(fit_db$BETA, c(2, 3), mean)
} else {
	# 2-D amen-style fit
	pm <- apply(fit_db$BETA, 2, mean)
}

# or just use coef(), which returns a [p, T] matrix for dynamic fits
# and a length-p vector for static fits -- no shape sniffing required
coef(fit_db)

## ----dynamic-beta-full-model, message=FALSE-----------------------------------
set.seed(2026)
n_f <- 25; T_f <- 5
beta_t_true <- seq(-0.5, 1.0, length.out = T_f)
a_true <- rnorm(n_f, 0, 0.5)       # sender heterogeneity (variance 0.25)
b_true <- rnorm(n_f, 0, 0.5)       # receiver heterogeneity
u_true <- rnorm(n_f, 0, 0.7)       # 1-D latent positions
v_true <- rnorm(n_f, 0, 0.7)

X_f <- lapply(seq_len(T_f), function(t)
	array(matrix(rnorm(n_f * n_f), n_f, n_f), c(n_f, n_f, 1),
	      dimnames = list(NULL, NULL, "x1")))
Y_f <- lapply(seq_len(T_f), function(t) {
	eta <- -0.3 + beta_t_true[t] * X_f[[t]][, , 1] +
		outer(a_true, b_true, "+") + outer(u_true, v_true)
	Yt <- matrix(rbinom(n_f * n_f, 1, pnorm(eta)), n_f, n_f)
	diag(Yt) <- NA
	rownames(Yt) <- colnames(Yt) <- paste0("a", seq_len(n_f))
	Yt
})
names(Y_f) <- paste0("t", seq_len(T_f))

fit_full <- lame(Y_f, Xdyad = X_f, family = "binary", R = 1,
                 rvar = TRUE, cvar = TRUE, dcor = TRUE,
                 dynamic_beta = "dyad",
                 nscan = 600, burn = 150, odens = 5, verbose = FALSE)

# the time-varying coefficient is recovered with the latent structure present
round(rbind(recovered = coef(fit_full)["x1_dyad", ],
            truth     = beta_t_true), 2)

# and the additive variance components stay pinned to the simulated truth
round(t(apply(fit_full$VC, 2, quantile,
              c(0.025, 0.5, 0.975))[, c("va", "vb")]), 2)

## ----dynamic-beta-kinds-------------------------------------------------------
kinds <- c("rw1", "rw2", "matern32")
paths <- sapply(kinds, function(k) {
	fit_k <- lame(Y_db, Xdyad = X_db_arr, family = "normal", R = 0,
	              nscan = 200, burn = 50, odens = 5,
	              dynamic_beta = "dyad", dynamic_beta_kind = k,
	              prior = if (k == "matern32")
	              	list(matern32_length_scale = 2) else NULL,
	              verbose = FALSE)
	coef(fit_k)["X1_dyad", ]
})
round(rbind(ar1 = coef(fit_db)["X1_dyad", ], t(paths),
            truth = beta_t_true), 3)

## ----dynamic-beta-time-index--------------------------------------------------
# observations at t = 1, 2, 4, 8 (a doubling gap structure)
fit_gaps <- lame(Y_db, Xdyad = X_db_arr, family = "normal", R = 0,
                 nscan = 200, burn = 50, odens = 5,
                 dynamic_beta = "dyad",
                 time_index = c(1, 2, 4, 8, 16)[seq_len(T_db)],
                 verbose = FALSE)
fit_gaps$time_index
round(coef(fit_gaps), 3)

## ----dynamic-beta-pool--------------------------------------------------------
fit_pool <- lame(Y_db, Xdyad = X_db_arr, family = "normal", R = 0,
                 nscan = 200, burn = 50, odens = 5,
                 dynamic_beta = c("intercept", "dyad"),
                 dynamic_beta_pool = "both",
                 verbose = FALSE)
fit_nopool <- lame(Y_db, Xdyad = X_db_arr, family = "normal", R = 0,
                   nscan = 200, burn = 50, odens = 5,
                   dynamic_beta = c("intercept", "dyad"),
                   dynamic_beta_pool = "none",
                   verbose = FALSE)

# posterior-mean AR(1) persistence per block, pooled vs independent
round(rbind(pooled      = colMeans(fit_pool$RHO_BETA),
            independent = colMeans(fit_nopool$RHO_BETA)), 3)

## ----dynamic-beta-per-actor---------------------------------------------------
fit_pa <- lame(Y_db, Xdyad = X_db_arr, family = "normal", R = 0,
               dynamic_beta_per_actor = "row",
               per_actor_covariate_idx = 1L,
               keep_per_actor = "summary",
               nscan = 200, burn = 50, odens = 5, verbose = FALSE)

# n_actors x T posterior-mean deviations, labelled by actor and period
round(fit_pa$theta_actor_mean[1:3, ], 2)

# the per-period sum-to-zero constraint holds exactly
round(colSums(fit_pa$theta_actor_mean), 10)

## ----change-point-demo--------------------------------------------------------
detect_change_point(fit_db, threshold_bf = 5)

## ----gof-temporal-demo--------------------------------------------------------
gof_temporal(fit_db, stat = "mean", n_rep = 200, seed = 1)

## ----prior-summary-demo-------------------------------------------------------
prior_summary(fit_db)

## ----rhat-dynamic-demo, eval = requireNamespace("posterior", quietly = TRUE)----
fit_list_db <- lame_parallel(
	Y_db, Xdyad = X_db_arr,
	family = "normal", R = 0,
	nscan = 200, burn = 50, odens = 5,
	dynamic_beta = "dyad",
	n_chains = 4, cores = 1,           # cores = 1 keeps the vignette portable
	combine_method = "list",            # one fit per chain
	verbose = FALSE
)

# multivariate (Brooks-Gelman) Rhat per coefficient on the length-T path,
# alongside the max univariate split-Rhat over the T per-period scalars.
rhat_dynamic_beta(fit_list_db)

# the pooled fit also routes through posterior::as_draws() with per-period
# coefficients flattened to "coef[t]" variable names, so summarise_draws()
# gives you ess_bulk / ess_tail per (coef, period).
fit_pool_db <- combine_ame_chains(fit_list_db)
posterior::summarise_draws(posterior::as_draws(fit_pool_db))

## ----mh-counters-demo---------------------------------------------------------
# per-block failure counts and acceptance proxies for the just-fit chain
str(fit_db$mh_counters)

## ----autoplot-example, fig.width = 6, fig.height = 4, fig.alt="Ribbon plot of the posterior path of the dyadic coefficient X1_dyad across five time periods, with the median as a line and the 95 percent credible interval as a shaded band, rising steadily from about minus 0.5 to 1."----
library(ggplot2)
autoplot(fit_db, coefs = "X1_dyad", probs = c(0.025, 0.5, 0.975)) +
	labs(title = "Posterior coefficient path",
	     subtitle = "Median line, 95% credible interval ribbon",
	     x = "Time period", y = "Coefficient value")

## ----save-log-lik-eval, eval = requireNamespace("loo", quietly = TRUE)--------
fit_ll <- lame(Y_db, Xdyad = X_db_arr, family = "normal", R = 0,
               nscan = 2000, burn = 200, odens = 5,
               dynamic_beta = "dyad",
               save_log_lik = TRUE,
               verbose = FALSE)

dim(fit_ll$log_lik)        # [n_stored, n_observations]
loo_db <- loo::loo(fit_ll)
print(loo_db)

## ----message=FALSE------------------------------------------------------------
library(lame)
set.seed(6886)

n_dyn <- 25
n_per <- 5
R_dyn <- 2

# evolving latent positions with AR(1) dynamics
U <- matrix(rnorm(n_dyn * R_dyn, 0, 1.5), n_dyn, R_dyn)
V <- matrix(rnorm(n_dyn * R_dyn, 0, 1.5), n_dyn, R_dyn)

Y_dyn <- list()
for(t in 1:n_per) {
	if(t > 1) {
		U <- 0.85 * U + matrix(rnorm(n_dyn * R_dyn, 0, 0.3), n_dyn, R_dyn)
		V <- 0.85 * V + matrix(rnorm(n_dyn * R_dyn, 0, 0.3), n_dyn, R_dyn)
	}
	eta <- -1 + U %*% t(V)
	Y_t <- matrix(rbinom(n_dyn * n_dyn, 1, pnorm(eta)), n_dyn, n_dyn)
	diag(Y_t) <- NA
	rownames(Y_t) <- colnames(Y_t) <- paste0("A", 1:n_dyn)
	Y_dyn[[t]] <- Y_t
}

# iterations are small so the vignette builds quickly
# use burn >= 1000 and nscan >= 5000 for real analyses.
fit_dyn_real <- lame(Y_dyn, R = 2,
	dynamic_uv = TRUE, dynamic_ab = TRUE,
	family = "binary",
	burn = 100, nscan = 1000, odens = 10,
	verbose = FALSE, plot = FALSE)

summary(fit_dyn_real)

## ----uv-traj-real, fig.width=7, fig.height=5, fig.alt="Trajectory plot of each actor's path through the 2D latent space across five periods under genuine AR(1) dynamics, showing coherent and gradual movement."----
uv_plot(fit_dyn_real, plot_type = "trajectory")

## ----uv-traj-highlight, fig.width=7, fig.height=5, fig.alt="Same trajectory plot as above but only four named actors (A1, A5, A10, A15) are coloured with the Okabe-Ito palette; remaining actors are drawn in light grey so the comparison reads at a glance."----
uv_plot(fit_dyn_real, plot_type = "trajectory",
        highlight = c("A1", "A5", "A10", "A15"))

## ----trace-dyn-real, fig.width=7, fig.height=5, fig.alt="MCMC trace and density plots for the sampled parameters of the dynamic-effects fit on truly correlated data, excluding the identification-fixed error variance; well-mixed traces indicate the sampler ran long enough."----
trace_plot(fit_dyn_real, exclude = "Error Variance")

## -----------------------------------------------------------------------------
set.seed(6886)
n <- 30
n_periods <- 5

# independent binary networks with no temporal structure
Y_list <- list()
for(t in 1:n_periods) {
	Y_t <- matrix(rbinom(n*n, 1, 0.2), n, n)
	diag(Y_t) <- NA
	rownames(Y_t) <- colnames(Y_t) <- paste0("Actor", 1:n)
	Y_list[[t]] <- Y_t
}

# fit with a custom prior on rho_uv to illustrate prior sensitivity
prior_custom <- list(
	rho_uv_mean = 0.95,    # expect very slow change in latent positions
	rho_uv_sd = 0.05       # tight prior
)

fit_null <- lame(
	Y = Y_list,
	R = 2,
	dynamic_uv = TRUE,
	dynamic_ab = TRUE,
	family = "binary",
	prior = prior_custom,
	burn = 100,
	nscan = 1000,
	odens = 10,
	# save per-draw latent trajectories for latent_positions() /
	# procrustes_align(per_draw = TRUE) below
	posterior_opts = list(save_UV = TRUE),
	verbose = FALSE,
	plot = FALSE
)

summary(fit_null)

## -----------------------------------------------------------------------------
cat("Independent data rho_uv:", round(mean(fit_null$rho_uv), 3), "\n")
cat("Correlated data rho_uv:", round(mean(fit_dyn_real$rho_uv), 3), "\n")

## ----uv-traj-null, fig.width=7, fig.height=5, fig.alt="Trajectory plot of latent positions for the null model with four actors highlighted in colour against the remaining 26 in light grey; the highlighted paths jump directionlessly and the whole cloud clumps near the origin, unlike the smooth dispersed paths of the genuine-dynamics fit."----
uv_plot(fit_null, plot_type = "trajectory",
        highlight = c("Actor1", "Actor5", "Actor10", "Actor15"))

## ----ab-traj-null, fig.width=7, fig.height=5, fig.alt="Trajectory plot of sender additive effects for the null model fit on independent networks; per-actor lines fluctuate without persistent trend, reflecting noise rather than temporal signal."----
ab_plot(fit_null, effect = "sender", plot_type = "trajectory")

## -----------------------------------------------------------------------------
# four small fits to compare specifications. nscan is intentionally
# small here so the vignette stays under CRAN's vignette compute
# budget; for a real comparison use nscan >= 2000 per fit and run them
# under `ame_parallel()` with 4 chains.
fit_static <- lame(Y_list, R = 2, family = "binary",
									burn = 50, nscan = 300, odens = 5,
									verbose = FALSE, plot = FALSE)

fit_uv <- lame(Y_list, R = 2, dynamic_uv = TRUE, family = "binary",
							burn = 50, nscan = 300, odens = 5,
							verbose = FALSE, plot = FALSE)

fit_ab <- lame(Y_list, R = 2, dynamic_ab = TRUE, family = "binary",
							burn = 50, nscan = 300, odens = 5,
							verbose = FALSE, plot = FALSE)

fit_full <- lame(Y_list, R = 2, dynamic_uv = TRUE, dynamic_ab = TRUE,
								family = "binary",
								burn = 50, nscan = 300, odens = 5,
								verbose = FALSE, plot = FALSE)

## -----------------------------------------------------------------------------
# compute p-values for each model and each GOF statistic
compute_pvals <- function(fit) {
	gof <- fit$GOF
	sapply(names(gof), function(stat) {
		mat <- gof[[stat]]
		obs <- mat[, 1]           # first column = observed
		sims <- mat[, -1]         # remaining = posterior predictive
		mean(colMeans(sims) >= mean(obs))
	})
}

gof_comparison <- rbind(
	Static     = compute_pvals(fit_static),
	Dynamic_UV = compute_pvals(fit_uv),
	Dynamic_AB = compute_pvals(fit_ab),
	Full       = compute_pvals(fit_full)
)
round(gof_comparison, 3)

## -----------------------------------------------------------------------------
# align latent positions across time
aligned <- procrustes_align(fit_null)
str(aligned$U)  # 3D array: actors x dimensions x time

## -----------------------------------------------------------------------------
lp <- latent_positions(fit_null, align = TRUE)
head(lp)

