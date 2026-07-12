## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>",
	fig.align = "center"
)

## ----load-data----------------------------------------------------------------
library(lame)
library(ggplot2)
set.seed(6886)

# load the Add Health friendship network
data(addhealthc3)

# convert valued network to binary (any nomination = friendship)
Y <- (addhealthc3$Y > 0) * 1
X_nodes <- addhealthc3$X

n <- nrow(Y)
cat("Students:", n, "\n")
cat("Friendships:", sum(Y, na.rm = TRUE), "\n")
cat("Network density:", round(mean(Y, na.rm = TRUE), 3), "\n")

## ----viz-network, fig.width=10, fig.height=4, fig.alt="Two faceted histograms showing the distribution of nominations sent and nominations received per student in the Add Health friendship network, illustrating heterogeneity in sociality and popularity."----
out_degree <- rowSums(Y, na.rm = TRUE)
in_degree <- colSums(Y, na.rm = TRUE)

degree_df <- data.frame(
	Degree = c(out_degree, in_degree),
	Type = rep(c("Nominations sent", "Nominations received"), each = n)
)

ggplot(degree_df, aes(x = Degree)) +
	geom_histogram(binwidth = 1) +
	facet_wrap(~Type, ncol = 2) +
	labs(x = "Number of Ties", y = "Count") +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.ticks = element_blank(),
		legend.position = "top",
		strip.background = element_rect(fill = "black", color = "black"),
		strip.text = element_text(color = "white", hjust = 0)
	)

## ----create-covariates--------------------------------------------------------
# `outer(x, x, FUN)` builds an n x n matrix whose (i,j) entry is FUN(x[i], x[j]).
# So `outer(female, female, "==")` returns TRUE where students i and j share the
# same gender. Multiplying by 1 converts TRUE/FALSE to 1/0.
#
# homophily indicators: 1 if same, 0 if different
same_female <- outer(X_nodes[,"female"], X_nodes[,"female"], "==") * 1
same_race <- outer(X_nodes[,"race"], X_nodes[,"race"], "==") * 1
same_grade <- outer(X_nodes[,"grade"], X_nodes[,"grade"], "==") * 1
# absolute grade difference; same shape, but a continuous covariate
grade_diff <- abs(outer(X_nodes[,"grade"], X_nodes[,"grade"], "-"))

# pack into a 3D array (n x n x p)
Xdyad <- array(NA, dim = c(n, n, 4))
Xdyad[,,1] <- same_female
Xdyad[,,2] <- same_race
Xdyad[,,3] <- same_grade
Xdyad[,,4] <- grade_diff
dimnames(Xdyad)[[3]] <- c('same_female', 'same_race', 'same_grade', 'grade_diff')
for(k in 1:4) diag(Xdyad[,,k]) <- NA

# nodal covariates (sender and receiver characteristics)
Xrow <- X_nodes[, c("female", "grade")]
Xcol <- X_nodes[, c("female", "grade")]

## ----fit-model, message=FALSE, warning=FALSE----------------------------------
fit <- ame(Y,
					Xdyad = Xdyad,
					Xrow = Xrow,
					Xcol = Xcol,
					R = 2,              # 2D latent space
					family = "binary",  # probit model for 0/1 data
					rvar = TRUE,        # sender random effects
					cvar = TRUE,        # receiver random effects
					dcor = TRUE,        # dyadic correlation (reciprocity)
					burn = 500,         # burn-in (lengthen for a final run)
					nscan = 2000,       # post-burn-in iterations
					odens = 25,         # thinning
					verbose = FALSE,
					gof = TRUE)

## ----print-summary------------------------------------------------------------
summary(fit)

## ----trace-plots, fig.width=10, fig.height=8, fig.alt="MCMC trace and density plots for each regression coefficient; stable traces fluctuate around a steady mean and the density plots are smooth and unimodal."----
trace_plot(fit, params = "beta", ncol = 3)

## ----trace-rho, fig.width=8, fig.height=5, fig.alt="MCMC trace and density plot for the dyadic correlation (rho); the trace moves slowly relative to the regression coefficients, illustrating the slow mixing discussed in the text."----
trace_plot(fit, params = "variance", include = "Dyadic Correlation")

## ----as-draws-rhat-ess, eval = requireNamespace("posterior", quietly = TRUE), message = FALSE, warning = FALSE----
library(posterior)
draws <- posterior::as_draws(fit)              # draws_array [iter, chain, var]
posterior::summarise_draws(draws)              # mean, sd, q5, q95, rhat,
                                               # ess_bulk, ess_tail per param

## ----ame-parallel-rhat, eval = requireNamespace("posterior", quietly = TRUE)----
# reduced chains/iterations for a fast build; use n_chains = 4 and
# nscan >= 2000 for a convergence check you will report.
fit_mc <- ame_parallel(
	Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
	family = "binary", R = 2,
	burn = 200, nscan = 800, odens = 10,
	n_chains = 2, cores = 1,            # cores = 1 keeps the chunk portable
	combine_method = "pool",            # pooled fit; chain_indicator is preserved
	verbose = FALSE
)
posterior::summarise_draws(posterior::as_draws(fit_mc))

## ----mh-counters--------------------------------------------------------------
str(fit$mh_counters)

## ----prior-summary-cross------------------------------------------------------
prior_summary(fit)

## ----loo-cross, eval = requireNamespace("loo", quietly = TRUE), message = FALSE----
# reduced iterations for a fast build; use nscan >= 2000 for a reported result.
fit_ll <- ame(Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
              family = "binary", R = 2,
              burn = 200, nscan = 800, odens = 10,
              save_log_lik = TRUE, verbose = FALSE)

# alternative model: drop the dyadic covariate
fit_ll_null <- ame(Y, Xrow = Xrow, Xcol = Xcol,
                   family = "binary", R = 2,
                   burn = 200, nscan = 800, odens = 10,
                   save_log_lik = TRUE, verbose = FALSE)

loo_dyad <- loo::loo(fit_ll)
loo_null <- loo::loo(fit_ll_null)

loo_dyad   # read the Pareto-k table before any elpd number

loo::loo_compare(list(with_dyad = loo_dyad, no_dyad = loo_null))

## ----gof-plot, fig.width=10, fig.height=6, fig.alt="Goodness-of-fit panels comparing observed network statistics (dashed orange vertical line) against the posterior predictive histograms (grey) for sender and receiver degree heterogeneity, dyadic dependence, triadic dependence, and transitivity. Observed and predicted differ on both colour and linetype so the cue survives grayscale printing and colour-blind viewing."----
gof_plot(fit)

## ----custom-gof---------------------------------------------------------------
# classic clustering coefficient: transitive triples / two-paths
trans_ratio <- function(Y) {
	Yb <- Y; Yb[is.na(Yb)] <- 0            # NA -> 0 for matrix multiply
	YY <- Yb %*% Yb
	triangles <- sum(YY * Yb)              # transitive triples
	two_path  <- sum(YY) - sum(diag(YY))   # two-paths (potential triangles)
	c(trans_ratio = triangles / max(two_path, 1))
}

gof_custom <- gof(fit, custom_gof = trans_ratio, nsim = 50, verbose = FALSE)
obs <- gof_custom[1, "trans_ratio"]
sim <- gof_custom[-1, "trans_ratio"]

# posterior-predictive p-value (right tail): near 0 = under-predicted
c(observed = round(obs, 3), sim_mean = round(mean(sim), 3),
  pp_p_right = round(mean(sim >= obs), 3))

## ----latent-space, fig.width=8, fig.height=6, fig.alt="Biplot of the estimated 2D latent space: triangles mark each student's sender position and circles mark receiver positions on the two latent dimensions; students plotted near each other have similar friendship patterns beyond what the covariates explain."----
uv_plot(fit, layout = "biplot", show.edges = FALSE, label.nodes = FALSE)

## ----individual-effects, fig.width=10, fig.height=8, fig.alt="Two stacked lollipop charts showing the posterior-mean sender effect (sociality) and receiver effect (popularity) for each student as a point connected by a stem to a dashed zero line, sorted from the most negative effect on the left to the most positive on the right; positive values sit above the line and negative below."----
p1 <- ab_plot(fit, effect = "sender", sorted = TRUE,
							title = "Sender Effects (Sociality)")
p2 <- ab_plot(fit, effect = "receiver", sorted = TRUE,
							title = "Receiver Effects (Popularity)")

library(patchwork)
p1 / p2

## ----predictions--------------------------------------------------------------
pred_resp <- predict(fit, type = "response")

# how well does the model classify?
Y_vec <- as.vector(Y)
pred_vec <- as.vector(pred_resp)
keep <- !is.na(Y_vec)

# simple classification at threshold 0.5
pred_binary <- (pred_vec[keep] > 0.5) * 1
confusion <- table(Actual = Y_vec[keep], Predicted = pred_binary)
knitr::kable(confusion, caption = "Confusion Matrix (threshold = 0.5)")

accuracy <- sum(diag(confusion)) / sum(confusion)
baseline <- max(mean(Y_vec[keep]), 1 - mean(Y_vec[keep]))  # always-predict-modal-class
cat("\nAccuracy:", round(accuracy, 3),
    " | Baseline (predict modal class):", round(baseline, 3),
    " | Lift:", round(accuracy - baseline, 3), "\n")

## ----auc-eval, eval = requireNamespace("pROC", quietly = TRUE) && requireNamespace("precrec", quietly = TRUE)----
# pROC: ROC curve and AUROC
# install.packages("pROC")
roc_obj <- pROC::roc(response = Y_vec[keep], predictor = pred_vec[keep],
                     quiet = TRUE)
auroc   <- as.numeric(pROC::auc(roc_obj))

# precrec: both ROC and Precision-Recall AUCs (the latter is more
# informative when the positive class is rare, as in friendship networks)
# install.packages("precrec")
ev <- precrec::evalmod(scores = pred_vec[keep], labels = Y_vec[keep])
precrec::auc(ev)   # returns AUROC and PR-AUC side-by-side

cat("AUROC:", round(auroc, 3), "\n",
    "Baseline AUROC (random ranker):", 0.5, "\n",
    "Baseline PR-AUC (density of Y):", round(mean(Y_vec[keep]), 3), "\n")

## ----heldout------------------------------------------------------------------
# 80/20 stratified mask of off-diagonal dyads
set.seed(6886)
off_diag    <- which(row(Y) != col(Y) & !is.na(Y))
pos_idx     <- off_diag[Y[off_diag] == 1]
neg_idx     <- off_diag[Y[off_diag] == 0]
test_idx    <- c(sample(pos_idx, round(0.2 * length(pos_idx))),
                 sample(neg_idx, round(0.2 * length(neg_idx))))

Y_train               <- Y
Y_train[test_idx]     <- NA      # mask test dyads from the likelihood

fit_train <- ame(Y_train, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
                 R = 2, family = "binary",
                 rvar = TRUE, cvar = TRUE, dcor = TRUE,
                 burn = 200, nscan = 800, odens = 10,   # reduced for build speed
                 verbose = FALSE, gof = FALSE)

pred_train <- predict(fit_train, type = "response")

# `evaluate_heldout()` scores predictions on a logical TEST mask:
# columns n_eval, auroc, auprc, brier, logloss
test_mask           <- matrix(FALSE, nrow(Y), ncol(Y),
                              dimnames = dimnames(Y))
test_mask[test_idx] <- TRUE
evaluate_heldout(y_obs = Y, y_pred = pred_train,
                 mask = test_mask, family = "binary")
cat("Test dyads:", length(test_idx),
    " | positives:", sum(Y[test_idx] == 1), "\n")

## ----simulation, fig.width=8, fig.height=4, fig.alt="Posterior predictive check histogram of simulated network densities (grey bars) with a dashed Okabe-Ito orange vertical line marking the observed density; the encoding is dual on colour and linetype so the cue survives greyscale printing and colour-blind viewing. Overlap indicates the model reproduces the observed density."----
sims <- simulate(fit, nsim = 100)

# compare simulated vs observed density
sim_densities <- sapply(sims$Y, function(y) mean(y, na.rm = TRUE))
obs_density <- mean(Y, na.rm = TRUE)

ggplot(data.frame(density = sim_densities), aes(x = density)) +
	geom_histogram(bins = 20, fill = "grey60") +
	# dual-encode the observed density as colour AND linetype so the cue
	# survives greyscale printing and colour-blind viewing (Okabe-Ito orange).
	geom_vline(xintercept = obs_density,
	           color = "#D55E00", linewidth = 1, linetype = "dashed") +
	labs(title = "Posterior Predictive Check: Network Density",
			subtitle = "Dashed orange line = observed; grey histogram = simulated from model",
			x = "Network Density (Mean Tie Probability)", y = "Replicate Count") +
	theme_bw() +
	theme(panel.border = element_blank(), axis.ticks = element_blank(),
	      legend.position = "top")

