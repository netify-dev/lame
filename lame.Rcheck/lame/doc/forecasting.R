## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>",
	fig.width = 7,
	fig.height = 5,
	eval = TRUE
)

## ----fit----------------------------------------------------------------------
library(lame)
set.seed(2026)

n_fc <- 30
T_fc <- 6
rho_true   <- 0.7
beta_bar   <- 0.5
beta_t_true <- numeric(T_fc)
beta_t_true[1] <- 0.9
for (t in 2:T_fc) {
	beta_t_true[t] <- beta_bar + rho_true * (beta_t_true[t - 1] - beta_bar) +
		rnorm(1, 0, 0.15)
}

Xdyad <- lapply(seq_len(T_fc), function(t) {
	x <- matrix(rnorm(n_fc * n_fc), n_fc, n_fc)
	array(x, dim = c(n_fc, n_fc, 1),
	      dimnames = list(NULL, NULL, "trade"))
})
a_fc <- rnorm(n_fc, 0, 0.4)
b_fc <- rnorm(n_fc, 0, 0.4)
Y <- lapply(seq_len(T_fc), function(t) {
	eta <- -0.5 + beta_t_true[t] * Xdyad[[t]][, , 1] + outer(a_fc, b_fc, "+")
	Yt  <- matrix(rbinom(n_fc * n_fc, 1, pnorm(eta)), n_fc, n_fc)
	diag(Yt) <- NA
	rownames(Yt) <- colnames(Yt) <- sprintf("a%02d", seq_len(n_fc))
	Yt
})
names(Y) <- paste0("t", seq_len(T_fc))

round(beta_t_true, 3)                                   # the truth we recover
sapply(Y, function(y) round(mean(y, na.rm = TRUE), 3))  # per-period density

fit <- lame(
	Y, Xdyad = Xdyad,
	family = "binary", R = 0,
	dynamic_beta = "dyad",      # the trade coefficient is AR(1)
	dynamic_beta_kind = "ar1",  # mean-reverting; switch to "rw1" for drift
	nscan = 800, burn = 200, odens = 5,
	verbose = FALSE
)

# the recovered per-period coefficient tracks the simulated path
coef_path <- coef(fit)
trade_row <- grep("^trade[._]dyad$|trade", rownames(coef_path), value = TRUE)[1]
round(coef_path[trade_row, ], 3)

## ----forecast-----------------------------------------------------------------
set.seed(1)   # predict(h=) propagates the AR(1) state stochastically; seed for reproducibility
# 3-step-ahead forecast on the link scale (linear predictor)
fc_link <- predict(fit, h = 3, type = "link")
length(fc_link)         # 3 (one matrix per future period)
dim(fc_link[[1]])       # n x n  (30 x 30 here; unipartite)

# response-scale forecast applies the family inverse link per draw
fc_resp <- predict(fit, h = 3, type = "response")
# binary: each entry is a posterior-mean predicted probability
range(fc_resp[[1]], na.rm = TRUE)

# by_draw = TRUE returns the full [n, n, h, n_draws] array
fc_full <- predict(fit, h = 3, type = "response", by_draw = TRUE)
dim(fc_full)            # n x n x 3 x n_draws

# interval = "credible" returns a list of length-3 (lower / median / upper)
# matrices per period at the requested quantiles (default 95%)
fc_ci <- predict(fit, h = 3, type = "response", interval = "credible")
str(fc_ci[[1]])         # list($lower, $median, $upper), each n x n

## ----forecast-plot, fig.width = 7, fig.height = 3.4, fig.alt="Left: heatmap of predicted tie probabilities among the 30 actors for the first forecast period, darker cells indicating higher predicted probability. Right: the across-dyad average predicted probability at horizons h = 1, 2, 3 with a shaded 95 percent credible band that widens with the horizon."----
library(ggplot2)
library(patchwork)

# left: forecast heatmap for the first future period
P1 <- fc_resp[[1]]
rn <- rownames(P1); if (is.null(rn)) rn <- sprintf("a%02d", seq_len(nrow(P1)))
cn <- colnames(P1); if (is.null(cn)) cn <- sprintf("a%02d", seq_len(ncol(P1)))
heat_df <- data.frame(
	sender   = factor(rn[row(P1)], levels = rev(rn)),
	receiver = factor(cn[col(P1)], levels = cn),
	prob     = as.vector(P1)
)
p_heat <- ggplot(heat_df, aes(receiver, sender, fill = prob)) +
	geom_tile() +
	scale_fill_viridis_c(limits = c(0, 1), name = "P(tie)") +
	labs(title = "Forecast for t7", x = "Receiver", y = "Sender") +
	theme_bw(base_size = 9) +
	theme(panel.border = element_blank(),
	      axis.text = element_blank(), axis.ticks = element_blank(),
	      panel.grid = element_blank(), legend.position = "top")

# right: across-dyad average predicted probability + average 95% band
horizon_df <- data.frame(
	period = factor(paste0("t", T_fc + seq_along(fc_ci)),
	                levels = paste0("t", T_fc + seq_along(fc_ci))),
	median = sapply(fc_ci, function(m) mean(m$median, na.rm = TRUE)),
	lower  = sapply(fc_ci, function(m) mean(m$lower,  na.rm = TRUE)),
	upper  = sapply(fc_ci, function(m) mean(m$upper,  na.rm = TRUE))
)
p_horizon <- ggplot(horizon_df, aes(period, median, group = 1)) +
	geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
	geom_line() + geom_point(size = 2) +
	labs(title = "Average forecast by horizon",
	     x = "Forecast period", y = "Mean P(tie)") +
	theme_bw(base_size = 9) +
	theme(panel.border = element_blank(), axis.ticks = element_blank(),
	      legend.position = "top")

p_heat + p_horizon

## ----counterfactual-----------------------------------------------------------
set.seed(1)   # forecast draws are stochastic; seed so the printed delta summary reproduces
# baseline scenario: freeze the trade covariate at the last observed period
last_X <- Xdyad[[length(Xdyad)]]
X_future <- list(last_X, last_X, last_X)
fc_cf <- predict(fit, h = 3, type = "response", newdata = X_future)

# counterfactual: shift every dyad's trade covariate up by one unit
# going forward. Because `trade` is a standardised (mean-zero) covariate,
# an additive +1 shift -- "one standard deviation more trade for every
# pair" -- is the interpretable counterfactual; a multiplicative scaling
# would be meaningless on a centred covariate (see caution 2 below).
X_future_up <- lapply(X_future, function(x) {
	x[, , 1] <- x[, , 1] + 1
	x
})
fc_up <- predict(fit, h = 3, type = "response", newdata = X_future_up)

# delta is the per-dyad change in posterior-mean cooperation probability
# at h = 3, attributable solely to the +1 trade shift
delta_h3 <- fc_up[[3]] - fc_cf[[3]]
summary(as.vector(delta_h3))

## ----autoplot, fig.width = 7, fig.height = 4, fig.alt="In-sample posterior coefficient paths from autoplot.lame: each panel shows the posterior median line and a 95 percent credible interval ribbon across the fitted training periods (t1-t6) for one regression coefficient. This plots the estimated in-sample path, not a forecast."----
library(ggplot2)
autoplot(fit, probs = c(0.025, 0.5, 0.975)) +
	labs(title = "Posterior coefficient paths",
	     subtitle = "median + 95% interval",
	     x = "Time period", y = "Coefficient value")

## ----sparse-fail--------------------------------------------------------------
set.seed(2026)
n_sp <- 20
T_sp <- 4
beta_drift <- c(0.4, 0.9, 1.4, 1.9)   # drifts up; never mean-reverts
X_sp <- lapply(seq_len(T_sp), function(t) {
	x <- matrix(rnorm(n_sp * n_sp), n_sp, n_sp)
	array(x, dim = c(n_sp, n_sp, 1),
	      dimnames = list(NULL, NULL, "sanction_risk"))
})
Y_sp <- lapply(seq_len(T_sp), function(t) {
	eta <- -2 + beta_drift[t] * X_sp[[t]][, , 1]
	Yt  <- matrix(rbinom(n_sp * n_sp, 1, pnorm(eta)), n_sp, n_sp)
	diag(Yt) <- NA
	rownames(Yt) <- colnames(Yt) <- sprintf("s%02d", seq_len(n_sp))
	Yt
})
names(Y_sp) <- paste0("t", seq_len(T_sp))
sapply(Y_sp, function(y) round(mean(y, na.rm = TRUE), 3))  # sparse early on

fit_sparse <- lame(
	Y_sp, Xdyad = X_sp,
	family = "binary", R = 0,
	dynamic_beta = "dyad", dynamic_beta_kind = "ar1",
	nscan = 400, burn = 100, odens = 5,
	verbose = FALSE
)
round(quantile(fit_sparse$RHO_BETA, c(0.5, 0.975)), 3)

## ----sparse-fail-warn---------------------------------------------------------
set.seed(1)
fc_sparse <- predict(fit_sparse, h = 6, type = "response",
                     interval = "credible")

## ----sparse-fail-width, fig.width = 6, fig.height = 3.2, fig.alt="Line plot of the across-dyad mean 95 percent forecast interval width at horizons one through six for two fits. The mean-reverting Step 1 fit's width is flat from roughly horizon three on, the saturation an AR(1) below the unit root predicts, while the near-unit-root sparse fit's width is still rising at horizon six, roughly double its one-step value."----
set.seed(1)
fc_healthy <- predict(fit, h = 6, type = "response",
                      interval = "credible")
width_by_h <- function(fc) {
	sapply(fc, function(m) mean(m$upper - m$lower, na.rm = TRUE))
}
w_healthy <- width_by_h(fc_healthy)
w_sparse  <- width_by_h(fc_sparse)
round(rbind(mean_reverting = w_healthy, near_unit_root = w_sparse), 3)

# at h = 3, how many of the dyads the covariate actually moves (top
# decile of |sanction_risk|) get an interval spanning nearly all of
# [0, 1]?
x_last <- X_sp[[T_sp]][, , 1]
hi_x   <- abs(x_last) >= quantile(abs(x_last), 0.9)
deg_h3 <- mean(fc_sparse[[3]]$lower[hi_x] < 0.05 &
               fc_sparse[[3]]$upper[hi_x] > 0.95, na.rm = TRUE)
round(deg_h3, 2)

width_df <- data.frame(
	h     = rep(seq_along(w_healthy), 2),
	fit   = rep(c("mean-reverting (Step 1 fit)",
	              "near-unit-root (sparse fit)"), each = length(w_healthy)),
	width = c(w_healthy, w_sparse)
)
ggplot(width_df, aes(h, width, linetype = fit)) +
	geom_line() + geom_point(size = 1.8) +
	labs(title = "Forecast interval width by horizon",
	     x = "Forecast horizon h", y = "Mean 95% interval width") +
	theme_bw(base_size = 9) +
	theme(panel.border = element_blank(), axis.ticks = element_blank(),
	      legend.position = "top", legend.title = element_blank())

## ----lfo----------------------------------------------------------------------
# refit on Y[1:(t-1)] for each t in periods and score Y[[t]].
# Skip the very first leave-out (would leave a 1-period training window
# which is too short for dynamic_beta); use the last 2 origins.
set.seed(1)  # the internal h=1 forecast draws inside lfo() are stochastic; seed for reproducibility
T_fit <- length(fit$YPM)
lfo_periods <- tail(seq_len(T_fit), max(1L, min(2L, T_fit - 2L)))
lfo_res <- lfo(fit, periods = lfo_periods, refit = TRUE,
               nscan = 200, burn = 50, odens = 5, verbose = FALSE)
print(lfo_res)

## ----pit----------------------------------------------------------------------
# fit on the first five periods, hold out the sixth
fit_train <- lame(
	Y[1:5], Xdyad = Xdyad[1:5],
	family = "binary", R = 0,
	dynamic_beta = "dyad",
	nscan = 400, burn = 100, odens = 5,
	verbose = FALSE
)
set.seed(1)   # the randomised PIT draws are stochastic; seed for reproducibility
pit <- forecast_pit(fit_train, y_future = list(Y[[6]]))
print(pit)

## ----pit-plot, fig.height = 3.6, fig.alt="PIT calibration plot for the held-out period: the empirical distribution of probability-integral-transform values compared with the Uniform(0,1) reference. Bars or a curve tracking the reference indicate calibration; a U-shape signals an over-confident (too-narrow) forecast and a central hump signals an under-confident (too-wide) one."----
plot(pit)

## ----loo-compare--------------------------------------------------------------
# reduced iterations so the vignette builds quickly; use nscan >= 2000,
# burn >= 500 for a comparison you will report.
fit_static <- lame(
	Y, Xdyad = Xdyad, family = "binary", R = 0,
	nscan = 600, burn = 150, odens = 5,
	save_log_lik = TRUE,
	verbose = FALSE
)
fit_dyn <- lame(
	Y, Xdyad = Xdyad, family = "binary", R = 0,
	dynamic_beta = "dyad",
	nscan = 600, burn = 150, odens = 5,
	save_log_lik = TRUE,
	verbose = FALSE
)

# rank the two fits: the row with elpd_diff = 0 is the best,
# subsequent rows show elpd_diff and its SE relative to it. Pass a NAMED
# list so the rows are labelled by model name rather than "model1"/"model2".
cmp <- loo::loo_compare(list(fit_static = loo::loo(fit_static),
                             fit_dyn    = loo::loo(fit_dyn)))
cmp

## ----loo-compare-read, echo = FALSE-------------------------------------------
best     <- rownames(cmp)[1]
runnerup <- rownames(cmp)[2]
diff_val <- abs(round(cmp[2, "elpd_diff"], 1))
se_val   <- round(cmp[2, "se_diff"], 1)
ratio    <- round(diff_val / se_val, 1)

## ----chunked, eval = FALSE----------------------------------------------------
# fit_chk <- lame(
# 	Y, Xdyad = Xdyad, family = "binary", R = 0,
# 	dynamic_beta = "dyad",
# 	nscan = 2000, burn = 500, odens = 5,
# 	save_log_lik = "chunked",          # per-column-chunk binary files
# 	log_lik_chunk_size = 5000L,
# 	verbose = FALSE
# )
# loo::loo(fit_chk)                       # reads the chunks transparently

