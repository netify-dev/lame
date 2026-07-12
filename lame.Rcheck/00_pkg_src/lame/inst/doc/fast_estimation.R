## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>"
)
library(lame)

## ----benchmark_self, eval = FALSE---------------------------------------------
# set.seed(1)
# n <- 100
# Y <- matrix(rnorm(n*n), n, n); diag(Y) <- NA
# rownames(Y) <- colnames(Y) <- paste0("a", sprintf("%03d", 1:n))
# 
# t_als <- system.time(ame_als(Y, R = 2, family = "normal", verbose = FALSE))
# t_mcmc <- system.time(ame(Y, R = 2, family = "normal", burn = 500,
#                           nscan = 4000, odens = 25, verbose = FALSE, plot = FALSE))
# cat("ALS:  ", round(t_als["elapsed"], 2), "s\n",
#     "MCMC: ", round(t_mcmc["elapsed"], 2), "s\n",
#     "ratio:", round(t_mcmc["elapsed"] / t_als["elapsed"], 1), "x\n")

## -----------------------------------------------------------------------------
set.seed(1)
n  <- 40
a  <- rnorm(n, 0, 0.5); b <- rnorm(n, 0, 0.5)
Xd <- matrix(rnorm(n * n), n, n)
Y  <- 0.5 + 0.8 * Xd + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
diag(Y) <- NA

fit <- ame_als(Y, Xdyad = Xd, R = 1, family = "normal", verbose = FALSE)
fit

## -----------------------------------------------------------------------------
devs <- sapply(0:3, function(r)
	ame_als(Y, Xdyad = Xd, R = r, family = "normal",
	            verbose = FALSE)$deviance)
data.frame(R = 0:3, deviance = round(devs, 1))

## -----------------------------------------------------------------------------
set.seed(2)
Yb <- 1 * (Y > median(Y, na.rm = TRUE)); diag(Yb) <- NA
fb <- ame_als(Yb, Xdyad = Xd, R = 0, family = "binary", verbose = FALSE)
coef(fb)

## -----------------------------------------------------------------------------
bt <- ame_als_bootstrap(fb, R = 100, type = "parametric",
                            seed = 1, verbose = FALSE)
summary(bt)

## ----als-vs-mcmc, fig.width = 7, fig.height = 3.5, fig.alt="Side-by-side coefficient plot comparing the ALS bootstrap and MCMC posterior on the same binary network. Each coefficient (intercept and the single dyadic covariate) has two horizontal intervals: one for the ALS bootstrap, one for the MCMC posterior. A dashed vertical line at zero marks the no-effect reference; intervals on the same side of zero and with overlapping range indicate the two estimators agree."----
# wrap Xd into a 3-D array with an explicit slice name so the ALS and MCMC
# paths produce identical coefficient names (the auto-naming default differs
# between the two engines)
Xd_arr <- array(Xd, dim = c(nrow(Xd), ncol(Xd), 1),
                dimnames = list(NULL, NULL, "Xd"))

# refit the binary ALS fit with the bootstrap attached so tidy() returns
# bootstrap-based intervals
fb_b <- ame_als(Yb, Xdyad = Xd_arr, R = 0, family = "binary",
                verbose = FALSE, bootstrap = 100, bootstrap_seed = 1)

# fit the same data with MCMC using a short chain for a quick comparison
fit_mcmc <- ame(Yb, Xdyad = Xd_arr, family = "binary", R = 0,
                burn = 200, nscan = 1000, odens = 5,
                verbose = FALSE, plot = FALSE, gof = FALSE,
                seed = 1)

# collect the common coefficient interval columns without relying on
# vignette-time s3 dispatch through broom
coef_interval_table <- function(fit) {
    est <- coef(fit)
    ci <- confint(fit)
    terms <- intersect(names(est), rownames(ci))
    data.frame(
        term = terms,
        estimate = unname(est[terms]),
        conf.low = unname(ci[terms, 1]),
        conf.high = unname(ci[terms, 2]),
        row.names = NULL
    )
}

als_tdy  <- coef_interval_table(fb_b)
mcmc_tdy <- coef_interval_table(fit_mcmc)
df <- rbind(
    cbind(estimator = "ALS bootstrap",
          als_tdy[, c("term", "estimate", "conf.low", "conf.high")]),
    cbind(estimator = "MCMC posterior",
          mcmc_tdy[, c("term", "estimate", "conf.low", "conf.high")])
)

library(ggplot2)
ggplot(df, aes(x = estimate, y = term,
                colour = estimator, shape = estimator)) +
	geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
	geom_pointrange(aes(xmin = conf.low, xmax = conf.high),
	                position = position_dodge(width = 0.45),
	                size = 0.5) +
	scale_colour_manual(values = c("ALS bootstrap"  = "#0072B2",
	                                "MCMC posterior" = "#D55E00")) +
	scale_shape_manual(values  = c("ALS bootstrap"  = 16,
	                                "MCMC posterior" = 17)) +
	labs(x = "coefficient (probit scale)", y = NULL,
	     colour = NULL, shape = NULL) +
	theme_bw() +
	theme(panel.border    = element_blank(),
	      axis.ticks      = element_blank(),
	      legend.position = "top")

## -----------------------------------------------------------------------------
set.seed(3)
Yl <- replicate(5, {
	# intercept 0.4 + sender/receiver structure + noise (no dyadic
	# covariate is passed below, so only the intercept is recoverable)
	m <- 0.4 + outer(a, b, "+") + matrix(rnorm(n * n), n, n)
	diag(m) <- NA
	m
}, simplify = FALSE)
lf <- lame_als(Yl, R = 1, family = "normal", verbose = FALSE)
coef(lf)

