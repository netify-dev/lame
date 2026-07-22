pkgname <- "lame"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "lame-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('lame')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("YX_bin")
### * YX_bin

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: YX_bin
### Title: binary relational data and covariates
### Aliases: YX_bin
### Keywords: datasets

### ** Examples


data(YX_bin)
gof_stats(YX_bin$Y) 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("YX_bin", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("YX_bin_long")
### * YX_bin_long

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: YX_bin_long
### Title: synthetic longitudinal binary relational data (latent-scale)
### Aliases: YX_bin_long
### Keywords: datasets

### ** Examples


data(YX_bin_long)
# threshold latent z to 0/1 before computing binary GOF stats
Yt <- 1 * (YX_bin_long$Y[, , 1] > 0)
diag(Yt) <- NA
gof_stats(Yt)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("YX_bin_long", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("YX_cbin")
### * YX_cbin

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: YX_cbin
### Title: Censored binary nomination data and covariates
### Aliases: YX_cbin
### Keywords: datasets

### ** Examples


data(YX_cbin)
gof_stats(YX_cbin$Y) 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("YX_cbin", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("YX_frn")
### * YX_frn

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: YX_frn
### Title: Fixed rank nomination data and covariates
### Aliases: YX_frn
### Keywords: datasets

### ** Examples


data(YX_frn)
gof_stats(YX_frn$Y) 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("YX_frn", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("YX_nrm")
### * YX_nrm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: YX_nrm
### Title: normal relational data and covariates
### Aliases: YX_nrm
### Keywords: datasets

### ** Examples


data(YX_nrm)
gof_stats(YX_nrm$Y)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("YX_nrm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("YX_ord")
### * YX_ord

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: YX_ord
### Title: ordinal relational data and covariates
### Aliases: YX_ord
### Keywords: datasets

### ** Examples


data(YX_ord)
gof_stats(YX_ord$Y)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("YX_ord", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ab_plot")
### * ab_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ab_plot
### Title: Visualize sender and receiver random effects
### Aliases: ab_plot

### ** Examples

## No test: 
# Fit an AME model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Visualize sender effects
ab_plot(fit, effect = "sender")

# Visualize receiver effects without sorting
ab_plot(fit, effect = "receiver", sorted = FALSE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ab_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("als_dynamic_beta")
### * als_dynamic_beta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: als_dynamic_beta
### Title: Penalised ALS time-varying coefficient estimate
### Aliases: als_dynamic_beta

### ** Examples

set.seed(1)
n <- 10; Tn <- 4
X <- replicate(Tn, array(rnorm(n*n*2), c(n, n, 2)), simplify = FALSE)
beta_true <- rbind(seq(-1, 1, length.out = Tn), seq(0.5, -0.5, length.out = Tn))
Y <- vector("list", Tn)
for (t in seq_len(Tn)) {
  Yt <- X[[t]][, , 1] * beta_true[1, t] + X[[t]][, , 2] * beta_true[2, t] +
        matrix(rnorm(n*n, 0, 0.2), n, n)
  diag(Yt) <- NA
  Y[[t]] <- Yt
}
fit_als <- als_dynamic_beta(Y, X, lambda = 0)     # per-period LS
fit_smooth <- als_dynamic_beta(Y, X, lambda = 10) # smoother path




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("als_dynamic_beta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ame")
### * ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ame
### Title: AME model fitting routine
### Aliases: ame

### ** Examples

## No test: 
data(YX_bin)
fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 10, nscan = 100, odens = 1,
           family = "binary", verbose = FALSE)
summary(fit)
# Note: you should run the Markov chain much longer in practice
## End(No test)
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ame_als")
### * ame_als

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ame_als
### Title: Fast (MCMC-free) AME estimation for a cross-sectional network
### Aliases: ame_als

### ** Examples

Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
coef(fit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ame_als", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ame_als_bootstrap")
### * ame_als_bootstrap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ame_als_bootstrap
### Title: Bootstrap uncertainty for the fast AME estimator
### Aliases: ame_als_bootstrap boot_ame

### ** Examples

## No test: 
Y <- replicate(6, { m <- matrix(rnorm(225), 15, 15); diag(m) <- NA; m },
               simplify = FALSE)
fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
# post-hoc bootstrap of an existing fit:
bt <- ame_als_bootstrap(fit, R = 50, type = "block",
                            seed = 1, verbose = FALSE)
# equivalent one-shot call:
# fit_b <- lame_als(Y, R = 1, family = "normal", verbose = FALSE,
#                   bootstrap = 50, bootstrap_type = "block",
#                   bootstrap_seed = 1)
print(bt)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ame_als_bootstrap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ame_als_refit")
### * ame_als_refit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ame_als_refit
### Title: Refit a fast AME model with a warm start
### Aliases: ame_als_refit

### ** Examples

Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
refit <- ame_als_refit(fit, verbose = FALSE)
coef(refit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ame_als_refit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ame_options")
### * ame_options

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ame_options
### Title: AME model fitting options
### Aliases: ame_options

### ** Examples

## No test: 
# Configure options
opts <- ame_options(verbose = TRUE, odens = 25)
opts
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ame_options", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ame_parallel")
### * ame_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ame_parallel
### Title: Run AME model with multiple parallel chains
### Aliases: ame_parallel

### ** Examples

## No test: 
# Run 2 chains sequentially
data(YX_nrm)
fit_parallel <- ame_parallel(YX_nrm$Y, Xdyad = YX_nrm$X,
                             n_chains = 2, cores = 1,
                             nscan = 100, burn = 10, odens = 1,
                             verbose = FALSE)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ame_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("as_lame_y")
### * as_lame_y

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: as_lame_y
### Title: Convert a graph object to a lame-ready adjacency matrix
### Aliases: as_lame_y

### ** Examples

# plain matrices are returned with the diagonal set to NA
m <- matrix(rbinom(25, 1, 0.4), 5, 5,
            dimnames = list(letters[1:5], letters[1:5]))
as_lame_y(m)

if (requireNamespace("igraph", quietly = TRUE)) {
  g <- igraph::sample_gnp(8, 0.3)
  Y <- as_lame_y(g)
  dim(Y)
}



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("as_lame_y", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("autoplot.lame")
### * autoplot.lame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: autoplot.lame
### Title: Ribbon plot of time-varying coefficients (or coefplot for static
###   fits)
### Aliases: autoplot.lame autoplot.ame

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            dynamic_beta = "dyad",
            nscan = 60, burn = 15, odens = 5, verbose = FALSE)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  autoplot(fit)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("autoplot.lame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("detect_change_point")
### * detect_change_point

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: detect_change_point
### Title: Detect potential change points in a dynamic_beta posterior path
### Aliases: detect_change_point

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            dynamic_beta = "dyad",
            nscan = 60, burn = 15, odens = 5, verbose = FALSE)
detect_change_point(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("detect_change_point", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dynamic_beta_prior_summary")
### * dynamic_beta_prior_summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dynamic_beta_prior_summary
### Title: Summarise the implied prior on a time-varying coefficient path
### Aliases: dynamic_beta_prior_summary

### ** Examples

## No test: 
# Default prior: AR(1) with rho_mean = 0.8, rho_sd = 0.15
s <- dynamic_beta_prior_summary(n_periods = 10, ndraws = 2000, seed = 1)
s$summary
# what's the probability that consecutive beta_t differ by more than 1?
s$prob_max_diff_gt_threshold

# Tighter prior on innovation variance
s_tight <- dynamic_beta_prior_summary(sigma_scale = 0.1, seed = 1)
s_tight$summary
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dynamic_beta_prior_summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("el2sm")
### * el2sm

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: el2sm
### Title: Edgelist to sociomatrix
### Aliases: el2sm

### ** Examples


Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
E<-sm2el(Y) 
el2sm(E) - Y 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("el2sm", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("evaluate_heldout")
### * evaluate_heldout

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: evaluate_heldout
### Title: Held-out predictive evaluation for an ame / lame fit
### Aliases: evaluate_heldout

### ** Examples

## No test: 
set.seed(1)
n <- 25; Y <- matrix(rbinom(n*n, 1, 0.3), n, n); diag(Y) <- NA
rownames(Y) <- colnames(Y) <- paste0("a", sprintf("%02d", 1:n))
# mask 20% of dyads
mask <- matrix(FALSE, n, n)
obs_idx <- which(!is.na(Y))
set.seed(1)
mask[sample(obs_idx, floor(0.2 * length(obs_idx)))] <- TRUE
Y_train <- Y; Y_train[mask] <- NA
fit <- ame(Y_train, R = 0, family = "binary",
           burn = 15, nscan = 60, odens = 5, verbose = FALSE, plot = FALSE)
y_pred <- predict(fit, type = "response")
evaluate_heldout(Y, y_pred, mask, family = "binary")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("evaluate_heldout", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("forecast_pit")
### * forecast_pit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: forecast_pit
### Title: Probability-integral-transform calibration check for h-step
###   forecasts
### Aliases: forecast_pit

### ** Examples

## No test: 
data(YX_bin_list)
Y_train <- YX_bin_list$Y[1:3]
Y_test <- YX_bin_list$Y[4]
X_train <- YX_bin_list$X[1:3]
fit <- lame(Y_train, Xdyad = X_train,
            family = "binary", R = 0,
            dynamic_beta = "dyad",
            nscan = 60, burn = 15, odens = 5, verbose = FALSE)
pit <- forecast_pit(fit, y_future = Y_test)
pit$ks_p
plot(pit)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("forecast_pit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glance.ame")
### * glance.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glance.ame
### Title: Glance method for fitted 'ame' / 'lame' objects
### Aliases: glance.ame glance.lame

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            nscan = 100, burn = 20, odens = 5, verbose = FALSE)
glance(fit)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glance.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glance.ame_als")
### * glance.ame_als

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glance.ame_als
### Title: Glance method for fitted 'ame_als' / 'lame_als' objects
### Aliases: glance.ame_als glance.lame_als

### ** Examples

## No test: 
data(YX_bin_list)
Y1 <- 1 * (YX_bin_list$Y[[1]] > 0); diag(Y1) <- NA
fit <- ame_als(Y = Y1, Xdyad = YX_bin_list$X[[1]],
               family = "binary", R = 1, verbose = FALSE)
glance(fit)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glance.ame_als", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gof")
### * gof

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gof
### Title: Compute GOF statistics from saved posterior samples
### Aliases: gof

### ** Examples

## No test: 
# Run model without GOF during fitting
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2, gof = FALSE,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Compute GOF post-hoc
gof_result <- gof(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gof", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gof_plot")
### * gof_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gof_plot
### Title: Visualize goodness-of-fit statistics for AME and LAME models
### Aliases: gof_plot

### ** Examples

## No test: 
# Fit an AME model
data(YX_nrm)
fit_ame <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, gof = TRUE,
               nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Basic GOF plot
gof_plot(fit_ame)

# Plot only degree-related statistics
gof_plot(fit_ame, statistics = c("sd.row", "sd.col"))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gof_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gof_stats")
### * gof_stats

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gof_stats
### Title: Goodness of fit statistics
### Aliases: gof_stats

### ** Examples


data(YX_nrm) 

# Auto-detect unipartite
gof_stats(YX_nrm$Y) 

# Bipartite (rectangular) networks are auto-detected; mode can also
# be given explicitly (required when a bipartite network is square)
Y_bip <- matrix(rnorm(120), 10, 12)
gof_stats(Y_bip, mode = "bipartite")

# Custom GOF function
my_stat <- function(Y) { c(my_measure = sum(Y > 0, na.rm = TRUE)) }
gof_stats(YX_nrm$Y, custom_gof = my_stat)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gof_stats", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gof_stats_bipartite")
### * gof_stats_bipartite

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gof_stats_bipartite
### Title: Goodness of fit statistics for bipartite networks
### Aliases: gof_stats_bipartite

### ** Examples

## No test: 
# Create a random bipartite network
Y <- matrix(rnorm(10*12), 10, 12)

# Calculate GOF statistics
gof_stats_bipartite(Y)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gof_stats_bipartite", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gof_stats_unipartite")
### * gof_stats_unipartite

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gof_stats_unipartite
### Title: Goodness of fit statistics for unipartite networks
### Aliases: gof_stats_unipartite

### ** Examples


# Create a random unipartite network
Y <- matrix(rnorm(100), 10, 10)
diag(Y) <- NA
gof_stats_unipartite(Y)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gof_stats_unipartite", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gof_temporal")
### * gof_temporal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gof_temporal
### Title: Posterior-predictive temporal-trend test
### Aliases: gof_temporal

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            dynamic_beta = "dyad",
            nscan = 60, burn = 15, odens = 5, verbose = FALSE)
gof_temporal(fit, stat = "density", n_rep = 50)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gof_temporal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame-package")
### * lame-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame-package
### Title: Longitudinal Additive and Multiplicative Effects Models for
###   Networks
### Aliases: lame-package
### Keywords: internal

### ** Examples

## No test: 
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, burn = 10, nscan = 100, odens = 1,
           family = "normal", verbose = FALSE)
summary(fit)
## End(No test)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame")
### * lame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame
### Title: AME model fitting routine for longitudinal relational data
### Aliases: lame

### ** Examples


data(YX_bin_list)
fit<-lame(YX_bin_list$Y,YX_bin_list$X,burn=5,nscan=5,odens=1,family="binary")
# you should run the Markov chain much longer than this

## No test: 
## Time-varying regression coefficients (dynamic_beta).
## Make every dyadic coefficient evolve as an AR(1):
fit_dyn <- lame(YX_bin_list$Y, YX_bin_list$X,
                family = "binary", R = 0,
                nscan = 60, burn = 15, odens = 5,
                dynamic_beta = "dyad")
dim(fit_dyn$BETA)        # [n_stored, p, T] -- 3-D when dynamic
coef(fit_dyn)            # [p, T] posterior-mean coefficient paths
confint(fit_dyn)         # per-period 95% credible intervals
summary(fit_dyn)         # prints a "Dynamic coefficients per period" block
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame_als")
### * lame_als

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame_als
### Title: Fast (MCMC-free) AME estimation for a longitudinal network
### Aliases: lame_als

### ** Examples

Y <- replicate(4, { m <- matrix(rnorm(400), 20, 20); diag(m) <- NA; m },
               simplify = FALSE)
fit <- lame_als(Y, R = 1, family = "normal", verbose = FALSE)
coef(fit)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame_als", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame_multi")
### * lame_multi

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame_multi
### Title: Multi-panel lame() with shared coefficients
### Aliases: lame_multi

### ** Examples

## No test: 
data(YX_bin_list)
fit_multi <- lame_multi(
  Y_list = list(YX_bin_list$Y, YX_bin_list$Y),
  Xdyad_list = list(YX_bin_list$X, YX_bin_list$X),
  family = "binary", R = 0,
  nscan = 100, burn = 25, odens = 5, verbose = FALSE)
dim(fit_multi$beta_shared)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame_multi", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame_parallel")
### * lame_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame_parallel
### Title: Run LAME (longitudinal AME) with multiple parallel chains
### Aliases: lame_parallel

### ** Examples

## No test: 
data(YX_bin_list)
fit_pll <- lame_parallel(YX_bin_list$Y, Xdyad = YX_bin_list$X,
                         family = "binary", n_chains = 2, cores = 1,
                         nscan = 50, burn = 10, odens = 5,
                         dynamic_beta = "dyad", verbose = FALSE)
dim(fit_pll$BETA)        # [iter * n_chains, p, T]
fit_pll$chain_indicator  # length iter * n_chains
if (requireNamespace("posterior", quietly = TRUE)) {
  posterior::summarise_draws(posterior::as_draws(fit_pll))
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame_resume")
### * lame_resume

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame_resume
### Title: Resume a 'lame()' MCMC run from a checkpoint
### Aliases: lame_resume

### ** Examples

## No test: 
ck <- tempfile(fileext = ".rds")
data(YX_bin_list)
# short run with very-aggressive max_seconds to force early termination
fit1 <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
             nscan = 5000, burn = 50, odens = 5,
             checkpoint_path = ck, checkpoint_every = 50L,
             max_seconds = 0.5, verbose = FALSE)
if (isTRUE(fit1$terminated_early)) {
  fit2 <- lame_resume(ck, nscan_more = 200)
  dim(fit2$BETA)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame_resume", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lame_snap_als")
### * lame_snap_als

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lame_snap_als
### Title: Fast approximate dynamic snap-shift AME estimator
### Aliases: lame_snap_als

### ** Examples

set.seed(1)
Y_bip <- lapply(seq_len(3), function(t) {
  m <- matrix(rnorm(6 * 5), 6, 5)
  rownames(m) <- paste0("r", seq_len(6))
  colnames(m) <- paste0("c", seq_len(5))
  m
})
fit_bip <- lame_snap_als(Y_bip, R = 1, mode = "bipartite",
                         max_iter = 3, verbose = FALSE)
dim(fit_bip$snap_prob_v)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lame_snap_als", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("latent_positions")
### * latent_positions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: latent_positions
### Title: Extract latent positions as a tidy data frame
### Aliases: latent_positions latent_positions.ame latent_positions.lame
###   latent_positions.ame_als

### ** Examples

## No test: 
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           burn = 5, nscan = 5, odens = 1, verbose = FALSE)
lp <- latent_positions(fit)
head(lp)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("latent_positions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ldZgbme")
### * ldZgbme

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ldZgbme
### Title: log density for GBME models
### Aliases: ldZgbme

### ** Examples

## For (overdispersed) Poisson regression, use
llYZ<-function(y,z){ dpois(y,z,log=TRUE) }




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ldZgbme", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lfo")
### * lfo

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lfo
### Title: Exact rolling-origin leave-future-out cross-validation
### Aliases: lfo

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            dynamic_beta = "dyad",
            nscan = 60, burn = 15, odens = 5, verbose = FALSE)
lfo_res <- lfo(fit, periods = 4L, refit = TRUE,
               nscan = 50, burn = 10, odens = 5, verbose = FALSE)
print(lfo_res)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lfo", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loo.ame")
### * loo.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loo.ame
### Title: Approximate leave-one-out cross-validation for AME / LAME fits
### Aliases: loo.ame loo.lame loo.ame_als

### ** Examples

## No test: 
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 0,
           nscan = 60, burn = 15, odens = 5,
           save_log_lik = TRUE, verbose = FALSE)
if (requireNamespace("loo", quietly = TRUE)) {
  loo_res <- loo::loo(fit)
  print(loo_res)
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loo.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("nodematch")
### * nodematch

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: nodematch
### Title: ERGM-style covariate helpers for ame() / lame()
### Aliases: nodematch absdiff nodefactor

### ** Examples

## No test: 
n <- 12
grp <- sample(letters[1:3], n, replace = TRUE)
age <- rnorm(n, 40, 10)
# build a single same-group homophily covariate
Xdyad <- array(nodematch(grp), dim = c(n, n, 1),
               dimnames = list(NULL, NULL, "same_group"))
# combine homophily + age difference
Xdyad <- array(c(nodematch(grp), absdiff(age)),
               dim = c(n, n, 2),
               dimnames = list(NULL, NULL, c("same_group", "age_diff")))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("nodematch", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("per_actor_slopes")
### * per_actor_slopes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: per_actor_slopes
### Title: Post-MCMC per-actor time-varying slopes
### Aliases: per_actor_slopes

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            nscan = 100, burn = 25, odens = 5, verbose = FALSE)
# post-MCMC per-row-actor slopes on the first dyadic covariate
pas <- per_actor_slopes(fit, kind = "row", lambda = 1)
dim(pas$slopes)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("per_actor_slopes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.ame")
### * plot.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.ame
### Title: Simple diagnostic plot for AME model fit
### Aliases: plot.ame

### ** Examples

## No test: 
# Fit an AME model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2, gof = TRUE,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Quick diagnostic plot (shows trace plots)
plot(fit)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.lame")
### * plot.lame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.lame
### Title: Plot diagnostics for a LAME model fit
### Aliases: plot.lame

### ** Examples

## No test: 
# Create simple longitudinal network data
set.seed(6886)
n <- 10
nms <- paste0("n", 1:n)
Y_list <- list(
  matrix(rnorm(n * n), n, n, dimnames = list(nms, nms)),
  matrix(rnorm(n * n), n, n, dimnames = list(nms, nms))
)
diag(Y_list[[1]]) <- diag(Y_list[[2]]) <- NA
fit <- lame(Y_list, family = "normal",
            nscan = 50, burn = 10, odens = 1, verbose = FALSE, plot = FALSE)

# default combined plot
plot(fit)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.lame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("predict.ame")
### * predict.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.ame
### Title: Predict method for AME models
### Aliases: predict.ame

### ** Examples

## No test: 
# Fit model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Point predictions
Y_pred <- predict(fit)

# Predictions on link scale
Y_link <- predict(fit, type = "link")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.ame")
### * print.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.ame
### Title: Print method for AME model objects
### Aliases: print.ame

### ** Examples

## No test: 
# Fit model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Display summary
print(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("procrustes_align")
### * procrustes_align

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: procrustes_align
### Title: Procrustes alignment of latent positions across time
### Aliases: procrustes_align

### ** Examples

## No test: 
data(YX_bin_list)
# YX_bin_list$Y stores latent-scale values; threshold to 0/1 first
Y_bin <- lapply(YX_bin_list$Y, function(y) 1 * (y > 0))
fit <- lame(Y_bin, Xdyad = YX_bin_list$X, R = 2,
            family = "binary", dynamic_uv = TRUE,
            burn = 5, nscan = 5, odens = 1,
            verbose = FALSE)
aligned <- procrustes_align(fit)
str(aligned$U)  # aligned 3D array [n, R, T]
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("procrustes_align", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rUV_sym_fc")
### * rUV_sym_fc

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rUV_sym_fc
### Title: Gibbs sampling of U and V
### Aliases: rUV_sym_fc

### ** Examples


U0<-matrix(rnorm(30,2),30,2) ; V0<-U0%*%diag(c(3,-2)) 
E<- U0%*%t(V0) + matrix(rnorm(30^2),30,30) 
rUV_sym_fc 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rUV_sym_fc", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rhat_dynamic_beta")
### * rhat_dynamic_beta

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rhat_dynamic_beta
### Title: Multivariate split-R-hat for dynamic_beta coefficient paths
### Aliases: rhat_dynamic_beta

### ** Examples

## No test: 
data(YX_bin_list)
# note: pass seed = to lame() -- an external set.seed() does not vary
# the sampler, so it would produce identical chains
fit_list <- lapply(c(1L, 2L), function(s) {
  lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
       dynamic_beta = "dyad", seed = s,
       nscan = 60, burn = 15, odens = 5, verbose = FALSE)
})
rhat_dynamic_beta(fit_list)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rhat_dynamic_beta", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rwish")
### * rwish

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rwish
### Title: Simulation from a Wishart distribution
### Aliases: rwish

### ** Examples


## The expectation is S0*nu

S0<-rwish(diag(3)) 

SS<-matrix(0,3,3) 
for(s in 1:1000) { SS<-SS+rwish(S0,5) }

SS/s 

S0*5





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rwish", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate.ame")
### * simulate.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate.ame
### Title: Simulate networks from a fitted AME model
### Aliases: simulate.ame

### ** Examples

## No test: 
# Fit a model
data(YX_bin)
fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 10, nscan = 100, odens = 1,
           family = "binary", verbose = FALSE)

# Simulate 10 networks from posterior
sims <- simulate(fit, nsim = 10)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate.lame")
### * simulate.lame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate.lame
### Title: Simulate longitudinal networks from a fitted LAME model
### Aliases: simulate.lame

### ** Examples

## No test: 
# Create simple longitudinal network data
set.seed(1)
n <- 10
nms <- paste0("n", 1:n)
Y_list <- list(
  matrix(rnorm(n * n), n, n, dimnames = list(nms, nms)),
  matrix(rnorm(n * n), n, n, dimnames = list(nms, nms))
)
diag(Y_list[[1]]) <- diag(Y_list[[2]]) <- NA
fit <- lame(Y_list, family = "normal",
            nscan = 50, burn = 10, odens = 1, verbose = FALSE, plot = FALSE)

# Simulate 10 network trajectories from posterior
sims <- simulate(fit, nsim = 10)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate.lame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_posterior")
### * simulate_posterior

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_posterior
### Title: Simulate posterior distributions from fitted AME model
### Aliases: simulate_posterior

### ** Examples

## No test: 
# Fit a model with multiplicative effects
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Get posterior samples of regression coefficients
beta_post <- simulate_posterior(fit, "beta", n_samples = 50)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_posterior", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sm2el")
### * sm2el

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sm2el
### Title: Sociomatrix to edgelist
### Aliases: sm2el

### ** Examples


Y<-matrix(rpois(10*10,.5),10,10) ; diag(Y)<-NA
E<-sm2el(Y) 
el2sm(E) - Y 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sm2el", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tidy.ame")
### * tidy.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tidy.ame
### Title: Tidy method for fitted 'ame' / 'lame' objects
### Aliases: tidy.ame tidy.lame

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary", R = 0,
            nscan = 100, burn = 20, odens = 5, verbose = FALSE)
tidy(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tidy.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tidy.ame_als")
### * tidy.ame_als

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tidy.ame_als
### Title: Tidy method for fitted 'ame_als' / 'lame_als' objects
### Aliases: tidy.ame_als tidy.lame_als

### ** Examples

## No test: 
data(YX_bin_list)
Y1 <- 1 * (YX_bin_list$Y[[1]] > 0); diag(Y1) <- NA
fit <- ame_als(Y = Y1, Xdyad = YX_bin_list$X[[1]],
               family = "binary", R = 1, verbose = FALSE)
tidy(fit)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tidy.ame_als", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("trace_plot")
### * trace_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: trace_plot
### Title: MCMC trace plots and density plots for AME/LAME model parameters
### Aliases: trace_plot

### ** Examples

## No test: 
# Fit an AME model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Basic trace plots for all parameters
trace_plot(fit)

# Only regression coefficients
trace_plot(fit, params = "beta")

# Only variance components
trace_plot(fit, params = "variance")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("trace_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("update.ame")
### * update.ame

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: update.ame
### Title: Update an AME / LAME fit
### Aliases: update.ame update.lame

### ** Examples

## No test: 
data(YX_bin_list)
fit <- lame(YX_bin_list$Y, YX_bin_list$X, family = "binary",
            burn = 5, nscan = 20, odens = 1, verbose = FALSE)
# toggle dynamic_beta on without rewriting the whole call
fit_dyn <- update(fit, dynamic_beta = "dyad")
dim(fit_dyn$BETA)  # 3-D now
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("update.ame", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("update.ame_als")
### * update.ame_als

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: update.ame_als
### Title: Update an 'ame_als' / 'lame_als' fit
### Aliases: update.ame_als update.lame_als

### ** Examples

## No test: 
Y <- matrix(rnorm(400), 20, 20); diag(Y) <- NA
fit <- ame_als(Y, R = 1, family = "normal", verbose = FALSE)
# warm-start refit with the same arguments (fast)
fit_w <- update(fit, max_iter = 50)
# cold-start refit with a different R (full restart)
fit_R2 <- update(fit, R = 2)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("update.ame_als", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("uv_plot")
### * uv_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: uv_plot
### Title: Visualize multiplicative effects (latent factors) from AME
###   models
### Aliases: uv_plot

### ** Examples

## No test: 
# Fit an AME model with multiplicative effects
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Basic visualization
uv_plot(fit)

# Use biplot layout
uv_plot(fit, layout = "biplot")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("uv_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("vignette_data")
### * vignette_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: vignette_data
### Title: TIES sanctions data for vignettes
### Aliases: vignette_data
### Keywords: datasets

### ** Examples

data(vignette_data)
cat("Countries:", nrow(Y[[1]]), "\n")
cat("Time periods:", length(Y), "\n")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("vignette_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
