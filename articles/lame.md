# Getting Started with lame

## What Is lame?

Networks of trade, conflict, alliance, friendship, and sanction share a
common inferential difficulty: the tie between two actors is rarely
independent of the ties around it. A more central actor sends more,
attracts more, and shifts the incentives of every actor sharing a
neighbour. The `lame` package fits **L**ongitudinal **A**dditive and
**M**ultiplicative **E**ffects models to this kind of data. It estimates
covariate effects on the directed tie scale while accounting for sender
heterogeneity, receiver heterogeneity, reciprocity, and the residual
relational structure that covariates alone cannot capture.

The model decomposes each tie into three pieces. Every actor has a
tendency to send ties ($`a_i`$) and to receive ties ($`b_j`$),
reflecting how active and how popular they are. Each actor also occupies
a latent position ($`u_i`$ as a sender, $`v_j`$ as a receiver) so that
actors with similar positions tend to form similar tie patterns.
Covariates shift the baseline probability of ties on top of those
actor-level pieces. When observations span multiple periods, any of
these components can evolve over time via an AR(1) process – short for
“autoregressive of order 1”, meaning each period’s value is a slightly
noisy copy of the previous period’s value (so positions drift smoothly
rather than jumping around).

## Coming from a tidy edgelist? Use `netify::to_amen(lame = TRUE)`

If your data already lives in a long-format edgelist (one row per sender
/ receiver / period), the cleanest way into `lame` is via the
[`netify`](https://netify-dev.github.io/netify/) package, which the
`lame` authors maintain as part of the same `netify-verse`. The
`netify::to_amen()` bridge produces the exact list shape (`Y`, `Xdyad`,
`Xrow`, `Xcol`) that
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) and
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) consume,
with `lame = TRUE` preserving the per-period dyadic-covariate axis
required for longitudinal models:

``` r

library(netify)
library(lame)

# 1. tidy edgelist -> netify object (one row per sender x receiver x year)
netlet <- netify::netify(
    edges_long,
    actor1 = "from", actor2 = "to", time = "year",
    weight = "tie", dyad_vars = "distance",
    symmetric = FALSE, mode = "unipartite",
    output_format = "longit_array"
)

# 2. netify -> lame-ready list (lame = TRUE keeps the per-period Xdyad axis)
amen_in <- netify::to_amen(netlet, lame = TRUE)

# 3. drop straight into lame()
fit <- lame(
    Y     = amen_in$Y,
    Xdyad = amen_in$Xdyad,
    Xrow  = amen_in$Xrow,
    Xcol  = amen_in$Xcol,
    family = "binary", R = 2,
    burn = 500, nscan = 4000, odens = 10,
    verbose = FALSE
)
```

### Pure base-R long-format CSV -\> 3-D array (no `netify` dependency)

If you would rather not pull in `netify`, the same long-format edgelist
can be reshaped to the `Y` / `Xdyad` lists
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) expects
with base R only. The pattern below is copy-pasteable: write a tiny CSV
to a tempfile, read it back, and pivot to the per-period arrays. The
reshape block is fully self-contained and evaluates here so you can
inspect the resulting list shapes; the
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) call at
the bottom is held out only so the vignette builds quickly.

``` r

set.seed(6886)
# write a tiny long-format CSV: one row per (sender, receiver, year)
edges_long <- data.frame(
    from     = rep(c("A", "B", "C", "D"), each = 8),
    to       = rep(rep(c("A", "B", "C", "D"), each = 2), 4),
    year     = rep(2020:2021, 16),
    tie      = sample(0:1, 32, replace = TRUE),
    distance = round(runif(32), 2)
)
csv_path <- tempfile(fileext = ".csv")
write.csv(edges_long, csv_path, row.names = FALSE)

# read back and reshape to a list-of-matrices keyed by year
dat       <- read.csv(csv_path, stringsAsFactors = FALSE)
actors    <- sort(unique(c(dat$from, dat$to)))
periods   <- sort(unique(dat$year))
n         <- length(actors)

empty_mat <- function() {
    matrix(NA_real_, n, n, dimnames = list(actors, actors))
}

Y_list <- lapply(periods, function(yr) {
    m  <- empty_mat()
    d  <- dat[dat$year == yr, , drop = FALSE]
    # vectorised fill: pair (from, to) -> tie
    idx <- cbind(match(d$from, actors), match(d$to, actors))
    m[idx] <- d$tie
    diag(m) <- NA  # self-ties undefined in a unipartite network
    m
})
names(Y_list) <- as.character(periods)

# Xdyad: one [n, n, 1] array per period, third-dim name = covariate name
Xdyad_list <- lapply(periods, function(yr) {
    m  <- matrix(0, n, n, dimnames = list(actors, actors))
    d  <- dat[dat$year == yr, , drop = FALSE]
    idx <- cbind(match(d$from, actors), match(d$to, actors))
    m[idx] <- d$distance
    array(m, dim = c(n, n, 1),
          dimnames = list(actors, actors, "distance"))
})
names(Xdyad_list) <- as.character(periods)

# inspect the shape contract lame() consumes
str(Y_list, max.level = 1)
#> List of 2
#>  $ 2020: num [1:4, 1:4] NA 1 1 1 1 NA 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ 2021: num [1:4, 1:4] NA 1 1 0 1 NA 0 0 0 1 ...
#>   ..- attr(*, "dimnames")=List of 2
str(Xdyad_list[[1]])
#>  num [1:4, 1:4, 1] 0.14 0.54 0.05 0.64 0.6 0.62 0.21 0.37 0.96 0.37 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ : chr [1:4] "A" "B" "C" "D"
#>   ..$ : chr [1:4] "A" "B" "C" "D"
#>   ..$ : chr "distance"
```

``` r

# the fit itself is identical to any other lame() call; this chunk is
# held out only so the vignette builds quickly. with the toy 4-actor
# data above, R = 0 and a small chain just confirms the pipe runs
# end-to-end; for substantive inference use your own larger panel.
fit <- lame(
    Y      = Y_list,
    Xdyad  = Xdyad_list,
    family = "binary", R = 0,
    burn = 200, nscan = 2000, odens = 10,
    verbose = FALSE
)
```

Using `netify` is not required. Anyone who can build a list of square
matrices, with optional 3-D covariate arrays, is fine. The
`netify::to_amen()` bridge is the path of least resistance from a tidy
data frame, and the base-R recipe above shows the underlying shape
contract for users who prefer to assemble the inputs themselves. If you
want a *worked* end-to-end example without writing your own CSV,
`data("YX_bin_long")` ships a 50-actor, 4-period synthetic panel as a
`[50, 50, 4]` $`Y`$ array plus a `[50, 50, 3, 4]` $`X`$ array. Despite
the dataset name, the shipped $`Y`$ holds the *latent* probit-scale
predictor rather than 0/1 ties (`range(YX_bin_long$Y, na.rm = TRUE)`
spans roughly $`-24`$ to $`20`$; the bare
[`range()`](https://rdrr.io/r/base/range.html) returns `NA` because the
diagonal is `NA`); threshold at zero to recover the binary tie indicator
and slice into the list-of-matrices shape
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md)
consumes:

``` r

data("YX_bin_long")
Y_bin <- lapply(seq_len(dim(YX_bin_long$Y)[3]), function(t) {
    Yt <- (YX_bin_long$Y[, , t] > 0) * 1   # threshold latent z to 0/1
    diag(Yt) <- NA
    Yt
})
```

The broom / ggplot2 / tidybayes round-trip on the resulting fit is
demonstrated in the [overview
vignette](https://netify-dev.github.io/lame/articles/lame-overview.md).

## A 5-Minute Example

Let’s fit a model to a small longitudinal binary network with one dyadic
covariate (a bilateral similarity score that drives tie formation). The
truth is `intercept = -1.0` and the covariate effect is `beta = 0.7`.

``` r

library(lame)
set.seed(6886)

# simulate 3 time periods of a 25-node directed network with a
# dyadic similarity covariate driving the ties.
n <- 25; T_periods <- 3
true_intercept <- -1.0
true_beta      <- 0.7

# zero-padded names sort the same way alphabetically as positionally,
# so `lame()`'s internal alphabetic actor sort leaves the row order
# unchanged. with names like "N1","N10","N2",... the sort reorders the
# actors (N1, N10, N2, ...); the fit is still correct as long as Y and
# the X arrays share the same names, but the stored output comes back in
# the sorted order, not your input order.
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
                            # (note: `lame()` draws live diagnostics when plot = TRUE;
                            #  `ame()` accepts plot = for signature parity but ignores it
                            #  -- the single-period sampler has no live plotting)
)

summary(fit)
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(similarity) + a[i] + b[j] + rho*e[ji] + U[i,1:2] %*% V[j,1:2], family = 'binary'"
#> 
#> Time periods: 3 
#> Family: binary 
#> Mode: unipartite 
#> 
#> Note: STATIC fit pooled across 3 time periods --
#>   U, V, a, b are time-invariant; per-period predictions vary
#>   only through per-period covariates. For time-varying effects,
#>   refit with dynamic_uv = TRUE and/or dynamic_ab = TRUE.
#> 
#> Regression coefficients:
#> ------------------------
#>                 Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept         -1.107    0.086   -12.8       0   -1.245   -0.915 ***
#> similarity_dyad    0.755     0.05  15.071       0    0.654    0.841 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.080    0.028
#> cab   -0.004    0.019
#> vb     0.074    0.022
#> rho   -0.135    0.093
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
```

The intercept lands near the true `-1.0` and `similarity_dyad` lands
near the true `0.7`, the basic sanity check that the sampler is
recovering known parameters. The `_dyad` suffix on the coefficient name
flags that the covariate is a pair-level quantity; sender-only and
receiver-only covariates carry `_row` and `_col` suffixes respectively.
The remaining variance components describe residual actor-level and
dyad-level structure that the covariate alone does not explain:

- `va`: variance of the sender random effects, capturing how much actors
  differ in their overall tendency to send ties.
- `vb`: variance of the receiver random effects, capturing how much they
  differ in their tendency to receive.
- `cab`: covariance between sender and receiver effects within an actor,
  indicating whether prolific senders are also frequent receivers.
- `rho`: within-dyad residual correlation, capturing reciprocity beyond
  what the additive effects and covariates explain.

**A note on MCMC settings.** The three key parameters are `burn`
(iterations discarded as burn-in), `nscan` (post-burn-in iterations kept
for inference), and `odens` (thinning: keep every `odens`-th sample to
reduce autocorrelation). With the settings above, we store
`nscan / odens` = 1000 / 10 = 100 posterior samples (kept small for
vignette build time). For a final run, aim for at least 1000 stored
samples with adequate effective sample sizes
(e.g. `burn = 1000, nscan = 25000, odens = 25`).

**A note on actor naming.**
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) sorts
actors alphabetically by name to align them across time periods and
against any covariate arrays. If your names sort in a different order
than they appear positionally (e.g. `"N1", "N10", "N11", ..., "N2"`),
pass `Xdyad` / `Xrow` / `Xcol` with matching dimnames so the package can
align them, or use names that sort correctly (zero-padded like
`"N01" ... "N25"`, or non-numeric like `"Alice", "Bob"`). The package
will inherit `Y`’s names onto unnamed covariate arrays when it can, but
explicit dimnames are safer.

## What Can You Do with a Fitted Model?

`lame` objects work with all the standard R methods you’d expect:

``` r

# regression coefficients (posterior means)
coef(fit)
#>       intercept similarity_dyad 
#>      -1.1065454       0.7550219

# 95% credible intervals
confint(fit)
#>                       2.5%      97.5%
#> intercept       -1.2452112 -0.9150417
#> similarity_dyad  0.6544237  0.8412236

# broom-style one-row-per-coefficient frame; ships with lame so it works
# without broom installed and dispatches through broom::tidy(fit) when
# broom is loaded. glance(fit) gives the one-row model summary that
# modelsummary uses for its lower panel. See the overview vignette for
# the full modelsummary / tidybayes / autoplot round-trip.
tidy(fit)
#> # A tibble: 2 × 7
#>   term            estimate std.error statistic p.value conf.low conf.high
#>   <chr>              <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 intercept         -1.11     0.0865     -12.8       0   -1.25     -0.915
#> 2 similarity_dyad    0.755    0.0501      15.1       0    0.654     0.841
glance(fit)
#> # A tibble: 1 × 13
#>    nobs n_actors n_row_actors n_col_actors n_periods n_stored family mode      R
#>   <int>    <int>        <int>        <int>     <int>    <int> <chr>  <chr> <int>
#> 1  1800       25           NA           NA         3      100 binary unip…     2
#> # ℹ 4 more variables: dynamic_uv <lgl>, dynamic_ab <lgl>, dynamic_beta <lgl>,
#> #   elpd_loo <dbl>

# predicted probabilities for every dyad at every time point.
# for a `lame()` fit (panel data), `predict()` returns a *list of length T*,
# where T is the number of time periods. Each element is an n x n matrix of
# posterior-mean predicted tie probabilities (between 0 and 1, since
# family = "binary"). For a single-period `ame()` fit, predict() returns
# the n x n matrix directly, not wrapped in a list.
Y_hat <- predict(fit, type = "response")
length(Y_hat)                      # 3, one matrix per period
#> [1] 3
dim(Y_hat[[1]])                    # 25 x 25
#> [1] 25 25
cat("Predicted probability range:",
        round(range(unlist(Y_hat), na.rm = TRUE), 3), "\n")
#> Predicted probability range: 0 0.97

# residuals: same list-of-matrices shape, observed minus predicted
resid_list <- residuals(fit)
cat("Residual SD:", round(sd(unlist(resid_list), na.rm = TRUE), 3), "\n")
#> Residual SD: 0.345
```

The two regression coefficients (intercept and `similarity_dyad`) should
land close to the simulation truth. Predicted probabilities span a wide
range because dyads with high similarity have an intercept-plus-effect
linear predictor much larger than dyads with low similarity. Use
`confint(fit)` for the 95% credible intervals on the coefficients.

`type = "response"` gives predictions back on the natural scale of the
outcome (probabilities for `family = "binary"`, expected counts for
`"poisson"`, expected values for `"normal"`). `type = "link"` gives the
underlying *linear predictor* (the probit-scale value before the
$`\Phi(\cdot)`$ transform for binary), which is what you want if you
plan to compute marginal effects by hand. For h-step-ahead forecasts on
a longitudinal fit, use `predict(fit, h = K)` (see the [forecasting
vignette](https://netify-dev.github.io/lame/articles/forecasting.md)).

## Checking Your Model

Two diagnostics matter for any AME fit. The first is whether the MCMC
sampler explored the posterior adequately. The trace plots should mix
freely around a stable mean, without trends or stuck regions, and the
marginal densities should be unimodal.

``` r

trace_plot(fit, params = "beta")
```

![MCMC trace and density plots for the regression coefficients
(intercept and similarity_dyad); well-mixed traces and unimodal
densities indicate adequate
convergence.](lame_files/figure-html/trace-plot-1.png)

For this fit the two regression parameters, the intercept and
`similarity_dyad`, mix around the simulation truths of $`-1.0`$ and
$`0.7`$.

The second diagnostic is whether the model reproduces the structural
features of the observed network. Posterior-predictive goodness-of-fit
plots simulate networks from the fitted model and compare their
structural statistics to those of the observed data. The observed series
is dual-encoded as an Okabe-Ito orange (`#D55E00`) solid line with
points; the posterior-predictive median is grey-dashed and the 95%
credible interval is the grey ribbon. The colour-plus-linetype encoding
stays legible in greyscale and for colour-blind readers.

``` r

gof_plot(fit)
```

![Longitudinal goodness-of-fit panels: each facet plots one network
statistic (sender/receiver degree heterogeneity, dyadic dependence,
triadic dependence, transitivity) across time. The observed series is a
solid orange line with points (Okabe-Ito \#D55E00); the
posterior-predictive median is a dark dashed line; the 95 percent
credible interval is a grey ribbon. The dual colour-plus-linetype
encoding survives greyscale and colour-blind
viewing.](lame_files/figure-html/gof-plot-1.png)

Panel titles use human-readable names (Sender Degree Heterogeneity,
Transitivity, and so on) rather than the underlying column codes such as
`sd.rowmean` or `trans.dep`; the internal codes are what `fit$GOF`
returns directly. Because we generated this data purely from
`eta = -1 + 0.7 * X` (no planted degree heterogeneity, reciprocity, or
clustering), the model reproduces almost every structural statistic:
across the five statistics and three periods, fourteen of the fifteen
observed values sit inside their 95% posterior-predictive bands. The
lone marginal excursion (transitivity at one period) sits just past the
band edge on a statistic whose values are all near zero,
i.e. Monte-Carlo noise rather than a structural misfit. On real
friendship or trade data the picture is different: at least one panel,
typically Transitivity or Sender Degree Heterogeneity, drifts clearly
outside the band. That misfit is informative: it tells the analyst which
higher-order features AME absorbs through its actor and latent structure
and which remain unexplained. It is a substantive finding about the
data, not grounds to discard the fit.

## Visualizing Network Structure

The latent space shows which actors play similar roles in the network:

``` r

uv_plot(fit)
```

![Circular-layout latent-space plot: each actor's sender position is a
triangle and receiver position a circle on concentric rings; actors
placed near each other share similar tie
patterns.](lame_files/figure-html/uv-plot-1.png)

And the additive effects show who’s unusually active or popular:

``` r

ab_plot(fit, effect = "sender")
```

![Lollipop plot of sender effects, sorted by posterior mean; each actor
is a stem from zero to its point estimate, with positive values marking
unusually active senders.](lame_files/figure-html/ab-plot-1.png)

## Cross-Sectional Models

If you have a single network (not a time series), use
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) instead of
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md):

``` r

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
#>       intercept similarity_dyad 
#>      -1.1772019       0.7421418
```

The output and methods are the same. The difference is that
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) pools
information across time periods, giving you more precise estimates when
the structure is stable.

## Dynamic Effects

When you believe the network structure is changing over time (alliances
shifting, friendships evolving), you can let the latent positions and
additive effects drift through time. The drift is modelled as an **AR(1)
process**: each period’s value is a noisy copy of the previous period’s
value, with a persistence parameter $`\rho \in (-1, 1)`$ controlling how
strongly past predicts present ($`\rho`$ near 1 = slow drift, $`\rho`$
near 0 = each period almost independent).

``` r

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
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(similarity) + a[i] + b[j] + rho*e[ji] + U[i,1:2] %*% V[j,1:2], family = 'binary'"
#> 
#> Time periods: 3 
#> Family: binary 
#> Mode: unipartite 
#> Dynamic latent positions: enabled (rho_uv = 0.685 )
#> Dynamic additive effects: enabled (rho_ab = 0.48 )
#> 
#> Regression coefficients:
#> ------------------------
#>                 Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept         -1.096    0.196  -5.591       0   -1.489   -0.711 ***
#> similarity_dyad    0.745    0.048  15.563       0    0.656    0.826 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.104    0.063
#> cab    0.041    0.061
#> vb     0.108    0.061
#> rho    0.204    0.157
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
```

The model estimates how persistent the latent positions are
($`\rho_{uv}`$) and how persistent the additive effects are
($`\rho_{ab}`$). Values near 1 indicate slow evolution; values near 0
indicate near-independent positions across periods. The simulated data
here use a single latent structure across all three periods, varying
only the dyadic similarity covariate, so the persistence estimates
reflect a mix of the prior and whatever stable-position signal the chain
extracts from three time points. The [dynamic effects
vignette](https://netify-dev.github.io/lame/articles/dynamic_effects.md)
provides a longer-panel example in which temporal evolution is genuinely
present in the simulation and the persistence parameters are
informative.

## Bipartite Networks

For two-mode networks (students and courses, countries and treaties),
pass a rectangular matrix and set `mode = "bipartite"`. Below we
simulate a 15 row × 10 column network where a dyadic similarity score
again drives the ties.

``` r

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
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(similarity) + a[i] + b[j] + U[i,1:2] %*% G %*% V[j,1:2]', family = 'binary'"
#> 
#> Time periods: 3 
#> Family: binary 
#> Mode: bipartite 
#> 
#> Note: STATIC fit pooled across 3 time periods --
#>   U, V, a, b are time-invariant; per-period predictions vary
#>   only through per-period covariates. For time-varying effects,
#>   refit with dynamic_uv = TRUE and/or dynamic_ab = TRUE.
#> 
#> Regression coefficients:
#> ------------------------
#>                 Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept         -0.896    0.382  -2.346   0.019   -1.587   -0.246   *
#> similarity_dyad    0.655    0.092   7.099       0    0.485    0.834 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.438    0.136
#> cab    0.000    0.000
#> vb     0.563    0.210
#> rho    0.000    0.000
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
#>   Note: bipartite model (rho fixed to 0, cab fixed to 0)
```

Notice that the bipartite summary reports `cab = 0.000` and
`rho = 0.000`. These two parameters describe relationships *between* the
sender and receiver of a tie, which only makes sense when the same set
of actors can play both roles. In a bipartite network the rows and
columns are different kinds of entities (rows = students, columns =
courses; rows = donors, columns = candidates), so a row actor is never
also a column actor. There is no within-actor sender / receiver
covariance and no reciprocity from $`i \to j`$ to $`j \to i`$ to
estimate, and the package fixes both at zero. The `similarity_dyad`
coefficient should land near the simulation truth of `0.6` (with this
seed the point estimate is about 0.66, a little above the truth, and the
95% interval comfortably covers 0.6). The [bipartite
vignette](https://netify-dev.github.io/lame/articles/bipartite.md) walks
through the full workflow.

## Supported Data Types

| Family | Data | Example |
|----|----|----|
| `"normal"` | Continuous | Trade volumes, survey ratings |
| `"binary"` | 0/1 | Friendships, alliances, sanctions |
| `"ordinal"` | Ordered categories | Conflict intensity (none/threat/action) |
| `"poisson"` | Counts | Number of co-sponsored bills |
| `"tobit"` | Non-negative continuous | Trade with zeros |
| `"cbin"` | Censored binary | Friendships with nomination limits |
| `"frn"` | Fixed rank nomination | “Name your top 5 friends” |
| `"rrl"` | Row-ranked likelihood | Ranked preferences |

## Power-User Features

Five features that show up in advanced workflows but rarely in a first
model.

**1. The `R > n/3` warning.** Both
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) and
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) warn
when the latent-space rank `R` exceeds `floor(n/3)` (or
`floor(min(nA, nB)/3)` in bipartite mode). Past that cap, the
multiplicative effects start absorbing structure that belongs to the
additive effects, and the latent positions stop being interpretable. The
warning is advisory rather than an error, so users with a substantive
reason to push higher can override it. In practice `R = 2` or `R = 3` is
the right default.

**2. Long runs that may not finish: `max_seconds` + `checkpoint_path` +
`resume_from`.** When you fit on a shared machine or a job queue with a
wall-clock limit, pass `max_seconds = S` and `checkpoint_path = "X.rds"`
to [`lame()`](https://netify-dev.github.io/lame/reference/lame.md). The
sampler writes a checkpoint every `checkpoint_every = 100L` iterations
(override per fit) and terminates cleanly when the wall clock hits `S`,
leaving `fit$terminated_early = TRUE`. Pick up exactly where you left
off by passing `resume_from = "X.rds"` to
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) (the
legacy `lame_resume("X.rds", nscan_more = 500)` is still supported and
produces an identical fit):

``` r

ck <- tempfile(fileext = ".rds")
fit1 <- lame(
    Y = Y_list, Xdyad = X_list, R = 2, family = "binary",
    nscan = 5000, burn = 200, odens = 25,
    max_seconds = 30,          # stop after 30 s wall clock
    checkpoint_path = ck,      # write state here
    checkpoint_every = 100L,
    verbose = FALSE
)
if (isTRUE(fit1$terminated_early)) {
    # consolidated entry point; `nscan` on the resume call is
    # treated as `nscan_more` (additional stored draws).
    fit2 <- lame(resume_from = ck, nscan = 2000)
    # equivalent legacy form:
    # fit2 <- lame_resume(ck, nscan_more = 2000)
}
```

The continuation reseeds the RNG from the checkpoint and re-invokes
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) with the
original call arguments (with `burn = 0` so the continuation does not
re-pay burn-in). The chain restarts from default prior-driven start
values rather than the checkpoint’s posterior means, so allow a short
re-equilibration before treating the continuation samples as drawn from
the same neighbourhood. (Inspect `fit2$BETA[1:50, ]` to see the warm-up
tail.)

Two caveats for iterative resume cycles. First, `resume_from` does not
currently honour user overrides for `checkpoint_path` or `max_seconds`:
both are stripped before re-evaluating the saved call, so a chain of
resumes cannot itself be checkpointed or time-budgeted through this
argument (open both with the original
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) call if
you need a long chain of resumes). Second, the continuation runs with
`burn = 0`; `verbose = TRUE` works on the resume call (it shows the
sampling progress bar without a burn-in phase).

**3. K-panel joint posterior:
[`lame_multi()`](https://netify-dev.github.io/lame/reference/lame_multi.md).**
When you have K parallel networks (K country panels, K experimental
conditions, K schools) and want one shared coefficient posterior, fit
each panel independently with
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) and pool
the per-panel beta posteriors into a precision-weighted shared
posterior:

``` r

# Each Y_panel* is itself a list (or 3-D array) of T per-panel networks
# -- i.e. the input to lame() for that panel on its own. Likewise each
# X_panel* is a list of T dyadic covariate arrays. Passing a bare matrix
# per panel will error from lame().
fit_multi <- lame_multi(
    Y_list     = list(Y_panel1, Y_panel2, Y_panel3),
    Xdyad_list = list(X_panel1, X_panel2, X_panel3),
    family = "binary", R = 2,
    nscan = 2000, burn = 200, odens = 20,
    verbose = FALSE
)
str(fit_multi$beta_shared)        # pooled beta mean (per-period if dynamic)
fit_multi$beta_deviations         # per-panel deviations from shared mean
```

This is exact under conditional independence of the panels given beta,
that is, when the panels share covariate effects but have panel-specific
actor effects, latent positions, and variance components. A single joint
MCMC across panels (full hierarchical shrinkage) requires a sampler
refactor;
[`lame_multi()`](https://netify-dev.github.io/lame/reference/lame_multi.md)
is the practical wrapper for the standard assumption.

**4. Memory-conscious
[`loo()`](https://netify-dev.github.io/lame/reference/loo.md):
`save_log_lik = "chunked"`.** For LOO or WAIC on long runs,
`save_log_lik = TRUE` stores a `[n_stored, n_obs]` log-likelihood matrix
in RAM. When that matrix would exceed memory, `save_log_lik = "chunked"`
streams per-column-chunk binary files (named `loglik_chunk_001.bin`,
`loglik_chunk_002.bin`, …) to the directory in `log_lik_path` (default:
a fresh subdirectory under
[`tempdir()`](https://rdrr.io/r/base/tempfile.html)), and `loo(fit)`
reads them back transparently. The on-disk layout is row-major within
each chunk so the readback is one `readBin` per chunk. Override the
per-file column width with `log_lik_chunk_size` (default `10000L`). See
`vignettes/dynamic_effects.Rmd` for a worked example.

**5. Multi-chain diagnostics:
[`lame_parallel()`](https://netify-dev.github.io/lame/reference/lame_parallel.md)
and
[`rhat_dynamic_beta()`](https://netify-dev.github.io/lame/reference/rhat_dynamic_beta.md).**
Single-chain $`\hat R`$ is a within-chain split statistic only. For an
honest between-chain $`\hat R`$, run K chains with different seeds via
`lame_parallel(..., n_chains = 4, combine_method = "list")`, then call
[`rhat_dynamic_beta()`](https://netify-dev.github.io/lame/reference/rhat_dynamic_beta.md)
on the returned list to get both the multivariate (Brooks-Gelman)
$`\hat R`$ on the length-T coefficient path and the max univariate
split-$`\hat R`$ over the T per-period scalars.
[`lame_parallel()`](https://netify-dev.github.io/lame/reference/lame_parallel.md)
(a thin wrapper around
[`ame_parallel()`](https://netify-dev.github.io/lame/reference/ame_parallel.md)
with `fitter = "lame"`) will warn and clamp `cores` if you ask for more
than the machine’s logical-core count, the most common silent foot-gun
on shared workers. The pooled-chain fit (`combine_method = "pool"`) also
routes through
[`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html)
so
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
gives per-`coef[t]` `ess_bulk` / `ess_tail`. The full demo lives in
`vignettes/dynamic_effects.Rmd`.

**6. Held-out predictive scoring:
[`evaluate_heldout()`](https://netify-dev.github.io/lame/reference/evaluate_heldout.md).**
For honest predictive validation, mask a random subset of dyads to `NA`
in `Y`, refit (the sampler imputes `NA` cells via data augmentation),
then score the masked cells with
`evaluate_heldout(y_obs, y_pred, mask, family = "binary")`. The helper
picks the right rule from the family: AUROC / PR-AUC / Brier / log-loss
for `binary` / `cbin`, RMSE / MAE / mean log-density for `normal` /
`tobit`, and RMSE + Poisson mean log-density for `poisson`. AUROC /
PR-AUC require (or AUROC alone via ); both are optional and the helper
degrades gracefully when neither is installed. See the [cross-sectional
vignette](https://netify-dev.github.io/lame/articles/cross_sec_ame.md)
for a worked end-to-end example.

## Where to Go Next

| I want to… | Read this |
|----|----|
| Understand the full workflow with real data | [lame overview](https://netify-dev.github.io/lame/articles/lame-overview.md) |
| Learn about cross-sectional models in depth | [Your first AME model](https://netify-dev.github.io/lame/articles/cross_sec_ame.md) |
| Model two-mode networks | [Bipartite networks](https://netify-dev.github.io/lame/articles/bipartite.md) |
| Let the network structure evolve over time | [Dynamic effects](https://netify-dev.github.io/lame/articles/dynamic_effects.md) |

## References

Hoff, PD (2021). Additive and Multiplicative Effects Network Models.
*Statistical Science*, 36, 34–50.

Minhas, S., Dorff, C., Gallop, M. B., Foster, M., Liu, H., Tellez, J., &
Ward, M. D. (2022). Taking dyads seriously. *Political Science Research
and Methods*, 10(4), 703–721.
