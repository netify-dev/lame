# lame Overview

## Package Overview

The `lame` package provides tools for fitting **L**ongitudinal
**A**dditive and **M**ultiplicative **E**ffects models to network data
observed over time. If you study relationships between actors (countries
trading with each other, legislators co-sponsoring bills, students
forming friendships across semesters) and you have repeated observations
of those relationships, `lame` is designed for you.

The modeling approach builds on the Additive and Multiplicative Effects
(AME) framework developed by Peter Hoff, whose
[`amen`](https://pdhoff.github.io/amen/) package
([GitHub](https://github.com/pdhoff/amen)) provides the foundational
implementation for cross-sectional network analysis. The `lame` package
extends this framework in several directions:

- **Longitudinal and replicated data**: While `amen` focuses on single
  networks (with
  [`ame()`](https://netify-dev.github.io/lame/reference/ame.md)) or
  basic replicated designs (with `ame_rep()`), `lame` is built from the
  ground up for panel data, networks observed across multiple time
  periods, potentially with actors entering and leaving the sample.
- **Dynamic effects**: Actors’ positions in latent space and their
  baseline activity levels can evolve over time via AR(1) processes
  (`dynamic_uv` and `dynamic_ab`), letting you model how network
  structure changes rather than assuming it is fixed.
- **C++ acceleration**: Core sampling routines are implemented in
  Rcpp/RcppArmadillo for faster computation, which matters when you have
  many time periods or actors.
- **Bipartite networks**: Explicit support for two-mode networks (e.g.,
  countries and treaties, users and products) with separate latent
  spaces for row and column nodes.
- **ggplot2-based diagnostics**: Built-in visualization functions
  (`trace_plot`, `gof_plot`, `ab_plot`, `uv_plot`) return standard
  ggplot2 figures you can theme and extend.
- **Standard S3 methods**:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html), and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) work as you
  would expect from any R model object.

The package is part of the `netify-verse`, and it works with the
[`netify`](https://netify-dev.github.io/netify/) package for data
preparation.

### Migrating from `amen`

Both `amen` and `lame` export an
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) function.
If you load both packages in the same session, `lame` prints a
`.onAttach` startup message so you know unqualified `ame(...)` calls
dispatch to whichever package was attached last; use `lame::ame(...)` or
`amen::ame(...)` to be explicit. The table below maps the breaking
differences for old `amen` scripts:

| Topic | `amen::ame()` | [`lame::ame()`](https://netify-dev.github.io/lame/reference/ame.md) / [`lame::lame()`](https://netify-dev.github.io/lame/reference/lame.md) | Old code breaks? |
|----|----|----|----|
| Family strings | `"nrm"`, `"bin"`, `"ord"`, `"cbin"`, `"frn"`, `"rrl"`, `"tobit"` | `"normal"`, `"binary"`, `"ordinal"`, `"cbin"`, `"frn"`, `"rrl"`, `"tobit"`, `"poisson"` | No for [`ame()`](https://netify-dev.github.io/lame/reference/ame.md) / [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) (short forms `"nrm"`, `"bin"`, `"ord"`, `"pois"`, `"tob"` are accepted as aliases via [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html) – two informational lines per call, recommending the full name; old scripts run unchanged). For [`ame_als()`](https://netify-dev.github.io/lame/reference/ame_als.md) / [`lame_als()`](https://netify-dev.github.io/lame/reference/lame_als.md) (MCMC-free estimator) the same short aliases (`"nrm"`, `"bin"`, `"ord"`, `"pois"`, `"tob"`) are accepted with an informational message, but the estimator *supports* only `"normal"`, `"binary"`, and `"poisson"`. The censoring / rank families `"ordinal"`, `"tobit"`, `"cbin"`, `"frn"`, and `"rrl"` are rejected with a direct error pointing to [`ame()`](https://netify-dev.github.io/lame/reference/ame.md) / [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) (MCMC), which is the posterior-sampling path for those families. |
| Bipartite rank | single `R` | `R` (square) or separate `R_row` / `R_col` (bipartite); single `R` still accepted | No (old `R` calls work; `R_row` / `R_col` is new) |
| Verbosity flag | `print = TRUE/FALSE` | `verbose = TRUE/FALSE` | No (deprecation warning) |
| `fit$BETA` shape (static) | `[n_stored, p]` matrix | `[n_stored, p]` matrix | No (identical) |
| `fit$BETA` shape (`dynamic_beta`) | not available | `[n_stored, p, T]` 3-D array | New feature, but watch out: under the new shape `apply(fit$BETA, 2, mean)` still returns a length-`p` vector that *looks* like the old amen result, while silently averaging across both iterations and time. Use `coef(fit)` (returns `[p, T]` for dynamic, length-`p` for static) or branch on `length(dim(fit$BETA))`. See [Dynamic Effects](https://netify-dev.github.io/lame/articles/dynamic_effects.md). |
| Longitudinal data | `ame_rep()` for replicates | [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) with a list of `T` matrices and optional `dynamic_uv` / `dynamic_ab` / `dynamic_beta` | Refit via [`lame()`](https://netify-dev.github.io/lame/reference/lame.md) |
| GOF output | `fit$GOF` matrix | `fit$GOF` (cross-sectional) or list-of-matrices per period (longitudinal) | Yes for longitudinal scripts |
| S3 methods | [`summary()`](https://rdrr.io/r/base/summary.html), [`plot()`](https://rdrr.io/r/graphics/plot.default.html) | adds [`coef()`](https://rdrr.io/r/stats/coef.html), [`confint()`](https://rdrr.io/r/stats/confint.html), [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`predict()`](https://rdrr.io/r/stats/predict.html), [`simulate()`](https://rdrr.io/r/stats/simulate.html), [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html), [`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md) | No (additive) |
| `loo` / `waic` | not available | `save_log_lik = TRUE` plus `loo::loo(fit)` / `loo::waic(fit)` S3 methods; see [Dynamic Effects](https://netify-dev.github.io/lame/articles/dynamic_effects.md) | N/A (new feature) |
| [`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) | not available | [`tidy.ame()`](https://netify-dev.github.io/lame/reference/tidy.ame.md) / [`tidy.lame()`](https://netify-dev.github.io/lame/reference/tidy.ame.md) registered against [`generics::tidy`](https://generics.r-lib.org/reference/tidy.html) | N/A (new feature) |
| [`broom::glance()`](https://generics.r-lib.org/reference/glance.html) | not available | [`glance.ame()`](https://netify-dev.github.io/lame/reference/glance.ame.md) / [`glance.lame()`](https://netify-dev.github.io/lame/reference/glance.ame.md) registered against [`generics::glance`](https://generics.r-lib.org/reference/glance.html); returns `nobs`, `n_actors`, `n_periods`, `n_stored`, `family`, `mode`, `R`, `dynamic_uv`, `dynamic_ab`, `dynamic_beta`, `elpd_loo`. Pairs with [`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md) for `modelsummary::modelsummary(fit)`. | N/A (new feature) |

The most common porting fix is `family = "nrm"` → `family = "normal"`
and `family = "bin"` → `family = "binary"`. Everything else either still
works unchanged or surfaces a deprecation warning rather than an error.

## Application: Dutch College Friendships

To see how `lame` works in practice, we analyze the **Dutch college
friendship network** (van de Bunt, van Duijn, and Snijders 1999), a
panel of 32 students who arrived as strangers and rated each other’s
friendship at seven time points. (We drop the cold-start first wave,
when the students barely knew each other and almost no ties exist; six
waves remain.) The substantive question is whether students sort by
gender, smoking status, and academic program, and whether residual
clustering remains once those covariates are absorbed.

The dataset includes:

- **Y**: directed friendship ratings on a -1 to 4 scale (-1 a
  negative/troubled relationship, 0 no or uncertain tie, 1-4 increasing
  friendship). We binarize any positive rating to a tie (`(Y > 0) * 1`),
  so the rare -1 ratings fold in with 0 as “no tie”.
- **X**: three node-level attributes: `male`, `smoker`, `program`.
- We construct three dyadic **homophily indicators** from the node
  attributes: `same_male`, `same_smoker`, `same_program`.

``` r

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
```

The friendship network’s average density rises across the panel as
students settle in, though the trajectory is not monotone (density
bounces between roughly 0.40 and 0.63 across the six waves):

``` r

sapply(Y, function(y) round(mean(y, na.rm = TRUE), 2))
#>   t2   t3   t4   t5   t6   t7 
#> 0.40 0.52 0.47 0.49 0.63 0.51
```

### Fitting the Model

We fit a binary probit AME model with sender random effects (`rvar`),
receiver random effects (`cvar`), dyadic correlation (`dcor`), and a
one-dimensional latent space (`R = 1`). The random effects capture
sender heterogeneity (some students simply nominate more friends) and
receiver heterogeneity (some students are nominated more often), while
the latent space picks up residual clustering: students who befriend
similar people for reasons that the observed covariates do not capture.
We use `R = 1` because with only 32 students a two-dimensional latent
space competes with the three dyadic homophily indicators for the same
signal and the chain mixes poorly on the dyadic coefficients; `R = 1`
resolves that cleanly and the substantive results are nearly identical.

``` r

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
    burn = 500,              # burn-in
    nscan = 4000,            # post-burn-in iterations
    odens = 10,              # thinning -> 400 stored draws
    verbose = FALSE,
    plot = FALSE
)
```

### Interpreting the Results

``` r

summary(fit)
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(same_male, same_smoker, same_program) + row(male, smoker, program) + col(male, smoker, program) + a[i] + b[j] + rho*e[ji] + U[i,1:1] %*% V[j,1:1], family = 'binary'"
#> 
#> Time periods: 6 
#> Family: binary 
#> Mode: unipartite 
#> 
#> Note: STATIC fit pooled across 6 time periods --
#>   U, V, a, b are time-invariant; per-period predictions vary
#>   only through per-period covariates. For time-varying effects,
#>   refit with dynamic_uv = TRUE and/or dynamic_ab = TRUE.
#> 
#> Regression coefficients:
#> ------------------------
#>                   Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept           -2.869    1.208  -2.376   0.018   -5.241   -0.572   *
#> male_row             0.939    0.417   2.251   0.024    0.212    1.761   *
#> smoker_row          -0.951     0.34  -2.799   0.005   -1.583   -0.311  **
#> program_row          0.497    0.231   2.153   0.031    0.059    0.947   *
#> male_col             0.991    0.267   3.708       0    0.446    1.516 ***
#> smoker_col          -0.062    0.217  -0.286   0.775   -0.525    0.341    
#> program_col          0.149    0.145   1.028   0.304   -0.111    0.437    
#> same_male_dyad       0.446    0.058   7.654       0    0.333    0.566 ***
#> same_smoker_dyad     0.143    0.056   2.551   0.011    0.037    0.244   *
#> same_program_dyad    0.577    0.053  10.948       0    0.472    0.669 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: stars are a visual hint from posterior mean / SD only; for inference use the credible intervals.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.797    0.217
#> cab    0.291    0.107
#> vb     0.268    0.078
#> rho    0.321    0.042
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
```

A lot lands at once, so let’s walk through it.

**Intercept.** The probit intercept is the latent linear predictor when
*every covariate is zero* and every random effect sits at its mean. That
is not the same as the average student: the covariates here are
uncentered, and `program` in particular is coded 2-4 and is never zero,
so “all covariates zero” extrapolates well outside the observed data.
The familiar `qnorm(density)` rule of thumb describes the predictor at
the *average* covariate profile, not at the all-zero corner. The two
reconcile exactly once you add back the mean covariate contribution:

``` r

bhat        <- colMeans(fit$BETA)
slice_means <- apply(fit$X[[1]], 3, mean, na.rm = TRUE)   # mean of each design column
centroid_lp <- sum(bhat * slice_means[names(bhat)])       # predictor at the average profile
c(intercept     = unname(bhat["intercept"]),
  centroid_lp   = centroid_lp,
  qnorm_density = qnorm(mean(unlist(Y), na.rm = TRUE)))
#>     intercept   centroid_lp qnorm_density 
#>   -2.86919980   -0.08283855    0.00842291
```

The intercept itself is strongly negative, but the predictor at the
average covariate profile (`centroid_lp`) sits close to `qnorm(density)`
– the uncentered covariates account for nearly the entire gap, leaving
only a small shift from the random-effect variance. The intercept
absorbs the mean covariate contribution; it is not a baseline
probability. If you want the intercept to read directly as a baseline,
center the covariates before fitting. Either way, read predicted
probabilities on the response scale off
`predict(fit, type = "response")`.

**Homophily.** The three dyadic coefficients (`same_male`,
`same_smoker`, `same_program`) measure whether students of the same type
are more likely to be friends, holding the sender, receiver, and latent
positions fixed. Positive values mean homophily; credible intervals that
exclude zero indicate that the direction of the association is
well-identified. Program homophily tends to be the strongest of the
three in this network, consistent with the substantive expectation that
students bond around shared coursework.

**Nodal covariates.** `male_row` / `male_col` measure whether male
students are unusually active senders or popular receivers, and
similarly for `smoker` and `program`. With only 32 students, the
credible intervals on individual nodal effects are wide; treat the sign
as suggestive, the magnitude with caution.

**Variance components.** `va` (sender variance) and `vb` (receiver
variance) quantify how much students differ in sociability and
popularity beyond what the covariates explain. The dyadic correlation
`rho` captures reciprocity: values away from zero indicate that ties
tend to be mutual. In friendship networks `rho > 0` is the rule rather
than the exception.

### Checking Convergence

The model is estimated via MCMC, so trace plots and posterior-density
plots are useful checks. The `trace_plot` function shows the sampled
values over iterations (top panels) and the corresponding posterior
densities (bottom panels) for each parameter.

``` r

trace_plot(fit, params = "beta")
```

![MCMC trace plots (top) and marginal posterior density plots (bottom)
for each regression coefficient of the Dutch college fit; well-mixed
fuzzy-caterpillar traces around a stable mean and smooth unimodal
densities together indicate the chain has
converged.](lame-overview_files/figure-html/overview-trace-1.png)

Trace plots should bounce around a stable mean without long-term trends
or sticky regions, and density plots should be smooth and unimodal.
Numerically, common MCMC checks are **split-$`\hat R`$ \< 1.01** for
monitored parameters and **bulk / tail ESS $`\ge`$ 400 per chain**.

With 400 stored samples, the regression coefficients mix well at this
seed. The variance components mix more slowly: in the `summarise_draws`
table below, `va` carries the highest $`\hat R`$ and the lowest ESS.
Which parameter mixes worst is seed-dependent in a short chain, so treat
this fit as a compact worked example rather than a final diagnostic run.

### Bayesian-ecosystem diagnostics: `posterior::as_draws()`

`lame` registers
[`as_draws()`](https://netify-dev.github.io/lame/reference/as_draws.md)
methods so any fit drops straight into the Stan-era diagnostic ecosystem
(`posterior`, `bayesplot`, `tidybayes`). The reshape lifts the `BETA`
and `VC` blocks into a `draws_array` indexed by iteration / chain /
variable, which means
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
gives you per-parameter posterior mean, median, sd, MAD, 5%/95%
quantiles, **split-$`\hat R`$**, and **bulk / tail ESS** in one call:

``` r

library(posterior)
#> This is posterior version 1.7.0
#> 
#> Attaching package: 'posterior'
#> The following object is masked from 'package:lame':
#> 
#>     as_draws
#> The following objects are masked from 'package:stats':
#> 
#>     mad, sd, var
#> The following objects are masked from 'package:base':
#> 
#>     %in%, match
draws <- posterior::as_draws(fit)              # draws_array [iter, chain, var]
posterior::summarise_draws(draws)              # rhat, ess_bulk, ess_tail per param
#> # A tibble: 15 × 10
#>    variable            mean  median     sd    mad      q5    q95   rhat ess_bulk
#>    <chr>              <dbl>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>    <dbl>
#>  1 intercept        -2.87   -2.81   1.21   1.22   -4.88   -0.870  1.00    378.  
#>  2 male_row          0.939   0.923  0.417  0.443   0.319   1.66   1.000   462.  
#>  3 smoker_row       -0.951  -0.918  0.340  0.317  -1.53   -0.418  1.00    348.  
#>  4 program_row       0.497   0.507  0.231  0.217   0.114   0.885  0.998   359.  
#>  5 male_col          0.991   0.987  0.267  0.234   0.560   1.44   0.998   460.  
#>  6 smoker_col       -0.0620 -0.0658 0.217  0.218  -0.397   0.276  1.00    306.  
#>  7 program_col       0.149   0.158  0.145  0.143  -0.0797  0.387  1.00    368.  
#>  8 same_male_dyad    0.446   0.445  0.0583 0.0542  0.351   0.546  1.00    333.  
#>  9 same_smoker_dyad  0.143   0.147  0.0561 0.0540  0.0461  0.236  1.00     75.7 
#> 10 same_program_dy…  0.577   0.580  0.0527 0.0507  0.488   0.659  1.00    200.  
#> 11 va                0.797   0.761  0.217  0.198   0.514   1.23   1.10      8.02
#> 12 cab               0.291   0.277  0.107  0.0939  0.145   0.477  1.05     25.3 
#> 13 vb                0.268   0.255  0.0784 0.0730  0.160   0.419  1.01    201.  
#> 14 rho               0.321   0.324  0.0424 0.0441  0.250   0.385  1.000   114.  
#> 15 ve                1       1      0      0       1       1     NA        NA   
#> # ℹ 1 more variable: ess_tail <dbl>
```

A single-chain fit like this one has `chain = 1`, so $`\hat R`$ reduces
to within-chain split-$`\hat R`$, which is informative about
within-chain non-stationarity but not about between-chain disagreement.
For honest $`\hat R`$ across independently initialised chains, see the
multi-chain section in the dynamic-effects vignette.

### Goodness of Fit

The `gof_plot` function compares observed network statistics to their
posterior predictive distributions. This tells us whether the model can
reproduce key structural features of the data, not just the dyad-level
relationships, but emergent properties like how much actors vary in
their degree (heterogeneity), whether friends-of-friends tend to be
friends (transitivity), and whether ties tend to be reciprocated (dyadic
dependence).

``` r

gof_plot(fit)
```

![Longitudinal goodness-of-fit panels for the Dutch college fit: each
facet plots one network statistic (Sender Degree Heterogeneity, Receiver
Degree Heterogeneity, Dyadic Dependence, Triadic Dependence,
Transitivity) across the six time periods. The observed series is a
solid orange line with filled orange points; the posterior-predictive
median is a dashed dark line; the 95 percent credible interval is a grey
ribbon. The three elements differ on both colour and linetype, so the
comparison survives grayscale printing and colour-blind viewing.
Observed series falling outside the band flag structural features the
model is missing.](lame-overview_files/figure-html/overview-gof-1.png)

For a longitudinal model, `gof_plot` shows these statistics across time.
The figure encodes three series on each panel using both colour and
linetype redundantly: the **observed series is a solid orange line with
filled orange points**, the **posterior-predictive median is a dashed
dark line**, and the **95% credible interval is a grey ribbon** (the
width is controlled by the `credible.level` argument). When the observed
series falls inside the ribbon, the model is reproducing that aspect of
the network; when it falls outside, the model is missing something.
Falling outside is itself a substantive finding about the
data-generating process, not a failure to report.

The panel-strip headers are the human-readable labels
(`Sender Degree Heterogeneity`, etc.). `fit$GOF` is a named list whose
elements are `sd.rowmean`, `sd.colmean`, `dyad.dep`, `cycle.dep` (cyclic
triads), and `trans.dep`; user-facing aliases for
`gof_plot(..., statistics = ...)` are `sd.row`, `sd.col`, `dyad.dep`,
`triad.dep`, `trans.dep`.

Reading the figure honestly: the Dutch college fit reproduces the
average level of reciprocity and most of the cyclic-triad structure but
under-predicts sender degree heterogeneity, transitivity, and (at the
panel ends) receiver heterogeneity. The observed Sender Degree
Heterogeneity series sits above the upper edge of the 95% band at five
of the six periods (inside the band only at $`t = 1`$); some students
nominate many more friends than the additive sender effect $`a_i`$ plus
latent $`u_i`$ can absorb. Transitivity is under-predicted at every
period – the observed value sits above the upper edge of the band at
$`t = 1`$ through $`t = 6`$ – so triangle closure is systematically too
low under the latent space. Receiver Degree Heterogeneity is above the
band at $`t = 1, 2`$, inside at $`t = 3`$ and $`t = 5`$, and below the
band at $`t = 4`$ and $`t = 6`$. The Triadic Dependence panel (cyclic
triads, `cycle.dep`) sits inside the band at five of the six periods
(outside only at $`t = 1`$). Dyadic Dependence (reciprocity) is the most
variable: above the band at $`t = 1`$ and $`5`$ (and marginally at
$`t = 2`$), below at $`t = 3`$ and $`t = 6`$, inside only at $`t = 4`$.

This is the standard signature of a latent-space model on a friendship
network. AME assumes conditional dyadic independence: given the additive
and multiplicative effects, dyads are independent, so triangle closure
has to be soaked up by clustering of the latent positions $`u_i`$. ERGMs
handle triangle closure directly via terms like `gwesp` or `triangle`.
The two model classes are complementary. ERGM is the right answer for
triangle-closure questions; AME is the right answer when the analyst
wants stable, non-degenerate estimates with explicit actor heterogeneity
and a coherent posterior to forecast or simulate from.

#### Posterior-predictive temporal-trend test: `gof_temporal()`

The figure read above is qualitative. For a scalar posterior-predictive
p-value on the temporal *trend* in a chosen network statistic,
[`gof_temporal()`](https://netify-dev.github.io/lame/reference/gof_temporal.md)
complements
[`gof_plot()`](https://netify-dev.github.io/lame/reference/gof_plot.md).
It computes the chosen statistic at each observed period, fits an OLS
slope on period index, draws posterior-predictive replicate panels via
`simulate(fit)`, recomputes the slope in each replicate, and returns a
two-sided Gelman-style posterior-predictive p-value defined as
`2 * min(p_up, 1 - p_up)`, so it runs from 0 (observed slope in the
extreme tail) up to 1 (observed slope dead-centre). A static fit on
truly trending data yields `p_pp` near 0 (the observed slope is
unreachable); a fit that captures the trend yields a `p_pp` well away
from 0 (the observed slope is central in the predictive distribution).

``` r

# density rises across the panel on net (the network gets denser as
# students settle in, though the trajectory is not monotone), so it is the
# natural temporal-trend target.
gt_density <- gof_temporal(fit, stat = "density", n_rep = 100, seed = 6886)
gt_density   # prints stat, observed slope, n_rep, and the p_pp
#> 
#> ── Temporal-trend posterior-predictive check ──
#> 
#> • Statistic: "density"
#> • Observed slope (per-period): 0.025
#> • Replicates: 100
#> • Posterior-predictive p-value (two-sided): 0
#> Observed temporal trend is incompatible with the fitted model.

# reciprocity is the natural second target for a directed network.
gt_recip <- gof_temporal(fit, stat = "reciprocity", n_rep = 100, seed = 6886)
gt_recip
#> 
#> ── Temporal-trend posterior-predictive check ──
#> 
#> • Statistic: "reciprocity"
#> • Observed slope (per-period): -0.0306
#> • Replicates: 100
#> • Posterior-predictive p-value (two-sided): 0
#> Observed temporal trend is incompatible with the fitted model.
```

A `p_pp` of 0 on `density` says the static fit cannot reproduce the
empirical net densification slope (about +0.025 per period in our run);
that is the expected diagnostic for a fit without `dynamic_beta` or
time-varying covariates and is consistent with the per-period mass shift
in
[`gof_plot()`](https://netify-dev.github.io/lame/reference/gof_plot.md).
The reciprocity check tells the same story in the opposite direction:
the static fit also rejects the observed reciprocity slope (`p_pp` near
0), because reciprocity actually drifts *downward* across the panel
(slope about $`-0.03`$ per period) and the static $`\rho`$ cannot track
that drift. Both checks are pointing at the same structural mismatch:
model parameters that are constant across periods cannot reproduce
period-to-period drift in either marginal quantity. The companion
diagnostic
[`detect_change_point()`](https://netify-dev.github.io/lame/reference/detect_change_point.md)
reports a heuristic tail-ratio score comparing the posterior path of
each `dynamic_beta` coefficient to its AR(1) prior null. It requires a
fit with `dynamic_beta` active and is therefore not applicable to the
static fit shown here; the worked example lives in the [Dynamic
Effects](https://netify-dev.github.io/lame/articles/dynamic_effects.md)
vignette.

### Latent Space

The multiplicative effects capture patterns of association that go
beyond what the covariates and additive effects explain. The `uv_plot`
function visualizes these latent positions. Because we fit `R = 1`, the
latent space is one-dimensional: each student’s sender position
(triangle) and receiver position (circle) is a single number, and
`uv_plot` lays them out along one axis (the two reference circles are
drawn only for scale). Students who land near each other share
friendship-formation patterns the covariates and additive effects miss;
sign separates the two latent groups. With `R = 2` the same call would
spread the positions across a genuine two-dimensional plane.

``` r

uv_plot(fit) +
    ggtitle("Dutch college friendship: multiplicative effects") +
    theme_bw() +
    theme(
        panel.border = element_blank(),
        legend.position = "none",
        axis.text  = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
    )
```

![One-dimensional latent-space layout for the Dutch college fit (R = 1):
triangles mark each student's sender position and circles the receiver
position, both falling along a single horizontal axis because the latent
space has one dimension, with two reference circles drawn only for scale
and marker size proportional to magnitude; students near each other
share friendship-formation patterns beyond what the covariates and
additive effects
explain.](lame-overview_files/figure-html/overview-uv-1.png)

With `R = 1` the latent space is one-dimensional, so the positions fall
along a single axis and the two circles serve only as a reference scale;
marker size encodes magnitude. Students at the extremes – and those
whose sender (triangle) and receiver (circle) positions diverge – have
the most distinctive friendship-formation patterns that the covariates
and additive effects do not explain on their own. With `R = 2` the same
call spreads the positions across a genuine plane.

## Extracting Latent Positions

If you need the estimated latent positions in a tidy format (for custom
plots, merging with external data, or exporting to other tools), the
[`latent_positions()`](https://netify-dev.github.io/lame/reference/latent_positions.md)
function returns a data frame:

``` r

lp <- latent_positions(fit)
#> ℹ `posterior_sd` is "NA" because U/V samples were not saved.
#> ℹ To get posterior SDs, refit with `posterior_opts = posterior_options(save_UV
#>   = TRUE)`.
#> This message is displayed once per session.
head(lp)
#>   actor dimension time      value posterior_sd type
#> 1   S01         1    1 -1.3214083           NA    U
#> 2   S02         1    1 -0.1260811           NA    U
#> 3   S03         1    1  0.2875313           NA    U
#> 4   S04         1    1  0.3985233           NA    U
#> 5   S05         1    1  0.5024891           NA    U
#> 6   S06         1    1  0.8535405           NA    U
```

Each row gives one actor’s position on one latent dimension at one time
point; a `type` column distinguishes the sender position (`"U"`) from
the receiver position (`"V"`), so a directed fit returns two rows per
actor per dimension. For a static fit like this one the positions do not
change over time, so
[`latent_positions()`](https://netify-dev.github.io/lame/reference/latent_positions.md)
returns a single slice (labelled `time = 1`) that applies to every
period rather than repeating it. For dynamic models
(`dynamic_uv = TRUE`), positions evolve over time and
[`latent_positions()`](https://netify-dev.github.io/lame/reference/latent_positions.md)
returns each time point separately (see the [dynamic effects
vignette](https://netify-dev.github.io/lame/articles/dynamic_effects.md)).

When using dynamic latent positions, the latent space is only identified
up to rotation at each time point. The
[`procrustes_align()`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
function aligns each time period’s positions to the previous one,
removing arbitrary rotations so that trajectories across time are
interpretable:

``` r

# procrustes alignment is useful for dynamic models.
# for a static model, positions are already aligned.
aligned <- procrustes_align(fit)
#> ℹ Latent positions are static (single time point). No alignment needed.
str(aligned$U)
#>  num [1:32, 1] -1.321 -0.126 0.288 0.399 0.502 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:32] "S01" "S02" "S03" "S04" ...
#>   ..$ : NULL
```

## Tidy-data round-trip via netify

`lame` integrates with the `netify` package via the
`netify::to_amen(netlet, lame = TRUE)` bridge: feed a long-format
edgelist (or any tidy network representation) through `netify::netify()`
and then pass the result of `to_amen()` straight into
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) or
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md). The
`lame = TRUE` flag preserves the per-period dyadic-covariate axis
(required for longitudinal models). Once the fit is built,
[`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md),
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html),
and
[`prediction_draws_long()`](https://netify-dev.github.io/lame/reference/prediction_draws_long.md)
give you broom-/ggplot-/marginaleffects-ready output that composes
naturally with the rest of the tidyverse.

### Entry: tidy edgelist → `lame()`

``` r

library(netify)

# 1. tidy edgelist -> netify
netlet <- netify::netify(
    edges_long,
    actor1 = "from", actor2 = "to", time = "year",
    weight = "tie", dyad_vars = "distance",
    symmetric = FALSE, mode = "unipartite",
    output_format = "longit_array"
)

# 2. netify -> lame-ready list (lame = TRUE for longitudinal)
amen_in <- netify::to_amen(netlet, lame = TRUE)

# 3. fit (same Dutch-college style fit we built above)
fit <- lame(
    Y      = amen_in$Y,
    Xdyad  = amen_in$Xdyad,
    Xrow   = amen_in$Xrow,
    Xcol   = amen_in$Xcol,
    family = "binary", R = 1,
    burn = 500, nscan = 4000, odens = 10,
    verbose = FALSE
)
```

### Exit: `tidy()`, `autoplot()`, `prediction_draws_long()`

The rest of the round-trip runs against the `fit` object built from the
Dutch-college data above (no `netify` install needed):

``` r

# tidy(fit) -> term / estimate / std.error / statistic / p.value /
# conf.low / conf.high; one row per coefficient (or per
# coefficient x period for dynamic_beta fits). The generic ships with
# lame so this runs without broom installed; when broom is on the path
# `broom::tidy(fit)` dispatches to the same method via
# generics::tidy registration.
tidy(fit)
#> # A tibble: 10 × 7
#>    term              estimate std.error statistic  p.value conf.low conf.high
#>    <chr>                <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#>  1 intercept          -2.87      1.21      -2.38  1.75e- 2  -5.24      -0.572
#>  2 male_row            0.939     0.417      2.25  2.44e- 2   0.212      1.76 
#>  3 smoker_row         -0.951     0.340     -2.80  5.12e- 3  -1.58      -0.311
#>  4 program_row         0.497     0.231      2.15  3.13e- 2   0.0592     0.947
#>  5 male_col            0.991     0.267      3.71  2.09e- 4   0.446      1.52 
#>  6 smoker_col         -0.0620    0.217     -0.286 7.75e- 1  -0.525      0.341
#>  7 program_col         0.149     0.145      1.03  3.04e- 1  -0.111      0.437
#>  8 same_male_dyad      0.446     0.0583     7.65  1.95e-14   0.333      0.566
#>  9 same_smoker_dyad    0.143     0.0561     2.55  1.07e- 2   0.0369     0.244
#> 10 same_program_dyad   0.577     0.0527    10.9   0          0.472      0.669
```

### Coefficient tables: `glance()` + `tidy()` -\> `modelsummary` / `gt` / `kableExtra`

`lame` ships both
[`tidy.lame()`](https://netify-dev.github.io/lame/reference/tidy.ame.md)
and
[`glance.lame()`](https://netify-dev.github.io/lame/reference/glance.ame.md),
registered against
[`generics::tidy`](https://generics.r-lib.org/reference/tidy.html) /
[`generics::glance`](https://generics.r-lib.org/reference/glance.html)
at load time, so any broom-aware table package consumes a `lame` fit
directly. The most compact route is `modelsummary::modelsummary(fit)`,
which uses
[`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md) for the
coefficient panel and
[`glance()`](https://netify-dev.github.io/lame/reference/glance.md) for
the goodness-of-fit footer:

``` r

# glance() returns the one-row model summary modelsummary uses for its
# lower panel: nobs, n_actors, n_periods, n_stored, family, mode, R,
# dynamic_uv / dynamic_ab / dynamic_beta, elpd_loo (NA unless save_log_lik=TRUE).
broom::glance(fit)
#> # A tibble: 1 × 13
#>    nobs n_actors n_row_actors n_col_actors n_periods n_stored family mode      R
#>   <int>    <int>        <int>        <int>     <int>    <int> <chr>  <chr> <int>
#> 1  5952       32           NA           NA         6      400 binary unip…     1
#> # ℹ 4 more variables: dynamic_uv <lgl>, dynamic_ab <lgl>, dynamic_beta <lgl>,
#> #   elpd_loo <dbl>

# pass the fit straight to modelsummary -- tidy() drives the upper
# (coefficients) panel, glance() drives the lower (GOF) panel.
modelsummary::modelsummary(
    list("Dutch college" = fit),
    statistic = "conf.int",
    gof_map = c("nobs", "n_actors", "n_periods", "n_stored",
                "family", "R", "dynamic_uv", "dynamic_ab")
)
```

|                   | Dutch college      |
|-------------------|--------------------|
| intercept         | -2.869             |
|                   | \[-5.241, -0.572\] |
| male_row          | 0.939              |
|                   | \[0.212, 1.761\]   |
| smoker_row        | -0.951             |
|                   | \[-1.583, -0.311\] |
| program_row       | 0.497              |
|                   | \[0.059, 0.947\]   |
| male_col          | 0.991              |
|                   | \[0.446, 1.516\]   |
| smoker_col        | -0.062             |
|                   | \[-0.525, 0.341\]  |
| program_col       | 0.149              |
|                   | \[-0.111, 0.437\]  |
| same_male_dyad    | 0.446              |
|                   | \[0.333, 0.566\]   |
| same_smoker_dyad  | 0.143              |
|                   | \[0.037, 0.244\]   |
| same_program_dyad | 0.577              |
|                   | \[0.472, 0.669\]   |
| Num.Obs.          | 5952               |
| n_actors          | 32                 |
| n_periods         | 6                  |
| n_stored          | 400                |
| family            | binary             |
| R                 | 1                  |
| dynamic_uv        | FALSE              |
| dynamic_ab        | FALSE              |

Multiple fits compose into a side-by-side table the same way they would
for [`lm()`](https://rdrr.io/r/stats/lm.html) /
[`glm()`](https://rdrr.io/r/stats/glm.html) – pass a named list, one
column per fit. Below we compare the full specification against a
no-latent-space variant (`R = 0`) to see how much of the coefficient
story is being soaked up by multiplicative effects:

``` r

fit_noUV <- lame(
    Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
    family = "binary", rvar = TRUE, cvar = TRUE, dcor = TRUE,
    R = 0,                          # no latent space
    symmetric = FALSE,
    burn = 500, nscan = 4000, odens = 10,
    verbose = FALSE, plot = FALSE
)

modelsummary::modelsummary(
    list("AME (R = 1)" = fit, "Additive only (R = 0)" = fit_noUV),
    statistic = "conf.int",
    gof_map = c("nobs", "n_actors", "n_periods", "n_stored",
                "family", "R", "dynamic_uv", "dynamic_ab")
)
```

|                   | AME (R = 1)        | Additive only (R = 0) |
|-------------------|--------------------|-----------------------|
| intercept         | -2.869             | -2.693                |
|                   | \[-5.241, -0.572\] | \[-5.172, -0.408\]    |
| male_row          | 0.939              | 0.915                 |
|                   | \[0.212, 1.761\]   | \[0.133, 1.698\]      |
| smoker_row        | -0.951             | -0.853                |
|                   | \[-1.583, -0.311\] | \[-1.508, -0.199\]    |
| program_row       | 0.497              | 0.498                 |
|                   | \[0.059, 0.947\]   | \[0.040, 0.981\]      |
| male_col          | 0.991              | 0.973                 |
|                   | \[0.446, 1.516\]   | \[0.501, 1.457\]      |
| smoker_col        | -0.062             | -0.094                |
|                   | \[-0.525, 0.341\]  | \[-0.499, 0.347\]     |
| program_col       | 0.149              | 0.104                 |
|                   | \[-0.111, 0.437\]  | \[-0.161, 0.387\]     |
| same_male_dyad    | 0.446              | 0.412                 |
|                   | \[0.333, 0.566\]   | \[0.298, 0.521\]      |
| same_smoker_dyad  | 0.143              | 0.191                 |
|                   | \[0.037, 0.244\]   | \[0.108, 0.275\]      |
| same_program_dyad | 0.577              | 0.594                 |
|                   | \[0.472, 0.669\]   | \[0.501, 0.684\]      |
| Num.Obs.          | 5952               | 5952                  |
| n_actors          | 32                 | 32                    |
| n_periods         | 6                  | 6                     |
| n_stored          | 400                | 400                   |
| family            | binary             | binary                |
| R                 | 1                  | 0                     |
| dynamic_uv        | FALSE              | FALSE                 |
| dynamic_ab        | FALSE              | FALSE                 |

Reading down the `R = 1` column versus the `R = 0` column shows how the
latent space can redistribute signal: when residual clustering is real,
homophily coefficients move toward zero in the AME column because
structure that was loaded onto `same_program` (etc.) is instead absorbed
by the latent positions. On this small panel (32 students, short chains)
the movement is modest and not uniform across the three homophily terms,
with widely overlapping credible intervals, so treat the comparison as
illustrating the mechanism rather than as a precise decomposition; the
direction is what to look for, not the exact shift.

**Side-by-side MCMC and ALS in one table.**
[`tidy.ame_als()`](https://netify-dev.github.io/lame/reference/tidy.ame_als.md)
and
[`glance.ame_als()`](https://netify-dev.github.io/lame/reference/glance.ame_als.md)
are registered alongside the MCMC methods, so
`modelsummary(list("MCMC" = fit, "ALS" = fit_als))` works directly. The
ALS row reports sandwich-based SEs when no bootstrap is attached and
bootstrap-based SEs when `ame_als(..., bootstrap = K)` was used; the
`se_source` column on `tidy(fit_als)` records which one supplied the
interval.

If you prefer `gt` directly, the same
[`tidy()`](https://netify-dev.github.io/lame/reference/tidy.md) frame
composes into any formatter; every column follows the broom contract so
`kableExtra::kbl()`, `flextable::flextable()`, and `gt::gt()` all work
without a custom shim:

``` r

library(dplyr)
library(gt)

broom::tidy(fit) |>
    mutate(`Estimate (95% CI)` = sprintf("%.3f (%.3f, %.3f)",
                                         estimate, conf.low, conf.high)) |>
    select(term, `Estimate (95% CI)`, p.value) |>
    gt::gt() |>
    gt::fmt_number(columns = p.value, decimals = 3) |>
    gt::tab_header(title = "Posterior coefficient summary",
                   subtitle = "lame::lame() fit, Dutch college panel")
```

``` r

# autoplot.lame returns a horizontal coefplot for static fits and a
# ribbon-per-period for dynamic_beta fits, so the same call works for
# both. ggplot layers compose on top.
autoplot(fit) +
    ggtitle("Dutch college friendship: posterior coefficients")
```

![Horizontal coefficient plot returned by autoplot.lame: each row is a
regression coefficient, with a point at the posterior mean and a
horizontal segment marking the 95 percent credible interval; a dashed
reference line at zero indicates where an effect is indistinguishable
from null.](lame-overview_files/figure-html/autoplot-fit-1.png)

The same generic dispatches on `which = "uv"` (latent positions) and
`which = "ab"` (sender / receiver random effects), so the entire
random-effects diagnostic suite is reachable from a single broom-style
verb:

``` r

# which = "uv": delegates to uv_plot(); trajectory for dynamic_uv,
# snapshot for static -- returns a ggplot in both cases.
autoplot(fit, which = "uv")
```

![Two separate ggplots returned by autoplot.lame, rendered one above the
other: a one-dimensional latent-space layout of posterior-mean U/V
positions, and a lollipop plot of sender random
effects.](lame-overview_files/figure-html/autoplot-uv-ab-1.png)

``` r


# which = "ab": delegates to ab_plot() (sender side by default);
# call ab_plot(fit, effect = "receiver") for the column-side plot.
autoplot(fit, which = "ab")
```

![Two separate ggplots returned by autoplot.lame, rendered one above the
other: a one-dimensional latent-space layout of posterior-mean U/V
positions, and a lollipop plot of sender random
effects.](lame-overview_files/figure-html/autoplot-uv-ab-2.png)

``` r

# prediction_draws_long() returns a long-format data frame with
# .chain / .iteration / .draw / period / period_label / i / j /
# actor_i / actor_j / .value -- ready for tidybayes / marginaleffects.
pdl <- prediction_draws_long(fit, type = "response", n_draws = 50)
# the diagonal (i == j) is a self-tie, undefined in a friendship network, but
# the linear predictor is still numerically defined there -- drop it so the
# per-dyad summaries below cover only real (off-diagonal) pairs.
pdl <- pdl[pdl$i != pdl$j, ]
head(pdl, 3)
#> # A tibble: 3 × 10
#>   .chain .iteration .draw period period_label     i     j actor_i actor_j
#>    <int>      <int> <int>  <int> <chr>        <int> <int> <chr>   <chr>  
#> 1      1          1     1      1 t2               2     1 S02     S01    
#> 2      1          1     1      1 t2               3     1 S03     S01    
#> 3      1          1     1      1 t2               4     1 S04     S01    
#> # ℹ 1 more variable: .value <dbl>
dim(pdl)
#> [1] 297600     10
```

### Worked tidybayes-style summary

Because
[`prediction_draws_long()`](https://netify-dev.github.io/lame/reference/prediction_draws_long.md)
already uses the `.draw` / `.value` column convention, you can pipe
straight into a `tidybayes::spread_draws()`- style workflow, or, since
the columns are already in the right shape, go directly to `dplyr` for
per-dyad posterior summaries:

``` r

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
pdl |>
    group_by(actor_i, actor_j, period_label) |>
    summarise(p_hat = mean(.value),
              lo    = quantile(.value, 0.025),
              hi    = quantile(.value, 0.975),
              .groups = "drop") |>
    head(3)
#> # A tibble: 3 × 6
#>   actor_i actor_j period_label  p_hat            lo    hi
#>   <chr>   <chr>   <chr>         <dbl>         <dbl> <dbl>
#> 1 S01     S02     t2           0.0199 0.00000000237 0.138
#> 2 S01     S02     t3           0.0199 0.00000000237 0.138
#> 3 S01     S02     t4           0.0199 0.00000000237 0.138
```

Because the `.chain` / `.iteration` / `.draw` / `.value` columns already
follow the `tidybayes` / `marginaleffects` convention, the data frame
composes directly with `ggdist::stat_pointinterval()`,
`tidybayes::point_interval()`, and `marginaleffects::posterior_draws()`
without an intermediate `spread_draws()` call (that generic dispatches
on draws / brms / Stan model objects, not on a long data frame – the
point of returning the long shape directly is to skip that step).

## Model comparison: `loo::loo_compare()` as a one-liner

When two specifications are on the table – say a static-coefficient
baseline against a `dynamic_beta = "dyad"` variant – the standard
out-of-sample yardstick is the expected log pointwise predictive density
(`elpd_loo`) returned by
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html). Both fits
need `save_log_lik = TRUE` so
[`loo()`](https://netify-dev.github.io/lame/reference/loo.md) can read
the pointwise log-likelihood matrix off `fit$log_lik`;
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
then ranks them in a single call:

``` r

fit_static <- lame(Y, Xdyad = Xdyad, family = "binary", R = 2,
                   burn = 500, nscan = 4000, odens = 10,
                   save_log_lik = TRUE, verbose = FALSE)
fit_dyn    <- lame(Y, Xdyad = Xdyad, family = "binary", R = 2,
                   dynamic_beta = "dyad",
                   burn = 500, nscan = 4000, odens = 10,
                   save_log_lik = TRUE, verbose = FALSE)

cmp <- loo::loo_compare(loo::loo(fit_static), loo::loo(fit_dyn))
cmp
#>            elpd_diff se_diff
#> fit_dyn     0.0       0.0
#> fit_static -8.4       3.1
```

**How to read the output.** The top row is the preferred model and
always carries `elpd_diff = 0`. Each subsequent row reports the
*difference* in expected log pointwise predictive density relative to
the top model, together with the standard error of that difference. In
the printout above, `fit_dyn` beats `fit_static` by 8.4 elpd units with
an SE of 3.1, so $`|\text{elpd\_diff}| / \text{se\_diff} \approx 2.7`$ –
the dynamic-beta specification is preferred at roughly the “2-sigma”
threshold. A natural way to phrase that comparison:

> *“The time-varying-coefficient model is preferred over the static
> specification by 8.4 elpd units (SE 3.1) on leave-one-out
> cross-validation, a `|`elpd_diff`|`/SE ratio of 2.7. This is
> consistent with non-trivial drift in the dyadic coefficient over the
> panel window; the static specification is rejected at conventional
> thresholds.”*

If `|elpd_diff|` were closer to its SE (say `-2.0 (SE 3.0)`), the two
specifications would be predictively indistinguishable and the simpler
static model would be the defensible choice. The numbers in the printout
above are illustrative; the preferred model and the gap size depend on
the data, the chain length, and how well AME captures the structural
features of the network. A
[`loo()`](https://netify-dev.github.io/lame/reference/loo.md) run on a
short Gibbs chain (the kind a replicator typically uses while iterating
on a specification) will also commonly flag a few percent of dyads with
Pareto $`\hat k > 0.7`$; treat that warning as a sign that the
importance-sampling approximation is straining and either lengthen the
chain or apply
[`loo::loo_moment_match()`](https://mc-stan.org/loo/reference/loo_moment_match.html)
before relying on the `elpd_diff` for inference.

Because `family = "binary"` is one of the six families with an exact
closed-form pointwise log-likelihood (`normal`, `binary`, `cbin`,
`tobit`, `poisson`, `ordinal`), the absolute `elpd_loo` values here are
directly comparable to a
[`loo()`](https://netify-dev.github.io/lame/reference/loo.md) result
from a probit-link `brms` fit on the same data, modulo $`\hat R`$ and
ESS being satisfactory on the `lame` Gibbs sampler. For `frn` and `rrl`,
see the `log_lik_method = "observed_ghk"` caveat in [Dynamic
Effects](https://netify-dev.github.io/lame/articles/dynamic_effects.html#loo-waic-via-save_log_lik-true).

## Reproducibility: `sessionInfo()` and `renv.lock`

MCMC posteriors depend on the RNG state, the C++ ABI of `RcppArmadillo`,
and minor internal changes in `lame` itself across releases. A `lame`
analysis that a reviewer or your future self can reproduce captures both
the seed and the exact package versions you fit under. The conventional
pattern at the top of any replication script:

``` r

set.seed(6886)
fit <- lame(...)
sessionInfo()                 # paste into supplementary materials
```

For long-running fits, save the fitted object together with the session
record:

``` r

saveRDS(list(fit = fit, session = sessionInfo()),
        file = "lame_fit_replication.rds")
```

For long-term reproducibility, pin `lame` to a specific commit hash (the
package is pre-CRAN at the time of writing):

``` r

remotes::install_github("netify-dev/lame@<commit-sha>")
```

or commit an `renv.lock` to your project repo so reviewers can recreate
the entire toolchain with `renv::restore()`. Random-seed + package
version is the floor; commit hash + `renv.lock` is the ceiling.

``` r

# this vignette's own session record -- copy the pattern verbatim
# into your replication archive
sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.2.1     posterior_1.7.0 ggplot2_4.0.3   lame_1.2.0     
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6          tensorA_0.36.2.1      bayestestR_0.18.1    
#>  [4] xfun_0.58             bslib_0.11.0          insight_1.5.1        
#>  [7] ggrepel_0.9.8         lattice_0.22-9        vctrs_0.7.3          
#> [10] tools_4.6.0           generics_0.1.4        datawizard_1.3.1     
#> [13] parallel_4.6.0        tibble_3.3.1          pkgconfig_2.0.3      
#> [16] tinytable_0.16.0      data.table_1.18.4     checkmate_2.3.4      
#> [19] RColorBrewer_1.1-3    S7_0.2.2              desc_1.4.3           
#> [22] distributional_0.7.0  lifecycle_1.0.5       compiler_4.6.0       
#> [25] farver_2.1.2          stringr_1.6.0         textshaping_1.0.5    
#> [28] ggforce_0.5.0         htmltools_0.5.9       sass_0.4.10          
#> [31] yaml_2.3.12           pillar_1.11.1         pkgdown_2.2.0        
#> [34] jquerylib_0.1.4       tidyr_1.3.2           MASS_7.3-65          
#> [37] cachem_1.1.0          abind_1.4-8           network_1.20.0       
#> [40] tidyselect_1.2.1      digest_0.6.39         performance_0.17.0   
#> [43] stringi_1.8.7         reshape2_1.4.5        purrr_1.2.2          
#> [46] labeling_0.4.3        polyclip_1.10-7       fastmap_1.2.0        
#> [49] grid_4.6.0            cli_3.6.6             magrittr_2.0.5       
#> [52] patchwork_1.3.2       loo_2.9.0             utf8_1.2.6           
#> [55] broom_1.0.13          withr_3.0.2           scales_1.4.0         
#> [58] backports_1.5.1       rmarkdown_2.31        matrixStats_1.5.0    
#> [61] modelsummary_2.6.0    ragg_1.5.2            coda_0.19-4.1        
#> [64] evaluate_1.0.5        knitr_1.51            parameters_0.29.1    
#> [67] rlang_1.2.0           Rcpp_1.1.1-1.1        glue_1.8.1           
#> [70] tweenr_2.0.3          jsonlite_2.0.0        R6_2.6.1             
#> [73] plyr_1.8.9            statnet.common_4.13.0 tables_0.9.33        
#> [76] systemfonts_1.3.2     fs_2.1.0
```

## What’s Next?

This vignette covered the core workflow: fitting, convergence, GOF, and
visualization. For more specialized topics:

- **Single networks (no time series)?** See [Your First AME
  Model](https://netify-dev.github.io/lame/articles/cross_sec_ame.md)
  for a detailed cross-sectional walkthrough
- **Two types of nodes?** See [Bipartite
  Networks](https://netify-dev.github.io/lame/articles/bipartite.md)
- **Evolving network structure?** See [Dynamic
  Effects](https://netify-dev.github.io/lame/articles/dynamic_effects.md)
- **Just want the quick version?** See [Getting
  Started](https://netify-dev.github.io/lame/articles/lame.md)

**References**:

Hoff, PD (2021) Additive and Multiplicative Effects Network Models.
Statistical Science 36, 34–50.

van de Bunt, G. G., van Duijn, M. A. J., & Snijders, T. A. B. (1999).
Friendship networks through time: An actor-oriented dynamic statistical
network model. *Computational & Mathematical Organization Theory*, 5(2),
167–192.

Minhas, S., Dorff, C., Gallop, M. B., Foster, M., Liu, H., Tellez, J., &
Ward, M. D. (2022). Taking dyads seriously. Political Science Research
and Methods, 10(4), 703–721. <doi:10.1017/psrm.2021.56>
