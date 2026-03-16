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
  (`trace_plot`, `gof_plot`, `ab_plot`, `uv_plot`) produce
  publication-ready ggplot2 figures.
- **Standard S3 methods**:
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html), and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) work as you
  would expect from any R model object.

The package is part of the `netify-verse`, so it works seamlessly with
the [`netify`](https://netify-dev.github.io/netify/) package for data
preparation.

## Application: International Sanctions

To see how `lame` works in practice, we analyze data from the Threat and
Imposition of Sanctions (TIES) dataset (Morgan et al. 2014). The data
cover 35 countries across four years (1993–1995 and 2000). The core
question is straightforward: *which country-level and dyad-level factors
predict whether one country sanctions another, and are there latent
patterns of sanctioning behavior beyond what these covariates explain?*

The dataset includes:

- **Y**: A list of directed binary adjacency matrices (one per year).
  Entry $y_{ij,t} = 1$ means country $i$ imposed sanctions on country
  $j$ in year $t$.
- **Xdyad**: Dyadic covariates (geographic distance and the number of
  shared intergovernmental organizations (IGOs)) that capture
  relationship-level predictors.
- **Xrow / Xcol**: Sender and receiver covariates (GDP and population)
  capturing country-level characteristics. We use the same variables for
  both senders and receivers here, but different covariates can be
  specified for each role.

### Fitting the Model

We fit a binary probit AME model with sender random effects (`rvar`),
receiver random effects (`cvar`), dyadic correlation (`dcor`), and a
2-dimensional latent space (`R = 2`). The random effects capture the
fact that some countries are simply more active sanctioners (sender
effects) or more frequent targets (receiver effects), while the latent
space picks up residual patterns, groups of countries that sanction
similar targets or are targeted by the same senders.

``` r
# Fit the AME model for binary data
# Note: burn and nscan are kept small for fast vignette building.
# For real analyses, use burn >= 1000 and nscan >= 5000.
fit <- lame(
  Y = Y,
  Xdyad = Xdyad,           # dyadic covariates
  Xrow = Xrow,             # sender covariates
  Xcol = Xcol,             # receiver covariates
  family = "binary",       # Binary probit model
  rvar = TRUE,             # sender random effects
  cvar = TRUE,             # receiver random effects
  dcor = TRUE,             # Dyadic correlation
  R = 2,                   # Multiplicative effects dimension
  symmetric = FALSE,       # Directed network
  burn = 100,              # Burn-in iterations
  nscan = 500,             # Post-burn-in iterations
  odens = 25,              # Output density (thinning)
  verbose = FALSE,           # Suppress iteration output
  plot = FALSE             # Suppress real-time plots
)
```

### Interpreting the Results

``` r
summary(fit)

=== Longitudinal AME Model Summary ===

Call:
NULL

Time periods: 4 
Family: binary 
Mode: unipartite 

Regression coefficients:
------------------------
                 Estimate StdError z_value p_value CI_lower CI_upper  
intercept          -5.534    2.328  -2.377   0.017  -10.097    -0.97 *
log_gdp.row        -0.152    0.226  -0.672   0.502   -0.595    0.291  
log_pop.row         0.161    0.291   0.552   0.581    -0.41    0.732  
log_gdp.col        -0.029     0.11  -0.262   0.793   -0.245    0.187  
log_pop.col         -0.06    0.151  -0.399    0.69   -0.357    0.236  
distance.dyad        0.04     0.04   0.986   0.324   -0.039    0.118  
shared_igos.dyad    0.017    0.008   2.249   0.024    0.002    0.032 *
---
Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Variance components:
-------------------
    Estimate StdError
va     6.170    2.126
cab    2.414    1.066
vb     1.002    0.550
rho    0.417    0.138
ve     1.000    0.000
  (va = sender, cab = sender-receiver covariance, vb = receiver,
   rho = dyadic correlation, ve = residual variance)
```

The regression coefficients tell us about the structural drivers of
sanctions. Shared IGO membership is positively associated with
sanctions, as countries embedded in the same international institutions
are more likely to sanction each other, probably because these
institutions provide both the information and the institutional
mechanisms to coordinate economic pressure. Distance has a negative
coefficient, meaning geographically closer countries are more likely to
sanction one another.

The variance components at the bottom of the summary capture
network-level heterogeneity. Sender variance (`va`) and receiver
variance (`vb`) indicate how much countries differ in their baseline
propensity to impose or receive sanctions. The dyadic correlation
(`rho`) captures reciprocity, the tendency for sanctions to be mutual.

### Checking Convergence

Because `lame` uses MCMC, we need to verify that the sampler has
converged before trusting the results. The `trace_plot` function shows
the sampled values over iterations (left panels) and the corresponding
posterior densities (right panels) for each parameter.

``` r
trace_plot(fit, params = "beta")
```

![MCMC trace plots and density plots showing convergence and posterior
distributions of regression
coefficients](lame-overview_files/figure-html/unnamed-chunk-4-1.png)

What to look for: the trace plots should look like “fuzzy caterpillars,”
bouncing around a stable mean without long-term trends or getting stuck
in one region. The density plots should be smooth and unimodal. With
only 20 post-burn-in samples (due to the small `nscan` we used here),
the chains have not fully converged. In a real analysis, you would run
much longer chains and check that increasing the run length doesn’t
substantially change the posterior summaries.

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

![Goodness of fit plots comparing observed network statistics to
posterior predictive distributions across time
periods](lame-overview_files/figure-html/unnamed-chunk-5-1.png)

For a longitudinal model, `gof_plot` shows these statistics across time.
The observed values (points) are overlaid on the posterior predictive
intervals (shaded bands). When the observed values fall inside the
bands, the model is capturing that aspect of the network well. When they
fall outside, it signals that the model is missing something. For
instance, if the model consistently under-predicts transitivity, it
might suggest a higher-dimensional latent space is needed.

### Latent Space

The multiplicative effects capture patterns of association that go
beyond what the covariates and additive effects explain. The `uv_plot`
function visualizes these latent positions. Countries positioned near
each other in the sender space (triangles) tend to sanction similar
targets; countries near each other in the receiver space (circles) tend
to be sanctioned by the same senders. Countries placed on opposite sides
of the plot tend to have dissimilar sanctioning patterns.

``` r

# Network plot showing multiplicative effects
uv_plot(fit) +
  ggtitle("Sanctions Network - Multiplicative Effects") +
  theme_minimal() +
  theme(
    legend.position='none',
    axis.text=element_blank(),
    axis.title=element_blank()
  )
```

![Network visualization showing multiplicative effects latent space
positions of countries in the sanctions
network](lame-overview_files/figure-html/unnamed-chunk-6-1.png)

The circular layout places sender positions (triangles) on the outer
ring and receiver positions (circles) on an inner ring, scaled by
magnitude. Countries with larger effects are placed further from the
center, indicating they play a more distinctive role in the network.
Countries clustered together (for example, Western allies like the USA
and Canada) share similar sanctioning profiles.

## Extracting Latent Positions

If you need the estimated latent positions in a tidy format (for custom
plots, merging with external data, or exporting to other tools), the
[`latent_positions()`](https://netify-dev.github.io/lame/reference/latent_positions.md)
function returns a data frame:

``` r
lp <- latent_positions(fit)
head(lp)
  actor dimension time      value posterior_sd type
1   AFG         1 <NA> -0.5253670           NA    U
2   ALB         1 <NA> -1.0708329           NA    U
3   ARG         1 <NA>  0.3959614           NA    U
4   AUS         1 <NA>  1.2064454           NA    U
5   BEL         1 <NA> -0.3992633           NA    U
6   BOL         1 <NA> -2.1687274           NA    U
```

Each row gives one actor’s position on one latent dimension at one time
point. For dynamic models (`dynamic_uv = TRUE`), the `time` column
tracks how positions evolve. You can optionally apply Procrustes
alignment to remove rotational indeterminacy across time periods:

``` r
aligned <- procrustes_align(fit)
 [36mℹ [39m Latent positions are static (single time point). No alignment needed.
str(aligned$U)  # aligned 3D array
 num [1:35, 1:2] -0.525 -1.071 0.396 1.206 -0.399 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:35] "AFG" "ALB" "ARG" "AUS" ...
  ..$ : NULL
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

Morgan, T. Clifton, Bapat, N., & Kobayashi, Y. (2014). Threat and
Imposition of Economic Sanctions 1945-2005. Conflict Management and
Peace Science 31(5): 541-558.

Minhas, S., Dorff, C., Gallop, M. B., Foster, M., Liu, H., Tellez, J., &
Ward, M. D. (2022). Taking dyads seriously. Political Science Research
and Methods, 10(4), 703–721. <doi:10.1017/psrm.2021.56>
