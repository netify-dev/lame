# Getting Started with lame

## What Is lame?

The `lame` package fits **L**ongitudinal **A**dditive and
**M**ultiplicative **E**ffects models to network data. If you have a
network (who trades with whom, who is friends with whom, who sanctions
whom) observed at one or more time points, `lame` lets you estimate
covariate effects while properly accounting for the dependencies that
make network data special.

The core idea: every actor has a tendency to send ties ($a_{i}$),
receive ties ($b_{j}$), and a latent position ($u_{i}$, $v_{j}$) that
captures who connects with whom. Covariates shift the baseline
probability of ties. And all of this can optionally evolve over time.

## A 5-Minute Example

Let’s fit a model to a small longitudinal binary network in under 20
lines of code.

``` r
library(lame)
set.seed(6886)

# Simulate 3 time periods of a 15-node directed network
n <- 15; T_periods <- 3
Y_list <- lapply(1:T_periods, function(t) {
  Y <- matrix(rbinom(n * n, 1, 0.15), n, n)
  diag(Y) <- NA
  rownames(Y) <- colnames(Y) <- paste0("N", 1:n)
  Y
})

# Fit a longitudinal AME model
fit <- lame(
  Y = Y_list,
  R = 2,                # 2D latent space
  family = "binary",    # probit for 0/1 networks
  burn = 100,           # burn-in (use 1000+ for real work)
  nscan = 500,          # post-burn-in samples (use 5000+)
  odens = 25,           # thinning
  verbose = FALSE,
  plot = FALSE
)

summary(fit)
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> NULL
#> 
#> Time periods: 3 
#> Family: binary 
#> Mode: unipartite 
#> 
#> Regression coefficients:
#> ------------------------
#>           Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept   -1.314    0.183  -7.199       0   -1.672   -0.956 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.117    0.050
#> cab    0.001    0.041
#> vb     0.099    0.033
#> rho    0.313    0.183
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
```

That’s it. You now have posterior estimates for regression coefficients,
variance components, and latent positions pooled across all time
periods.

## What Can You Do with a Fitted Model?

`lame` objects work with all the standard R methods you’d expect:

``` r
# Regression coefficients (posterior means)
coef(fit)
#> intercept 
#>  -1.31421

# 95% credible intervals
confint(fit)
#>                2.5%     97.5%
#> intercept -1.632783 -1.026529

# Predicted probabilities for every dyad at every time point
Y_hat <- predict(fit, type = "response")
cat("Predicted probability range:",
    round(range(unlist(Y_hat), na.rm = TRUE), 3), "\n")
#> Predicted probability range: 0.014 0.625

# Residuals
resid_list <- residuals(fit)
cat("Residual SD:", round(sd(unlist(resid_list), na.rm = TRUE), 3), "\n")
#> Residual SD: 0.35
```

## Checking Your Model

Two essential diagnostics:

**1. Did the MCMC converge?** The trace plots should look like fuzzy
caterpillars, not random walks with trends.

``` r
trace_plot(fit, params = "beta")
```

![](lame_files/figure-html/trace-plot-1.png)

**2. Does the model fit the data?** The GOF plots compare observed
network statistics to what the model predicts. Observed values should
fall within the simulated distributions.

``` r
gof_plot(fit)
```

![](lame_files/figure-html/gof-plot-1.png)

## Visualizing Network Structure

The latent space shows which actors play similar roles in the network:

``` r
uv_plot(fit)
```

![](lame_files/figure-html/uv-plot-1.png)

And the additive effects show who’s unusually active or popular:

``` r
ab_plot(fit, effect = "sender")
```

![](lame_files/figure-html/ab-plot-1.png)

## Cross-Sectional Models

If you have a single network (not a time series), use
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md) instead of
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md):

``` r
fit_cs <- ame(
  Y = Y_list[[1]],     # just one time period
  R = 2,
  family = "binary",
  burn = 100,
  nscan = 500,
  odens = 25,
  verbose = FALSE
)

coef(fit_cs)
#> intercept 
#>  -1.21051
```

The output and methods are the same. The difference is that
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) pools
information across time periods, giving you more precise estimates when
the structure is stable.

## Dynamic Effects

When you believe the network structure is changing over time (alliances
shifting, friendships evolving), you can let the latent positions and
additive effects drift via AR(1) processes:

``` r
fit_dyn <- lame(
  Y = Y_list,
  R = 2,
  dynamic_ab = TRUE,    # time-varying sociality/popularity
  dynamic_uv = TRUE,    # time-varying latent positions
  family = "binary",
  burn = 100,
  nscan = 500,
  odens = 25,
  verbose = FALSE,
  plot = FALSE
)

summary(fit_dyn)
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> NULL
#> 
#> Time periods: 3 
#> Family: binary 
#> Mode: unipartite 
#> Dynamic latent positions: enabled (rho_uv = 0.322 )
#> Dynamic additive effects: enabled (rho_ab = 0.294 )
#> 
#> Regression coefficients:
#> ------------------------
#>           Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept   -1.199    0.158  -7.572       0   -1.509   -0.889 ***
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.114    0.058
#> cab    0.042    0.050
#> vb     0.123    0.051
#> rho    0.378    0.120
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
```

The model estimates how persistent the latent positions are
($\rho_{uv}$): values near 1 mean slow evolution, values near 0 mean
rapid change. See the [dynamic effects
vignette](https://netify-dev.github.io/lame/articles/dynamic_effects.md)
for a deeper dive.

## Bipartite Networks

For two-mode networks (students and courses, countries and treaties),
pass a rectangular matrix and set `mode = "bipartite"`:

``` r
nA <- 10; nB <- 8
Y_bip <- lapply(1:3, function(t) {
  Y <- matrix(rbinom(nA * nB, 1, 0.2), nA, nB)
  rownames(Y) <- paste0("R", 1:nA)
  colnames(Y) <- paste0("C", 1:nB)
  Y
})

fit_bip <- lame(
  Y = Y_bip,
  mode = "bipartite",
  R = 2,
  family = "binary",
  burn = 100, nscan = 500, odens = 25,
  verbose = FALSE, plot = FALSE
)

summary(fit_bip)
#> 
#> === Longitudinal AME Model Summary ===
#> 
#> Call:
#> NULL
#> 
#> Time periods: 3 
#> Family: binary 
#> Mode: bipartite 
#> 
#> Regression coefficients:
#> ------------------------
#>           Estimate StdError z_value p_value CI_lower CI_upper  
#> intercept   -0.895    0.438  -2.045   0.041   -1.753   -0.037 *
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.589    0.224
#> cab    0.000    0.000
#> vb     0.770    0.322
#> rho    0.000    0.000
#> ve     1.000    0.000
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
```

See the [bipartite
vignette](https://netify-dev.github.io/lame/articles/bipartite.md) for a
full walkthrough.

## Supported Data Types

| Family      | Data                    | Example                                 |
|-------------|-------------------------|-----------------------------------------|
| `"normal"`  | Continuous              | Trade volumes, survey ratings           |
| `"binary"`  | 0/1                     | Friendships, alliances, sanctions       |
| `"ordinal"` | Ordered categories      | Conflict intensity (none/threat/action) |
| `"poisson"` | Counts                  | Number of co-sponsored bills            |
| `"tobit"`   | Non-negative continuous | Trade with zeros                        |
| `"cbin"`    | Censored binary         | Friendships with nomination limits      |
| `"frn"`     | Fixed rank nomination   | “Name your top 5 friends”               |
| `"rrl"`     | Row-ranked likelihood   | Ranked preferences                      |

## Where to Go Next

| I want to…                                  | Read this                                                                           |
|---------------------------------------------|-------------------------------------------------------------------------------------|
| Understand the full workflow with real data | [lame overview](https://netify-dev.github.io/lame/articles/lame-overview.md)        |
| Learn about cross-sectional models in depth | [Your first AME model](https://netify-dev.github.io/lame/articles/cross_sec_ame.md) |
| Model two-mode networks                     | [Bipartite networks](https://netify-dev.github.io/lame/articles/bipartite.md)       |
| Let the network structure evolve over time  | [Dynamic effects](https://netify-dev.github.io/lame/articles/dynamic_effects.md)    |

## References

Hoff, PD (2021). Additive and Multiplicative Effects Network Models.
*Statistical Science*, 36, 34–50.

Minhas, S., Dorff, C., Gallop, M. B., Foster, M., Liu, H., Tellez, J., &
Ward, M. D. (2022). Taking dyads seriously. *Political Science Research
and Methods*, 10(4), 703–721.
