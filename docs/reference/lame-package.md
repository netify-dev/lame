# Longitudinal Additive and Multiplicative Effects Models for Networks

An extension of the 'amen' package specifically designed for
longitudinal network analysis. This package provides tools for fitting
Additive and Multiplicative Effects (AME) models to
longitudinal/replicated relational data with several extensions: (1)
handling of changing actor compositions across time periods, allowing
for networks with different sets of actors at each time point; (2)
significant performance improvements through C++ implementations via
Rcpp and RcppArmadillo; (3) specialized functions for temporal network
dynamics. The package supports eight data types: normal (nrm), binary
(bin), ordinal (ord), Poisson count (poisson), tobit/censored continuous
(tobit), censored binary (cbin), fixed-rank nomination (frn), and
row-ranked (rrl). Based on the AME framework originally developed by
Hoff (2009) and Hoff, Fosdick, Volfovsky and Stovel (2013).

## Details

**Estimators.** The package offers two estimation routes:

- [`ame`](https://netify-dev.github.io/lame/reference/ame.md) /
  [`lame`](https://netify-dev.github.io/lame/reference/lame.md) — the
  Bayesian MCMC estimators, for calibrated posterior inference
  (cross-sectional and longitudinal respectively).

- [`ame_als`](https://netify-dev.github.io/lame/reference/ame_als.md) /
  [`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
  — a fast, MCMC-free point estimator by iterative block coordinate
  descent, with bootstrap uncertainty via
  [`ame_als_bootstrap`](https://netify-dev.github.io/lame/reference/ame_als_bootstrap.md).
  Use it for rapid model exploration, rank selection and starting
  values; use the MCMC estimators for final inference.

|          |         |
|----------|---------|
| Package: | lame    |
| Type:    | Package |
| Version: | 1.0.0   |
| Date:    | 2026    |
| License: | MIT     |

## Author

Shahryar Minhas, Tosin Salau, Cassy Dorff

Maintainer: Shahryar Minhas <minhassh@msu.edu>

## Examples

``` r
# \donttest{
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, burn = 10, nscan = 100, odens = 1,
           family = "normal", verbose = FALSE)
summary(fit)
#> 
#> Regression coefficients:
#>                 pmean   psd  z-stat p-val
#> intercept_dyad -2.709 0.115 -23.527 0.000
#> rgpa_dyad       0.219 0.101   2.165 0.030
#> rsmoke_dyad     0.263 0.120   2.186 0.029
#> cgpa_dyad       0.164 0.039   4.165 0.000
#> csmoke_dyad     0.112 0.056   1.995 0.046
#> igrade_dyad     1.129 0.037  30.452 0.000
#> igpa_dyad       0.055 0.016   3.413 0.001
#> ismoke_dyad     0.029 0.025   1.166 0.243
#> 
#> Variance parameters:
#>     pmean   psd
#> va  0.853 0.120
#> cab 0.094 0.034
#> vb  0.120 0.019
#> rho 0.920 0.002
#> ve  1.053 0.020
# }


```
