# Longitudinal Additive and Multiplicative Effects Models for Networks

An extension of the 'amen' package specifically designed for
longitudinal network analysis. This package provides tools for fitting
Additive and Multiplicative Effects (AME) models to
longitudinal/replicated relational data with several key enhancements:
(1) robust handling of changing actor compositions across time periods,
allowing for networks with different sets of actors at each time point;
(2) significant performance improvements through C++ implementations via
Rcpp and RcppArmadillo; (3) specialized functions for temporal network
dynamics. The package supports eight data types: normal (nrm), binary
(bin), ordinal (ord), Poisson count (poisson), tobit/censored continuous
(tobit), censored binary (cbin), fixed-rank nomination (frn), and
row-ranked (rrl). Based on the AME framework originally developed by
Hoff (2009) and Hoff, Fosdick, Volfovsky and Stovel (2013).

## Details

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
#> === AME Model Summary ===
#> 
#> Call:
#> [1] "Y ~ dyad(intercept, rgpa, rsmoke, cgpa, csmoke, igrade, igpa, ismoke) + a[i] + b[j] + rho*e[ji], family = 'normal'"
#> 
#> Regression coefficients:
#> ------------------------
#>                Estimate StdError z_value p_value CI_lower CI_upper    
#> intercept_dyad   -2.709    0.115 -23.527       0   -2.901   -2.515 ***
#> rgpa_dyad         0.219    0.101   2.165    0.03    0.025    0.444   *
#> rsmoke_dyad       0.263     0.12   2.186   0.029    0.051    0.515   *
#> cgpa_dyad         0.164    0.039   4.165       0    0.091    0.241 ***
#> csmoke_dyad       0.112    0.056   1.995   0.046        0    0.205   *
#> igrade_dyad       1.129    0.037  30.452       0    1.059    1.197 ***
#> igpa_dyad         0.055    0.016   3.413   0.001     0.02    0.081 ***
#> ismoke_dyad       0.029    0.025   1.166   0.243   -0.018    0.075    
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Note: p-values are approximate (posterior mean / SD); use credible intervals for inference.
#> 
#> Variance components:
#> -------------------
#>     Estimate StdError
#> va     0.853    0.120
#> cab    0.094    0.034
#> vb     0.120    0.019
#> rho    0.920    0.002
#> ve     1.053    0.020
#>   (va = sender, cab = sender-receiver covariance, vb = receiver,
#>    rho = dyadic correlation, ve = residual variance)
# }


```
