# Longitudinal Additive and Multiplicative Effects Models for Networks

An extension of the 'amen' package specifically designed for
longitudinal network analysis. This package provides tools for fitting
Additive and Multiplicative Effects (AME) models to
longitudinal/replicated relational data with several key enhancements:
(1) robust handling of changing actor compositions across time periods,
allowing for networks with different sets of actors at each time point;
(2) significant performance improvements through C++ implementations via
Rcpp and RcppArmadillo; (3) specialized functions for temporal network
dynamics. The package supports various data types including
binary/network data (bin), normal relational data (nrm), ordinal
relational data (ord), censored binary data (cbin), fixed-rank
nomination schemes (frn), and row-ranked data (rrl). Based on the AME
framework originally developed by Hoff (2009) and Hoff, Fosdick,
Volfovsky and Stovel (2013).

## Details

|          |            |
|----------|------------|
| Package: | lame       |
| Type:    | Package    |
| Version: | 0.0.0.9000 |
| Date:    | 2025       |
| License: | MIT        |

## Author

Shahryar Minhas, Tosin Salau, Cassy Dorff

Maintainer: Shahryar Minhas <sminhas@example.com>

## Examples

``` r
# \donttest{
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, burn = 10, nscan = 100, odens = 1,
           family = "normal", print = FALSE)
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
#> intercept_dyad   -2.709    0.115 -23.527       0   -2.934   -2.483 ***
#> rgpa_dyad         0.219    0.101   2.165    0.03    0.021    0.418   *
#> rsmoke_dyad       0.263     0.12   2.186   0.029    0.027    0.499   *
#> cgpa_dyad         0.164    0.039   4.165       0    0.087    0.242 ***
#> csmoke_dyad       0.112    0.056   1.995   0.046    0.002    0.222   *
#> igrade_dyad       1.129    0.037  30.452       0    1.056    1.202 ***
#> igpa_dyad         0.055    0.016   3.413   0.001    0.023    0.086 ***
#> ismoke_dyad       0.029    0.025   1.166   0.243    -0.02    0.079    
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
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
