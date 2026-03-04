# Print method for AME model objects

Displays a formatted summary of a fitted AME (Additive and
Multiplicative Effects) model. This method provides a concise overview
of model structure, parameter estimates, and goodness-of-fit statistics
without generating new data.

## Usage

``` r
# S3 method for class 'ame'
print(x, ...)
```

## Arguments

- x:

  an object of class "ame" from fitting an AME model

- ...:

  additional arguments (not currently used)

## Value

Invisibly returns the input object (for method chaining)

## Details

The print method displays:

- Model type (unipartite/bipartite, symmetric/asymmetric)

- Network dimensions and observation count

- Family and link function used

- Number of MCMC iterations

- Parameter counts for regression coefficients and latent factors

- Basic convergence diagnostics if available

Unlike `simulate`, this method only formats existing results for display
and does not perform any new computations or data generation.

## See also

[`summary.ame`](https://netify-dev.github.io/lame/reference/summary.ame.md)
for detailed summaries,
[`simulate.ame`](https://netify-dev.github.io/lame/reference/simulate.ame.md)
for generating new networks,
[`predict.ame`](https://netify-dev.github.io/lame/reference/predict.ame.md)
for predictions

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Fit model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, R = 2,
           nscan = 100, burn = 10, odens = 1, print = FALSE)

# Display summary
print(fit)
#> 
#> Additive and Multiplicative Effects (AME) Model
#> ================================================
#> 
#> Network dimensions:  100 x 100 
#> MCMC iterations:  100 
#> Family:  normal 
#> Mode:  unipartite 
#> Symmetric:  FALSE 
#> 
#> Number of parameters:
#>   Regression coefficients:  8 
#>   Row/sender effects: enabled
#>   Column/receiver effects: enabled
#>   Dyadic correlation: enabled
#>   Multiplicative effects dimension:  2 
#> 
#> Use summary(object) for detailed results
# }
```
