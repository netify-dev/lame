# Extract fitted values from AME model

Returns the posterior mean of the network on the response scale (YPM).
For binary models, these are predicted probabilities between 0 and 1.
For normal models, these are predicted continuous values. For Poisson
models, these are predicted counts.

## Usage

``` r
# S3 method for class 'ame'
fitted(object, ...)
```

## Arguments

- object:

  Fitted AME model object (class "ame").

- ...:

  Additional arguments (not used).

## Value

An n x n matrix (unipartite) or nA x nB matrix (bipartite) of fitted
values on the response scale. Diagonal entries are `NA` for unipartite
networks.

## See also

[`predict.ame`](https://netify-dev.github.io/lame/reference/predict.ame.md)
for predictions with type control,
[`residuals.ame`](https://netify-dev.github.io/lame/reference/residuals.ame.md)
for residuals

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
