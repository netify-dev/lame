# Reconstruct EZ and UVPM matrices from AME model output

Helper functions to reconstruct EZ (linear predictor) and UVPM
(posterior mean of UV) matrices that are no longer stored in the model
output to save memory.

Note: EZ returns the LINEAR PREDICTOR (\\\eta\\), not the response:

- For Gaussian: \\EZ = \eta = \mu\\ (identity link)

- For Poisson: \\EZ = \eta = \log(\lambda)\\ (can be negative)

- For Binary: \\EZ = \eta\\ = probit inverse of p (can be any real
  value) Use YPM for predictions on the response scale.

## Usage

``` r
reconstruct_EZ(fit, X = NULL)

reconstruct_UVPM(fit)
```

## Arguments

- fit:

  Fitted AME model object

- X:

  Covariate array (optional, will use fit\$X if available)

## Value

Reconstructed matrix

## Details

These matrices are not stored by default to save memory, but can be
reconstructed when needed for diagnostics or analysis.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau
