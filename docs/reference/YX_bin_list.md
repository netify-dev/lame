# Synthetic longitudinal binary relational data, list-form (latent-scale)

The list-form sibling of
[`YX_bin_long`](https://netify-dev.github.io/lame/reference/YX_bin_long.md):
same synthetic 4-period 50-actor panel reshaped from arrays to lists.
The shipped `Y` entries hold the *latent probit-scale predictor* (range
roughly -20 to 17), not 0/1 ties. To recover the binary tie indicator
the dataset name implies, threshold at zero, e.g.
`Y_bin <- lapply(YX_bin_list$Y, function(Yt) 1 * (Yt > 0))`. If passed
unthresholded to `lame(..., family = "binary")` the fit will warn and
silently apply the same `Y > 0` threshold.

## Format

A list with two elements:

- Y:

  List of 4 numeric matrices, each `[50, 50]`, with `NA` on the
  diagonal. Latent probit-scale predictor; threshold at 0 for the binary
  indicator.

- X:

  List of 4 numeric arrays, each `[50, 50, 3]` of dyadic covariates.
