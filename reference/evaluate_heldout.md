# Held-out predictive evaluation for an ame / lame fit

Computes family-appropriate held-out predictive scores given a fit and a
logical mask of cells to score. The function does not split the data for
you: it expects you to have either (a) refit on a training subset and
now want to score the held-out cells of the same matrix, or (b) have
predicted probabilities you want scored. Both AUROC and PR-AUC are
reported when applicable, alongside Brier and mean log density.

## Usage

``` r
evaluate_heldout(y_obs, y_pred, mask, family = "binary")
```

## Arguments

- y_obs:

  Observed outcomes. For longitudinal fits, a list of per-period
  matrices; for cross-sectional, a single matrix. Cells not in `mask`
  are ignored.

- y_pred:

  Predicted probabilities / means on the response scale. Same shape as
  `y_obs`.

- mask:

  Logical mask of the same shape as `y_obs` marking cells to score
  (`TRUE` = include in evaluation).

- family:

  Family string; used to pick the scoring rule. Defaults to `"binary"`.

## Value

A one-row data frame with columns appropriate to the family: `n_eval`,
plus `auroc` + `auprc` + `brier` + `logloss` (binary / cbin), or
`rmse` + `mae` + `mean_logdens` (normal / tobit), or `mean_logdens` +
`rmse` (poisson).

## Details

**Workflow.** The typical pattern is: mask a random sample of dyads to
`NA` in `Y`, refit
([`lame()`](https://netify-dev.github.io/lame/reference/lame.md) handles
`NA` internally via data augmentation), call
`predict(fit, type = "response")`, then pass that prediction alongside
the original `Y` and the held-out mask to this function. See the
examples.

**Dependencies.** AUROC / PR-AUC use precrec when available; if not
installed, only Brier and mean log-density are computed and a one-line
note is emitted.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 25; Y <- matrix(rbinom(n*n, 1, 0.3), n, n); diag(Y) <- NA
rownames(Y) <- colnames(Y) <- paste0("a", sprintf("%02d", 1:n))
# mask 20% of dyads
mask <- matrix(FALSE, n, n)
obs_idx <- which(!is.na(Y))
set.seed(1)
mask[sample(obs_idx, floor(0.2 * length(obs_idx)))] <- TRUE
Y_train <- Y; Y_train[mask] <- NA
fit <- ame(Y_train, R = 0, family = "binary",
           burn = 50, nscan = 200, odens = 5, verbose = FALSE, plot = FALSE)
y_pred <- predict(fit, type = "response")
evaluate_heldout(Y, y_pred, mask, family = "binary")
#>   n_eval     auroc     auprc    brier   logloss
#> 1    120 0.5659341 0.2869632 0.190125 0.5682987
# }
```
