# Realign a newdata covariate array to the fit's internal actor order

Reorders the rows and columns of each `n x m x p` slice by their
dimnames so they line up with the fit's stored (sorted) additive /
multiplicative effects. Only acts when `newdata` carries dimnames whose
name set exactly matches the fit's; otherwise it falls back to the
supplied positional order unchanged (the historical behaviour, correct
when the user already passes actors in the fit's order or supplies no
names).

## Usage

``` r
.reorder_newdata_actors(X, row_order, col_order)
```

## Arguments

- X:

  an `n x m x p` array (already promoted from a matrix).

- row_order, col_order:

  canonical actor names from
  [`.fit_actor_order`](https://netify-dev.github.io/lame/reference/dot-fit_actor_order.md).

## Value

`X` with rows / columns permuted into canonical order, or unchanged.
