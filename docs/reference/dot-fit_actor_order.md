# The fit's internal (canonical) actor order

[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) sorts
actors alphabetically when it ingests a list of `Y` matrices (via
`list_to_array`), so the stored additive / multiplicative effects and
design slices are in SORTED actor order, which need not match the order
a user later supplies in `newdata`. This returns the row and column
actor names the fit was estimated in, tried from the most authoritative
source down, so `newdata` can be realigned by name before prediction.

## Usage

``` r
.fit_actor_order(object)
```

## Arguments

- object:

  a fitted `lame` (or `ame`) object.

## Value

a list with `rows` and `cols` character vectors (or `NULL` when the
order cannot be recovered).
