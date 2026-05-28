# ERGM-style covariate helpers for ame() / lame()

`nodematch(x)` returns an `n x n` matrix with `1` where `x[i] == x[j]`
and `0` otherwise – the AME analogue of ERGM's `nodematch` term,
suitable to drop into `Xdyad` as a single homophily covariate.
`absdiff(x)` returns the absolute difference `|x[i] - x[j]|`, the AME
analogue of ERGM's `absdiff`. `nodefactor(x)` returns the dyadic matrix
of `x[i] + x[j]` for a numeric `x` (or a list of per-level binary dyadic
indicators for a factor / character `x`).

## Usage

``` r
nodematch(x, na_diag = TRUE)

absdiff(x, na_diag = TRUE)

nodefactor(x, na_diag = TRUE)
```

## Arguments

- x:

  a numeric, factor, or character vector of length `n`.

- na_diag:

  logical: set the diagonal to `NA`? Defaults to `TRUE`.

## Value

An `n x n` numeric matrix (or, for a factor `x` passed to `nodefactor`,
a named list of such matrices).

## Details

All helpers return a matrix *or* list of matrices that is ready to wrap
into an `n x n x p` `Xdyad` array via `simplify2array` or
[`array()`](https://rdrr.io/r/base/array.html).

## Examples

``` r
# \donttest{
n <- 12
grp <- sample(letters[1:3], n, replace = TRUE)
age <- rnorm(n, 40, 10)
# build a single same-group homophily covariate
Xdyad <- array(nodematch(grp), dim = c(n, n, 1),
               dimnames = list(NULL, NULL, "same_group"))
# combine homophily + age difference
Xdyad <- array(c(nodematch(grp), absdiff(age)),
               dim = c(n, n, 2),
               dimnames = list(NULL, NULL, c("same_group", "age_diff")))
# }
```
