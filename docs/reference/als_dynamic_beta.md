# Penalised ALS time-varying coefficient estimate

Computes a fast point estimate of a time-varying coefficient vector
\\\beta_t\\ for a longitudinal network model by solving the
first-difference-penalised least-squares problem. No MCMC, no
multiplicative effects, no random effects — purely a regression point
estimate with a smoothing penalty on the coefficient path. Useful for
rapid exploration and for generating starting values for a full
`lame(dynamic_beta = ...)` fit.

## Usage

``` r
als_dynamic_beta(Y, Xdyad, lambda = 0, intercept = TRUE)
```

## Arguments

- Y:

  A list of \\T\\ response matrices (or a 3-D array with third dimension
  time).

- Xdyad:

  A list of \\T\\ dyadic covariate arrays (each \\n \times n \times p\\,
  or a single \\n \times n\\ matrix if \\p=1\\).

- lambda:

  Non-negative smoothing parameter. Default `0` (no penalty, per-period
  OLS).

- intercept:

  Logical; include an intercept (default `TRUE`).

## Value

A list with `beta` (a \\p \times T\\ matrix of estimates), `lambda`,
`residual_ss` (sum of squared residuals across periods), `intercept`,
and `call`. Class `"als_dynamic_beta"`.

## Details

\$\$\min\_{\beta\_{1:T}} \sum_t \\y_t - X_t \beta_t\\^2 + \lambda
\sum\_{t=2}^{T} \\\beta_t - \beta\_{t-1}\\^2\$\$

At \\\lambda = 0\\ this is \\T\\ independent per-period OLS fits. As
\\\lambda \to \infty\\ the solution converges to a single pooled
\\\beta\\ (constant over time).

## Relation to other entry points

`als_dynamic_beta` is a *regression-only smoother*; the returned object
carries time-varying \\\beta_t\\ alone, no additive (\\a, b\\) or
multiplicative (\\U, V\\) AME components. It is therefore not a special
case of
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
(which always estimates the full AME decomposition) and not a special
case of
[`lame`](https://netify-dev.github.io/lame/reference/lame.md)`(dynamic_beta = TRUE)`
(which is Bayesian with an AR(1) / RW1 / RW2 / Matern 3/2 state-space
prior). Use it when you want a fast, point-only, deterministic smoother
for the coefficient path; use
[`lame`](https://netify-dev.github.io/lame/reference/lame.md)`(dynamic_beta = ...)`
when you need full posterior uncertainty and additive / multiplicative
effects.

## See also

[`lame`](https://netify-dev.github.io/lame/reference/lame.md) for the
full Bayesian dynamic-\\\beta\\ fit (AR(1) / RW1 / RW2 / Matern 3/2);
[`lame_als`](https://netify-dev.github.io/lame/reference/lame_als.md)
for the longitudinal AME point estimator (static \\\beta\\, full \\a, b,
U, V\\).

## Examples

``` r
set.seed(1)
n <- 10; T <- 4
X <- replicate(T, array(rnorm(n*n*2), c(n, n, 2)), simplify = FALSE)
beta_true <- rbind(seq(-1, 1, length.out = T), seq(0.5, -0.5, length.out = T))
Y <- vector("list", T)
for (t in seq_len(T)) {
  Yt <- X[[t]][, , 1] * beta_true[1, t] + X[[t]][, , 2] * beta_true[2, t] +
        matrix(rnorm(n*n, 0, 0.2), n, n)
  diag(Yt) <- NA
  Y[[t]] <- Yt
}
fit_als <- als_dynamic_beta(Y, X, lambda = 0)     # per-period LS
fit_smooth <- als_dynamic_beta(Y, X, lambda = 10) # smoother path
```
