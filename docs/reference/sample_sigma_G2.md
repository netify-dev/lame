# Inverse-gamma posterior draw for sigma_G^2 given the FFBS path

Conjugate IG update on `sigma_G^2` given the AR(1) innovations
`eta_t = g_t - rho * g_{t-1}`, t = 2..N, plus the stationary density at
t = 1. Default prior is IG(2, 1) on `sigma_G^2`.

## Usage

``` r
sample_sigma_G2(
  vecG_path,
  rho_G,
  prior_shape = 2,
  prior_rate = 1,
  s2_obs = 1,
  v_cap_mult = 4
)
```

## Arguments

- s2_obs:

  observation-noise variance (1 for probit/binary).

- v_cap_mult:

  cap on the stationary G-state variance in units of `s2_obs` (default
  4).

## Details

Scale identification. The model identifies only the product \\U_t G_t
V_t'\\, not \\G_t\\ alone, so the overall scale of `vec(G_t)` is free:
left unchecked the chain finds a degenerate mode where `sigma_G^2` (and
hence `G_t`) inflates by orders of magnitude while `U,V` collapse to
compensate, leaving the linear predictor unchanged but the reported
`G_cube` meaningless. Clamping `rho_G` alone does not bound this because
the stationary state variance is `sigma_G^2 / (1 - rho_G^2)`. We
therefore cap the implied stationary variance at `v_cap_mult * s2_obs`
(a few observation-noise units), which forces the multiplicative scale
onto `U,V` – where the regularising N(0, s2) prior pins it – and keeps
`G_t` on a scale comparable to a static-G fit. The cap is loose enough
never to bind on a genuinely small-variation `G_t`.
