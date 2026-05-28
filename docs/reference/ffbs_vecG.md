# Forward-filter / backward-sample for vec(G_t) under AR(1) state prior

One Carter-Kohn sweep over the per-period observations of the form
`vec(E_t) = H_t g_t + eps_t`, with state transition
`g_t = rho * g_{t-1} + eta_t`.

## Usage

``` r
ffbs_vecG(E_cube, U_cube, V_cube, s2, rho_G, sigma_G2)
```

## Arguments

- E_cube:

  nA x nB x T residual cube (Z minus base - a - b - UV')

- U_cube:

  nA x RA x T latent row factor cube

- V_cube:

  nB x RB x T latent column factor cube

- s2:

  scalar observation variance

- rho_G:

  AR(1) coefficient in `(-1, 1)` (use `rho = 1` for the RW1 limit; the
  forward variance is clamped via a tiny floor).

- sigma_G2:

  state innovation variance

## Value

list with `G_cube` (RA x RB x T) and `vecG_path` (p x T draws of the
vectorised state).
