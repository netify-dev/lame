# Rotation-drift diagnostic for the canonical (U_t, G_t, V_t) trio

Reports a scalar ratio comparing the variance of the *raw* G_t entries
(per element across t) to the variance of the *canonical* G_t entries
(after per-period SVD). When the raw / canonical variance ratio is large
(\>= 5), the apparent G_t time-variation is dominated by rotation drift
in U_t, V_t rather than real temporal change. The canonicalisation
removes that rotation and the canonical entries should be the
user-facing summary.

## Usage

``` r
.rotation_drift_diagnostic(G_cube_raw, U_cube, V_cube)
```

## Arguments

- G_cube_raw:

  RA x RB x T raw (FFBS-sampled) G cube

- U_cube:

  nA x RA x T raw U cube

- V_cube:

  nB x RB x T raw V cube

## Value

list with `ratio` (numeric), `var_raw` (RA x RB), `var_canonical` (RA x
RB), `flag` (logical: ratio \>= 5)
