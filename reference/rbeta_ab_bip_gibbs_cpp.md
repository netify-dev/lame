# Full bipartite Gibbs update for beta, a, b

Single C++ function replacing the bipartite beta/a/b update block in
lame.R

## Usage

``` r
rbeta_ab_bip_gibbs_cpp(
  Z,
  Xlist,
  UV_eff,
  a_current,
  b_current,
  s2,
  g_prior,
  va,
  vb,
  rvar,
  cvar
)
```

## Arguments

- Z:

  3D array (nA x nB x T)

- Xlist:

  List of T arrays (nA x nB x p)

- UV_eff:

  nA x nB matrix (U*G*V' or U\*V')

- a_current:

  Current row effects (length nA)

- b_current:

  Current column effects (length nB)

- s2:

  Dyadic variance

- g_prior:

  G-prior parameter

- va:

  Row effect variance (diagonal element of Sab)

- vb:

  Column effect variance (diagonal element of Sab)

- rvar:

  Whether to update row effects

- cvar:

  Whether to update column effects

## Value

List with beta, a, b
