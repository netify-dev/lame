# TIES sanctions data for vignettes

Longitudinal directed binary network of international economic sanctions
among 35 countries across four years (1993, 1994, 1995, 2000), derived
from the Threat and Imposition of Sanctions (TIES) dataset (Morgan et
al. 2014). Entry \\y\_{ij,t} = 1\\ means country \\i\\ imposed sanctions
on country \\j\\ in year \\t\\. Network density is approximately 2\\

## Usage

``` r
data(vignette_data)
```

## Format

Four objects:

- Y:

  List of 4 binary adjacency matrices (35 x 35), one per year.

- Xdyad:

  List of 4 arrays (35 x 35 x 2) of dyadic covariates: `distance`
  (geographic) and `shared_igos` (shared IGO memberships).

- Xrow:

  List of 4 matrices (35 x 2) of sender covariates: `log_gdp` and
  `log_pop`.

- Xcol:

  List of 4 matrices (35 x 2) of receiver covariates: `log_gdp` and
  `log_pop`.

## Details

When loaded via `data(vignette_data)`, the following objects are placed
in the calling environment: `Y`, `Xdyad`, `Xrow`, `Xcol`.

## References

Morgan, T. Clifton, Bapat, N., & Kobayashi, Y. (2014). Threat and
Imposition of Economic Sanctions 1945–2005. *Conflict Management and
Peace Science*, 31(5), 541–558.

## Examples

``` r
data(vignette_data)
cat("Countries:", nrow(Y[[1]]), "\n")
#> Countries: 35 
cat("Time periods:", length(Y), "\n")
#> Time periods: 4 
```
