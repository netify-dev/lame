# Recovery grid

`tests/testthat/helper_recovery_grid.R` simulates networks from the model's
own generative process with known parameters, fits them with the estimator a
cell names (`ame`, `lame`, `ame_als`, `lame(method = "als")`, the dynamic ALS
estimator), scores recovery against the realized truth, and checks the
post-processing methods against exact identities. One row per fit.

Conventions the scores rest on: symmetric data carry one latent per unordered
pair; variance components are scored against the realized sample variance of
the drawn effects; first moments use `fit$YPM`; with a multiplicative term the
effects are scored against `a + rowMeans(U G V')` and the intercept against
the identified level.

`test_recovery_grid.R` runs a fast subset on every check. The full grid runs
from the source tree:

```r
source("tests/testthat/helper_recovery_grid.R")
library(lame)
for (g in c("G1_static_uni", "G2_bipartite", "G3_dynamic_uni", "G4_als", "G5_kinds", "G6_dynrho")) {
	df <- grid_run(grid_cells(g), seeds = 1:2, out_csv = paste0(g, ".csv"))
	print(grid_systematic(df))
}
```

`grid_run()` is resumable (finished cells are skipped), `grid_flags()` marks
rows whose scores fall outside the triage thresholds, `grid_systematic()`
lists the flags that fire in every seed of a cell, and `grid_repro(cell, seed)`
refits one cell for diagnosis. Flags are thresholds for triage, not verdicts:
read the numbers.
