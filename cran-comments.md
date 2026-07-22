# cran-comments

## Resubmission of lame, now 1.3.3 (new package)

Thank you for the review. Every point raised has been addressed:

* **Single quotes around function names.** Removed; the Description now
  writes `ame()` and `lame()` unquoted. Package names (`Rcpp`,
  `RcppArmadillo`) remain single-quoted.
* **Redundant "Tools for".** Removed from the start of the Description.
* **References in the form authors (year) <doi:...>.** Added: the
  Description now cites Sewell and Chen (2015) <doi:10.1080/01621459.2014.988214>
  and Durante and Dunson (2014) <doi:10.1093/biomet/asu040>.
* **`T`/`F` instead of `TRUE`/`FALSE`.** The three flagged functions
  (`init_bipartite_startvals()`, `init_dynamic_ab_cpp()`,
  `init_dynamic_positions()`) used `T` as the argument holding the number of
  time points; it is renamed to `Tn` throughout, in the R sources, the C++
  sources, and the regenerated `RcppExports`. A further sweep of the sources
  found and fixed four uses of bare `T` as a value (`diag = T`,
  `byrow = T`) and one `T` variable in an example.
* **`options(warn = -1)`.** Removed entirely from `R/lame.R`. The call was a
  legacy guard against warnings from an older sampler; we verified across
  the symmetric binary, normal, ordinal, high-rank, goodness-of-fit and
  dynamic-UV paths that no warnings are emitted, so no replacement
  `suppressWarnings()` was needed and genuine warnings now reach the user.
* **Writing to `.GlobalEnv`.** The fitting functions set a seed internally
  and previously left the resulting `.Random.seed` in the user's workspace
  when none had existed. All nine sites now restore the caller's
  `.Random.seed` on exit and remove it if the session had none, so a fit
  leaves the global environment byte-identical to how it found it. This is
  verified in both directions in the test suite.
* **Authors, contributors and copyright holders.** Several samplers and
  helpers in this package are derived from Peter Hoff's `amen` package.
  Peter Hoff has been added to `Authors@R` with `ctb` and `cph` roles, and
  a new `inst/COPYRIGHTS` file records the derivation and lists the affected
  routines. Because `amen` is distributed under GPL-3, **the license of this
  package has been changed from MIT to GPL-3** so that the upstream license
  is preserved rather than weakened; the former MIT `LICENSE` file has been
  removed accordingly.

The version was bumped to 1.3.3; it also folds in bug fixes made since the
previous submission (see NEWS.md). The vignettes were compacted so the suite
rebuilds quickly on CRAN's builders.

## Test environments

* local: Ubuntu 24.04 (WSL2), R 4.3.3 — `R CMD check --as-cran`
* (please also see win-builder / R-hub results if requested)

## R CMD check results

0 errors | 0 warnings | 4 notes

* NOTE: New submission.
* NOTE: installed size (sub-directory `libs`). The local check reports a
  28.1Mb installed package, including 21.2Mb in `libs`. The package ships
  C++ MCMC samplers built with Rcpp/RcppArmadillo; all compiled code is
  required at runtime. CRAN's Windows pretest reports the installed package
  size as OK.
* NOTE (local only): future file timestamps -- the WSL2 build environment
  cannot verify the current time; not expected on CRAN's builders.
* NOTE (local only): non-portable compiler flag `-mno-omit-leaf-frame-pointer`
  comes from the local Ubuntu r-base build's default Makeconf flags, not from
  the package's Makevars; not expected on CRAN's builders.

## Comments

* The `netify` (>= 1.5.3) dependency was published on CRAN on 2026-06-30.
* Long-running MCMC estimation is exercised in `\donttest{}` examples with
  deliberately short chains, and the full test suite gates its long MCMC
  tests behind `skip_on_cran()`; CRAN-visible test time is about 30s.
* Vignettes use compact live examples. Rank-selection and PSIS-LOO examples
  retain their runnable code but are not evaluated during the build because
  those workflows require longer converged chains. Rebuilding all vignette
  outputs takes about 75 seconds on the local check environment.
* `ame()`/`lame()` and the ALS fitters accept an optional `seed` argument;
  the user's `.Random.seed` is saved and restored on exit (and removed if
  the session had none) so fitting a model never perturbs the caller's RNG
  stream or leaves anything behind in the global environment.
