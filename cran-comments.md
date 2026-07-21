# cran-comments

## Resubmission of lame, now 1.3.2 (new package)

This resubmission addresses the Rd markup issues reported by CRAN's incoming
checks on the previous submission (1.3.1): the affected man pages
(`glance.ame.Rd`, `procrustes_align.Rd`, `rUV_dynamic_t_fc_cpp.Rd`,
`rZ_ord_sym_fc.Rd`, `sample_rho_uv.Rd`) have been corrected and the package
re-documented. The version was bumped to 1.3.2; it also folds in bug fixes
made since that submission (see NEWS.md). The vignettes were compacted so
the suite rebuilds quickly on CRAN's builders.

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
  tests behind `skip_on_cran()`; CRAN-visible test time is ~20s.
* Vignettes use compact live examples. Rank-selection and PSIS-LOO examples
  retain their runnable code but are not evaluated during the build because
  those workflows require longer converged chains. Rebuilding all vignette
  outputs takes about 45 seconds on the local check environment.
* `ame()`/`lame()` and the ALS fitters accept an optional `seed` argument;
  the user's `.Random.seed` is saved and restored on exit so fitting a
  model never perturbs the caller's RNG stream.
