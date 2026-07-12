# cran-comments

## Submission of lame 1.3.1 (new submission)

## Test environments

* local: Ubuntu 24.04 (WSL2), R 4.3.3 — `R CMD check --as-cran`
* (please also see win-builder / R-hub results if requested)

## R CMD check results

0 errors | 0 warnings | 4 notes

* NOTE: New submission.
* NOTE: installed size (sub-directory `libs`). The package ships extensive
  C++ MCMC samplers built with Rcpp/RcppArmadillo; the templated Armadillo
  code produces a large shared object (~30Mb unstripped locally; ~1Mb once
  the installer strips debug symbols, as CRAN's builders do -- the remaining
  size is R code, data, and vignettes). All compiled code is required at
  runtime.
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
* `ame()`/`lame()` and the ALS fitters accept an optional `seed` argument;
  the user's `.Random.seed` is saved and restored on exit so fitting a
  model never perturbs the caller's RNG stream.
