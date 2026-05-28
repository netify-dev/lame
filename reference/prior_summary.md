# Print the priors used by an AME / LAME / ame_als fit

Bayesian-hygiene helper modelled on `rstanarm::prior_summary()`. Prints
the priors that were actually in effect for a fit, with defaults filled
in. For `ame_als` fits it just states that the estimator is a
frequentist point estimator and there are no priors.

## Usage

``` r
prior_summary(object, ...)

# Default S3 method
prior_summary(object, ...)

# S3 method for class 'ame'
prior_summary(object, ...)

# S3 method for class 'lame'
prior_summary(object, ...)

# S3 method for class 'ame_als'
prior_summary(object, ...)
```

## Arguments

- object:

  a fitted `ame`, `lame`, or `ame_als` object.

- ...:

  ignored.

## Value

`object`, invisibly.
