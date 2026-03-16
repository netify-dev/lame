# MCMC trace plots and density plots for AME/LAME model parameters

Creates diagnostic plots for Markov Chain Monte Carlo (MCMC) samples
from AME or LAME models. Displays trace plots to assess convergence and
mixing, alongside density plots to visualize posterior distributions.

## Usage

``` r
trace_plot(
  fit,
  params = c("all", "beta", "variance"),
  include = NULL,
  exclude = NULL,
  ncol = 3,
  nrow = NULL,
  burn.in = 0,
  thin = 1,
  title = NULL
)
```

## Arguments

- fit:

  An object of class "ame" or "lame" containing MCMC samples

- params:

  Character vector specifying which parameters to plot: "beta" for
  regression coefficients, "variance" for variance components, or "all"
  (default) for both

- include:

  Character vector of specific parameter names to include

- exclude:

  Character vector of specific parameter names to exclude

- ncol:

  Number of columns for plot layout (default 3)

- nrow:

  Number of rows for plot layout (default NULL, determined
  automatically)

- burn.in:

  Number of initial iterations to exclude as burn-in when calculating
  statistics (default 0, assumes burn-in already removed)

- thin:

  Thinning interval for display (default 1, no thinning)

- title:

  Optional title for the plot

## Value

A ggplot2 object that can be further customized

## Details

This function produces two types of diagnostic plots:

- Trace plots:

  Show the evolution of parameter values across MCMC iterations. Good
  mixing is indicated by rapid exploration of the parameter space with
  no trends or stuck periods.

- Density plots:

  Show the posterior distribution of parameters. Multiple modes may
  indicate identification issues or convergence problems.

The plots help diagnose:

- Convergence: Has the chain reached the stationary distribution?

- Mixing: Is the chain exploring the parameter space efficiently?

- Autocorrelation: Are successive samples highly correlated?

Parameters displayed include:

- Regression coefficients (beta)

- Variance components (va, vb, cab, rho, ve)

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Fit an AME model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Basic trace plots for all parameters
trace_plot(fit)


# Only regression coefficients
trace_plot(fit, params = "beta")


# Only variance components
trace_plot(fit, params = "variance")

# }
```
