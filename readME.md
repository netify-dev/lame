# **lame** <img src="man/figures/lame_hex.png" align="right" alt="hex" width="200px">

> **L**ongitudinal **A**dditive and **M**ultiplicative **E**ffects Models for Networks

## Overview

The `lame` package extends the Additive and Multiplicative Effects (AME) framework for network analysis, providing support for both cross-sectional and longitudinal networks. The package includes two main functions: `ame` for cross-sectional network analysis and `lame` for longitudinal network analysis with dynamic effects that capture temporal heterogeneity through autoregressive processes. Both functions support unipartite (square) and bipartite (rectangular) network structures, making them suitable for analyzing various types of relational data.

## Installation

You have two options for installing `lame`:

### 🔧 Option 1: Install from GitHub (requires build tools)

> ⚠️ **Requires R build tools:**
> - **macOS**: Xcode Command Line Tools
> - **Windows**: Rtools
> - **Linux**: build-essential and related packages

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("netify-dev/lame", dependencies = TRUE)
```

### 📦 Option 2: Install from CRAN (coming soon ... maybe ... probably)

```r
# Once available on CRAN
install.packages("lame")
```

#### First, install dependencies (if needed):

```r
# Install required packages if you don't already have them
deps <- c("Rcpp", "RcppArmadillo", "ggplot2", "plyr", "reshape2", 
          "ggrepel", "network", "gridExtra", "coda", "cli", 
          "rlang", "lifecycle", "purrr")

# Check which packages are not installed
missing_deps <- deps[!deps %in% installed.packages()[,"Package"]]

# Install missing packages
if(length(missing_deps) > 0) {
  install.packages(missing_deps, repos='https://cloud.r-project.org/')
}
```

## Quick Start

```r
library(lame)

# Load example data
data("vignette_data")

# Fit a basic longitudinal AME model
fit <- lame(
  Y = Y,                    # List of T network matrices
  Xdyad = Xdyad,           # Dyadic covariates
  Xrow = Xrow,             # Sender covariates
  Xcol = Xcol,             # Receiver covariates
  family = "binary",       # Network type
  R = 2,                   # Latent dimensions
  burn = 100,              # Burn-in iterations
  nscan = 500              # Post-burn samples
)

# Summary and diagnostics
summary(fit)
```

### Dynamic Effects Models

The key innovation in `lame` is the ability to model time-varying network effects:

```r
# Fit model with dynamic effects
fit_dynamic <- lame(
  Y = Y,
  Xdyad = Xdyad,
  Xrow = Xrow,
  Xcol = Xcol,
  family = "binary",
  dynamic_ab = TRUE,       # Time-varying sender/receiver effects
  dynamic_uv = TRUE,       # Time-varying latent positions
  R = 2,
  burn = 1000,            # Longer burn-in for dynamic models
  nscan = 5000,           # More samples for temporal correlation
  prior = list(           # Customize temporal persistence
    rho_uv_mean = 0.9,    # High persistence for latent factors
    rho_ab_mean = 0.8     # Moderate persistence for additive effects
  )
)
```

## Visualization

Plotting functions:

```r
# Visualize additive effects (sender/receiver)
ab_plot(fit, effect = "sender")                    # Static effects
ab_plot(fit_dynamic, plot_type = "trajectory")     # Dynamic effects over time
ab_plot(fit_dynamic, time_point = "average")       # Time-averaged dynamic effects

# Visualize multiplicative effects (latent factors)
uv_plot(fit)                                       # Static latent positions
uv_plot(fit_dynamic, plot_type = "trajectory")     # Dynamic trajectories
uv_plot(fit_dynamic, plot_type = "faceted")        # Multiple time points

# Model diagnostics
trace_plot(list(BETA = fit$BETA, VC = fit$VC))    # MCMC convergence
gof_plot(fit)                                       # Goodness-of-fit
```

## Key Features

### 📈 Dynamic Effects Modeling

- **Time-varying latent positions** (`dynamic_uv`): Captures evolving community structure and homophily patterns via AR(1) processes
- **Time-varying heterogeneity** (`dynamic_ab`): Models changing individual activity levels and popularity over time
- **Bipartite network support**: Handles rectangular adjacency matrices with separate latent factors for row and column nodes
- **Flexible priors**: Customizable temporal persistence (ρ) and innovation variance (σ²) parameters
- **Theoretical foundation**: Based on Sewell & Chen (2015, JASA) and Durante & Dunson (2014, Biometrika)

## Documentation

- **Getting Started**: `vignette("lame-overview")`
- **Dynamic Effects**: `vignette("dynamic_effects")`
- **Bipartite Networks**: `vignette("bipartite")`

## Citation

If you use `lame` in your research, please cite:

```bibtex
@Manual{lame2024,
  title = {lame: Longitudinal Additive and Multiplicative Effects Models for Networks},
  author = {Shahryar Minhas and Tosin Salau and Cassy Dorff},
  year = {2024},
  note = {R package version 0.10},
  url = {https://github.com/netify-dev/lame},
}
```

The dynamic effects implementation draws on:

- Sewell, D. K., & Chen, Y. (2015). Latent space models for dynamic networks. *Journal of the American Statistical Association*, 110(512), 1646-1657.
- Durante, D., & Dunson, D. B. (2014). Nonparametric Bayes dynamic modeling of relational data. *Biometrika*, 101(4), 883-898.

## Contributors

- **Shahryar Minhas** (Michigan State University)
- **Cassy Dorff** (Vanderbilt University)
- **Tosin Salau** (Michigan State University)

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/netify-dev/lame/issues)
- 💬 **Questions**: [GitHub Discussions](https://github.com/netify-dev/lame/discussions)
- 📧 **Contact**: <minhassh@msu.edu>
