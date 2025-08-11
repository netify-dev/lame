# lame

## Overview <img src="https://github.com/netify-dev/lame/blob/main/man/figures/lame_hex.png" align = "right" alt="hex" width="200px">

The `lame` R package provides a framework for modeling longitudinal network data via the Additive and Multiplicative Effects Model (`amen`). Key features include:

- **Dynamic effects modeling**: Both additive (sender/receiver) and multiplicative (latent factor) effects can evolve over time following AR(1) processes, capturing temporal heterogeneity in network structure
- **Temporal handling**: Manages networks with changing actor compositions across time periods
- **Performance**: Computations implemented in C++ via Rcpp/RcppArmadillo

## Installation

    if (!require(devtools)) {
        install.packages("devtools")
      }
      library(devtools)
      install_github("netify-dev/lame")

## Usage

See our `lame-overview` vignette for detailed examples and documentation. To get started, supply data to the `lame` function. For example, we use the code below:

```{r}
    library(lame)
    
    # load data
    data("vignette_data")
    
    # run model with static effects
    fit <- lame(
      Y = Y,                    
      Xdyad = Xdyad,           # dyadic covariates
      Xrow = Xrow,             # sender covariates
      Xcol = Xcol,             # receiver covariates 
      family = "binary",       # Binary model
      rvar = TRUE,             # sender random effects
      cvar = TRUE,             # receiver random effects
      dcor = TRUE,             # Dyadic correlation
      R = 2,                   # Multiplicative effects dimension
      symmetric = FALSE,       # Directed networks
      burn = 100,              # Burn-in iterations
      nscan = 500,             # Post-burn-in iterations
      odens = 25,              # Output density
      print = FALSE,           # Suppress iteration output
      plot = FALSE             # Suppress real-time plots
    )
    
    # run model with dynamic effects (AR(1) temporal evolution)
    fit_dynamic <- lame(
      Y = Y,                    
      Xdyad = Xdyad,           
      Xrow = Xrow,             
      Xcol = Xcol,             
      family = "binary",       
      dynamic_ab = TRUE,       # Time-varying sender/receiver effects
      dynamic_uv = TRUE,       # Time-varying latent positions
      R = 2,                   
      symmetric = FALSE,       
      burn = 1000,             # Longer burn-in recommended for dynamic models
      nscan = 5000,            # More samples for temporal correlation
      odens = 25,              
      print = FALSE,           
      plot = FALSE,
      prior = list(            # Customize temporal persistence
        rho_uv_mean = 0.9,     # High persistence for latent factors
        rho_ab_mean = 0.8      # Moderate persistence for additive effects
      )
    )
    
    # summarize model results
    summary(fit)
    
    # analyze model convergence
    trace_plot(list(BETA = fit$BETA, VC = fit$VC))
```

## Visualization

The package provides unified plotting functions that automatically detect whether effects are static or dynamic:

```{r}
    # Visualize additive effects (sender/receiver)
    ab_plot(fit, effect = "sender")                    # Static effects
    ab_plot(fit_dynamic, plot_type = "trajectory")     # Dynamic effects over time
    ab_plot(fit_dynamic, time_point = "average")       # Time-averaged dynamic effects
    
    # Visualize multiplicative effects (latent factors)
    uv_plot(fit)                                       # Static latent positions
    uv_plot(fit_dynamic, plot_type = "trajectory")     # Dynamic trajectories
    uv_plot(fit_dynamic, plot_type = "faceted")        # Multiple time points
```

## Key Improvements over amen

### Dynamic Effects Modeling

- **Time-varying latent positions** (`dynamic_uv`): Captures evolving community structure, homophily patterns, and network reconfiguration via AR(1) processes on multiplicative effects
- **Time-varying heterogeneity** (`dynamic_ab`): Models changing individual activity levels and popularity through AR(1) evolution of additive effects
- **Priors**: Customizable temporal persistence (ρ) and innovation variance (σ²) parameters with defaults
- **Theoretical foundation**: Based on Sewell & Chen (2015, JASA) and Durante & Dunson (2014, Biometrika)

### Performance & Implementation

- **C++ acceleration**: Functions implemented in C++ via Rcpp/RcppArmadillo
- **Memory efficiency**: ~70% reduction in memory usage for dynamic effects
- **Sampling**: Forward-filtering backward-sampling (FFBS) for temporal updates
- **Processing**: Vectorized operations for actor updates

### User Interface

- **Temporal flexibility**: Handles networks with actors entering/exiting over time
- **CLI**: Console output with progress bars via cli, crayon, and progress packages
- **Plotting**: Various additional plotting functions
- **Diagnostics**: Convergence monitoring 

## Documentation

For information on the dynamic effects implementation, see:

- `vignette("dynamic_effects")` - Guide to temporal modeling
- `?lame` - Function documentation with mathematical details
- Package website with examples and tutorials (coming soon)

## Citation

If you use the dynamic effects features in your research, please cite:

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

- Sewell, D. K., & Chen, Y. (2015). Latent space models for dynamic networks. *JASA*, 110(512), 1646-1657.
- Durante, D., & Dunson, D. B. (2014). Nonparametric Bayes dynamic modeling of relational data. *Biometrika*, 101(4), 883-898.

## Contributors

- **Cassy Dorff** (Vanderbilt University)
- **Shahryar Minhas** (Michigan State University)
- **Tosin Salau** (Michigan State University)