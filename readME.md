# lame

## Overview <img src="https://github.com/netify-dev/lame/blob/main/man/figures/lame_hex.png" align = "right" alt="hex" width="200px">

The `lame` R package provides easy way to model longitudinal network data via the Additive and Multiplicative Effects Model (`amen`). It handles networks with changing actor compositions, as well as C++ integration. This package enables researchers to efficiently model networks with temporal dependencies.

## Installation
    if (!require(devtools)) {
        install.packages("devtools")
      }
      library(devtools)
      install_github("netify-dev/lame")
      
## Usage
See our `lame-overview` vignette for detailed examples and documentation. To get started, supply data to the `lame` function. For example, we use the code below:

    library(lame)
    
    # load data
    data("vignette_data")
    
    # run model
    fit <- lame(
      Y = Y,                    
      Xdyad = Xdyad,           # dyadic covariates
      Xrow = Xrow,             # sender covariates
      Xcol = Xcol,             # receiver covariates 
      family = "bin",          # Binary  model
      rvar = TRUE,             # sender random effects
      cvar = TRUE,             # receiver random effects
      dcor = TRUE,            # Dyadic correlation
      R = 2,                   # Multiplicative effects dimension
      symmetric = FALSE,       # Directed networks
      burn = 100,              # Burn-in iterations
      nscan = 500,             # Post-burn-in iterations
      odens = 25,              # Output density
      print = FALSE,           # Suppress iteration output
      plot = FALSE             # Suppress real-time plots
    )
    
    # summarize model results
    summary(fit)
    
    # analyze model convergence
    paramPlot(fit$BETA)

## Contributors 
Cassy Dorff (Vanderbilt University), Shahryar Minhas (Michigan State University), Tosin Salau (Michigan State University)
