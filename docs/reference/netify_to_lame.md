# Using `lame` with `netify` data

For data prepared with the `netify` package, load `netify` and convert
to `lame`-ready input with `netify::to_amen(netlet, lame = TRUE)`. Pass
the result's `Y`, `Xdyad`, `Xrow`, `Xcol` components straight to
[`lame()`](https://netify-dev.github.io/lame/reference/lame.md) or
[`ame()`](https://netify-dev.github.io/lame/reference/ame.md).

## Details

    library(netify); library(lame)
    amen_in <- netify::to_amen(netlet, lame = TRUE)
    fit <- lame(
        Y      = amen_in$Y,
        Xdyad  = amen_in$Xdyad,
        family = "binary", R = 2,
        nscan  = 1000, burn = 250, odens = 5
    )

For longitudinal models, set `lame = TRUE` on the converter so the
dyadic covariates carry the period axis.
