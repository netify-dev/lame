# Visualize sender and receiver random effects

Creates a visualization of the additive sender (row) and receiver
(column) random effects from an AME or LAME model. Automatically detects
whether effects are static or dynamic and provides appropriate
visualization options.

## Usage

``` r
ab_plot(
  fit,
  effect = c("sender", "receiver"),
  sorted = TRUE,
  labels = NULL,
  title = NULL,
  time_point = NULL,
  plot_type = c("snapshot", "trajectory", "faceted", "ribbon"),
  show_actors = NULL
)
```

## Arguments

- fit:

  An object of class "ame" or "lame" from fitting an AME model

- effect:

  Character string specifying which effect to plot: "sender" (default)
  or "receiver"

- sorted:

  Logical; if TRUE (default), actors are sorted by effect magnitude

- labels:

  Logical; if TRUE, actor labels are shown on x-axis (default TRUE for n
  \<= 50 actors)

- title:

  Optional title for the plot

- time_point:

  For dynamic effects, which time point to plot (default: last). Can be
  numeric index, "all" for faceted plot, or "average" for time-averaged

- plot_type:

  For dynamic effects: "snapshot" (single time), "trajectory" (evolution
  over time), "faceted" (grid of time points), or "ribbon" (confidence
  bands over time). For static effects, this parameter is ignored.

- show_actors:

  Character vector of specific actors to highlight (for dynamic
  trajectory/ribbon plots)

## Value

A ggplot2 object that can be further customized

## Details

The additive effects in AME models represent:

- Sender effects (a):

  Actor-specific tendencies to form outgoing ties. Positive values
  indicate actors who send more ties than expected; negative values
  indicate actors who send fewer ties.

- Receiver effects (b):

  Actor-specific tendencies to receive incoming ties. Positive values
  indicate actors who receive more ties than expected; negative values
  indicate actors who receive fewer ties.

For static effects, the plot displays these effects as a dot plot with
vertical lines extending from zero to each effect estimate.

For dynamic effects (when fit contains a_dynamic/b_dynamic), additional
options are available to visualize how effects evolve over time.

## Author

Cassy Dorff, Shahryar Minhas, Tosin Salau

## Examples

``` r
# \donttest{
# Fit an AME model
data(YX_nrm)
fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X,
           nscan = 100, burn = 10, odens = 1, verbose = FALSE)

# Visualize sender effects
ab_plot(fit, effect = "sender")


# Visualize receiver effects without sorting
ab_plot(fit, effect = "receiver", sorted = FALSE)

# }
```
