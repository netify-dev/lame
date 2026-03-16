# Dynamic Effects in Longitudinal AME Models

## Introduction

Networks change. Countries that were close allies a decade ago may have
drifted apart. A legislator’s co-sponsorship patterns shift as they gain
seniority or switch committee assignments. A user’s purchasing habits
evolve as their tastes change.

The standard AME model assumes that each actor has a fixed latent
position and fixed additive effects across all time periods. This is a
reasonable starting point (pooling across time gives you more data and
more precise estimates), but it can miss important dynamics. The `lame`
package provides two mechanisms for letting the model capture temporal
change: `dynamic_uv` for time-varying latent positions and `dynamic_ab`
for time-varying additive effects.

This vignette explains what these options do, when to use them, and how
to interpret the results.

## What Are Dynamic Effects?

### The Static Baseline

In the standard AME model for longitudinal data, the tie between actors
$i$ and $j$ at time $t$ is:

$$y_{ij,t} = \beta\prime x_{ij,t} + a_{i} + b_{j} + u_{i}\prime v_{j} + \epsilon_{ij,t}$$

The covariates ($x_{ij,t}$) can vary over time, but everything else (the
sender effect $a_{i}$, the receiver effect $b_{j}$, and the latent
positions $u_{i}$ and $v_{j}$) is constant. This means the model assumes
that a country’s tendency to sanction others, or its position in the
latent “sanctioning space,” is the same in 1993 as in 2000.

### Dynamic Latent Positions (`dynamic_uv = TRUE`)

When you set `dynamic_uv = TRUE`, each actor’s latent position evolves
over time according to an AR(1) process:

$$U_{i,k,t} = \rho_{uv}\, U_{i,k,t - 1} + \epsilon_{i,k,t}$$$$V_{j,k,t} = \rho_{uv}\, V_{j,k,t - 1} + \eta_{j,k,t}$$

The autoregressive parameter $\rho_{uv}$ controls how persistent the
positions are. A value close to 1 means positions change slowly, as last
year’s position is a strong predictor of this year’s. A value close to 0
means positions are essentially re-drawn each period. In practice,
$\rho_{uv}$ is estimated from the data, so you don’t need to choose it
yourself.

This is useful when you believe the underlying community structure is
evolving: alliances shift, social groups re-form, trading blocs realign.

### Dynamic Additive Effects (`dynamic_ab = TRUE`)

When you set `dynamic_ab = TRUE`, the sender and receiver effects evolve
over time:

$$a_{i,t} = \rho_{ab}\, a_{i,t - 1} + \epsilon_{i,t}$$$$b_{j,t} = \rho_{ab}\, b_{j,t - 1} + \eta_{j,t}$$

This captures changes in actors’ overall activity levels. A country
might become a more active sanctioner after a change in government. A
student might become more or less socially active across semesters.
These are changes in how much an actor participates, not in who they
connect with (that’s what `dynamic_uv` captures).

You can use either option alone or combine them. In practice,
`dynamic_ab` is often a good place to start, since changes in overall
activity are common and relatively easy to estimate. `dynamic_uv` adds
more flexibility but also more parameters.

## Prior Specifications

The AR(1) coefficients and innovation variances have sensible default
priors:

- $\rho_{uv} \sim \text{TruncNormal}(0.9,0.1,0,1)$, centered on high
  persistence
- $\rho_{ab} \sim \text{TruncNormal}(0.8,0.15,0,1)$, slightly less
  persistent
- $\sigma_{uv}^{2},\sigma_{ab}^{2} \sim \text{InverseGamma}(2,1)$

These can be customized via the `prior` argument if you have strong
beliefs about the rate of change:

``` r
prior_custom <- list(
  rho_uv_mean = 0.95,    # expect very slow change in latent positions
  rho_uv_sd = 0.05,      # tight prior
  sigma_uv_shape = 3,
  sigma_uv_scale = 2
)
```

## Practical Example

We simulate a simple longitudinal network and fit models with and
without dynamic effects to see the difference.

``` r
library(lame)

set.seed(6886)
n <- 25  # actors
n_periods <- 5

# Generate sparse binary networks
Y_list <- list()
for(t in 1:n_periods) {
  Y_t <- matrix(rbinom(n*n, 1, 0.1), n, n)
  diag(Y_t) <- NA
  rownames(Y_t) <- colnames(Y_t) <- paste0("Actor", 1:n)
  Y_list[[t]] <- Y_t
}

# Fit model with both dynamic effects
# Note: iterations are reduced for vignette speed.
# Use burn >= 1000 and nscan >= 5000 for real analyses.
fit_dynamic <- lame(
  Y = Y_list,
  R = 2,
  dynamic_uv = TRUE,
  dynamic_ab = TRUE,
  family = "binary",
  burn = 100,
  nscan = 500,
  odens = 25,
  verbose = FALSE,
  plot = FALSE
)
#> Note: Square matrix assumed to be unipartite. Use mode='bipartite' for square bipartite networks.
```

### Visualizing Dynamic Effects

The trajectory plot shows how each actor’s latent position moves through
the 2-dimensional space over time. Each line traces one actor’s path
from the first to the last period. Actors whose lines are short stayed
in roughly the same position; actors with long, wandering lines
experienced substantial shifts in their network role.

``` r
uv_plot(fit_dynamic, plot_type = "trajectory")
```

![](dynamic_effects_files/figure-html/unnamed-chunk-3-1.png)

The `ab_plot` with `plot_type = "trajectory"` shows how each actor’s
additive effect changes over time. Rising lines indicate actors who
became more active (or more popular) over the observation window;
falling lines indicate the opposite. This is particularly useful for
identifying actors who experienced a structural break in their network
behavior.

``` r
ab_plot(fit_dynamic, effect = "sender", plot_type = "trajectory")
#> ℹ Showing top 5 and bottom 5 actors by average effect
#> → Use `show_actors` to specify actors to display
```

![](dynamic_effects_files/figure-html/unnamed-chunk-4-1.png)

As always, check convergence. Dynamic models have more parameters, so
adequate mixing is harder to achieve and requires longer chains.

``` r
trace_plot(fit_dynamic)
```

![](dynamic_effects_files/figure-html/unnamed-chunk-5-1.png)

### Comparing Model Specifications

A natural workflow is to fit several specifications and compare them:

``` r
# Static model (no dynamic effects)
fit_static <- lame(Y_list, R = 2, family = "binary",
                   burn = 100, nscan = 400, odens = 20,
                   verbose = FALSE, plot = FALSE)

# Dynamic latent positions only
fit_uv <- lame(Y_list, R = 2, dynamic_uv = TRUE, family = "binary",
               burn = 100, nscan = 400, odens = 20,
               verbose = FALSE, plot = FALSE)

# Dynamic additive effects only
fit_ab <- lame(Y_list, R = 2, dynamic_ab = TRUE, family = "binary",
               burn = 100, nscan = 400, odens = 20,
               verbose = FALSE, plot = FALSE)

# Full dynamic
fit_full <- lame(Y_list, R = 2, dynamic_uv = TRUE, dynamic_ab = TRUE,
                 family = "binary",
                 burn = 100, nscan = 400, odens = 20,
                 verbose = FALSE, plot = FALSE)
```

Compare the GOF statistics across models. The model that best reproduces
the observed network statistics (heterogeneity, transitivity, dyadic
dependence) without unnecessary complexity is usually the best choice.

``` r
if(!is.null(fit_static$GOF) && !is.null(fit_full$GOF)) {
  cat("Static model GOF (sample):\n")
  print(head(fit_static$GOF))
  cat("\nFull dynamic model GOF (sample):\n")
  print(head(fit_full$GOF))
}
#> Static model GOF (sample):
#> $sd.rowmean
#>             obs          1          2          3          4          5
#> [1,] 0.04659041 0.05576140 0.05874238 0.06140575 0.05450382 0.09871170
#> [2,] 0.06641285 0.06997142 0.07756288 0.08380931 0.05331666 0.06705222
#> [3,] 0.05033223 0.06997142 0.05033223 0.04369592 0.05964338 0.07302967
#> [4,] 0.06332982 0.06294972 0.06031031 0.05416026 0.05656854 0.05805744
#> [5,] 0.05941941 0.06248200 0.07240626 0.05656854 0.07578918 0.06780364
#>               6          7          8          9         10         11
#> [1,] 0.06645299 0.05276362 0.04827698 0.05479659 0.06974238 0.04760952
#> [2,] 0.06862458 0.06057502 0.06916647 0.05425864 0.06332982 0.07578918
#> [3,] 0.07529498 0.07578918 0.06008882 0.04898979 0.05381450 0.05537749
#> [4,] 0.06273755 0.06760671 0.05787343 0.06740920 0.07393691 0.06031031
#> [5,] 0.04085748 0.05523284 0.08333067 0.07393691 0.08089499 0.05874238
#>              12         13         14         15         16         17
#> [1,] 0.05331666 0.05291503 0.06740920 0.06218253 0.05689757 0.07718376
#> [2,] 0.06544209 0.08174758 0.06000000 0.07129282 0.05537749 0.06997142
#> [3,] 0.06218253 0.06760671 0.06605049 0.07185170 0.05479659 0.05381450
#> [4,] 0.06248200 0.05741080 0.05713143 0.06183850 0.05833238 0.05856051
#> [5,] 0.05810336 0.06531973 0.05736433 0.07016172 0.05828665 0.05276362
#>              18         19         20
#> [1,] 0.05571355 0.07440430 0.04898979
#> [2,] 0.07571878 0.05986652 0.05713143
#> [3,] 0.07422488 0.07302967 0.07454752
#> [4,] 0.06218253 0.05689757 0.07050296
#> [5,] 0.07118052 0.05356616 0.06544209
#> 
#> $sd.colmean
#>             obs          1          2          3          4          5
#> [1,] 0.06955094 0.06035451 0.06311894 0.05919459 0.05070174 0.04332820
#> [2,] 0.05547372 0.05741080 0.05901412 0.06601010 0.06764614 0.08142072
#> [3,] 0.06429101 0.06503332 0.05773503 0.05205126 0.04855238 0.08000000
#> [4,] 0.05043808 0.05503938 0.06140575 0.07023769 0.05656854 0.05200000
#> [5,] 0.04393935 0.05919459 0.06358197 0.06110101 0.06740920 0.08039900
#>               6          7          8          9         10         11
#> [1,] 0.06337192 0.07012370 0.04963869 0.06881860 0.05351635 0.05773503
#> [2,] 0.04052982 0.08207314 0.06819580 0.07031358 0.04630335 0.06540133
#> [3,] 0.05717808 0.06540133 0.06332982 0.04760952 0.07547185 0.06831301
#> [4,] 0.07529498 0.06858571 0.05670979 0.06226824 0.05773503 0.05571355
#> [5,] 0.04400000 0.06096994 0.07310267 0.06633250 0.08412689 0.07012370
#>              12         13         14         15         16         17
#> [1,] 0.04942334 0.06733003 0.08647158 0.06000000 0.06458070 0.07543651
#> [2,] 0.07313914 0.05787343 0.06218253 0.05901412 0.05656854 0.07807689
#> [3,] 0.05416026 0.06858571 0.07458329 0.08223543 0.04969239 0.05503938
#> [4,] 0.07418895 0.07892613 0.06374951 0.05619015 0.06379133 0.05503938
#> [5,] 0.06252466 0.06110101 0.05964338 0.06209670 0.05225578 0.07472617
#>              18         19         20
#> [1,] 0.04936936 0.06685307 0.06429101
#> [2,] 0.06531973 0.05642694 0.07436845
#> [3,] 0.06564551 0.05416026 0.04855238
#> [4,] 0.06831301 0.05805744 0.06031031
#> [5,] 0.06928203 0.04969239 0.06939741
#> 
#> $dyad.dep
#>               obs           1           2           3             4
#> [1,]  0.019269395 -0.08131488 0.008055339  0.13695707  2.206228e-02
#> [2,] -0.009636809  0.04974761 0.099728860  0.03073677 -5.142499e-02
#> [3,]  0.021291209  0.04974761 0.078171091  0.03103972  1.081274e-01
#> [4,]  0.146174863  0.04974761 0.019269395 -0.06932153 -2.018887e-16
#> [5,] -0.003450161  0.09772784 0.089561985  0.06692407  3.717591e-02
#>                 5          6           7            8            9          10
#> [1,] -0.009636809 0.14054890 0.051336261  0.092084751 0.0154744255 0.003507653
#> [2,] -0.019766758 0.05672532 0.042929560  0.008055339 0.0008429252 0.025178253
#> [3,]  0.075757576 0.01547443 0.037175910  0.073508894 0.0231259968 0.080538169
#> [4,]  0.136957068 0.01547443 0.019269395 -0.002692947 0.0008429252 0.115044248
#> [5,]  0.049639476 0.02601295 0.008055339 -0.044451872 0.0044247788 0.073508894
#>                 11          12          13            14           15
#> [1,] -1.946482e-16  0.05672532 0.158653846  0.0008429252  0.227272727
#> [2,] -4.445187e-02  0.06442411 0.068927978  0.0181992337 -0.112099644
#> [3,]  1.427432e-01 -0.01304945 0.006229799 -0.0642826735  0.113973145
#> [4,]  1.926939e-02  0.05849862 0.013668218 -0.0461501674  0.108282646
#> [5,]  5.133626e-02  0.07227207 0.023125997  0.0675877278 -0.006846557
#>               16           17          18           19           20
#> [1,]  0.01926939  0.030736767  0.07095916  0.042929560 -0.069321534
#> [2,]  0.01819923  0.047103194  0.10287081  0.008055339  0.124649860
#> [3,] -0.07388316 -0.053201733  0.07227207  0.021291209 -0.032514831
#> [4,]  0.13130096  0.087757708 -0.01304945 -0.019959829 -0.098418278
#> [5,]  0.16734987  0.008055339 -0.05661882 -0.023134421 -0.002692947
#> 
#> $cycle.dep
#>              obs             1            2             3            4
#> [1,] 0.004000134  0.0009460156  0.008101547 -0.0051947200 -0.013828250
#> [2,] 0.021198238  0.0255683605 -0.018621227  0.0087498469 -0.005484807
#> [3,] 0.002352264 -0.0108900015  0.009409970  0.0009646276 -0.008276965
#> [4,] 0.006223826  0.0089952373  0.007683338  0.0265516983 -0.001790639
#> [5,] 0.006189746 -0.0118244863  0.012019258  0.0028009721  0.009436416
#>                 5            6            7            8             9
#> [1,]  0.017513991 -0.020749485 -0.014606455  0.001772131 -0.0018047684
#> [2,]  0.009494850 -0.018688376  0.017700353 -0.026890506  0.0282077891
#> [3,] -0.018531061 -0.001830455 -0.035246722 -0.003425668 -0.0009210738
#> [4,] -0.001511517 -0.011347488 -0.001511517  0.004261137  0.0166002566
#> [5,] -0.004380404 -0.001319378  0.034413021 -0.007208952  0.0271510595
#>                10           11            12           13           14
#> [1,]  0.016857230 -0.017138977 -0.0019615453  0.024977095  0.020516333
#> [2,]  0.006912160 -0.003810647 -0.0053466624  0.000589853 -0.025551703
#> [3,] -0.035570531 -0.004358251  0.0255157815  0.011403254 -0.012019132
#> [4,] -0.032485374  0.035530987  0.0001985422 -0.007158847 -0.009241299
#> [5,]  0.009084672  0.023223107  0.0195900289 -0.009075946  0.004837166
#>                15           16            17            18           19
#> [1,] -0.005135113  0.008301590 -0.0001026881 -0.0156038561  0.004586374
#> [2,]  0.020719552 -0.016058536  0.0131795710  0.0159952321 -0.003178302
#> [3,] -0.011319684  0.008810302  0.0055513704 -0.0011332106  0.009786137
#> [4,] -0.022624671 -0.005490852 -0.0287553516 -0.0003411681 -0.010969457
#> [5,]  0.019064862 -0.014574196  0.0105659490 -0.0016624258 -0.006929324
#>                20
#> [1,] -0.004974697
#> [2,] -0.007717614
#> [3,]  0.010442032
#> [4,]  0.022652929
#> [5,]  0.008643732
#> 
#> $trans.dep
#>              obs           1            2           3           4           5
#> [1,] -0.03405087 -0.01725131 -0.015921262 -0.02439561 -0.01790421 -0.01257219
#> [2,] -0.02658332 -0.02213112 -0.036695001 -0.01839281 -0.02778143 -0.03720588
#> [3,] -0.02865812 -0.02199730 -0.008510928 -0.02102015 -0.04129031 -0.01012138
#> [4,] -0.02792264 -0.02067580 -0.020594018 -0.03995741 -0.02622008 -0.02148851
#> [5,] -0.03012973 -0.03257758 -0.029288661 -0.04035354 -0.02066110 -0.04298702
#>                6            7           8           9          10          11
#> [1,] -0.03013803 -0.023489311 -0.02458471 -0.03134053 -0.03089647 -0.00908110
#> [2,] -0.03266918 -0.009742238 -0.03737341 -0.02292007 -0.04269136 -0.03696153
#> [3,] -0.03039011 -0.032159196 -0.03006828 -0.02911492 -0.03697473 -0.03039543
#> [4,] -0.02529981 -0.041496198 -0.02472288 -0.02136849 -0.03302480 -0.01707744
#> [5,] -0.02057081 -0.030282635 -0.02545055 -0.04445262 -0.02781321 -0.02200062
#>               12          13          14          15          16          17
#> [1,] -0.01246107 -0.01222818 -0.04026493 -0.03415967 -0.02038793 -0.04622922
#> [2,] -0.02104329 -0.03827604 -0.03052009 -0.02027961 -0.03735162 -0.03767165
#> [3,] -0.02589286 -0.01979562 -0.02528432 -0.02751107 -0.01690447 -0.03127169
#> [4,] -0.03319584 -0.01564613 -0.01126292 -0.03214361 -0.02891739 -0.02347771
#> [5,] -0.02187535 -0.01610756 -0.02600661 -0.02684629 -0.04581472 -0.02869283
#>               18          19          20
#> [1,] -0.02612325 -0.04022310 -0.02069794
#> [2,] -0.03762923 -0.04281308 -0.03501448
#> [3,] -0.02943985 -0.01831534 -0.02481649
#> [4,] -0.03609199 -0.01838410 -0.02550496
#> [5,] -0.01493937 -0.02426377 -0.03944244
#> 
#> 
#> Full dynamic model GOF (sample):
#> $sd.rowmean
#>             obs          1          2          3          4          5
#> [1,] 0.04659041 0.05736433 0.05896892 0.04636090 0.05499091 0.05642694
#> [2,] 0.06641285 0.05406169 0.09318083 0.05887841 0.08554141 0.08460102
#> [3,] 0.05033223 0.07364781 0.06831301 0.04898979 0.06641285 0.06740920
#> [4,] 0.06332982 0.07050296 0.06858571 0.07571878 0.08830251 0.08223543
#> [5,] 0.05941941 0.05773503 0.05805744 0.04000000 0.06661331 0.08000000
#>               6          7          8          9         10         11
#> [1,] 0.07016172 0.04963869 0.02026491 0.09729680 0.05225578 0.08304216
#> [2,] 0.06877984 0.06429101 0.03525148 0.07418895 0.03555278 0.07564831
#> [3,] 0.07735632 0.07393691 0.04687572 0.09443163 0.04543127 0.05623759
#> [4,] 0.07436845 0.04898979 0.05163978 0.07125541 0.04270831 0.05547372
#> [5,] 0.06601010 0.04601449 0.02612789 0.08092795 0.03946306 0.03843609
#>              12         13         14         15         16         17
#> [1,] 0.11530828 0.06862458 0.07635007 0.04543127 0.09572182 0.05656854
#> [2,] 0.09965273 0.05782733 0.07635007 0.06337192 0.08349052 0.06166577
#> [3,] 0.07110556 0.06075086 0.08009994 0.04163332 0.07292005 0.08320256
#> [4,] 0.06661331 0.06916647 0.03448671 0.04749737 0.07236942 0.05600000
#> [5,] 0.08979978 0.05075431 0.05887841 0.04543127 0.05986652 0.08082904
#>              18         19         20
#> [1,] 0.05225578 0.06008882 0.05163978
#> [2,] 0.05600000 0.08326664 0.03525148
#> [3,] 0.03815757 0.06544209 0.04472136
#> [4,] 0.04548993 0.07302967 0.04630335
#> [5,] 0.02508652 0.05764258 0.03158058
#> 
#> $sd.colmean
#>             obs          1          2          3          4          5
#> [1,] 0.06955094 0.08845338 0.05896892 0.06843001 0.05964338 0.08708616
#> [2,] 0.05547372 0.06725077 0.08092795 0.04320494 0.06621178 0.05964338
#> [3,] 0.06429101 0.09429033 0.08246211 0.04163332 0.08252676 0.06839103
#> [4,] 0.05043808 0.05070174 0.05919459 0.05163978 0.06269503 0.07368401
#> [5,] 0.04393935 0.06928203 0.05070174 0.04472136 0.05919459 0.06928203
#>               6          7          8          9         10         11
#> [1,] 0.08634813 0.06374951 0.03282276 0.08793937 0.04079216 0.07547185
#> [2,] 0.07346655 0.06218253 0.03709447 0.09748846 0.03158058 0.07203703
#> [3,] 0.08475848 0.06928203 0.03555278 0.07821338 0.03912374 0.06803920
#> [4,] 0.06374951 0.06110101 0.04163332 0.08569714 0.04270831 0.05547372
#> [5,] 0.06075086 0.07200000 0.02856571 0.08726970 0.04715930 0.07922962
#>              12         13         14         15         16         17
#> [1,] 0.09360912 0.07943131 0.05968808 0.04393935 0.09846827 0.05291503
#> [2,] 0.07346655 0.06226824 0.06079474 0.05787343 0.05070174 0.05600000
#> [3,] 0.08787870 0.08460102 0.05552177 0.03651484 0.07012370 0.06823489
#> [4,] 0.08428523 0.05986652 0.06725077 0.04150502 0.08662563 0.05600000
#> [5,] 0.06053098 0.04519587 0.04898979 0.05351635 0.06205374 0.07302967
#>              18         19         20
#> [1,] 0.03362539 0.07310267 0.05033223
#> [2,] 0.05717808 0.06633250 0.04805552
#> [3,] 0.04460194 0.06544209 0.04898979
#> [4,] 0.05230679 0.07211103 0.04484046
#> [5,] 0.03600000 0.05991105 0.04687572
#> 
#> $dyad.dep
#>               obs           1          2           3          4          5
#> [1,]  0.019269395  0.02704806 0.09753216  0.05198566 0.15785326 0.04548721
#> [2,] -0.009636809  0.11900926 0.17999831  0.06174334 0.08521711 0.23204433
#> [3,]  0.021291209  0.12607413 0.10606061  0.06174334 0.16443850 0.16443850
#> [4,]  0.146174863 -0.02683461 0.13695707 -0.04347826 0.16115425 0.20178799
#> [5,] -0.003450161  0.18261563 0.11985605 -0.04166667 0.05849862 0.22566372
#>              6          7           8         9        10        11        12
#> [1,] 0.2089786 0.25274988  0.12321721 0.3102797 0.2776422 0.2352557 0.3982670
#> [2,] 0.2678688 0.11504425  0.52887080 0.5001272 0.0368563 0.2149730 0.3360054
#> [3,] 0.1125448 0.11123853  0.13338880 0.3680763 0.0368563 0.2358734 0.2501781
#> [4,] 0.2029234 0.19299451  0.45833333 0.3580203 0.1683059 0.3733289 0.3989662
#> [5,] 0.1451059 0.08640996 -0.01957586 0.3122566 0.2608798 0.1992536 0.3255506
#>                13         14          15         16         17          18
#> [1,] 0.2779284710 0.18084310  0.03685630 0.10779835 0.10287081 -0.03993344
#> [2,] 0.0008429252 0.08775771  0.32753519 0.10397769 0.09269212 -0.02393511
#> [3,] 0.1260741343 0.14054890  0.12500000 0.03727665 0.05387632 -0.03820598
#> [4,] 0.2546213476 0.03510538 -0.04515050 0.02332011 0.12590905  0.01190978
#> [5,] 0.0107119193 0.17391304 -0.04340568 0.01745541 0.19299451 -0.02796053
#>              19         20
#> [1,] 0.10984188 0.14529915
#> [2,] 0.07575758 0.02787748
#> [3,] 0.06892798 0.03846154
#> [4,] 0.04129794 0.05678174
#> [5,] 0.11900926 0.11711827
#> 
#> $cycle.dep
#>              obs            1            2            3            4
#> [1,] 0.004000134 -0.004863698 -0.010484886 -0.021129984  0.018973901
#> [2,] 0.021198238 -0.017274191 -0.005873411  0.010925973  0.009328688
#> [3,] 0.002352264 -0.028981326 -0.003795519  0.002227218 -0.005053393
#> [4,] 0.006223826 -0.021882546  0.004000134  0.010871739 -0.007014093
#> [5,] 0.006189746 -0.026868396 -0.024463582 -0.002375735 -0.003839827
#>                 5            6            7            8            9
#> [1,] -0.008536230 -0.039303728 -0.008285261 -0.006319956 -0.031421404
#> [2,] -0.015537100 -0.016013806 -0.009290098 -0.018414400 -0.001847216
#> [3,] -0.012795810  0.005121326 -0.012259145  0.012855259 -0.005000272
#> [4,]  0.008669044  0.001299266 -0.018818113 -0.020702835 -0.007237220
#> [5,] -0.004375336  0.010796637  0.013911535 -0.004497010 -0.015735101
#>                10           11           12           13            14
#> [1,] -0.003976914  0.007972305 -0.041457924 -0.016025802 -0.0002344697
#> [2,] -0.013985926 -0.037935778 -0.005536017  0.008369461 -0.0172758605
#> [3,]  0.007071495 -0.004166583 -0.029591872 -0.037431037  0.0365569627
#> [4,]  0.012392833 -0.005757249 -0.049784475  0.018714811  0.0198213194
#> [5,]  0.012672939 -0.024970320 -0.008342543  0.002984650 -0.0094648083
#>                15           16           17            18           19
#> [1,] -0.001975112 -0.018128598 -0.002134195  0.0003731671 -0.002534585
#> [2,] -0.011133850 -0.033652564 -0.028814645  0.0084323291 -0.006474708
#> [3,] -0.012557457 -0.005626893  0.002051449  0.0001126648 -0.017307656
#> [4,]  0.014364305  0.017188891  0.006429071  0.0035152534 -0.009469906
#> [5,]  0.004030296  0.007977255 -0.030992427 -0.0015654845 -0.020329116
#>                20
#> [1,] -0.013969712
#> [2,]  0.003129132
#> [3,] -0.009789250
#> [4,] -0.016958608
#> [5,]  0.015078705
#> 
#> $trans.dep
#>              obs           1            2            3           4            5
#> [1,] -0.03405087 -0.01061236 -0.021724260 -0.017865778 -0.03260994 -0.030969234
#> [2,] -0.02658332 -0.03694693 -0.044817618 -0.034753000 -0.03164718 -0.038850210
#> [3,] -0.02865812 -0.03296481 -0.008037569 -0.027104818 -0.04159892 -0.021102777
#> [4,] -0.02792264 -0.03245686 -0.011065922 -0.010488031 -0.02633641 -0.003111548
#> [5,] -0.03012973 -0.03095706 -0.026873769 -0.007805987 -0.04010184 -0.024793572
#>                6           7           8           9            10          11
#> [1,] -0.04255652 -0.02656948 -0.01426403 -0.04868229 -0.0293330781 -0.04475317
#> [2,] -0.05054492 -0.01242675 -0.02487811 -0.03806017 -0.0059530524 -0.04797191
#> [3,] -0.03308699 -0.02729886 -0.02795989 -0.03181522 -0.0216466701 -0.02575960
#> [4,] -0.03805555 -0.01813578 -0.02647248 -0.03324686 -0.0067048086 -0.01712892
#> [5,] -0.04215480 -0.01692475 -0.01398437 -0.04468652  0.0009406761 -0.04749371
#>               12          13          14           15          16           17
#> [1,] -0.07415339 -0.03666564 -0.02796993 -0.026664008 -0.03025114 -0.028845334
#> [2,] -0.04006404 -0.02959147 -0.02296751 -0.018021200 -0.04610224 -0.029388322
#> [3,] -0.06497210 -0.02797427 -0.01721970 -0.016969537 -0.03071219 -0.029924585
#> [4,] -0.08435361 -0.02350339 -0.02408716 -0.005985249 -0.02527722 -0.005199903
#> [5,] -0.03329674 -0.02908207 -0.03632440 -0.021313036 -0.01503230 -0.024312715
#>                18          19           20
#> [1,] -0.025195687 -0.03453933 -0.038286069
#> [2,] -0.009011132 -0.03081068 -0.030140114
#> [3,] -0.013557125 -0.01512401 -0.003065672
#> [4,] -0.011282700 -0.03158633 -0.028775805
#> [5,] -0.015283647 -0.03594802 -0.016975798
```

### When to Use Dynamic Effects

**Use `dynamic_ab`** when you suspect actors’ overall activity levels
change over time. This is common in many settings: countries go through
isolationist vs. interventionist phases, users churn in and out of
platforms, students become more or less engaged across semesters.

**Use `dynamic_uv`** when you suspect the underlying community structure
is shifting. This is a stronger claim, not just that actors are more or
less active, but that the pattern of *who connects with whom* is
changing. Examples include political realignment, market disruption, or
generational turnover in a social network.

**Use both** when you believe both types of change are happening. This
is the most flexible specification but also the most data-hungry. With
short panels (few time periods) or sparse networks, the dynamic
parameters may not be well-identified, and you might be better off with
the static model.

**Stick with static** when you have few time periods, when the network
structure is genuinely stable, or when you primarily care about the
covariate effects (which are the same across specifications) rather than
the latent structure.

## Dealing with Rotational Indeterminacy

One subtlety of dynamic latent space models: the latent space is only
identified up to rotation at each time point. This means that even if
actor positions are evolving smoothly, the raw estimated $U_{t}$
matrices might appear to “jump” between time periods due to arbitrary
rotations.

The
[`procrustes_align()`](https://netify-dev.github.io/lame/reference/procrustes_align.md)
function solves this by applying Procrustes rotation to align each time
period’s positions to the previous one:

``` r
# Align latent positions across time
aligned <- procrustes_align(fit_dynamic)
str(aligned$U)  # 3D array: actors x dimensions x time
#>  num [1:25, 1:2, 1:5] -0.3423 -0.0519 -0.1383 -0.0499 -0.0464 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ : chr [1:25] "Actor1" "Actor2" "Actor3" "Actor4" ...
#>   ..$ : NULL
#>   ..$ : NULL
```

You can also extract aligned positions as a tidy data frame using
[`latent_positions()`](https://netify-dev.github.io/lame/reference/latent_positions.md)
with `align = TRUE` (the default for `lame` objects):

``` r
lp <- latent_positions(fit_dynamic, align = TRUE)
head(lp)
#>    actor dimension time       value posterior_sd type
#> 1 Actor1         1    1 -0.34234484           NA    U
#> 2 Actor2         1    1 -0.05188764           NA    U
#> 3 Actor3         1    1 -0.13831092           NA    U
#> 4 Actor4         1    1 -0.04987621           NA    U
#> 5 Actor5         1    1 -0.04644628           NA    U
#> 6 Actor6         1    1  0.19183921           NA    U
```

This data frame is ready for `ggplot2` if you want to build custom
trajectory plots (e.g., highlighting specific actors or overlaying
external events).

## Implementation Notes

The dynamic effects are implemented in C++ via Rcpp and RcppArmadillo.
The AR(1) updates use Gibbs sampling with block updates across all
actors at each time point, which is efficient for large networks.

## References

1.  **Hoff, PD (2021)**. Additive and Multiplicative Effects Network
    Models. *Statistical Science* 36, 34–50.

2.  **Sewell, D. K., & Chen, Y. (2015)**. Latent space models for
    dynamic networks. *Journal of the American Statistical Association*,
    110(512), 1646-1657.

3.  **Durante, D., & Dunson, D. B. (2014)**. Nonparametric Bayes dynamic
    modeling of relational data. *Biometrika*, 101(4), 883-898.
