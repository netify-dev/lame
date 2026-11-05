#' AME model fitting routine
#'
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to cross-sectional relational data of various types.
#' This function supports both unipartite (square) and bipartite (rectangular)
#' networks. For longitudinal networks, use the \code{lame} function. The
#' additive and multiplicative effects framework is due to Hoff (2005, 2021).
#'
#' @details
#' This command provides posterior inference for parameters in AME models of
#' cross-sectional relational data, assuming one of six possible data types/models.
#' The function supports both unipartite networks (square adjacency matrices) and
#' bipartite networks (rectangular adjacency matrices with distinct row and column
#' node sets) for single time point analysis.
#'
#' \strong{Model.} For a dyad \eqn{(i, j)} the AME linear predictor is
#'
#' \deqn{\eta_{ij} = \beta_0 + x_{ij}'\beta + a_i + b_j + u_i' D v_j,}
#'
#' where \eqn{\beta} are regression coefficients on dyadic / nodal
#' covariates, \eqn{a_i} is a row (sender) random effect, \eqn{b_j} is
#' a column (receiver) random effect, and \eqn{u_i' D v_j} is the
#' multiplicative latent-factor term (rank \code{R}). The observation
#' model is \eqn{Y_{ij} \sim F(\eta_{ij}, \theta)} with \eqn{F}
#' specified by \code{family} (Gaussian for \code{"normal"}, probit
#' for \code{"binary"} / \code{"cbin"}, etc.). For unipartite
#' (\code{symmetric = FALSE}) the residual error has dyad-level
#' correlation \eqn{\rho} between \eqn{(i, j)} and \eqn{(j, i)}; for
#' bipartite, dyad correlation is fixed at 0.
#'
#' \strong{Priors (in brief).} \eqn{\beta} has a Zellner-style g-prior
#' (\code{g}); \eqn{(a_i, b_i)} are jointly Normal with covariance
#' \eqn{\Sigma_{ab}} (Inverse-Wishart prior \code{Sab0 / eta0});
#' \eqn{u_i, v_j} are independent Normal with covariance
#' \eqn{\Sigma_{uv}} (Inverse-Wishart prior with scale \code{kappa0 * Suv0}
#' and \code{kappa0} degrees of freedom);
#' the dyad-correlation \eqn{\rho} has the prior \eqn{2\,\mathrm{Beta}(2, 2) - 1}
#' on \eqn{(-1, 1)} (density proportional to \eqn{1-\rho^2}, the prior
#' \code{\link{lame}} uses), updated with Metropolis steps. See
#' \code{prior_summary(fit)} for the priors actually used.
#'
#' \strong{Choosing R.} The multiplicative rank \code{R} controls the
#' dimensionality of latent homophily / heterogeneity not explained
#' by covariates and additive effects. \code{R = 0} fits an
#' additive-only social-relations model; \code{R = 1} or \code{2} is
#' typical for small / medium networks. A rank above \code{floor(n/3)} is
#' fitted like any other (on 18-actor panels with a rank-2 truth, ranks 6 to 17
#' ran without numerical failure; the chain mixed more slowly and the fitted
#' values lost accuracy gradually as the rank grew past the truth), so compare
#' ranks by their fit rather than by a cap. Latent factors capture
#' unobserved structure (clusters, hub patterns, transitive triangles
#' the covariates miss) and are accessed at \code{fit$U}, \code{fit$V}.
#'
#' \strong{Identifiability.} The latent factor term \eqn{u_i' D v_j}
#' is invariant to rotation and reflection of \eqn{U, V}; the package
#' canonicalises with an SVD so successive draws are interpretable.
#' For unipartite fits with both resolved actor margins enabled, reporting
#' centers static latent positions and transfers the removed row, column and
#' grand-mean components to the actor effects and intercept. A symmetric shared
#' actor absorbs half the grand constant in each margin when no intercept is
#' fitted. Restricted unipartite models retain sampled coordinates in their
#' stored factor draws. This reporting map
#' preserves each draw's predictor and leaves the sampler state unchanged.
#' For visual stability across posterior summaries see
#' \code{\link{procrustes_align}} and \code{\link{latent_positions}}.
#'
#' \strong{Theoretical Foundation:}
#'
#' The AME model decomposes network structure into several components:
#' \deqn{y_{ij} = \beta'x_{ij} + a_i + b_j + u_i'v_j + \epsilon_{ij}}
#' where:
#' \itemize{
#'   \item \eqn{\beta'x_{ij}}: Fixed effects of dyadic/nodal covariates
#'   \item \eqn{a_i}: Additive sender (row) effect for node i
#'   \item \eqn{b_j}: Additive receiver (column) effect for node j
#'   \item \eqn{u_i'v_j}: Multiplicative interaction between latent factors
#'   \item \eqn{\epsilon_{ij}}: Dyadic error term (may be correlated)
#' }
#'
#' This specification generalizes the social relations model (Warner et al. 1979)
#' and latent space models (Hoff et al. 2002) within a unified framework.
#'
#' \strong{Prior Distributions:}
#'
#' The model uses conjugate and semi-conjugate priors where possible:
#' \itemize{
#'   \item Regression coefficients: \eqn{\beta \sim N(0, g\sigma^2(X'X)^{-1})} (g-prior)
#'   \item Additive effects: \eqn{(a_i, b_i)' \sim N(0, \Sigma_{ab})} jointly
#'   \item Covariance: \eqn{\Sigma_{ab} \sim IW(\eta_0, \eta_0 S_{ab0})} (inverse-Wishart)
#'   \item Multiplicative effects: Hierarchical shrinkage via \eqn{\eta_0}
#'   \item Dyadic correlation: \eqn{2\,\mathrm{Beta}(2, 2) - 1} prior on
#'         \eqn{(-1, 1)}, density \eqn{p(\rho) \propto 1-\rho^2} (the prior of
#'         \code{\link{lame}}), with Metropolis updates
#' }
#'
#' The inverse-Wishart prior on \eqn{\Sigma_{ab}} allows learning correlation between
#' sender and receiver effects, capturing reciprocity patterns.
#'
#' Note on the additive-effects variance scale: for a unipartite fit \code{Sab0}
#' defaults to a data-scaled matrix rather than a fixed \code{diag(2)}. For the
#' \code{"normal"} and \code{"poisson"} families the default is
#' \code{Sab0 = diag(2) * vscale}, where \code{vscale} is an empirical-Bayes
#' moment estimate of the sender/receiver variance: the mean of the variances of
#' the centred row means and column means of \code{Y} (of \code{log1p(Y)} for
#' \code{"poisson"}). The \code{"binary"}, \code{"ordinal"}, \code{"cbin"} and
#' \code{"frn"} families apply the same idea to a probit-moment residual, but
#' only when \code{start_vals} is not supplied and the design has at least one
#' column; otherwise they fall back to \code{Sab0 = diag(2)} with
#' \code{eta0 = 4}. A bipartite fit scales \code{Sab0} the same way (one scale
#' for both sides).
#' Call \code{prior_summary()} on a fitted object to see the prior actually used.
#'
#' The data scaling matters because the inverse-Wishart prior contributes
#' pseudo-data on the scale of \code{eta0 * Sab0} no matter what units \code{Y}
#' is in, so a fixed \code{Sab0 = diag(2)} pulls the additive-effects variances
#' \code{va}/\code{vb} upward whenever the true sender/receiver variance is well
#' below 1, most visibly at small \eqn{n}. The data-scaled default sets the
#' prior on the scale of \code{Y}; \code{prior_summary()} prints it.
#'
#' To use a fixed identity scale instead, pass
#' \code{prior = list(Sab0 = diag(2), eta0 = 4)}. Changing this prior can
#' affect other jointly estimated quantities as well as \code{va},
#' \code{vb} and \code{cab}. Matching these prior settings alone does not
#' establish equivalence to another package's posterior.
#'
#' \strong{Multiplicative Effects (Latent Factors):}
#'
#' When R > 0, the model includes R-dimensional latent factors:
#' \itemize{
#'   \item Asymmetric case: \eqn{u_i, v_j \in \mathbb{R}^R} with \eqn{u_i'v_j} interaction
#'   \item Symmetric case: \eqn{u_i = v_i} with eigendecomposition \eqn{ULU'}
#'   \item Captures homophily, transitivity, and community structure
#'   \item R chosen via model selection or set to 2-3 for visualization
#' }
#'
#' \strong{Estimation Algorithm:}
#'
#' The model uses a Gibbs sampler with the following updates:
#' 1. Sample latent Z given parameters (data augmentation for non-normal families)
#' 2. Update regression coefficients \eqn{\beta} via g-prior conjugate update
#' 3. Update additive effects (a,b) jointly with \eqn{\beta}
#' 4. Update covariance \eqn{\Sigma_{ab}} from inverse-Wishart
#' 5. Update multiplicative effects U,V via Gibbs or Metropolis-Hastings
#' 6. Update dyadic correlation \eqn{\rho} via Metropolis-Hastings
#' 7. Update variance \eqn{\sigma^2} (for continuous families)
#'
#' \strong{Standard Model Types:}
#'
#' The following data types/models are available:
#'
#' "normal": A normal AME model (identity link: \eqn{E[Y] = \eta}).
#'
#' "binary": A binary probit AME model (probit link: \eqn{P(Y=1) = \Phi(\eta)}).
#'
#' "ordinal": An ordinal probit AME model (cumulative probit link). An intercept
#' is not identifiable with data-induced cutpoints. Explicit cutpoints anchored
#' at zero support an intercept.
#'
#' "cbin": An AME model for censored binary data (probit link with censoring). The value of
#' 'odmax' specifies the maximum number of links each row may have.
#'
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#'
#' "poisson": An overdispersed Poisson-lognormal model for counts (negative
#' or non-integer observed values are an error):
#' \eqn{Y \sim \mathrm{Poisson}(e \exp(z))} with
#' \eqn{z \sim N(\eta, \sigma^2)}. Given \eqn{\eta} and \eqn{\sigma^2}, a
#' cell's mean count is \eqn{e \exp(\eta + \sigma^2/2)}. \code{YPM},
#' \code{fitted()} and \code{predict(type = "response")} report it at the
#' posterior mean \eqn{\eta} of the linear predictor and the posterior mean
#' \eqn{\sigma^2} of the residual variance, the form the fast estimators
#' report. Predictive quantiles and event probabilities come from posterior
#' predictive draws (\code{simulate()}, \code{predict(type = "distribution")}).
#' Zero exposure gives a structural zero count (cross-sectional fits use unit
#' exposure).
#'
#' @param Y For unipartite: an n x n square relational matrix. For bipartite:
#' an nA x nB rectangular relational matrix where nA is the number of row
#' nodes and nB is the number of column nodes. A cross-sectional
#' \pkg{netify} object is also accepted and converted with
#' \code{netify::to_lame()}; when \code{family} or \code{mode} is omitted,
#' the value inferred by \pkg{netify} is used. See family below for data
#' types.
#' @param Xdyad For unipartite: an n x n x pd array of dyadic covariates (e.g. distance, shared group). For
#' bipartite: an nA x nB x pd array. A 2-D matrix (single dyadic covariate) must be wrapped as
#' \code{array(x, dim = c(n, n, 1))}; \code{Inf}/\code{NaN} entries are rejected.
#' @param Xrow For unipartite: an n x pr matrix of \strong{sender} (row) covariates (e.g. sender's age, group).
#' For bipartite: an nA x pr matrix of row-node covariates. A \code{data.frame} is accepted and coerced
#' to numeric matrix internally.
#' @param Xcol For unipartite: an n x pc matrix of \strong{receiver} (column) covariates. For bipartite:
#' an nB x pc matrix of column-node covariates. A \code{data.frame} is accepted and coerced internally.
#' @param rvar logical: fit row random effects (asymmetric case)?
#' @param cvar logical: fit column random effects (asymmetric case)?
#' @param dcor logical: fit a dyadic correlation (asymmetric case)? Note: not used for bipartite networks. With \code{method = "als"} the argument goes to \code{\link{ame_als}}, whose directed binary probit, Poisson or normal fit then has the correlation in its likelihood.
#' @param nvar logical: fit nodal random effects (symmetric case)?
#' @param R integer: dimension of the multiplicative effects (can be zero). For bipartite networks,
#' this is used as the default for both R_row and R_col if they are not specified.
#' A binary fit whose multiplicative term reproduces the observed ties
#' exactly is separated: the size of the term is then set by its prior, the
#' fit says so, and the signs and the fitted probabilities are what the data
#' determine; a smaller \code{R} (down to zero for a panel whose ties barely
#' change over time) reproduces them without saturating.
#' @param R_row integer: for bipartite networks, dimension of row node multiplicative effects (defaults to R)
#' @param R_col integer: for bipartite networks, dimension of column node multiplicative effects (defaults to R)
#' @param mode character: either "unipartite" (default) for square networks or "bipartite" for rectangular networks.
#' Not all combinations of \code{family} and \code{mode} are supported -- see the table under \strong{Supported
#' family x mode combinations} below.
#' @param family character: one of "normal","binary","ordinal","cbin","frn","poisson". See the
#' \strong{Supported family x mode combinations} table below for which combinations are valid; see \strong{Details}
#' for the model assumptions behind each family.
#' @param intercept logical: fit model with an intercept?
#' @param symmetric logical: Is the sociomatrix symmetric by design?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' families. A cap below a row's observed out-degree contradicts \code{Y} and
#' is an error.
#'   A cap that most rows' observed outdegrees reach erases the cross-row
#'   spread in outdegree that the sender effects rest on: a row at its cap
#'   shows only that its level is high enough to fill the cap. The sampler
#'   pools the prior scale of the row side with the column side once
#'   `r round(100 * .lame_odmax_sat_threshold)` percent of the observed rows
#'   (`r round(100 * .lame_odmax_sat_threshold_single)` percent in a single
#'   period) nominate the full \code{odmax}, counting a row with hidden cells
#'   by what a full cap thins to at its observed share;
#'   \code{fit$odmax_saturated_share} records the share and
#'   \code{fit$odmax_sat_threshold} the threshold. \code{va}, \code{vb} and
#'   \code{cab} are reported at every share: on simulated panels of 24 to 30
#'   actors their 95 percent intervals covered the generating values in 0.88
#'   to 1.00 of the panels at every share, with and without a fifth of the
#'   cells hidden and at generating sender sds of 0.25, 0.5 and 1 (at 0.25
#'   the posterior mean of \code{va} runs up to 2.3 times the generating
#'   variance while the interval still covers). \code{va} and the level of the rows follow the design. A design of two periods
#'   or more with a cap share of 0.5 or more and a posterior-mean \code{va} from
#'   0.5 to 1.5 reports the intercept and the row-covariate coefficients as not
#'   estimated: the rows below the cap carry too little of the spread of the
#'   senders and the posterior of \code{va} sits near its prior. In simulated
#'   panels (24 actors, 4 and 6 periods, one chain at the package defaults,
#'   1440 fits at generating sender sds of 0.5, 1 and 1.5) the intercept was
#'   covered in 0.89 of the 312 panels in that region (0.76 at sd 0.5 and 1.5,
#'   0.97 at sd 1) and in 0.96 of the 1128 outside it. \code{va} is reported in
#'   the region: its interval covered 0.91 of the 312 panels (0.79 at sd 0.5,
#'   1.00 at sd 1, 0.72 at sd 1.5) and 0.99 outside it. The posterior mean
#'   of \code{va} ran 0.18 and 0.23 above the generating variance at sd 0.5
#'   (shares 0.3 and 0.7), 0.08 and 0.00 above it at sd 1 and 0.02 and 0.18
#'   below it at sd 1.5. In a single period (\code{ame()} and
#'   one-period \code{lame()} fits) a row at its cap leaves the level of the
#'   rows unbounded above: with every row at the cap the intercept interval
#'   covered 0.64 of the panels, and with cells hidden and a share from 0.5 it
#'   covered 0.84, against 0.98 or more elsewhere. The intercept and any
#'   row-covariate coefficient are reported as not estimated there, with the
#'   draws kept in \code{fit$BETA}. On 200 simulated 170 by 100 bipartite panels
#'   (cap 10, shares 0.3 and 0.7) the intercept was covered in 1.00 of the
#'   panels, \code{vb} in 0.98 to 1.00 and \code{va} in 0.90 to 1.00, its
#'   posterior mean 0.09 to 0.15 below the generating variance; nothing is withheld there.
#'   When the effects prior is supplied (\code{prior$Sab0}, \code{eta0} or
#'   \code{etaab}), the fit names the prior in the summary
#'   (\code{fit$prior_supplied} records it); a prior on the effects does not
#'   fix the level of a row, so the intercept stays withheld in those designs, and
#'   \code{va} is reported under that prior.
#'   A cap the design itself imposes cannot be relaxed: read \code{APM} /
#'   \code{BPM} as relative sender and receiver levels there.
#' @param prior a list containing hyperparameters for the prior distributions.
#' Available options and their defaults:
#' \describe{
#'   \item{Sab0}{Prior scale matrix for the additive-effects covariance. A 2x2
#'         matrix where Sab0\\[1,1\\] is the prior variance for row effects,
#'         Sab0\\[2,2\\] is the prior variance for column effects, and off-diagonals
#'         control correlation between row and column effects. For a unipartite
#'         fit this defaults to \code{diag(2)} scaled to the observed
#'         sender/receiver heterogeneity in \code{Y} (unconditionally for the
#'         continuous families, conditionally for the discrete ones); a
#'         bipartite fit is scaled the same way, with one scale for both
#'         sides (a single period is too little to read the sender side of a
#'         capped \code{"cbin"} or \code{"frn"} panel on its own, so the
#'         scale is read over every row, not only the rows below the cap as
#'         for panels; see \code{Sab0} in \code{\link{lame}}). See the note on the
#'         additive-effects variance scale in Details. Pass \code{diag(2)}
#'         explicitly for a fixed unit-scale prior.}
#'   \item{eta0}{Prior degrees of freedom for the additive-effects covariance
#'         \eqn{\Sigma_{ab}} (default: round(4 + 3 \\* n/100) for a continuous-family
#'         fit, unipartite or bipartite (n is then the number of row actors,
#'         as in \code{\link{lame}}), and round(4 \\* vdfmlt) for a
#'         discrete-family fit, where n is the number of actors and vdfmlt
#'         is a probit-moment variance multiplier estimated from \code{Y}).
#'         Higher values impose stronger shrinkage of the row/column effects
#'         toward the prior scale. The multiplicative-effects degrees of
#'         freedom are controlled by \code{kappa0}, not \code{eta0}.}
#'   \item{etaab}{Prior degrees of freedom for covariance of additive effects
#'         (default: 4 + 3 \\* n/100). Controls shrinkage of row/column random effects.
#'         Larger values shrink effects toward zero. \emph{Used by the bipartite
#'         and longitudinal paths; ignored for a unipartite cross-sectional
#'         \code{ame()} fit, which controls the additive prior through
#'         \code{Sab0} and \code{eta0}.}}
#'   \item{s20}{Fixed positive prior scale for the dyadic residual variance
#'         \eqn{v_e}. The prior has one degree of freedom. By default, normal
#'         fits use the variance of the observed responses, computed once before
#'         sampling; constant responses use their mean square, or 1 if all zero.
#'         Poisson fits use 1 on the latent log-rate scale. The resolved scale is
#'         reported as \code{fit$s2_prior_scale}. It is never recomputed from the
#'         current residuals. A supplied value overrides this default. In a
#'         bipartite cross-sectional \code{ame()} fit, \code{s20} also
#'         multiplies the coefficient prior variance (\code{g * s20}); its
#'         default there is the variance of the response for the normal family
#'         (so that the prior is the same in any unit of the response) and 1 on
#'         the latent scale of the other families; elsewhere that prior is
#'         controlled separately by \code{g}.}
#'   \item{s2u0}{Prior variance for multiplicative effects (default: 1).
#'         \emph{Used by the bipartite and longitudinal paths.}}
#'   \item{Suv0}{Inverse-Wishart (inverse) scale matrix for the
#'         multiplicative-effects covariance. The prior scale is
#'         \code{kappa0 * Suv0} and the prior mean is
#'         \code{kappa0 * Suv0 / (kappa0 - 2 * R - 1)}. When not supplied it
#'         is \code{s * diag(2 * R)}, \eqn{s} the scale of the additive effects
#'         estimated from the data, the one the default \code{Sab0} reads: for
#'         \code{"binary"}, the variance of the actor means of the
#'         standardised residuals \eqn{(y - \bar y) / \phi(\Phi^{-1}(\bar y))}
#'         of the less frequent value; for \code{"ordinal"}, \code{"cbin"} and
#'         \code{"frn"}, the same probit moments of the ranks; for
#'         \code{"poisson"}, the variance of the centred log row and column
#'         rates; for \code{"normal"}, the square root of the mean variance of
#'         the row and column means of the centred response (in the unit of the
#'         response, as the product \eqn{u_i'v_j} is). Every period of a
#'         \code{\link{lame}} panel enters, and a single period gives the value of
#'         \code{ame()}.
#'         \code{\link{lame}} has the same default. A bipartite network has no covariance of this kind
#'         (the rows of its factors have a unit-free N(0, I) prior and the interaction matrix carries the
#'         scale), and supplying \code{Suv0} or \code{kappa0} for one is an error of class
#'         \code{"lame_prior_unused"}.}
#'   \item{kappa0}{Prior degrees of freedom for the multiplicative-effects
#'         covariance \eqn{\Sigma_{uv}}. The default when R > 0 is
#'         \code{round((2 * R + 2) * m)}: \eqn{2R + 2}, the smallest value for
#'         which the prior has a finite mean, times
#'         \eqn{m = 0.25 / \overline{p(1 - p)}}, the inverse Bernoulli
#'         information of the cells relative to a cell at \eqn{p = 1/2}, with
#'         \eqn{p} the probit fit of
#'         the additive moments, for the \code{"binary"}, \code{"cbin"} and
#'         \code{"frn"} families (\eqn{m \ge 1}), and \eqn{m = 1} for the
#'         others. A sparse network has little information per cell about its
#'         factors, and the prior then holds \eqn{\Sigma_{uv}} near the data
#'         scale \eqn{s I} with more pseudo-observations; \eqn{m} depends on
#'         the density, not on the number of actors. Higher values impose
#'         stronger shrinkage of the latent factors toward the prior scale.}
#'   \item{actor_move}{Logical, default \code{TRUE}. For \code{"poisson"} and
#'         \code{"ordinal"} fits (unipartite, directed or symmetric, and
#'         bipartite): translate each actor's level together with the latent
#'         values of its row and column (an exact draw: the Gaussian density of
#'         the latent array is unchanged, so the conditional is the level prior
#'         times the Poisson factor, or the interval the ranks or cutpoints
#'         allow), draw the location of the latent array (ordinal), and, for the
#'         Poisson family, translate the coefficients and move the dyadic
#'         variance and correlation together with the latent residuals and
#'         redraw the positions with the latent values of their rows and columns. The
#'         augmented sweep moves an actor with no events, and the intercept and
#'         the spread of the levels with it, by steps of order
#'         \code{sqrt(s2 / number of cells)}, which takes thousands of sweeps on a
#'         sparse panel; the updates leave the posterior unchanged. For binary fits (unipartite,
#'         directed or symmetric, with any additive margin in the model, a fit
#'         with none drawing the positions alone, and bipartite; complete or with
#'         missing cells and actors):
#'         redraw each actor's level and static position from its exact
#'         conditional with the latent array integrated out (an elliptical
#'         slice step). The augmented sweep moves an actor whose cells hold no
#'         events by a step of order \code{1 / sqrt(number of cells)} however
#'         wide its posterior is, so such an actor relaxes slowly; the update
#'         leaves the posterior unchanged. A bipartite fit of any family with an
#'         intercept also draws, exactly, along the ridge that the intercept
#'         shares with the means of its (uncentred) additive effects (a shift of
#'         the intercept and the opposite shift of the levels changes no
#'         predictor, so the draw has no likelihood term). For directed
#'         \code{"cbin"} and \code{"frn"} fits (unipartite or bipartite), each
#'         receiver's level is translated together with its latent column, an
#'         exact draw from the prior of the level cut to the interval the order
#'         constraints of the column leave (the augmented sweep moves a receiver
#'         nobody nominates by the same slow step), the rows free above translate
#'         with the intercept against the rows with room, and a directed fit with
#'         \code{dcor = TRUE} updates the dyadic correlation with the
#'         non-nominations integrated out (see the Details of
#'         \code{\link{lame}}). \code{FALSE} omits them; the work is recorded in
#'         \code{fit$actor_move}.}
#'   \item{chain_extensions}{Number of times a chain that is still moving at
#'   the end of the run is continued before the call is refused (default 3;
#'   0 refuses at the first flag). See the convergence policy.}
#'   \item{actor_move_every}{Number of sweeps between visits of that update
#'         (default 10, and 2 for the Poisson translations; 1 for the translations of \code{"cbin"} and
#'         \code{"frn"} fits, whose reciprocity update runs on every tenth
#'         visit, and for the draw along the intercept ridge of bipartite fits
#'         of the other families).}
#'   \item{uv_g_scale_move}{Logical, default \code{TRUE}. For bipartite fits with
#'         positions: each sweep also draws, exactly, the factors by which the
#'         columns of the positions and the rows and columns of the interaction
#'         matrix are rescaled without changing their product, because the unit
#'         priors anchor the split of scale between them weakly and the
#'         conditional updates cross it slowly. \code{FALSE} omits it; its work is
#'         recorded in \code{fit$uv_g_scale_move}.}
#'   \item{uv_split_move}{Logical, default \code{TRUE}. For directed fits with latent factors
#'         (\code{R > 0}): draw the split of the scale between \eqn{U} and \eqn{V}
#'         exactly (\eqn{U \to cU}, \eqn{V \to V/c} leaves the product and the likelihood
#'         unchanged; with the covariance of the rows integrated out \eqn{\log c} has the
#'         matrix-t density along the orbit and is drawn by a slice step). The conditional
#'         draws moved the split by a few percent per sweep against a posterior sd of
#'         0.2 to 0.5. Recorded in \code{fit$uv_split_move}.}
#'   \item{sparse_position_move}{Logical, default \code{TRUE}. For directed \code{"cbin"},
#'         \code{"frn"} and \code{"ordinal"} fits with latent factors (\code{R > 0}): each
#'         sweep also draws the position of every sender and receiver together with the
#'         latent values of its row or column, exactly. Moving \eqn{u_i} changes the
#'         predictor of row \eqn{i} by a different amount in each column; carrying the
#'         latent values along changes no residual, so the conditional of the position is
#'         its prior times the indicator that the shifted row is still a latent row the
#'         data allow, drawn by an elliptical slice step. Actors the data bound on one
#'         side only (no nominations, no nominators, only the lowest category) move across
#'         their wide posterior in a few sweeps instead of by random walk. Not applied to
#'         \code{"cbin"} / \code{"frn"} panels with hidden cells. Recorded in
#'         \code{fit$sparse_position_move}.
#'         An ordinal fit that runs it visits the level translations of
#'         \code{actor_move} every sweep unless \code{actor_move_every} says otherwise.}
#'   \item{sparse_position_every}{Sweeps between visits of \code{sparse_position_move}
#'         (default 1: every sweep).}
#'   \item{cov_ridge_move}{Logical, default \code{TRUE}. For bipartite fits with a row
#'         or column covariate (a design slice that is constant along the columns or
#'         the rows): each sweep also draws, exactly, the amount by which the
#'         coefficient and the additive effects of that margin trade against one
#'         another without changing any cell's predictor. The two conditional
#'         draws cross this ridge slowly (a row coefficient had 18 to 44 effective
#'         draws of 1,600 in a nomination panel). Directed unipartite fits of \code{ame()} and
#'         \code{lame()} (static or dynamic additive effects, static coefficients) take the same
#'         draw. \code{FALSE} omits it.}
#'   \item{ordinal_scale_move}{Logical, default \code{TRUE}. For
#'         \code{"ordinal"} fits: each sweep also draws, exactly, the scale of the
#'         whole state (the coefficients, the additive effects and the latent
#'         array together, with the free cutpoints of
#'         \code{ordinal_cutpoints = "explicit"}) and, with the default
#'         data-induced cutpoints, the shift and the scale of the latent array
#'         alone, along transformations that keep every rank. The windows the ranks
#'         leave the latent values move the implied cutpoints by steps of
#'         \code{1 / number of cells}, while the scale of the configuration has
#'         a posterior of width \code{1 / sqrt(number of cells)} and the
#'         coefficients and the dyadic correlation follow it: on 150 actors a
#'         dyadic coefficient had 8 effective draws of 600 and the correlation 95,
#'         now 600 and 400. \code{FALSE} omits it; its work is recorded in
#'         \code{fit$ordinal_scale_move}. See \code{\link{lame}}.}
#'   \item{binary_scale_move}{Logical, default \code{TRUE}. For \code{"binary"}
#'         fits: every \code{binary_scale_every} sweeps the coefficients and the
#'         additive effects are multiplied by a factor drawn exactly from its
#'         conditional density, with the latent array integrated out. The
#'         augmented sweep crosses that direction by steps of
#'         \code{1 / sqrt(number of cells)}, and when a covariate separates the
#'         ties the posterior along it is tens of times wider: on 30 actors,
#'         after 8,000 sweeps, the separating coefficient of a directed panel had
#'         12 effective draws of 2,000, now 460, and that of a bipartite panel,
#'         which had drifted over those sweeps (14.0 and 15.6 for two seeds),
#'         now has a posterior mean of 59 at 2,000, 8,000 and 32,000 sweeps. \code{FALSE} omits it; its
#'         work is recorded in \code{fit$binary_scale_move}. See
#'         \code{\link{lame}}.}
#'   \item{binary_scale_every}{Positive whole number, default 5: the sweeps
#'         between the draws of that factor.}
#' }
#' Common usage: prior = list(Sab0 = diag(c(2, 2)), eta0 = 10) for moderate
#' shrinkage, or prior = list(Sab0 = diag(c(0.5, 0.5))) for tighter control.
#' For a unipartite cross-sectional \code{ame()} fit the prior is controlled by
#' \code{Sab0}, \code{eta0}, \code{Suv0} and \code{g}.
#' @param g optional \strong{scalar} regression prior scale, independent of
#' the residual variance. The unipartite Gaussian block uses precision
#' \code{(XtX + 1e-6 * diag(diag(XtX))) / g}; an all-ones first column
#' additionally receives intercept variance \code{100 * g / XtX[1,1]}.
#' If not specified, \code{g = n_observed_cells * var(Y)} for normal data
#' and \code{g = n_observed_cells} otherwise, where observed cells exclude
#' the diagonal (before covariate-missingness propagation). Per-coefficient (vector) \code{g} is not supported by the
#' unipartite path. A finite positive \code{g} is required when regression
#' coefficients are fitted; invalid scales are not replaced or clipped.
#' \emph{Note:} \code{g} is a top-level argument to \code{ame()},
#' \strong{not} an element of \code{prior = list(...)}; passing
#' \code{prior = list(g = 0.1)} is a no-op (warned about).
#' @param seed random seed for the MCMC sampler (default 6886). The sampler
#'   is seeded internally with this value, so results are reproducible by
#'   default and an external \code{set.seed()} call has \emph{no} effect on
#'   the chain -- pass a different \code{seed} here to vary the draws (e.g.
#'   when running multiple chains). The caller's \code{.Random.seed} is
#'   restored on exit, so fitting never perturbs your RNG stream.
#' @param nscan number of iterations of the Markov chain (beyond burn-in)
#' @param burn burn in for the Markov chain. Typical use is \code{burn} far
#' smaller than \code{nscan}; \code{burn > nscan} is allowed but warns.
#' @param odens output density (thinning interval) for the Markov chain.
#' \code{nscan / odens} samples are stored.
#' @param plot accepted for signature parity with \code{\link{lame}};
#'   ignored. \code{ame()} runs a single-period MCMC and does not draw
#'   live trace panels. Default \code{FALSE}.
#' @param verbose logical: print progress while running? Default TRUE.
#' @param print Deprecated. Use \code{verbose} instead.
#' @param gof logical: calculate goodness of fit statistics? Setting to TRUE
#'   adds approximately 2-5% to runtime. For faster sampling without GOF overhead,
#'   set gof=FALSE and use gof() after model fitting.
#' @param custom_gof optional function or list of named functions for computing
#'   custom goodness-of-fit statistics. Each function must accept a single matrix Y
#'   as input and return a numeric vector. If a single function is provided, it
#'   should return a named vector. If a list of functions is provided, each function
#'   should return a single value and will be named according to the list names.
#'   Custom statistics will be computed in addition to default statistics.
#'   Example: custom_gof = function(Y) c(density = mean(Y > 0, na.rm = TRUE))
#' @param start_vals List from previous model run containing parameter starting values for new MCMC
#' @param periodic_save logical: indicating whether to periodically save MCMC results
#' @param out_file character vector indicating name and path in which file should be stored if periodic_save is selected. For example, on an Apple OS out_file="~/Desktop/ameFit.rda".
#' @param save_interval quantile interval indicating when to save during the post-burn-in period.
#' @param model.name optional string for model selection output
#' @param posterior_opts optional list of posterior draw-storage options,
#'   usually built with \code{\link{posterior_options}}. Recognised names
#'   (unknown names trigger a warning): \code{save_UV}, \code{save_UV_draws},
#'   \code{save_ab}, \code{thin_UV}, \code{thin_ab}.
#'   \code{save_UV = TRUE} stores per-iteration latent-position draws on
#'   \code{fit$U_samples} / \code{fit$V_samples} (\code{[actor, dim, draw]}
#'   arrays; bipartite fits also store the interaction-matrix draws on
#'   \code{fit$G_samples}); \code{save_UV_draws} (the \code{\link{lame}}
#'   spelling) is accepted as an alias. \code{save_ab = TRUE} stores
#'   additive-effect draws on \code{fit$a_samples} / \code{fit$b_samples}
#'   (a symmetric fit has one effect per actor: the two hold the same draws).
#'   \code{thin_UV} / \code{thin_ab} thin the stored draws. Default
#'   \code{NULL}.
#' @param n_chains integer: the older name of \code{chains}.
#' @param cores integer or \code{NULL}: workers for the chains. The chains run
#'   forked on Unix, at most one worker per chain; \code{NULL} uses
#'   \code{getOption("lame.cores")} if set and otherwise all cores but one.
#'   Windows runs the chains one after the other.
#' @param chains integer: the number of independently seeded chains (default 2,
#'   or \code{getOption("lame.chains")}). Chain \eqn{k} runs with seed
#'   \code{seed + 1000 * (k - 1)}. With \code{chains > 1} the chains are checked,
#'   compared and pooled (see the convergence policy); \code{chains = 1} runs one
#'   chain and the drift check only.
#' @param use_sparse_matrices logical: use sparse matrix storage for large networks? (default: FALSE).
#'   Recommended only for truly sparse networks (< 10% non-zero entries).
#' @param method character: \code{"mcmc"} (default, the Bayesian MCMC fit) or
#'   \code{"als"} (the fast, MCMC-free iterative block coordinate descent point
#'   estimator). When \code{method = "als"}, MCMC-specific arguments (\code{nscan},
#'   \code{burn}, \code{odens}, ...) are warned about and ignored; the call
#'   forwards to \code{\link{ame_als}}, which uses \code{prior$s20}.
#' @param bootstrap integer (only used when \code{method = "als"}): number of
#'   bootstrap replicates. \code{0} (default) skips the bootstrap; \code{N >= 30}
#'   runs \code{N} replicates and attaches the result so that \code{\link{confint}}
#'   returns bootstrap intervals (a smaller positive \code{N} is an error).
#' @param bootstrap_type character (only used when \code{method = "als"}):
#'   \code{"parametric"} (default) or \code{"block"}.
#' @param bootstrap_block_length integer: block length for the block bootstrap.
#' @param bootstrap_seed optional integer seed for the bootstrap.
#' @param multistart character (only used when \code{method = "als"} and
#'   \code{R > 0}): \code{"none"} (default), \code{"cheap"} (4 starts), or
#'   \code{"full"} (8 starts). Forwarded to \code{\link{ame_als}}, which
#'   reruns the block-coordinate fit from several starts and keeps the one
#'   with the lowest objective.
#' @param bias_correct character or \code{NULL} (only used when
#'   \code{method = "als"}): how the ALS point estimate
#'   handles the finite-sample bias of a fixed-effects fit, forwarded to
#'   \code{\link{ame_als}}. \code{NULL} (default) fits a normal, binary or
#'   count network by the empirical-Bayes estimator (\code{"eb"}; Laplace standard
#'   errors from \code{vcov()}), directed, bipartite or undirected, with or
#'   without a latent rank (see \code{bias_correct} in
#'   \code{\link{ame_als}}). \code{"none"} returns the
#'   fixed-effects estimate at any size; \code{"jackknife"} and
#'   \code{"bootstrap"} require that correction of it.
#' @param ... reserved for future use. Passing \code{lame()}-only
#'   arguments (e.g. \code{dynamic_beta}, \code{period_exposure}) here
#'   triggers a clean abort directing you to \code{\link{lame}}; passing
#'   any other unrecognised name warns so typos are visible rather than
#'   silently dropped.
#' @param save_log_lik logical: attach a per-iteration matrix of matched
#'   observed marginal scores. LOO and WAIC additionally require a verified
#'   observed likelihood factorization; directed residual dyad correlation
#'   prevents cellwise factorization. Default \code{FALSE}.
#' @param log_lik_method \code{"observed_exact"} (default) supports normal,
#'   binary, Poisson-lognormal and explicit-cutpoint ordinal marginal scores.
#'   Poisson integration uses an adaptive log-integral with concavity tail
#'   bounds. Capped
#'   nominations, fixed ranks and data-induced ordinal models are refused.
#'   \code{"observed_approx"} retains the existing approximate diagnostic
#'   scores for these families without claiming an observed likelihood.
#'   AME does not save latent Z draws for augmented scoring. Inspect
#'   \code{fit$log_lik_contract}; approximate diagnostics cannot be used
#'   for observed-data LOO or WAIC.
#' @param ordinal_cutpoints character: cutpoint convention for
#'   \code{family = "ordinal"}. \code{"data_induced"} (default) uses the
#'   data-induced cutpoints; \code{"explicit"} samples ordered cutpoints with the first fixed at zero and unit latent residual
#'   variance. Independent directed and symmetric models use the Cowles (1996)
#'   marginal Metropolis-Hastings update. Directed models with \code{dcor = TRUE}
#'   use cutpoint Gibbs updates conditional on the correlated latent observations.
#'   Explicit cutpoints are supported for unipartite models; predictions retain
#'   the original numeric category labels and use each saved cutpoint draw.
#'   Ignored for other families.
#'
#' @section Supported family x mode combinations:
#' Every \code{family} is supported under
#' both \code{mode}s. The bipartite Z-samplers live in \code{R/rZ_bipartite.R}
#' and dispatch per family; see also the inline comment in \code{R/lame.R}
#' (the rectangular samplers live in \code{R/rZ_bipartite.R}).
#'
#' \tabular{lll}{
#'   \strong{family}  \tab \strong{unipartite} \tab \strong{bipartite} \cr
#'   normal           \tab yes                 \tab yes                \cr
#'   binary           \tab yes                 \tab yes                \cr
#'   ordinal          \tab yes                 \tab yes                \cr
#'   cbin             \tab yes                 \tab yes                \cr
#'   frn              \tab yes                 \tab yes                \cr
#'   poisson          \tab yes                 \tab yes                \cr
#' }
#'
#' Symmetric (\code{symmetric = TRUE}) fits require a symmetric \code{Y} and
#' symmetric \code{Xdyad} slices: the regression term of a symmetric fit is
#' symmetric, so an asymmetric slice is an error (class
#' \code{lame_symmetric_xdyad}) that names the slice. Replace it by
#' \code{(X + t(X)) / 2}. A value missing in one triangle of a slice is not an
#' asymmetry: it is filled from the other triangle (see \emph{Missing
#' covariates}). \code{\link{lame}} applies the same checks to every period,
#' and the ALS routes (\code{method = "als"}, \code{\link{ame_als}}) refuse an
#' asymmetric slice the same way.
#' \code{family = "ordinal"} with \code{symmetric = TRUE} is supported via the
#' unit-variance unordered-dyad samplers for either cutpoint convention,
#' which mirror upper-triangle draws to the lower
#' triangle so \eqn{Z = t(Z)} holds at every sweep. The rank-nomination
#' families (\code{"cbin"}, \code{"frn"}) model each actor's own ranked
#' nominations and are directed by construction, so \code{symmetric = TRUE}
#' is rejected for them.
#'
#' \strong{Symmetric input with one triangle missing.} When
#' \code{symmetric = TRUE} and one triangle of \code{Y} is fully \code{NA}
#' (the user stored only the lower or upper triangle), the symmetry check
#' \code{any(is.finite(Y - t(Y)))} is \code{FALSE}, and the call stops with an
#' error instead of mirroring the populated triangle. Symmetrize \code{Y}
#' explicitly (for example \code{Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]})
#' before calling; \code{lame()} applies the same check to every period.
#'
#' @section Data preparation:
#' \code{ame()} accepts a matrix directly. For long-format edgelists,
#' bipartite data, and covariates, use \pkg{netify} to build the network object
#' and pass that object as \code{Y}. \code{ame()} will call
#' \code{netify::to_lame()} internally. If you already have an \pkg{igraph} or
#' \pkg{network} object and only need a plain adjacency matrix,
#' \code{\link{as_lame_y}} is still available as a small convenience helper.
#'
#' Directed normal MCMC models with estimated reciprocity use the prior
#' \eqn{2\,\mathrm{Beta}(2, 2) - 1} on the full open interval \eqn{(-1,1)},
#' with density proportional to \eqn{1 - \rho^2}. If every observed reciprocal
#' response pair is exactly equal, the posterior can be improper at
#' \eqn{\rho = 1}. A conservative structural check counts reciprocal pairs and
#' bounds the antisymmetric actor/regression design dimension. With
#' multiplicative rank \code{R}, it also includes a conservative \code{n * R}
#' factor bound. The check refuses cases where the pair count exceeds that
#' combined bound by four or more (the posterior is improper from that excess
#' under this prior; under the arc-sine prior of versions before 1.4.0 it was
#' improper from an excess of one); it does not certify propriety of other
#' designs. For an undirected network,
#' use \code{symmetric = TRUE}; to fix reciprocity at zero in a directed
#' model, use \code{dcor = FALSE}. An exactly symmetric \code{Y} fitted with
#' \code{symmetric = FALSE} is allowed: the two cells of a dyad are equal, so
#' the residuals of reciprocal cells are equal and the posterior of the
#' dyadic correlation concentrates near one (with \code{dcor = FALSE} the two
#' cells count as two observations). A message says so and points at
#' \code{symmetric = TRUE}, which fits the undirected model.
#'
#' For an undirected/symmetric network, pass \code{symmetric = TRUE}; for a
#' rectangular two-mode network (students x courses, donors x candidates), pass
#' \code{mode = "bipartite"} or build the \pkg{netify} object with
#' \code{mode = "bipartite"}. \code{Xrow} and \code{Xcol} accept either a
#' numeric matrix or a data.frame (coerced internally); \code{Xdyad} must be a
#' 3-D array \code{n x n x p} of numeric covariates with no
#' \code{Inf}/\code{NaN}.
#'
#' @section Missing covariates:
#' A cell whose dyadic covariate, sender (row) covariate or receiver (column)
#' covariate is missing is an unobserved tie, for unipartite and bipartite
#' networks alike: it is set to \code{NA} in \code{Y} before the design is
#' built, the sampler augments it like any other missing cell, and the
#' covariate coefficients are estimated from the complete cells. A message
#' gives the number of such cells. A missing value is never modelled as 0,
#' which would bias the coefficient of a covariate that is not centred. In an
#' undirected network (\code{symmetric = TRUE}) a dyad has one covariate value:
#' a hole in one triangle of a dyadic covariate is filled from the other
#' triangle (a message gives the number of values), so the slice stays
#' symmetric and the dyad stays observed, and a dyad missing in both triangles
#' is unobserved. A missing node covariate removes the actor's row and column.
#' A covariate that is missing everywhere stops the fit. \code{\link{lame}}
#' (per period) and the ALS routes (\code{\link{ame_als}}) apply the same rule
#' and give the same messages.
#'
#' @section Collinear covariates:
#' A regression design whose columns are exactly aliased over the observed
#' cells (a covariate that duplicates another, or one that is a linear
#' combination of others, counting the intercept and the sender and receiver
#' covariates) is an error of class \code{lame_design_rank}. It names the
#' aliased columns and the combination each equals: the data cannot separate
#' their coefficients, so the posterior along that direction would be the
#' prior's alone. A column is aliased when it is a combination of the columns
#' before it to a relative 1e-8 (as \code{lm()} reports them); a design that is
#' only nearly collinear is fitted.
#'
#' @section When to use AME vs ERGM:
#' AME and ERGM are complementary tools for binary network analysis,
#' not direct substitutes. \strong{ERGM} is a class of exponential-family
#' models built around explicit network statistics (counts of edges,
#' mutual ties, triangles, geometrically-weighted shared partners, ...).
#' You write the statistics you think matter, ERGM gives you their
#' coefficients. ERGM excels when you have a substantive theory about
#' which configurations drive tie formation.
#'
#' \strong{AME} models latent homophily / heterogeneity directly via
#' sender, receiver, and multiplicative latent-factor effects. You don't
#' enumerate triadic terms; the multiplicative-effects rank \code{R}
#' captures higher-order structure (clustering, transitivity, hub
#' patterns) implicitly. AME excels when (a) you have dyadic / nodal
#' covariates whose effects you want to interpret cleanly without ERGM
#' degeneracy, (b) higher-order structure is "nuisance" that you want
#' to absorb but not parameterise, or (c) you need a posterior
#' distribution over predictions for forecasting or imputation.
#'
#' Practical guidance: if your research question is "do nodes that
#' share attribute X tend to form triangles together?", reach for
#' ERGM's \code{gwesp}. If your research question is "controlling for
#' unobserved sender / receiver heterogeneity and latent clustering,
#' what is the effect of dyadic covariate X?", reach for AME. The
#' \code{R = 0} additive-only case is the social relations model
#' (Warner, Kenny, Stoto 1979); \code{R >= 1} adds latent space.
#'
#' @section ERGM to AME translation:
#' For users coming from \code{statnet::ergm}, the rough analogues are:
#' \tabular{ll}{
#'   \strong{ERGM term}            \tab \strong{AME analogue} \cr
#'   \code{edges}                  \tab \code{intercept} (probit link, not logit) \cr
#'   \code{nodecov("x")}           \tab \code{Xrow = x} or \code{Xcol = x} \cr
#'   \code{nodematch("g")}         \tab dyadic covariate via \code{\link{nodematch}(g)} into \code{Xdyad} \cr
#'   \code{nodefactor("g")}        \tab dyadic covariate via \code{\link{nodefactor}(g)} into \code{Xdyad} (drop one level) \cr
#'   \code{absdiff("z")}           \tab dyadic covariate via \code{\link{absdiff}(z)} into \code{Xdyad} \cr
#'   \code{mutual}                 \tab \code{dcor = TRUE} -> the \code{rho} parameter (probit-scale, not log-odds; not numerically comparable to ERGM's \code{mutual}) \cr
#'   \code{gwesp} / transitivity   \tab \code{R >= 1} multiplicative latent factors (not the same statistic) \cr
#'   sender activity heterogeneity \tab \code{rvar = TRUE}, gives \code{a_i}, \code{va} \cr
#'   receiver popularity heterog.  \tab \code{cvar = TRUE}, gives \code{b_j}, \code{vb} \cr
#' }
#' For \code{family = "binary"} the link is probit, so the intercept is on the
#' probit scale; do not compare it to an ERGM \code{edges} estimate by simple
#' arithmetic.
#'
#' @section Migration from amen:
#' Both \code{amen} and \code{lame} export \code{ame()}; loading both
#' packages fires a startup warning telling you to call
#' \code{lame::ame(...)} or \code{amen::ame(...)} explicitly.
#'
#' For default cross-sectional calls, \code{lame::ame()} follows the
#' \code{amen::ame()} interface: the \code{fit$BETA} slot is a
#' 2-D \code{[n_stored, p]} matrix in both packages, so scripts that
#' call \code{colMeans(fit$BETA)} or \code{apply(fit$BETA, 2, mean)}
#' continue to work unchanged. \code{lame} additionally accepts
#' \code{family = "binary"} (which \code{amen} 1.4.5 no longer
#' accepts; \code{amen} requires \code{"bin"}). The \code{print}
#' argument is deprecated in favour of \code{verbose}; calls that pass
#' \code{print = ...} still work but warn.
#'
#' The cross-sectional path has no \code{dynamic_beta} option (an
#' AR(1) prior on a single-period coefficient is unidentified), so the
#' \code{BETA} 3-D shape that \code{lame()} can produce never arises
#' from \code{ame()}. See \code{\link{lame}} for the longitudinal path
#' and the silent-aggregation hazard with 2-D \code{apply(fit$BETA, 2, mean)}
#' scripts under \code{dynamic_beta = TRUE}.
#'
#' @section Notes on priors:
#' \itemize{
#'   \item The regression-coefficient prior is a Zellner g-prior
#'         (\code{beta ~ N(0, g * sigma^2 * (XtX)^-1)} where \code{XtX} is
#'         the design cross-product). \code{g} is a \emph{top-level}
#'         argument of \code{ame()}, not an entry in \code{prior = list(...)}
#'         (a common slip).
#'   \item There is no per-coefficient prior knob. If you need student_t,
#'         horseshoe, or to centre a slope away from zero, AME does not
#'         currently expose it; the g-prior structure is the only knob.
#'   \item Unknown names in \code{prior = list(...)} are warned about (a typo
#'         like \code{Sab = ...} instead of \code{Sab0} would otherwise be
#'         silently dropped).
#' }
#'
#' @section Chains and convergence policy:
#' \code{ame()} runs \code{chains} chains (default 2; \code{options(lame.chains = k)}
#' sets another default). Chain \eqn{k} is seeded with \code{seed + 1000 * (k - 1)},
#' so \code{chains = 1} is the first chain of a two-chain fit. On Unix the chains
#' run forked, on at most \code{cores} workers (all cores but one by default); on
#' Windows they run one after the other. The draws do not depend on the number of
#' workers. The result is the pooled fit of \code{\link{combine_ame_chains}}, with
#' \code{fit$chains}, \code{fit$chain_seeds}, \code{fit$rhat} (the R-hat of the
#' monitored set), the state each chain stopped in, and the \code{chain_indicator}
#' that \code{as_draws()} and \code{trace_plot()} use. Two chains cost twice the CPU
#' time and about the same wall time on a multicore machine.
#'
#' At the end of the post-burn run every chain is checked, and the chains are
#' compared with each other. Within a chain, the first and last thirds of the
#' stored draws of the intercept (or first coefficient), \code{va}, \code{vb},
#' \code{ve} and \code{rho} are compared with batch-means standard errors (ten
#' batches per third); a quantity is flagged when the two means differ by more than
#' 3 standard errors and by more than half a posterior standard deviation of the
#' whole series. Between chains, the split rank-normalised R-hat (Vehtari et al.,
#' 2021) of the monitored set is compared with 1.05: every regression coefficient,
#' \code{va}, \code{vb}, \code{cab}, \code{ve}, \code{rho}, and summaries of the
#' additive and latent structure recorded at every stored draw (the rms of
#' \code{a} and \code{b}, the largest actor effect, the rms of the latent product
#' \eqn{UV'}, the norm of \code{G} of a bipartite fit and the split of scale between
#' \code{U} and \code{V} of a directed fit). The dyadic products are not
#' monitored. Quantities a fit reports as not estimated are left out, and a run
#' with fewer than 60 stored draws per chain is not checked.
#' The gate is computed on the draws the returned fit holds (the last block of
#' every chain), and \code{summary()} and \code{\link{compute_mcmc_diagnostics}}
#' report the same statistic of the same draws: the Rhat column of
#' \code{summary()} is \code{fit$rhat} for every monitored coefficient and
#' variance component.
#'
#' When a chain drifts or an R-hat exceeds 1.05, every chain is extended by a block
#' of \code{nscan} sweeps from its own final state (the seeds advance by the sweeps
#' already run); the flagged block counts as extra burn-in, its draws are
#' discarded, and both checks run again. Up to \code{prior$chain_extensions}
#' extensions are made (default 3, so at most four times \code{nscan} post-burn
#' sweeps per chain). The returned fit has the same number of stored draws as a fit
#' that needed no extension, and \code{fit$auto_extend} records the extensions, the
#' sweeps discarded and the quantities that triggered each. If the last block is
#' still flagged, \code{ame()} stops with an error of class
#' \code{lame_chain_unsettled}: its \code{fits} field holds the chains, \code{fit}
#' the pooled fit of the last block, \code{rhat} the monitored table and
#' \code{moving} the offending quantities, and \code{lame_resume(e$fit, nscan =
#' nscan)} continues every chain (\code{lame_resume(e$fits[[k]], nscan = nscan)}
#' one). \code{prior = list(chain_extensions = 0)} refuses at the first flag;
#' \code{chains = 1} runs one chain with the drift check only. There is no mode
#' that returns a flagged fit. The measured false-extension rate and power are in
#' \code{\link{lame}}.
#'
#' @return
#' \strong{Posterior Samples (full MCMC chains):}
#' \item{BETA}{Regression coefficients (\code{nscan/odens} x p matrix; one
#' row per stored draw)}
#' \item{VC}{Variance components (\code{nscan/odens} x k matrix)}
#' \item{GOF}{Goodness-of-fit statistics ((\code{nscan/odens} + 1) x 5 matrix;
#' three columns, \code{sd.rowmean}, \code{sd.colmean} and \code{four.cycles},
#' for \code{mode = "bipartite"}).
#' First row contains observed values, remaining rows contain posterior predictive samples.
#' See \code{\link{gof}} for post-hoc computation and \code{\link{gof_plot}} for visualization.}
#'
#' \strong{Posterior Means (averaged over chain):}
#' \item{APM}{Additive row/sender effects (n-vector)}
#' \item{APM_sd, BPM_sd}{posterior standard deviation of each additive effect
#'   (static effects), used by \code{simulate()} to draw the effects}
#' \item{BPM}{Additive column/receiver effects (m-vector); NULL for symmetric networks}
#' \item{U}{Multiplicative row/sender factors (n xR matrix)}
#' \item{V}{Multiplicative column/receiver factors (m xR matrix); NULL for symmetric networks}
#' \item{L}{Eigenvalue matrix (R xR diagonal); symmetric networks only}
#' \item{YPM}{Stored response-scale fitted value: the average of the outcomes
#' replicated at each stored draw, except for Poisson fits, where it is the mean
#' count \eqn{\exp(\eta + \sigma^2 / 2)} at the posterior means of the linear
#' predictor and of the residual variance; see the Poisson model details.}
#'
#' \strong{Metadata:}
#' \item{family}{Model family (normal, binary, etc.)}
#' \item{mode}{Network mode (unipartite or bipartite)}
#' \item{sampler_health}{For MCMC fits, counts of numerical failures,
#' separated from valid prior-bound rejections. A numerical failure (a caught
#' collapsed reciprocity or utility-refresh error in unipartite fits, a failed
#' U, V or G draw in bipartite ones) aborts the fit with an error of class
#' \code{lame_sampler_failure} that carries the counts in its \code{health}
#' field, so a returned fit has none. A proposal outside the prior's box
#' (\code{prior$uv_max_abs}, \code{prior$G_max}) is a Metropolis rejection of
#' the model: that conditional keeps its state, and the rejections are counted
#' here. U and V proposal counts in bipartite fits are per proposed actor row
#' and G counts per Gaussian interaction-matrix draw.}
#' \item{sampler_health_iterations}{Number of burn-in plus production sweeps
#' included in the MCMC health counts. U and V proposal counts in bipartite
#' fits can exceed this number because each sweep visits multiple actor rows.}
#' \item{uv_g_scale_move}{For bipartite fits with positions, the number of sweeps
#' and of rescaled coordinates of the draw along the scale ridge of the positions
#' and the interaction matrix (see \code{uv_g_scale_move} in \code{prior}); absent
#' when the update did not run}
#' \item{ordinal_scale_move}{For \code{"ordinal"} fits, the sweeps the scale moves
#' visited (\code{visits}) and the sweeps at which the scale of the whole state
#' moved (\code{orbit_moved}); absent when the moves did not run}
#' \item{poisson_period_move}{For \code{"poisson"} fits, the sweeps the scale move
#' visited (\code{visits}) and its proposals and accepts (\code{scale_attempts},
#' \code{scale_accepts}, and for the correlation alone \code{scale_rho_attempts},
#' \code{scale_rho_accepts}); absent when the move did not run}
#' \item{binary_scale_move}{For \code{"binary"} fits, the draws of the scale of the linear
#' predictor (\code{visits}), the evaluations of its density they used (\code{evals}) and the
#' draws that moved the state (\code{moved}); absent when the draw did not run}
#' \item{actor_move}{For Poisson and ordinal fits, the work of the translations of
#' the levels with the latent array: \code{visits}, the level translations drawn
#' (\code{level_calls}), the likelihood evaluations of the Poisson ones
#' (\code{level_evals}), \code{stuck} and \code{skipped} steps, the draws of the
#' location of the latent array (\code{location_calls}, ordinal), the
#' coefficient translations (\code{coef_calls}, \code{coef_evals},
#' \code{coef_stuck}; Poisson), the position steps (\code{uv_calls},
#' \code{uv_evals}, \code{uv_stuck}; Poisson with positions) and the attempts and accepts of the scale move of
#' the dyadic variance and correlation (\code{scale_attempts},
#' \code{scale_accepts}, \code{rho_attempts}, \code{rho_accepts}; Poisson);
#' for binary fits, the work of the collapsed
#' per-actor update: \code{visits}, the slice steps of the level translations
#' (\code{level_calls}) and of the static positions (\code{uv_calls}), the
#' likelihood evaluations they used (\code{level_evals}, \code{uv_evals}), the
#' steps that failed to terminate (\code{stuck}) and, for bipartite fits, the
#' draws along the ridge of the intercept and the additive-effect means
#' (\code{ridge_calls}); for \code{"cbin"} and \code{"frn"} fits, the work of
#' the receiver-level translation: \code{visits}, the
#' receivers moved (\code{col_calls}) and the receivers left to the sampler
#' because the current state lay outside the interval the translation needs
#' (\code{col_skips}), the sweeps in which the level update moved a set of rows
#' (\code{level_calls}) and the rows it moved (\code{level_rows}), the proposals
#' and accepts of the reciprocity update (\code{rho_attempts},
#' \code{rho_accepts}) and, for a bipartite fit, \code{ridge_calls}; for the
#' other families of a bipartite fit with an intercept, \code{visits} and the
#' draws along the ridge of the intercept and the level means
#' (\code{ridge_calls}); absent when the updates did not run (see
#' \code{actor_move} in \code{prior})}
#' \item{symmetric}{Logical indicating if network is symmetric}
#' \item{R}{Dimension of multiplicative effects}
#' \item{prior_supplied}{\code{TRUE} when \code{prior} carried \code{Sab0},
#'   \code{eta0} or \code{etaab}: where the data alone leave the effect
#'   variances to the prior (a saturated nomination cap, see \code{odmax}),
#'   the reporting methods then show them under the supplied prior instead of
#'   as not estimated}
#'
#' \strong{Optional Posterior Samples (if requested via posterior_options):}
#' \item{U_samples}{Samples of U (n x R x iterations array)}
#' \item{V_samples}{Samples of V (m x R x iterations array)}
#' \item{a_samples}{Samples of row effects (n x iterations matrix)}
#' \item{b_samples}{Samples of column effects (m x iterations matrix); in a
#'   symmetric fit the same draws as \code{a_samples}}
#'
#' \strong{Note on the latent-scale matrices:}
#' The posterior-mean multiplicative product is stored on the fit
#' (\code{UVPM}, or \code{ULUPM} for symmetric fits); EZ (the expected
#' latent network) is not stored, to save memory. Accessors:
#' \itemize{
#'   \item \code{reconstruct_EZ(fit)} - Returns linear predictor (link scale, not response scale)
#'   \item \code{reconstruct_UVPM(fit)} - Returns the stored posterior-mean
#'     multiplicative product (\code{UVPM} / \code{ULUPM}) when present,
#'     otherwise U\%*\%t(V) or U\%*\%L\%*\%t(U)
#' }
#'
#' \strong{Generating posterior distributions:}
#' Use \code{simulate_posterior(fit, component="UV")} to generate posterior samples
#' for components where only means are stored, or use \code{posterior_options()}
#' during model fitting to save full posterior samples.
#' \item{model.name}{Name of the model (if provided)}
#' @seealso \code{\link{lame}} for longitudinal models,
#'   \code{\link{gof}} for post-hoc goodness-of-fit computation,
#'   \code{\link{gof_plot}} for visualizing GOF results,
#'   \code{\link{latent_positions}} for extracting latent positions as a tidy data frame,
#'   \code{\link{procrustes_align}} for Procrustes alignment of latent positions,
#'   \code{\link{summary.ame}} for model summaries,
#'   \code{\link{coef.ame}} for coefficient extraction
#' @author Cassy Dorff, Shahryar Minhas, Tosin Salau
#' @examples
#' \donttest{
#' data(YX_nrm)
#' # the default runs two chains and compares them on the regression
#' # coefficients and variance components
#' fit <- ame(YX_nrm$Y, Xdyad = YX_nrm$X, burn = 100, nscan = 600, odens = 10,
#'            verbose = FALSE)
#' summary(fit)
#' fit$rhat  # between-chain split R-hat; every value is below 1.05
#' # the defaults (burn = 500, nscan = 10000, odens = 25) are the length of a final run
#' }
#'
#' @export ame
ame <- function(
	Y,
	Xdyad = NULL,
	Xrow = NULL,
	Xcol = NULL,
	rvar = TRUE,
	cvar = TRUE,
	dcor = !symmetric,
	nvar = TRUE,
	R = 0,
	R_row = NULL,
	R_col = NULL,
	mode = c("unipartite", "bipartite"),
	family = "normal",
	intercept = !(family == "ordinal"),
	symmetric = FALSE,
	odmax = rep(max(apply(Y > 0, 1, sum, na.rm = TRUE)), nrow(Y)),
	prior = list(),
	g = NA,
	seed = 6886,
	nscan = 10000,
	burn = 500,
	odens = 25,
	verbose = TRUE,
	gof = TRUE,
	custom_gof = NULL,
	plot = FALSE,
	start_vals = NULL,
	periodic_save = FALSE,
	out_file = NULL,
	save_interval = 0.25,
	posterior_opts = NULL,
	n_chains = 1,
	cores = NULL,
	chains = .lame_chains_default(),
	use_sparse_matrices = FALSE,
	method = c("mcmc", "als"),
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	multistart = c("none", "cheap", "full"),
	save_log_lik = FALSE,
	ordinal_cutpoints = c("data_induced", "explicit"),
	print,
	...,
	# after `...` so a typo such as `model = "bin"` cannot partial-match into it
	log_lik_method = c("observed_exact", "observed_approx"),
	model.name = NULL,
	bias_correct = NULL
) {
	family_missing <- missing(family)
	mode_missing <- missing(mode)
	symmetric_missing <- missing(symmetric)
	netify_input <- .lame_prepare_netify_input(
		Y = Y,
		Xdyad = Xdyad,
		Xrow = Xrow,
		Xcol = Xcol,
		mode = mode,
		family = family,
		mode_missing = mode_missing,
		family_missing = family_missing,
		symmetric = symmetric,
		symmetric_missing = symmetric_missing,
		longitudinal = FALSE,
		caller = "ame",
		fit_method = if (identical(method[1L], "als")) "als" else "gibbs"
	)
	if (isTRUE(netify_input$converted)) {
		Y <- netify_input$Y
		Xdyad <- netify_input$Xdyad
		Xrow <- netify_input$Xrow
		Xcol <- netify_input$Xcol
		mode <- netify_input$mode
		family <- netify_input$family
		symmetric <- netify_input$symmetric
	}

	# resolve ordinal_cutpoints early so the non-ordinal warning fires once
	log_lik_method <- match.arg(log_lik_method)
	ordinal_cutpoints <- match.arg(ordinal_cutpoints)
	if (identical(ordinal_cutpoints, "explicit") && family != "ordinal") {
		cli::cli_warn(c(
			"{.arg ordinal_cutpoints} = {.val explicit} only takes effect when {.code family = \"ordinal\"}.",
			"i" = "Argument ignored for {.code family = \"{family}\"}."
		))
		ordinal_cutpoints <- "data_induced"
	}
	# accepted for signature parity with lame(); ame() has no live plotting
	invisible(plot)

	# lame()-only arguments have no meaning on a single cross-section
	dot_names <- names(list(...))
	lame_only <- c(
		"dynamic_beta", "dynamic_beta_kind", "dynamic_beta_pool",
		"dynamic_beta_per_actor", "dynamic_ab", "dynamic_uv",
		"dynamic_G", "time_index", "period_exposure",
		"max_seconds", "checkpoint_path", "checkpoint_every",
		"log_lik_path", "log_lik_chunk_size"
	)
	bad <- intersect(dot_names, lame_only)
	if (length(bad) > 0L) {
		cli::cli_abort(c(
			"{.arg {bad}} {?is/are} longitudinal-only argument{?s} and not valid for {.fn ame} (cross-sectional).",
			"i" = "Use {.fn lame} with a list of {.val T >= 2} matrices to enable dynamic / time-varying effects."
		))
	}
	# remaining unrecognised arguments (usually typos) get one warning
	if (length(dot_names) > 0L) {
		unknown <- setdiff(dot_names, lame_only)
		if (length(unknown) > 0L) {
			cli::cli_warn(c(
				"Unrecognised argument{?s} passed to {.fn ame}: {.arg {unknown}}.",
				"i" = "Check {.fn ?ame} for the supported signature; typos here are silently ignored."
			))
		}
	}

	# deprecated print argument; warn only when verbose was not also supplied
	if (!missing(print)) {
		if (missing(verbose)) {
			cli::cli_warn(c(
				"The {.arg print} argument is deprecated.",
				"i" = "Use {.arg verbose} instead."
			))
			verbose <- print
		}
	}

	# ALS fits go to the point estimator
	method <- match.arg(method)
	if (identical(method, "als")) {
		mode_arg <- match.arg(mode)
		# ordinal is not an ALS family; the ALS driver reports that, so skip the intercept guard
		if (!isTRUE(intercept) &&
				!identical(family, "ordinal") && !identical(family, "ord")) {
			cli::cli_abort(c(
				"{.fn ame} {.arg method} = {.val als} does not currently support {.code intercept = FALSE}.",
				"i" = "Use {.code method = \"mcmc\"} for a no-intercept fit, or keep the ALS intercept and center covariates as needed."
			))
		}
		# MCMC-only arguments the user supplied
		user_args <- setdiff(names(match.call()), "")
		forwarded <- c(
			"Y", "Xdyad", "Xrow", "Xcol", "R", "R_row", "R_col",
			"family", "mode", "symmetric",
			"bootstrap", "bootstrap_type",
			"bootstrap_block_length", "bootstrap_seed",
			"multistart", "bias_correct", "dcor",
			"verbose", "seed", "method", "intercept", "odmax", "prior",
			# cosmetic MCMC arguments are absorbed silently
			"model.name"
		)
		dropped <- setdiff(user_args, forwarded)
		if (length(dropped) > 0L) {
			cli::cli_warn(c(
				"{.fn ame} {.arg method} = {.val als}: ignored {length(dropped)} arguments that apply only to the MCMC path: {.arg {dropped}}.",
				"i" = "These do not affect the ALS fit. Use {.code method = \"mcmc\"} (the default) if you need them."
			))
		}
		if (length(family) == 1L && family %in% c("poisson", "pois")) .lame_check_poisson_counts(Y)
		# a bipartite fit uses one rank on the static ALS route
		R_als <- R
		if (!is.null(R_row) || !is.null(R_col)) {
			rr <- R_row %||% R; rc <- R_col %||% R
			if (rr != rc) {
				cli::cli_warn(c(
					"{.arg R_row} ({rr}) != {.arg R_col} ({rc}), but the ALS estimator uses a single latent dimension.",
					"i" = "Using {.code R = {max(rr, rc)}} for the ALS fit (the larger of the two)."
				))
			}
			R_als <- max(rr, rc, R)
		}
		return(ame_als(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			R = R_als,
			family = family,
			mode = mode_arg,
			symmetric = symmetric,
			bootstrap = bootstrap,
			bootstrap_type = bootstrap_type,
			bootstrap_block_length = bootstrap_block_length,
			bootstrap_seed = bootstrap_seed,
			multistart = match.arg(multistart),
			bias_correct = bias_correct,
			verbose = verbose,
			seed = seed,
			dcor = dcor,
			prior = prior
		))
	}

	# match call and mode
	mc <- match.call()
	mode <- match.arg(mode)

	# a chain of a multi-chain fit is a call with chains = 1 and a marker on the
	# prior; a call with chains > 1 runs the chains and pools them
	.blk <- attr(prior, "lame_chain_block", exact = TRUE)
	if (!is.null(.blk)) {
		attr(prior, "lame_chain_block") <- NULL
		mc <- .blk$mc
	} else {
		if (!missing(n_chains) && n_chains > 1) {
			if (!missing(chains) && !identical(as.integer(chains), as.integer(n_chains))) {
				cli::cli_abort(
					"{.arg chains} and {.arg n_chains} disagree; {.arg n_chains} is the older name of {.arg chains}.", class = "lame_chains_arg"
				)
			}
			chains <- n_chains
		}
		chains <- .lame_chains_check(chains)
		if (chains > 1L) return(.lame_chains_from_ame(environment(), mc, chains))
	}

	# warn on unknown posterior_opts keys and accept lame()'s save_UV_draws as an alias for save_UV
	if (!is.null(posterior_opts)) {
		if (!is.list(posterior_opts)) {
			cli::cli_abort(c(
				"{.arg posterior_opts} must be a list.",
				"i" = "Build it with {.fn posterior_options}, e.g. {.code posterior_opts = posterior_options(save_UV = TRUE)}."
			))
		}
		known_posterior_opts <- c(
			"save_UV", "save_UV_draws", "save_ab",
			"thin_UV", "thin_ab"
		)
		bad_opts <- setdiff(names(posterior_opts), known_posterior_opts)
		if (length(bad_opts) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg posterior_opts} list element{?s}: {.val {bad_opts}}.",
				"i" = "Recognised names: {.val {known_posterior_opts}}."
			))
		}
		if (isTRUE(posterior_opts$save_UV_draws)) {
			posterior_opts$save_UV <- TRUE
		}
	}
	# gof(), posterior simulation and the pointwise log-likelihood need the effect
	# draws aligned with BETA; for one network they are small, so keep them by default
	if (is.null(posterior_opts)) {
		uv_rank <- max(c(R, R_row, R_col, 0), na.rm = TRUE)
		posterior_opts <- list(save_UV = uv_rank > 0, save_ab = TRUE, thin_UV = 1L, thin_ab = 1L)
	}

	# the log-likelihood is computed after sampling from the saved draws, so it
	# needs per-iteration U, V, a and b storage
	want_log_lik <- isTRUE(save_log_lik) ||
		identical(save_log_lik, "chunked")
	if (want_log_lik) {
		if (identical(save_log_lik, "chunked")) {
			cli::cli_warn(c(
				"{.arg save_log_lik} = {.val chunked} is not yet supported on the cross-sectional {.fn ame} path.",
				"i" = "Falling back to in-memory storage."
			))
		}
		if (is.null(posterior_opts)) posterior_opts <- list()
		posterior_opts$save_UV <- TRUE
		posterior_opts$save_ab <- TRUE
		# the log-likelihood pairs BETA draw s with effect draw s, so effect storage must not be thinned
		if (!is.null(posterior_opts$thin_UV) && posterior_opts$thin_UV != 1L ||
				!is.null(posterior_opts$thin_ab) && posterior_opts$thin_ab != 1L) {
			cli::cli_warn(c(
				"{.arg save_log_lik} requires unthinned effect draws; overriding {.code thin_UV}/{.code thin_ab} to 1.",
				"i" = "The pointwise log-likelihood pairs each stored BETA draw with its own U/V and a/b draw."
			))
		}
		posterior_opts$thin_UV <- 1L
		posterior_opts$thin_ab <- 1L
	}
	# dynamic_beta is available in lame() only

	# seed locally and restore the global RNG stream on exit
	.had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
	.old_seed <- if (.had_seed) get(".Random.seed", envir = globalenv()) else NULL
	on.exit({
		if (.had_seed) {
			assign(".Random.seed", .old_seed, envir = globalenv())
		} else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
			rm(".Random.seed", envir = globalenv())
		}
	}, add = TRUE)

	# valid families, with the short amen spellings as aliases
	valid_families <- c(
		"normal", "binary", "ordinal",
		"cbin", "frn", "poisson"
	)
	amen_aliases <- c(
		nrm = "normal", bin = "binary", ord = "ordinal",
		pois = "poisson"
	)
	if (length(family) == 1L && is.character(family) &&
			family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		# announce the alias once per session
		if (!isTRUE(.lame_state$family_alias_msg_shown)) {
			cli::cli_inform(c(
				"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
				"i" = "Update your script to {.code family = \"{new_family}\"} when convenient.",
				"i" = "This note is shown once per session."
			))
			.lame_state$family_alias_msg_shown <- TRUE
		}
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
			!family %in% valid_families) {
		cli::cli_abort(c(
			"{.arg family} = {.val {family}} is not a recognised family.",
			"i" = "Choose one of: {.val {valid_families}}."
		))
	}
	.lame_explicit_ordinal_support(
		family,
		ordinal_cutpoints,
		identical(mode, "bipartite"),
		symmetric,
		dcor
	)
	.lame_validate_score_request(family, ordinal_cutpoints, log_lik_method, want_log_lik)
	# cbin and frn model each actor's ranked nominations row by row, so they are directed
	if (isTRUE(symmetric) && !identical(mode, "bipartite") &&
			family %in% c("cbin", "frn")) {
		cli::cli_abort(c(
			"{.arg symmetric} = TRUE is not supported for {.code family = \"{family}\"}.",
			"i" = "The {.val cbin} and {.val frn} families model each actor's own ranked nominations, so they are directed by construction.",
			"i" = "Fit with {.code symmetric = FALSE}, or recode the ties and use {.val binary} / {.val ordinal} for a symmetric fit."
		))
	}

	# bipartite ordinal, cbin, frn and Poisson use the rectangular samplers in R/rZ_bipartite.R

	# validate Y and the key MCMC arguments up front
	if (!is.matrix(Y) && !is.array(Y)) {
		cli::cli_abort(c(
			"{.arg Y} must be a matrix.",
			"i" = "For a numeric adjacency {.cls data.frame} (e.g. read in with {.fn read.csv}), convert it with {.code as.matrix(Y)}.",
			"i" = "For a long-format edgelist, build the network with {.fn netify::netify} and pass that object directly.",
			"i" = "{.fn as_lame_y} converts {.pkg igraph}/{.pkg network} objects to a lame-ready matrix."
		))
	}
	if (!is.numeric(Y) && !is.logical(Y)) {
		cli::cli_abort(c(
			"{.arg Y} must be numeric.",
			"i" = "Character input usually means an edgelist rather than an adjacency matrix; build the network with {.fn netify::netify} instead."
		))
	}
	if (any(is.infinite(Y) | is.nan(Y))) {
		cli::cli_abort(c(
			"{.arg Y} contains infinite or NaN values.",
			"i" = "Only finite values and {.val NA} (missing) are allowed."
		))
	}
	# symmetric = TRUE requires a symmetric Y
	if (isTRUE(symmetric) && nrow(Y) == ncol(Y)) {
		off <- Y - t(Y)
		if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
				"i" = "Symmetrize {.arg Y} (e.g. {.code (Y + t(Y))/2}) or set {.arg symmetric} = FALSE."
			))
		}
		# one triangle entirely NA passes the symmetry check vacuously, so refuse it
		ut <- upper.tri(Y); lt <- lower.tri(Y)
		ut_obs <- sum(is.finite(Y[ut]))
		lt_obs <- sum(is.finite(Y[lt]))
		if ((ut_obs == 0L) != (lt_obs == 0L)) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but only one triangle of {.arg Y} is observed (the other is entirely {.val NA}).",
				"i" = "The symmetric sampler would silently mirror the populated triangle.",
				"i" = "Symmetrize {.arg Y} explicitly with {.code Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]} (or vice versa) before fitting."
			))
		}
		# symmetric = TRUE with an asymmetric Xdyad slice gives a non-symmetric predictor,
		# so it is refused; a hole in one triangle is filled from the other in ame_unipartite()
		if (!is.null(Xdyad)) {
			Xd_check <- if (is.matrix(Xdyad)) array(Xdyad, c(dim(Xdyad), 1L)) else Xdyad
			if (length(dim(Xd_check)) == 3L) .lame_check_symmetric_xdyad(Xd_check)
		}
	}
	obs_y <- Y[is.finite(Y)]
	if (length(obs_y) == 0L) {
		cli::cli_abort("{.arg Y} has no observed (non-missing, finite) values.")
	}
	if (nrow(Y) < 3L || ncol(Y) < 3L) {
		cli::cli_abort(c(
			"{.arg Y} is too small ({nrow(Y)}x{ncol(Y)}) to fit an AME model.",
			"i" = "At least 3 actors per mode are required."
		))
	}
	# a constant response cannot be fit under any family
	eff_y <- if (family %in% c("binary", "cbin")) 1 * (obs_y > 0) else obs_y
	if (length(unique(eff_y)) == 1L) {
		hint <- if (family %in% c("binary", "cbin")) {
			c(
				"i" = "All observed cells threshold to {.val {unique(eff_y)}}.",
				"i" = "If your edgelist lists only the ties that occurred, absent dyads should enter as 0s (in {.pkg netify}, revisit {.arg missing_to_zero}), not be missing."
			)
		} else {
			c("i" = "Every observed cell equals {.val {unique(eff_y)}}.")
		}
		cli::cli_abort(c(
			"{.arg Y} is constant; there is no relational variation to model.",
			hint
		))
	}
	if (length(R) != 1L || anyNA(R) || R < 0 ||
			!isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	if (R >= min(nrow(Y), ncol(Y))) {
		cli::cli_abort("{.arg R} = {R} must be smaller than the network dimension ({min(nrow(Y), ncol(Y))}).")
	}
	# bipartite latent dimensions: U is nA x R_row and V is nB x R_col
	for (rr in c("R_row", "R_col")) {
		rv <- get(rr)
		if (!is.null(rv)) {
			if (length(rv) != 1L || anyNA(rv) || rv < 0 ||
					!isTRUE(all.equal(rv, round(rv)))) {
				cli::cli_abort("{.arg {rr}} must be a single non-negative integer.")
			}
			rmax_bip <- if (identical(rr, "R_row")) nrow(Y) else ncol(Y)
			if (rv >= rmax_bip) {
				cli::cli_abort("{.arg {rr}} = {rv} must be smaller than the network dimension ({rmax_bip}).")
			}
		}
	}
	if (length(nscan) != 1L || anyNA(nscan) || nscan < 1L) {
		cli::cli_abort("{.arg nscan} must be a positive integer (post-burn-in iterations).")
	}
	if (length(burn) != 1L || anyNA(burn) || burn < 0L) {
		cli::cli_abort("{.arg burn} must be a non-negative integer.")
	}
	if (length(odens) != 1L || anyNA(odens) || odens < 1L || odens > nscan) {
		cli::cli_abort(c(
			"{.arg odens} must be a positive integer no greater than {.arg nscan} = {nscan}.",
			"i" = "Otherwise no posterior draws would be stored."
		))
	}
	if (family == "poisson") .lame_check_poisson_counts(obs_y)
	# a nomination cap below a row's observed out-degree contradicts Y
	if (family %in% c("cbin", "frn") && !is.null(odmax) && !anyNA(odmax) &&
			length(odmax) %in% c(1L, nrow(Y))) {
		Y_od <- Y
		if (!identical(mode, "bipartite")) diag(Y_od) <- NA
		.lame_check_odmax_outdegree(rowSums(Y_od > 0, na.rm = TRUE), odmax, rownames(Y))
	}
	# binary and cbin threshold observed values at zero
	if (family %in% c("binary", "cbin")) {
		uy <- unique(obs_y)
		if (any(!uy %in% c(0, 1))) {
			cli::cli_warn(c(
				"{.arg family} = {.val {family}} but {.arg Y} contains values other than 0/1.",
				"i" = "{.arg Y} will be thresholded to {.code 1 * (Y > 0)}; if you meant counts, use {.val poisson}, or {.val ordinal}/{.val normal} as appropriate."
			))
		}
	}
	# binary-looking Y under family = "normal" warns but does not block
	if (family == "normal") {
		uy <- unique(obs_y)
		if (length(uy) <= 2L && all(uy %in% c(0, 1))) {
			cli::cli_inform(c(
				"i" = "{.arg Y} looks binary (values in {.val {{0, 1}}}) but {.arg family} = {.val normal}.",
				"i" = "If this is a 0/1 network, pass {.code family = \"binary\"}."
			))
		}
	}
	# below 2 percent density few mutual dyads are expected; the note says what that means for rho
	if (family %in% c("binary", "cbin") && isTRUE(dcor) &&
			mode == "unipartite" && !symmetric) {
		dens <- mean(obs_y > 0)
		if (is.finite(dens) && dens < 0.02) {
			cli::cli_warn(c(
				"Network density is {round(100 * dens, 2)}%; the dyadic correlation {.code rho} is informed by the joint pattern of ties within each dyad, and with few mutual dyads mainly by the absence of mutual ties.",
				"i" = .lame_rho_few_mutual_text
			))
		}
	}
	# check Xdyad's actor dimensions before the separation check indexes its slices against Y
	if (!is.null(Xdyad)) {
		dimd_early <- dim(Xdyad)
		if (!is.null(dimd_early) && length(dimd_early) == 3L &&
				(dimd_early[1] != nrow(Y) || dimd_early[2] != ncol(Y))) {
			cli::cli_abort(c(
				"{.arg Xdyad} dimensions ({dimd_early[1]} x {dimd_early[2]} x {dimd_early[3]}) do not match {.arg Y} ({nrow(Y)} x {ncol(Y)}).",
				"i" = "The dyadic covariate array's first two dimensions must be the actor dimensions of {.arg Y}."
			))
		}
	}
	# complete separation by a single dyadic covariate: the data do not bound its
	# coefficient, so the posterior is set by the prior
	if (family %in% c("binary", "cbin") && !is.null(Xdyad) &&
			length(dim(Xdyad)) == 3L) {
		yb <- 1 * (Y > 0)
		sep_k <- integer(0)
		sep_above <- logical(0)
		for (k in seq_len(dim(Xdyad)[3L])) {
			xk <- Xdyad[, , k]
			ok <- !is.na(yb) & !is.na(xk)
			x1 <- xk[ok & yb == 1]
			x0 <- xk[ok & yb == 0]
			if (length(x1) > 0L && length(x0) > 0L &&
					(min(x1) > max(x0) || max(x1) < min(x0))) {
				sep_k <- c(sep_k, k)
				sep_above <- c(sep_above, min(x1) > max(x0))
			}
		}
		if (length(sep_k) > 0L) {
			sep_sd <- .ame_separation_prior_sd(Y, Xdyad, Xrow, Xcol, intercept, mode, g, prior, sep_k)
			for (j in seq_along(sep_k)) {
				nm_k <- dimnames(Xdyad)[[3L]][sep_k[j]] %||% as.character(sep_k[j])
				dir_word <- if (sep_above[j]) "above" else "below"
				sd_txt <- if (is.finite(sep_sd[j])) paste0(
					", a normal with sd ",
					signif(sep_sd[j], 2L),
					","
				) else ""
				cli::cli_warn(c(
					"{.arg Xdyad} covariate {.val {nm_k}} separates the ties from the non-ties: every tie lies {dir_word} every non-tie on it.",
					"i" = "The likelihood increases with its coefficient without a bound, so the data do not determine it: its posterior is its prior{sd_txt} carried out to where the dyads nearest the boundary are classified. Read its sign and the fitted probabilities, not its size.",
					"i" = if (identical(family, "binary")) "The sampler draws the scale of the coefficients and additive effects exactly ({.code prior$binary_scale_move}), so the posterior does not depend on the number of sweeps."
					else "The censored-binary sampler has no exact draw of the scale and crosses that direction slowly (3 to 10 effective draws of 2,000 stored on a 30-actor panel)."
				))
			}
		}
	}
	# an exactly symmetric Y under the directed model pushes rho towards one; point at symmetric = TRUE
	if (!symmetric && mode == "unipartite" && nrow(Y) == ncol(Y) &&
			isTRUE(all.equal(unname(Y), unname(t(Y))))) {
		.lame_inform_symmetric_y(dcor, family = family)
	}
	if (burn > nscan) {
		cli::cli_warn(c(
			"{.arg burn} = {burn} is greater than {.arg nscan} = {nscan}.",
			"i" = "Most iterations will be discarded; consider {.arg burn} <= {.arg nscan}."
		))
	}
	# fully unobserved actors (all-NA rows or columns) are fit from the prior alone
	if (mode == "unipartite") {
		drow <- apply(Y, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y, 2, function(c) all(is.na(c) | !is.finite(c)))
		# the diagonal does not count towards an all-NA row
		drow_ex <- vapply(
			seq_len(nrow(Y)),
			function(i)
				all(is.na(Y[i, -i]) | !is.finite(Y[i, -i])),
			logical(1)
		)
		dcol_ex <- vapply(
			seq_len(ncol(Y)),
			function(j)
				all(is.na(Y[-j, j]) | !is.finite(Y[-j, j])),
			logical(1)
		)
		if (any(drow_ex) || any(dcol_ex)) {
			cli::cli_inform(c(
				"i" = "{.arg Y} has {sum(drow_ex)} fully-unobserved row{?s} and {sum(dcol_ex)} fully-unobserved column{?s}.",
				"i" = "Their additive effects (APM/BPM) will be estimated entirely from the prior."
			))
		}
	} else {
		drow <- apply(Y, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y, 2, function(c) all(is.na(c) | !is.finite(c)))
		# the bipartite binary Z sampler needs observed data in every row and column
		if ((any(drow) || any(dcol)) && family %in% c("binary")) {
			cli::cli_abort(c(
				"{.arg Y} (bipartite, {.val {family}}) has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
				"x" = "The bipartite {.val {family}} sampler cannot proceed with rows/columns that have no observed data.",
				"i" = "Drop the empty actor(s) from {.arg Y} (and matching covariates) before fitting, or use {.val normal} family."
			))
		}
		if (any(drow) || any(dcol)) {
			cli::cli_inform(c(
				"i" = "{.arg Y} has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
				"i" = "Their additive effects (APM/BPM) will be estimated entirely from the prior."
			))
		}
	}
	# covariate finiteness and shape checks
	.check_cov_finite <- function(x, nm) {
		if (is.null(x)) return(invisible())
		v <- if (is.data.frame(x)) as.matrix(x) else x
		if (!is.numeric(v) && !is.logical(v)) {
			cli::cli_abort("{.arg {nm}} must be numeric.")
		}
		if (any(is.infinite(v) | is.nan(v))) {
			cli::cli_abort(c(
				"{.arg {nm}} contains infinite or NaN values.",
				"i" = "Only finite values and {.val NA} are allowed."
			))
		}
		invisible()
	}
	.check_cov_finite(Xrow, "Xrow")
	.check_cov_finite(Xcol, "Xcol")
	.check_cov_finite(Xdyad, "Xdyad")
	# covariates are aligned to Y by actor name below, after Xdyad becomes a 3-D array

	# a 2-D Xdyad matrix is wrapped into a 3-D array with a short message
	if (!is.null(Xdyad)) {
		dimd <- dim(Xdyad)
		if (is.null(dimd) || length(dimd) == 2L) {
			if (is.null(dimd)) {
				cli::cli_abort(c(
					"{.arg Xdyad} must be a 2-D matrix or 3-D array, not a {.cls {class(Xdyad)[1]}}."
				))
			}
			if (dimd[1] != nrow(Y) || dimd[2] != ncol(Y)) {
				cli::cli_abort(c(
					"{.arg Xdyad} dimensions ({dimd[1]} x {dimd[2]}) do not match {.arg Y} ({nrow(Y)} x {ncol(Y)})."
				))
			}
			xname <- deparse(substitute(Xdyad))[1]
			# flag only a real variable name, not an inline expression
			if (nchar(xname) > 0L && nchar(xname) < 40L &&
					!grepl("[(){}]", xname)) {
				cli::cli_inform(c(
					"i" = "{.arg Xdyad} was a 2-D matrix; wrapping as {.code array(., dim = c(n, n, 1))} with slice name {.val {xname}}."
				))
			}
			Xdyad <- array(
				Xdyad,
				dim = c(dimd[1], dimd[2], 1L),
				dimnames = list(rownames(Xdyad), colnames(Xdyad), xname)
			)
			dimd <- dim(Xdyad)
		} else if (length(dimd) != 3L) {
			cli::cli_abort(c(
				"{.arg Xdyad} must be a 3-dimensional array ({nrow(Y)} x {ncol(Y)} x p).",
				"i" = "Got a {length(dimd)}-D object; wrap a single covariate as {.code array(x, dim = c(n, n, 1))}."
			))
		}
		if (dimd[1] != nrow(Y) || dimd[2] != ncol(Y)) {
			cli::cli_abort(c(
				"{.arg Xdyad} dimensions ({dimd[1]} x {dimd[2]} x {dimd[3]}) do not match {.arg Y} ({nrow(Y)} x {ncol(Y)})."
			))
		}
	}
	# nodal covariates need one row per actor
	if (!is.null(Xrow) && !is.null(dim(Xrow)) && nrow(Xrow) != nrow(Y)) {
		cli::cli_abort(c(
			"{.arg Xrow} has {nrow(Xrow)} rows but {.arg Y} has {nrow(Y)} row actors.",
			"i" = "{.arg Xrow} needs exactly one row per row actor of {.arg Y}."
		))
	}
	if (!is.null(Xcol) && !is.null(dim(Xcol)) && nrow(Xcol) != ncol(Y)) {
		cli::cli_abort(c(
			"{.arg Xcol} has {nrow(Xcol)} rows but {.arg Y} has {ncol(Y)} column actors.",
			"i" = "{.arg Xcol} needs exactly one row per column actor of {.arg Y}."
		))
	}
	# align covariates to Y by actor name when both are named: permuted inputs are
	# reordered, mismatched actors abort, unnamed inputs stay positional
	.align_idx <- function(have, want, what) {
		if (is.null(have) || is.null(want)) return(NULL)
		if (length(have) != length(want) || identical(have, want)) return(NULL)
		# default sequential rownames ("1", "2", ...) are not real names
		if (identical(have, as.character(seq_along(have)))) return(NULL)
		if (setequal(have, want)) return(match(want, have))
		cli::cli_abort(c(
			"{.arg {what}} carries actor names that do not match {.arg Y}'s.",
			"i" = "Covariates are aligned to {.arg Y} by name; rename {.arg {what}} to match, or drop its dimnames to force positional alignment."
		))
	}
	if (!is.null(Xdyad) && length(dim(Xdyad)) == 3L) {
		ir <- .align_idx(dimnames(Xdyad)[[1]], rownames(Y), "Xdyad")
		jc <- .align_idx(dimnames(Xdyad)[[2]], colnames(Y), "Xdyad")
		if (!is.null(ir) || !is.null(jc)) {
			Xdyad <- Xdyad[ir %||% seq_len(dim(Xdyad)[1L]),
				jc %||% seq_len(dim(Xdyad)[2L]), , drop = FALSE]
			cli::cli_inform(c(
				"i" = "Reordered {.arg Xdyad} to match {.arg Y}'s actor order (aligned by dimnames)."
			))
		}
	}
	if (!is.null(Xrow) && !is.null(dim(Xrow))) {
		ir <- .align_idx(rownames(Xrow), rownames(Y), "Xrow")
		if (!is.null(ir)) {
			Xrow <- Xrow[ir, , drop = FALSE]
			cli::cli_inform(c(
				"i" = "Reordered {.arg Xrow} to match {.arg Y}'s actor order (aligned by rownames)."
			))
		}
	}
	if (!is.null(Xcol) && !is.null(dim(Xcol))) {
		ir <- .align_idx(rownames(Xcol), colnames(Y), "Xcol")
		if (!is.null(ir)) {
			Xcol <- Xcol[ir, , drop = FALSE]
			cli::cli_inform(c(
				"i" = "Reordered {.arg Xcol} to match {.arg Y}'s column actor order (aligned by rownames)."
			))
		}
	}
	# an entirely-NA covariate cannot be identified
	.check_all_na_cov <- function(x, nm) {
		if (is.null(x)) return(invisible())
		if (length(dim(x)) == 3L) {
			bad <- which(apply(x, 3L, function(s) all(is.na(s))))
		} else if (!is.null(dim(x))) {
			bad <- which(apply(x, 2L, function(cc) all(is.na(cc))))
		} else {
			return(invisible())
		}
		if (length(bad) > 0L) {
			cli::cli_abort(c(
				"{.arg {nm}} covariate{?s} {.val {bad}} {?is/are} entirely NA.",
				"i" = "An all-NA covariate cannot be identified; drop it or supply values."
			))
		}
		invisible()
	}
	.check_all_na_cov(Xdyad, "Xdyad")
	.check_all_na_cov(Xrow, "Xrow")
	.check_all_na_cov(Xcol, "Xcol")
	# warn on unrecognised prior-list keys
	if (length(prior) > 0L && !is.null(names(prior))) {
		known_prior <- c(
			"Sab0", "eta0", "etaab", "Suv0", "kappa0",
			"s20", "s2u0",
			"rho_uv_mean", "rho_uv_sd",
			"rho_ab_mean", "rho_ab_sd",
			"uv_split_move",
			"cov_ridge_move",
			"sparse_position_move",
			"sparse_position_every",
			"actor_move", "actor_move_every", "uv_g_scale_move", "ordinal_scale_move",
			"binary_scale_move", "binary_scale_every", "chain_extensions"
		)
		bad <- setdiff(names(prior), known_prior)
		if (length(bad) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg prior} list element{?s}: {.val {bad}}.",
				"i" = "Recognised names: {.val {known_prior}}.",
				"i" = "Note: {.code g} is a top-level argument to {.fn ame}, not a {.arg prior} entry."
			))
		}
		non_num <- intersect(
			names(prior)[!vapply(prior, is.numeric, logical(1))],
			setdiff(known_prior, c("actor_move", "uv_g_scale_move", "ordinal_scale_move", "binary_scale_move"))
		)
		non_num <- setdiff(non_num, c("uv_split_move", "sparse_position_move", "cov_ridge_move"))
		if (length(non_num) > 0L) {
			cli::cli_abort(c(
				"{.arg prior} element{?s} {.val {non_num}} must be numeric.",
				"i" = "Got class{?es}: {.cls {vapply(prior[non_num], function(x) class(x)[1], character(1))}}."
			))
		}
		.lame_chain_check_prior(prior)
		if (!is.null(prior[["actor_move"]]) &&
				(!is.logical(prior[["actor_move"]]) || length(prior[["actor_move"]]) != 1L ||
						is.na(prior[["actor_move"]])))
			cli::cli_abort("prior$actor_move must be TRUE or FALSE.")
		if (!is.null(prior[["cov_ridge_move"]]) &&
				(!is.logical(prior[["cov_ridge_move"]]) || length(prior[["cov_ridge_move"]]) != 1L ||
						is.na(prior[["cov_ridge_move"]])))
			cli::cli_abort("prior$cov_ridge_move must be TRUE or FALSE.")
		if (!is.null(prior[["uv_split_move"]]) &&
				(!is.logical(prior[["uv_split_move"]]) || length(prior[["uv_split_move"]]) != 1L ||
						is.na(prior[["uv_split_move"]])))
			cli::cli_abort("prior$uv_split_move must be TRUE or FALSE.")
		if (!is.null(prior[["sparse_position_every"]]) &&
				(!is.numeric(prior[["sparse_position_every"]]) || length(prior[["sparse_position_every"]]) != 1L ||
						!is.finite(prior[["sparse_position_every"]]) || prior[["sparse_position_every"]] < 1))
			cli::cli_abort("prior$sparse_position_every must be a positive integer.")
		if (!is.null(prior[["sparse_position_move"]]) &&
				(!is.logical(prior[["sparse_position_move"]]) || length(prior[["sparse_position_move"]]) != 1L ||
						is.na(prior[["sparse_position_move"]])))
			cli::cli_abort("prior$sparse_position_move must be TRUE or FALSE.")
		if (!is.null(prior[["uv_g_scale_move"]]) &&
				(!is.logical(prior[["uv_g_scale_move"]]) || length(prior[["uv_g_scale_move"]]) != 1L ||
						is.na(prior[["uv_g_scale_move"]])))
			cli::cli_abort("prior$uv_g_scale_move must be TRUE or FALSE.")
		for (nm_ in c("ordinal_scale_move", "binary_scale_move")) {
			if (!is.null(prior[[nm_]]) && (!is.logical(prior[[nm_]]) || length(prior[[nm_]]) != 1L || is.na(prior[[nm_]])))
				cli::cli_abort(paste0("prior$", nm_, " must be TRUE or FALSE."))
		}
		if (!is.null(prior[["actor_move_every"]]) &&
				(!is.numeric(prior[["actor_move_every"]]) || length(prior[["actor_move_every"]]) != 1L ||
						!is.finite(prior[["actor_move_every"]]) || prior[["actor_move_every"]] < 1 ||
						prior[["actor_move_every"]] != floor(prior[["actor_move_every"]])))
			cli::cli_abort("prior$actor_move_every must be a positive whole number.")
		if (!is.null(prior[["binary_scale_every"]]) &&
				(!is.numeric(prior[["binary_scale_every"]]) || length(prior[["binary_scale_every"]]) != 1L ||
						!is.finite(prior[["binary_scale_every"]]) || prior[["binary_scale_every"]] < 1 ||
						prior[["binary_scale_every"]] != floor(prior[["binary_scale_every"]])))
			cli::cli_abort("prior$binary_scale_every must be a positive whole number.")
	}
	if (isTRUE(periodic_save) && is.null(out_file)) {
		cli::cli_abort(c(
			"{.arg periodic_save} = TRUE requires an {.arg out_file} path.",
			"i" = "Supply {.arg out_file} or set {.arg periodic_save} = FALSE."
		))
	}

	# share of rows at the nomination cap (cbin / frn), read by the intercept
	# reporting rule (.lame_cap_level_unreported)
	.odmax_sat_share <- .lame_odmax_saturated_share(Y, odmax, family)
	.cells_hidden <- family %in% c("cbin", "frn") && is.matrix(Y) &&
		any(.lame_nom_hidden_cells(Y, identical(mode[1L], "bipartite")))
	# a supplied effects prior (Sab0, eta0 or etaab) means the effect variances are reported under it
	.ab_prior_supplied <- .lame_ab_prior_given(prior)
	for (st_ in .lame_ab_statements(list(
		family = family,
		mode = mode[1L],
		odmax_saturated_share = .odmax_sat_share,
		cells_hidden = .cells_hidden,
		n_time = 1L,
		prior_supplied = .ab_prior_supplied
	))) {
		cli::cli_inform(c("i" = "{st_}"))
	}

	# unipartite dispatch
	if (mode == "unipartite") {
		if (nrow(Y) != ncol(Y)) {
			cli::cli_abort(c(
				"Matrix Y must be square for unipartite networks.",
				"i" = "Y has dimensions {nrow(Y)} x {ncol(Y)}.",
				"i" = "Use mode='bipartite' for rectangular matrices."
			))
		}
		# symmetric ordinal uses one unit-variance utility per unordered dyad, mirrored each sweep

		fit <- ame_unipartite(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			rvar = rvar,
			cvar = cvar,
			dcor = dcor,
			nvar = nvar,
			R = R,
			family = family,
			intercept = intercept,
			symmetric = symmetric,
			odmax = odmax,
			prior = prior,
			g = g,
			seed = seed,
			nscan = nscan,
			burn = burn,
			odens = odens,
			verbose = verbose,
			gof = gof,
			custom_gof = custom_gof,
			start_vals = start_vals,
			periodic_save = periodic_save,
			out_file = out_file,
			save_interval = save_interval,
			model.name = model.name,
			posterior_opts = posterior_opts,
			use_sparse_matrices = use_sparse_matrices,
			ordinal_cutpoints = ordinal_cutpoints
		)
		if (!isTRUE(.blk$nopolicy)) fit <- .ame_chain_policy(fit, environment(), ame_unipartite)
		fit$call <- mc

		# bipartite dispatch
	} else if (mode == "bipartite") {
		if (symmetric) {
			cli::cli_warn("Symmetric option ignored for bipartite networks.")
			symmetric <- FALSE
		}

		if (is.null(R_row)) R_row <- R
		if (is.null(R_col)) R_col <- R

		fit <- ame_bipartite(
			Y = Y,
			Xdyad = Xdyad,
			Xrow = Xrow,
			Xcol = Xcol,
			rvar = rvar,
			cvar = cvar,
			R_row = R_row,
			R_col = R_col,
			family = family,
			intercept = intercept,
			odmax = odmax,
			prior = prior,
			g = g,
			seed = seed,
			nscan = nscan,
			burn = burn,
			odens = odens,
			verbose = verbose,
			gof = gof,
			custom_gof = custom_gof,
			start_vals = start_vals,
			periodic_save = periodic_save,
			out_file = out_file,
			save_interval = save_interval,
			model.name = model.name,
			posterior_opts = posterior_opts,
			use_sparse_matrices = use_sparse_matrices
		)
		if (!isTRUE(.blk$nopolicy)) fit <- .ame_chain_policy(fit, environment(), ame_bipartite)
		fit$call <- mc
	}
	if (!isTRUE(.blk$record)) fit$chain_summary <- NULL
	# nomination cap share and effects-prior flag for the reporting methods
	if (is.finite(.odmax_sat_share)) {
		fit$odmax_saturated_share <- .odmax_sat_share
		fit$odmax_sat_threshold <- .lame_odmax_sat_threshold_single
	}
	fit$cells_hidden <- .cells_hidden
	fit$prior_supplied <- .ab_prior_supplied
	fit$coefficient_prior <- if (identical(mode, "bipartite")) {
		list(
			kind = "independent_fixed", description = "Regression coefficients have independent N(0, g * s20) priors, with fixed s20 = prior$s20 (default 1 on the latent scale of the binary, ordinal, cbin, frn and Poisson families, the variance of the response for the normal family)."
		)
	} else {
		list(
			kind = "single_design", description = "Regression coefficients have a zero-mean design-based Gaussian g-prior, independent of residual s2, with a proportional numerical ridge and additional intercept variance; see rbeta_ab_fc()."
		)
	}
	# cbin / frn with hidden cells: effective sample size of the intercept (.lame_nom_hidden_ess)
	if (family %in% c("cbin", "frn") && is.matrix(Y)) {
		fit$hidden_cell_ess <- .lame_nom_hidden_ess(
			fit,
			any(.lame_nom_hidden_cells(Y, identical(mode, "bipartite")))
		)
	}

	# coefficient-row indices of independently thinned saved blocks; the symmetric sampler uses nvar for both margins
	fit$rvar <- if (symmetric && mode == "unipartite") nvar else rvar
	fit$cvar <- if (symmetric && mode == "unipartite") nvar else cvar
	if (!is.null(fit$a_samples) || !is.null(fit$b_samples))
		fit$ab_draw_index <- seq_len(ncol(fit$a_samples %||% fit$b_samples)) * (posterior_opts$thin_ab %||% 1L)
	if (!is.null(fit$U_samples))
		fit$uv_draw_index <- seq_len(dim(fit$U_samples)[3L]) * (posterior_opts$thin_UV %||% 1L)

	# observed marginal scores, or explicitly approximate ones; LOO/WAIC eligibility is recorded separately
	if (want_log_lik) {
		# cells left unobserved for a missing covariate are not scored
		Y_scored <- Y
		if (is.matrix(fit$Y) && identical(dim(fit$Y), dim(Y))) Y_scored[is.na(fit$Y)] <- NA
		fit$log_lik <- .ame_compute_log_lik_post_hoc(fit, Y = Y_scored, family = family, method = log_lik_method)
		fit$log_lik_method <- log_lik_method
		fit$log_lik_index <- .lame_score_index(Y_scored, identical(mode, "bipartite"), symmetric)
		fit$log_lik_contract <- .lame_score_contract(
			family,
			fit$log_lik_method,
			ordinal_cutpoints,
			identical(mode, "bipartite"),
			symmetric,
			dcor,
			fit$log_lik_index,
			if ("rho" %in% colnames(fit$VC)) fit$VC[, "rho"] else NULL
		)
	}

	return(fit)
}

# per-iteration observed marginal scores from matched saved draws of every active block
#' @noRd
.ame_compute_log_lik_post_hoc <- function(fit, Y, family, method = "observed_exact") {
	.lame_validate_score_request(family, fit$ordinal_cutpoints, method, TRUE)
	if (is.null(fit$BETA) || is.null(fit$VC)) {
		cli::cli_warn(c("Cannot compute log_lik: BETA / VC missing from fit."))
		return(NULL)
	}
	n_iter <- nrow(fit$BETA)
	# a symmetric relation is scored once per unordered dyad
	bip <- identical(fit$mode, "bipartite") || (is.null(fit$mode) && nrow(Y) != ncol(Y))
	if (family %in% c("binary", "cbin")) Y <- 1 * (Y > 0)
	index <- .lame_score_index(Y, bip, isTRUE(fit$symmetric))
	obs_idx <- index[, 1L] + (index[, 2L] - 1L) * nrow(Y)
	n_obs <- length(obs_idx)
	if (n_obs == 0L) return(NULL)
	n_a <- nrow(Y); n_b <- ncol(Y)
	ve_col <- match("ve", colnames(fit$VC))
	if (is.na(ve_col)) cli::cli_abort("Post-hoc scoring requires named observation variance draws.")
	beta_names <- colnames(fit$BETA)
	has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"
	# active latent factors need matched saved U and V draws
	have_U <- !is.null(fit$U_samples) && !is.null(fit$V_samples)
	# Xdyad is stored on fit$X when ame() saves the call arguments; otherwise rebuilt from the call
	X <- fit[["X"]]
	ll <- matrix(NA_real_, nrow = n_iter, ncol = n_obs)
	# linear predictor for draw s
	compute_eta_draw <- function(s) {
		beta_s <- if (length(dim(fit$BETA)) == 3L) fit$BETA[s, , 1] else fit$BETA[s, ]
		eta_s <- matrix(0, n_a, n_b)
		if (is.null(X)) {
			if (has_intercept) eta_s <- eta_s + beta_s[1]
			if (length(beta_s) > as.integer(has_intercept)) cli::cli_abort("Post-hoc scoring requires the fitted design array.")
		} else {
			Xarr <- if (length(dim(X)) == 3L) X else array(X, c(n_a, n_b, 1L))
			p_x <- dim(Xarr)[3L]
			p_b <- length(beta_s)
			if (p_b == p_x) {
				# BETA columns match the design slices 1:1; slice 1 is the intercept when present
				for (k in seq_len(p_x)) {
					eta_s <- eta_s + beta_s[k] * Xarr[, , k]
				}
			} else {
				# symmetric fits store one coefficient per nodal covariate while the design keeps
				# row and column slices: [intercept | row | col | dyad] vs [intercept | node | dyad]
				ic <- as.integer(has_intercept)
				pr <- p_x - p_b
				pd <- p_b - ic - pr
				if (pr >= 0 && pd >= 0) {
					if (ic) eta_s <- eta_s + beta_s[1] * Xarr[, , 1]
					for (j in seq_len(pr)) {
						eta_s <- eta_s + beta_s[ic + j] *
							(Xarr[, , ic + j] + Xarr[, , ic + pr + j])
					}
					for (k in seq_len(pd)) {
						eta_s <- eta_s + beta_s[ic + pr + k] *
							Xarr[, , ic + 2L * pr + k]
					}
				} else cli::cli_abort("Saved coefficients cannot be mapped to the fitted design array.")
			}
		}
		# every active block needs the draw that belongs to BETA[s, ]
		match_draw <- function(ids, count, block) {
			if (is.null(ids)) ids <- if (count == n_iter) seq_len(count) else integer(0)
			k <- match(s, ids)
			if (is.na(k) || k > count) cli::cli_abort("Post-hoc scoring requires matched {block} draws; refit with {.code save_log_lik = TRUE}.")
			k
		}
		actor_draw <- function(name, active, size) {
			if (!active) return(rep(0, size))
			x <- fit[[name]]
			if (is.null(x) || nrow(x) != size) cli::cli_abort("Post-hoc scoring requires saved draws for each active actor margin.")
			k <- match_draw(fit$ab_draw_index, ncol(x), "actor")
			x[, k]
		}
		ractive <- isTRUE(fit$rvar) || (is.null(fit$rvar) && any(fit$APM != 0))
		cactive <- isTRUE(fit$cvar) || (is.null(fit$cvar) && any(fit$BPM != 0))
		a_s <- actor_draw("a_samples", ractive, n_a)
		b_s <- if (isTRUE(fit$symmetric)) a_s else actor_draw("b_samples", cactive, n_b)
		eta_s <- eta_s + outer(a_s, b_s, "+")
		rank_active <- if (bip && !is.null(fit$R_row) && !is.null(fit$R_col)) {
			fit$R_row > 0L && fit$R_col > 0L
		} else if (!is.null(fit$R)) fit$R > 0L else
			!is.null(fit[["U"]]) && length(dim(fit[["U"]])) >= 2L && dim(fit[["U"]])[2L] > 0L
		if (rank_active) {
			if (!have_U) cli::cli_abort("Post-hoc scoring requires matched saved U/V draws; refit with save_log_lik = TRUE.")
			k <- match_draw(fit$uv_draw_index, dim(fit$U_samples)[3L], "factor")
			if (dim(fit$V_samples)[3L] < k) cli::cli_abort("Post-hoc scoring requires matched saved V draws.")
			U_s <- matrix(fit$U_samples[, , k], n_a, dim(fit$U_samples)[2L])
			V_s <- matrix(fit$V_samples[, , k], n_b, dim(fit$V_samples)[2L])
			if (bip) {
				if (is.null(fit$G_samples) || dim(fit$G_samples)[3L] < k)
					cli::cli_abort("Post-hoc bipartite scoring requires matched saved interaction-matrix draws.")
				G_s <- matrix(fit$G_samples[, , k], dim(fit$G_samples)[1L], dim(fit$G_samples)[2L])
				eta_s <- eta_s + U_s %*% G_s %*% t(V_s)
			} else eta_s <- eta_s + tcrossprod(U_s, V_s)
		}
		if (isTRUE(fit$symmetric)) eta_s <- (eta_s + t(eta_s)) / 2
		eta_s
	}
	# per-iteration pointwise log density by family
	y_obs <- Y[obs_idx]
	for (s in seq_len(n_iter)) {
		eta_s <- compute_eta_draw(s)
		ve_s <- fit$VC[s, ve_col]
		eta_obs <- eta_s[obs_idx]

		if (family == "normal") {
			ve_s <- fit$VC[s, ve_col]
			if (!is.finite(ve_s) || ve_s <= 0) {
				cli::cli_abort("Normal log-likelihood requires a finite positive residual variance.")
			}
			ll[s, ] <- stats::dnorm(
				y_obs,
				mean = eta_obs,
				sd = sqrt(ve_s),
				log = TRUE
			)
		} else if (family == "binary" || family == "cbin") {
			# cbin uses the pointwise marginal approximation
			ll[s, ] <- .binary_pointwise_loglik(y_obs, eta_obs)
		} else if (family == "poisson") {
			if (!"ve" %in% colnames(fit$VC))
				cli::cli_abort("Poisson post-hoc scoring requires named observation variance draws.")
			ll[s, ] <- .lame_pois_marginal_loglik(round(y_obs), eta_obs, fit$VC[s, "ve"])
		} else if (family == "ordinal" && identical(fit$ordinal_cutpoints, "explicit")) {
			levels <- fit$ordinal_levels
			if (length(levels) > 2L && (is.null(fit$ALPHA) || nrow(fit$ALPHA) != n_iter || ncol(fit$ALPHA) != length(levels) - 2L))
				cli::cli_abort("Explicit ordinal scoring requires ALPHA rows matched to each coefficient draw.")
			alpha <- if (length(levels) == 2L) 0 else c(0, fit$ALPHA[s, ])
			ll[s, ] <- .lame_explicit_ordinal_loglik(y_obs, eta_obs, alpha, levels)
		} else if (family == "ordinal") {
			ll[s, ] <- .ordinal_pointwise_loglik(y_obs, eta_obs)
		} else {
			# frn falls back to the augmented-Z normal density
			if (s == 1L) {
				cli::cli_warn(c(
					"!" = "The requested rank score is an approximate Gaussian response proxy, not an observed likelihood or a saved-Z density.",
					"i" = "This diagnostic cannot be used for observed-data LOO or WAIC."
				))
			}
			ll[s, ] <- stats::dnorm(
				y_obs,
				mean = eta_obs,
				sd = sqrt(ve_s),
				log = TRUE
			)
		}
	}
	ll
}
