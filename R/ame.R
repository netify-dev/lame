#' AME model fitting routine
#' 
#' An MCMC routine providing a fit to an additive and multiplicative effects
#' (AME) regression model to cross-sectional relational data of various types. 
#' This function supports both unipartite (square) and bipartite (rectangular) 
#' networks. For longitudinal networks, use the \code{lame} function. Original
#' implementation by Peter Hoff.
#' 
#' @details
#' This command provides posterior inference for parameters in AME models of
#' cross-sectional relational data, assuming one of eight possible data types/models.
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
#' \eqn{\Sigma_{uv}} (Inverse-Wishart prior \code{Suv0 / kappa0});
#' the dyad-correlation \eqn{\rho} has a Normal prior on
#' \eqn{\mathrm{atanh}(\rho)}. See \code{prior_summary(fit)} for the
#' priors actually used.
#'
#' \strong{Choosing R.} The multiplicative rank \code{R} controls the
#' dimensionality of latent homophily / heterogeneity not explained
#' by covariates and additive effects. \code{R = 0} fits an
#' additive-only social-relations model; \code{R = 1} or \code{2} is
#' typical for small / medium networks; \code{R > floor(n/3)} is rarely
#' identifiable and will issue a warning. Latent factors capture
#' unobserved structure (clusters, hub patterns, transitive triangles
#' the covariates miss) and are accessed at \code{fit$U}, \code{fit$V}.
#'
#' \strong{Identifiability.} The latent factor term \eqn{u_i' D v_j}
#' is invariant to rotation and reflection of \eqn{U, V}; the package
#' canonicalises with an SVD so successive draws are interpretable.
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
#'   \item Dyadic correlation: \eqn{\rho \sim Uniform(-1, 1)} with Metropolis updates
#' }
#' 
#' The inverse-Wishart prior on \eqn{\Sigma_{ab}} allows learning correlation between
#' sender and receiver effects, capturing reciprocity patterns.
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
#' "tobit": A tobit AME model for censored continuous data. Values are censored
#' at zero, appropriate for non-negative continuous relational data (identity link with censoring).
#' 
#' "binary": A binary probit AME model (probit link: \eqn{P(Y=1) = \Phi(\eta)}).
#' 
#' "ordinal": An ordinal probit AME model (cumulative probit link). An intercept is not 
#' identifiable in this model.
#' 
#' "cbin": An AME model for censored binary data (probit link with censoring). The value of 
#' 'odmax' specifies the maximum number of links each row may have.
#' 
#' "frn": An AME model for fixed rank nomination networks. A higher value of
#' the rank indicates a stronger relationship. The value of 'odmax' specifies
#' the maximum number of links each row may have.
#' 
#' "rrl": An AME model based on the row ranks. This is appropriate if the
#' relationships across rows are not directly comparable in terms of scale. An
#' intercept, row random effects and row regression effects are not estimable
#' for this model.
#' 
#' "poisson": An overdispersed Poisson AME model for count data (log link: \eqn{E[Y] = \exp(\eta)}).
#' The linear predictor \eqn{\eta} represents \eqn{\log(\lambda)} where \eqn{\lambda} is the expected count.
#' 
#' @param Y For unipartite: an n x n square relational matrix. For bipartite: an nA x nB 
#' rectangular relational matrix where nA is the number of row nodes and nB is the 
#' number of column nodes. See family below for various data types.
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
#' @param dcor logical: fit a dyadic correlation (asymmetric case)? Note: not used for bipartite networks.
#' @param nvar logical: fit nodal random effects (symmetric case)?
#' @param R integer: dimension of the multiplicative effects (can be zero). For bipartite networks, 
#' this is used as the default for both R_row and R_col if they are not specified.
#' @param R_row integer: for bipartite networks, dimension of row node multiplicative effects (defaults to R)
#' @param R_col integer: for bipartite networks, dimension of column node multiplicative effects (defaults to R)
#' @param mode character: either "unipartite" (default) for square networks or "bipartite" for rectangular networks.
#' Not all combinations of \code{family} and \code{mode} are supported -- see the table under \strong{Supported
#' family x mode combinations} below.
#' @param family character: one of "normal","tobit","binary","ordinal","cbin","frn","rrl","poisson". See the
#' \strong{Supported family x mode combinations} table below for which combinations are valid; see \strong{Details}
#' for the model assumptions behind each family.
#' @param intercept logical: fit model with an intercept? 
#' @param symmetric logical: Is the sociomatrix symmetric by design?
#' @param odmax a scalar integer or vector of length n giving the maximum
#' number of nominations that each node may make - used for "frn" and "cbin"
#' families
#' @param prior a list containing hyperparameters for the prior distributions.
#' Available options and their defaults:
#' \describe{
#'   \item{Sab0}{Prior covariance matrix for additive effects (default: diag(2)). 
#'         A 2x2 matrix where Sab0\\[1,1\\] is the prior variance for row effects,
#'         Sab0\\[2,2\\] is the prior variance for column effects, and off-diagonals
#'         control correlation between row and column effects.}
#'   \item{eta0}{Prior degrees of freedom for covariance of multiplicative effects 
#'         (default: 4 + 3 \\* n/100, where n is the number of actors). Higher values 
#'         impose stronger shrinkage on the latent factors. Common values: 4-10 for 
#'         weak shrinkage, 10-20 for moderate, >20 for strong shrinkage.}
#'   \item{etaab}{Prior degrees of freedom for covariance of additive effects
#'         (default: 4 + 3 \\* n/100). Controls shrinkage of row/column random effects.
#'         Larger values shrink effects toward zero. \emph{Used by the bipartite
#'         and longitudinal paths; ignored for a unipartite cross-sectional
#'         \code{ame()} fit, which controls the additive prior through
#'         \code{Sab0} and \code{eta0}.}}
#'   \item{s20}{Prior variance for regression coefficients (default: 1).
#'         Larger values allow for larger coefficient values. \emph{Used by the
#'         bipartite and longitudinal paths; for a unipartite \code{ame()} fit
#'         the regression prior is controlled by \code{g}.}}
#'   \item{s2u0}{Prior variance for multiplicative effects (default: 1).
#'         \emph{Used by the bipartite and longitudinal paths.}}
#'   \item{Suv0}{Prior covariance for multiplicative effects (default: identity matrix).}
#' }
#' Common usage: prior = list(Sab0 = diag(c(2, 2)), eta0 = 10) for moderate
#' shrinkage, or prior = list(Sab0 = diag(c(0.5, 0.5))) for tighter control.
#' For a unipartite cross-sectional \code{ame()} fit the prior is controlled by
#' \code{Sab0}, \code{eta0}, \code{Suv0} and \code{g}.
#' @param g optional \strong{scalar} for the Zellner g-prior on regression
#' coefficients (\code{beta ~ N(0, g * sigma^2 * solve(XtX))} where
#' \code{XtX} is the design cross-product). If not specified, defaults are:
#' for \code{normal} family, \code{g = n * var(Y)}; for \code{tobit},
#' \code{g = n * var(Y) * 4}; for other families, \code{g = n} (number of
#' non-missing dyads). Per-coefficient (vector) \code{g} is NOT currently
#' supported by the unipartite path -- pass a scalar.
#' \emph{Note:} \code{g} is a top-level argument to \code{ame()},
#' \strong{not} an element of \code{prior = list(...)}; passing
#' \code{prior = list(g = 0.1)} is a no-op (warned about).
#' @param seed random seed
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
#' @param save_interval quantile interval indicating when to save during post burn-in phase.
#' @param model.name optional string for model selection output
#' @param posterior_opts optional list of posterior sampling options
#' @param n_chains integer: number of MCMC chains to run (default: 1)
#' @param cores integer: number of cores for parallel chains (default: 1)
#' @param use_sparse_matrices logical: use sparse matrix storage for large networks? (default: FALSE).
#'   Recommended only for truly sparse networks (< 10% non-zero entries).
#' @param method character: \code{"mcmc"} (default, the Bayesian MCMC fit) or
#'   \code{"als"} (the fast, MCMC-free iterative block coordinate descent point
#'   estimator). When \code{method = "als"}, MCMC-specific arguments (\code{nscan},
#'   \code{burn}, \code{odens}, \code{prior}, ...) are silently ignored and the
#'   call forwards to \code{\link{ame_als}}.
#' @param bootstrap integer (only used when \code{method = "als"}): number of
#'   bootstrap replicates. \code{0} (default) skips the bootstrap; \code{N > 0}
#'   runs \code{N} replicates and attaches the result so that \code{\link{confint}}
#'   returns bootstrap intervals.
#' @param bootstrap_type character (only used when \code{method = "als"}):
#'   \code{"parametric"} (default) or \code{"block"}.
#' @param bootstrap_block_length integer: block length for the block bootstrap.
#' @param bootstrap_seed optional integer seed for the bootstrap.
#' @param ... reserved for future use. Passing \code{lame()}-only
#'   arguments (e.g. \code{dynamic_beta}, \code{period_exposure}) here
#'   triggers a clean abort directing you to \code{\link{lame}}; passing
#'   any other unrecognised name warns so typos are visible rather than
#'   silently dropped.
#' @param save_log_lik logical: when \code{TRUE}, attach a per-iteration
#'   pointwise log-likelihood matrix \code{fit$log_lik} (an
#'   \code{n_iter x n_obs} array) using the augmented-Z normal
#'   approximation. Required for \code{loo::loo(fit)} and
#'   \code{loo::waic(fit)}. The augmented-Z approximation is the same
#'   semantics used by \code{lame()}; for non-normal families it is the
#'   predictive criterion on the latent-Z scale, not the family-specific
#'   Y density. Default \code{FALSE} (no log_lik storage; default fit is
#'   byte-identical to previous releases).
#' @param ordinal_cutpoints character: cutpoint convention for
#'   \code{family = "ordinal"}. \code{"data_induced"} (default) uses the
#'   data-induced cutpoints; \code{"explicit"} samples explicit cutpoints via
#'   a Cowles (1996) Metropolis-Hastings update. Ignored for other families.
#'
#' @section Supported family x mode combinations:
#' Every \code{family} is supported under
#' both \code{mode}s. The bipartite Z-samplers live in \code{R/rZ_bipartite.R}
#' and dispatch per family; see also the inline comment in \code{R/lame.R}
#' (\dQuote{first-class samplers as of the bipartite-families round}).
#'
#' \tabular{lll}{
#'   \strong{family}  \tab \strong{unipartite} \tab \strong{bipartite} \cr
#'   normal           \tab yes                 \tab yes                \cr
#'   binary           \tab yes                 \tab yes                \cr
#'   tobit            \tab yes                 \tab yes                \cr
#'   ordinal          \tab yes                 \tab yes                \cr
#'   cbin             \tab yes                 \tab yes                \cr
#'   frn              \tab yes                 \tab yes                \cr
#'   poisson          \tab yes                 \tab yes                \cr
#'   rrl              \tab yes                 \tab yes                \cr
#' }
#'
#' Symmetric (\code{symmetric = TRUE}) fits require a symmetric \code{Y}.
#' \code{family = "ordinal"} with \code{symmetric = TRUE} is supported via the
#' dedicated sampler in \code{R/rZ_ord_sym_fc.R}, which uses the
#' symmetric-doubled precision and mirrors upper-triangle draws to the lower
#' triangle so \eqn{Z = t(Z)} holds at every sweep. The previous \dQuote{ordinal
#' + symmetric is unsupported} restriction has been lifted. \code{rrl} +
#' bipartite is first-class: the rectangular row-rank Z-sampler handles the
#' full per-row permutation likelihood, and \code{ame()} recovers regression
#' coefficients comparably to the unipartite rrl path on similarly-sized data.
#'
#' \strong{Symmetric input with one triangle missing.} When
#' \code{symmetric = TRUE} and one triangle of \code{Y} is fully \code{NA}
#' (the user stored only the lower or upper triangle), the symmetry validator
#' \code{any(is.finite(Y - t(Y)))} evaluates \code{FALSE} and the call
#' proceeds: the sampler then treats the populated triangle as the
#' symmetric data and mirrors it. This is usually intended, but if the
#' upper / lower triangles were meant to differ, the model is silently
#' fitting half the data. Audit \code{anyNA(Y[upper.tri(Y)]) !=
#' anyNA(Y[lower.tri(Y)])} before calling if you are unsure.
#'
#' @section Converting a network from another R package:
#' \code{ame()} expects \code{Y} as a (potentially named) matrix, not a graph
#' object. The most common conversions:
#' \itemize{
#'   \item \pkg{igraph}: \code{Y <- igraph::as_adjacency_matrix(g, sparse = FALSE);
#'         diag(Y) <- NA}
#'   \item \pkg{network} / \pkg{statnet}: \code{Y <- as.matrix(net); diag(Y) <- NA}
#'   \item \pkg{tibble} / long-format \pkg{tidyverse}: \code{Y <-
#'         tidyr::pivot_wider(df, names_from = receiver, values_from = tie) |>
#'         tibble::column_to_rownames("sender") |> as.matrix()}; or use the
#'         convenience helper \code{\link{as_lame_y}}.
#' }
#' For an undirected/symmetric network, pass \code{symmetric = TRUE}; for a
#' rectangular two-mode network (students x courses, donors x candidates), pass
#' \code{mode = "bipartite"}. \code{Xrow} and \code{Xcol} accept either a numeric
#' matrix or a data.frame (coerced internally); \code{Xdyad} must be a 3-D array
#' \code{n x n x p} of numeric covariates with no \code{Inf}/\code{NaN}.
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
#'   \code{mutual}                 \tab \code{dcor = TRUE} -> the \code{rho} parameter (probit-scale, NOT log-odds; not numerically comparable to ERGM's \code{mutual}) \cr
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
#' For default cross-sectional calls, \code{lame::ame()} is a drop-in
#' replacement for \code{amen::ame()}: the \code{fit$BETA} slot is a
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
#' and the silent-aggregation hazard with old \code{apply(fit$BETA, 2, mean)}
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
#' @return
#' \strong{Posterior Samples (full MCMC chains):}
#' \item{BETA}{Regression coefficients (nscan xp matrix)}
#' \item{VC}{Variance components (nscan xk matrix)}
#' \item{GOF}{Goodness-of-fit statistics (nscan x4 matrix).
#' First row contains observed values, remaining rows contain posterior predictive samples.
#' See \code{\link{gof}} for post-hoc computation and \code{\link{gof_plot}} for visualization.}
#' 
#' \strong{Posterior Means (averaged over chain):}
#' \item{APM}{Additive row/sender effects (n-vector)}
#' \item{BPM}{Additive column/receiver effects (m-vector); NULL for symmetric networks}
#' \item{U}{Multiplicative row/sender factors (n xR matrix)}
#' \item{V}{Multiplicative column/receiver factors (m xR matrix); NULL for symmetric networks}
#' \item{L}{Eigenvalue matrix (R xR diagonal); symmetric networks only}
#' \item{YPM}{Posterior mean of Y on response scale (for predictions and imputing missing values)}
#' 
#' \strong{Metadata:}
#' \item{family}{Model family (normal, binary, etc.)}
#' \item{mode}{Network mode (unipartite or bipartite)}
#' \item{symmetric}{Logical indicating if network is symmetric}
#' \item{R}{Dimension of multiplicative effects}
#' 
#' \strong{Optional Posterior Samples (if requested via posterior_options):}
#' \item{U_samples}{Samples of U (n xR xiterations array)}
#' \item{V_samples}{Samples of V (m xR xiterations array)}
#' \item{a_samples}{Samples of row effects (n xiterations matrix)}
#' \item{b_samples}{Samples of column effects (m xiterations matrix)}
#' 
#' \strong{Note on reconstructing removed matrices:}
#' To save memory, EZ (expected latent network) and UVPM/ULUPM (multiplicative products) 
#' are not stored but can be reconstructed using:
#' \itemize{
#'   \item \code{reconstruct_EZ(fit)} - Returns linear predictor (link scale, not response scale)
#'   \item \code{reconstruct_UVPM(fit)} - Returns U\%*\%t(V) or U\%*\%L\%*\%t(U)
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
#' data(YX_bin)
#' fit <- ame(YX_bin$Y, Xdyad = YX_bin$X, burn = 10, nscan = 100, odens = 1,
#'            family = "binary", verbose = FALSE)
#' summary(fit)
#' # Note: you should run the Markov chain much longer in practice
#' }
#'  
#' @export ame
ame<-function (
	Y,Xdyad=NULL, Xrow=NULL, Xcol=NULL,
	rvar = !(family=="rrl") , cvar = TRUE, dcor = !symmetric,
	nvar=TRUE,
	R = 0, R_row = NULL, R_col = NULL,
	mode = c("unipartite", "bipartite"),
	family="normal",
	intercept=!is.element(family,c("rrl","ordinal")),
	symmetric=FALSE,
	odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
	prior=list(), g=NA,
	seed = 6886, nscan = 10000, burn = 500, odens = 25,
	verbose = TRUE, gof=TRUE, custom_gof=NULL,
	plot = FALSE,
	start_vals=NULL, periodic_save=FALSE, out_file=NULL,
	save_interval=0.25, model.name=NULL,
	posterior_opts = NULL, n_chains = 1, cores = 1,
	use_sparse_matrices = FALSE,
	method = c("mcmc", "als"),
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	save_log_lik = FALSE,
	ordinal_cutpoints = c("data_induced", "explicit"),
	print, ...
	){
	# resolve ordinal_cutpoints early so the dispatch below carries the
	# canonical scalar value and any non-ordinal-family warning fires once.
	ordinal_cutpoints <- match.arg(ordinal_cutpoints)
	if (identical(ordinal_cutpoints, "explicit") && family != "ordinal") {
		cli::cli_warn(c(
			"{.arg ordinal_cutpoints} = {.val explicit} only takes effect when {.code family = \"ordinal\"}.",
			"i" = "Argument ignored for {.code family = \"{family}\"}."))
		ordinal_cutpoints <- "data_induced"
	}
	# plot is accepted for signature parity with lame(). ame() runs as a
	# single-period mcmc so live trace plotting is not implemented; the
	# argument is silently ignored.
	invisible(plot)

	# detect lame()-only args passed into ame() and abort with a clear
	# explanation rather than R's primitive "unused argument" message.
	# These all assume a longitudinal panel (T >= 2) and have no meaning
	# on a single cross-section.
	dot_names <- names(list(...))
	lame_only <- c("dynamic_beta", "dynamic_beta_kind", "dynamic_beta_pool",
	               "dynamic_beta_per_actor", "dynamic_ab", "dynamic_uv",
	               "dynamic_G", "time_index", "period_exposure",
	               "max_seconds", "checkpoint_path", "checkpoint_every",
	               "log_lik_path", "log_lik_chunk_size", "log_lik_method")
	bad <- intersect(dot_names, lame_only)
	if (length(bad) > 0L) {
		cli::cli_abort(c(
			"{.arg {bad}} {?is/are} longitudinal-only argument{?s} and not valid for {.fn ame} (cross-sectional).",
			"i" = "Use {.fn lame} with a list of {.val T >= 2} matrices to enable dynamic / time-varying effects."))
	}
	# any remaining unrecognised args: emit a single warning rather than
	# silently dropping them. (R's primitive `unused argument` would
	# have aborted; we accept via `...` to give the cleaner message above,
	# so we have to catch typos ourselves.)
	if (length(dot_names) > 0L) {
		unknown <- setdiff(dot_names, lame_only)
		if (length(unknown) > 0L) {
			cli::cli_warn(c(
				"Unrecognised argument{?s} passed to {.fn ame}: {.arg {unknown}}.",
				"i" = "Check {.fn ?ame} for the supported signature; typos here are silently ignored."))
		}
	}

	# handle the deprecated print argument. emit the warning only when
	# the user did not also supply verbose (otherwise their preferred
	# arg already takes effect and the noise just spams multi-chain
	# wrappers like ame_parallel that pass both for back-compat).
	if (!missing(print)) {
		if (missing(verbose)) {
			cli::cli_warn(c(
				"The {.arg print} argument is deprecated.",
				"i" = "Use {.arg verbose} instead."
			))
			verbose <- print
		}
	}

	####
	# single-front-door dispatch: method = "als" forwards to ame_als()
	# (the fast, MCMC-free point estimator). MCMC-specific args don't apply;
	# warn loudly when the user supplied any of them so they know what was
	# ignored (rather than silently believing their MCMC settings took).
	####
	method <- match.arg(method)
	if (identical(method, "als")) {
		mode_arg <- match.arg(mode)
		# detect MCMC-only args the user explicitly supplied
		user_args <- setdiff(names(match.call()), "")
		forwarded <- c("Y", "Xdyad", "Xrow", "Xcol", "R", "R_row", "R_col",
		               "family", "mode", "symmetric",
		               "bootstrap", "bootstrap_type",
		               "bootstrap_block_length", "bootstrap_seed",
		               "verbose", "seed", "method", "intercept", "odmax",
		               # silently absorb cosmetic-only MCMC args
		               "model.name")
		dropped <- setdiff(user_args, forwarded)
		if (length(dropped) > 0L) {
			cli::cli_warn(c(
				"{.fn ame} {.arg method} = {.val als}: ignored {length(dropped)} argument{?s} that apply only to the MCMC path: {.arg {dropped}}.",
				"i" = "These do not affect the ALS fit. Use {.code method = \"mcmc\"} (the default) if you need them."))
		}
		# ame_als() has a single latent dimension R; the bipartite MCMC path
		# accepts asymmetric R_row / R_col. collapse to a single R for the
		# ALS path (max of the two so nothing is truncated) and warn when
		# they differ since ALS can't honour them separately.
		R_als <- R
		if (!is.null(R_row) || !is.null(R_col)) {
			rr <- R_row %||% R; rc <- R_col %||% R
			if (rr != rc) {
				cli::cli_warn(c(
					"{.arg R_row} ({rr}) != {.arg R_col} ({rc}), but the ALS estimator uses a single latent dimension.",
					"i" = "Using {.code R = {max(rr, rc)}} for the ALS fit (the larger of the two)."))
			}
			R_als <- max(rr, rc, R)
		}
		return(ame_als(
			Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
			R = R_als, family = family, mode = mode_arg, symmetric = symmetric,
			bootstrap = bootstrap, bootstrap_type = bootstrap_type,
			bootstrap_block_length = bootstrap_block_length,
			bootstrap_seed = bootstrap_seed,
			verbose = verbose, seed = seed))
	}

	####
	# match call and mode
	####
	mc <- match.call()
	mode <- match.arg(mode)

	# save_log_lik on the cross-sectional path: the underlying samplers do
	# not store a per-iteration log_lik matrix during MCMC, so we compute
	# log_lik post-hoc from saved per-iteration parameters. That requires
	# per-iteration U / V / a / b storage which is opt-in via
	# posterior_opts; force it on here so the log_lik is faithful.
	want_log_lik <- isTRUE(save_log_lik) ||
	                identical(save_log_lik, "chunked")
	if (want_log_lik) {
		if (identical(save_log_lik, "chunked")) {
			cli::cli_warn(c(
				"{.arg save_log_lik} = {.val chunked} is not yet supported on the cross-sectional {.fn ame} path.",
				"i" = "Falling back to in-memory storage."))
		}
		if (is.null(posterior_opts)) posterior_opts <- list()
		# field names match what ame_unipartite/ame_bipartite check:
		# `save_UV` (uppercase) for U/V storage, `save_ab` for a/b storage
		posterior_opts$save_UV <- TRUE
		posterior_opts$save_ab <- TRUE
		if (is.null(posterior_opts$thin_UV)) posterior_opts$thin_UV <- 1L
		if (is.null(posterior_opts$thin_ab)) posterior_opts$thin_ab <- 1L
	}
	# basic identifiability check on the multiplicative rank R. with R too
	# close to the network size the latent factors can absorb the entire
	# residual structure and lose interpretability; warn loudly so the user
	# can lower it.
	if (!is.null(R) && length(R) == 1L && is.finite(R) && R > 0L) {
		n_use <- if (identical(mode, "bipartite")) {
			min(nrow(Y), ncol(Y))
		} else {
			nrow(Y)
		}
		n_third <- max(1L, floor(n_use / 3))
		if (length(n_use) == 1L && is.finite(n_use) && R > n_third) {
			cli::cli_warn(c(
				"!" = "R = {.val {R}} is large relative to n = {.val {n_use}} (recommended cap ~ n/3 = {.val {n_third}}).",
				"i" = "High-rank latent factors can absorb structure that belongs to the additive effects, hurting interpretability.",
				"i" = "Consider {.code R = 2} or {.code R = 3} unless you have a specific reason."))
		}
	}
	# (note: dynamic_beta is intentionally NOT a formal argument here -- the
	# cross-sectional path has only one time period, so an AR(1) prior on beta
	# is unidentified. Use lame() with >= 2 periods to fit dynamic
	# coefficients.)

	# seed locally: restore the global RNG stream on exit so a downstream
	# random draw is not silently perturbed by having fit a model
	if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		.old_seed <- get(".Random.seed", envir = globalenv())
		on.exit(assign(".Random.seed", .old_seed, envir = globalenv()),
		        add = TRUE)
	}

	# validate family; accept short amen spellings as aliases
	valid_families <- c("normal", "tobit", "binary", "ordinal",
	                    "cbin", "frn", "rrl", "poisson")
	amen_aliases <- c(nrm = "normal", bin = "binary", ord = "ordinal",
	                  pois = "poisson", tob = "tobit")
	if (length(family) == 1L && is.character(family) &&
	    family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		cli::cli_inform(c(
			"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
			"i" = "Update your script to {.code family = \"{new_family}\"} when convenient."))
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
	    !family %in% valid_families) {
		cli::cli_abort(c(
			"{.arg family} = {.val {family}} is not a recognised family.",
			"i" = "Choose one of: {.val {valid_families}}."))
	}

	####
	# bipartite ordinal / cbin / frn / poisson / rrl are first-class; the
	# rectangular Z samplers in R/rZ_bipartite.R back the ame_bipartite()
	# dispatch. no guard here.

	####
	# validate Y and the key MCMC arguments up front so degenerate input
	# surfaces as a clear error rather than a deep-stack failure.
	####
	if (!is.matrix(Y) && !is.array(Y)) {
		cli::cli_abort("{.arg Y} must be a matrix.")
	}
	if (!is.numeric(Y) && !is.logical(Y)) {
		cli::cli_abort("{.arg Y} must be numeric.")
	}
	if (any(is.infinite(Y) | is.nan(Y))) {
		cli::cli_abort(c(
			"{.arg Y} contains infinite or NaN values.",
			"i" = "Only finite values and {.val NA} (missing) are allowed."))
	}
	# symmetric=TRUE must be accompanied by a symmetric Y. Otherwise the
	# model silently averages the upper and lower triangles and discards
	# half the observed asymmetry, no warning emitted.
	if (isTRUE(symmetric) && nrow(Y) == ncol(Y)) {
		off <- Y - t(Y)
		if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
				"i" = "Symmetrize {.arg Y} (e.g. {.code (Y + t(Y))/2}) or set {.arg symmetric} = FALSE."))
		}
		# when one full triangle is NA the symmetry check above passes
		# vacuously (Y - t(Y) is all-NA off-diagonal) and the sampler would
		# silently mirror the populated triangle, almost certainly user error.
		ut <- upper.tri(Y); lt <- lower.tri(Y)
		ut_obs <- sum(is.finite(Y[ut]))
		lt_obs <- sum(is.finite(Y[lt]))
		if ((ut_obs == 0L) != (lt_obs == 0L)) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but only one triangle of {.arg Y} is observed (the other is entirely {.val NA}).",
				"i" = "The symmetric sampler would silently mirror the populated triangle.",
				"i" = "Symmetrize {.arg Y} explicitly with {.code Y[lower.tri(Y)] <- t(Y)[lower.tri(Y)]} (or vice versa) before fitting."))
		}
		# warn when symmetric=TRUE is paired with an asymmetric Xdyad slice.
		# The linear predictor eta_ij = beta'X_ij + a_i + b_j + u_i'G v_j is
		# structurally symmetric on the (a, U, V) side (a == b averaged after
		# each sweep; rUV_sym_fc parameterises V = U L for L diagonal) but the
		# regression term inherits whatever symmetry the input X has. The
		# ordinal symmetric sampler (rZ_ord_sym_fc) averages (EZ[i,j] +
		# EZ[j,i])/2 internally to compensate, so non-ordinal symmetric fits
		# carry the most risk of silent misspecification when Xdyad is
		# asymmetric.
		if (!is.null(Xdyad)) {
			Xd_check <- if (is.matrix(Xdyad)) array(Xdyad, c(dim(Xdyad), 1L)) else Xdyad
			if (length(dim(Xd_check)) == 3L) {
				asym_slices <- integer(0)
				for (k in seq_len(dim(Xd_check)[3L])) {
					Xk <- Xd_check[, , k]
					off_k <- Xk - t(Xk)
					if (any(is.finite(off_k)) && max(abs(off_k), na.rm = TRUE) > 1e-8) {
						asym_slices <- c(asym_slices, k)
					}
				}
				if (length(asym_slices) > 0L) {
					slice_nms <- dimnames(Xd_check)[[3]] %||% paste0("k=", asym_slices)
					cli::cli_warn(c(
						"{.arg symmetric} = TRUE, but {length(asym_slices)} {.arg Xdyad} slice{?s} {.val {slice_nms[asym_slices]}} {?is/are} not symmetric.",
						"i" = "The regression term {.code beta'X} will not be structurally symmetric. The sampler relies on the dyadic-correlation coupling (or the explicit symmetric ordinal Z step) to approximate symmetry. For exact structural symmetry, pre-symmetrise: {.code (Xdyad[, , k] + t(Xdyad[, , k]))/2}."))
				}
			}
		}
	}
	obs_y <- Y[is.finite(Y)]
	if (length(obs_y) == 0L) {
		cli::cli_abort("{.arg Y} has no observed (non-missing, finite) values.")
	}
	if (nrow(Y) < 3L || ncol(Y) < 3L) {
		cli::cli_abort(c(
			"{.arg Y} is too small ({nrow(Y)}x{ncol(Y)}) to fit an AME model.",
			"i" = "At least 3 actors per mode are required."))
	}
	# a constant continuous response cannot be fit (the residual variance is
	# zero); a constant discrete response is degenerate but the latent-variable
	# samplers handle it, so it is left to proceed
	if (family %in% c("normal", "tobit") && length(unique(obs_y)) == 1L) {
		cli::cli_abort(c(
			"{.arg Y} is constant; there is no relational variation to model."))
	}
	if (length(R) != 1L || anyNA(R) || R < 0 ||
	    !isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	if (R >= min(nrow(Y), ncol(Y))) {
		cli::cli_abort("{.arg R} = {R} must be smaller than the network dimension ({min(nrow(Y), ncol(Y))}).")
	}
	# bipartite latent dimensions get the same non-negative-integer / range
	# checks as R (otherwise R_row = 1.7 silently recycles a corrupt factor)
	rmax_bip <- min(nrow(Y), ncol(Y))
	for (rr in c("R_row", "R_col")) {
		rv <- get(rr)
		if (!is.null(rv)) {
			if (length(rv) != 1L || anyNA(rv) || rv < 0 ||
			    !isTRUE(all.equal(rv, round(rv)))) {
				cli::cli_abort("{.arg {rr}} must be a single non-negative integer.")
			}
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
			"i" = "Otherwise no posterior draws would be stored."))
	}
	if (family == "poisson" && any(obs_y < 0)) {
		cli::cli_abort(c(
			"{.arg family} = {.val poisson} requires non-negative counts.",
			"i" = "{.arg Y} has negative observed values."))
	}
	# binary/cbin silently threshold via 1*(Y>0); warn so a novice fitting
	# count data with family="binary" sees that their counts were collapsed
	if (family %in% c("binary", "cbin")) {
		uy <- unique(obs_y)
		if (any(!uy %in% c(0, 1))) {
			cli::cli_warn(c(
				"{.arg family} = {.val {family}} but {.arg Y} contains values other than 0/1.",
				"i" = "{.arg Y} will be thresholded to {.code 1 * (Y > 0)}; if you meant counts, use {.val poisson}, or {.val ordinal}/{.val normal} as appropriate."))
		}
	}
	# binary-looking Y with family="normal" is a common novice slip; nudge but
	# do not block (the user may genuinely want a linear-probability model)
	if (family == "normal") {
		uy <- unique(obs_y)
		if (length(uy) <= 2L && all(uy %in% c(0, 1))) {
			cli::cli_inform(c(
				"i" = "{.arg Y} looks binary (values in {.val {{0, 1}}}) but {.arg family} = {.val normal}.",
				"i" = "If this is a 0/1 network, pass {.code family = \"binary\"}."))
		}
	}
	# warn when burn dwarfs nscan -- typical use is burn << nscan; burn > nscan
	# wastes most of the chain
	if (burn > nscan) {
		cli::cli_warn(c(
			"{.arg burn} = {burn} is greater than {.arg nscan} = {nscan}.",
			"i" = "Most iterations will be discarded; consider {.arg burn} <= {.arg nscan}."))
	}
	# fully unobserved actors (all-NA rows / cols) silently get APM/BPM purely
	# from the prior; surface this so an upstream data-pipeline bug is visible
	if (mode == "unipartite") {
		drow <- apply(Y, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y, 2, function(c) all(is.na(c) | !is.finite(c)))
		# drop the diagonal NA from the "all NA" judgement
		drow_ex <- vapply(seq_len(nrow(Y)), function(i)
			all(is.na(Y[i, -i]) | !is.finite(Y[i, -i])), logical(1))
		dcol_ex <- vapply(seq_len(ncol(Y)), function(j)
			all(is.na(Y[-j, j]) | !is.finite(Y[-j, j])), logical(1))
		if (any(drow_ex) || any(dcol_ex)) {
			cli::cli_inform(c(
				"i" = "{.arg Y} has {sum(drow_ex)} fully-unobserved row{?s} and {sum(dcol_ex)} fully-unobserved column{?s}.",
				"i" = "Their additive effects (APM/BPM) will be estimated entirely from the prior."))
		}
	} else {
		drow <- apply(Y, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y, 2, function(c) all(is.na(c) | !is.finite(c)))
		# the bipartite-binary / tobit Z sampler cannot proceed when a row
		# or column has no observed data; abort with a clear remediation
		# hint rather than letting the sampler crash deeper in the stack.
		if ((any(drow) || any(dcol)) && family %in% c("binary", "tobit")) {
			cli::cli_abort(c(
				"{.arg Y} (bipartite, {.val {family}}) has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
				"x" = "The bipartite {.val {family}} sampler cannot proceed with rows/columns that have no observed data.",
				"i" = "Drop the empty actor(s) from {.arg Y} (and matching covariates) before fitting, or use {.val normal} family."))
		}
		if (any(drow) || any(dcol)) {
			cli::cli_inform(c(
				"i" = "{.arg Y} has {sum(drow)} fully-unobserved row{?s} and {sum(dcol)} fully-unobserved column{?s}.",
				"i" = "Their additive effects (APM/BPM) will be estimated entirely from the prior."))
		}
	}
	# covariate finiteness and shape checks; surface bad inputs early
	# instead of letting them blow up deep in the sampler.
	.check_cov_finite <- function(x, nm) {
		if (is.null(x)) return(invisible())
		v <- if (is.data.frame(x)) as.matrix(x) else x
		if (!is.numeric(v) && !is.logical(v)) {
			cli::cli_abort("{.arg {nm}} must be numeric.")
		}
		if (any(is.infinite(v) | is.nan(v))) {
			cli::cli_abort(c(
				"{.arg {nm}} contains infinite or NaN values.",
				"i" = "Only finite values and {.val NA} are allowed."))
		}
		invisible()
	}
	.check_cov_finite(Xrow, "Xrow")
	.check_cov_finite(Xcol, "Xcol")
	.check_cov_finite(Xdyad, "Xdyad")
	# warn when Y and Xdyad both carry dimnames but in different orders --
	# the package aligns positionally, so a name mismatch would silently
	# misalign the whole design and can flip coefficient signs.
	if (!is.null(Xdyad) && length(dim(Xdyad)) == 3L &&
	    !is.null(rownames(Y)) && !is.null(dimnames(Xdyad)[[1]])) {
		yr <- rownames(Y); xr <- dimnames(Xdyad)[[1]]
		yc <- colnames(Y); xc <- dimnames(Xdyad)[[2]]
		if (length(yr) == length(xr) && !identical(yr, xr)) {
			if (setequal(yr, xr)) {
				cli::cli_warn(c(
					"{.arg Y} and {.arg Xdyad} have the same row actors but in different orders.",
					"i" = "Reorder {.arg Xdyad} to match {.arg Y}'s rownames or the design will be silently misaligned.",
					"i" = "Try {.code Xdyad <- Xdyad[rownames(Y), colnames(Y), , drop = FALSE]}."))
			} else {
				cli::cli_warn(c(
					"{.arg Y} and {.arg Xdyad} have different row actors -- alignment is positional.",
					"i" = "If actors should match by name, reorder {.arg Xdyad} so its rownames equal {.code rownames(Y)}."))
			}
		}
		if (!is.null(yc) && !is.null(xc) &&
		    length(yc) == length(xc) && !identical(yc, xc) && setequal(yc, xc)) {
			cli::cli_warn(c(
				"{.arg Y} and {.arg Xdyad} have the same column actors but in different orders.",
				"i" = "Reorder columns of {.arg Xdyad} to match {.code colnames(Y)}."))
		}
	}

	# Xdyad should be a 3-D array n x n x p (or nA x nB x p for bipartite).
	# a 2-D n x n matrix is a common novice slip; coerce with a message so
	# the variable name is preserved and the user knows the wrap happened.
	if (!is.null(Xdyad)) {
		dimd <- dim(Xdyad)
		if (is.null(dimd) || length(dimd) == 2L) {
			if (is.null(dimd)) {
				cli::cli_abort(c(
					"{.arg Xdyad} must be a 2-D matrix or 3-D array, not a {.cls {class(Xdyad)[1]}}."))
			}
			if (dimd[1] != nrow(Y) || dimd[2] != ncol(Y)) {
				cli::cli_abort(c(
					"{.arg Xdyad} dimensions ({dimd[1]} x {dimd[2]}) do not match {.arg Y} ({nrow(Y)} x {ncol(Y)})."))
			}
			xname <- deparse(substitute(Xdyad))[1]
			# only flag when this looks like a real variable, not an inline expr
			if (nchar(xname) > 0L && nchar(xname) < 40L &&
			    !grepl("[(){}]", xname)) {
				cli::cli_inform(c(
					"i" = "{.arg Xdyad} was a 2-D matrix; wrapping as {.code array(., dim = c(n, n, 1))} with slice name {.val {xname}}."))
			}
			Xdyad <- array(Xdyad, dim = c(dimd[1], dimd[2], 1L),
			               dimnames = list(rownames(Xdyad), colnames(Xdyad), xname))
			dimd <- dim(Xdyad)
		} else if (length(dimd) != 3L) {
			cli::cli_abort(c(
				"{.arg Xdyad} must be a 3-dimensional array ({nrow(Y)} x {ncol(Y)} x p).",
				"i" = "Got a {length(dimd)}-D object; wrap a single covariate as {.code array(x, dim = c(n, n, 1))}."))
		}
		if (dimd[1] != nrow(Y) || dimd[2] != ncol(Y)) {
			cli::cli_abort(c(
				"{.arg Xdyad} dimensions ({dimd[1]} x {dimd[2]} x {dimd[3]}) do not match {.arg Y} ({nrow(Y)} x {ncol(Y)})."))
		}
	}
	# warn on prior-list keys we don't recognise (silent typos otherwise mean
	# the user's "custom prior" was actually the default)
	if (length(prior) > 0L && !is.null(names(prior))) {
		known_prior <- c("Sab0", "eta0", "etaab", "Suv0", "kappa0",
		                 "s20", "s2u0",
		                 "rho_uv_mean", "rho_uv_sd",
		                 "rho_ab_mean", "rho_ab_sd")
		bad <- setdiff(names(prior), known_prior)
		if (length(bad) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg prior} list element{?s}: {.val {bad}}.",
				"i" = "Recognised names: {.val {known_prior}}.",
				"i" = "Note: {.code g} is a top-level argument to {.fn ame}, not a {.arg prior} entry."))
		}
	}
	if (isTRUE(periodic_save) && is.null(out_file)) {
		cli::cli_abort(c(
			"{.arg periodic_save} = TRUE requires an {.arg out_file} path.",
			"i" = "Supply {.arg out_file} or set {.arg periodic_save} = FALSE."))
	}

	####
	# parallel chains dispatch
	####
	if(n_chains > 1) {
		return(ame_parallel(Y = Y, n_chains = n_chains, cores = cores,
											 combine_method = "pool",
											 Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
											 rvar = rvar, cvar = cvar, dcor = dcor, nvar = nvar,
											 R = R, R_row = R_row, R_col = R_col,
											 mode = mode, family = family, 
											 intercept = intercept, symmetric = symmetric,
											 odmax = odmax, prior = prior, g = g,
											 seed = seed, nscan = nscan, burn = burn, odens = odens,
											 verbose = verbose, gof = gof,
											 start_vals = start_vals, periodic_save = periodic_save,
											 out_file = out_file, save_interval = save_interval,
											 model.name = model.name, posterior_opts = posterior_opts))
	}

	####
	# unipartite dispatch
	####
	if(mode == "unipartite") {
		if(nrow(Y) != ncol(Y)) {
			cli::cli_abort(c(
				"Matrix Y must be square for unipartite networks.",
				"i" = "Y has dimensions {nrow(Y)} x {ncol(Y)}.",
				"i" = "Use mode='bipartite' for rectangular matrices."
			))
		}
		# symmetric ordinal is first-class. the sampler in R/rZ_ord_sym_fc.R
		# uses the symmetric-doubled precision (variance 1/2) and mirrors
		# upper-triangle draws to the lower triangle so Z = t(Z) holds at
		# every sweep.
		
		fit <- ame_unipartite(
			Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol,
			rvar=rvar, cvar=cvar, dcor=dcor, nvar=nvar, R=R,
			family=family, intercept=intercept, symmetric=symmetric,
			odmax=odmax, prior=prior, g=g,
			seed=seed, nscan=nscan, burn=burn, odens=odens,
			verbose=verbose, gof=gof, custom_gof=custom_gof,
			start_vals=start_vals, periodic_save=periodic_save,
			out_file=out_file, save_interval=save_interval,
			model.name=model.name, posterior_opts=posterior_opts,
			use_sparse_matrices=use_sparse_matrices,
			ordinal_cutpoints=ordinal_cutpoints
		)
		fit$call <- mc

	####
	# bipartite dispatch
	####
	} else if(mode == "bipartite") {
		if(symmetric) {
			cli::cli_warn("Symmetric option ignored for bipartite networks.")
			symmetric <- FALSE
		}
		
		if(is.null(R_row)) R_row <- R
		if(is.null(R_col)) R_col <- R
		
		fit <- ame_bipartite(
			Y=Y, Xdyad=Xdyad, Xrow=Xrow, Xcol=Xcol,
			rvar=rvar, cvar=cvar, R_row=R_row, R_col=R_col,
			family=family, intercept=intercept,
			odmax=odmax, prior=prior, g=g,
			seed=seed, nscan=nscan, burn=burn, odens=odens,
			verbose=verbose, gof=gof, custom_gof=custom_gof,
			start_vals=start_vals, periodic_save=periodic_save,
			out_file=out_file, save_interval=save_interval,
			model.name=model.name, posterior_opts=posterior_opts,
			use_sparse_matrices=use_sparse_matrices
		)
		fit$call <- mc
	}

	# post-hoc log_lik computation when save_log_lik = TRUE. uses the
	# augmented-Z normal approximation (same semantics as lame()): for
	# non-normal families this is the predictive criterion on the latent
	# Z scale, not the family-specific Y density. matches the documented
	# behaviour of loo.ame and is sufficient for relative model comparison.
	if (want_log_lik) {
		fit$log_lik <- .ame_compute_log_lik_post_hoc(fit, Y = Y, family = family)
		fit$log_lik_method <- if (family %in% c("frn", "rrl")) "augmented" else "observed_exact"
	}

	return(fit)
}

# post-hoc per-iteration log_lik for ame() fits. relies on
# posterior_opts$store_uv / store_ab having stored U_samples / V_samples /
# a_samples / b_samples; falls back to posterior means when storage was
# disabled (rare: only triggered for chunked path on bipartite, which
# warns elsewhere).
#' @noRd
.ame_compute_log_lik_post_hoc <- function(fit, Y, family) {
	if (is.null(fit$BETA) || is.null(fit$VC)) {
		cli::cli_warn(c("Cannot compute log_lik: BETA / VC missing from fit."))
		return(NULL)
	}
	n_iter <- nrow(fit$BETA)
	# observations: off-diagonal finite cells of Y
	if (is.array(Y) && length(dim(Y)) == 2L) {
		if (nrow(Y) == ncol(Y)) diag(Y) <- NA
	}
	if (identical(family, "binary"))  Y <- 1 * (Y > 0)
	if (identical(family, "cbin"))    Y <- 1 * (Y > 0)
	if (identical(family, "tobit"))   Y[Y < 0] <- 0
	obs_idx <- which(!is.na(Y) & is.finite(Y))
	n_obs <- length(obs_idx)
	if (n_obs == 0L) return(NULL)
	# pre-compute per-cell design components that don't change per iter
	n_a <- nrow(Y); n_b <- ncol(Y)
	# locate sigma^2 column on VC: last column is "ve" by convention.
	ve_col <- ncol(fit$VC)
	# beta_names: tells us which beta columns are dyadic vs intercept / nodal
	beta_names <- colnames(fit$BETA)
	has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"
	# pull out per-iteration U / V / a / b if available; else fall back to
	# posterior-mean U / V / APM / BPM (the latter understates per-draw spread
	# but at least lets loo() run on a fit where the samples were not stored)
	have_U <- !is.null(fit$U_samples) && !is.null(fit$V_samples)
	have_ab <- !is.null(fit$a_samples) && !is.null(fit$b_samples)
	# build a per-cell Xdyad x BETA matrix once. Xdyad is stored on fit$X
	# when ame() saves the call args; rebuild from the call if missing.
	X <- fit$X
	ll <- matrix(NA_real_, nrow = n_iter, ncol = n_obs)
	# linear predictor for draw s
	compute_eta_draw <- function(s) {
		beta_s <- if (length(dim(fit$BETA)) == 3L) fit$BETA[s, , 1] else fit$BETA[s, ]
		eta_s <- matrix(0, n_a, n_b)
		if (has_intercept) eta_s <- eta_s + beta_s[1]
		# dyadic covariates
		if (!is.null(X)) {
			# X is [n_a, n_b, p_dyad]; map beta_s[non-intercept dyad slots] -> X
			is_dyad <- if (is.null(beta_names)) {
				c(FALSE, rep(TRUE, length(beta_s) - 1L))
			} else {
				grepl("_dyad$|\\.dyad$", beta_names)
			}
			dyad_b <- beta_s[is_dyad]
			p_dyad <- if (length(dim(X)) == 3L) dim(X)[3L] else 1L
			for (k in seq_len(min(p_dyad, length(dyad_b)))) {
				eta_s <- eta_s + dyad_b[k] *
					(if (length(dim(X)) == 3L) X[, , k] else X)
			}
		}
		# additive effects
		if (have_ab) {
			ab_idx <- min(s, dim(fit$a_samples)[2L])
			a_s <- fit$a_samples[, ab_idx]
			b_s <- if (!is.null(fit$b_samples)) fit$b_samples[, ab_idx] else a_s
		} else {
			a_s <- fit$APM %||% rep(0, n_a)
			b_s <- fit$BPM %||% rep(0, n_b)
		}
		eta_s <- eta_s + outer(a_s, b_s, "+")
		# multiplicative effects
		if (have_U) {
			uv_idx <- min(s, dim(fit$U_samples)[3L])
			U_s <- fit$U_samples[, , uv_idx]
			V_s <- fit$V_samples[, , uv_idx]
			eta_s <- eta_s + U_s %*% t(V_s)
		} else if (!is.null(fit$U) && !is.null(fit$V) &&
		           NCOL(fit$U) > 0L && NCOL(fit$V) > 0L) {
			eta_s <- eta_s + fit$U %*% t(fit$V)
		}
		eta_s
	}
	# per-iteration family-specific pointwise log-Y density
	y_obs <- Y[obs_idx]
	for (s in seq_len(n_iter)) {
		eta_s <- compute_eta_draw(s)
		ve_s <- max(fit$VC[s, ve_col], .Machine$double.eps)
		eta_obs <- eta_s[obs_idx]

		if (family == "normal") {
			ll[s, ] <- stats::dnorm(y_obs, mean = eta_obs,
			                        sd = sqrt(ve_s), log = TRUE)
		} else if (family == "binary" || family == "cbin") {
			p_hat <- stats::pnorm(eta_obs)
			p_hat <- pmin(pmax(p_hat, 1e-12), 1 - 1e-12)
			ll[s, ] <- y_obs * log(p_hat) + (1 - y_obs) * log(1 - p_hat)
		} else if (family == "tobit") {
			sigma_t <- sqrt(ve_s)
			is_zero <- y_obs <= 0
			ll[s, is_zero]  <- stats::pnorm(0, mean = eta_obs[is_zero],
			                                sd = sigma_t, log.p = TRUE)
			ll[s, !is_zero] <- stats::dnorm(y_obs[!is_zero],
			                                mean = eta_obs[!is_zero],
			                                sd = sigma_t, log = TRUE)
		} else if (family == "poisson") {
			lam <- pmin(exp(eta_obs), 1e8)
			ll[s, ] <- stats::dpois(round(y_obs), lambda = lam, log = TRUE)
		} else if (family == "ordinal") {
			ll[s, ] <- .ordinal_pointwise_loglik(y_obs, eta_obs)
		} else {
			# augmented-Z normal fallback for frn / rrl
			if (s == 1L) {
				cli::cli_warn(c(
					"!" = "Pointwise log-lik for family {.val {family}} uses the augmented-Z normal approximation.",
					"i" = "Treat {.fn loo::loo} output as latent-space relative comparison only; the exact marginal needs GHK and is not implemented for the cross-sectional path."))
			}
			ll[s, ] <- stats::dnorm(y_obs, mean = eta_obs,
			                        sd = sqrt(ve_s), log = TRUE)
		}
	}
	ll
}