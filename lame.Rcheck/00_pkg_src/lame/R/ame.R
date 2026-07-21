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
#' \eqn{\Sigma_{uv}} (Inverse-Wishart prior with scale \code{kappa0 * Suv0}
#' and \code{kappa0} degrees of freedom);
#' the dyad-correlation \eqn{\rho} has an arc-sine prior on \eqn{(-1, 1)}
#' (density proportional to \eqn{(1-\rho^2)^{-1/2}}), updated with
#' Metropolis steps. See \code{prior_summary(fit)} for the priors
#' actually used.
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
#'   \item Dyadic correlation: arc-sine prior on \eqn{(-1, 1)}, density
#'         \eqn{p(\rho) \propto (1-\rho^2)^{-1/2}}, with Metropolis updates
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
#' \code{eta0 = 4}. A bipartite fit always uses \code{Sab0 = diag(2)}.
#' Call \code{prior_summary()} on a fitted object to see the prior actually used.
#'
#' The data scaling matters because the inverse-Wishart prior contributes
#' pseudo-data on the scale of \code{eta0 * Sab0} no matter what units \code{Y}
#' is in, so a fixed \code{Sab0 = diag(2)} pulls the additive-effects variances
#' \code{va}/\code{vb} upward whenever the true sender/receiver variance is well
#' below 1, most visibly at small \eqn{n}. In Social Relations Model simulations
#' the data-scaled default recovers small true variances with modest bias where
#' a fixed \code{diag(2)} prior can nearly double them; the data-scaled prior is
#' mildly more conservative in the opposite regime, when the true variance is
#' large relative to the residual scale.
#'
#' This default differs from the \code{amen} package, whose
#' \code{ame(family = "nrm")} leaves \code{Sab0 = diag(2)} and \code{eta0 = 4}.
#' To reproduce \code{amen}'s additive-effects posterior, pass
#' \code{prior = list(Sab0 = diag(2), eta0 = 4)}. The choice affects only
#' \code{va}, \code{vb} and (weakly) \code{cab}; \code{ve} and \code{rho} are
#' unchanged.
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
#' "poisson": An overdispersed Poisson AME model for count data:
#' \eqn{Y \sim \mathrm{Poisson}(\exp(z))} with \eqn{z \sim N(\eta, \sigma^2)},
#' a lognormal-mixed Poisson. The conditional mean given the latent \eqn{z} is
#' \eqn{\exp(z)}; the marginal mean is \eqn{\exp(\eta + \sigma^2/2)}, not
#' \eqn{\exp(\eta)}.
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
#' @param dcor logical: fit a dyadic correlation (asymmetric case)? Note: not used for bipartite networks.
#' @param nvar logical: fit nodal random effects (symmetric case)?
#' @param R integer: dimension of the multiplicative effects (can be zero). For bipartite networks, 
#' this is used as the default for both R_row and R_col if they are not specified.
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
#' families
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
#'         bipartite fit defaults to \code{diag(2)}. See the note on the
#'         additive-effects variance scale in Details. Pass \code{diag(2)}
#'         explicitly for a fixed unit-scale prior.}
#'   \item{eta0}{Prior degrees of freedom for the additive-effects covariance
#'         \eqn{\Sigma_{ab}} (default: round(4 + 3 \\* n/100) for a unipartite
#'         continuous-family fit, round(4 \\* vdfmlt) for a unipartite
#'         discrete-family fit, and 4 + 3 \\* (nA + nB) / 200 for a bipartite
#'         fit -- except bipartite \code{"binary"}, which also uses
#'         round(4 \\* vdfmlt) -- where n is the number of actors and vdfmlt
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
#'   \item{s20}{Prior variance for regression coefficients (default: 1).
#'         Larger values allow for larger coefficient values. \emph{Used by the
#'         bipartite and longitudinal paths; for a unipartite \code{ame()} fit
#'         the regression prior is controlled by \code{g}.}}
#'   \item{s2u0}{Prior variance for multiplicative effects (default: 1).
#'         \emph{Used by the bipartite and longitudinal paths.}}
#'   \item{Suv0}{Inverse-Wishart (inverse) scale matrix for the
#'         multiplicative-effects covariance. The prior scale is
#'         \code{kappa0 * Suv0} and the prior mean is
#'         \code{kappa0 * Suv0 / (kappa0 - 2 * R - 1)}. When not supplied it
#'         defaults to \code{diag(2 * R)} times a scale estimated from the
#'         data.}
#'   \item{kappa0}{Prior degrees of freedom for the multiplicative-effects
#'         covariance \eqn{\Sigma_{uv}} (default when R > 0: 2 \\* R + 2, times
#'         the same probit-moment multiplier as \code{eta0} for a unipartite
#'         discrete-family fit). Higher values impose stronger shrinkage of the
#'         latent factors toward the prior scale \code{kappa0 * Suv0}.}
#' }
#' Common usage: prior = list(Sab0 = diag(c(2, 2)), eta0 = 10) for moderate
#' shrinkage, or prior = list(Sab0 = diag(c(0.5, 0.5))) for tighter control.
#' For a unipartite cross-sectional \code{ame()} fit the prior is controlled by
#' \code{Sab0}, \code{eta0}, \code{Suv0} and \code{g}.
#' @param g optional \strong{scalar} for the Zellner g-prior on regression
#' coefficients (\code{beta ~ N(0, g * sigma^2 * solve(XtX))} where
#' \code{XtX} is the design cross-product). If not specified, defaults are:
#' for \code{normal} family, \code{g = n * var(Y)}; for other families,
#' \code{g = n} (number of
#' non-missing dyads). Per-coefficient (vector) \code{g} is not currently
#' supported by the unipartite path -- pass a scalar.
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
#'   additive-effect draws on \code{fit$a_samples} / \code{fit$b_samples}.
#'   \code{thin_UV} / \code{thin_ab} thin the stored draws. Default
#'   \code{NULL}.
#' @param n_chains integer: number of MCMC chains to run (default: 1)
#' @param cores integer: number of cores for parallel chains (default: 1)
#' @param use_sparse_matrices logical: use sparse matrix storage for large networks? (default: FALSE).
#'   Recommended only for truly sparse networks (< 10% non-zero entries).
#' @param method character: \code{"mcmc"} (default, the Bayesian MCMC fit) or
#'   \code{"als"} (the fast, MCMC-free iterative block coordinate descent point
#'   estimator). When \code{method = "als"}, MCMC-specific arguments (\code{nscan},
#'   \code{burn}, \code{odens}, \code{prior}, ...) are warned about and ignored;
#'   the call forwards to \code{\link{ame_als}}.
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
#'   \code{n_iter x n_obs} array), required for \code{loo(fit)} and
#'   \code{waic(fit)}. For families normal, binary, cbin, poisson,
#'   and ordinal this is the exact observed-data log density of Y
#'   (\code{fit$log_lik_method = "observed_exact"}); only frn falls
#'   back to the augmented-Z normal approximation on the latent scale
#'   (with a one-time warning). Default \code{FALSE} (no log_lik storage;
#'   default fit is byte-identical to previous releases).
#' @param ordinal_cutpoints character: cutpoint convention for
#'   \code{family = "ordinal"}. \code{"data_induced"} (default) uses the
#'   data-induced cutpoints; \code{"explicit"} samples explicit cutpoints via
#'   a Cowles (1996) Metropolis-Hastings update. Ignored for other families.
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
#' Symmetric (\code{symmetric = TRUE}) fits require a symmetric \code{Y}.
#' \code{family = "ordinal"} with \code{symmetric = TRUE} is supported via the
#' dedicated sampler in \code{R/rZ_ord_sym_fc.R}, which uses the
#' symmetric-doubled precision and mirrors upper-triangle draws to the lower
#' triangle so \eqn{Z = t(Z)} holds at every sweep.
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
#' @section Data preparation:
#' \code{ame()} accepts a matrix directly. For long-format edgelists,
#' bipartite data, and covariates, use \pkg{netify} to build the network object
#' and pass that object as \code{Y}. \code{ame()} will call
#' \code{netify::to_lame()} internally. If you already have an \pkg{igraph} or
#' \pkg{network} object and only need a plain adjacency matrix,
#' \code{\link{as_lame_y}} is still available as a small convenience helper.
#'
#' For an undirected/symmetric network, pass \code{symmetric = TRUE}; for a
#' rectangular two-mode network (students x courses, donors x candidates), pass
#' \code{mode = "bipartite"} or build the \pkg{netify} object with
#' \code{mode = "bipartite"}. \code{Xrow} and \code{Xcol} accept either a
#' numeric matrix or a data.frame (coerced internally); \code{Xdyad} must be a
#' 3-D array \code{n x n x p} of numeric covariates with no
#' \code{Inf}/\code{NaN}.
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
#' @return
#' \strong{Posterior Samples (full MCMC chains):}
#' \item{BETA}{Regression coefficients (\code{nscan/odens} x p matrix; one
#' row per stored draw)}
#' \item{VC}{Variance components (\code{nscan/odens} x k matrix)}
#' \item{GOF}{Goodness-of-fit statistics ((\code{nscan/odens} + 1) x 5 matrix).
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
	rvar = TRUE, cvar = TRUE, dcor = !symmetric,
	nvar=TRUE,
	R = 0, R_row = NULL, R_col = NULL,
	mode = c("unipartite", "bipartite"),
	family="normal",
	intercept=!(family == "ordinal"),
	symmetric=FALSE,
	odmax=rep(max(apply(Y>0,1,sum,na.rm=TRUE)),nrow(Y)),
	prior=list(), g=NA,
	seed = 6886, nscan = 10000, burn = 500, odens = 25,
	verbose = TRUE, gof=TRUE, custom_gof=NULL,
	plot = FALSE,
	start_vals=NULL, periodic_save=FALSE, out_file=NULL,
	save_interval=0.25,
	posterior_opts = NULL, n_chains = 1, cores = 1,
	use_sparse_matrices = FALSE,
	method = c("mcmc", "als"),
	bootstrap = 0L,
	bootstrap_type = c("parametric", "block"),
	bootstrap_block_length = 1L,
	bootstrap_seed = NULL,
	save_log_lik = FALSE,
	ordinal_cutpoints = c("data_induced", "explicit"),
	print, ...,
	# after `...` so a typo like `model = "bin"` cannot partial-match
	# into it and silently fit the default family
	model.name=NULL
	){
	family_missing <- missing(family)
	mode_missing <- missing(mode)
	symmetric_missing <- missing(symmetric)
	netify_input <- .lame_prepare_netify_input(
		Y = Y, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol,
		mode = mode, family = family,
		mode_missing = mode_missing,
		family_missing = family_missing,
		symmetric = symmetric,
		symmetric_missing = symmetric_missing,
		longitudinal = FALSE,
		caller = "ame",
		fit_method = if (identical(method[1L], "als")) "als" else "gibbs")
	if (isTRUE(netify_input$converted)) {
		Y <- netify_input$Y
		Xdyad <- netify_input$Xdyad
		Xrow <- netify_input$Xrow
		Xcol <- netify_input$Xcol
		mode <- netify_input$mode
		family <- netify_input$family
		symmetric <- netify_input$symmetric
	}

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
	# explanation rather than r's primitive "unused argument" message.
	# these all assume a longitudinal panel and have no meaning
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
	# silently dropping them. (r's primitive `unused argument` would
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
	# route als fits through the point estimator
	####
	method <- match.arg(method)
	if (identical(method, "als")) {
		mode_arg <- match.arg(mode)
		# skip the intercept guard for the ordinal family (and its "ord"
		# alias): its default sets intercept = FALSE internally, and the ALS
		# driver aborts with the real reason (ordinal is not an ALS family),
		# so blaming intercept here would misdirect the user.
		if (!isTRUE(intercept) &&
		    !identical(family, "ordinal") && !identical(family, "ord")) {
			cli::cli_abort(c(
				"{.fn ame} {.arg method} = {.val als} does not currently support {.code intercept = FALSE}.",
				"i" = "Use {.code method = \"mcmc\"} for a no-intercept fit, or keep the ALS intercept and center covariates as needed."))
		}
		# detect mcmc-only args the user explicitly supplied
		user_args <- setdiff(names(match.call()), "")
		forwarded <- c("Y", "Xdyad", "Xrow", "Xcol", "R", "R_row", "R_col",
		               "family", "mode", "symmetric",
		               "bootstrap", "bootstrap_type",
		               "bootstrap_block_length", "bootstrap_seed",
		               "verbose", "seed", "method", "intercept", "odmax",
		               # silently absorb cosmetic-only mcmc args
		               "model.name")
		dropped <- setdiff(user_args, forwarded)
		if (length(dropped) > 0L) {
			cli::cli_warn(c(
				"{.fn ame} {.arg method} = {.val als}: ignored {length(dropped)} arguments that apply only to the MCMC path: {.arg {dropped}}.",
				"i" = "These do not affect the ALS fit. Use {.code method = \"mcmc\"} (the default) if you need them."))
		}
		# collapse bipartite ranks for the static als route
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

	# posterior_opts: warn on unknown keys (usually a typo, and silently
	# dropping a requested storage option is worse than noise), and accept
	# the longitudinal save_UV_draws spelling (what lame() uses) as an
	# alias for save_UV so the same option list behaves consistently
	# across ame() and lame().
	if (!is.null(posterior_opts)) {
		if (!is.list(posterior_opts)) {
			cli::cli_abort(c(
				"{.arg posterior_opts} must be a list.",
				"i" = "Build it with {.fn posterior_options}, e.g. {.code posterior_opts = posterior_options(save_UV = TRUE)}."))
		}
		known_posterior_opts <- c("save_UV", "save_UV_draws", "save_ab",
		                          "thin_UV", "thin_ab")
		bad_opts <- setdiff(names(posterior_opts), known_posterior_opts)
		if (length(bad_opts) > 0L) {
			cli::cli_warn(c(
				"Unknown {.arg posterior_opts} list element{?s}: {.val {bad_opts}}.",
				"i" = "Recognised names: {.val {known_posterior_opts}}."))
		}
		if (isTRUE(posterior_opts$save_UV_draws)) {
			posterior_opts$save_UV <- TRUE
		}
	}

	# save_log_lik on the cross-sectional path: the underlying samplers do
	# not store a per-iteration log_lik matrix during mcmc, so we compute
	# log_lik post-hoc from saved per-iteration parameters. that requires
	# per-iteration u / v / a / b storage which is opt-in via
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
		# `save_uv` (uppercase) for u/v storage, `save_ab` for a/b storage
		posterior_opts$save_UV <- TRUE
		posterior_opts$save_ab <- TRUE
		# the pointwise log-lik pairs BETA draw s with a_samples/U_samples
		# draw s, so thinned effect storage would silently mis-pair draws
		# and corrupt WAIC/LOO -- force 1:1 storage
		if (!is.null(posterior_opts$thin_UV) && posterior_opts$thin_UV != 1L ||
		    !is.null(posterior_opts$thin_ab) && posterior_opts$thin_ab != 1L) {
			cli::cli_warn(c(
				"{.arg save_log_lik} requires unthinned effect draws; overriding {.code thin_UV}/{.code thin_ab} to 1.",
				"i" = "The pointwise log-likelihood pairs each stored BETA draw with its own U/V and a/b draw."))
		}
		posterior_opts$thin_UV <- 1L
		posterior_opts$thin_ab <- 1L
	}
	# basic identifiability check on the multiplicative rank r. with r too
	# close to the network size the latent factors can absorb the entire
	# residual structure and lose interpretability; warn loudly so the user
	# can lower it.
	rank_check_max <- suppressWarnings(max(
		if (!is.null(R) && length(R) == 1L && is.finite(R)) R else 0L,
		if (!is.null(R_row) && length(R_row) == 1L && is.finite(R_row)) R_row else 0L,
		if (!is.null(R_col) && length(R_col) == 1L && is.finite(R_col)) R_col else 0L))
	if (rank_check_max > 0L) {
		n_use <- if (identical(mode, "bipartite")) {
			min(nrow(Y), ncol(Y))
		} else {
			nrow(Y)
		}
		n_third <- max(1L, floor(n_use / 3))
		if (length(n_use) == 1L && is.finite(n_use) && rank_check_max > n_third) {
			cli::cli_warn(c(
				"!" = "R = {.val {rank_check_max}} is large relative to n = {.val {n_use}} (recommended cap ~ n/3 = {.val {n_third}}).",
				"i" = "High-rank latent factors can absorb structure that belongs to the additive effects, hurting interpretability.",
				"i" = "Consider {.code R = {n_third}} or lower unless you have a specific reason."))
		}
	}
	# dynamic_beta is intentionally absent from ame(); use lame() for dynamic
	# coefficients

	# seed locally: restore the global rng stream on exit so a downstream
	# random draw is not silently perturbed by having fit a model
	if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
		.old_seed <- get(".Random.seed", envir = globalenv())
		on.exit(assign(".Random.seed", .old_seed, envir = globalenv()),
		        add = TRUE)
	}

	# validate family; accept short amen spellings as aliases
	valid_families <- c("normal", "binary", "ordinal",
	                    "cbin", "frn", "poisson")
	amen_aliases <- c(nrm = "normal", bin = "binary", ord = "ordinal",
	                  pois = "poisson")
	if (length(family) == 1L && is.character(family) &&
	    family %in% names(amen_aliases)) {
		new_family <- unname(amen_aliases[family])
		# announce the alias once per session; a porting script calls ame() in a
		# loop and the reminder is identical on every call
		if (!isTRUE(.lame_state$family_alias_msg_shown)) {
			cli::cli_inform(c(
				"i" = "{.arg family} = {.val {family}} (amen-style) accepted as alias for {.val {new_family}}.",
				"i" = "Update your script to {.code family = \"{new_family}\"} when convenient.",
				"i" = "This note is shown once per session."))
			.lame_state$family_alias_msg_shown <- TRUE
		}
		family <- new_family
	}
	if (length(family) != 1L || !is.character(family) ||
	    !family %in% valid_families) {
		cli::cli_abort(c(
			"{.arg family} = {.val {family}} is not a recognised family.",
			"i" = "Choose one of: {.val {valid_families}}."))
	}

	####
	# bipartite ordinal / cbin / frn / poisson share rectangular samplers;
	# rectangular z samplers in r/rz_bipartite.r back the ame_bipartite()
	# dispatch. no guard here.

	####
	# validate y and the key mcmc arguments up front so degenerate input
	# surfaces as a clear error rather than a deep-stack failure.
	####
	if (!is.matrix(Y) && !is.array(Y)) {
		cli::cli_abort(c(
			"{.arg Y} must be a matrix.",
			"i" = "For a numeric adjacency {.cls data.frame} (e.g. read in with {.fn read.csv}), convert it with {.code as.matrix(Y)}.",
			"i" = "For a long-format edgelist, build the network with {.fn netify::netify} and pass that object directly.",
			"i" = "{.fn as_lame_y} converts {.pkg igraph}/{.pkg network} objects to a lame-ready matrix."))
	}
	if (!is.numeric(Y) && !is.logical(Y)) {
		cli::cli_abort(c(
			"{.arg Y} must be numeric.",
			"i" = "Character input usually means an edgelist rather than an adjacency matrix; build the network with {.fn netify::netify} instead."))
	}
	if (any(is.infinite(Y) | is.nan(Y))) {
		cli::cli_abort(c(
			"{.arg Y} contains infinite or NaN values.",
			"i" = "Only finite values and {.val NA} (missing) are allowed."))
	}
	# symmetric=true must be accompanied by a symmetric y. otherwise the
	# model silently averages the upper and lower triangles and discards
	# half the observed asymmetry, no warning emitted.
	if (isTRUE(symmetric) && nrow(Y) == ncol(Y)) {
		off <- Y - t(Y)
		if (any(is.finite(off)) && max(abs(off), na.rm = TRUE) > 1e-8) {
			cli::cli_abort(c(
				"{.arg symmetric} = TRUE, but {.arg Y} is not symmetric.",
				"i" = "Symmetrize {.arg Y} (e.g. {.code (Y + t(Y))/2}) or set {.arg symmetric} = FALSE."))
		}
		# when one full triangle is na the symmetry check above passes
		# vacuously (y - t(y) is all-na off-diagonal) and the sampler would
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
		# warn when symmetric=true is paired with an asymmetric xdyad slice.
		# the linear predictor eta_ij = beta'x_ij + a_i + b_j + u_i'g v_j is
		# structurally symmetric on the (a, u, v) side (a == b averaged after
		# each sweep; ruv_sym_fc parameterises v = u l for l diagonal) but the
		# regression term inherits whatever symmetry the input x has. the
		# ordinal symmetric sampler (rz_ord_sym_fc) averages (ez[i,j] +
		# ez[j,i])/2 internally to compensate, so non-ordinal symmetric fits
		# carry the most risk of silent misspecification when xdyad is
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
	# a constant response cannot be fit under any family: continuous
	# families have zero residual variance, and a one-category discrete
	# response makes the latent covariance degenerate (chol failure deep
	# in the sampler)
	eff_y <- if (family %in% c("binary", "cbin")) 1 * (obs_y > 0) else obs_y
	if (length(unique(eff_y)) == 1L) {
		hint <- if (family %in% c("binary", "cbin")) {
			c("i" = "All observed cells threshold to {.val {unique(eff_y)}}.",
			  "i" = "If your edgelist lists only the ties that occurred, absent dyads should enter as 0s (in {.pkg netify}, revisit {.arg missing_to_zero}), not be missing.")
		} else {
			c("i" = "Every observed cell equals {.val {unique(eff_y)}}.")
		}
		cli::cli_abort(c(
			"{.arg Y} is constant; there is no relational variation to model.",
			hint))
	}
	if (length(R) != 1L || anyNA(R) || R < 0 ||
	    !isTRUE(all.equal(R, round(R)))) {
		cli::cli_abort("{.arg R} must be a single non-negative integer.")
	}
	if (R >= min(nrow(Y), ncol(Y))) {
		cli::cli_abort("{.arg R} = {R} must be smaller than the network dimension ({min(nrow(Y), ncol(Y))}).")
	}
	# bipartite latent dimensions get the same non-negative-integer / range
	# checks as r (otherwise r_row = 1.7 silently recycles a corrupt factor)
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
	if (family == "poisson" && any(obs_y %% 1 != 0)) {
		cli::cli_warn(c(
			"{.arg family} = {.val poisson} but {.arg Y} contains non-integer values.",
			"i" = "Values will be treated as counts; for continuous outcomes use {.val normal}."))
	}
	# binary and cbin threshold observed values at zero
	if (family %in% c("binary", "cbin")) {
		uy <- unique(obs_y)
		if (any(!uy %in% c(0, 1))) {
			cli::cli_warn(c(
				"{.arg family} = {.val {family}} but {.arg Y} contains values other than 0/1.",
				"i" = "{.arg Y} will be thresholded to {.code 1 * (Y > 0)}; if you meant counts, use {.val poisson}, or {.val ordinal}/{.val normal} as appropriate."))
		}
	}
	# warn on binary-looking y when family = "normal", but do not block
	if (family == "normal") {
		uy <- unique(obs_y)
		if (length(uy) <= 2L && all(uy %in% c(0, 1))) {
			cli::cli_inform(c(
				"i" = "{.arg Y} looks binary (values in {.val {{0, 1}}}) but {.arg family} = {.val normal}.",
				"i" = "If this is a 0/1 network, pass {.code family = \"binary\"}."))
		}
	}
	# rho is weakly identified in very sparse binary networks: replicate
	# fits of the same sparse DGP give wildly different rho estimates
	if (family %in% c("binary", "cbin") && isTRUE(dcor) &&
	    mode == "unipartite" && !symmetric) {
		dens <- mean(obs_y > 0)
		if (is.finite(dens) && dens < 0.02) {
			cli::cli_warn(c(
				"Network density is {round(100 * dens, 2)}%; the dyadic correlation {.code rho} is weakly identified below ~2% density.",
				"i" = "Treat the reported {.code rho} with caution, or fit with {.code dcor = FALSE}."))
		}
	}
	# validate Xdyad's actor dimensions up front, before any family-specific
	# scan indexes its slices against Y. the binary/cbin separation check
	# below compares Xdyad[, , k] with Y element-wise and would otherwise fail
	# with a bare "non-conformable arrays" error (rather than a clear message)
	# when a 3-D Xdyad's actor dimension does not match Y. the full 2-D/3-D
	# validation still runs later; this only hoists the actor-dimension guard.
	if (!is.null(Xdyad)) {
		dimd_early <- dim(Xdyad)
		if (!is.null(dimd_early) && length(dimd_early) == 3L &&
		    (dimd_early[1] != nrow(Y) || dimd_early[2] != ncol(Y))) {
			cli::cli_abort(c(
				"{.arg Xdyad} dimensions ({dimd_early[1]} x {dimd_early[2]} x {dimd_early[3]}) do not match {.arg Y} ({nrow(Y)} x {ncol(Y)}).",
				"i" = "The dyadic covariate array's first two dimensions must be the actor dimensions of {.arg Y}."))
		}
	}
	# complete separation by a single dyadic covariate: the coefficient
	# posterior random-walks upward and its mean depends on chain length,
	# exactly as in unpenalized logistic regression
	if (family %in% c("binary", "cbin") && !is.null(Xdyad) &&
	    length(dim(Xdyad)) == 3L) {
		yb <- 1 * (Y > 0)
		for (k in seq_len(dim(Xdyad)[3L])) {
			xk <- Xdyad[, , k]
			ok <- !is.na(yb) & !is.na(xk)
			x1 <- xk[ok & yb == 1]
			x0 <- xk[ok & yb == 0]
			if (length(x1) > 0L && length(x0) > 0L &&
			    (min(x1) > max(x0) || max(x1) < min(x0))) {
				nm_k <- dimnames(Xdyad)[[3L]][k] %||% as.character(k)
				cli::cli_warn(c(
					"{.arg Xdyad} covariate {.val {nm_k}} perfectly separates ties from non-ties.",
					"i" = "Its coefficient is only weakly identified through the prior: the chain drifts and the posterior mean depends on chain length.",
					"i" = "Consider a stronger prior (smaller {.arg g}) or removing the covariate."))
			}
		}
	}
	# an exactly symmetric Y under the directed default pushes the dyadic
	# correlation to the boundary; point at symmetric = TRUE instead of
	# letting the fit run degenerate
	if (!symmetric && mode == "unipartite" && nrow(Y) == ncol(Y) &&
	    isTRUE(all.equal(unname(Y), unname(t(Y))))) {
		cli::cli_warn(c(
			"{.arg Y} is exactly symmetric but {.arg symmetric} = FALSE.",
			"i" = "Fitting a directed model to undirected data drives the dyadic correlation to its boundary; pass {.code symmetric = TRUE} for the undirected model."))
	}
	# warn when burn dwarfs nscan -- typical use is burn << nscan; burn > nscan
	# wastes most of the chain
	if (burn > nscan) {
		cli::cli_warn(c(
			"{.arg burn} = {burn} is greater than {.arg nscan} = {nscan}.",
			"i" = "Most iterations will be discarded; consider {.arg burn} <= {.arg nscan}."))
	}
	# fully unobserved actors (all-na rows / cols) silently get apm/bpm purely
	# from the prior; surface this so an upstream data-pipeline bug is visible
	if (mode == "unipartite") {
		drow <- apply(Y, 1, function(r) all(is.na(r) | !is.finite(r)))
		dcol <- apply(Y, 2, function(c) all(is.na(c) | !is.finite(c)))
		# drop the diagonal na from the "all na" judgement
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
		# the bipartite-binary z sampler cannot proceed when a row
		# or column has no observed data; abort with a clear remediation
		# hint rather than letting the sampler crash deeper in the stack.
		if ((any(drow) || any(dcol)) && family %in% c("binary")) {
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
	# covariates are aligned to Y by actor name when both sides carry
	# names (matching lame()'s semantics); the realignment itself happens
	# below, after Xdyad has been normalised to a 3-D array.

	# xdyad should be a 3d array; wrap a 2d matrix with a short message
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
	# nodal covariates must have one row per actor; catch the off-by-one
	# early instead of failing inside design_array with a base error
	if (!is.null(Xrow) && !is.null(dim(Xrow)) && nrow(Xrow) != nrow(Y)) {
		cli::cli_abort(c(
			"{.arg Xrow} has {nrow(Xrow)} rows but {.arg Y} has {nrow(Y)} row actors.",
			"i" = "{.arg Xrow} needs exactly one row per row actor of {.arg Y}."))
	}
	if (!is.null(Xcol) && !is.null(dim(Xcol)) && nrow(Xcol) != ncol(Y)) {
		cli::cli_abort(c(
			"{.arg Xcol} has {nrow(Xcol)} rows but {.arg Y} has {ncol(Y)} column actors.",
			"i" = "{.arg Xcol} needs exactly one row per column actor of {.arg Y}."))
	}
	# align covariates to Y by actor name when both sides are named.
	# named-but-permuted inputs are reordered (with a message); named
	# inputs whose actors differ from Y's abort rather than silently
	# fitting a misaligned design. unnamed inputs stay positional.
	.align_idx <- function(have, want, what) {
		if (is.null(have) || is.null(want)) return(NULL)
		if (length(have) != length(want) || identical(have, want)) return(NULL)
		# default sequential rownames ("1", "2", ...) are not real names
		if (identical(have, as.character(seq_along(have)))) return(NULL)
		if (setequal(have, want)) return(match(want, have))
		cli::cli_abort(c(
			"{.arg {what}} carries actor names that do not match {.arg Y}'s.",
			"i" = "Covariates are aligned to {.arg Y} by name; rename {.arg {what}} to match, or drop its dimnames to force positional alignment."))
	}
	if (!is.null(Xdyad) && length(dim(Xdyad)) == 3L) {
		ir <- .align_idx(dimnames(Xdyad)[[1]], rownames(Y), "Xdyad")
		jc <- .align_idx(dimnames(Xdyad)[[2]], colnames(Y), "Xdyad")
		if (!is.null(ir) || !is.null(jc)) {
			Xdyad <- Xdyad[ir %||% seq_len(dim(Xdyad)[1L]),
			               jc %||% seq_len(dim(Xdyad)[2L]), , drop = FALSE]
			cli::cli_inform(c(
				"i" = "Reordered {.arg Xdyad} to match {.arg Y}'s actor order (aligned by dimnames)."))
		}
	}
	if (!is.null(Xrow) && !is.null(dim(Xrow))) {
		ir <- .align_idx(rownames(Xrow), rownames(Y), "Xrow")
		if (!is.null(ir)) {
			Xrow <- Xrow[ir, , drop = FALSE]
			cli::cli_inform(c(
				"i" = "Reordered {.arg Xrow} to match {.arg Y}'s actor order (aligned by rownames)."))
		}
	}
	if (!is.null(Xcol) && !is.null(dim(Xcol))) {
		ir <- .align_idx(rownames(Xcol), colnames(Y), "Xcol")
		if (!is.null(ir)) {
			Xcol <- Xcol[ir, , drop = FALSE]
			cli::cli_inform(c(
				"i" = "Reordered {.arg Xcol} to match {.arg Y}'s column actor order (aligned by rownames)."))
		}
	}
	# an entirely-NA covariate can never be identified; catch it here (the
	# bipartite path in particular has no missing-covariate augmentation,
	# so an all-NA slice would otherwise produce a silent random walk on
	# its coefficient)
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
				"i" = "An all-NA covariate cannot be identified; drop it or supply values."))
		}
		invisible()
	}
	.check_all_na_cov(Xdyad, "Xdyad")
	.check_all_na_cov(Xrow, "Xrow")
	.check_all_na_cov(Xcol, "Xcol")
	if (mode == "bipartite" &&
	    (anyNA(Xdyad) || anyNA(Xrow) || anyNA(Xcol))) {
		cli::cli_warn(c(
			"Missing covariate values are not data-augmented on the bipartite path.",
			"i" = "Coefficients touching NA design cells can be weakly identified; impute or drop the affected covariate values."))
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
		non_num <- intersect(names(prior)[!vapply(prior, is.numeric, logical(1))],
		                     known_prior)
		if (length(non_num) > 0L) {
			cli::cli_abort(c(
				"{.arg prior} element{?s} {.val {non_num}} must be numeric.",
				"i" = "Got class{?es}: {.cls {vapply(prior[non_num], function(x) class(x)[1], character(1))}}."))
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
		# symmetric ordinal uses the sampler in r/rz_ord_sym_fc.r
		# uses the symmetric-doubled precision (variance 1/2) and mirrors
		# upper-triangle draws to the lower triangle so z = t(z) holds at
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

	# post-hoc log_lik computation when save_log_lik = true. exact
	# observed-data family densities for normal/binary/cbin/poisson/
	# ordinal; augmented-z normal approximation only for frn. matches
	# the documented behaviour of loo.ame.
	if (want_log_lik) {
		fit$log_lik <- .ame_compute_log_lik_post_hoc(fit, Y = Y, family = family)
		fit$log_lik_method <- if (family %in% c("frn")) "augmented" else "observed_exact"
	}

	# surface a slow-mixing rho: short chains can leave too few effective
	# draws to trust its point estimate or interval
	if (verbose && !is.null(fit$VC) && is.matrix(fit$VC) &&
	    "rho" %in% colnames(fit$VC) && nrow(fit$VC) > 10L &&
	    requireNamespace("coda", quietly = TRUE)) {
		rho_ess <- tryCatch(
			as.numeric(coda::effectiveSize(coda::mcmc(fit$VC[, "rho"]))),
			error = function(e) NA_real_)
		if (is.finite(rho_ess) && rho_ess < 50) {
			cli::cli_inform(c(
				"i" = "Effective sample size for {.code rho} is only {round(rho_ess)}; its estimate and interval are unstable.",
				"i" = "Increase {.arg nscan} if inference on the dyadic correlation matters."))
		}
	}

	return(fit)
}

# post-hoc per-iteration log_lik for ame() fits. relies on
# posterior_opts$store_uv / store_ab having stored u_samples / v_samples /
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
	# observations: off-diagonal finite cells of y
	if (is.array(Y) && length(dim(Y)) == 2L) {
		if (nrow(Y) == ncol(Y)) diag(Y) <- NA
	}
	if (identical(family, "binary"))  Y <- 1 * (Y > 0)
	if (identical(family, "cbin"))    Y <- 1 * (Y > 0)
	obs_idx <- which(!is.na(Y) & is.finite(Y))
	n_obs <- length(obs_idx)
	if (n_obs == 0L) return(NULL)
	# pre-compute per-cell design components that don't change per iter
	n_a <- nrow(Y); n_b <- ncol(Y)
	# locate sigma^2 column on vc: last column is "ve" by convention.
	ve_col <- ncol(fit$VC)
	# beta_names: tells us which beta columns are dyadic vs intercept / nodal
	beta_names <- colnames(fit$BETA)
	has_intercept <- !is.null(beta_names) && beta_names[1] == "intercept"
	# pull out per-iteration u / v / a / b if available; else fall back to
	# posterior-mean u / v / apm / bpm (the latter understates per-draw spread
	# but at least lets loo() run on a fit where the samples were not stored)
	have_U <- !is.null(fit$U_samples) && !is.null(fit$V_samples)
	have_ab <- !is.null(fit$a_samples) && !is.null(fit$b_samples)
	# build a per-cell xdyad x beta matrix once. xdyad is stored on fit$x
	# when ame() saves the call args; rebuild from the call if missing.
	X <- fit$X
	ll <- matrix(NA_real_, nrow = n_iter, ncol = n_obs)
	# linear predictor for draw s
	compute_eta_draw <- function(s) {
		beta_s <- if (length(dim(fit$BETA)) == 3L) fit$BETA[s, , 1] else fit$BETA[s, ]
		eta_s <- matrix(0, n_a, n_b)
		if (is.null(X)) {
			if (has_intercept) eta_s <- eta_s + beta_s[1]
		} else {
			Xarr <- if (length(dim(X)) == 3L) X else array(X, c(n_a, n_b, 1L))
			p_x <- dim(Xarr)[3L]
			p_b <- length(beta_s)
			if (p_b == p_x) {
				# BETA columns correspond 1:1 to the design-array slices
				# (colnames(BETA) <- dimnames(X)[[3]]); slice 1 is the
				# all-ones intercept when an intercept is present
				for (k in seq_len(p_x)) {
					eta_s <- eta_s + beta_s[k] * Xarr[, , k]
				}
			} else {
				# symmetric fits store one coefficient per nodal covariate
				# while the design array keeps separate row/col slices:
				# slices [intercept | row (pr) | col (pr) | dyad (pd)],
				# betas  [intercept | node (pr)          | dyad (pd)]
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
				} else if (has_intercept) {
					eta_s <- eta_s + beta_s[1]
				}
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
	# per-iteration family-specific pointwise log-y density
	y_obs <- Y[obs_idx]
	# gauss-hermite rule for the poisson marginal (golub-welsch, base R;
	# nodes/weights for integrating against exp(-x^2))
	gh <- NULL
	if (family == "poisson") {
		k_gh <- 20L
		i_gh <- seq_len(k_gh - 1L)
		J_gh <- matrix(0, k_gh, k_gh)
		J_gh[cbind(i_gh, i_gh + 1L)] <- sqrt(i_gh / 2)
		J_gh[cbind(i_gh + 1L, i_gh)] <- sqrt(i_gh / 2)
		e_gh <- eigen(J_gh, symmetric = TRUE)
		gh <- list(nodes = e_gh$values, weights = e_gh$vectors[1, ]^2)
	}
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
		} else if (family == "poisson") {
			# marginal likelihood of the fitted overdispersed model,
			# p(y | eta, s2) = int Pois(y | e^z) N(z | eta, s2) dz,
			# by adaptive 20-node gauss-hermite: the rule is recentred per
			# observation at the laplace mode z* of the integrand, because
			# a rule centred at eta places all nodes within ~5.4*sqrt(s2)
			# of eta and understates poorly fit tail cells by tens of
			# nats. the plug-in dpois(y, exp(eta)) is the likelihood of a
			# different (equidispersed) model and grossly understates the
			# fitted model's log-likelihood.
			y_r <- round(y_obs)
			# newton for the mode of l(z) = y z - e^z - (z - eta)^2/(2 s2);
			# the init sits at/above the mode for y >= 1 at any eta, so
			# newton on the strictly concave l is monotonically convergent
			z_m <- log(y_r + pmax(eta_obs, 0) / ve_s + 1e-8)
			for (it_n in seq_len(50L)) {
				step_n <- (y_r - exp(z_m) - (z_m - eta_obs) / ve_s) /
					(-exp(z_m) - 1 / ve_s)
				z_m <- z_m - step_n
				if (max(abs(step_n)) < 1e-10) break
			}
			sd_star <- 1 / sqrt(exp(z_m) + 1 / ve_s)
			# log p = log( sqrt(2) sd* sum_k sqrt(pi) wbar_k e^{x_k^2}
			#              Pois(y | e^{z_k}) N(z_k | eta, s2) )
			acc <- matrix(-Inf, length(y_obs), k_gh)
			for (k_n in seq_len(k_gh)) {
				z_k <- z_m + sqrt(2) * sd_star * gh$nodes[k_n]
				acc[, k_n] <- stats::dpois(y_r, exp(z_k), log = TRUE) +
					stats::dnorm(z_k, eta_obs, sqrt(ve_s), log = TRUE) +
					gh$nodes[k_n]^2 + log(gh$weights[k_n]) + 0.5 * log(pi)
			}
			M_gh <- apply(acc, 1, max)
			# guard: if exp(z_k) overflowed at every node the row is all
			# -Inf and the log-sum-exp would return NaN; report -Inf
			ll[s, ] <- ifelse(is.finite(M_gh),
				M_gh + log(rowSums(exp(acc - M_gh))) +
					0.5 * log(2) + log(sd_star),
				-Inf)
		} else if (family == "ordinal") {
			ll[s, ] <- .ordinal_pointwise_loglik(y_obs, eta_obs)
		} else {
			# augmented-z normal fallback for frn
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
