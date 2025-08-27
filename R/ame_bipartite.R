#' AME model fitting for bipartite (rectangular) networks
#' 
#' Internal function that handles AME model fitting for bipartite networks.
#' This function contains the core logic for rectangular adjacency matrices
#' with distinct row and column node sets.
#' 
#' @inheritParams ame
#' @param R_row integer: dimension of row node multiplicative effects
#' @param R_col integer: dimension of column node multiplicative effects
#' @return An ame object with posterior samples and estimates
#' @keywords internal
#' @noRd
ame_bipartite <- function(
  Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
  rvar = !(family=="rrl"), cvar = TRUE, R_row = 0, R_col = 0,
  family="normal", intercept=!is.element(family,c("rrl","ordinal")),
  odmax=NULL, prior=list(), g=NA,
  seed = 6886, nscan = 10000, burn = 500, odens = 25,
  print = TRUE, gof=TRUE, custom_gof=NULL,
  start_vals=NULL, periodic_save=FALSE, out_file=NULL,
  save_interval=0.25, model.name=NULL,
  posterior_opts = NULL, use_sparse_matrices = FALSE
  ){
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Get dimensions
  nA <- nrow(Y)  # Number of row nodes
  nB <- ncol(Y)  # Number of column nodes
  
  # Set default node names if not provided
  if(is.null(rownames(Y))) {
    rownames(Y) <- paste0("rowNode", 1:nA)
  }
  if(is.null(colnames(Y))) {
    colnames(Y) <- paste0("colNode", 1:nB)
  }
  
  # Note: bipartite networks don't have diagonals to set to NA
  
  # Initialize storage for MCMC samples
  n_eff <- floor((nscan - burn) / odens)
  
  # Process prior parameters
  Sab0 <- prior$Sab0
  if(is.null(Sab0)) Sab0 <- diag(2)
  eta0 <- prior$eta0
  if(is.null(eta0)) eta0 <- 4 + 3 * (nA + nB) / 200
  etaab <- prior$etaab
  if(is.null(etaab)) etaab <- eta0
  s20 <- prior$s20
  if(is.null(s20)) s20 <- 1
  s2u0 <- prior$s2u0
  if(is.null(s2u0)) s2u0 <- 1
  Suv0 <- prior$Suv0
  if(is.null(Suv0)) Suv0 <- diag(max(R_row, R_col))
  
  # Process covariates
  if(!is.null(Xdyad)) {
    Xdyad <- as.array(Xdyad)
    if(length(dim(Xdyad)) == 2) {
      Xdyad <- array(Xdyad, c(nA, nB, 1))
    }
    p_dyad <- dim(Xdyad)[3]
  } else {
    p_dyad <- 0
  }
  
  if(!is.null(Xrow)) {
    Xrow <- as.matrix(Xrow)
    p_row <- ncol(Xrow)
  } else {
    p_row <- 0
  }
  
  if(!is.null(Xcol)) {
    Xcol <- as.matrix(Xcol)
    p_col <- ncol(Xcol)
  } else {
    p_col <- 0
  }
  
  # Total number of regression coefficients
  p <- p_dyad + p_row + p_col + as.numeric(intercept)
  
  # Initialize parameters
  if(is.null(start_vals)) {
    # Initialize Z
    if(family == "normal") {
      Z <- Y
    } else if(family == "tobit") {
      Z <- Y
      Z[Y < 0] <- 0  # Tobit censoring at 0
      if(sum(Y > 0, na.rm = TRUE) > 0) {
        Z[Y == 0] <- -mean(Y[Y > 0], na.rm = TRUE)
      }
    } else if(family == "binary") {
      Z <- matrix(rnorm(nA * nB), nA, nB)
      Z[Y == 1] <- abs(Z[Y == 1])
      Z[Y == 0] <- -abs(Z[Y == 0])
    } else if(family == "ordinal") {
      # Convert to z-scores
      Z <- matrix(qnorm((Y + 0.5) / (max(Y, na.rm = TRUE) + 1)), nA, nB)
    } else if(family == "cbin") {
      Z <- matrix(rnorm(nA * nB), nA, nB)
      Z[Y == 1] <- abs(Z[Y == 1])
      Z[Y == 0] <- -abs(Z[Y == 0])
    } else if(family == "frn") {
      # Fixed rank nomination
      Z <- matrix(rnorm(nA * nB), nA, nB)
      for(i in 1:nA) {
        if(any(Y[i,] > 0, na.rm = TRUE)) {
          ranked <- which(Y[i,] > 0)
          Z[i, ranked] <- sort(abs(rnorm(length(ranked))), decreasing = TRUE)
        }
      }
    } else if(family == "rrl") {
      # Row rank likelihood - use row-wise z-scores
      Z <- matrix(0, nA, nB)
      for(i in 1:nA) {
        row_vals <- Y[i,]
        if(any(!is.na(row_vals))) {
          Z[i,] <- qnorm((rank(row_vals, na.last = "keep") - 0.5) / sum(!is.na(row_vals)))
        }
      }
    } else if(family == "poisson") {
      Z <- log(Y + 1)
    } else {
      Z <- Y  # Default
    }
    
    # Initialize parameters
    beta <- rep(0, p)
    a <- rep(0, nA) 
    b <- rep(0, nB)
    Sab <- Sab0
    s2 <- 1
    
    # Initialize multiplicative effects
    if(R_row > 0) {
      U <- matrix(rnorm(nA * R_row, 0, 0.1), nA, R_row)
    } else {
      U <- NULL
    }
    
    if(R_col > 0) {
      V <- matrix(rnorm(nB * R_col, 0, 0.1), nB, R_col)
    } else {
      V <- NULL
    }
    
    # Initialize interaction matrix G
    if(R_row > 0 && R_col > 0) {
      G <- diag(min(R_row, R_col), R_row, R_col)
    } else {
      G <- NULL
    }
  } else {
    # Use provided starting values
    Z <- start_vals$Z
    beta <- start_vals$beta
    a <- start_vals$a
    b <- start_vals$b
    Sab <- start_vals$Sab
    s2 <- start_vals$s2
    U <- start_vals$U
    V <- start_vals$V
    G <- start_vals$G
  }
  
  # Storage for posterior samples - pre-allocated for efficiency
  BETA <- matrix(NA_real_, n_eff, p)
  VC <- matrix(NA_real_, n_eff, 5)  # Store variance components
  sample_idx <- 0  # Track current sample index
  
  # Initialize posterior sample storage if requested
  U_samples <- NULL
  V_samples <- NULL
  a_samples <- NULL
  b_samples <- NULL
  
  if(!is.null(posterior_opts)) {
    if(posterior_opts$save_UV && R_row > 0 && R_col > 0) {
      # Calculate thinned sample size
      n_UV_samples <- ceiling(n_eff / posterior_opts$thin_UV)
      U_samples <- array(NA_real_, dim = c(nA, R_row, n_UV_samples))
      V_samples <- array(NA_real_, dim = c(nB, R_col, n_UV_samples))
    }
    if(posterior_opts$save_ab && (rvar || cvar)) {
      # Calculate thinned sample size
      n_ab_samples <- ceiling(n_eff / posterior_opts$thin_ab)
      if(rvar) a_samples <- matrix(NA_real_, nA, n_ab_samples)
      if(cvar) b_samples <- matrix(NA_real_, nB, n_ab_samples)
    }
  }
  
  # Track sample indices for thinned storage
  UV_sample_idx <- 0
  ab_sample_idx <- 0
  APS <- rep(0, nA)  # Additive row effects posterior sum
  BPS <- rep(0, nB)  # Additive column effects posterior sum
  YPS <- array(0, dim = c(nA, nB))  # Posterior predictive sum
  # Note: EZ not accumulated - can be reconstructed from components if needed
  
  if(gof) {
    # Determine number of GOF statistics
    n_base_stats <- 3  # sd.rowmean, sd.colmean, four.cycles
    
    # Process custom GOF functions
    custom_gof_names <- NULL
    n_custom_stats <- 0
    if(!is.null(custom_gof)) {
      if(is.function(custom_gof)) {
        # Single function - get example output to determine size
        tryCatch({
          test_output <- custom_gof(Y)
          n_custom_stats <- length(test_output)
          custom_gof_names <- names(test_output)
          if(is.null(custom_gof_names)) {
            custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
          }
        }, error = function(e) {
          warning("Custom GOF function failed on initial test. It will be skipped.")
          custom_gof <- NULL
        })
      } else if(is.list(custom_gof)) {
        # List of functions
        n_custom_stats <- length(custom_gof)
        custom_gof_names <- names(custom_gof)
        if(is.null(custom_gof_names)) {
          custom_gof_names <- paste0("custom", seq_len(n_custom_stats))
        }
      }
    }
    
    # Initialize GOF matrix
    n_total_stats <- n_base_stats + n_custom_stats
    GOF <- matrix(NA, n_eff, n_total_stats)
    
    # Set column names
    base_names <- c("sd.rowmean", "sd.colmean", "four.cycles")
    all_names <- c(base_names, custom_gof_names)
    colnames(GOF) <- all_names
  } else {
    GOF <- NULL
    custom_gof <- NULL  # Don't compute custom GOF if gof=FALSE
  }
  
  # Construct design matrix for regression
  if(p > 0) {
    X <- array(0, c(nA, nB, p))
    idx <- 1
    
    if(intercept) {
      X[,,idx] <- 1
      idx <- idx + 1
    }
    
    if(p_dyad > 0) {
      for(k in 1:p_dyad) {
        X[,,idx] <- Xdyad[,,k]
        idx <- idx + 1
      }
    }
    
    if(p_row > 0) {
      for(k in 1:p_row) {
        X[,,idx] <- Xrow[,k]
        idx <- idx + 1
      }
    }
    
    if(p_col > 0) {
      for(k in 1:p_col) {
        for(j in 1:nB) {
          X[,j,idx] <- Xcol[j,k]
        }
        idx <- idx + 1
      }
    }
  } else {
    X <- NULL
  }
  
  # MCMC iterations
  if(print) {
    cli::cli_h2("Running MCMC for bipartite network")
    cli::cli_text("Dimensions: {.val {nA}} x {.val {nB}}")
    cli::cli_text("R_row = {.val {R_row}}, R_col = {.val {R_col}}")
    cli::cli_text("Settings: Burn-in = {.val {burn}}, MCMC = {.val {nscan}}, Thin = {.val {odens}}")
    
    # Estimate time impact of GOF
    if(gof) {
      n_gof_stats <- 3
      if(!is.null(custom_gof)) {
        if(is.function(custom_gof)) {
          n_gof_stats <- n_gof_stats + length(custom_gof(Y))
        } else if(is.list(custom_gof)) {
          n_gof_stats <- n_gof_stats + length(custom_gof)
        }
      }
    }
    
    # Start timing
    start_time <- Sys.time()
    
    if(burn > 0) {
      cli::cli_progress_bar("Burn-in", total = burn, clear = TRUE)
    }
  }
  
  # Pre-select family-specific functions to avoid string comparisons
  update_Z_fn <- switch(family,
    normal = function(Z, EZ, s2, Y) Y,
    tobit = function(Z, EZ, s2, Y) {
      Z_new <- EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      Z_new[Y == 0] <- pmin(Z_new[Y == 0], 0)
      Z_new[Y > 0] <- Y[Y > 0]
      Z_new
    },
    binary = function(Z, EZ, s2, Y) {
      Z_new <- EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      Z_new[Y == 0] <- pmin(Z_new[Y == 0], 0)
      Z_new[Y == 1] <- pmax(Z_new[Y == 1], 0)
      Z_new
    },
    poisson = function(Z, EZ, s2, Y) {
      lambda <- exp(pmin(EZ, 10))
      # Vectorized update for non-missing values
      not_na <- !is.na(Y)
      n_update <- sum(not_na)
      if(n_update > 0) {
        Z_new <- Z
        Z_new[not_na] <- log(Y[not_na] + rgamma(n_update, 1, lambda[not_na]))
        Z_new
      } else {
        Z
      }
    },
    ordinal = function(Z, EZ, s2, Y) {
      EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
    },
    cbin = function(Z, EZ, s2, Y) {
      Z_new <- EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      Z_new[Y == 0] <- pmin(Z_new[Y == 0], 0)
      Z_new[Y == 1] <- pmax(Z_new[Y == 1], 0)
      Z_new
    },
    frn = function(Z, EZ, s2, Y) {
      EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
    },
    rrl = function(Z, EZ, s2, Y) {
      EZ + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
    },
    stop(paste("Unknown family:", family))
  )
  
  needs_s2 <- family %in% c("normal", "tobit", "poisson")
  
  # Initialize cache components for EZ computation
  Xbeta_cache <- matrix(0, nA, nB)
  if(p > 0) {
    for(k in 1:p) Xbeta_cache <- Xbeta_cache + beta[k] * X[,,k]
  }
  
  a_mat_cache <- if(rvar) matrix(a, nA, nB, byrow = FALSE) else matrix(0, nA, nB)
  b_mat_cache <- if(cvar) matrix(b, nA, nB, byrow = TRUE) else matrix(0, nA, nB)
  UV_cache <- if(!is.null(U) && !is.null(V) && !is.null(G)) {
    U %*% G %*% t(V)
  } else {
    matrix(0, nA, nB)
  }
  
  EZ_cache <- Xbeta_cache + a_mat_cache + b_mat_cache + UV_cache
  
  # Pre-compute X statistics once if needed
  if(p > 0) {
    XtX <- matrix(0, p, p)
    for(k1 in 1:p) {
      for(k2 in k1:p) {
        XtX[k1,k2] <- sum(X[,,k1] * X[,,k2], na.rm = TRUE)
        if(k1 != k2) XtX[k2,k1] <- XtX[k1,k2]
      }
    }
  }
  
  iter_save <- 0
  
  for(iter in 1:(nscan + burn)) {
    # Use cached EZ
    EZ_curr <- EZ_cache
    
    # Update Z using pre-selected function
    Z <- update_Z_fn(Z, EZ_curr, s2, Y)
    
    # Update regression coefficients
    if(p > 0) {
      # Compute residuals
      resid <- Z
      if(rvar) resid <- resid - a
      if(cvar) {
        for(j in 1:nB) {
          resid[,j] <- resid[,j] - b[j]
        }
      }
      if(!is.null(U) && !is.null(V) && !is.null(G)) {
        resid <- resid - UV_cache
      }
      
      # Use g-prior update
      XtX <- matrix(0, p, p)
      Xty <- rep(0, p)
      
      for(k1 in 1:p) {
        for(k2 in 1:p) {
          XtX[k1,k2] <- sum(X[,,k1] * X[,,k2], na.rm = TRUE)
        }
        Xty[k1] <- sum(X[,,k1] * resid, na.rm = TRUE)
      }
      
      if(is.na(g)) g <- nA * nB
      
      V_beta <- solve(XtX / s2 + diag(p) / (g * s20))
      m_beta <- V_beta %*% (Xty / s2)
      
      beta <- m_beta + t(chol(V_beta)) %*% rnorm(p)
    }
    
    # Update additive effects
    if(rvar) {
      # Update row effects
      resid <- Z
      if(p > 0) {
        for(k in 1:p) {
          resid <- resid - beta[k] * X[,,k]
        }
      }
      if(cvar) {
        for(j in 1:nB) {
          resid[,j] <- resid[,j] - b[j]
        }
      }
      if(!is.null(U) && !is.null(V) && !is.null(G)) {
        resid <- resid - UV_cache
      }
      
      for(i in 1:nA) {
        prec_a <- nB / s2 + 1 / Sab[1,1]
        mean_a <- sum(resid[i,], na.rm = TRUE) / s2 / prec_a
        a[i] <- rnorm(1, mean_a, sqrt(1/prec_a))
      }
    }
    
    if(cvar) {
      # Update column effects
      resid <- Z
      if(p > 0) {
        for(k in 1:p) {
          resid <- resid - beta[k] * X[,,k]
        }
      }
      if(rvar) resid <- resid - a
      if(!is.null(U) && !is.null(V) && !is.null(G)) {
        resid <- resid - UV_cache
      }
      
      for(j in 1:nB) {
        prec_b <- nA / s2 + 1 / Sab[2,2]
        mean_b <- sum(resid[,j], na.rm = TRUE) / s2 / prec_b
        b[j] <- rnorm(1, mean_b, sqrt(1/prec_b))
      }
    }
    
    # Update covariance of additive effects
    if(rvar || cvar) {
      SS <- matrix(0, 2, 2)
      if(rvar) SS[1,1] <- sum(a^2)
      if(cvar) SS[2,2] <- sum(b^2)
      
      # Sample from inverse Wishart
      nu_post <- etaab + max(nA * rvar, nB * cvar)
      S_post <- etaab * Sab0 + SS
      
      # Simple approximation for 2x2 inverse Wishart
      Sab <- solve(matrix(rWishart(1, nu_post, solve(S_post)), 2, 2))
    }
    
    # Update multiplicative effects
    if(R_row > 0 && R_col > 0) {
      # Update U
      resid <- Z
      if(p > 0) {
        for(k in 1:p) {
          resid <- resid - beta[k] * X[,,k]
        }
      }
      if(rvar) resid <- resid - a
      if(cvar) {
        for(j in 1:nB) {
          resid[,j] <- resid[,j] - b[j]
        }
      }
      
      # Simple random walk update for U
      for(i in 1:nA) {
        U_prop <- U
        U_prop[i,] <- U[i,] + rnorm(R_row, 0, 0.1)
        
        # Compute log ratio
        ll_curr <- -sum((resid[i,] - U[i,] %*% G %*% t(V))^2, na.rm = TRUE) / (2*s2)
        ll_prop <- -sum((resid[i,] - U_prop[i,] %*% G %*% t(V))^2, na.rm = TRUE) / (2*s2)
        
        if(log(runif(1)) < ll_prop - ll_curr) {
          U <- U_prop
        }
      }
      
      # Update V similarly
      for(j in 1:nB) {
        V_prop <- V
        V_prop[j,] <- V[j,] + rnorm(R_col, 0, 0.1)
        
        ll_curr <- -sum((resid[,j] - U %*% G %*% V[j,])^2, na.rm = TRUE) / (2*s2)
        ll_prop <- -sum((resid[,j] - U %*% G %*% V_prop[j,])^2, na.rm = TRUE) / (2*s2)
        
        if(log(runif(1)) < ll_prop - ll_curr) {
          V <- V_prop
        }
      }
      
      # Update UV cache after modifications
      UV_cache <- U %*% G %*% t(V)
    }
    
    # Update error variance
    if(family == "normal") {
      resid <- Z
      if(p > 0) {
        for(k in 1:p) {
          resid <- resid - beta[k] * X[,,k]
        }
      }
      if(rvar) resid <- resid - a
      if(cvar) {
        for(j in 1:nB) {
          resid[,j] <- resid[,j] - b[j]
        }
      }
      if(!is.null(U) && !is.null(V) && !is.null(G)) {
        resid <- resid - UV_cache
      }
      
      shape <- (nA * nB) / 2
      rate <- sum(resid^2, na.rm = TRUE) / 2
      s2 <- 1 / rgamma(1, shape, rate)
    }
    
    # Update EZ cache for next iteration
    EZ_cache <- Xbeta_cache + a_mat_cache + b_mat_cache + UV_cache
    
    # Store samples after burn-in
    if(iter > burn && (iter - burn) %% odens == 0) {
      iter_save <- iter_save + 1
      
      if(iter_save <= n_eff) {
        BETA[iter_save,] <- beta
        VC[iter_save,] <- c(Sab[1,1], Sab[1,2], Sab[2,2], s2, 0)  # No rho for bipartite
        
        if(rvar) APS <- APS + a
        if(cvar) BPS <- BPS + b
        if(!is.null(U) && !is.null(V) && !is.null(G)) {
          # UVPS not accumulated
        }
        
        # Save posterior samples if requested
        if(!is.null(posterior_opts)) {
          # Save U/V samples (thinned)
          if(posterior_opts$save_UV && !is.null(U) && !is.null(V)) {
            if(iter_save %% posterior_opts$thin_UV == 0) {
              UV_sample_idx <- UV_sample_idx + 1
              if(UV_sample_idx <= dim(U_samples)[3]) {
                U_samples[,,UV_sample_idx] <- U
                V_samples[,,UV_sample_idx] <- V
              }
            }
          }
          
          # Save a/b samples (thinned)
          if(posterior_opts$save_ab) {
            if(iter_save %% posterior_opts$thin_ab == 0) {
              ab_sample_idx <- ab_sample_idx + 1
              if(rvar && !is.null(a_samples) && ab_sample_idx <= ncol(a_samples)) {
                a_samples[,ab_sample_idx] <- a
              }
              if(cvar && !is.null(b_samples) && ab_sample_idx <= ncol(b_samples)) {
                b_samples[,ab_sample_idx] <- b
              }
            }
          }
        }
      }
      
      # Compute expected value
      EZ_curr <- matrix(0, nA, nB)
      if(p > 0) {
        for(k in 1:p) {
          EZ_curr <- EZ_curr + beta[k] * X[,,k]
        }
      }
      if(rvar) EZ_curr <- EZ_curr + a
      if(cvar) {
        for(j in 1:nB) {
          EZ_curr[,j] <- EZ_curr[,j] + b[j]
        }
      }
      if(!is.null(U) && !is.null(V) && !is.null(G)) {
        EZ_curr <- EZ_curr + U %*% G %*% t(V)
      }
      
      # Generate posterior predictive sample
      if(family == "normal") {
        Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      } else if(family == "tobit") {
        Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
        Y_sim[Y_sim < 0] <- 0  # Censor at 0
      } else if(family == "binary") {
        Y_sim <- 1 * (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB) > 0)
      } else if(family == "ordinal") {
        # Simplified ordinal simulation
        Y_sim <- round(EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB))
        Y_sim[Y_sim < 0] <- 0
        if(any(!is.na(Y))) {
          max_cat <- max(Y, na.rm = TRUE)
          Y_sim[Y_sim > max_cat] <- max_cat
        }
      } else if(family == "cbin") {
        Y_sim <- 1 * (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB) > 0)
        # Apply row constraints if odmax specified
        if(!is.null(odmax)) {
          for(i in 1:nA) {
            if(sum(Y_sim[i,]) > odmax[i]) {
              # Keep only top odmax[i] values
              thresh <- sort(EZ_curr[i,], decreasing = TRUE)[odmax[i]]
              Y_sim[i, EZ_curr[i,] < thresh] <- 0
            }
          }
        }
      } else if(family == "frn") {
        # Fixed rank nomination
        Y_sim <- matrix(0, nA, nB)
        if(!is.null(odmax)) {
          for(i in 1:nA) {
            # Rank the values and assign ranks to top odmax[i]
            ranks <- rank(-EZ_curr[i,], ties.method = "random")
            Y_sim[i, ranks <= odmax[i]] <- ranks[ranks <= odmax[i]]
          }
        }
      } else if(family == "rrl") {
        # Row rank likelihood
        Y_sim <- matrix(0, nA, nB)
        for(i in 1:nA) {
          Y_sim[i,] <- rank(EZ_curr[i,] + rnorm(nB, 0, sqrt(s2)))
        }
      } else if(family == "poisson") {
        Y_sim <- matrix(rpois(nA * nB, exp(EZ_curr)), nA, nB)
      } else {
        Y_sim <- EZ_curr
      }
      
      if(iter_save <= n_eff) {
        YPS <- YPS + Y_sim
      }
      
      # Compute GOF statistics
      if(gof && iter_save <= n_eff) {
        # Use bipartite-specific GOF function for base stats
        gof_vals <- gof_stats_bipartite(Y)
        GOF[iter_save, 1:3] <- gof_vals
        
        # Compute custom GOF statistics if provided
        if(!is.null(custom_gof)) {
          if(is.function(custom_gof)) {
            # Single function returning possibly multiple values
            tryCatch({
              custom_vals <- custom_gof(Y)
              GOF[iter_save, (4):(3 + length(custom_vals))] <- custom_vals
            }, error = function(e) {
              # Leave as NA if function fails
            })
          } else if(is.list(custom_gof)) {
            # List of functions each returning single value
            for(i in seq_along(custom_gof)) {
              if(is.function(custom_gof[[i]])) {
                tryCatch({
                  GOF[iter_save, 3 + i] <- custom_gof[[i]](Y)
                }, error = function(e) {
                  # Leave as NA if function fails
                })
              }
            }
          }
        }
      }
    }
    
    # Update progress bars
    if(print) {
      if(iter <= burn) {
        cli::cli_progress_update(id = NULL)
        if(iter == burn) {
          cli::cli_progress_done()
          cli::cli_progress_bar("Sampling", total = nscan, clear = TRUE)
        }
      } else if(iter > burn) {
        cli::cli_progress_update(id = NULL)
      }
    }
  }
  
  # Close progress bar if printing
  if(print) {
    cli::cli_progress_done()
    
    # Report timing
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    cli::cli_alert_success("MCMC complete in {.val {round(elapsed, 1)}} seconds")
  }
  
  # Compute posterior means
  APM <- APS / n_eff
  BPM <- BPS / n_eff
  YPM <- YPS / n_eff
  
  # Prepare output
  fit <- list(
    BETA = BETA,
    VC = VC,
    APM = APM,
    BPM = BPM,
    U = U,
    V = V,
    G = G,
    YPM = YPM,
    GOF = GOF,
    X = X,
    Y = Y,
    family = family,
    R_row = R_row,
    R_col = R_col,
    rvar = rvar,
    cvar = cvar,
    symmetric = FALSE,  # Bipartite networks are never symmetric
    mode = "bipartite",
    nA = nA,
    nB = nB
  )
  
  # Add posterior samples if saved
  if(!is.null(posterior_opts)) {
    if(!is.null(U_samples)) {
      # Trim to actual number of samples saved
      actual_UV_samples <- min(UV_sample_idx, dim(U_samples)[3])
      if(actual_UV_samples > 0) {
        fit$U_samples <- U_samples[,,1:actual_UV_samples, drop=FALSE]
        fit$V_samples <- V_samples[,,1:actual_UV_samples, drop=FALSE]
      }
    }
    if(!is.null(a_samples)) {
      actual_ab_samples <- min(ab_sample_idx, ncol(a_samples))
      if(actual_ab_samples > 0 && rvar) {
        fit$a_samples <- a_samples[,1:actual_ab_samples, drop=FALSE]
      }
    }
    if(!is.null(b_samples)) {
      actual_ab_samples <- min(ab_sample_idx, ncol(b_samples))
      if(actual_ab_samples > 0 && cvar) {
        fit$b_samples <- b_samples[,1:actual_ab_samples, drop=FALSE]
      }
    }
  }
  
  class(fit) <- "ame"
  
  # Apply sparse matrix optimization if requested  
  if(use_sparse_matrices) {
    fit <- compact_ame(fit, use_sparse_matrices = use_sparse_matrices)
  }
  
  return(fit)
}