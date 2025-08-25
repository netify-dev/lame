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
ame_bipartite <- function(Y, Xdyad=NULL, Xrow=NULL, Xcol=NULL, 
                        rvar = !(family=="rrl"), cvar = TRUE, R_row = 0, R_col = 0,
                        family="normal", intercept=!is.element(family,c("rrl","ordinal")),
                        odmax=NULL, prior=list(), g=NA,
                        seed = 6886, nscan = 10000, burn = 500, odens = 25,
                        plot=FALSE, print = TRUE, gof=TRUE, 
                        start_vals=NULL, periodic_save=FALSE, out_file=NULL,
                        save_interval=0.25, model.name=NULL) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Get dimensions
  nA <- nrow(Y)  # Number of row nodes
  nB <- ncol(Y)  # Number of column nodes
  
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
  
  # Storage for posterior samples
  BETA <- matrix(NA, n_eff, p)
  VC <- matrix(NA, n_eff, 5)  # Store variance components
  if(!is.null(U) && !is.null(V) && !is.null(G)) {
    UVPS <- U %*% G %*% t(V) * 0  # Initialize UV product sum
  } else {
    UVPS <- matrix(0, nA, nB)
  }
  APS <- rep(0, nA)  # Additive row effects posterior sum
  BPS <- rep(0, nB)  # Additive column effects posterior sum
  YPS <- array(0, dim = c(nA, nB))  # Posterior predictive sum
  EZ <- array(0, dim = c(nA, nB))  # Expected value sum
  
  if(gof) {
    GOF <- matrix(NA, n_eff, 4)  # Store GOF statistics
  } else {
    GOF <- NULL
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
    cat("Running MCMC for bipartite network...\n")
    cat("Dimensions: ", nA, "x", nB, "\n")
    cat("R_row =", R_row, ", R_col =", R_col, "\n")
  }
  
  iter_save <- 0
  
  for(iter in 1:nscan) {
    # Compute current expected value
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
    
    # Update Z given other parameters based on family
    if(family == "normal") {
      # For normal, Z = Y (observed)
      Z <- Y
    } else if(family == "tobit") {
      # Truncated normal for censored data
      Z <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      Z[Y == 0] <- pmin(Z[Y == 0], 0)  # Censored at 0
      Z[Y > 0] <- Y[Y > 0]  # Observed values
    } else if(family == "binary") {
      # Truncated normal
      Z <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      Z[Y == 0] <- pmin(Z[Y == 0], 0)
      Z[Y == 1] <- pmax(Z[Y == 1], 0)
    } else if(family == "ordinal") {
      # Ordinal probit - simplified version
      # Would need threshold parameters in full implementation
      Z <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
    } else if(family == "cbin") {
      # Censored binary - similar to binary but with row constraints
      Z <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      Z[Y == 0] <- pmin(Z[Y == 0], 0)
      Z[Y == 1] <- pmax(Z[Y == 1], 0)
    } else if(family == "frn") {
      # Fixed rank nomination - simplified
      Z <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
    } else if(family == "rrl") {
      # Row rank likelihood - simplified
      Z <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
    } else if(family == "poisson") {
      # Use auxiliary variable approach
      for(i in 1:nA) {
        for(j in 1:nB) {
          if(!is.na(Y[i,j])) {
            lambda <- exp(EZ_curr[i,j])
            # Use data augmentation
            Z[i,j] <- log(Y[i,j] + rgamma(1, 1, lambda))
          }
        }
      }
    }
    
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
        resid <- resid - U %*% G %*% t(V)
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
        resid <- resid - U %*% G %*% t(V)
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
        resid <- resid - U %*% G %*% t(V)
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
    
    # Update multiplicative effects (simplified)
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
        resid <- resid - U %*% G %*% t(V)
      }
      
      shape <- (nA * nB) / 2
      rate <- sum(resid^2, na.rm = TRUE) / 2
      s2 <- 1 / rgamma(1, shape, rate)
    }
    
    # Store samples after burn-in
    if(iter > burn && (iter - burn) %% odens == 0) {
      iter_save <- iter_save + 1
      
      BETA[iter_save,] <- beta
      VC[iter_save,] <- c(Sab[1,1], Sab[1,2], Sab[2,2], s2, 0)  # No rho for bipartite
      
      if(rvar) APS <- APS + a
      if(cvar) BPS <- BPS + b
      if(!is.null(U) && !is.null(V) && !is.null(G)) {
        UVPS <- UVPS + U %*% G %*% t(V)
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
      
      EZ <- EZ + EZ_curr
      
      # Generate posterior predictive sample
      if(family == "normal") {
        Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
      } else if(family == "tobit") {
        Y_sim <- EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)
        Y_sim[Y_sim < 0] <- 0  # Censor at 0
      } else if(family == "binary") {
        Y_sim <- (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)) > 0
      } else if(family == "ordinal") {
        # Simplified ordinal simulation
        Y_sim <- round(EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB))
        Y_sim[Y_sim < 0] <- 0
        if(any(!is.na(Y))) {
          max_cat <- max(Y, na.rm = TRUE)
          Y_sim[Y_sim > max_cat] <- max_cat
        }
      } else if(family == "cbin") {
        Y_sim <- (EZ_curr + matrix(rnorm(nA * nB, 0, sqrt(s2)), nA, nB)) > 0
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
      
      YPS <- YPS + Y_sim
      
      # Compute GOF statistics
      if(gof) {
        # Simple GOF: density, mean, variance, correlation
        GOF[iter_save, 1] <- mean(Y > 0, na.rm = TRUE)  # Observed density
        GOF[iter_save, 2] <- mean(Y_sim > 0, na.rm = TRUE)  # Simulated density
        GOF[iter_save, 3] <- mean(Y, na.rm = TRUE)  # Observed mean
        GOF[iter_save, 4] <- mean(Y_sim, na.rm = TRUE)  # Simulated mean
      }
    }
    
    # Progress
    if(print && iter %% 1000 == 0) {
      cat("Iteration", iter, "of", nscan, "\n")
    }
  }
  
  # Compute posterior means
  APM <- APS / n_eff
  BPM <- BPS / n_eff
  UVPM <- UVPS / n_eff
  EZ <- EZ / n_eff
  YPM <- YPS / n_eff
  
  # Transform for count models
  if(family == "poisson") {
    EZ <- exp(EZ)
    EZ[EZ > 1e6] <- 1e6
  }
  
  # Prepare output
  fit <- list(
    BETA = BETA,
    VC = VC,
    APM = APM,
    BPM = BPM,
    U = U,
    V = V,
    G = G,
    UVPM = UVPM,
    EZ = EZ,
    YPM = YPM,
    GOF = GOF,
    X = X,
    Y = Y,
    family = family,
    R_row = R_row,
    R_col = R_col,
    rvar = rvar,
    cvar = cvar,
    mode = "bipartite",
    nA = nA,
    nB = nB
  )
  
  class(fit) <- "ame"
  
  return(fit)
}