#' Bipartite network helper functions
#' 
#' @name bipartite_helpers
#' @keywords internal

#' Initialize bipartite starting values
#' @keywords internal
init_bipartite_startvals <- function(Y, family, nA, nB, RA, RB, T, 
                                    Xlist = NULL, odmax = NULL) {
  # Initialize Z
  Z <- array(dim = c(nA, nB, T))
  
  for(t in 1:T) {
    Yt <- Y[,,t]
    
    if(family == "normal") {
      Z[,,t] <- Yt
    } else if(family == "binary") {
      # Use z-scores then center between classes
      Z[,,t] <- matrix(qnorm((Yt + 0.5)/2), nA, nB)
      
      # Find split point
      if(sum(Yt == 0, na.rm = TRUE) > 0 && sum(Yt == 1, na.rm = TRUE) > 0) {
        z0_max <- max(Z[,,t][Yt == 0], na.rm = TRUE)
        z1_min <- min(Z[,,t][Yt == 1], na.rm = TRUE)
        z_split <- (z0_max + z1_min) / 2
        Z[,,t] <- Z[,,t] - z_split
      }
    } else if(family == "poisson") {
      Z[,,t] <- log(Yt + 1)
    } else if(family %in% c("ordinal", "tobit")) {
      # Similar to unipartite
      Z[,,t] <- Yt
    }
  }
  
  # Initialize parameters
  a <- matrix(0, nA, T)
  b <- matrix(0, nB, T)
  
  # Initialize U and V
  if(RA > 0) {
    U <- array(rnorm(nA * RA * T, 0, 0.1), c(nA, RA, T))
  } else {
    U <- array(0, c(nA, 1, T))  # Minimal dimension
  }
  
  if(RB > 0) {
    V <- array(rnorm(nB * RB * T, 0, 0.1), c(nB, RB, T))
  } else {
    V <- array(0, c(nB, 1, T))  # Minimal dimension
  }
  
  # Initialize G as rectangular diagonal
  if(RA > 0 && RB > 0) {
    min_rank <- min(RA, RB)
    G <- matrix(0, RA, RB)
    diag(G)[1:min_rank] <- 1
  } else {
    G <- matrix(0, max(1,RA), max(1,RB))  # Minimal size
  }
  
  # Beta coefficients
  if(!is.null(Xlist)) {
    p <- dim(Xlist[[1]])[3]
    beta <- rep(0, p)
  } else {
    beta <- numeric(0)
  }
  
  # Variance parameters
  sigma2_a <- 1
  sigma2_b <- 1
  s2 <- 1
  
  # No rho for bipartite
  rho <- NULL
  
  list(
    Z = Z,
    a = a,
    b = b,
    U = U,
    V = V,
    G = G,
    beta = beta,
    sigma2_a = sigma2_a,
    sigma2_b = sigma2_b,
    s2 = s2,
    rho = rho
  )
}

#' Sample additive effects for bipartite networks
#' @keywords internal
sample_ab_bipartite <- function(Z, EZ_without_ab, sigma2_a, sigma2_b, s2) {
  dims <- dim(Z)
  nA <- dims[1]
  nB <- dims[2] 
  T <- dims[3]
  
  # Residuals without additive effects
  R <- Z - EZ_without_ab
  
  # Sample row effects (a)
  a <- matrix(0, nA, T)
  prec_a <- T / s2 + 1 / sigma2_a
  
  for(i in 1:nA) {
    for(t in 1:T) {
      mean_a <- sum(R[i,,t], na.rm = TRUE) / (nB * s2)
      mean_a <- mean_a / prec_a
      a[i,t] <- rnorm(1, mean_a, sqrt(1/prec_a))
    }
  }
  
  # Update residuals
  for(t in 1:T) {
    R[,,t] <- R[,,t] - a[,t]
  }
  
  # Sample column effects (b)
  b <- matrix(0, nB, T)
  prec_b <- T / s2 + 1 / sigma2_b
  
  for(j in 1:nB) {
    for(t in 1:T) {
      mean_b <- sum(R[,j,t], na.rm = TRUE) / (nA * s2)
      mean_b <- mean_b / prec_b
      b[j,t] <- rnorm(1, mean_b, sqrt(1/prec_b))
    }
  }
  
  list(a = a, b = b)
}

#' Update variance parameters for bipartite
#' @keywords internal
update_variances_bipartite <- function(a, b, eta0_a = 2, eta0_b = 2,
                                      Sab0_aa = 1, Sab0_bb = 1) {
  nA <- nrow(a)
  nB <- nrow(b)
  T <- ncol(a)
  
  # Sample row variance
  shape_a <- (eta0_a + nA * T) / 2
  rate_a <- (eta0_a * Sab0_aa + sum(a^2)) / 2
  sigma2_a <- 1 / rgamma(1, shape_a, rate_a)
  
  # Sample column variance  
  shape_b <- (eta0_b + nB * T) / 2
  rate_b <- (eta0_b * Sab0_bb + sum(b^2)) / 2
  sigma2_b <- 1 / rgamma(1, shape_b, rate_b)
  
  list(sigma2_a = sigma2_a, sigma2_b = sigma2_b)
}

#' Compute GOF statistics for bipartite networks
#' @keywords internal
compute_gof_bipartite <- function(Y_obs, Y_sim, family) {
  T <- dim(Y_obs)[3]
  
  gof_stats <- list()
  
  # Density
  gof_stats$density_obs <- apply(Y_obs, 3, function(y) mean(y > 0, na.rm = TRUE))
  gof_stats$density_sim <- apply(Y_sim, 3, function(y) mean(y > 0, na.rm = TRUE))
  
  # Degrees
  deg_bip <- compute_degrees_bip_cpp(Y_obs)
  gof_stats$row_degree_mean_obs <- colMeans(deg_bip$row_degrees)
  gof_stats$col_degree_mean_obs <- colMeans(deg_bip$col_degrees)
  
  deg_bip_sim <- compute_degrees_bip_cpp(Y_sim)
  gof_stats$row_degree_mean_sim <- colMeans(deg_bip_sim$row_degrees)
  gof_stats$col_degree_mean_sim <- colMeans(deg_bip_sim$col_degrees)
  
  # Four-cycles
  gof_stats$four_cycles_obs <- as.numeric(count_four_cycles_bip_cpp(Y_obs))
  gof_stats$four_cycles_sim <- as.numeric(count_four_cycles_bip_cpp(Y_sim))
  
  # Clustering
  gof_stats$clustering_obs <- as.numeric(clustering_coef_bip_cpp(Y_obs))
  gof_stats$clustering_sim <- as.numeric(clustering_coef_bip_cpp(Y_sim))
  
  gof_stats
}

#' Convert bipartite list data to array format
#' @keywords internal
list_to_array_bipartite <- function(rowActorSet, colActorSet, Y_list, 
                                   Xdyad = NULL, Xrow = NULL, Xcol = NULL) {
  nA <- length(rowActorSet)
  nB <- length(colActorSet)
  T <- length(Y_list)
  
  # Initialize arrays
  Y_array <- array(NA, c(nA, nB, T))
  
  # Fill Y array
  for(t in 1:T) {
    Yt <- Y_list[[t]]
    row_idx <- match(rownames(Yt), rowActorSet)
    col_idx <- match(colnames(Yt), colActorSet)
    
    # Only fill observed dyads
    for(i in seq_along(row_idx)) {
      for(j in seq_along(col_idx)) {
        if(!is.na(row_idx[i]) && !is.na(col_idx[j])) {
          Y_array[row_idx[i], col_idx[j], t] <- Yt[i, j]
        }
      }
    }
  }
  
  # Coerce Xdyad to 3D arrays if needed (use the same helpers from list_to_array.R)
  if(!is.null(Xdyad)) {
    Xdyad <- .coerce_Xdyad_3d(Xdyad)
  }
  
  # Handle covariates similarly
  # ... (implementation depends on your covariate structure)
  
  list(Y = Y_array, Xdyad = Xdyad, Xrow = Xrow, Xcol = Xcol)
}