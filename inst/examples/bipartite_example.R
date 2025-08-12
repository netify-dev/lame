# Example of using bipartite network support in lame package

library(lame)

# Simulate a simple bipartite network
set.seed(6886)

# Dimensions
nA <- 10  # Number of row nodes (e.g., users)
nB <- 15  # Number of column nodes (e.g., items)
T <- 5    # Number of time periods

# Generate synthetic bipartite network data
Y_list <- list()
for(t in 1:T) {
  # Create adjacency matrix with some structure
  Y_t <- matrix(0, nA, nB)
  
  # Add some baseline connectivity
  baseline_prob <- 0.2
  Y_t <- matrix(rbinom(nA * nB, 1, baseline_prob), nA, nB)
  
  # Add some row effects (some users more active)
  row_effects <- rnorm(nA, 0, 0.5)
  for(i in 1:nA) {
    Y_t[i,] <- Y_t[i,] | (runif(nB) < plogis(row_effects[i]))
  }
  
  # Add some column effects (some items more popular)
  col_effects <- rnorm(nB, 0, 0.5)
  for(j in 1:nB) {
    Y_t[,j] <- Y_t[,j] | (runif(nA) < plogis(col_effects[j]))
  }
  
  # Add row and column names
  rownames(Y_t) <- paste0("User", 1:nA)
  colnames(Y_t) <- paste0("Item", 1:nB)
  
  Y_list[[t]] <- Y_t
}
names(Y_list) <- paste0("T", 1:T)

# Display basic network statistics
cat("Bipartite Network Statistics:\n")
cat("Dimensions:", nA, "x", nB, "over", T, "time periods\n")
cat("Density by time period:\n")
densities <- sapply(Y_list, function(y) mean(y))
print(round(densities, 3))

# Fit bipartite LAME model
cat("\nFitting bipartite LAME model...\n")
cat("(Note: Use more iterations in practice)\n")

# This would be the call to fit the model:
# fit_bip <- lame(
#   Y = Y_list,
#   mode = "bipartite",
#   R_row = 2,  # 2 latent dimensions for row nodes
#   R_col = 3,  # 3 latent dimensions for column nodes
#   family = "binary",
#   nscan = 1000,
#   burn = 100,
#   odens = 10
# )

# Compute GOF statistics using the helper functions
Y_array <- array(0, c(nA, nB, T))
for(t in 1:T) {
  Y_array[,,t] <- Y_list[[t]]
}

# Four-cycle counts
if(exists("count_four_cycles_bip_cpp")) {
  cycles <- lame:::count_four_cycles_bip_cpp(Y_array)
  cat("\nFour-cycle counts by time period:\n")
  print(as.numeric(cycles))
}

# Degree distributions
if(exists("compute_degrees_bip_cpp")) {
  degrees <- lame:::compute_degrees_bip_cpp(Y_array)
  cat("\nMean row degrees (users) by time:\n")
  print(round(colMeans(degrees$row_degrees), 2))
  cat("\nMean column degrees (items) by time:\n")
  print(round(colMeans(degrees$col_degrees), 2))
}

cat("\nBipartite network example complete.\n")
cat("The package now supports rectangular networks with different dimensions.\n")