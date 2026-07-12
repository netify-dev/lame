####
# bipartite network example
library(lame)
set.seed(6886)
####

####
# dimensions
nA = 10
nB = 15
T = 5

# generate bipartite network data
Y_list = list()
for(t in 1:T) {
	Y_t = matrix(rbinom(nA * nB, 1, 0.2), nA, nB)

	# row effects
	row_effects = rnorm(nA, 0, 0.5)
	for(i in 1:nA) {
		Y_t[i,] = Y_t[i,] | (runif(nB) < plogis(row_effects[i]))
	}

	# col effects
	col_effects = rnorm(nB, 0, 0.5)
	for(j in 1:nB) {
		Y_t[,j] = Y_t[,j] | (runif(nA) < plogis(col_effects[j]))
	}

	rownames(Y_t) = paste0("User", 1:nA)
	colnames(Y_t) = paste0("Item", 1:nB)
	Y_list[[t]] = Y_t
}
names(Y_list) = paste0("T", 1:T)
####

####
# network statistics
cat("Bipartite Network Statistics:\n")
cat("Dimensions:", nA, "x", nB, "over", T, "time periods\n")
cat("Density by time period:\n")
densities = sapply(Y_list, function(y) mean(y))
print(round(densities, 3))
####

####
# gof statistics
Y_array = array(0, c(nA, nB, T))
for(t in 1:T) {
	Y_array[,,t] = Y_list[[t]]
}

# four-cycle counts
if(exists("count_four_cycles_bip_cpp")) {
	cycles = lame:::count_four_cycles_bip_cpp(Y_array)
	cat("\nFour-cycle counts by time period:\n")
	print(as.numeric(cycles))
}

# degree distributions
if(exists("compute_degrees_bip_cpp")) {
	degrees = lame:::compute_degrees_bip_cpp(Y_array)
	cat("\nMean row degrees (users) by time:\n")
	print(round(colMeans(degrees$row_degrees), 2))
	cat("\nMean column degrees (items) by time:\n")
	print(round(colMeans(degrees$col_degrees), 2))
}

cat("\nBipartite network example complete.\n")
####
