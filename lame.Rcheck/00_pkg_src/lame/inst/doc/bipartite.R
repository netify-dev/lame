## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
	collapse = TRUE,
	comment = "#>",
	fig.width = 7,
	fig.height = 5,
	dpi = 96,
	out.width = "100%",
	fig.align = "center"
)

## ----setup--------------------------------------------------------------------
library(lame)
library(ggplot2)
set.seed(6886)

## ----tibble_bridge, eval = FALSE----------------------------------------------
# rx_net <- netify::netify(
# 	rx,
# 	actor1 = "patient_id", actor2 = "drug",
# 	weight = "prescribed",
# 	mode = "bipartite",
# 	symmetric = FALSE,
# 	missing_to_zero = FALSE
# )
# 
# fit <- ame(rx_net, family = "binary", mode = "bipartite", R = 2)

## ----simulate_cross_sectional-------------------------------------------------
# simulate a student-course enrollment network
n_students <- 30
n_courses <- 20

# true latent positions (unobserved in practice)
U_true <- matrix(rnorm(n_students * 2, 0, 0.8), n_students, 2)
V_true <- matrix(rnorm(n_courses * 2, 0, 0.8), n_courses, 2)

# interaction matrix: dimension 1 has positive affinity,
# dimension 2 has negative (students high on dim 2 avoid courses high on dim 2)
G_true <- matrix(c(1, 0.5, 0.5, -1), 2, 2)

# generate enrollment probabilities and binary outcomes
# negative intercept keeps enrollment rate realistic (~25-30%)
eta <- -0.8 + U_true %*% G_true %*% t(V_true)
prob <- pnorm(eta)
Y_bipartite <- matrix(rbinom(n_students * n_courses, 1, prob),
											n_students, n_courses)

rownames(Y_bipartite) <- paste0("Student", 1:n_students)
colnames(Y_bipartite) <- paste0("Course", 1:n_courses)

cat("Network dimensions:", dim(Y_bipartite), "\n")
cat("Enrollment rate:", round(mean(Y_bipartite), 2), "\n")

## ----fit_cross_sectional, message=FALSE, warning=FALSE------------------------
# burn and nscan are small so the vignette builds quickly
# for real analyses, use burn >= 1000 and nscan >= 5000.
fit_cross <- ame(
	Y = Y_bipartite,
	mode = "bipartite",
	R_row = 2,              # latent dimensions for students
	R_col = 2,              # latent dimensions for courses
	family = "binary",
	burn = 100,
	nscan = 500,
	odens = 5,
	verbose = FALSE,
	# save thinned U/V draws so latent_positions() can report posterior SDs
	posterior_opts = posterior_options(save_UV = TRUE)
)

summary(fit_cross)

## ----visualize_cross_sectional, fig.height=5, fig.alt="Bipartite latent-space biplot placing students and courses in the same 2D coordinate system, with proximity indicating a higher predicted enrollment probability."----
uv_plot(fit_cross, layout = "biplot") +
	ggtitle("Bipartite Latent Space: Students and Courses")

## ----interaction_heatmap, fig.height=6, fig.alt="Heatmap of the fitted multiplicative structure U G V' across student-course pairs using a diverging colour-blind-safe ColorBrewer RdBu palette with the scale capped at the 99th percentile of absolute fitted values: red tiles mark positive latent affinity (predicted enrollment above baseline), blue tiles mark negative latent affinity (predicted below baseline), and near-white tiles mark a posterior estimate near zero."----
G_mult <- fit_cross$U %*% fit_cross$G %*% t(fit_cross$V)
G_df <- data.frame(
	Student = rownames(G_mult)[row(G_mult)],
	Course = colnames(G_mult)[col(G_mult)],
	Fitted = as.vector(G_mult)
)
# order axes numerically (Student1, Student2, ...) instead of the default
# alphabetical factor order (Student1, Student10, Student11, ...)
G_df$Student <- factor(G_df$Student, levels = paste0("Student", n_students:1))
G_df$Course <- factor(G_df$Course, levels = paste0("Course", 1:n_courses))
# rdbu reversed (direction = -1): positive => red, negative => blue,
# zero => white. By default scale_fill_distiller() centres white at the
# midpoint of the DATA range, not at zero -- so we pass symmetric limits
# c(-mx, mx) to force the white midpoint onto zero, which is the right
# choice for a signed effect. We cap the fill at the 99th percentile of
# |UGV'| (clamping the few more extreme cells to the cap) so a single
# outlying cell cannot wash the rest of the map to near-white. RdBu is
# colour-blind-safe per ColorBrewer's protanopia / deuteranopia checks.
mx <- quantile(abs(G_df$Fitted), 0.99, na.rm = TRUE)
G_df$Fitted_capped <- pmax(pmin(G_df$Fitted, mx), -mx)
ggplot(G_df, aes(x = Course, y = Student, fill = Fitted_capped)) +
	geom_tile() +
	scale_fill_distiller(palette = "RdBu", direction = -1,
			limits = c(-mx, mx)) +
	labs(title = "Fitted Multiplicative Structure (UGV')",
			subtitle = "Red = positive, blue = negative; colour scale capped at the 99th percentile",
			x = "Course", y = "Student",
			fill = "U G V'") +
	theme_bw() +
	theme(panel.border = element_blank(),
			axis.ticks = element_blank(),
			legend.position = "top",
			axis.text.x = element_text(angle = 45, hjust = 1),
			axis.text.y = element_text(size = 6))

## ----simulate_longitudinal----------------------------------------------------
n_periods <- 5
n_users <- 20
n_items <- 15

true_rho <- 0.9
G_long <- matrix(c(1, 0.3, 0.3, -0.8), 2, 2)

U_t <- vector("list", n_periods)
V_t <- vector("list", n_periods)
U_t[[1]] <- matrix(rnorm(n_users * 2), n_users, 2)
V_t[[1]] <- matrix(rnorm(n_items * 2), n_items, 2)
for(t in 2:n_periods) {
	U_t[[t]] <- true_rho * U_t[[t-1]] +
		sqrt(1 - true_rho^2) * matrix(rnorm(n_users * 2), n_users, 2)
	V_t[[t]] <- true_rho * V_t[[t-1]] +
		sqrt(1 - true_rho^2) * matrix(rnorm(n_items * 2), n_items, 2)
}

Y_list <- list()
for(t in 1:n_periods) {
	eta_t <- U_t[[t]] %*% G_long %*% t(V_t[[t]])
	# continuous-valued interactions (normal family): a numerically stable
	# outcome for demonstrating the dynamic bipartite sampler end to end.
	Y_list[[t]] <- eta_t +
		matrix(rnorm(n_users * n_items, 0, 0.5), n_users, n_items)
	rownames(Y_list[[t]]) <- paste0("User", 1:n_users)
	colnames(Y_list[[t]]) <- paste0("Item", 1:n_items)
}
names(Y_list) <- paste0("T", 1:n_periods)

cat("Time periods:", length(Y_list), "\n")
cat("Dimensions per period:", dim(Y_list[[1]]), "\n")
cat("Average edge weight:", round(mean(sapply(Y_list, mean)), 2), "\n")
cat("True latent AR(1) coefficient:", true_rho, "\n")

## ----fit_longitudinal_static, message=FALSE-----------------------------------
# iterations are small so the vignette builds quickly; use burn >= 1000
# and nscan >= 5000 for real analyses.
fit_static <- lame(
	Y = Y_list,
	mode = "bipartite",
	R_row = 2,
	R_col = 2,
	family = "normal",
	dynamic_uv = FALSE,
	dynamic_ab = FALSE,
	burn = 100,
	nscan = 500,
	odens = 5,
	verbose = FALSE,
	plot = FALSE
)

summary(fit_static)

## ----fit_longitudinal_dynamic, message=FALSE----------------------------------
# same reduced iterations as above for vignette speed.
fit_dynamic <- lame(
	Y = Y_list,
	mode = "bipartite",
	R_row = 2,
	R_col = 2,
	family = "normal",
	dynamic_uv = TRUE,
	dynamic_ab = TRUE,
	burn = 100,
	nscan = 500,
	odens = 5,
	verbose = FALSE,
	plot = FALSE
)

summary(fit_dynamic)

## ----dynamic_rho_recovery-----------------------------------------------------
cat(sprintf("True rho_uv = %.2f | posterior mean = %.3f | 95%% CI = [%.3f, %.3f]\n",
						true_rho,
						mean(fit_dynamic$rho_uv),
						quantile(fit_dynamic$rho_uv, 0.025),
						quantile(fit_dynamic$rho_uv, 0.975)))

## ----fit_dynamic_G, message=FALSE---------------------------------------------
fit_dynG <- lame(
	Y = Y_list,
	mode = "bipartite",
	R_row = 2,
	R_col = 2,
	family = "normal",
	dynamic_uv = TRUE,
	dynamic_ab = TRUE,
	dynamic_G = TRUE,
	burn = 100,
	nscan = 500,
	odens = 5,
	verbose = FALSE,
	plot = FALSE
)

# fit$G_cube_post_mean has dimensions R_row x R_col x n_periods
if (!is.null(fit_dynG$G_cube_post_mean)) {
	cat("G_cube_post_mean dimensions:", dim(fit_dynG$G_cube_post_mean), "\n")
	cat("Frobenius norm of posterior-mean G_t per period:\n")
	print(round(apply(fit_dynG$G_cube_post_mean, 3,
	                  function(g) sqrt(sum(g^2))), 3))
	cat("Rotation-drift ratio (>= 5 flags latent-rotation domination):",
	    round(fit_dynG$G_rotation_drift$ratio, 2), "\n")
} else {
	cat("G_cube was not attached (no bipartite + RA>0 + RB>0 case).\n")
}

## ----visualize_longitudinal, fig.height=5, fig.alt="Trajectory plot tracing each actor's latent position from the first to the last period; paths drift gradually and stay clustered rather than scattering."----
uv_plot(fit_dynamic, plot_type = "trajectory", label.nodes = FALSE) +
	guides(color = "none") +
	ggtitle("Dynamic latent trajectories")

## ----visualize_longitudinal_highlight, fig.height=5, fig.alt="Trajectory plot with three highlighted actors (User1, User5, Item3) coloured and named in the legend while all other trajectories are greyed out for context."----
uv_plot(fit_dynamic, plot_type = "trajectory",
	highlight = c("User1", "User5", "Item3"),
	label.nodes = FALSE) +
	ggtitle("Highlighted actor trajectories")

## ----trace_dynamic, fig.height=6, fig.alt="MCMC trace plots and posterior densities for the sampled dynamic bipartite model parameters (intercept and variance components); well-mixed caterpillar traces indicate adequate sampling."----
trace_plot(fit_dynamic,
	exclude = c("Dyadic Correlation", "Sender-Receiver Covariance"))

## ----model_comparison, fig.height=4, fig.alt="Side-by-side goodness-of-fit panels comparing posterior predictive coverage of bipartite sender and receiver degree heterogeneity for the static versus dynamic model."----
# gof plots for each model; four-cycles is omitted because it is
# constant on this dense continuous network (see text above)
gof_plot(fit_static, statistics = c("sd.row", "sd.col"))
gof_plot(fit_dynamic, statistics = c("sd.row", "sd.col"))

## ----model_comparison_boxplot, fig.height=4, fig.alt="Boxplots of posterior predictive bipartite statistics (row-mean SD and column-mean SD) faceted by static and dynamic models, with horizontal reference lines at observed values."----
gof_static <- fit_static$GOF
gof_dynamic <- fit_dynamic$GOF

# extract posterior predictive samples (exclude column 1, which is observed)
# each element is a matrix [n_time x n_mcmc], so colMeans averages across time
extract_ppc <- function(gof, stat) {
	mat <- gof[[stat]]
	colMeans(mat[, -1, drop = FALSE])
}

# observed values (column 1, averaged across time periods)
obs_vals <- c(
	mean(gof_static$sd.rowmean[, 1]),
	mean(gof_static$sd.colmean[, 1])
)

gof_df <- data.frame(
	value = c(
		extract_ppc(gof_static, "sd.rowmean"),
		extract_ppc(gof_dynamic, "sd.rowmean"),
		extract_ppc(gof_static, "sd.colmean"),
		extract_ppc(gof_dynamic, "sd.colmean")
	),
	model = rep(rep(c("Static", "Dynamic"),
							each = ncol(gof_static$sd.rowmean) - 1), 2),
	statistic = rep(c("Row Mean SD", "Column Mean SD"),
									each = 2 * (ncol(gof_static$sd.rowmean) - 1))
)

obs_df <- data.frame(
	statistic = c("Row Mean SD", "Column Mean SD"),
	observed = obs_vals
)

ggplot(gof_df, aes(x = model, y = value)) +
	geom_boxplot() +
	geom_hline(data = obs_df, aes(yintercept = observed),
						linetype = 2) +
	facet_wrap(~statistic, scales = "free_y") +
	labs(title = "GOF Comparison: Static vs Dynamic",
			subtitle = "Dashed Line = Observed Value",
			x = "Model", y = "Simulated Statistic") +
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.ticks = element_blank(),
		legend.position = "top",
		strip.background = element_rect(fill = "black", color = "black"),
		strip.text = element_text(color = "white", hjust = 0)
	)

## ----dimension_selection, message=FALSE, warning=FALSE------------------------
# include R=0 as a no-latent-space baseline
dims_to_test <- list(
	c(0, 0),
	c(1, 1),
	c(2, 2)
)

# compute observed GOF statistics from the data
obs_gof <- gof_stats(Y_bipartite, mode = "bipartite")

gof_results <- list()
for(i in seq_along(dims_to_test)) {
	fit_temp <- ame(
		Y = Y_bipartite,
		mode = "bipartite",
		R_row = dims_to_test[[i]][1],
		R_col = dims_to_test[[i]][2],
		family = "binary",
		burn = 200,
		nscan = 2500,
		odens = 5,
		verbose = FALSE
	)
	# posterior predictive tail probability: proportion of simulated
	# statistics at or above the observed value (one-sided)
	gof_results[[i]] <- c(
		R_row = dims_to_test[[i]][1],
		R_col = dims_to_test[[i]][2],
		pval_rowmean = mean(fit_temp$GOF[, "sd.rowmean"] >= obs_gof["sd.rowmean"]),
		pval_colmean = mean(fit_temp$GOF[, "sd.colmean"] >= obs_gof["sd.colmean"]),
		pval_fourcycles = mean(fit_temp$GOF[, "four.cycles"] >= obs_gof["four.cycles"])
	)
}

do.call(rbind, gof_results)

## -----------------------------------------------------------------------------
lp <- latent_positions(fit_cross)
head(lp)

# filter to just the row nodes (students)
lp_students <- lp[lp$type == "U", ]
cat("Student positions:", nrow(lp_students), "rows\n")

