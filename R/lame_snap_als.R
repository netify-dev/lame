# the snap kind of dynamic_uv has no ALS route: its drift scale would be read
# from paths the penalty has already shrunk. called by lame() for snap with
# method = "als"
.lame_snap_als_refuse <- function(
	shape = c("directed", "bipartite", "symmetric")
) {
	shape <- match.arg(shape)
	reason <- switch(
		shape,
		directed = "On directed snap panels a snap ALS fit reads its drift scale from paths the penalty has already shrunk, so it either tracks the latent product less closely than the {.val ar1} kind or misses resets.",
		bipartite = "On bipartite snap panels a snap ALS fit reads its drift scale from paths the penalty has already shrunk, so it tracks the latent product less closely than the {.val ar1} kind.",
		symmetric = "On symmetric snap panels a snap ALS fit reads its drift scale from paths the penalty has already shrunk, so it tracks the latent product less closely than the {.val ar1} kind."
	)
	use <- "Use {.code dynamic_uv_kind = \"snap\"} with {.code method = \"mcmc\"}, or {.code dynamic_uv_kind = \"ar1\"} with {.code method = \"als\"}."
	cli::cli_abort(
		c(
			"{.code dynamic_uv_kind = \"snap\"} with {.code method = \"als\"} is not available.",
			"i" = reason,
			"i" = use
		),
		call = parent.frame()
	)
}
