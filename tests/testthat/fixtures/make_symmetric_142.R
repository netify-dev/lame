# symmetric fits saved by lame 1.4.2 for test_symmetric_node_saved.R; run with lame 1.4.2
# installed: Rscript make_symmetric_142.R <output .rds>
suppressMessages(library(lame))
stopifnot(identical(format(packageVersion("lame")), "1.4.2"))
set.seed(3); n = 10; nm = sprintf("a%02d", 1:n)
X = matrix(rnorm(n * 2), n, 2, dimnames = list(nm, c("x", "z")))
G = matrix(rnorm(n), n, 1, dimnames = list(nm, "growth"))
sym = function(m) { m = (m + t(m)) / 2; diag(m) = NA; dimnames(m) = list(nm, nm); m }
Y = sym(matrix(rnorm(n * n), n) + outer(X[, 1] + G[, 1], X[, 1] + G[, 1], "+"))
Xd = array((function(a) (a + t(a)) / 2)(matrix(rnorm(n * n), n)), c(n, n, 1), list(nm, nm, "d"))
o = list(Xd = Xd)
o$mcmc = ame(Y, Xdyad = Xd, Xrow = G, symmetric = TRUE, R = 0, burn = 20, nscan = 30, verbose = FALSE, plot = FALSE, gof = FALSE)
o$mcmc_pred = predict(o$mcmc)
YS = lapply(1:3, function(t) sym(Y + matrix(rnorm(n * n), n)))
o$lame_db = lame(YS, Xrow = rep(list(G), 3), symmetric = TRUE, R = 0, dynamic_beta = TRUE, burn = 20, nscan = 60, odens = 3, verbose = FALSE, plot = FALSE, gof = FALSE)
o$als = suppressWarnings(ame_als(Y, Xrow = X[, 1, drop = FALSE], Xcol = X[, 2, drop = FALSE], symmetric = TRUE, R = 0, verbose = FALSE))
o$als_refit = coef(suppressWarnings(ame_als_refit(o$als)))
o$als_pred = predict(o$als)
saveRDS(o, commandArgs(TRUE)[1], compress = "xz")
