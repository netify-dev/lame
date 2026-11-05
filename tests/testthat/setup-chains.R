# the suite fits single chains unless a test asks for more: the package default is two chains
old_lame_chains = options(lame.chains = 1L)
withr::defer(options(old_lame_chains), teardown_env())
