# CRAN's two-core policy: run the OpenMP engine single-threaded here.
library(cusna)
old_threads <- cusna_set_threads(1L)
have_rsiena <- requireNamespace("RSiena", quietly = TRUE)

library(RSiena)

net <- sienaDependent(array(unlist(s50$waves), dim = c(50, 50, 3)))
dat <- sienaDataCreate(net)
eff <- includeEffects(getEffects(dat), transTrip)

alg <- sienaAlgorithmCreate(projname = NULL, nsub = 2, n3 = 300, seed = 2026)
alg$FRAN <- cusna_fran(
  waves        = s50$waves,                         # the same 0/1 wave matrices
  effect_names = c("density", "recip", "transTrip"),# in the RSiena effects order
  conditional  = TRUE)                              # must match `cond`

ans <- siena07(alg, data = dat, effects = eff, batch = TRUE,
               verbose = FALSE, silent = TRUE, useCluster = FALSE)
ans                              # a normal sienaFit
attr(alg$FRAN, "counter")$i      # simulation calls made by siena07()

cusna_set_threads(old_threads)

