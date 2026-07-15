library(cusna)
library(RSiena)

# ... assemble `dat` (sienaData) and `eff` (sienaEffects) as usual ...

alg <- sienaAlgorithmCreate(projname = NULL, cond = FALSE)

alg$FRAN <- cusna_fran(
  waves        = list(w1, w2, w3),                  # the same 0/1 wave matrices
  effect_names = c("density", "recip", "transTrip"),# in the RSiena effects order
  conditional  = FALSE)                             # must match `cond`

ans <- siena07(alg, data = dat, effects = eff, useCluster = FALSE)
ans   # a normal sienaFit: estimates, standard errors, convergence as usual

