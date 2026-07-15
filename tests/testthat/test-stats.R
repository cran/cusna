# A5 smoke tests: the native host statistics agree with hand-computed values.
# The heavy machine-zero / RSiena validation lives in native/test (the C ABI is
# the same); these just prove the R -> cpp11 -> native path is wired correctly.

test_that("abi/cuda introspection works", {
  expect_true(cusna_abi_version() >= 1L)
  expect_type(cusna_has_cuda(), "logical")
})

test_that("cusna_set_threads caps the OpenMP thread count", {
  prev <- cusna_set_threads(1L)
  expect_type(prev, "integer")
  after <- cusna_set_threads(prev)
  # after == 1L when OpenMP is compiled in (we just set it to 1 above);
  # 0L on a build without OpenMP support (the call is then a no-op).
  expect_true(after == 1L || after == 0L)
})

test_that("network stats reproduce hand-computed density/recip/transTrip", {
  # 3 actors, ties 1->2, 2->1, 2->3 (1-indexed)
  a <- matrix(0L, 3, 3)
  a[1, 2] <- 1L; a[2, 1] <- 1L; a[2, 3] <- 1L
  k <- 3L; n <- 3L
  out <- cusna_network_stats(
    a, type_ids = c(0L, 1L, 2L),          # density, recip, transTrip
    p1 = rep(0, k), p2 = rep(0, k),
    vmat = matrix(0, k, n), vmiss = matrix(0L, k, n),
    dyn_flags = rep(0L, k), vb = rep(0, n),
    od0 = as.integer(rowSums(a)), id0 = as.integer(colSums(a)),
    wtab = matrix(0, k, n + 1L), symflag = 0L)
  expect_equal(out[1], 3)   # density (3 ties)
  expect_equal(out[2], 2)   # recip (1<->2, both directions)
  expect_equal(out[3], 0)   # transTrip (no transitive triple)
})

test_that("behavior linear/quad match the centered sums", {
  # 4 actors, behavior values, complete star predictor network
  z <- matrix(c(1L, 2L, 3L, 2L), nrow = 1)   # one chain
  xnet <- matrix(0L, 4, 4)
  xnet[1, 2] <- 1L; xnet[1, 3] <- 1L; xnet[1, 4] <- 1L
  bkdist <- rep(1L, 4)
  zbar <- 2.0
  out <- cusna_behavior_stats(z, bty = c(0L, 1L), xnet = xnet,
                              bkdist = bkdist, zbar = zbar,
                              brange = 2.0, bsimmean = 0.0)
  u <- c(1, 2, 3, 2) - zbar                    # -1, 0, 1, 0
  expect_equal(out[1, 1], sum(u))              # linear
  expect_equal(out[1, 2], sum(u^2))            # quad
})
