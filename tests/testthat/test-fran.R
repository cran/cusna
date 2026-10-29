# cusna_fran() must honour RSiena's simstats0c contract: siena07() sums `fra`
# over periods with colSums() and reads `sc` one row per period, so both are
# periods x parameters matrices, and `ntim0` holds one duration per period.
# cusna 0.1.0 returned `fra` as the period-summed vector and `sc` for the last
# period only, so siena07() stopped at its first call (in colSums(zz$fra)).

skip_if_not_installed("jsonlite")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "estimate", "estimate.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)
as_mat <- function(rows) do.call(rbind, lapply(rows, as.numeric))
waves <- lapply(1:3, function(i) as_mat(g$inputs$s50_waves[[i]]))
effs <- c("density", "recip", "transTrip")

test_that("fra and sc are periods x parameters matrices (conditional)", {
  fran <- cusna_fran(waves, effs)
  set.seed(1)
  r <- fran(list(theta = c(-2.5, 2.4, 0.6), cconditional = TRUE), NULL)
  expect_true(r$OK)
  expect_true(is.matrix(r$fra))
  expect_equal(dim(r$fra), c(2L, 3L))
  expect_true(is.matrix(r$sc))
  expect_equal(dim(r$sc), c(2L, 3L))
  expect_length(r$ntim0, 2L)
  expect_true(all(r$ntim0 > 0))
  expect_length(colSums(r$fra), 3L)           # what siena07() computes first
})

test_that("unconditional rows carry each period's rate statistic in its own column", {
  fran <- cusna_fran(waves, c("density", "recip"), conditional = FALSE)
  set.seed(2)
  r <- fran(list(theta = c(6, 5, -2.5, 2.4), cconditional = FALSE), NULL)
  expect_equal(dim(r$fra), c(2L, 4L))
  expect_equal(dim(r$sc), c(2L, 4L))
  expect_equal(r$fra[1, 2], 0)
  expect_equal(r$fra[2, 1], 0)
  expect_gt(r$fra[1, 1], 0)
  expect_gt(r$fra[2, 2], 0)
})

test_that("fra and sc rows are the per-period pieces of the native simulator", {
  fran <- cusna_fran(waves, effs)
  theta <- c(-2.5, 2.4, 0.6)
  set.seed(3)
  r <- fran(list(theta = theta), NULL)
  est <- .mom_estimator(.build_saom_data(waves),
                        .build_specs(lapply(effs, cusna_effect)), conditional = TRUE)
  set.seed(3)
  s <- .mom_simulate(est, theta, 1L)
  for (pp in 1:2) {
    expect_equal(r$fra[pp, ], s$pieces[[pp]]$sp[1, ])
    expect_equal(r$sc[pp, ], s$pieces[[pp]]$scp[1, ])
  }
  expect_equal(colSums(r$fra), s$stats[1, ])
  expect_equal(r$ntim0, s$ntim[1, ])
})

test_that("the closure counts its calls and leaves the RNG stream to RSiena", {
  fran <- cusna_fran(waves, effs)
  set.seed(4)
  r1 <- fran(list(theta = c(-2.5, 2.4, 0.6)), NULL)
  r2 <- fran(list(theta = c(-2.5, 2.4, 0.6)), NULL)
  expect_equal(attr(fran, "counter")$i, 2L)
  expect_false(identical(r1$fra, r2$fra))     # successive calls draw new seeds
  set.seed(4)
  r1b <- fran(list(theta = c(-2.5, 2.4, 0.6)), NULL)
  expect_identical(r1$fra, r1b$fra)           # ... reproducibly
})

test_that("a theta that does not match the effects is rejected with a clear message", {
  fran <- cusna_fran(waves, c("density", "recip"))
  expect_error(fran(list(theta = c(1, 2, 3)), NULL), "3 parameters")
  expect_error(fran(list(theta = c(-2, 1), cconditional = FALSE), NULL),
               "conditional")
})

test_that("siena07() runs end to end on the cusna FRAN backend", {
  skip_if_not_installed("RSiena")
  skip_on_cran()
  net <- RSiena::sienaDependent(array(unlist(waves), dim = c(50, 50, 3)))
  dat <- RSiena::sienaDataCreate(net)
  eff <- RSiena::includeEffects(RSiena::getEffects(dat), "transTrip",
                                character = TRUE, verbose = FALSE)
  alg <- RSiena::sienaAlgorithmCreate(projname = NULL, nsub = 2, n3 = 200,
                                      seed = 11, cond = TRUE)
  alg$FRAN <- cusna_fran(waves, effs)
  utils::capture.output(
    ans <- RSiena::siena07(alg, data = dat, effects = eff, batch = TRUE,
                           verbose = FALSE, silent = TRUE, useCluster = FALSE))
  expect_s3_class(ans, "sienaFit")
  expect_length(ans$theta, 3L)
  expect_true(all(is.finite(ans$theta)))
  expect_true(all(is.finite(ans$se)))
  # the s50 reference estimates are about (-2.7, 2.4, 0.63); a short run lands near
  expect_true(all(abs(ans$theta - c(-2.7, 2.4, 0.63)) < 4 * ans$se + 0.2))
  expect_gt(attr(alg$FRAN, "counter")$i, 200L)
})
