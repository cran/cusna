# R4 acceptance (determinism): the cpp11 wrapper cusna_simulate_period_cpu_
# marshals cusna_saom_problem identically to the reference ctypes transcription.
# sim.json.gz was frozen from one native cusna_saom_simulate_period_cpu call per
# config (baseline / covariate / coevolution, fixed seeds); the wrapper must
# reproduce it bit-for-bit (integer trajectory) / to float32 (scores).

skip_if_not_installed("jsonlite")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "sim", "sim.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)
as_mat <- function(rows) do.call(rbind, lapply(rows, as.numeric))
waves <- lapply(1:3, function(i) as_mat(g$inputs$s50_waves[[i]]))
alc <- as_mat(g$inputs$alcohol)
rm2 <- function(x) as.vector(t(as.matrix(x)))

run_wrap <- function(cfg, data, effects, beh_effects = NULL) {
  n <- data$n; k <- length(effects); M <- cfg$M; p <- cfg$m + 1L
  beh <- data$beh
  z_start <- if (!is.null(beh)) .beh_values_less_missing_starts(beh, p) else NULL
  pe <- .prepare_effects_v2(effects, n, p, beh, z_start)
  wtab <- .effect_wtabs(effects, n)
  kb <- if (!is.null(beh_effects)) length(beh_effects) else 0L
  bty <- if (kb > 0) vapply(beh_effects, function(e) e$type_id, integer(1)) else integer(0)
  cusna_simulate_period_cpu_(n, k, M, cfg$seed, as.integer(rm2(data$imputed[[p]])),
    cfg$lam, as.numeric(cfg$theta), as.integer(pe$type_ids), rep(0L, k),
    as.integer(pe$dyn_flags), as.integer(pe$dyc_idx), rm2(pe$vmat), rm2(wtab),
    as.numeric(pe$v_dyn), as.numeric(pe$p1), as.numeric(pe$p2),
    as.integer(rm2(pe$vmiss)), 0L, 0L, 0L, if (isTRUE(data$symmetric)) 1L else 0L,
    kb, if (kb > 0) cfg$lamb else 0, if (kb > 0) as.numeric(cfg$thetab) else numeric(0),
    as.integer(bty), if (kb > 0) as.integer(beh$imputed[, p]) else integer(0),
    if (kb > 0) as.integer(.beh_kmask_dist(beh, p)) else integer(0),
    if (kb > 0) beh$zmin else 0L, if (kb > 0) beh$zmax else 0L,
    if (kb > 0) beh$overall_mean else 0, if (kb > 0) beh$range else 0,
    if (kb > 0) beh$sim_mean else 0, 0L)
}
E <- .make_effect_spec

test_that("wrapper reproduces the native baseline simulation bit-for-bit", {
  d <- .build_saom_data(waves)
  r <- run_wrap(g$baseline, d, lapply(c("density","recip","transTrip"), E))
  M <- g$baseline$M
  expect_identical(as.integer(r$dist), as.integer(unlist(g$baseline$dist)))
  expect_identical(as.integer(r$dist_stat), as.integer(unlist(g$baseline$dist_stat)))
  expect_identical(as.integer(r$steps), as.integer(unlist(g$baseline$steps)))
  expect_equal(matrix(r$stats, M, byrow = TRUE), as_mat(g$baseline$stats), tolerance = 1e-9)
  expect_equal(matrix(r$scores, M, byrow = TRUE), as_mat(g$baseline$scores), tolerance = 2e-4)
})

test_that("wrapper reproduces the native covariate simulation", {
  d <- .build_saom_data(waves)
  effc <- c(lapply(c("density","recip","transTrip"), E),
            list(E("egoX", alc[,1]), E("altX", alc[,1]), E("simX", alc[,1])))
  r <- run_wrap(g$covariate, d, effc)
  M <- g$covariate$M
  expect_identical(as.integer(r$dist), as.integer(unlist(g$covariate$dist)))
  expect_identical(as.integer(r$steps), as.integer(unlist(g$covariate$steps)))
  expect_equal(matrix(r$stats, M, byrow = TRUE), as_mat(g$covariate$stats), tolerance = 1e-9)
})

test_that("wrapper reproduces the native coevolution simulation (net + behavior)", {
  d <- .build_saom_data(waves, behavior = alc)
  effv <- c(lapply(c("density","recip","transTrip"), E),
            list(E("egoX", dyn=TRUE), E("altX", dyn=TRUE), E("simX", dyn=TRUE)))
  bef <- lapply(c("linear","quad","avSim"), .make_beh_effect_spec)
  r <- run_wrap(g$coevol, d, effv, bef)
  M <- g$coevol$M
  expect_identical(as.integer(r$dist), as.integer(unlist(g$coevol$dist)))
  expect_identical(as.integer(r$bdist), as.integer(unlist(g$coevol$bdist)))
  expect_equal(matrix(r$z_end, M, byrow = TRUE), as_mat(g$coevol$z_end), tolerance = 1e-9)
  expect_equal(matrix(r$stats, M, byrow = TRUE), as_mat(g$coevol$stats), tolerance = 1e-9)
})
