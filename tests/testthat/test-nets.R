# mom_estimate(return_nets = TRUE) keeps the phase-3 end networks, as RSiena's
# siena07(returnDeps = TRUE) does. cusna 0.1.0 always called the simulator with
# want_ends = 0, so no simulated network could be recovered from a fit.

skip_if_not_installed("jsonlite")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "estimate", "estimate.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)
as_mat <- function(rows) do.call(rbind, lapply(rows, as.numeric))
waves <- lapply(1:3, function(i) as_mat(g$inputs$s50_waves[[i]]))
alc <- as_mat(g$inputs$alcohol)

# Recompute the evaluation statistics of one end network with the validated host
# statistic (the same call .eval_targets makes); s50 has no missing ties, so the
# statistic network is the end network itself.
net_stats <- function(a, effects, d, pp) {
  pe <- .prepare_effects_v2(effects, d$n, pp, NULL, NULL)
  x0 <- d$imputed[[pp]]
  cusna_network_stats(a, pe$type_ids, pe$p1, pe$p2, pe$vmat, pe$vmiss,
                      pe$dyn_flags, pe$v_dyn, as.integer(rowSums(x0)),
                      as.integer(colSums(x0)), .effect_wtabs(effects, d$n))
}

test_that("the period simulator returns end networks that match its statistics", {
  d <- .build_saom_data(waves)
  # egoX is orientation-sensitive, so a transposed network would fail here
  effects <- .build_specs(list(cusna_effect("density"), cusna_effect("recip"),
                               cusna_effect("egoX", covariate = alc[, 1])))
  est <- .mom_estimator(d, effects, conditional = TRUE)
  set.seed(5)
  s <- .mom_simulate(est, c(-2.5, 2.4, 0.1), 12L, keep_nets = TRUE)
  expect_length(s$nets, 2L)
  for (pp in 1:2) {
    a <- s$nets[[pp]]
    expect_equal(dim(a), c(12L, 50L, 50L))
    expect_true(is.integer(a))
    for (c in 1:12) {
      m <- a[c, , ]
      expect_true(all(diag(m) == 0L))
      # conditional simulation stops exactly at the observed distance
      expect_equal(sum(m != d$imputed[[pp]]), est$dists[pp])
      expect_equal(net_stats(m, effects, d, pp), s$pieces[[pp]]$sp[c, ],
                   tolerance = 1e-12)
    }
  }
  # networks are only materialised on request
  set.seed(5)
  expect_null(.mom_simulate(est, c(-2.5, 2.4, 0.1), 4L)$nets)
})

test_that("mom_estimate(return_nets = TRUE) returns the phase-3 networks", {
  ctl <- mom_control(n1 = 60, nsub = 1, n2 = 6, batch2 = 30, nD = 60, n3 = 40,
                     seed = 3)
  effects <- list(cusna_effect("density"), cusna_effect("recip"))
  fit <- mom_estimate(saom_data(waves), effects, conditional = TRUE,
                      control = ctl, return_nets = TRUE)
  expect_length(fit$sims, 2L)
  for (pp in 1:2) {
    expect_equal(dim(fit$sims[[pp]]), c(40L, 50L, 50L))
    expect_true(all(fit$sims[[pp]] %in% c(0L, 1L)))
  }
  # they are the phase-3 draws: their statistics average to the reported means
  ties <- function(a) apply(a, 1, sum)
  recip <- function(a) apply(a, 1, function(m) sum(m * t(m)))
  expect_equal(mean(ties(fit$sims[[1]]) + ties(fit$sims[[2]])),
               unname(fit$sim_means[1]))
  expect_equal(mean(recip(fit$sims[[1]]) + recip(fit$sims[[2]])),
               unname(fit$sim_means[2]))
  # and feed the GoF primitive directly
  indeg <- cusna_gof_distribution(fit$sims[[2]], kind = 0L, levls = 0:8)
  expect_equal(dim(indeg), c(40L, 9L))
  expect_output(print(fit), "phase-3 networks")

  fit0 <- mom_estimate(saom_data(waves), effects, conditional = TRUE,
                       control = ctl)
  expect_null(fit0$sims)
  expect_error(mom_estimate(saom_data(waves), effects, return_nets = NA),
               "return_nets")
})

test_that("co-evolution fits also return the simulated behavior", {
  ctl <- mom_control(n1 = 60, nsub = 1, n2 = 6, batch2 = 30, nD = 60, n3 = 20,
                     seed = 4)
  fit <- mom_estimate(saom_data(waves, behavior = alc),
                      effects = list(cusna_effect("density"), cusna_effect("recip")),
                      beh_effects = list(cusna_beh_effect("linear")),
                      control = ctl, return_nets = TRUE)
  expect_length(fit$sims, 2L)
  expect_length(fit$sims_behavior, 2L)
  expect_equal(dim(fit$sims_behavior[[1]]), c(20L, 50L))
  expect_true(all(fit$sims_behavior[[2]] %in% 1:5))
})

test_that("saom_simulate() simulates at a fixed theta and returns the networks", {
  d <- saom_data(waves)
  effects <- list(cusna_effect("density"), cusna_effect("recip"))
  sim <- saom_simulate(d, effects, theta = c(-2.6, 2.4), nsim = 25,
                       conditional = TRUE, seed = 9)
  expect_s3_class(sim, "cusna_sim")
  expect_equal(dim(sim$stats), c(25L, 2L))
  expect_equal(colnames(sim$stats), c("density", "recip"))
  expect_equal(dim(sim$durations), c(25L, 2L))
  expect_equal(dim(sim$nets[[1]]), c(25L, 50L, 50L))
  ties <- apply(sim$nets[[1]], 1, sum) + apply(sim$nets[[2]], 1, sum)
  expect_equal(ties, unname(sim$stats[, 1]))
  expect_output(print(sim), "observed")
  # reproducible from the seed, and identical for any number of threads
  old <- cusna_set_threads(1L)
  sim1 <- saom_simulate(d, effects, theta = c(-2.6, 2.4), nsim = 25,
                        conditional = TRUE, seed = 9)
  cusna_set_threads(old)
  expect_identical(sim1$stats, sim$stats)
  expect_identical(sim1$nets, sim$nets)
  # the theta layout is checked (unconditional: two rates + two effects)
  expect_error(saom_simulate(d, effects, theta = c(1, 2, 3)), "4 values")
  expect_null(saom_simulate(d, effects, c(-2.6, 2.4), nsim = 3, conditional = TRUE,
                            return_nets = FALSE)$nets)
})
