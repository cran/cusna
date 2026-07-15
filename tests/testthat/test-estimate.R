# R4 acceptance (convergence): the R Robbins-Monro estimator (.mom_estimator +
# .mom_run) lands within SE of the Python cpu.py estimates. estimate.json.gz
# was frozen from MomEstimator.estimate on s50 with a fast schedule; the R run
# uses a different RNG stream (so within-SE, not bit-exact) and the float32 ABI
# scores. This test drives many native simulations -- opt-in (slow).

skip_if_not_installed("jsonlite")
skip_if(Sys.getenv("CUSNA_RUN_ESTIMATE_TESTS") == "",
        "slow: set CUSNA_RUN_ESTIMATE_TESTS=1 to run the RM estimation tests")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "estimate", "estimate.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)
as_mat <- function(rows) do.call(rbind, lapply(rows, as.numeric))
waves <- lapply(1:3, function(i) as_mat(g$inputs$s50_waves[[i]]))
alc <- as_mat(g$inputs$alcohol)
S <- g$settings
ctl <- list(firstg = S$firstg, nsub = as.integer(S$nsub), n2 = as.numeric(unlist(S$n2)),
            batch2 = as.integer(S$batch2), n1 = as.integer(S$n1), n3 = as.integer(S$n3),
            nD = as.integer(S$nD), diagonalize = S$diagonalize, seed = as.integer(S$seed),
            maxdegree = 0L)
E <- .make_effect_spec

expect_within_se <- function(fit, gold) {
  th_r <- fit$theta; th_p <- as.numeric(unlist(gold$theta)); se <- as.numeric(unlist(gold$se))
  for (i in seq_along(th_p)) {
    expect_lt(abs(th_r[i] - th_p[i]), 1.0 * se[i] + 0.1)
  }
}

test_that("R RM estimate matches Python within SE (baseline, unconditional)", {
  d <- .build_saom_data(waves)
  est <- .mom_estimator(d, lapply(c("density","recip","transTrip"), E), control = ctl)
  expect_within_se(.mom_run(est), g$configs$baseline)
})

test_that("R RM estimate matches Python within SE (baseline, conditional)", {
  d <- .build_saom_data(waves)
  est <- .mom_estimator(d, lapply(c("density","recip","transTrip"), E),
                        conditional = TRUE, control = ctl)
  expect_within_se(.mom_run(est), g$configs$baseline_cond)
})

test_that("R RM estimate matches Python within SE (covariate)", {
  d <- .build_saom_data(waves)
  effc <- c(lapply(c("density","recip","transTrip"), E),
            list(E("egoX", alc[,1]), E("altX", alc[,1]), E("simX", alc[,1])))
  est <- .mom_estimator(d, effc, control = ctl)
  expect_within_se(.mom_run(est), g$configs$covariate)
})

test_that("R RM estimate matches Python within SE (coevolution)", {
  d <- .build_saom_data(waves, behavior = alc)
  effv <- c(lapply(c("density","recip","transTrip"), E),
            list(E("egoX", dyn=TRUE), E("altX", dyn=TRUE), E("simX", dyn=TRUE)))
  bef <- lapply(c("linear","quad","avSim"), .make_beh_effect_spec)
  est <- .mom_estimator(d, effv, beh_effects = bef, control = ctl)
  expect_within_se(.mom_run(est), g$configs$coevol)
})
