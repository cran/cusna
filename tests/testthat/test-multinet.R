# R4b acceptance: the multi-network Robbins-Monro estimator (mom_estimate_
# multinet over the native cusna_saom_simulate_multi_cpu + host network/cross
# statistics) lands within SE of Python MultiMomEstimator. multinet.json.gz was
# frozen from a converging 2-network model on s50-derived panels. Slow -- opt-in.

skip_if_not_installed("jsonlite")
skip_if(Sys.getenv("CUSNA_RUN_ESTIMATE_TESTS") == "",
        "slow: set CUSNA_RUN_ESTIMATE_TESTS=1 to run the multi-network estimation")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "multinet", "multinet.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)
as_mat <- function(rows) do.call(rbind, lapply(rows, as.numeric))
netA <- lapply(g$netA_waves, as_mat); netB <- lapply(g$netB_waves, as_mat)
S <- g$settings
ctl <- mom_control(firstg = S$firstg, nsub = as.integer(S$nsub), n2 = as.numeric(unlist(S$n2)),
  batch2 = as.integer(S$batch2), n1 = as.integer(S$n1), n3 = as.integer(S$n3),
  nD = as.integer(S$nD), diagonalize = S$diagonalize, seed = as.integer(S$seed))

test_that("multi-network estimate matches Python within SE", {
  mn <- saom_multinet_data(list(netA, netB))
  ebn <- list(list(cusna_effect("density"), cusna_effect("recip")),
              list(cusna_effect("density"), cusna_effect("recip")))
  fit <- mom_estimate_multinet(mn, ebn, control = ctl)
  expect_s3_class(fit, "cusna_fit")
  expect_length(coef(fit), 6)                          # 2 nets x (rate + density + recip)
  th_p <- as.numeric(unlist(g$theta)); se <- as.numeric(unlist(g$se))
  for (i in seq_along(th_p)) {
    expect_lt(abs(coef(fit)[i] - th_p[i]), 1.0 * se[i] + 0.2)
  }
})
