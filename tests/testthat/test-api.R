# R5 acceptance: the public API (saom_data / cusna_effect & friends /
# mom_control / mom_estimate / cusna_fit methods) builds and validates without
# Python. Constructor/validation checks always run; the end-to-end fit (which
# drives the native simulator) is opt-in (slow).

skip_if_not_installed("jsonlite")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "estimate", "estimate.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)
as_mat <- function(rows) do.call(rbind, lapply(rows, as.numeric))
waves <- lapply(1:3, function(i) as_mat(g$inputs$s50_waves[[i]]))
alc <- as_mat(g$inputs$alcohol)

test_that("saom_data builds a cusna_data with the internal engine", {
  d <- saom_data(waves)
  expect_s3_class(d, "cusna_data")
  expect_equal(d$n, 50L)
  expect_equal(d$n_periods, 2L)
  expect_true(inherits(d$internal, "cusna_saom_internal"))
  expect_output(print(d), "50 actors")
})

test_that("effect constructors return typed spec objects", {
  expect_s3_class(cusna_effect("density"), "cusna_effect")
  expect_s3_class(cusna_effect("egoX", covariate = alc[, 1]), "cusna_effect")
  expect_s3_class(cusna_beh_effect("linear"), "cusna_beh_effect")
  expect_s3_class(cusna_rate_effect("outRate"), "cusna_rate_effect")
  expect_s3_class(cusna_interaction(c(1, 2)), "cusna_interaction")
  expect_error(cusna_interaction(1), "2 or 3")
  expect_error(cusna_effect(c("a", "b")), "single effect name")
})

test_that("mom_control carries the schedule", {
  ctl <- mom_control(n1 = 100, seed = 9)
  expect_s3_class(ctl, "cusna_control")
  expect_equal(ctl$n1, 100L)
  expect_equal(ctl$seed, 9L)
})

test_that("mom_estimate validates its inputs (P-11 messages)", {
  expect_error(mom_estimate(list(1), list(cusna_effect("density"))),
               "saom_data")
  expect_error(mom_estimate(saom_data(waves), list()), "non-empty")
  expect_error(mom_estimate(saom_data(waves), list(cusna_effect("density")),
                            backend = "gpu"), "gpu")
})

test_that("saom_multinet_data builds a multi-network container", {
  mn <- saom_multinet_data(list(waves, waves))
  expect_s3_class(mn, "cusna_multinet_data")
  expect_equal(mn$n_networks, 2L)
  expect_equal(mn$internal$D, 2L)
})

test_that("cusna_fran returns a working FRAN closure (simstats0c contract)", {
  fran <- cusna_fran(waves, c("density", "recip"))
  expect_true(is.function(fran))
  r <- fran(list(theta = c(-2, 1)), NULL)
  expect_true(r$OK)
  expect_equal(dim(r$fra), c(2L, 2L))  # periods x eval effects (conditional)
  expect_equal(dim(r$sc), c(2L, 2L))
  expect_length(r$ntim0, 2)            # one duration per period
})

test_that("mom_estimate produces a cusna_fit with working methods", {
  skip_if(Sys.getenv("CUSNA_RUN_ESTIMATE_TESTS") == "",
          "slow: set CUSNA_RUN_ESTIMATE_TESTS=1")
  d <- saom_data(waves)
  ctl <- mom_control(n1 = 150, nsub = 1, n2 = 8, batch2 = 100, nD = 150, n3 = 300, seed = 7)
  fit <- mom_estimate(d, list(cusna_effect("density"), cusna_effect("recip"),
                              cusna_effect("transTrip")), control = ctl)
  expect_s3_class(fit, "cusna_fit")
  expect_length(coef(fit), 5)                       # 2 rates + 3 eval
  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_s3_class(as.data.frame(fit), "data.frame")
  expect_named(coef(fit))
  # delta-method vcov diagonal == se^2
  expect_equal(unname(sqrt(diag(vcov(fit)))), unname(fit$se), tolerance = 1e-6)
  expect_output(print(fit), "convergence ratio")
  expect_output(summary(fit), "overall max. convergence ratio")
})
