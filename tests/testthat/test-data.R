# R1 acceptance: the pure-R data engine (.build_saom_data + .saom_* / .beh_*
# carriers) reproduces the Python data.py conventions bit-for-bit. The golden
# fixtures under fixtures/data/<name>.json.gz were frozen from data.py on s50,
# s50 (missing/structural injections), Knecht, COW alliances and a synthetic
# panel; they carry the raw inputs so the R side rebuilds from identical data.

skip_if_not_installed("jsonlite")

fix_dir <- test_path("fixtures", "data")

rd <- function(name) {
  jsonlite::fromJSON(gzfile(file.path(fix_dir, paste0(name, ".json.gz"))),
                     simplifyVector = TRUE, simplifyDataFrame = FALSE,
                     simplifyMatrix = FALSE)
}
as_mat <- function(rows) {
  if (is.matrix(rows)) return(rows)
  do.call(rbind, lapply(rows, as.numeric))
}
# bit-exact integer comparison (masks / imputed / targets), dim-agnostic
eq_int <- function(Rval, gold) {
  expect_equal(as.integer(as.vector(as.matrix(Rval) * 1L)),
               as.integer(as.vector(as_mat(gold))))
}
eq_ivec <- function(Rval, gold) {
  expect_equal(as.integer(Rval), as.integer(unlist(gold)))
}
eq_flt <- function(Rval, gold, tol = 1e-9) {
  expect_equal(as.numeric(Rval), as.numeric(unlist(gold)), tolerance = tol)
}

check_dataset <- function(name) {
  g <- rd(name)
  meta <- g$meta
  nw <- meta$n_waves
  P <- meta$n_periods

  waves <- lapply(seq_len(nw), function(m) as_mat(g$waves_input[[m]]))
  behavior <- if (!is.null(g$behavior_input)) as_mat(g$behavior_input) else NULL
  composition <- g$composition_input

  obj <- .build_saom_data(waves, behavior = behavior,
                          composition = composition, cc_option = meta$cc_option)

  expect_identical(as.logical(obj$symmetric), as.logical(meta$symmetric))
  expect_identical(obj$n_periods, as.integer(meta$n_periods))

  for (m in seq_len(nw)) {
    eq_int(obj$imputed[[m]], g$imputed[[m]])
    eq_int(obj$miss[[m]], g$miss[[m]])
    eq_int(obj$struct[[m]], g$struct[[m]])
    eq_int(obj$rawvals[[m]], g$rawvals[[m]])
  }
  if (!is.null(g$comp_active)) eq_int(obj$comp_active, g$comp_active)

  for (p in seq_len(P)) {
    gp <- g$periods[[p]]
    eq_int(.saom_kmask_dist(obj, p), gp$kmask_dist)
    expect_equal(.saom_observed_distance(obj, p), gp$observed_distance)
    expect_equal(.saom_rate_target(obj, p), gp$rate_target)
    expect_equal(.saom_rate_const(obj, p), gp$rate_const)
    eq_ivec(.saom_rate_target_rows(obj, p), gp$rate_target_rows)
    eq_ivec(.saom_rate_const_rows(obj, p), gp$rate_const_rows)
    eq_int(.saom_network_less_missing(obj, p), gp$network_less_missing)
    eq_int(.saom_network_less_missing_start(obj, p), gp$network_less_missing_start)
    kf <- .saom_stat_keep_force(obj, p)
    eq_int(kf$keep, gp$keep)
    eq_int(kf$force, gp$force)
    kfs <- .saom_stat_keep_force_struct(obj, p)
    eq_int(kfs$keep, gp$keepS)
    eq_int(kfs$force, gp$forceS)
    eq_int(.saom_target_network(obj, p), gp$target_network)
    eq_int(.saom_permitted(obj, p), gp$permitted)
  }

  if (!is.null(g$behavior)) {
    bg <- g$behavior
    beh <- obj$beh
    expect_equal(beh$zmin, bg$zmin)
    expect_equal(beh$zmax, bg$zmax)
    eq_flt(beh$overall_mean, bg$overall_mean, 1e-12)
    eq_flt(beh$range, bg$range, 1e-12)
    eq_flt(beh$sim_mean, bg$sim_mean, 1e-12)
    eq_int(beh$imputed, bg$imputed)
    eq_int(beh$miss, bg$miss)
    for (p in seq_len(P)) {
      bp <- bg$periods[[p]]
      eq_ivec(as.integer(.beh_kmask_dist(beh, p)), bp$kmask_dist)
      expect_equal(.beh_observed_distance(beh, p), bp$observed_distance)
      eq_flt(.beh_values_less_missing_starts(beh, p),
             bp$values_less_missing_starts, 1e-12)
      eq_flt(.beh_stat_values(beh, beh$imputed[, p + 1L], p),
             bp$stat_values_target, 1e-12)
    }
  }
}

test_that("R data engine is bit-identical to Python data.py (s50)", {
  check_dataset("s50")
})
test_that("R data engine is bit-identical to Python data.py (s50 + behavior)", {
  check_dataset("s50_beh")
})
test_that("R data engine reproduces missing masks (s50 missing injections)", {
  check_dataset("s50miss")
})
test_that("R data engine reproduces structural masks (s50 struct injections)", {
  check_dataset("s50str")
})
test_that("R data engine reproduces Knecht (structural + missing + coevolution)", {
  check_dataset("knecht")
})
test_that("R data engine reproduces COW (symmetric + composition, n=116)", {
  check_dataset("cow")
})
test_that("R data engine reproduces the synthetic panel", {
  check_dataset("synth")
})
