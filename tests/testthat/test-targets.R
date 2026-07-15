# R3 acceptance: the R moment-target assembly (.mom_targets over the native
# cusna_network_stats / cusna_behavior_stats carriers) reproduces the Python
# MomEstimator targets. The golden targets.json.gz was frozen from MomEstimator
# on baseline (uncond/cond), covariate, coevolution, rate-effect and the real
# Knecht / COW configs; it carries the raw inputs so R rebuilds identically.
# Since Python == RSiena to machine zero (A4.6), R == RSiena transitively.

skip_if_not_installed("jsonlite")

g <- jsonlite::fromJSON(gzfile(test_path("fixtures", "targets", "targets.json.gz")),
                        simplifyVector = TRUE, simplifyDataFrame = FALSE,
                        simplifyMatrix = FALSE)

as_mat <- function(rows) do.call(rbind, lapply(rows, function(r)
  vapply(r, function(x) if (is.null(x)) NA_real_ else as.numeric(x), double(1))))
alc_const <- as.numeric(unlist(g$inputs$alc_const))
alc_change <- as_mat(g$inputs$alc_change)
resolve_cov <- function(key) switch(key, alc_const = alc_const,
                                    alc_change = alc_change, stop("cov ", key))
build_eff <- function(r) .make_effect_spec(name = r$name,
  cov = if (!is.null(r$cov)) resolve_cov(r$cov) else NULL, parm = r$parm,
  dyn = isTRUE(r$dyn), etype = if (is.null(r$etype)) "eval" else r$etype,
  net_ref = r$net_ref)

eqm <- function(R, gold_rows, tol = 1e-7) {
  if (is.null(gold_rows)) return(invisible())
  go <- as_mat(gold_rows)
  expect_equal(matrix(as.numeric(R), nrow(go), ncol(go)), go, tolerance = tol)
}
eqv <- function(R, gold, tol = 1e-7) {
  if (is.null(gold)) return(invisible())
  expect_equal(as.numeric(unlist(R)), as.numeric(unlist(gold)), tolerance = tol)
}

check_config <- function(cfg) {
  nw <- length(cfg$waves_input)
  waves <- lapply(seq_len(nw), function(m) as_mat(cfg$waves_input[[m]]))
  behavior <- if (!is.null(cfg$behavior_input)) as_mat(cfg$behavior_input) else NULL
  data <- .build_saom_data(waves, behavior = behavior,
                           composition = cfg$composition_input)
  effects <- lapply(cfg$effects, build_eff)
  beh_effects <- if (!is.null(cfg$beh_effects))
    lapply(unlist(cfg$beh_effects), .make_beh_effect_spec) else NULL
  rate_effects <- if (!is.null(cfg$rate_effects)) lapply(cfg$rate_effects,
    function(r) .make_rate_effect_spec(r$name,
      if (!is.null(r$cov)) resolve_cov(r$cov) else NULL)) else NULL
  mt <- .mom_targets(data, effects, beh_effects, rate_effects, cfg$conditional)
  o <- cfg$out
  eqm(mt$eval, o$eval_targets)
  eqv(mt$basic_rate, o$rate_targets)
  eqm(mt$beh, o$beh_targets)
  eqv(mt$beh_rate, o$bdists)
  if (!is.null(mt$rate_eff)) eqv(colSums(mt$rate_eff), o$rate_eff_targets)
  eqv(mt$target_vector, o$target_vector)
}

for (cfg in g$configs) {
  local({
    cc <- cfg
    test_that(sprintf("moment targets match Python MomEstimator (%s)", cc$name), {
      check_config(cc)
    })
  })
}
