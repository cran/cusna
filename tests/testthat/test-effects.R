# R2 acceptance: the pure-R effect preprocessing (.make_effect_spec +
# .prepare_effects_v2 + .interaction_comp_idx / .effect_etypes / .effect_netrefs
# / .effect_wtabs / .typed_pair_matrix + .gwesp_weights) reproduces the Python
# stats.py descriptors. The golden fixture effects.json.gz was frozen from
# prepare_effects_v2 & friends over 16 configs covering every branch (structural,
# degree, actor covariates const/changing/missing, gwesp/gwdsp + integer coercion,
# dyadic const/changing, dependent-behavior dyn, interactions, cross-network,
# typed creation/endowment); wtab is the float64 host table (gwesp_weights).

skip_if_not_installed("jsonlite")

g <- jsonlite::fromJSON(
  gzfile(test_path("fixtures", "effects", "effects.json.gz")),
  simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
n <- g$n

as_mat <- function(rows) do.call(rbind, lapply(rows, function(r)
  vapply(r, function(x) if (is.null(x)) NA_real_ else as.numeric(x), double(1))))

inp <- g$inputs
v_const <- as.numeric(unlist(inp$v_const))
v_change <- as_mat(inp$v_change)
v_miss <- vapply(inp$v_miss, function(x) if (is.null(x)) NA_real_ else as.numeric(x),
                 double(1))
dyad_const <- as_mat(inp$dyad_const)
dyad_change <- array(0.0, c(n, n, 3))
for (m in 1:3) dyad_change[, , m] <- as_mat(g$dyad_change[[m]])
beh <- .build_behavior_data(v_change, n, 3L)

resolve_cov <- function(key) switch(key,
  v_const = v_const, v_change = v_change, v_miss = v_miss,
  dyad_const = dyad_const, dyad_change = dyad_change, stop("unknown cov ", key))

build_effect <- function(r) {
  if (isTRUE(r$interaction)) {
    return(.make_interaction_spec(as.integer(unlist(r$comps))))
  }
  cov <- if (!is.null(r$cov)) resolve_cov(r$cov) else NULL
  .make_effect_spec(name = r$name, cov = cov, parm = r$parm, dyn = isTRUE(r$dyn),
                    etype = if (is.null(r$etype)) "eval" else r$etype,
                    net_ref = r$net_ref)
}

eqi <- function(R, gold) expect_equal(as.integer(R), as.integer(unlist(gold)))
eqf <- function(R, gold, tol = 1e-9) {
  expect_equal(as.numeric(R), as.numeric(unlist(gold)), tolerance = tol)
}
eqm <- function(R, gold_rows, tol = 1e-9) {
  go <- as_mat(gold_rows)
  expect_equal(matrix(as.numeric(R), nrow(go), ncol(go)), go, tolerance = tol)
}

check_config <- function(cfg) {
  effects <- lapply(cfg$recipe, build_effect)
  p <- cfg$period + 1L
  use_beh <- isTRUE(cfg$use_beh)
  z_start <- if (use_beh) .beh_values_less_missing_starts(beh, p) else NULL
  pe <- .prepare_effects_v2(effects, n, p, if (use_beh) beh else NULL, z_start)
  eqi(pe$type_ids, cfg$type_ids)
  eqf(pe$p1, cfg$p1)
  eqf(pe$p2, cfg$p2)
  eqm(pe$vmat, cfg$vmat)
  eqi(pe$dyn_flags, cfg$dyn_flags)
  eqf(pe$v_dyn, cfg$v_dyn)
  eqi(pe$dyc_idx, cfg$dyc_idx)
  eqi(as.integer(t(pe$vmiss)), cfg$vmiss)              # Python ravel is row-major
  if (!is.null(cfg$dymats)) {
    for (j in seq_along(cfg$dymats)) eqm(pe$dymats[[j]], cfg$dymats[[j]])
  }
  eqm(.effect_wtabs(effects, n), cfg$wtab)
  eqi(as.integer(t(.interaction_comp_idx(effects))),
      as.integer(unlist(lapply(cfg$int_idx, function(r) as.integer(unlist(r))))))
  eqi(.effect_etypes(effects), cfg$etypes)
  eqi(.effect_netrefs(effects), cfg$netrefs)
  if (!is.null(cfg$typed_F)) {
    for (j in seq_along(effects)) {
      eqm(.typed_pair_matrix(effects[[j]], n, p), cfg$typed_F[[j]])
    }
  }
}

for (cfg in g$configs) {
  local({
    cc <- cfg
    test_that(sprintf("effect descriptors match Python (%s)", cc$name), {
      check_config(cc)
    })
  })
}

test_that("gwesp_weights matches Python for several parm", {
  for (parm in names(g$gwesp_weights)) {
    expect_equal(.gwesp_weights(as.numeric(parm), n),
                 as.numeric(unlist(g$gwesp_weights[[parm]])), tolerance = 1e-12)
  }
})
