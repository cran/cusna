# =============================================================================
# Moment targets -- pure-R assembly over the native statistic carriers
# (track A-RAPI, milestone R3). Applies the R1 period transforms to the observed
# end waves and calls the validated native cusna_network_stats /
# cusna_behavior_stats (bit-exact host float64), yielding the per-period
# evaluation / behavior / rate targets the Robbins-Monro loop (R4) matches.
#
# Bar (R3): eval/behavior targets == Python MomEstimator (same native carriers)
# and == RSiena bench targets to machine precision on baseline / covariate /
# coevolution / real data (Knecht, COW).
#
# typed (creation/endowment) and cross-network (multinet) targets need a cpp11
# wrapper over cusna_typed_stats / cusna_cross_stats (both already in the C ABI);
# added in R3b together with the R4 simulator glue.
# =============================================================================

# Per-period evaluation-effect targets (P x k): apply the statistic transform to
# the observed end wave, then call the native float64 host statistic.
.eval_targets <- function(data, effects, symflag = NULL) {
  n <- data$n
  P <- data$n_periods
  k <- length(effects)
  beh <- data$beh
  if (is.null(symflag)) symflag <- if (isTRUE(data$symmetric)) 1L else 0L
  wtab <- .effect_wtabs(effects, n)
  out <- matrix(0.0, P, k)
  for (p in seq_len(P)) {
    z_start <- if (!is.null(beh)) .beh_values_less_missing_starts(beh, p) else NULL
    pe <- .prepare_effects_v2(effects, n, p, beh, z_start)
    a_start <- data$imputed[[p]]
    a_tgt <- .saom_target_network(data, p)
    od0 <- as.integer(rowSums(a_start))
    id0 <- as.integer(colSums(a_start))
    out[p, ] <- cusna_network_stats(
      a_tgt, type_ids = pe$type_ids, p1 = pe$p1, p2 = pe$p2, vmat = pe$vmat,
      vmiss = pe$vmiss, dyn_flags = pe$dyn_flags, vb = pe$v_dyn,
      od0 = od0, id0 = id0, wtab = wtab, symflag = symflag)
  }
  out
}

# Per-period behavior-effect targets (P x kb) on the observed end behavior wave.
.beh_targets <- function(data, beh_effects) {
  beh <- data$beh
  if (is.null(beh)) stop("behavior targets need co-evolving behavior data.")
  P <- data$n_periods
  kb <- length(beh_effects)
  bty <- vapply(beh_effects, function(e) as.integer(e$type_id), integer(1))
  out <- matrix(0.0, P, kb)
  for (p in seq_len(P)) {
    z_end <- beh$imputed[, p + 1L]
    xnet <- .saom_network_less_missing_start(data, p)
    bkdist <- as.integer(.beh_kmask_dist(beh, p))
    out[p, ] <- cusna_behavior_stats(
      z_end, bty = bty, xnet = xnet, bkdist = bkdist,
      zbar = beh$overall_mean, brange = beh$range, bsimmean = beh$sim_mean)
  }
  out
}

# Basic-rate targets per period. Composition-change gate (mom.py): WITH
# composition RSiena feeds the C++ StatisticCalculator convention (rate_target,
# imputed values / asymmetric masks); WITHOUT it siena07 reports the observed
# Hamming distance (attr(depvar,"distance")). They agree on complete data and
# differ only when a network is missing at a period end (e.g. Knecht).
.basic_rate_targets <- function(data) {
  fn <- if (isTRUE(data$has_composition)) .saom_rate_target else .saom_observed_distance
  vapply(seq_len(data$n_periods),
         function(p) as.numeric(fn(data, p)), numeric(1))
}

# Non-basic rate-effect target per period: the per-ego tie differences of the
# period weighted by fixed per-actor weights (covariate for RateX, degrees of
# wave m less its own missings for outRate/inRate).
.rate_effect_targets <- function(data, rate_effects) {
  P <- data$n_periods
  nr <- length(rate_effects)
  out <- matrix(0.0, P, nr)
  for (p in seq_len(P)) {
    rows <- .saom_rate_target_rows(data, p)              # per-ego differences
    for (r in seq_len(nr)) {
      e <- rate_effects[[r]]
      w <- switch(as.character(e$kind),
        "0" = .rate_cov_values(e, p),                     # RateX (centered cov)
        "1" = rowSums(.saom_network_less_missing(data, p)),   # outRate: out-deg
        "2" = colSums(.saom_network_less_missing(data, p)))   # inRate:  in-deg
      out[p, r] <- sum(as.numeric(w) * rows)
    }
  }
  out
}

# Assemble the full moment-target bundle + the theta-ordered target vector.
# Layouts (matching mom.py / the R4 loop):
#   unconditional        : [basic rate x P, non-basic rate x nr, eval x k]
#   conditional          : [non-basic rate x nr, eval x k]  (basic rate is the
#                          stopping target, not a moment)
#   coevolution (uncond) : [net rate x P, beh rate x P, eval x k, beh x kb]
.mom_targets <- function(data, effects, beh_effects = NULL, rate_effects = NULL,
                         conditional = TRUE) {
  P <- data$n_periods
  eval_t <- .eval_targets(data, effects)
  basic_rate <- .basic_rate_targets(data)
  rate_eff <- if (!is.null(rate_effects)) .rate_effect_targets(data, rate_effects) else NULL
  beh_t <- if (!is.null(beh_effects)) .beh_targets(data, beh_effects) else NULL
  bh_rate <- if (!is.null(beh_t)) .beh_rate_targets(data) else NULL

  # column sums are the aggregated moment targets (RM matches the batch mean)
  eval_vec <- colSums(eval_t)
  parts <- list()
  if (!is.null(beh_t)) {
    # coevolution: net basic rate per period, behavior basic rate per period
    parts$net_rate <- basic_rate
    parts$beh_rate <- bh_rate
    parts$eval <- eval_vec
    parts$beh <- colSums(beh_t)
  } else {
    if (!conditional) parts$rate <- basic_rate
    if (!is.null(rate_eff)) parts$rate_eff <- colSums(rate_eff)
    parts$eval <- eval_vec
  }
  list(eval = eval_t, beh = beh_t, basic_rate = basic_rate,
       rate_eff = rate_eff, beh_rate = bh_rate,
       target_vector = unlist(parts, use.names = FALSE),
       layout = names(parts))
}

# Behavior basic-rate target per period (observed behavior distance).
.beh_rate_targets <- function(data) {
  beh <- data$beh
  vapply(seq_len(data$n_periods),
         function(p) as.numeric(.beh_observed_distance(beh, p)), numeric(1))
}
