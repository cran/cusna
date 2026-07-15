# =============================================================================
# siena07() acceleration backend (track A-RAPI, milestone R5). cusna_fran()
# returns a FRAN closure that plugs the native simulator into RSiena's
# unmodified estimator. The full RSiena bridge (mapping siena07's effects object
# and theta layout onto the cusna specs) is a follow-up; the closure here builds
# the native estimator context and simulates one batch per Robbins-Monro call.
# =============================================================================

#' Native FRAN backend for RSiena's siena07()
#'
#' Returns a FRAN closure that drives the native SAOM simulator inside an
#' unmodified `RSiena::siena07()` run:
#'
#' \preformatted{
#'   alg <- RSiena::sienaAlgorithmCreate(...)
#'   alg$FRAN <- cusna_fran(waves, effect_names)
#'   RSiena::siena07(alg, data = dat, effects = eff, useCluster = FALSE)
#' }
#'
#' @param waves A list of adjacency matrices (the panel).
#' @param effect_names Character vector of effect names, in the order RSiena's
#'   `theta` uses (after the basic rates).
#' @param covariate Optional covariate shared by the covariate effects.
#' @param conditional Conditional (RSiena default) or unconditional simulation.
#' @return A FRAN closure for `alg$FRAN`.
#' @export
cusna_fran <- function(waves, effect_names, covariate = NULL, conditional = TRUE) {
  data <- saom_data(waves)
  cov_effects <- c("egoX", "altX", "simX", "sameX", "higher", "egoXaltX",
                   "egoPlusAltX", "altRThresholdX")
  effects <- lapply(effect_names, function(nm)
    cusna_effect(nm, covariate = if (nm %in% cov_effects) covariate else NULL))
  est <- .mom_estimator(data$internal, .build_specs(effects),
                        conditional = conditional, control = unclass(mom_control()))
  counter <- new.env(parent = emptyenv()); counter$i <- 0L

  function(z, x, ...) {
    counter$i <- counter$i + 1L
    theta <- as.numeric(z$theta)
    set.seed(counter$i)
    s <- .mom_simulate(est, theta, 1L)
    list(fra = as.numeric(s$stats[1, ]),
         sc = matrix(s$pieces[[length(s$pieces)]]$scp[1, ], nrow = 1),
         ntim0 = as.numeric(s$ntim[1, ]), OK = TRUE, feasible = TRUE)
  }
}
