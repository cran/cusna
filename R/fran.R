# =============================================================================
# siena07() acceleration backend (track A-RAPI, milestone R5). cusna_fran()
# returns a FRAN closure that plugs the native simulator into RSiena's
# unmodified estimator. The full RSiena bridge (mapping siena07's effects object
# and theta layout onto the cusna specs) is a follow-up; the closure here builds
# the native estimator context and simulates one chain per Robbins-Monro call.
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
  conditional <- isTRUE(conditional)
  est <- .mom_estimator(data$internal, .build_specs(effects),
                        conditional = conditional, control = unclass(mom_control()))
  counter <- new.env(parent = emptyenv()); counter$i <- 0L

  fran <- function(z, x, ...) {
    theta <- as.numeric(z$theta)
    if (!is.null(z$cconditional) && !identical(isTRUE(z$cconditional), conditional)) {
      stop(sprintf("cusna_fran(conditional = %s) does not match the algorithm (cond = %s).",
                   conditional, isTRUE(z$cconditional)), call. = FALSE)
    }
    if (length(theta) != est$p) {
      stop(sprintf(paste0("siena07() passed %d parameters, but the cusna model has %d ",
                          "(%s); check `effect_names` against the RSiena effects object."),
                   length(theta), est$p, paste(est$param_names, collapse = ", ")),
           call. = FALSE)
    }
    counter$i <- counter$i + 1L
    # one chain; its per-period seeds come from R's RNG stream, which siena07()
    # seeds from the algorithm, so a run is reproducible under a fixed `seed`
    s <- .mom_simulate(est, theta, 1L)
    # simstats0c contract: siena07() forms colSums(fra) and reads sc row by row,
    # so both are periods x parameters, one row per period
    fra <- do.call(rbind, lapply(s$pieces, function(pc) pc$sp[1L, ]))
    sc <- do.call(rbind, lapply(s$pieces, function(pc) pc$scp[1L, ]))
    list(fra = fra, sc = sc, ntim0 = as.numeric(s$ntim[1L, ]),
         OK = TRUE, feasible = TRUE)
  }
  attr(fran, "counter") <- counter
  fran
}
