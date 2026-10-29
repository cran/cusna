# =============================================================================
# Forward simulation of a SAOM at given parameters -- the SAOM counterpart of
# ergm_simulate() / tergm_simulate() / alaam_simulate(). It builds the same
# estimator context as mom_estimate() (data masks, effect descriptors, targets)
# and runs one batch of chains per period on the native simulator, so the
# statistics are exactly the moment statistics the estimator matches.
# =============================================================================

#' Simulate a SAOM at given parameters
#'
#' @param data A [saom_data] object.
#' @param effects List of [cusna_effect] / [cusna_interaction] specifications.
#' @param theta Parameter vector in the layout of `coef()` of the matching
#'   fit, or a `cusna_fit` object.
#' @param nsim Number of independent chains per period.
#' @param beh_effects Optional list of [cusna_beh_effect] (co-evolution).
#' @param conditional Condition each period on the observed distance.
#' @param maxdegree Optional out-degree bound (0 = none).
#' @param seed Optional seed for R's random number generator, from which the
#'   per-period simulation seeds are drawn.
#' @param return_nets Return the simulated end networks.
#' @return A `cusna_sim` object.
#' @export
saom_simulate <- function(data, effects, theta, nsim = 100, beh_effects = NULL,
                          conditional = FALSE, maxdegree = 0, seed = NULL,
                          return_nets = TRUE) {
  if (!inherits(data, "cusna_data")) {
    stop("`data` must be created with saom_data().", call. = FALSE)
  }
  if (!is.list(effects) || !length(effects)) {
    stop("`effects` must be a non-empty list of cusna_effect objects.", call. = FALSE)
  }
  nsim <- suppressWarnings(as.integer(nsim))
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L) {
    stop("`nsim` must be a positive integer.", call. = FALSE)
  }
  if (inherits(theta, "cusna_fit")) theta <- stats::coef(theta)
  bspecs <- if (!is.null(beh_effects))
    lapply(beh_effects, function(e) .make_beh_effect_spec(e$name)) else NULL
  est <- .mom_estimator(data$internal, .build_specs(effects), beh_effects = bspecs,
                        conditional = conditional,
                        control = list(maxdegree = as.integer(maxdegree)))
  theta <- as.numeric(theta)
  if (length(theta) != est$p || anyNA(theta)) {
    stop(sprintf("`theta` must have %d values (%s), the layout of coef() for this model.",
                 est$p, paste(est$param_names, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  s <- .mom_simulate(est, theta, nsim, keep_nets = isTRUE(return_nets))
  stats <- s$stats
  colnames(stats) <- est$param_names
  structure(list(stats = stats,
                 targets = stats::setNames(as.numeric(est$targets), est$param_names),
                 theta = stats::setNames(theta, est$param_names),
                 durations = if (est$conditional) s$ntim else NULL,
                 nets = s$nets, behavior = s$z_ends,
                 nsim = nsim, n_periods = est$P, conditional = est$conditional),
            class = "cusna_sim")
}

#' @export
print.cusna_sim <- function(x, ...) {
  cat(sprintf("cusna SAOM simulation: %d chains x %d period%s (%s)\n", x$nsim,
              x$n_periods, if (x$n_periods == 1L) "" else "s",
              if (x$conditional) "conditional" else "unconditional"))
  tab <- data.frame(statistic = names(x$targets), observed = x$targets,
                    sim_mean = colMeans(x$stats),
                    sim_sd = apply(x$stats, 2, stats::sd), row.names = NULL)
  tab$t <- (tab$sim_mean - tab$observed) / tab$sim_sd
  tab[, -1] <- round(tab[, -1], 3)
  print(tab, row.names = FALSE)
  if (!is.null(x$nets)) cat("end networks kept: $nets (one chains x n x n array per period)\n")
  invisible(x)
}
