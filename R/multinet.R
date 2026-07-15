# =============================================================================
# Multi-network co-evolution (track A-RAPI, milestones R5 + R4b). saom_multinet_
# data() builds the container; mom_estimate_multinet() runs the unconditional
# per-variable Robbins-Monro over the native multi-network simulator. The data
# container ships in R5; the estimator body lands in R4b.
# =============================================================================

# Internal multi-network data: one internal SaomData per network, shared actor
# set and wave count (port of data.py MultiNetData).
.build_multinet_data <- function(waves_list, composition = NULL, cc_option = 1) {
  if (!is.list(waves_list) || length(waves_list) < 2L) {
    stop("`waves_list` needs >= 2 networks.", call. = FALSE)
  }
  nets <- lapply(waves_list, function(net)
    .build_saom_data(net, composition = composition, cc_option = cc_option))
  n0 <- nets[[1]]$n; p0 <- nets[[1]]$n_periods
  for (d in seq_along(nets)) {
    if (nets[[d]]$n != n0) stop(sprintf("network %d: actor set size differs", d), call. = FALSE)
    if (nets[[d]]$n_periods != p0) stop(sprintf("network %d: wave count differs", d), call. = FALSE)
  }
  list(nets = nets, D = length(nets), n = n0, n_periods = p0)
}

#' Multi-network co-evolution data
#'
#' Several dependent networks over one actor set, with cross-network effects.
#' Estimation is unconditional only (as in RSiena for more than one dependent
#' variable).
#'
#' @param waves_list A list of >= 2 networks, each itself a list of wave matrices.
#' @param composition Optional shared composition-change intervals (see
#'   [saom_data]).
#' @param cc_option The sienaCompositionChange option.
#' @return A `cusna_multinet_data` object.
#' @export
saom_multinet_data <- function(waves_list, composition = NULL, cc_option = 1) {
  if (!is.list(waves_list) || length(waves_list) < 2L) {
    stop("`waves_list` must be a list of >= 2 networks, each a list of wave ",
         "matrices.", call. = FALSE)
  }
  n_waves <- length(waves_list[[1]])
  for (net in waves_list) {
    if (!is.list(net) || length(net) < 2L) {
      stop("each network must be a list of >= 2 wave matrices.", call. = FALSE)
    }
    if (length(net) != n_waves) {
      stop("all networks must have the same number of waves.", call. = FALSE)
    }
  }
  internal <- .build_multinet_data(waves_list, composition, cc_option)
  structure(list(internal = internal, n_networks = length(waves_list),
    n_waves = n_waves, n_periods = n_waves - 1L), class = "cusna_multinet_data")
}

#' @export
print.cusna_multinet_data <- function(x, ...) {
  cat(sprintf("<cusna multi-network data> %d networks, %d waves (%d period%s)\n",
              x$n_networks, x$n_waves, x$n_periods,
              if (x$n_periods == 1L) "" else "s"))
  invisible(x)
}

#' Method-of-Moments estimation of a multi-network co-evolution model
#'
#' @param data A [saom_multinet_data] object.
#' @param effects_by_net A list with one effects list per network; cross-network
#'   effects carry `net_ref` = the 1-based source network index.
#' @param backend `"cpu"` or `"gpu"`.
#' @param control A [mom_control] list.
#' @param verbose Print the parameter trace.
#' @return A `cusna_fit` object.
#' @export
mom_estimate_multinet <- function(data, effects_by_net, backend = c("cpu", "gpu"),
                                  control = mom_control(), verbose = FALSE) {
  backend <- match.arg(backend)
  if (!inherits(data, "cusna_multinet_data")) {
    stop("`data` must be created with saom_multinet_data().", call. = FALSE)
  }
  if (!is.list(effects_by_net) || length(effects_by_net) != data$n_networks) {
    stop("`effects_by_net` must be a list with one effects list per network (",
         data$n_networks, ").", call. = FALSE)
  }
  if (backend == "gpu") {
    stop("the R estimator currently drives the CPU simulator; use backend = 'cpu'.",
         call. = FALSE)
  }
  specs_by_net <- lapply(effects_by_net, function(effs) {
    if (!is.list(effs)) stop("each element of `effects_by_net` must be a list of effects.",
                             call. = FALSE)
    .build_specs(effs)
  })
  ctl <- unclass(control)
  if (!exists(".multi_mom_estimator", mode = "function")) {
    stop("multi-network estimation is not available in this build.", call. = FALSE)
  }
  est <- .multi_mom_estimator(data$internal, specs_by_net, control = ctl)
  new_cusna_fit(.multi_mom_run(est, verbose = verbose, backend = backend))
}
