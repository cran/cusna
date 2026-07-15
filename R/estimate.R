# =============================================================================
# Method-of-Moments estimation -- the user-facing API (track A-RAPI, milestone
# R5). The interface mirrors the reticulate package (Rpkg/R/estimate.R); the
# bodies drive the native engine (R1-R4) instead of Python: model specs are
# translated to internal effect specs, the Robbins-Monro loop runs on the
# compiled C ABI simulator, and the result is an R-native cusna_fit object.
# =============================================================================

#' Control settings for the Robbins-Monro estimator
#'
#' @param firstg Initial gain of the Robbins-Monro updates.
#' @param nsub Number of phase-2 subphases.
#' @param n2 Integer vector (length `nsub`): phase-2 iterations per subphase.
#' @param batch2 Chains per phase-2 iteration.
#' @param n1,n3 Chains for the phase-1 derivative and the phase-3 estimate.
#' @param nD Chains for the phase-2 subphase derivatives.
#' @param diagonalize Weight in [0,1] mixing the full and diagonal derivative.
#' @param seed Random seed for the simulation stream.
#' @return A `cusna_control` list.
#' @export
mom_control <- function(firstg = 0.4, nsub = 4, n2 = c(40, 40, 60, 80),
                        batch2 = 256, n1 = 1000, n3 = 4000, nD = 1000,
                        diagonalize = 0.2, seed = 1234) {
  structure(list(firstg = firstg, nsub = as.integer(nsub), n2 = as.numeric(n2),
    batch2 = as.integer(batch2), n1 = as.integer(n1), n3 = as.integer(n3),
    nD = as.integer(nD), diagonalize = diagonalize, seed = as.integer(seed)),
    class = "cusna_control")
}

# Translate the R model specs into internal effect specs. Interaction component
# positions and cross-network net_ref are 1-based in R, 0-based in the engine.
.build_specs <- function(effects) {
  lapply(effects, function(e) {
    if (inherits(e, "cusna_interaction")) {
      .make_interaction_spec(as.integer(e$components) - 1L)
    } else if (inherits(e, "cusna_effect")) {
      .make_effect_spec(name = e$name, cov = e$covariate, parm = e$parameter,
        dyn = isTRUE(e$dyn), etype = e$type,
        net_ref = if (is.null(e$net_ref)) NULL else as.integer(e$net_ref) - 1L)
    } else {
      stop("`effects` must contain cusna_effect / cusna_interaction objects.",
           call. = FALSE)
    }
  })
}

#' Method-of-Moments estimation of a SAOM
#'
#' Fits a stochastic actor-oriented model by the RSiena-style Method of Moments,
#' driving the native engine (no Python). Returns a [cusna_fit] object.
#'
#' @param data A [saom_data] object.
#' @param effects List of [cusna_effect] / [cusna_interaction] specifications.
#' @param beh_effects Optional list of [cusna_beh_effect] (co-evolution).
#' @param rate_effects Optional list of [cusna_rate_effect] (non-basic rates).
#' @param conditional Condition on the observed distances (RSiena default for a
#'   single dependent variable). Unconditional adds the basic rates to the model.
#' @param backend `"cpu"` (the from-source CRAN build) or `"gpu"` (a CUDA build).
#' @param maxdegree Optional out-degree bound (0 = none).
#' @param control A [mom_control] list.
#' @param verbose Print the parameter trace per phase.
#' @return A `cusna_fit` object.
#' @export
mom_estimate <- function(data, effects, beh_effects = NULL, rate_effects = NULL,
                         conditional = FALSE, backend = c("cpu", "gpu"),
                         maxdegree = 0, control = mom_control(), verbose = FALSE) {
  backend <- match.arg(backend)
  if (!inherits(data, "cusna_data")) {
    stop("`data` must be created with saom_data().", call. = FALSE)
  }
  if (!is.list(effects) || !length(effects)) {
    stop("`effects` must be a non-empty list of cusna_effect objects.", call. = FALSE)
  }
  if (backend == "gpu" && !cusna_has_cuda()) {
    stop("the 'gpu' backend needs a CUDA build of cusna; use backend = 'cpu'.",
         call. = FALSE)
  }
  if (backend == "gpu") {
    stop("the R estimator currently drives the CPU simulator; use backend = 'cpu'.",
         call. = FALSE)
  }
  specs <- .build_specs(effects)
  bspecs <- if (!is.null(beh_effects))
    lapply(beh_effects, function(e) .make_beh_effect_spec(e$name)) else NULL
  rspecs <- if (!is.null(rate_effects))
    lapply(rate_effects, function(e) .make_rate_effect_spec(e$name, e$covariate)) else NULL
  ctl <- unclass(control); ctl$maxdegree <- as.integer(maxdegree)
  est <- .mom_estimator(data$internal, specs, beh_effects = bspecs,
                        rate_effects = rspecs, conditional = conditional, control = ctl)
  new_cusna_fit(.mom_run(est, verbose = verbose, backend = backend))
}

# Assemble the R-native fit object from the engine result.
new_cusna_fit <- function(res) {
  nm <- as.character(res$param_names)
  theta <- as.numeric(res$theta); se <- as.numeric(res$se); tconv <- as.numeric(res$tconv)
  names(theta) <- names(se) <- names(tconv) <- nm
  estimates <- data.frame(effect = nm, estimate = theta, se = se, t_conv = tconv,
                          stringsAsFactors = FALSE, row.names = NULL)
  rates <- NULL
  if (isTRUE(res$conditional) && !is.null(res$rate_est)) {
    rates <- data.frame(period = seq_along(res$rate_est), rate = as.numeric(res$rate_est),
                        sd_duration = as.numeric(res$rate_sd))
  }
  structure(list(coefficients = theta, se = se, tconv = tconv, estimates = estimates,
    rates = rates, tconv_max = as.numeric(res$tconv_max), n_sims = as.integer(res$n_sims),
    wall_time = as.numeric(res$wall_time), backend = as.character(res$backend),
    conditional = isTRUE(res$conditional), targets = as.numeric(res$targets),
    sim_means = as.numeric(res$sim_means), D = res$D, msf = res$msf),
    class = "cusna_fit")
}

#' @export
print.cusna_fit <- function(x, ...) {
  cat(sprintf("cusna MoM fit  (%s backend, %s)\n", x$backend,
              if (x$conditional) "conditional" else "unconditional"))
  est <- x$estimates
  est$estimate <- round(est$estimate, 4); est$se <- round(est$se, 4)
  est$t_conv <- round(est$t_conv, 3)
  print(est, row.names = FALSE)
  cat(sprintf("overall max. convergence ratio: %.4f\n%d simulations, %.2fs\n",
              x$tconv_max, x$n_sims, x$wall_time))
  invisible(x)
}

#' @export
summary.cusna_fit <- function(object, ...) {
  lines <- character(0)
  if (object$conditional && !is.null(object$rates)) {
    for (i in seq_len(nrow(object$rates))) {
      lines <- c(lines, sprintf("%-16s%10.4f  (sd of durations %.4f)",
        paste0("rate period ", object$rates$period[i]),
        object$rates$rate[i], object$rates$sd_duration[i]))
    }
  }
  lines <- c(lines, sprintf("%-16s%10s%9s%9s", "parameter", "estimate", "s.e.", "t-conv"))
  for (i in seq_along(object$coefficients)) {
    lines <- c(lines, sprintf("%-16s%10.4f%9.4f%9.3f", names(object$coefficients)[i],
      object$coefficients[i], object$se[i], object$tconv[i]))
  }
  lines <- c(lines, sprintf("overall max. convergence ratio: %.4f", object$tconv_max),
    sprintf("total simulations: %d, wall time: %.2fs, backend: %s, conditional: %s",
      object$n_sims, object$wall_time, object$backend, object$conditional))
  cat(paste(lines, collapse = "\n"), "\n")
  invisible(object)
}

#' @export
coef.cusna_fit <- function(object, ...) object$coefficients

# Delta-method covariance: Dinv %*% MSF %*% t(Dinv), the construction the engine
# uses for the reported standard errors.
#' @export
vcov.cusna_fit <- function(object, ...) {
  dinv <- .safe_solve(as.matrix(object$D))
  V <- dinv %*% as.matrix(object$msf) %*% t(dinv)
  nm <- names(object$coefficients)
  dimnames(V) <- list(nm, nm)
  V
}

#' @export
as.data.frame.cusna_fit <- function(x, ...) x$estimates
