# =============================================================================
# Effect registry + preprocessing -- pure-R port of the effect-side of
# prototype/cusna/saom/stats.py (track A-RAPI, milestone R2). Builds the float64
# descriptors the native C ABI consumes (type_ids/p1/p2/vmat/vmiss/dyn/v_dyn/
# dymats/dyc_idx + gwesp weight tables + interaction indices + typed F matrices),
# WITHOUT Python. Bar: descriptors == Python prepare_effects_v2 & friends.
#
# Period index convention (as in data.R): a 1-based period p in 1..P; Python's
# 0-based `period = p - 1`. Covariate/dyadic column selection clamps to the last
# available column (RSiena varCovar/varDyadCovar), matching min(period, cols-1).
# =============================================================================

# --- numeric effect ids (shared with the CUDA kernel / numba reference) ------
.EFF <- c(
  density = 0L, recip = 1L, transTrip = 2L, cycle3 = 3L, inPop = 4L,
  transTies = 5L, transRecTrip = 6L, inPopSqrt = 7L, outPop = 8L,
  outPopSqrt = 9L, egoX = 10L, altX = 11L, simX = 12L, sameX = 13L,
  inAct = 14L, inActSqrt = 15L, outAct = 16L, outActSqrt = 17L,
  outTrunc = 18L, reciAct = 19L, higher = 20L, egoXaltX = 21L, X = 22L,
  transTriads = 23L, gwesp = 24L, egoPlusAltX = 25L, gwespFB = 26L,
  gwespBF = 27L, gwespBB = 28L, gwdspFF = 29L, gwdspFB = 30L,
  altRThresholdX = 31L, interaction = 32L, crprod = 33L, crprodRecip = 34L,
  crprodMutual = 35L, from = 36L, to = 37L)

.STRUCTURAL <- c(density = 0L, recip = 1L, transTrip = 2L, cycle3 = 3L,
                 inPop = 4L, transTies = 5L, transRecTrip = 6L, inPopSqrt = 7L,
                 outPop = 8L, outPopSqrt = 9L, inAct = 14L, inActSqrt = 15L,
                 outAct = 16L, outActSqrt = 17L, outTrunc = 18L, reciAct = 19L,
                 transTriads = 23L, gwesp = 24L, gwespFF = 24L, gwespFB = 26L,
                 gwespBF = 27L, gwespBB = 28L, gwdspFF = 29L, gwdspFB = 30L)
.COVARIATE <- c(egoX = 10L, altX = 11L, simX = 12L, sameX = 13L, higher = 20L,
                egoXaltX = 21L, egoPlusAltX = 25L, altRThresholdX = 31L)
.CROSSNET <- c(crprod = 33L, crprodRecip = 34L, crprodMutual = 35L,
               from = 36L, to = 37L)
.DYADIC <- c(X = 22L)
.GW_FAMILY <- c(24L, 26L, 27L, 28L, 29L, 30L)       # gwesp/gwdsp ids
.CROSS_IDS <- c(33L, 34L, 35L, 36L, 37L)
.RATE_KINDS <- c(RateX = 0L, outRate = 1L, inRate = 2L)
.BEH_NAMES <- c(linear = 0L, quad = 1L, avAlt = 2L, avSim = 3L, totSim = 4L)
.ETYPES <- c(eval = 0L, creation = 1L, endow = 2L)
.EFF_SIMX <- 12L; .EFF_X <- 22L; .EFF_ALTRTHRX <- 31L
.EFF_EGOPLUSALTX <- 25L; .EFF_INTERACTION <- 32L
.DEG_P1 <- c(18L, 8L, 9L, 14L, 15L)                 # outTrunc/outPop*/inAct*

.DEFAULT_PARM <- c(outTrunc = 1.0, outPop = 1.0, inAct = 1.0, outPopSqrt = 0.0,
                   inActSqrt = 0.0, gwesp = 69.0, gwespFF = 69.0, gwespFB = 69.0,
                   gwespBF = 69.0, gwespBB = 69.0, gwdspFF = 69.0, gwdspFB = 69.0,
                   altRThresholdX = 0.0)
.INT_COERCE <- c("gwesp", "gwespFF", "gwespFB", "gwespBF", "gwespBB",
                 "gwdspFF", "gwdspFB", "altRThresholdX")

# RSiena GwespFunction cumulative weights: w[s]=e^a(1-(1-e^-a)^s), w[0]=0.
.gwesp_weights <- function(parm, n) {
  alpha <- 0.01 * as.numeric(parm)
  w <- numeric(n + 1L)
  factor <- 1.0 - exp(-alpha)
  pw <- 1.0
  for (s in seq_len(n)) {
    pw <- pw * factor
    w[s + 1L] <- exp(alpha) * (1.0 - pw)
  }
  w
}

# Build an internal effect spec (port of EffectSpec.__post_init__).
.make_effect_spec <- function(name, cov = NULL, parm = NULL, dyn = FALSE,
                              etype = "eval", net_ref = NULL) {
  if (!etype %in% names(.ETYPES)) stop(sprintf("unknown effect type %s", etype))
  spec <- list(name = name, cov = cov, parm = parm, dyn = isTRUE(dyn),
               etype = etype, net_ref = net_ref)
  if (name %in% names(.CROSSNET)) {
    if (is.null(net_ref)) {
      stop(sprintf("%s requires net_ref (the other dependent network index)",
                   name), call. = FALSE)
    }
    spec$type_id <- unname(.CROSSNET[[name]])
    spec$needs_cov <- FALSE
    if (is.null(spec$parm)) {
      spec$parm <- if (name %in% names(.DEFAULT_PARM)) unname(.DEFAULT_PARM[[name]]) else 1.0
    }
    return(spec)
  }
  if (name %in% names(.STRUCTURAL)) {
    spec$type_id <- unname(.STRUCTURAL[[name]])
    spec$needs_cov <- FALSE
  } else if (name %in% names(.COVARIATE)) {
    spec$type_id <- unname(.COVARIATE[[name]])
    spec$needs_cov <- !spec$dyn
    if (is.null(spec$cov) && spec$needs_cov) {
      stop(sprintf("%s requires a covariate", name), call. = FALSE)
    }
  } else if (name %in% names(.DYADIC)) {
    spec$type_id <- unname(.DYADIC[[name]])
    spec$needs_cov <- TRUE
    w <- if (is.null(spec$cov)) NULL else spec$cov
    d <- if (is.null(w)) NULL else dim(as.array(w))
    if (is.null(d) || !(length(d) %in% c(2L, 3L)) || d[1] != d[2]) {
      stop(sprintf("%s requires an (n,n) or (n,n,periods) dyadic covariate",
                   name), call. = FALSE)
    }
  } else {
    stop(sprintf("unknown effect %s", name), call. = FALSE)
  }
  if (is.null(spec$parm)) {
    spec$parm <- if (name %in% names(.DEFAULT_PARM)) unname(.DEFAULT_PARM[[name]]) else 1.0
  }
  # RSiena stores parm in an INTEGER column -> non-integral parameters are
  # truncated before reaching the C++ effect (gwesp(0.69)->0, threshold 0.5->0).
  if (name %in% .INT_COERCE) spec$parm <- as.numeric(trunc(spec$parm))
  spec
}

.make_beh_effect_spec <- function(name) {
  if (!name %in% names(.BEH_NAMES)) {
    stop(sprintf("unknown behavior effect %s", name), call. = FALSE)
  }
  list(name = name, type_id = unname(.BEH_NAMES[[name]]))
}

.make_rate_effect_spec <- function(name, cov = NULL) {
  if (!name %in% names(.RATE_KINDS)) {
    stop(sprintf("unknown rate effect %s", name), call. = FALSE)
  }
  kind <- unname(.RATE_KINDS[[name]])
  if (kind == 0L && is.null(cov)) stop("RateX requires a covariate", call. = FALSE)
  list(name = name, cov = cov, kind = kind)
}

.make_interaction_spec <- function(comps) {
  c0 <- as.integer(comps)
  if (!length(c0) %in% c(2L, 3L)) stop("interaction needs 2 or 3 components")
  list(name = "unspInt", comps = c0, type_id = .EFF_INTERACTION,
       needs_cov = FALSE, cov = NULL, parm = 1.0, dyn = FALSE, etype = "eval")
}

.is_interaction <- function(e) identical(e$type_id, .EFF_INTERACTION)

# Centered actor-covariate values active in period p (grand mean over all cols).
.effect_cov_values <- function(spec, p = 1L) {
  v <- spec$cov
  gmean <- mean(as.numeric(v), na.rm = TRUE)
  if (is.null(dim(v)) || length(dim(v)) == 1L) {
    return(as.numeric(v) - gmean)
  }
  col <- min(p, ncol(v))
  v[, col] - gmean
}

# (1/range, simMean) for simX: all covariate columns pooled, range over all.
.effect_sim_params <- function(spec) {
  v <- spec$cov
  rng <- max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  cols <- if (is.null(dim(v))) matrix(as.numeric(v), ncol = 1L) else
    matrix(as.numeric(v), nrow = nrow(v))
  tot <- 0.0
  cnt <- 0L
  for (m in seq_len(ncol(cols))) {
    x <- cols[!is.na(cols[, m]), m]
    if (length(x) > 1L) {
      d <- 1.0 - abs(outer(x, x, "-")) / rng
      tot <- tot + (sum(d) - length(x))
      cnt <- cnt + length(x) * (length(x) - 1L)
    }
  }
  c(1.0 / rng, if (cnt > 0L) tot / cnt else 0.0)
}

.rate_cov_values <- function(spec, p = 1L) {
  v <- spec$cov
  gmean <- mean(as.numeric(v), na.rm = TRUE)
  vm <- if (is.null(dim(v))) as.numeric(v) - gmean else v[, min(p, ncol(v))] - gmean
  ifelse(is.na(vm), 0.0, vm)
}

# Port of prepare_effects_v2: returns the float64 descriptor bundle.
.prepare_effects_v2 <- function(effects, n, p = 1L, beh = NULL, z_start = NULL) {
  k <- length(effects)
  type_ids <- vapply(effects, function(e) as.integer(e$type_id), integer(1))
  p1 <- numeric(k)
  p2 <- numeric(k)
  vmat <- matrix(0.0, k, n)
  vmiss <- matrix(0L, k, n)
  dyn_flags <- vapply(effects, function(e) if (isTRUE(e$dyn)) 1L else 0L, integer(1))
  dymats <- list()
  dyc_idx <- rep(-1L, k)
  off <- !diag(TRUE, n)
  for (idx in seq_len(k)) {
    e <- effects[[idx]]
    t <- as.integer(e$type_id)
    if (t == .EFF_SIMX) {
      if (isTRUE(e$dyn)) {
        if (!is.null(beh)) {
          p1[idx] <- 1.0 / beh$range
          p2[idx] <- beh$sim_mean
        }
      } else {
        sp <- .effect_sim_params(e)
        p1[idx] <- sp[1]; p2[idx] <- sp[2]
      }
    } else if (t %in% .DEG_P1) {
      p1[idx] <- as.numeric(e$parm)
    } else if (t == .EFF_ALTRTHRX) {
      p1[idx] <- as.numeric(trunc(e$parm))
    } else if (t %in% .CROSS_IDS) {
      p1[idx] <- as.numeric(e$parm)
    }
    if (t == .EFF_X) {
      w <- e$cov
      d <- dim(as.array(w))
      if (length(d) == 3L) {
        gmean <- mean(matrix(as.numeric(w), nrow = n * n)[as.vector(off), ],
                      na.rm = TRUE)
        wm <- w[, , min(p, d[3])] - gmean
      } else {
        wm <- w - mean(w[off])
      }
      wm <- as.matrix(wm)
      diag(wm) <- 0.0
      dymats[[length(dymats) + 1L]] <- wm
      dyc_idx[idx] <- length(dymats) - 1L         # 0-based index into dymats
    } else if (isTRUE(e$needs_cov) && !isTRUE(e$dyn)) {
      vv <- .effect_cov_values(e, p)
      if (length(vv) != n) stop(sprintf("%s: covariate length != n", e$name))
      miss <- is.na(vv)
      vmiss[idx, ] <- as.integer(miss)
      vmat[idx, ] <- ifelse(miss, 0.0, vv)
    } else if (isTRUE(e$dyn) && !is.null(beh)) {
      vmiss[idx, ] <- as.integer(beh$miss[, min(p, ncol(beh$miss))])
    }
  }
  if (any(dyn_flags == 1L) && !is.null(beh) && !is.null(z_start)) {
    v_dyn <- as.numeric(z_start) - beh$overall_mean
  } else {
    v_dyn <- numeric(n)
  }
  list(type_ids = type_ids, p1 = p1, p2 = p2, vmat = vmat, dyn_flags = dyn_flags,
       v_dyn = v_dyn, dymats = dymats, dyc_idx = dyc_idx, vmiss = vmiss)
}

# (k, 3) 0-based component indices for interactions (-1 padded), validated.
.interaction_comp_idx <- function(effects) {
  k <- length(effects)
  out <- matrix(-1L, k, 3L)
  for (e_i in seq_len(k)) {
    e <- effects[[e_i]]
    if (!.is_interaction(e)) next
    for (slot in seq_along(e$comps)) {
      ci <- e$comps[slot]                          # 0-based
      if (!(ci >= 0L && ci < (e_i - 1L))) {
        stop(sprintf("interaction component %d must precede position %d",
                     ci, e_i - 1L), call. = FALSE)
      }
      comp <- effects[[ci + 1L]]
      if (.is_interaction(comp)) stop("nested interactions not supported")
      if (!identical(comp$etype, "eval")) {
        stop("interaction components must be eval effects", call. = FALSE)
      }
      out[e_i, slot] <- as.integer(ci)
    }
  }
  out
}

# (k) 0 eval / 1 creation / 2 endowment (creation/endow only for egoPlusAltX).
.effect_etypes <- function(effects) {
  vapply(effects, function(e) {
    et <- unname(.ETYPES[[if (is.null(e$etype)) "eval" else e$etype]])
    if (et != 0L && !identical(as.integer(e$type_id), .EFF_EGOPLUSALTX)) {
      stop("creation/endowment types are supported for egoPlusAltX only",
           call. = FALSE)
    }
    et
  }, integer(1))
}

# (k) source-network index for cross-network effects (-1 otherwise).
.effect_netrefs <- function(effects) {
  vapply(effects, function(e) {
    r <- e$net_ref
    if (is.null(r)) -1L else as.integer(r)
  }, integer(1))
}

# (k, n+1) gwesp/gwdsp cumulative weight tables (float64; zero rows otherwise).
.effect_wtabs <- function(effects, n) {
  k <- length(effects)
  tabs <- matrix(0.0, k, n + 1L)
  for (i in seq_len(k)) {
    if (as.integer(effects[[i]]$type_id) %in% .GW_FAMILY) {
      tabs[i, ] <- .gwesp_weights(effects[[i]]$parm, n)
    }
  }
  tabs
}

# (n, n) f_ij pair matrix for a typed egoPlusAltX creation/endowment effect.
.typed_pair_matrix <- function(spec, n, p = 1L) {
  if (as.integer(spec$type_id) != .EFF_EGOPLUSALTX) {
    stop("typed_pair_matrix supports egoPlusAltX only", call. = FALSE)
  }
  v <- .effect_cov_values(spec, p)
  cov <- spec$cov
  miss <- is.na(as.numeric(cov))
  if (!is.null(dim(cov)) && length(dim(cov)) == 2L) {
    miss <- is.na(cov[, min(p, ncol(cov))])
  }
  vv <- ifelse(miss, 0.0, v)
  ok2 <- !(outer(miss, miss, "|"))
  f <- outer(vv, vv, "+") * ok2
  diag(f) <- 0.0
  f
}

# ============================================================================
# Public API (milestone R5): model-specification constructors. Lightweight R
# objects (no Python); mom_estimate() translates them into internal specs.
# ============================================================================

#' Network evaluation effect
#'
#' @param name Effect name (e.g. "density", "recip", "transTrip", "egoX",
#'   "gwesp", a cross-network "crprod"/"from"/...).
#' @param covariate Covariate values for covariate effects: `(n)` constant or
#'   `(n, n_periods)` changing (actor); `(n, n)` or `(n, n, periods)` (dyadic).
#' @param parameter Internal effect parameter (RSiena `parameter`; `NULL` =
#'   default).
#' @param dyn `TRUE` if the covariate is the co-evolving dependent behavior.
#' @param type `"eval"`, `"creation"` or `"endow"`.
#' @param net_ref For cross-network effects: the 1-based index of the source
#'   dependent network.
#' @return A `cusna_effect` object.
#' @export
cusna_effect <- function(name, covariate = NULL, parameter = NULL, dyn = FALSE,
                         type = c("eval", "creation", "endow"), net_ref = NULL) {
  type <- match.arg(type)
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single effect name (a string).", call. = FALSE)
  }
  structure(list(name = name, covariate = covariate, parameter = parameter,
    dyn = isTRUE(dyn), type = type, net_ref = net_ref), class = "cusna_effect")
}

#' @export
print.cusna_effect <- function(x, ...) {
  cat("<cusna effect>", x$name)
  if (!identical(x$type, "eval")) cat(" [", x$type, "]", sep = "")
  if (isTRUE(x$dyn)) cat(" (behavior-dynamic)")
  if (!is.null(x$covariate)) cat(" + covariate")
  if (!is.null(x$parameter)) cat("  parm =", x$parameter)
  if (!is.null(x$net_ref)) cat("  <- net", x$net_ref)
  cat("\n")
  invisible(x)
}

#' Behavior evaluation effect
#' @param name One of "linear", "quad", "avAlt", "avSim", "totSim".
#' @return A `cusna_beh_effect` object.
#' @export
cusna_beh_effect <- function(name) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single behavior-effect name (a string).", call. = FALSE)
  }
  structure(list(name = name), class = "cusna_beh_effect")
}

#' @export
print.cusna_beh_effect <- function(x, ...) {
  cat("<cusna behavior effect>", x$name, "\n"); invisible(x)
}

#' Non-basic rate effect
#' @param name One of "RateX", "outRate", "inRate".
#' @param covariate Covariate for "RateX".
#' @return A `cusna_rate_effect` object.
#' @export
cusna_rate_effect <- function(name, covariate = NULL) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be one of \"RateX\", \"outRate\", \"inRate\".", call. = FALSE)
  }
  structure(list(name = name, covariate = covariate), class = "cusna_rate_effect")
}

#' @export
print.cusna_rate_effect <- function(x, ...) {
  cat("<cusna rate effect>", x$name)
  if (!is.null(x$covariate)) cat(" + covariate")
  cat("\n"); invisible(x)
}

#' User interaction effect
#'
#' The product of 2 or 3 component effects that appear earlier in the effects
#' list. `components` are their 1-based positions in that list.
#' @param components Integer vector of 2 or 3 positions.
#' @return A `cusna_interaction` object.
#' @export
cusna_interaction <- function(components) {
  components <- as.integer(components)
  if (length(components) < 2L || length(components) > 3L || anyNA(components) ||
      any(components < 1L)) {
    stop("`components` must be 2 or 3 positive positions in the effects list.",
         call. = FALSE)
  }
  structure(list(components = components), class = "cusna_interaction")
}

#' @export
print.cusna_interaction <- function(x, ...) {
  cat("<cusna interaction> components", paste(x$components, collapse = " x "), "\n")
  invisible(x)
}
