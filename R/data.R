# =============================================================================
# SAOM data conventions carrier -- pure-R port of prototype/cusna/saom/data.py
# (track A-RAPI, milestone R1). This is the host-side "data engine": it applies
# RSiena's imputation and builds every per-period mask/target the native engine
# and the moment targets need, WITHOUT Python. The public saom_data() surface
# (milestone R5) is a thin wrapper over the internal builder below.
#
# Bit-identity bar (R1): the integer masks / imputed states / rate targets are
# bit-identical to data.py on s50, s50miss, s50str, Knecht, COW, synthetic;
# float64 behavior scalars (overallMean/range/simMean/valuesLessMissingStarts)
# match to machine precision (sum-order differences only).
#
# Convention notes mirror the Python docstrings; see data.py for the RSiena
# source references (initializeFRAN.r, StatisticCalculator, NetworkVariable.cpp).
#
# Period index convention: methods take a 1-based period `p` in 1..P; period p
# spans waves p (start) and p+1 (end). This maps to Python's 0-based m = p-1.
# =============================================================================

.NA_CODE <- 9L
.STRUCT_ZERO <- 10L
.STRUCT_ONE <- 11L

# Split one raw wave into (vals 0/1 uint8-like, miss logical, struct logical).
# Structural codes contribute their value (10 -> 0, 11 -> 1) and are flagged;
# code 9 (integer input) or NaN (floating input) are missing.
.saom_parse_wave <- function(w) {
  w <- as.matrix(w)
  n <- nrow(w)
  if (n != ncol(w)) stop("each wave must be a square matrix.", call. = FALSE)
  if (is.double(w) && any(is.nan(w))) {
    # floating input: NaN = missing, values truncated to int (numpy astype)
    miss <- is.nan(w)
    wi <- w
    wi[miss] <- 0
    wi <- matrix(as.integer(trunc(wi)), n, n)
  } else {
    wi <- matrix(as.integer(round(w)), n, n)
    miss <- wi == .NA_CODE
  }
  struct <- (wi == .STRUCT_ZERO) | (wi == .STRUCT_ONE)
  vals <- wi
  vals[wi == .STRUCT_ZERO] <- 0L
  vals[wi == .STRUCT_ONE] <- 1L
  vals[miss] <- 0L
  chk <- vals[!miss & !struct]
  if (!all(chk %in% c(0L, 1L))) {
    stop("tie values must be 0/1/9(NA)/10/11 or NaN", call. = FALSE)
  }
  storage.mode(vals) <- "integer"
  diag(vals) <- 0L
  diag(miss) <- FALSE
  diag(struct) <- FALSE
  list(vals = vals, miss = miss, struct = struct)
}

# Build the internal SAOM data object (list) from raw integer-coded waves,
# optional behavior (n x n_waves, NaN = missing) and optional composition
# (per-actor list of c(start, end) intervals, 1-based wave numbers).
.build_saom_data <- function(waves, behavior = NULL, composition = NULL,
                             cc_option = 1) {
  if (!is.list(waves) || length(waves) < 2L) {
    stop("`waves` must be a list of at least two adjacency matrices.",
         call. = FALSE)
  }
  n <- nrow(as.matrix(waves[[1]]))
  nw <- length(waves)
  parsed <- lapply(waves, .saom_parse_wave)
  for (p in parsed) {
    if (nrow(p$vals) != n) {
      stop("all waves must have the same node set.", call. = FALSE)
    }
  }
  vals <- lapply(parsed, `[[`, "vals")
  miss <- lapply(parsed, `[[`, "miss")
  struct <- lapply(parsed, `[[`, "struct")

  # RSiena imputation: wave 1 NA -> 0 (vals already 0), later waves NA -> the
  # previous imputed value (recursive carry-forward).
  imputed <- vector("list", nw)
  imputed[[1]] <- vals[[1]]
  for (m in 2:nw) {
    w <- vals[[m]]
    mm <- miss[[m]]
    w[mm] <- imputed[[m - 1]][mm]
    imputed[[m]] <- w
  }

  # symmetric panel (identical(m, t(m)) incl. NA/struct pattern, every wave)
  symmetric <- all(vapply(seq_len(nw), function(m) {
    isSymmetric(unname(vals[[m]])) &&
      all(miss[[m]] == t(miss[[m]])) &&
      all(struct[[m]] == t(struct[[m]]))
  }, logical(1)))

  obj <- list(
    n = n, n_waves = nw, n_periods = nw - 1L,
    vals = vals, miss = miss, struct = struct,
    imputed = imputed, rawvals = vals,        # rawvals == vals pre-composition
    symmetric = symmetric,
    has_composition = FALSE, comp_active = NULL, comp_intervals = NULL,
    cc_option = as.integer(cc_option),
    beh = NULL
  )

  if (!is.null(composition)) {
    obj <- .saom_apply_composition(obj, composition, cc_option)
  }
  if (!is.null(behavior)) {
    if (!is.null(composition)) {
      stop("composition change with a dependent behavior is not supported.",
           call. = FALSE)
    }
    obj$beh <- .build_behavior_data(behavior, n, nw)
  }
  class(obj) <- "cusna_saom_internal"
  obj
}

# Replicate RSiena's composition-change data treatment (sienaDataCreate.r +
# initializeFRAN.r unpackOneMode, integer-time events only -> pure activeStart
# gating). Modifies miss/rawvals/imputed and sets comp_active (P x n).
.saom_apply_composition <- function(obj, composition, cc_option) {
  n <- obj$n
  nw <- obj$n_waves
  cc <- as.integer(cc_option)
  if (!cc %in% c(1L, 2L, 3L)) stop("ccOption must be 1, 2 or 3.", call. = FALSE)
  if (length(composition) != n) {
    stop("composition needs one interval list per actor.", call. = FALSE)
  }

  intervals <- vector("list", n)
  active_start <- matrix(FALSE, n, nw)
  for (j in seq_len(n)) {
    flat <- as.numeric(unlist(composition[[j]]))
    if (length(flat) < 2L || length(flat) %% 2L != 0L) {
      stop(sprintf("actor %d: interval list must have an even length.", j),
           call. = FALSE)
    }
    if (!all(flat == round(flat))) {
      stop(sprintf("actor %d: only integer event times are supported.", j),
           call. = FALSE)
    }
    pairs <- vector("list", length(flat) %/% 2L)
    for (q in seq_len(length(flat) %/% 2L)) {
      s <- as.integer(flat[2L * q - 1L])
      e <- as.integer(flat[2L * q])
      if (!(1L <= s && s <= e && e <= nw)) {
        stop(sprintf("actor %d: interval (%d,%d) outside 1..%d", j, s, e, nw),
             call. = FALSE)
      }
      pairs[[q]] <- c(s, e)
      active_start[j, s:e] <- TRUE
    }
    intervals[[j]] <- pairs
  }

  # action codes for waves where the actor is not active at the start:
  # 1 before joining, 2 after leaving, 3 temporary gap.
  action <- matrix(0L, n, nw)
  for (j in seq_len(n)) {
    for (m in seq_len(nw)) {
      if (active_start[j, m]) next
      prec <- if (m > 1L) any(active_start[j, seq_len(m - 1L)]) else FALSE
      succ <- if (m < nw) any(active_start[j, (m + 1L):nw]) else FALSE
      action[j, m] <- if (!prec) 1L else if (!succ) 2L else 3L
    }
  }

  # data mods, per wave ascending; within a wave in the R order 1, 3, 2. The NA
  # pattern comes from the actor's ROW and is applied to row AND column.
  for (m in seq_len(nw)) {
    for (code in c(1L, 3L, 2L)) {
      js <- which(action[, m] == code)
      for (j in js) {
        row_na <- obj$miss[[m]][j, ]           # logical length n (copied)
        zero_branch <- (code == 1L && cc %in% c(1L, 2L)) ||
          (code == 3L && cc %in% c(1L, 2L)) || (code == 2L && cc == 1L)
        if (zero_branch) {
          obj$miss[[m]][j, row_na] <- FALSE
          obj$miss[[m]][row_na, j] <- FALSE
          obj$rawvals[[m]][j, row_na] <- 0L
          obj$rawvals[[m]][row_na, j] <- 0L
          if (code == 1L) {
            obj$imputed[[m]][j, row_na] <- 0L
            obj$imputed[[m]][row_na, j] <- 0L
          } else if (m == 1L) {
            obj$imputed[[m]][j, row_na] <- 0L
            obj$imputed[[m]][row_na, j] <- 0L
          } else {
            obj$imputed[[m]][j, row_na] <- obj$imputed[[m - 1L]][j, row_na]
            obj$imputed[[m]][row_na, j] <- obj$imputed[[m - 1L]][row_na, j]
          }
        } else {
          # whole row+column missing (keep diagonal non-missing)
          obj$miss[[m]][j, ] <- TRUE
          obj$miss[[m]][, j] <- TRUE
          obj$miss[[m]][j, j] <- FALSE
        }
      }
    }
  }

  P <- nw - 1L
  obj$has_composition <- TRUE
  obj$comp_intervals <- intervals
  # active DURING period m == active at the period's START wave
  obj$comp_active <- t(active_start[, seq_len(P), drop = FALSE])   # P x n
  obj
}

# Dependent behavior variable (n, n_waves), integer values; NaN = missing.
.build_behavior_data <- function(values, n, nw) {
  v <- matrix(as.numeric(as.matrix(values)), n, nw)
  miss <- is.na(v) | is.nan(v)
  obs <- v[!miss]
  # overallMean: mean of PER-WAVE means of observed values
  wave_means <- vapply(seq_len(nw), function(m) mean(v[!miss[, m], m]),
                       double(1))
  overall_mean <- mean(wave_means)
  zmin <- as.integer(min(obs))
  zmax <- as.integer(max(obs))
  range <- as.numeric(zmax - zmin)
  # simMean: pooled over ordered pairs within each wave EXCEPT the last;
  # range taken over all waves.
  tot <- 0.0
  cnt <- 0L
  if (nw >= 2L) {
    for (m in seq_len(nw - 1L)) {
      col <- v[, m]
      x <- col[!is.na(col)]
      if (length(x) > 1L) {
        d <- 1.0 - abs(outer(x, x, "-")) / range
        tot <- tot + (sum(d) - length(x))
        cnt <- cnt + length(x) * (length(x) - 1L)
      }
    }
  }
  sim_mean <- if (cnt > 0L) tot / cnt else 0.0

  # imputation: all-NA row -> mode; wave-1 backward fill; later carry-forward.
  imp <- v
  allna <- apply(miss, 1L, all)
  if (any(allna)) {
    tab <- table(obs)
    mode_val <- as.numeric(names(tab))[which.max(tab)]
    imp[allna, ] <- mode_val
  }
  first <- imp[, 1]
  for (m in 2:nw) first <- ifelse(is.na(first), imp[, m], first)
  imp[, 1] <- first
  for (m in 2:nw) imp[, m] <- ifelse(is.na(imp[, m]), imp[, m - 1L], imp[, m])
  if (anyNA(imp)) stop("behavior imputation failed.", call. = FALSE)
  imputed <- matrix(as.integer(round(imp)), n, nw)

  list(n = n, n_waves = nw, values = v, miss = miss, imputed = imputed,
       overall_mean = overall_mean, zmin = zmin, zmax = zmax,
       range = range, sim_mean = sim_mean)
}

# ------------------------------------------------------------ period-m masks
# All take the internal object and a 1-based period p (1..P).

.saom_kmask_dist <- function(obj, p) {
  k <- !(obj$miss[[p]] | obj$miss[[p + 1L]] |
           obj$struct[[p]] | obj$struct[[p + 1L]])
  diag(k) <- FALSE
  k
}

.saom_observed_distance <- function(obj, p) {
  excl <- obj$miss[[p]] | obj$miss[[p + 1L]] |
    obj$struct[[p]] | obj$struct[[p + 1L]]
  diff <- (obj$rawvals[[p]] != obj$rawvals[[p + 1L]]) & !excl
  diag(diff) <- FALSE
  sum(diff)
}

.saom_rate_target <- function(obj, p) {
  a <- (obj$imputed[[p]] != 0L) & !obj$miss[[p]]
  b <- (obj$imputed[[p + 1L]] != 0L) & !obj$miss[[p]] & !obj$miss[[p + 1L]]
  diff <- xor(a, b) & !obj$struct[[p]] & !obj$struct[[p + 1L]]
  diag(diff) <- FALSE
  sum(diff)
}

.saom_rate_const <- function(obj, p) {
  cc <- (obj$imputed[[p]] != 0L) & !obj$miss[[p]] & obj$miss[[p + 1L]] &
    !obj$struct[[p]] & !obj$struct[[p + 1L]]
  diag(cc) <- FALSE
  sum(cc)
}

.saom_rate_target_rows <- function(obj, p) {
  a <- (obj$imputed[[p]] != 0L) & !obj$miss[[p]]
  b <- (obj$imputed[[p + 1L]] != 0L) & !obj$miss[[p]] & !obj$miss[[p + 1L]]
  diff <- xor(a, b) & !obj$struct[[p]] & !obj$struct[[p + 1L]]
  diag(diff) <- FALSE
  as.integer(rowSums(diff))
}

.saom_rate_const_rows <- function(obj, p) {
  cc <- (obj$imputed[[p]] != 0L) & !obj$miss[[p]] & obj$miss[[p + 1L]] &
    !obj$struct[[p]] & !obj$struct[[p + 1L]]
  diag(cc) <- FALSE
  as.integer(rowSums(cc))
}

.saom_network_less_missing <- function(obj, p) {
  ((obj$imputed[[p]] != 0L) & !obj$miss[[p]]) * 1L
}

.saom_network_less_missing_start <- function(obj, p) {
  ((obj$imputed[[p]] != 0L) & !obj$miss[[p]] & !obj$miss[[p + 1L]]) * 1L
}

# (keep, force) so the statistic network is (A & keep) | force: subtract ties
# missing at either end, then write back structural values of wave p+1, then p.
.saom_stat_keep_force <- function(obj, p) {
  keep <- !(obj$miss[[p]] | obj$miss[[p + 1L]])
  force <- matrix(FALSE, obj$n, obj$n)
  s1 <- obj$struct[[p + 1L]]
  keep <- keep & !s1
  force[s1] <- (obj$imputed[[p + 1L]][s1] != 0L)
  s0 <- obj$struct[[p]]
  keep <- keep & !s0
  force[s0] <- (obj$imputed[[p]][s0] != 0L)
  diag(keep) <- FALSE
  diag(force) <- FALSE
  list(keep = keep, force = force)
}

# (keepS, forceS): ONLY the structural replacements of period p (no missing
# subtraction) -- the current* network of the endowment statistics.
.saom_stat_keep_force_struct <- function(obj, p) {
  keep <- matrix(TRUE, obj$n, obj$n)
  force <- matrix(FALSE, obj$n, obj$n)
  s1 <- obj$struct[[p + 1L]]
  keep <- keep & !s1
  force[s1] <- (obj$imputed[[p + 1L]][s1] != 0L)
  s0 <- obj$struct[[p]]
  keep <- keep & !s0
  force[s0] <- (obj$imputed[[p]][s0] != 0L)
  diag(keep) <- FALSE
  diag(force) <- FALSE
  list(keep = keep, force = force)
}

.saom_stat_network <- function(obj, a, p) {
  kf <- .saom_stat_keep_force(obj, p)
  (((a != 0L) & kf$keep) | kf$force) * 1L
}

.saom_target_network <- function(obj, p) {
  .saom_stat_network(obj, obj$imputed[[p + 1L]], p)
}

.saom_permitted <- function(obj, p) {
  pmat <- !obj$struct[[p]]
  diag(pmat) <- FALSE
  pmat
}

# ------------------------------------------------------------ behavior period masks

.beh_kmask_dist <- function(beh, p) {
  !(beh$miss[, p] | beh$miss[, p + 1L])
}

.beh_observed_distance <- function(beh, p) {
  d <- abs(beh$imputed[, p + 1L] - beh$imputed[, p])
  sum(d[.beh_kmask_dist(beh, p)])
}

.beh_values_less_missing_starts <- function(beh, p) {
  ifelse(beh$miss[, p], 0.0, as.numeric(beh$imputed[, p]))
}

.beh_stat_values <- function(beh, z_end, p) {
  u <- as.numeric(z_end) - beh$overall_mean
  u[!.beh_kmask_dist(beh, p)] <- 0.0
  u
}

# ============================================================================
# Public API (milestone R5): the user-facing SAOM data container. Wraps the
# internal data engine above; carries it in `$internal` for the estimator.
# ============================================================================

#' Longitudinal SAOM data
#'
#' Builds a SAOM data container from a panel of adjacency matrices, applying
#' RSiena's imputation and per-period masks (no Python). Missing ties are `NA`
#' or the code 9; structural zeros/ones are 10/11, read from the wave matrices.
#'
#' @param waves A list of at least two square adjacency matrices (the panel).
#' @param behavior Optional `n * n_waves` matrix of a co-evolving dependent
#'   behavior (`NA` for missing).
#' @param composition Optional per-actor composition-change intervals: a list of
#'   length `n`, each element a numeric vector `c(start, end)` (or several
#'   stacked pairs) of 1-based wave numbers the actor is present.
#' @param cc_option The sienaCompositionChange option (1, 2 or 3).
#' @return A `cusna_data` object.
#' @export
saom_data <- function(waves, behavior = NULL, composition = NULL, cc_option = 1) {
  if (!is.list(waves) || length(waves) < 2L) {
    stop("`waves` must be a list of at least two adjacency matrices.", call. = FALSE)
  }
  mats <- lapply(waves, function(w) {
    m <- as.matrix(w)
    if (nrow(m) != ncol(m)) stop("each wave must be a square matrix.", call. = FALSE)
    m
  })
  n <- nrow(mats[[1]])
  if (!all(vapply(mats, nrow, integer(1)) == n)) {
    stop("all waves must share the same actor set (same dimension).", call. = FALSE)
  }
  internal <- .build_saom_data(mats, behavior = behavior, composition = composition,
                               cc_option = cc_option)
  structure(list(internal = internal, n = n, n_waves = length(mats),
    n_periods = length(mats) - 1L, has_behavior = !is.null(behavior),
    has_composition = !is.null(composition)), class = "cusna_data")
}

#' @export
print.cusna_data <- function(x, ...) {
  cat(sprintf("<cusna SAOM data> %d actors, %d waves (%d period%s)\n", x$n, x$n_waves,
              x$n_periods, if (x$n_periods == 1L) "" else "s"))
  if (x$has_behavior) cat("  + co-evolving behavior\n")
  if (x$has_composition) cat("  + composition change\n")
  invisible(x)
}
