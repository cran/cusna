# =============================================================================
# Method-of-Moments Robbins-Monro estimator -- pure-R port of
# prototype/cusna/saom/mom.py (track A-RAPI, milestone R4). Drives the native
# period simulator (cusna_simulate_period_cpu_) batch-by-batch through the three
# RSiena phases; the derivative (score-covariance), guarded step, diagonalised
# inverse, standard errors and convergence ratios follow mom.py exactly.
#
# Regimes: unconditional [rate x P, rate-eff x R, eval x k]; conditional
# [eval x k] (internal rate 1, stop at the observed distance); coevolution
# [net rate x P, beh rate x P, eval x k, beh x kb]. Bar: estimates track Python
# cpu.py within SE (float32 ABI scores; DESIGN.md 3.1).
# =============================================================================

# Moore-Penrose pseudo-inverse via SVD (no MASS dependency).
.pinv <- function(m, tol = sqrt(.Machine$double.eps)) {
  s <- svd(m)
  keep <- s$d > max(tol * s$d[1], 0)
  s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
}
.safe_solve <- function(m) tryCatch(solve(m), error = function(e) .pinv(m))

# Build the estimator context: targets (R3), parameter layout, per-period
# float64 descriptors, control settings.
.mom_estimator <- function(data, effects, beh_effects = NULL, rate_effects = NULL,
                           conditional = FALSE, control = list()) {
  ctl <- utils::modifyList(list(firstg = 0.4, nsub = 4L, n2 = c(40, 40, 60, 80),
    batch2 = 256L, n1 = 1000L, n3 = 4000L, nD = 1000L, diagonalize = 0.2,
    seed = 1234L, maxdegree = 0L), control)
  n <- data$n; P <- data$n_periods; k <- length(effects)
  coevol <- !is.null(beh_effects) && length(beh_effects) > 0
  kb <- if (coevol) length(beh_effects) else 0L
  R <- if (!is.null(rate_effects)) length(rate_effects) else 0L
  if (R > 0) stop("rate effects are not yet supported in the R estimator (R4b).")
  beh <- data$beh
  symflag <- if (isTRUE(data$symmetric)) 1L else 0L

  mt <- .mom_targets(data, effects, beh_effects, rate_effects, conditional)
  targets <- mt$target_vector
  if (coevol) p <- 2L * P + k + kb else p <- if (conditional) k else P + R + k

  # parameter names, in the theta order (mom.py conventions)
  eff_names <- vapply(effects, function(e) e$name, character(1))
  if (coevol) {
    beh_names <- vapply(beh_effects, function(e) e$name, character(1))
    param_names <- c(sprintf("rate net (p%d)", seq_len(P)),
                     sprintf("rate beh (p%d)", seq_len(P)), eff_names,
                     paste("beh", beh_names))
  } else if (conditional) {
    param_names <- eff_names
  } else {
    rate_names <- if (R > 0) vapply(rate_effects, function(e) e$name, character(1)) else character(0)
    param_names <- c(sprintf("rate (period %d)", seq_len(P)), rate_names, eff_names)
  }
  # rate_consts: zero without composition change (s50 / complete data)
  rate_consts <- if (isTRUE(data$has_composition))
    vapply(seq_len(P), function(pp) as.numeric(.saom_rate_const(data, pp)), numeric(1))
    else rep(0.0, P)
  dists <- vapply(seq_len(P), function(pp) as.numeric(.saom_observed_distance(data, pp)), numeric(1))

  # per-period float64 descriptors
  bty <- if (coevol) vapply(beh_effects, function(e) as.integer(e$type_id), integer(1)) else integer(0)
  per <- vector("list", P)
  for (pp in seq_len(P)) {
    z_start <- if (coevol) .beh_values_less_missing_starts(beh, pp) else NULL
    pe <- .prepare_effects_v2(effects, n, pp, if (coevol) beh else NULL, z_start)
    d <- list(type_ids = pe$type_ids, p1 = pe$p1, p2 = pe$p2, vmat = pe$vmat,
              vmiss = pe$vmiss, dyn_flags = pe$dyn_flags, dyc_idx = pe$dyc_idx,
              v_dyn = pe$v_dyn, wtab = .effect_wtabs(effects, n),
              x0 = data$imputed[[pp]])
    if (coevol) {
      d$z0 <- as.integer(beh$imputed[, pp])
      d$bkdist <- as.integer(.beh_kmask_dist(beh, pp))
      d$xnet <- .saom_network_less_missing_start(data, pp)
    }
    per[[pp]] <- d
  }

  structure(list(data = data, effects = effects, beh_effects = beh_effects,
    conditional = conditional, coevol = coevol, n = n, P = P, k = k, kb = kb,
    R = R, p = p, symflag = symflag, targets = targets, param_names = param_names,
    rate_consts = rate_consts, dists = dists, bty = bty, per = per, beh = beh,
    control = ctl), class = "cusna_mom_est")
}

# row-major flatten for the C ABI
.rm2 <- function(m) as.vector(t(as.matrix(m)))

# One period, M chains, via the native simulator; returns the pieces the loop
# accumulates (network stats/scores, behavior dist/scores/stats, ntim).
.run_period_sim <- function(est, pp, lam, th, M, seed, lamb = 0, thb = numeric(0)) {
  d <- est$per[[pp]]; n <- est$n; k <- est$k; kb <- est$kb
  beh <- est$beh
  tgt <- if (est$conditional) as.integer(est$dists[pp]) else 0L
  r <- cusna_simulate_period_cpu_(
    n, k, M, seed, as.integer(.rm2(d$x0)), lam, as.numeric(th),
    as.integer(d$type_ids), rep(0L, k), as.integer(d$dyn_flags),
    as.integer(d$dyc_idx), .rm2(d$vmat), .rm2(d$wtab), as.numeric(d$v_dyn),
    as.numeric(d$p1), as.numeric(d$p2), as.integer(.rm2(d$vmiss)),
    as.integer(est$control$maxdegree), if (est$conditional) 1L else 0L, tgt,
    est$symflag, kb, lamb, as.numeric(thb), as.integer(est$bty),
    if (kb > 0) d$z0 else integer(0), if (kb > 0) d$bkdist else integer(0),
    if (kb > 0) beh$zmin else 0L, if (kb > 0) beh$zmax else 0L,
    if (kb > 0) beh$overall_mean else 0, if (kb > 0) beh$range else 0,
    if (kb > 0) beh$sim_mean else 0, 0L)
  stats <- matrix(r$stats, M, k, byrow = TRUE)
  scores <- matrix(r$scores, M, k + 1L, byrow = TRUE)
  out <- list(stats = stats, scores = scores, ntim = r$ntim,
              dist_stat = r$dist_stat)
  if (kb > 0) {
    z_end <- matrix(r$z_end, M, n, byrow = TRUE)
    bstats <- cusna_behavior_stats(z_end, est$bty, d$xnet, d$bkdist,
                                   beh$overall_mean, beh$range, beh$sim_mean)
    out$bstats <- bstats
    out$bdist <- r$bdist
    out$bscores <- matrix(r$bscores, M, kb + 1L, byrow = TRUE)
  }
  out
}

# Simulate all periods; returns (stats M x p, pieces list of (sp,scp), ntim M x P).
.mom_simulate <- function(est, theta, M) {
  P <- est$P; p <- est$p; k <- est$k; kb <- est$kb
  if (est$coevol) {
    lam <- theta[1:P]; lamb <- theta[(P + 1):(2 * P)]
    th <- theta[(2 * P + 1):(2 * P + k)]; thb <- theta[(2 * P + k + 1):p]; off <- 2L * P
  } else if (est$conditional) {
    lam <- rep(1, P); th <- theta; thb <- numeric(0); off <- 0L
  } else {
    lam <- theta[1:P]; th <- theta[(P + 1):p]; thb <- numeric(0); off <- P
  }
  stats_acc <- matrix(0, M, p); pieces <- vector("list", P)
  ntim <- matrix(0, M, P)
  for (pp in seq_len(P)) {
    seed <- sample.int(.Machine$integer.max, 1L)
    r <- .run_period_sim(est, pp, as.numeric(lam[pp]), th, M, seed,
                         lamb = if (est$coevol) as.numeric(lamb[pp]) else 0,
                         thb = thb)
    sp <- matrix(0, M, p); scp <- matrix(0, M, p)
    if (!est$conditional) {
      sp[, pp] <- r$dist_stat + est$rate_consts[pp]
      scp[, pp] <- r$scores[, 1]
    }
    if (est$coevol) {
      sp[, P + pp] <- r$bdist
      scp[, P + pp] <- r$bscores[, 1]
      sp[, (off + k + 1):p] <- r$bstats
      scp[, (off + k + 1):p] <- r$bscores[, 2:(kb + 1), drop = FALSE]
    }
    sp[, (off + 1):(off + k)] <- r$stats
    scp[, (off + 1):(off + k)] <- r$scores[, 2:(k + 1), drop = FALSE]
    ntim[, pp] <- r$ntim
    stats_acc <- stats_acc + sp
    pieces[[pp]] <- list(sp = sp, scp = scp)
  }
  list(stats = stats_acc, pieces = pieces, ntim = ntim)
}

# Score-function derivative, summed per period (derivativeFromScoresAndDeviations).
.mom_derivative <- function(pieces) {
  p <- ncol(pieces[[1]]$sp); M <- nrow(pieces[[1]]$sp)
  D <- matrix(0, p, p)
  for (pc in pieces) {
    s_c <- sweep(pc$sp, 2, colMeans(pc$sp))
    sc_c <- sweep(pc$scp, 2, colMeans(pc$scp))
    D <- D + crossprod(s_c, sc_c) / M
  }
  D
}

.mom_dinv <- function(D, diagonalize) {
  mix <- (1 - diagonalize) * D + diagonalize * diag(diag(D))
  .safe_solve(mix)
}

.mom_guard_step <- function(est, theta, step) {
  cap <- pmax(1.0, 0.5 * abs(theta))
  step <- pmin(pmax(step, -cap), cap)
  n_rates <- if (est$conditional) 0L else if (est$coevol) 2L * est$P else est$P
  if (n_rates > 0) for (m in seq_len(n_rates)) {
    if (theta[m] - step[m] < 0.1 * theta[m]) step[m] <- 0.5 * theta[m]
  }
  step
}

.mom_default_theta0 <- function(est) {
  ev0 <- numeric(est$k)
  for (idx in seq_along(est$effects)) if (est$effects[[idx]]$name == "density") ev0[idx] <- -1.0
  if (est$conditional) return(ev0)
  rate_t <- est$targets[1:est$P]
  lam0 <- 1.5 * rate_t / est$n + 1.0
  if (est$coevol) {
    bdists <- est$targets[(est$P + 1):(2 * est$P)]
    lamb0 <- 1.5 * bdists / est$n + 1.0
    return(c(lam0, lamb0, ev0, numeric(est$kb)))
  }
  c(lam0, numeric(est$R), ev0)
}

# Run the Robbins-Monro loop; returns the fit list (theta/se/tconv/D/msf/...).
.mom_run <- function(est, theta0 = NULL, verbose = FALSE, backend = "cpu") {
  ctl <- est$control
  t0 <- Sys.time()
  n_sims <- ctl$n1 + ctl$n3 +
    sum(vapply(seq_len(ctl$nsub), function(s) ctl$nD + ctl$n2[s] * ctl$batch2, numeric(1)))
  set.seed(ctl$seed)
  theta <- if (!is.null(theta0)) as.numeric(theta0) else .mom_default_theta0(est)

  # phase 1
  s <- .mom_simulate(est, theta, ctl$n1)
  D <- .mom_derivative(s$pieces); dinv <- .mom_dinv(D, ctl$diagonalize)
  dev <- colMeans(s$stats) - est$targets
  theta <- theta - .mom_guard_step(est, theta, 0.5 * ctl$firstg * (dinv %*% dev))
  theta <- as.numeric(theta)
  if (verbose) cat("[phase 1] theta =", round(theta, 3), "\n")

  # phase 2
  for (sub in seq_len(ctl$nsub)) {
    gain <- ctl$firstg * (0.5^(sub - 1))
    s <- .mom_simulate(est, theta, ctl$nD)
    D <- .mom_derivative(s$pieces); dinv <- .mom_dinv(D, ctl$diagonalize)
    thav <- numeric(est$p)
    for (it in seq_len(ctl$n2[sub])) {
      s <- .mom_simulate(est, theta, ctl$batch2)
      dev <- colMeans(s$stats) - est$targets
      theta <- theta - .mom_guard_step(est, theta, gain * (dinv %*% dev))
      theta <- as.numeric(theta)
      thav <- thav + theta
    }
    theta <- thav / ctl$n2[sub]
    if (verbose) cat(sprintf("[phase 2.%d] gain=%.3f theta = %s\n", sub, gain,
                             paste(round(theta, 3), collapse = " ")))
  }

  # phase 3
  s <- .mom_simulate(est, theta, ctl$n3)
  D <- .mom_derivative(s$pieces); dinv <- .safe_solve(D)
  msf <- stats::cov(s$stats)
  dev_mean <- colMeans(s$stats) - est$targets
  tconv <- dev_mean / sqrt(diag(msf))
  tconv_max <- as.numeric(sqrt(dev_mean %*% .safe_solve(msf) %*% dev_mean))
  covtheta <- dinv %*% msf %*% t(dinv)
  se <- sqrt(abs(diag(covtheta)))

  rate_est <- if (est$conditional) colMeans(s$ntim) else NULL
  list(theta = theta, se = se, tconv = tconv, tconv_max = tconv_max,
       param_names = est$param_names, targets = est$targets,
       sim_means = colMeans(s$stats), D = D, msf = msf,
       conditional = est$conditional, rate_est = rate_est,
       rate_sd = if (est$conditional) apply(s$ntim, 2, stats::sd) else NULL,
       n_sims = as.integer(n_sims), backend = backend,
       wall_time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

# =============================================================================
# Multi-network co-evolution estimator (R4b) -- port of mom.py MultiMomEstimator.
# Per-variable parameter blocks [rate_net_d x P, eval_net_d ...]; the eval / cross
# statistics are computed on the host from the D end networks (cusna_network_stats
# + cusna_cross_stats), matching cpu.py::simulate_batch_cpu_multi.
# =============================================================================

# Evaluation statistics of one network over a batch of MASKED end networks
# (list of M dense matrices): non-cross effects via cusna_network_stats (per
# chain), cross-network effects via cusna_cross_stats (batched over M).
.batch_net_stats <- function(a_masked, effects, net, cross_pred, p, symflag) {
  k <- length(effects); M <- length(a_masked); n <- net$n
  type_ids <- vapply(effects, function(e) as.integer(e$type_id), integer(1))
  is_cross <- type_ids %in% .CROSS_IDS
  nc_idx <- which(!is_cross); cr_idx <- which(is_cross)
  out <- matrix(0.0, M, k)
  if (length(nc_idx)) {
    pe <- .prepare_effects_v2(effects[nc_idx], n, p, NULL, NULL)
    wtab <- .effect_wtabs(effects[nc_idx], n)
    a_start <- net$imputed[[p]]
    od0 <- as.integer(rowSums(a_start)); id0 <- as.integer(colSums(a_start))
    for (c in seq_len(M)) {
      out[c, nc_idx] <- cusna_network_stats(a_masked[[c]], pe$type_ids, pe$p1,
        pe$p2, pe$vmat, pe$vmiss, pe$dyn_flags, pe$v_dyn, od0, id0, wtab, symflag)
    }
  }
  if (length(cr_idx)) {
    a_flat <- as.integer(unlist(lapply(a_masked, function(a) as.vector(t(a)))))
    for (ci in cr_idx) {
      e <- effects[[ci]]
      wp <- cross_pred[[e$net_ref + 1L]]                # net_ref is 0-based
      out[, ci] <- cusna_cross_stats_(a_flat, as.integer(.rm2(wp)), n, M,
                                      as.integer(e$type_id), as.numeric(e$parm))
    }
  }
  out
}

.multi_mom_estimator <- function(mn, specs_by_net, control = list()) {
  ctl <- utils::modifyList(list(firstg = 0.4, nsub = 4L, n2 = c(40, 40, 60, 80),
    batch2 = 256L, n1 = 1000L, n3 = 4000L, nD = 1000L, diagonalize = 0.2,
    seed = 1234L, maxdegree = 0L), control)
  D <- mn$D; P <- mn$n_periods; nets <- mn$nets; n <- mn$n
  effects <- do.call(c, specs_by_net)
  netof <- unlist(lapply(seq_len(D), function(d) rep(d - 1L, length(specs_by_net[[d]]))))
  k <- length(effects)

  # per-variable layout
  rate_pos <- matrix(0L, D, P); eff_pos <- integer(k); names_v <- character(0)
  pos <- 0L; eidx <- 0L
  for (d in seq_len(D)) {
    for (m in seq_len(P)) { rate_pos[d, m] <- pos; pos <- pos + 1L
      names_v <- c(names_v, sprintf("rate net%d (p%d)", d - 1L, m)) }
    for (e in specs_by_net[[d]]) { eidx <- eidx + 1L; eff_pos[eidx] <- pos; pos <- pos + 1L
      nr <- if (is.null(e$net_ref)) "" else sprintf("<-net%d", e$net_ref)
      names_v <- c(names_v, sprintf("net%d %s%s", d - 1L, e$name, nr)) }
  }
  p <- pos
  symf <- vapply(nets, function(nd) if (isTRUE(nd$symmetric)) 1L else 0L, integer(1))
  rate_targets <- t(vapply(nets, function(nd)
    vapply(seq_len(P), function(m) as.numeric(.saom_rate_target(nd, m)), numeric(1)), numeric(P)))
  rate_consts <- t(vapply(nets, function(nd)
    vapply(seq_len(P), function(m) as.numeric(.saom_rate_const(nd, m)), numeric(1)), numeric(P)))
  if (P == 1L) { rate_targets <- matrix(rate_targets, D, 1); rate_consts <- matrix(rate_consts, D, 1) }

  # per-period cross predictors (period-START, less missing at both ends)
  cross_pred <- lapply(seq_len(P), function(m)
    lapply(nets, function(nd) .saom_network_less_missing_start(nd, m)))

  # targets
  targets <- numeric(p)
  for (d in seq_len(D)) for (m in seq_len(P)) targets[rate_pos[d, m] + 1L] <- rate_targets[d, m]
  eidx <- 0L
  for (d in seq_len(D)) {
    eff_d <- specs_by_net[[d]]
    if (length(eff_d)) {
      for (m in seq_len(P)) {
        tg <- .batch_net_stats(list(.saom_target_network(nets[[d]], m)), eff_d,
                               nets[[d]], cross_pred[[m]], m, symf[d])
        for (q in seq_along(eff_d)) targets[eff_pos[eidx + q] + 1L] <-
          targets[eff_pos[eidx + q] + 1L] + tg[1, q]
      }
    }
    eidx <- eidx + length(eff_d)
  }

  structure(list(mn = mn, nets = nets, D = D, P = P, n = n, k = k, p = p,
    effects = effects, specs_by_net = specs_by_net, netof = netof,
    rate_pos = rate_pos, eff_pos = eff_pos, symf = symf, targets = targets,
    rate_targets = rate_targets, rate_consts = rate_consts, cross_pred = cross_pred,
    param_names = names_v, control = ctl), class = "cusna_multi_mom_est")
}

.multi_default_theta0 <- function(est) {
  th <- numeric(est$p)
  for (d in seq_len(est$D)) for (m in seq_len(est$P))
    th[est$rate_pos[d, m] + 1L] <- 1.5 * est$rate_targets[d, m] / est$n + 1.0
  eidx <- 0L
  for (d in seq_len(est$D)) for (e in est$specs_by_net[[d]]) {
    eidx <- eidx + 1L
    if (e$name == "density") th[est$eff_pos[eidx] + 1L] <- -1.0
  }
  th
}

.multi_guard_step <- function(est, theta, step) {
  cap <- pmax(1.0, 0.5 * abs(theta))
  step <- pmin(pmax(step, -cap), cap)
  for (d in seq_len(est$D)) for (m in seq_len(est$P)) {
    q <- est$rate_pos[d, m] + 1L
    if (theta[q] - step[q] < 0.1 * theta[q]) step[q] <- 0.5 * theta[q]
  }
  step
}

.multi_mom_simulate <- function(est, theta, M) {
  D <- est$D; P <- est$P; k <- est$k; p <- est$p; n <- est$n; nets <- est$nets
  th_flat <- theta[est$eff_pos + 1L]
  pe <- .prepare_effects_v2(est$effects, n, 1L, NULL, NULL)   # descriptors (period-invariant covariate here)
  wtab <- .effect_wtabs(est$effects, n)
  etypes <- .effect_etypes(est$effects)
  netref <- .effect_netrefs(est$effects)
  stats_acc <- matrix(0, M, p); pieces <- vector("list", P); ntim <- matrix(0, M, P)
  for (m in seq_len(P)) {
    lam_d <- vapply(seq_len(D), function(d) theta[est$rate_pos[d, m] + 1L], numeric(1))
    seed <- sample.int(.Machine$integer.max, 1L)
    x0 <- as.integer(unlist(lapply(nets, function(nd) .rm2(nd$imputed[[m]]))))
    perm <- as.integer(unlist(lapply(nets, function(nd) .rm2(.saom_permitted(nd, m)))))
    pe_m <- .prepare_effects_v2(est$effects, n, m, NULL, NULL)
    r <- cusna_simulate_multi_cpu_(n, k, M, D, seed, x0, as.numeric(lam_d),
      as.numeric(th_flat), as.integer(pe_m$type_ids), as.integer(etypes),
      as.integer(est$netof), as.integer(netref), .rm2(wtab), as.numeric(pe_m$p1),
      as.numeric(pe_m$p2), .rm2(pe_m$vmat), perm, as.integer(est$symf))
    rscores <- matrix(r$rscores, M, D, byrow = TRUE)
    scores <- matrix(r$scores, M, k, byrow = TRUE)
    a_ends <- r$a_ends                                    # length M*D*n*n row-major
    sp <- matrix(0, M, p); scp <- matrix(0, M, p)
    eidx <- 0L
    for (d in seq_len(D)) {
      x0d <- nets[[d]]$imputed[[m]]
      kdist <- .saom_kmask_dist(nets[[d]], m)
      kf <- .saom_stat_keep_force(nets[[d]], m)
      # slice per-chain end network of net d
      a_list <- vector("list", M); dist_stat <- numeric(M)
      for (c in seq_len(M)) {
        off <- ((c - 1L) * D + (d - 1L)) * n * n
        ad <- matrix(a_ends[(off + 1L):(off + n * n)], n, n, byrow = TRUE)
        dist_stat[c] <- sum((ad != x0d) & kdist)
        a_list[[c]] <- ((ad != 0L) & kf$keep) | kf$force      # masked statistic net
      }
      q <- est$rate_pos[d, m] + 1L
      sp[, q] <- dist_stat + est$rate_consts[d, m]
      scp[, q] <- rscores[, d]
      eff_d <- est$specs_by_net[[d]]
      if (length(eff_d)) {
        sd <- .batch_net_stats(a_list, eff_d, nets[[d]], est$cross_pred[[m]], m, est$symf[d])
        for (qq in seq_along(eff_d)) {
          sp[, est$eff_pos[eidx + qq] + 1L] <- sd[, qq]
          scp[, est$eff_pos[eidx + qq] + 1L] <- scores[, eidx + qq]
        }
      }
      eidx <- eidx + length(eff_d)
    }
    ntim[, m] <- r$ntim
    stats_acc <- stats_acc + sp
    pieces[[m]] <- list(sp = sp, scp = scp)
  }
  list(stats = stats_acc, pieces = pieces, ntim = ntim)
}

.multi_mom_run <- function(est, theta0 = NULL, verbose = FALSE, backend = "cpu") {
  ctl <- est$control; t0 <- Sys.time()
  n_sims <- ctl$n1 + ctl$n3 +
    sum(vapply(seq_len(ctl$nsub), function(s) ctl$nD + ctl$n2[s] * ctl$batch2, numeric(1)))
  set.seed(ctl$seed)
  theta <- if (!is.null(theta0)) as.numeric(theta0) else .multi_default_theta0(est)

  s <- .multi_mom_simulate(est, theta, ctl$n1)
  D <- .mom_derivative(s$pieces); dinv <- .mom_dinv(D, ctl$diagonalize)
  dev <- colMeans(s$stats) - est$targets
  theta <- as.numeric(theta - .multi_guard_step(est, theta, 0.5 * ctl$firstg * (dinv %*% dev)))
  if (verbose) cat("[phase 1] theta =", round(theta, 3), "\n")

  for (sub in seq_len(ctl$nsub)) {
    gain <- ctl$firstg * (0.5^(sub - 1))
    s <- .multi_mom_simulate(est, theta, ctl$nD)
    D <- .mom_derivative(s$pieces); dinv <- .mom_dinv(D, ctl$diagonalize)
    thav <- numeric(est$p)
    for (it in seq_len(ctl$n2[sub])) {
      s <- .multi_mom_simulate(est, theta, ctl$batch2)
      dev <- colMeans(s$stats) - est$targets
      theta <- as.numeric(theta - .multi_guard_step(est, theta, gain * (dinv %*% dev)))
      thav <- thav + theta
    }
    theta <- thav / ctl$n2[sub]
    if (verbose) cat(sprintf("[phase 2.%d] theta = %s\n", sub, paste(round(theta, 3), collapse = " ")))
  }

  s <- .multi_mom_simulate(est, theta, ctl$n3)
  D <- .mom_derivative(s$pieces); dinv <- .safe_solve(D)
  msf <- stats::cov(s$stats); dev_mean <- colMeans(s$stats) - est$targets
  tconv <- dev_mean / sqrt(diag(msf))
  tconv_max <- as.numeric(sqrt(dev_mean %*% .safe_solve(msf) %*% dev_mean))
  se <- sqrt(abs(diag(dinv %*% msf %*% t(dinv))))
  list(theta = theta, se = se, tconv = tconv, tconv_max = tconv_max,
       param_names = est$param_names, targets = est$targets, sim_means = colMeans(s$stats),
       D = D, msf = msf, conditional = FALSE, rate_est = NULL, rate_sd = NULL,
       n_sims = as.integer(n_sims), backend = backend,
       wall_time = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
