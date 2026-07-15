# ERGM building blocks (native): observed statistics, a batched TNT simulator and
# a maximum pseudo-likelihood estimator, over the native C ABI (cpp11) -- no
# Python. The interface mirrors the reticulate-backed track-B package; the bodies
# call the vendored CPU host code (cusna_ergm_*_), so they work in the CPU-only
# CRAN build. Full MCMC-MLE for ERGMs remains the domain of the 'ergm' package;
# cusna provides the simulator and MPLE (framing M-3).

# term-name -> type id registry (matches ergm_kernels.cu / tnt_gpu _TERM_IDS).
# ids 7 (memory) / 8 (delrecip) are the TEMPORAL terms: they are only meaningful
# with a lagged network (tergm_simulate) and are NOT cross-sectional statistics,
# so ergm_stats / ergm_simulate reject them (see .ergm_tergm_only_ids).
.ergm_term_ids <- c(edges = 0L, mutual = 1L, nodecov = 2L, nodematch = 3L,
                    absdiff = 4L, gwesp = 5L, gwdsp = 6L,
                    memory = 7L, delrecip = 8L)

# temporal-only term ids -- valid in tergm_simulate (with a lag), not elsewhere
.ergm_tergm_only_ids <- c(7L, 8L)

#' Construct an ERGM term
#'
#' @param name A single ERGM term name: one of `edges`, `mutual`, `nodecov`,
#'   `nodematch`, `absdiff`, `gwesp`, `gwdsp`, or the temporal (tergm-only)
#'   terms `memory` (lagged tie) and `delrecip` (lagged reciprocity).
#' @param parameter Numeric term parameter (e.g. the gwesp/gwdsp decay); ignored
#'   by terms that take none.
#' @return An `ergm_term` object.
#' @export
ergm_term <- function(name, parameter = 0) {
  if (!is.character(name) || length(name) != 1L) {
    stop("`name` must be a single ERGM term name (a string).", call. = FALSE)
  }
  if (is.na(.ergm_term_ids[name])) {
    stop("unknown ERGM term '", name, "'; one of: ",
         paste(names(.ergm_term_ids), collapse = ", "), call. = FALSE)
  }
  structure(list(name = name, parameter = parameter), class = "ergm_term")
}

#' @export
print.ergm_term <- function(x, ...) {
  cat("<ergm term>", x$name)
  if (!identical(x$parameter, 0) && !identical(x$parameter, 0L)) {
    cat(" (parm =", x$parameter, ")")
  }
  cat("\n")
  invisible(x)
}

# split a list of ergm_term objects into parallel id / parm / name vectors
.ergm_terms_split <- function(terms) {
  if (!is.list(terms) || !length(terms)) {
    stop("`terms` must be a non-empty list of ergm_term() objects.",
         call. = FALSE)
  }
  ids <- integer(length(terms))
  parm <- numeric(length(terms))
  nm <- character(length(terms))
  for (i in seq_along(terms)) {
    t <- terms[[i]]
    if (!inherits(t, "ergm_term")) {
      stop("each term must be an ergm_term() object.", call. = FALSE)
    }
    ids[i] <- .ergm_term_ids[[t$name]]
    parm[i] <- as.numeric(t$parameter)
    nm[i] <- t$name
  }
  list(ids = ids, parm = parm, names = nm)
}

#' Observed ERGM statistics of a network
#'
#' Native host statistics for the ERGM term list (a bit-exact port of the
#' reference `network_stats_ergm`).
#'
#' @param x Integer/logical n*n adjacency matrix (0/1, zero diagonal).
#' @param terms A non-empty list of [ergm_term()] objects.
#' @param attr Optional numeric node covariate of length n (for
#'   nodecov/nodematch/absdiff); `NULL` for none.
#' @param directed Logical: is the network directed (enables `mutual`)?
#' @return A named numeric vector of length `length(terms)`.
#' @export
# reject temporal-only terms where no lagged network is available
.ergm_reject_temporal <- function(ts, where) {
  bad <- ts$names[ts$ids %in% .ergm_tergm_only_ids]
  if (length(bad)) {
    stop("temporal term(s) ", paste(unique(bad), collapse = ", "),
         " need a lagged network; use tergm_simulate(), not ", where, ".",
         call. = FALSE)
  }
}

ergm_stats <- function(x, terms, attr = NULL, directed = FALSE) {
  n <- nrow(x)
  ts <- .ergm_terms_split(terms)
  .ergm_reject_temporal(ts, "ergm_stats()")
  av <- if (is.null(attr)) numeric(0) else as.double(attr)
  out <- cusna_ergm_stats_(as.integer(.rowmajor(x)), n, length(ts$ids),
                           as.integer(directed), ts$ids, ts$parm, av)
  stats <- as.numeric(out)
  names(stats) <- ts$names
  stats
}

#' Simulate networks from an ERGM (TNT sampler)
#'
#' Draws `nsim` independent networks with the native TNT Metropolis sampler (the
#' CPU mirror, validated distributionally against `ergm::simulate`) and returns
#' their sufficient statistics.
#'
#' @param x Integer/logical n*n start adjacency matrix (0/1, zero diagonal).
#' @param coef Numeric vector of natural parameters, one per term.
#' @param terms A non-empty list of [ergm_term()] objects.
#' @param nsim Number of independent chains / simulated networks.
#' @param attr Optional numeric node covariate of length n; `NULL` for none.
#' @param directed Logical: is the network directed?
#' @param burnin,interval Metropolis burn-in and (unused for one sample per chain)
#'   interval; the default burn-in should be raised to ~10x the dyad count for
#'   large networks so chains forget the observed start.
#' @param seed Integer RNG seed.
#' @param return_nets Logical: if `TRUE`, also return the simulated networks.
#' @return If `return_nets = FALSE` (default) a `nsim` x `length(terms)` matrix of
#'   sampled statistics (named columns). If `return_nets = TRUE`, a list with
#'   `stats` (that matrix) and `nets` (a list of `nsim` n*n integer adjacency
#'   matrices, one per simulated network).
#' @export
ergm_simulate <- function(x, coef, terms, nsim = 100, attr = NULL,
                          directed = FALSE, burnin = 16384L, interval = 1024L,
                          seed = 1234L, return_nets = FALSE) {
  n <- nrow(x)
  ts <- .ergm_terms_split(terms)
  .ergm_reject_temporal(ts, "ergm_simulate()")
  if (length(coef) != length(terms)) {
    stop("`coef` must have one value per term.", call. = FALSE)
  }
  av <- if (is.null(attr)) numeric(0) else as.double(attr)
  res <- cusna_ergm_simulate_cpu_(
    as.integer(.rowmajor(x)), n, length(ts$ids), as.integer(directed),
    as.double(coef), ts$ids, ts$parm, av, as.integer(nsim), as.double(seed),
    as.integer(burnin), as.integer(interval), 1L, as.integer(isTRUE(return_nets)))
  flat <- if (isTRUE(return_nets)) res$stats else res
  m <- matrix(flat, nrow = nsim, ncol = length(ts$ids), byrow = TRUE)
  colnames(m) <- ts$names
  if (!isTRUE(return_nets)) return(m)
  list(stats = m, nets = .ergm_unpack_nets(res$nets, nsim, n))
}

# unpack a C*n*n row-major integer vector into a list of nsim n*n R matrices
.ergm_unpack_nets <- function(flat, nsim, n) {
  lapply(seq_len(nsim), function(c) {
    off <- (c - 1L) * n * n
    matrix(flat[(off + 1L):(off + n * n)], n, n, byrow = TRUE)
  })
}

#' Maximum pseudo-likelihood estimate (directed edges/mutual/ttriple demo)
#'
#' A logistic regression of tie presence on the dyadic change statistics the
#' native engine computes (Strauss & Ikeda 1990). This is the MPLE demo for the
#' three directed terms edges / mutual / transitive-triple; full MCMC-MLE stays
#' with the 'ergm' package (framing M-3).
#'
#' @param x Integer/logical n*n directed adjacency matrix (0/1, zero diagonal).
#' @return A named numeric vector of pseudo-MLE coefficients
#'   (`edges`, `mutual`, `transitive`).
#' @export
ergm_mple <- function(x) {
  n <- nrow(x)
  cs <- cusna_ergm_mple_(as.integer(.rowmajor(x)), n)
  X <- matrix(cs$X, nrow = cs$nrow, ncol = cs$ncol, byrow = TRUE)
  colnames(X) <- c("edges", "mutual", "transitive")
  y <- as.integer(cs$y)
  fit <- stats::glm.fit(X, y, family = stats::binomial())
  coefs <- fit$coefficients
  names(coefs) <- colnames(X)
  coefs
}

# numerically-safe inverse: fall back to a tiny ridge if singular
.ergm_safe_solve <- function(M) {
  k <- nrow(M)
  tryCatch(solve(M),
           error = function(e)
             solve(M + diag(1e-8 * mean(abs(diag(M))) + 1e-12, k)))
}

# starting parameters: edges from the observed density (logit), rest 0
.ergm_mcmle_seed <- function(gobs, ts, n, directed) {
  ndyads <- if (directed) n * (n - 1) else n * (n - 1) / 2
  th <- numeric(length(ts$ids))
  ei <- which(ts$ids == 0L)                    # edges has type id 0
  if (length(ei)) {
    dens <- min(max(as.numeric(gobs[ei]) / ndyads, 1e-4), 1 - 1e-4)
    th[ei] <- stats::qlogis(dens)
  }
  th
}

#' Maximum likelihood estimate of an ERGM (MCMC-MLE)
#'
#' Fits an ERGM by MCMC maximum likelihood -- an external moment-matching loop
#' over the native TNT sampler. For the canonical exponential family the Fisher
#' information is the covariance of the sufficient statistics, so each iteration
#' simulates from the current parameters, forms the mean and covariance of the
#' sufficient statistics, and takes a damped Newton step
#' \eqn{\theta \leftarrow \theta + \gamma\,\Sigma^{-1}(g_{obs}-\bar g)} until the
#' observed statistics are matched. This closes the sampler+MPLE -> MLE gap; the
#' loop is pure R over the native engine (no Python, no C ABI change).
#'
#' @param x Integer/logical n*n observed adjacency (0/1, zero diagonal).
#' @param terms A non-empty list of [ergm_term()] objects.
#' @param attr Optional numeric node covariate of length n; `NULL` for none.
#' @param directed Logical: is the network directed?
#' @param coef0 Optional starting parameters (default: edges from density, rest 0).
#' @param nsim Networks simulated per iteration (must exceed the number of terms).
#' @param burnin,interval MCMC schedule per iteration; the chains start from the
#'   observed network, so a modest burn-in decorrelates near the MLE. Default
#'   burn-in scales with the dyad count.
#' @param maxit Maximum outer iterations.
#' @param tol Convergence: stop when every standardized deviation |t| is below it.
#' @param damp Step-length damping (Hummel-style partial stepping when far off).
#' @param seed Integer RNG seed.
#' @param verbose Print the per-iteration trajectory.
#' @return An `ergm_mcmle` object: `coef`, `se` (from the inverse Fisher
#'   information), `tstat`, `iterations`, `converged`, `sim_mean`, `obs`.
#' @export
ergm_mcmle <- function(x, terms, attr = NULL, directed = FALSE, coef0 = NULL,
                       nsim = 500L, burnin = NULL, interval = 1024L,
                       maxit = 25L, tol = 0.1, damp = 0.5, seed = 1234L,
                       verbose = FALSE) {
  ts <- .ergm_terms_split(terms)
  k <- length(ts$ids)
  n <- nrow(x)
  if (nsim <= k) stop("`nsim` must exceed the number of terms.", call. = FALSE)
  ndyads <- if (directed) n * (n - 1) else n * (n - 1) / 2
  if (is.null(burnin)) burnin <- as.integer(max(4096, 3 * ndyads))
  gobs <- ergm_stats(x, terms, attr = attr, directed = directed)
  theta <- if (is.null(coef0)) .ergm_mcmle_seed(gobs, ts, n, directed)
           else as.numeric(coef0)

  Sig <- diag(k)
  gbar <- gobs
  tstat <- rep(Inf, k)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    G <- ergm_simulate(x, theta, terms, nsim = nsim, attr = attr,
                       directed = directed, burnin = burnin, interval = interval,
                       seed = seed + 1000L * it)
    gbar <- colMeans(G)
    Sig <- stats::cov(G)
    dev <- as.numeric(gobs) - gbar
    sdg <- sqrt(pmax(diag(Sig), 1e-12))
    tstat <- dev / sdg
    if (max(abs(tstat)) < tol) { converged <- TRUE; break }
    Siginv <- .ergm_safe_solve(Sig)
    step <- as.numeric(Siginv %*% dev)
    # Hummel-style partial step: a big Mahalanobis gap (observed far outside the
    # simulated cloud) shrinks the step; near the solution take the full Newton.
    md <- sqrt(max(as.numeric(dev %*% Siginv %*% dev), 0))
    gamma <- max(min(1, damp * sqrt(k) / max(md, 1e-9)), 0.05)
    cap <- pmax(1.0, 0.5 * abs(theta))          # guard against runaway steps
    theta <- theta + pmin(pmax(gamma * step, -cap), cap)
    if (verbose) {
      cat(sprintf("[mcmle %2d] gamma=%.3f max|t|=%.3f  theta = %s\n",
                  it, gamma, max(abs(tstat)),
                  paste(sprintf("%.4f", theta), collapse = "  ")))
    }
  }
  se <- sqrt(pmax(diag(.ergm_safe_solve(Sig)), 0))
  structure(list(
    coef = stats::setNames(theta, ts$names),
    se = stats::setNames(se, ts$names),
    tstat = stats::setNames(as.numeric(tstat), ts$names),
    iterations = it, converged = converged,
    sim_mean = stats::setNames(as.numeric(gbar), ts$names),
    obs = stats::setNames(as.numeric(gobs), ts$names)),
    class = "ergm_mcmle")
}

#' @export
print.ergm_mcmle <- function(x, ...) {
  cat("ERGM MCMC-MLE (cusna native)\n")
  tab <- data.frame(estimate = round(x$coef, 4), s.e. = round(x$se, 4),
                    `t-conv` = round(x$tstat, 3), check.names = FALSE)
  print(tab)
  cat(sprintf("iterations: %d   converged: %s   max|t-conv|: %.3f\n",
              x$iterations, x$converged, max(abs(x$tstat))))
  invisible(x)
}

# off-diagonal dyad order (i outer, j inner, i != j) matching cusna_ergm_mple_
.tergm_dyad_order <- function(n) {
  g <- expand.grid(j = seq_len(n), i = seq_len(n))[, c("i", "j")]
  g[g$i != g$j, , drop = FALSE]
}

#' Temporal ERGM by bootstrap pseudo-likelihood (btergm style)
#'
#' Fits a temporal ERGM to a sequence of directed networks by pooled maximum
#' pseudo-likelihood -- the btergm estimator (Desmarais & Cranmer 2012). For every
#' transition \eqn{Y^{t-1} \to Y^t} the native engine builds the 0->1 change
#' statistics of every dyad of \eqn{Y^t} (edges / mutual / transitive triple, via
#' \code{\link{ergm_mple}}'s design matrix) and the temporal terms are read off the
#' lagged network (\code{memory} = \eqn{Y^{t-1}_{ij}}, \code{delrecip} =
#' \eqn{Y^{t-1}_{ji}}); the rows are pooled across transitions and fit by weighted
#' logistic regression. Confidence intervals come from a bootstrap over the pooled
#' cases. The design matrix is host C++ (no Python; works in the CPU-only build).
#'
#' @param nets A list of >= 2 directed n*n adjacency matrices (0/1, zero diagonal),
#'   in time order.
#' @param mutual,ttriple Include the cross-sectional mutual / transitive-triple
#'   terms (from the native change-statistic matrix).
#' @param memory,delrecip Include the temporal terms (lagged tie / lagged
#'   reciprocity).
#' @param R Bootstrap replications for the confidence intervals.
#' @param level Confidence level for the bootstrap intervals.
#' @param seed Integer RNG seed for the bootstrap.
#' @return A `tergm_mple` object: `coef`, `ci` (bootstrap interval), `boot` (the
#'   replicate coefficients), `nobs` (pooled dyads), `transitions`.
#' @export
tergm_mple <- function(nets, mutual = TRUE, ttriple = TRUE, memory = TRUE,
                       delrecip = TRUE, R = 200L, level = 0.95, seed = 1234L) {
  if (!is.list(nets) || length(nets) < 2L) {
    stop("`nets` must be a list of >= 2 adjacency matrices.", call. = FALSE)
  }
  n <- nrow(nets[[1]])
  ord <- .tergm_dyad_order(n)
  ii <- ord$i
  jj <- ord$j

  # one design/response block per transition (kept separate for the block
  # bootstrap, which resamples whole time steps -- the btergm scheme)
  Xblocks <- vector("list", length(nets) - 1L)
  yblocks <- vector("list", length(nets) - 1L)
  for (t in 2:length(nets)) {
    Yt <- nets[[t]]
    Ylag <- nets[[t - 1]]
    cs <- cusna_ergm_mple_(as.integer(.rowmajor(Yt)), n)
    Xc <- matrix(cs$X, nrow = cs$nrow, ncol = 3L, byrow = TRUE)  # edges,mutual,ttriple
    cols <- list(edges = Xc[, 1])
    if (mutual)   cols$mutual   <- Xc[, 2]
    if (ttriple)  cols$ttriple  <- Xc[, 3]
    if (memory)   cols$memory   <- Ylag[cbind(ii, jj)]           # lagged tie i->j
    if (delrecip) cols$delrecip <- Ylag[cbind(jj, ii)]           # lagged tie j->i
    Xblocks[[t - 1]] <- do.call(cbind, cols)
    yblocks[[t - 1]] <- as.integer(cs$y)
  }
  Xp <- do.call(rbind, Xblocks)
  yp <- unlist(yblocks)

  fit <- stats::glm.fit(Xp, yp, family = stats::binomial())
  coefs <- fit$coefficients
  names(coefs) <- colnames(Xp)

  # block bootstrap: resample whole transitions with replacement (btergm's
  # bootstrap resamples time steps, not dyads, to respect temporal dependence)
  set.seed(seed)
  nb <- length(Xblocks)
  B <- matrix(NA_real_, R, length(coefs))
  for (r in seq_len(R)) {
    sel <- sample.int(nb, nb, replace = TRUE)
    Xb <- do.call(rbind, Xblocks[sel])
    yb <- unlist(yblocks[sel])
    fr <- tryCatch(
      stats::glm.fit(Xb, yb, family = stats::binomial()),
      error = function(e) NULL)
    if (!is.null(fr)) B[r, ] <- fr$coefficients
  }
  a <- (1 - level) / 2
  ci <- t(apply(B, 2, stats::quantile, c(a, 1 - a), na.rm = TRUE))
  rownames(ci) <- names(coefs)
  colnames(ci) <- c("lower", "upper")

  structure(list(coef = coefs, ci = ci, boot = B, nobs = nrow(Xp),
                 transitions = length(nets) - 1L, level = level),
            class = "tergm_mple")
}

#' @export
print.tergm_mple <- function(x, ...) {
  cat(sprintf("Temporal ERGM (bootstrap MPLE, cusna native): %d transitions, %d pooled dyads\n",
              x$transitions, x$nobs))
  tab <- data.frame(estimate = round(x$coef, 4),
                    lower = round(x$ci[, 1], 4),
                    upper = round(x$ci[, 2], 4), check.names = FALSE)
  print(tab)
  cat(sprintf("(%.0f%% bootstrap interval, %d replications)\n",
              100 * x$level, nrow(x$boot)))
  invisible(x)
}

# start-network absolute statistics (s0) for a tergm term list: cross-sectional
# terms via ergm_stats(directed = TRUE), the temporal terms memory / delrecip by
# their closed forms over the lagged network. Returns an s0 vector aligned to ids.
.tergm_s0 <- function(x_start, lag, ids, parm, nms, attr) {
  s0 <- numeric(length(ids))
  cross <- !(ids %in% .ergm_tergm_only_ids)
  if (any(cross)) {
    ct <- lapply(which(cross),
                 function(i) ergm_term(nms[i], parameter = parm[i]))
    s0[cross] <- as.numeric(ergm_stats(x_start, ct, attr = attr, directed = TRUE))
  }
  for (i in which(!cross)) {
    s0[i] <- switch(as.character(ids[i]),
                    "7" = sum(x_start * lag),          # memory: lagged tie i->j
                    "8" = sum(x_start * t(lag)))       # delrecip: lagged tie j->i
  }
  s0
}

#' Simulate from a temporal ERGM (TNT sampler with a lagged network)
#'
#' Draws `nsim` directed networks from the conditional ERGM
#' \eqn{Y^t \mid Y^{t-1}} with the native CPU TNT sampler. The term list may mix
#' the cross-sectional ERGM terms (edges, mutual, nodecov, nodematch, absdiff)
#' with the temporal terms memory (lagged tie) and delrecip (lagged reciprocity).
#' Python-free (CPU engine). Returns the sampled sufficient statistics and,
#' optionally, the simulated networks.
#'
#' @param x_start Integer/logical n*n start adjacency (0/1, zero diagonal), the
#'   chain's starting state.
#' @param lag Integer/logical n*n lagged network \eqn{Y^{t-1}}.
#' @param coef Numeric parameters, one per active term (matching `terms`, or the
#'   legacy order edges, [mutual], [memory], [delrecip] when `terms` is `NULL`).
#' @param terms A list of [ergm_term()] objects (may include the temporal terms
#'   `memory` / `delrecip`). If `NULL` (default) the term list is assembled from
#'   the deprecated logical flags `mutual` / `memory` / `delrecip` for backward
#'   compatibility.
#' @param attr Optional numeric node covariate of length n (for the
#'   nodecov / nodematch / absdiff cross-sectional terms); `NULL` for none.
#' @param mutual,memory,delrecip Deprecated: include the corresponding term. Used
#'   only when `terms` is `NULL`; the assembled order is edges, mutual, memory,
#'   delrecip (as before).
#' @param nsim Number of independent chains / simulated networks.
#' @param burnin MCMC burn-in per chain (default scales with the dyad count).
#' @param seed Integer RNG seed.
#' @param return_nets Logical: if `TRUE`, also return the simulated networks.
#' @return If `return_nets = FALSE` (default) a `nsim` x k matrix of sampled
#'   statistics (named columns). If `return_nets = TRUE`, a list with `stats`
#'   (that matrix) and `nets` (a list of `nsim` n*n integer adjacency matrices).
#' @export
tergm_simulate <- function(x_start, lag, coef, terms = NULL, attr = NULL,
                           mutual = TRUE, memory = TRUE, delrecip = TRUE,
                           nsim = 100, burnin = NULL, seed = 1234L,
                           return_nets = FALSE) {
  n <- nrow(x_start)
  if (is.null(terms)) {
    # legacy path: assemble ids from the logical flags (order edges, mutual,
    # memory, delrecip), parm all zero, no attr -- byte-identical to <=0.1.0.
    ids <- 0L; nms <- "edges"; parm <- 0
    if (mutual)   { ids <- c(ids, 1L); nms <- c(nms, "mutual");   parm <- c(parm, 0) }
    if (memory)   { ids <- c(ids, 7L); nms <- c(nms, "memory");   parm <- c(parm, 0) }
    if (delrecip) { ids <- c(ids, 8L); nms <- c(nms, "delrecip"); parm <- c(parm, 0) }
    av <- numeric(0)
  } else {
    ts <- .ergm_terms_split(terms)
    ids <- ts$ids; nms <- ts$names; parm <- ts$parm
    av <- if (is.null(attr)) numeric(0) else as.double(attr)
  }
  k <- length(ids)
  if (length(coef) != k) {
    stop("`coef` must have one value per active term.", call. = FALSE)
  }
  s0 <- .tergm_s0(x_start, lag, ids, parm, nms,
                  attr = if (length(av)) av else NULL)
  ndyads <- n * (n - 1)
  if (is.null(burnin)) burnin <- as.integer(max(4096, 3 * ndyads))
  res <- cusna_ergm_tsim_cpu_(
    as.integer(.rowmajor(x_start)), as.integer(.rowmajor(lag)), n, k, 1L,
    as.double(coef), as.integer(ids), as.double(parm), as.double(av),
    as.double(s0), as.integer(nsim), as.double(seed), as.integer(burnin), 1024L,
    1L, as.integer(isTRUE(return_nets)))
  flat <- if (isTRUE(return_nets)) res$stats else res
  m <- matrix(flat, nrow = nsim, ncol = k, byrow = TRUE)
  colnames(m) <- nms
  if (!isTRUE(return_nets)) return(m)
  list(stats = m, nets = .ergm_unpack_nets(res$nets, nsim, n))
}

# ============================================================================
#  ALAAM (autologistic actor attribute model) -- attribute-DV extension (X1.5).
#  The dependent variable is a BINARY actor attribute at a FIXED network; the
#  contagion change statistic is the masked popcount of the network row against
#  the attribute vector (cusna_alaam_stats). Python-free (host engine, CPU-only).
# ============================================================================

# ALAAM contagion change statistics of an attribute at a fixed network. Returns
# a list(con_out, con_in): out/in-neighbours with the attribute (coincide when
# the network is undirected/symmetric = the neighbour count with y=1).
.alaam_contagion <- function(net, y) {
  n <- nrow(net)
  cusna_alaam_stats_(as.integer(.rowmajor(net)), as.integer(y != 0), n)
}

#' Autologistic actor attribute model by maximum pseudo-likelihood (ALAAM)
#'
#' Fits an ALAAM (Robins, Pattison & Elliott 2001) -- a model for a binary actor
#' attribute at a fixed, exogenous network -- by maximum pseudo-likelihood, i.e. a
#' logistic regression of the attribute on the native contagion change statistic
#' (the number of network neighbours that also carry the attribute, a masked
#' popcount) plus any actor covariates. The dependent variable is the attribute,
#' not a tie; the network is the same bit-kernel primitive applied to
#' (network-row AND attribute-vector). Best suited to an undirected network, where
#' the contagion is the neighbour count with `y = 1`; for a directed network
#' `contagion = "both"` includes out- and in-contagion as separate terms.
#'
#' @param y Binary actor attribute (0/1), length n -- the dependent variable.
#' @param net Integer/logical n*n adjacency of the fixed (exogenous) network.
#' @param cov Optional actor covariate(s): a length-n vector or an n*p matrix,
#'   entered as change statistics (the covariate value at each actor).
#' @param contagion Which contagion change statistic to include: `"out"`
#'   (out-neighbours with `y = 1`, the default and the neighbour count when
#'   undirected), `"in"`, `"both"`, or `FALSE` for none.
#' @param directed Logical: is `net` directed? (documentation only; the contagion
#'   columns are selected by `contagion`).
#' @return An `alaam_mple` object: `coef`, `se`, `z`, `pval` (from the logistic
#'   regression), `glm` (the fitted model), `n`.
#' @export
alaam_mple <- function(y, net, cov = NULL, contagion = "out", directed = FALSE) {
  n <- nrow(net)
  yv <- as.integer(y != 0)
  if (length(yv) != n) stop("`y` must have length nrow(net).", call. = FALSE)
  cs <- .alaam_contagion(net, yv)
  cols <- list()
  if (isTRUE(contagion) || identical(contagion, "out")) {
    cols$contagion <- cs$con_out
  } else if (identical(contagion, "in")) {
    cols$contagion <- cs$con_in
  } else if (identical(contagion, "both")) {
    cols$contagion_out <- cs$con_out
    cols$contagion_in  <- cs$con_in
  } else if (!identical(contagion, FALSE)) {
    stop("`contagion` must be one of \"out\", \"in\", \"both\", or FALSE.",
         call. = FALSE)
  }
  if (!is.null(cov)) {
    covm <- as.matrix(cov)
    cn <- colnames(covm)
    if (is.null(cn)) cn <- if (ncol(covm) == 1L) "cov" else paste0("cov", seq_len(ncol(covm)))
    for (j in seq_len(ncol(covm))) cols[[cn[j]]] <- as.numeric(covm[, j])
  }
  if (!length(cols)) stop("ALAAM needs at least one predictor (contagion or cov).",
                          call. = FALSE)
  df <- data.frame(c(list(y = yv), cols), check.names = FALSE)
  fit <- stats::glm(y ~ ., data = df, family = stats::binomial())
  sm <- stats::coef(summary(fit))
  co <- sm[, 1]; se <- sm[, 2]
  nm <- rownames(sm); nm[nm == "(Intercept)"] <- "density"
  names(co) <- nm; names(se) <- nm
  structure(list(coef = co, se = se, z = sm[, 3], pval = sm[, 4],
                 glm = fit, n = n),
            class = "alaam_mple")
}

#' @export
print.alaam_mple <- function(x, ...) {
  cat(sprintf("ALAAM (autologistic actor attribute, MPLE; cusna native): n = %d actors\n",
              x$n))
  tab <- data.frame(estimate = round(x$coef, 4), s.e. = round(x$se, 4),
                    z = round(x$z, 3), check.names = FALSE)
  print(tab)
  invisible(x)
}

#' Simulate the actor attribute from an ALAAM (Gibbs sampler)
#'
#' Draws `nsim` binary attribute vectors from an ALAAM at a fixed network with the
#' native Gibbs sampler: each full conditional is
#' \eqn{\mathrm{logit}\,P(y_i=1\mid y_{-i}) = \beta_0 + \beta_c c_i + \beta_x' x_i},
#' with the contagion \eqn{c_i} the number of neighbours with `y = 1` (the masked
#' popcount used by \code{\link{alaam_mple}}, so estimator and simulator share one
#' statistic). Python-free (host engine, CPU-only).
#'
#' @param net Integer/logical n*n adjacency of the fixed network.
#' @param coef Parameters in the order `density` (intercept), `contagion`, then
#'   one per covariate column of `cov`.
#' @param cov Optional actor covariate(s) (length-n vector or n*p matrix) matching
#'   the covariate coefficients in `coef`.
#' @param y0 Optional start attribute (length n); default all zero.
#' @param contagion `"out"` (default) uses out-neighbours; `"total"` adds
#'   in-neighbours (the joint-consistent contagion for a directed network).
#' @param nsim Number of independent chains.
#' @param sweeps Gibbs sweeps per chain (burn-in over the attribute vector).
#' @param seed Integer RNG seed.
#' @return A list with `y` (a `nsim` x n matrix of simulated attributes) and
#'   `stats` (a `nsim` x 2 matrix of `density` and `contagion`).
#' @export
alaam_simulate <- function(net, coef, cov = NULL, y0 = NULL,
                           contagion = "out", nsim = 100L, sweeps = 200L,
                           seed = 1234L) {
  n <- nrow(net)
  if (length(coef) < 2L) stop("`coef` needs at least density and contagion.",
                              call. = FALSE)
  beta0 <- as.double(coef[1]); beta_c <- as.double(coef[2])
  beta_x <- if (length(coef) > 2L) as.double(coef[-(1:2)]) else numeric(0)
  p <- length(beta_x)
  Xflat <- if (p > 0L) {
    covm <- matrix(as.numeric(as.matrix(cov)), nrow = n)
    if (ncol(covm) != p) stop("`cov` columns must match the covariate coefficients.",
                              call. = FALSE)
    as.double(.rowmajor(covm))
  } else numeric(0)
  if (is.null(y0)) y0 <- integer(n)
  add_in <- if (identical(contagion, "total")) 1L else 0L
  res <- cusna_alaam_gibbs_cpu_(
    as.integer(.rowmajor(net)), n, add_in, as.integer(y0 != 0), Xflat,
    as.integer(p), beta0, beta_c, beta_x, as.integer(nsim), as.double(seed),
    as.integer(sweeps))
  ymat <- matrix(res$y_ends, nrow = nsim, ncol = n, byrow = TRUE)
  smat <- matrix(res$stats, nrow = nsim, ncol = 2L, byrow = TRUE)
  colnames(smat) <- c("density", "contagion")
  list(y = ymat, stats = smat)
}

# covariate matrix + parameter names for the single-contagion ALAAM
.alaam_design <- function(net, cov, contagion) {
  n <- nrow(net)
  covm <- if (is.null(cov)) NULL else matrix(as.numeric(as.matrix(cov)), nrow = n)
  p <- if (is.null(covm)) 0L else ncol(covm)
  cn <- if (p > 0L) colnames(covm) else character(0)
  if (p > 0L && is.null(cn)) cn <- if (p == 1L) "cov" else paste0("cov", seq_len(p))
  halfc <- if (identical(contagion, "total")) 1.0 else 0.5
  list(covm = covm, p = p, names = c("density", "contagion", cn), halfc = halfc)
}

#' Maximum likelihood estimate of an ALAAM (MCMC-MLE)
#'
#' Fits an ALAAM by MCMC maximum likelihood -- an external moment-matching loop
#' over the native Gibbs sampler, the attribute-DV twin of \code{\link{ergm_mcmle}}.
#' The autologistic model is a canonical exponential family, so its Fisher
#' information is the covariance of the sufficient statistics (density, contagion,
#' covariate products); each iteration simulates the attribute at the current
#' parameters, forms their mean and covariance, and takes a damped Newton step
#' until the observed statistics are matched. Unlike \code{\link{alaam_mple}} (a
#' pseudo-likelihood), the fitted model reproduces the observed moments. Pure R
#' over the native engine (no Python).
#'
#' @param y Binary actor attribute (0/1), length n.
#' @param net Integer/logical n*n adjacency of the fixed network.
#' @param cov Optional actor covariate(s) (length-n vector or n*p matrix).
#' @param contagion `"out"` (undirected neighbour count, default) or `"total"`.
#' @param coef0 Optional starting parameters (default: the MPLE).
#' @param nsim Chains simulated per iteration (must exceed the parameter count).
#' @param sweeps Gibbs sweeps per chain (default scales with n).
#' @param maxit,tol,damp Outer iterations, convergence tolerance on |t|, and
#'   Hummel-style step damping.
#' @param seed Integer RNG seed.
#' @param verbose Print the per-iteration trajectory.
#' @return An `alaam_mcmle` object: `coef`, `se` (from the inverse Fisher
#'   information), `tstat`, `iterations`, `converged`, `sim_mean`, `obs`.
#' @export
alaam_mcmle <- function(y, net, cov = NULL, contagion = "out", coef0 = NULL,
                        nsim = 500L, sweeps = NULL, maxit = 30L, tol = 0.12,
                        damp = 0.5, seed = 1234L, verbose = FALSE) {
  n <- nrow(net)
  yv <- as.integer(y != 0)
  des <- .alaam_design(net, cov, contagion)
  covm <- des$covm; p <- des$p; nm <- des$names; halfc <- des$halfc
  k <- 2L + p
  if (nsim <= k) stop("`nsim` must exceed the number of parameters.", call. = FALSE)
  if (is.null(sweeps)) sweeps <- as.integer(max(300, 6 * n))

  # observed sufficient statistics: density, contagion, covariate products
  Ay <- as.numeric(net %*% yv)
  gobs <- c(sum(yv), halfc * sum(yv * Ay))
  if (p > 0L) gobs <- c(gobs, colSums(covm * yv))
  names(gobs) <- nm

  theta <- if (!is.null(coef0)) as.numeric(coef0)
           else unname(alaam_mple(yv, net, cov = cov, contagion = "out")$coef)

  Sig <- diag(k); gbar <- gobs; tstat <- rep(Inf, k); converged <- FALSE
  for (it in seq_len(maxit)) {
    sim <- alaam_simulate(net, theta, cov = covm, contagion = contagion,
                          nsim = nsim, sweeps = sweeps, seed = seed + 1000L * it)
    G <- sim$stats
    if (p > 0L) G <- cbind(G, sim$y %*% covm)
    colnames(G) <- nm
    gbar <- colMeans(G); Sig <- stats::cov(G)
    dev <- gobs - gbar
    sdg <- sqrt(pmax(diag(Sig), 1e-12)); tstat <- dev / sdg
    if (max(abs(tstat)) < tol) { converged <- TRUE; break }
    Siginv <- .ergm_safe_solve(Sig)
    step <- as.numeric(Siginv %*% dev)
    md <- sqrt(max(as.numeric(dev %*% Siginv %*% dev), 0))
    gamma <- max(min(1, damp * sqrt(k) / max(md, 1e-9)), 0.05)
    cap <- pmax(0.5, 0.5 * abs(theta))
    theta <- theta + pmin(pmax(gamma * step, -cap), cap)
    if (verbose) {
      cat(sprintf("[alaam mcmle %2d] gamma=%.3f max|t|=%.3f  theta = %s\n",
                  it, gamma, max(abs(tstat)),
                  paste(sprintf("%.4f", theta), collapse = "  ")))
    }
  }
  se <- sqrt(pmax(diag(.ergm_safe_solve(Sig)), 0))
  structure(list(
    coef = stats::setNames(theta, nm), se = stats::setNames(se, nm),
    tstat = stats::setNames(as.numeric(tstat), nm),
    iterations = it, converged = converged,
    sim_mean = stats::setNames(as.numeric(gbar), nm),
    obs = stats::setNames(as.numeric(gobs), nm)),
    class = "alaam_mcmle")
}

#' @export
print.alaam_mcmle <- function(x, ...) {
  cat("ALAAM MCMC-MLE (cusna native)\n")
  tab <- data.frame(estimate = round(x$coef, 4), s.e. = round(x$se, 4),
                    `t-conv` = round(x$tstat, 3), check.names = FALSE)
  print(tab)
  cat(sprintf("iterations: %d   converged: %s   max|t-conv|: %.3f\n",
              x$iterations, x$converged, max(abs(x$tstat))))
  invisible(x)
}

# ============================================================================
#  STERGM (separable temporal ERGM) by conditional MLE (X2).
#  Each transition Y^{t-1} -> Y^t splits into two conditionally independent
#  ERGMs on discord-constrained networks (Krivitsky & Handcock 2014):
#    formation    Y+ = Y^{t-1} | Y^t  (union), free dyads = absent at t-1;
#    dissolution  Y- = Y^{t-1} & Y^t  (intersection, persistence
#                 parameterisation), free dyads = present at t-1.
#  CMLE = ERGM MCMC-MLE twice (reuse of X0), each with the native TNT sampler
#  restricted to the free dyads (the permitted mask). Python-free (CPU engine).
# ============================================================================

# pooled MCMC-MLE for one separable phase (formation or dissolution). blocks is a
# list of per-transition list(x0, permitted); the moment equation sums observed
# and simulated sufficient statistics over transitions, the Fisher information is
# the summed covariance (transitions are conditionally independent).
.stergm_fit_phase <- function(blocks, gobs, seed_theta, ts, n, directed,
                              nsim, burnin, interval, maxit, tol, damp, seed,
                              verbose, label) {
  k <- length(ts$ids)
  theta <- seed_theta
  Sig <- diag(k); gbar <- gobs; tstat <- rep(Inf, k); converged <- FALSE
  for (it in seq_len(maxit)) {
    means <- matrix(0, length(blocks), k)
    Sig <- matrix(0, k, k)
    for (bi in seq_along(blocks)) {
      b <- blocks[[bi]]
      flat <- cusna_ergm_sim_masked_cpu_(
        as.integer(.rowmajor(b$x0)), as.integer(.rowmajor(b$permitted)), n, k,
        as.integer(directed), as.double(theta), ts$ids, ts$parm, as.double(b$attr),
        as.integer(nsim), as.double(seed + 1000L * it + bi),
        as.integer(burnin), as.integer(interval), 1L)
      G <- matrix(flat, nrow = nsim, ncol = k, byrow = TRUE)
      means[bi, ] <- colMeans(G)
      Sig <- Sig + stats::cov(G)                 # summed Fisher info (independent)
    }
    gbar <- colSums(means)
    dev <- gobs - gbar
    sdg <- sqrt(pmax(diag(Sig), 1e-12)); tstat <- dev / sdg
    if (max(abs(tstat)) < tol) { converged <- TRUE; break }
    Siginv <- .ergm_safe_solve(Sig)
    step <- as.numeric(Siginv %*% dev)
    md <- sqrt(max(as.numeric(dev %*% Siginv %*% dev), 0))
    gamma <- max(min(1, damp * sqrt(k) / max(md, 1e-9)), 0.05)
    cap <- pmax(1.0, 0.5 * abs(theta))
    theta <- theta + pmin(pmax(gamma * step, -cap), cap)
    if (verbose) {
      cat(sprintf("[stergm %-11s %2d] gamma=%.3f max|t|=%.3f  theta = %s\n",
                  label, it, gamma, max(abs(tstat)),
                  paste(sprintf("%.4f", theta), collapse = "  ")))
    }
  }
  se <- sqrt(pmax(diag(.ergm_safe_solve(Sig)), 0))
  list(coef = stats::setNames(theta, ts$names), se = stats::setNames(se, ts$names),
       tstat = stats::setNames(as.numeric(tstat), ts$names),
       iterations = it, converged = converged,
       sim_mean = stats::setNames(as.numeric(gbar), ts$names),
       obs = stats::setNames(as.numeric(gobs), ts$names))
}

# build the constrained networks, observed statistics and a seed for one phase
.stergm_phase_setup <- function(nets, terms, phase, attr, directed) {
  ts <- .ergm_terms_split(terms)
  k <- length(ts$ids); n <- nrow(nets[[1]])
  blocks <- vector("list", length(nets) - 1L)
  gobs <- numeric(k); on_free <- 0; n_free <- 0
  for (t in 2:length(nets)) {
    Ylag <- (nets[[t - 1]] != 0) * 1L
    Yt   <- (nets[[t]] != 0) * 1L
    if (phase == "formation") {
      x0 <- ((Ylag | Yt) != 0) * 1L               # union
      permitted <- (Ylag == 0) * 1L               # free = absent at t-1
    } else {
      x0 <- (Ylag & Yt) * 1L                       # intersection (persistence)
      permitted <- (Ylag != 0) * 1L               # free = present at t-1
    }
    diag(x0) <- 0L; diag(permitted) <- 0L
    s0 <- ergm_stats(x0, terms, attr = attr, directed = directed)
    gobs <- gobs + as.numeric(s0)
    fr <- permitted == 1L; diag(fr) <- FALSE
    on_free <- on_free + sum(Yt[fr]); n_free <- n_free + sum(fr)
    blocks[[t - 1]] <- list(x0 = x0, permitted = permitted,
                            attr = if (is.null(attr)) numeric(0) else as.double(attr))
  }
  # seed: edges from the pooled free-dyad on-rate (formation/persistence), rest 0
  theta <- numeric(k)
  ei <- which(ts$ids == 0L)
  if (length(ei) && n_free > 0) {
    rate <- min(max(on_free / n_free, 1e-3), 1 - 1e-3)
    theta[ei] <- stats::qlogis(rate)
  }
  list(ts = ts, blocks = blocks, gobs = gobs, seed_theta = theta, n = n)
}

#' Separable temporal ERGM by conditional MLE (STERGM CMLE)
#'
#' Fits a STERGM (Krivitsky & Handcock 2014) to a network sequence by conditional
#' maximum likelihood. Each transition \eqn{Y^{t-1} \to Y^t} separates into two
#' conditionally independent ERGMs on discord-constrained networks: a *formation*
#' model on the union \eqn{Y^{t-1}\cup Y^t} (only dyads absent at \eqn{t-1} are
#' free to toggle) and a *dissolution* model on the intersection
#' \eqn{Y^{t-1}\cap Y^t} (only dyads present at \eqn{t-1} are free), the latter in
#' the persistence parameterisation (positive coefficients favour tie
#' persistence, matching \code{tergm}'s \code{Persist()}). The CMLE is two ERGM
#' MCMC-MLE fits (a reuse of \code{\link{ergm_mcmle}}) with the native TNT sampler
#' restricted to the free dyads; statistics are pooled across transitions.
#' Python-free (CPU engine).
#'
#' @param nets A list of >= 2 n*n adjacency matrices (0/1, zero diagonal), in
#'   time order.
#' @param formation A non-empty list of [ergm_term()] objects for the formation
#'   model.
#' @param dissolution A list of [ergm_term()] objects for the dissolution
#'   (persistence) model; defaults to the formation terms.
#' @param attr Optional numeric node covariate of length n (for
#'   nodecov/nodematch/absdiff); `NULL` for none.
#' @param directed Logical: are the networks directed?
#' @param nsim Networks simulated per iteration per transition.
#' @param burnin,interval MCMC schedule (default burn-in scales with the dyad
#'   count).
#' @param maxit,tol,damp Outer iterations, convergence tolerance and step damping.
#' @param seed Integer RNG seed.
#' @param verbose Print the per-iteration trajectory.
#' @return A `stergm_cmle` object with `formation` and `dissolution` sub-fits
#'   (each with `coef`, `se`, `tstat`, ...) and `transitions`.
#' @export
stergm_cmle <- function(nets, formation, dissolution = formation, attr = NULL,
                        directed = TRUE, nsim = 500L, burnin = NULL,
                        interval = 1024L, maxit = 30L, tol = 0.1, damp = 0.5,
                        seed = 1234L, verbose = FALSE) {
  if (!is.list(nets) || length(nets) < 2L) {
    stop("`nets` must be a list of >= 2 adjacency matrices.", call. = FALSE)
  }
  n <- nrow(nets[[1]])
  ndyads <- if (directed) n * (n - 1) else n * (n - 1) / 2
  if (is.null(burnin)) burnin <- as.integer(max(4096, 2 * ndyads))

  fs <- .stergm_phase_setup(nets, formation, "formation", attr, directed)
  ds <- .stergm_phase_setup(nets, dissolution, "dissolution", attr, directed)

  form <- .stergm_fit_phase(fs$blocks, fs$gobs, fs$seed_theta, fs$ts, n, directed,
                            nsim, burnin, interval, maxit, tol, damp, seed,
                            verbose, "formation")
  diss <- .stergm_fit_phase(ds$blocks, ds$gobs, ds$seed_theta, ds$ts, n, directed,
                            nsim, burnin, interval, maxit, tol, damp,
                            seed + 500000L, verbose, "dissolution")

  structure(list(formation = form, dissolution = diss,
                 transitions = length(nets) - 1L),
            class = "stergm_cmle")
}

#' @export
print.stergm_cmle <- function(x, ...) {
  cat(sprintf("Separable temporal ERGM (CMLE, cusna native): %d transitions\n",
              x$transitions))
  fmt <- function(f, title) {
    cat("\n", title, " (converged: ", f$converged, ", max|t-conv|: ",
        sprintf("%.3f", max(abs(f$tstat))), ")\n", sep = "")
    print(data.frame(estimate = round(f$coef, 4), s.e. = round(f$se, 4),
                     check.names = FALSE))
  }
  fmt(x$formation, "Formation")
  fmt(x$dissolution, "Dissolution (persistence)")
  invisible(x)
}
