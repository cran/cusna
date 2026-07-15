#' Native cusna engine: ABI version and CUDA availability
#'
#' `cusna_abi_version()` returns the C ABI version compiled into the package.
#' `cusna_has_cuda()` returns `TRUE` when the package was built with the GPU path
#' (a CUDA toolkit was detected at configure time) and a usable device is
#' present; the CRAN CPU-only build always returns `FALSE`.
#' @return An integer (`cusna_abi_version`) or logical (`cusna_has_cuda`).
#' @export
cusna_abi_version <- function() cusna_abi_version_()

#' @rdname cusna_abi_version
#' @export
cusna_has_cuda <- function() cusna_has_cuda_() != 0L

#' Set the OpenMP thread count for the native CPU backend
#'
#' Explicitly caps the number of threads the CPU backend uses for its
#' OpenMP-parallel batch loops (SAOM simulation, ERGM sampling and MPLE).
#' Unlike `Sys.setenv(OMP_NUM_THREADS = ...)`, this takes effect immediately
#' via the OpenMP runtime API and does not depend on when the environment
#' variable would otherwise be read. Has no effect on a build without OpenMP
#' (e.g. the default macOS toolchain without `libomp`).
#'
#' @param n Requested thread count (coerced to an integer, at least 1).
#' @return Invisibly, the previous maximum thread count (`0L` if OpenMP is
#'   not compiled in).
#' @export
cusna_set_threads <- function(n) invisible(cusna_set_num_threads_(as.integer(n)))

# row-major flatten of an n*n (or r*c) matrix for the C ABI (R is column-major)
.rowmajor <- function(m) as.vector(t(m))

#' Host network evaluation statistics (RSiena conventions)
#'
#' Thin R wrapper over the native `cusna_network_stats` (a bit-exact port of the
#' reference `network_stats_ex`). Given the masked statistic network and the
#' float64 effect descriptors, it returns the `k` effect statistics; on public
#' datasets these reproduce the RSiena targets to machine precision.
#'
#' @param a Integer/logical n*n adjacency matrix (the masked statistic network,
#'   `(A & keep) | force`), 0/1, zero diagonal.
#' @param type_ids Integer vector (k) of effect type ids (see the kernel registry).
#' @param p1,p2 Numeric vectors (k) of effect parameters 1 and 2.
#' @param vmat Numeric k*n matrix of per-effect covariate values (rows = effects).
#' @param vmiss Integer/logical k*n matrix of covariate missing flags.
#' @param dyn_flags Integer vector (k): 1 if the effect covariate is behavior.
#' @param vb Numeric vector (n): centered start behavior (for dyn effects).
#' @param od0,id0 Integer vectors (n): start out/in degrees (parm<=0 degree effects).
#' @param wtab Numeric k*(n+1) matrix of gwesp/gwdsp weight tables (rows = effects).
#' @param symflag Integer: 1 for symmetric-network conventions (density halved).
#' @return Numeric vector of length `k`.
#' @export
cusna_network_stats <- function(a, type_ids, p1, p2, vmat, vmiss, dyn_flags, vb,
                                od0, id0, wtab, symflag = 0L) {
  n <- nrow(a); k <- length(type_ids)
  cusna_network_stats_(
    as.integer(.rowmajor(a)), n, k, as.integer(type_ids),
    as.double(p1), as.double(p2), as.double(.rowmajor(vmat)),
    as.integer(.rowmajor(vmiss)), as.integer(dyn_flags), as.double(vb),
    as.integer(od0), as.integer(id0), as.double(.rowmajor(wtab)),
    as.integer(symflag))
}

#' Host behavior evaluation statistics (linear/quad/avAlt/avSim/totSim)
#'
#' Wrapper over the native `cusna_behavior_stats` (port of `behavior_stats`,
#' RSiena convention 9). Returns the M*kb behavior statistics for a batch of end
#' behavior vectors.
#'
#' @param z Integer M*n matrix of end behavior values (rows = chains).
#' @param bty Integer vector (kb): behavior effect type ids
#'   (0 linear, 1 quad, 2 avAlt, 3 avSim, 4 totSim).
#' @param xnet Integer/logical n*n predictor network (start-of-period, less
#'   missing ties).
#' @param bkdist Integer vector (n): 1 if the actor counts (kmask_dist).
#' @param zbar,brange,bsimmean Numeric: overall mean, range and simMean.
#' @return Numeric M*kb matrix (rows = chains).
#' @export
cusna_behavior_stats <- function(z, bty, xnet, bkdist, zbar, brange, bsimmean) {
  z <- rbind(z); M <- nrow(z); n <- ncol(z); kb <- length(bty)
  out <- cusna_behavior_stats_(
    as.integer(.rowmajor(z)), n, M, kb, as.integer(bty),
    as.integer(.rowmajor(xnet)), as.integer(bkdist),
    as.double(zbar), as.double(brange), as.double(bsimmean))
  matrix(out, nrow = M, byrow = TRUE)
}

#' sienaGOF auxiliary distributions
#'
#' Wrapper over the native `cusna_gof_distribution`. For a batch of `M` end
#' networks, returns the cumulative indegree/outdegree counts (`kind` 0/1) or the
#' geodesic-distance pair counts plus the unreachable count (`kind` 2).
#'
#' @param a Integer M*n*n array of end networks (first dim = chains).
#' @param kind Integer: 0 indegree, 1 outdegree, 2 geodesic.
#' @param levls Integer vector of cut levels.
#' @return Numeric matrix: M rows, `length(levls)` (kind 0/1) or
#'   `length(levls)+1` (kind 2) columns.
#' @export
cusna_gof_distribution <- function(a, kind, levls) {
  d <- dim(a); M <- d[1]; n <- d[2]
  # row-major flatten over (chain, i, j)
  flat <- as.integer(aperm(a, c(3, 2, 1)))
  cols <- if (kind == 2L) length(levls) + 1L else length(levls)
  out <- cusna_gof_distribution_(flat, n, M, as.integer(kind),
                                 as.integer(levls), length(levls))
  matrix(out, nrow = M, byrow = TRUE)
}
