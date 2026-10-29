#' cusna: native simulation and estimation of network models
#'
#' A self-contained native engine (a C interface over C++ host logic) for
#' SAOM/RSiena and ERGM, callable from R without a Python runtime. The package
#' builds the CPU backend of the engine (C++17, OpenMP); a CUDA backend
#' implements the same interface and is built separately from source.
#'
#' Milestone A5 exposes the validated host statistics -- the RSiena Appendix B
#' convention carriers (\code{cusna_network_stats}, \code{cusna_behavior_stats},
#' \code{cusna_gof_distribution}) -- which reproduce RSiena targets on public
#' datasets to machine precision. The period simulator and the Method-of-Moments
#' estimator are exposed as the native engine matures (see the package roadmap in
#' \code{docs/packaging/track-A-native.md}).
#'
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
  # CRAN checks may use at most two cores; R CMD check --as-cran signals this
  # with _R_CHECK_LIMIT_CORES_, so cap the OpenMP team for examples, tests and
  # vignettes run under it.
  lim <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(lim) && !identical(tolower(lim), "false")) cusna_set_num_threads_(2L)
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  gpu <- tryCatch(cusna_has_cuda_() != 0L, error = function(e) FALSE)
  packageStartupMessage(
    sprintf("cusna native engine (ABI %d) -- %s build",
            cusna_abi_version_(), if (gpu) "GPU" else "CPU-only"))
}
