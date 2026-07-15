#' cusna: native GPU-accelerated stochastic actor-oriented models
#'
#' A self-contained native engine (C ABI over CUDA kernels and C++ host logic)
#' for SAOM/RSiena and ERGM, callable from R without a Python runtime. Modelled
#' on the \pkg{torch} package: the CRAN build is CPU-only from source, and the
#' GPU path is compiled from source when a CUDA toolkit is detected at configure
#' time (\code{cusna_has_cuda()} reports which build you have).
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
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  gpu <- tryCatch(cusna_has_cuda_() != 0L, error = function(e) FALSE)
  packageStartupMessage(
    sprintf("cusna native engine (ABI %d) -- %s build",
            cusna_abi_version_(), if (gpu) "GPU" else "CPU-only"))
}
