\name{cusna-package}
\alias{cusna-package}
\alias{cusna}
\docType{package}
\title{cusna: native GPU-accelerated stochastic actor-oriented models}
\description{
A self-contained native engine (C ABI over CUDA kernels and C++ host logic) for
SAOM/RSiena and ERGM, callable from R without a Python runtime. Modelled on the
\pkg{torch} package: the CRAN build is CPU-only from source, and the GPU path is
compiled from source when a CUDA toolkit is detected at configure time.
}
\details{
The package covers four model families, all on the same compiled engine:

\itemize{
\item \strong{SAOM} (the RSiena model family): \code{\link{saom_data}} +
  \code{\link{cusna_effect}} + \code{\link{mom_estimate}} run the full
  Method-of-Moments (Robbins--Monro) estimator in R over the native
  simulator, including behavior co-evolution, composition change, and
  multi-network models (\code{\link{mom_estimate_multinet}}). The data
  preparation and moment targets are validated bit-for-bit against the
  reference implementation; the estimates agree within simulation standard
  errors.
\item \strong{ERGM}: \code{\link{ergm_simulate}} (native TNT sampler),
  \code{\link{ergm_stats}}, \code{\link{ergm_mple}}, and the
  \code{\link{ergm_mcmle}} maximum-likelihood estimator.
\item \strong{Temporal ERGM}: \code{\link{tergm_mple}} (pooled MPLE with
  block bootstrap), \code{\link{tergm_simulate}}, and the separable
  \code{\link{stergm_cmle}} (formation/persistence CMLE).
\item \strong{ALAAM} (autologistic actor attribute models):
  \code{\link{alaam_mple}}, \code{\link{alaam_mcmle}},
  \code{\link{alaam_simulate}}.
}

The validated host statistics -- the RSiena Appendix B convention carriers
(\code{\link{cusna_network_stats}}, \code{\link{cusna_behavior_stats}},
\code{\link{cusna_gof_distribution}}) -- reproduce RSiena targets on public
datasets to machine precision and remain available as low-level entry points.
}
