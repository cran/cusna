\name{cusna-package}
\alias{cusna-package}
\alias{cusna}
\docType{package}
\title{cusna: native simulation and estimation of network models}
\description{
A self-contained native engine (a C interface over C++ host logic) for
SAOM/RSiena and ERGM, callable from R without a Python runtime. The package
builds the CPU backend of the engine (C++17, multi-threaded with OpenMP; see
\code{\link{cusna_set_threads}}); a CUDA backend implements the same interface
and is built separately from source.
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

The \code{\link{s50}} panel (three waves of friendship among 50 girls, with
alcohol use) is included for the examples.
}
