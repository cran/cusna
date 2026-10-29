\name{cusna_fran}
\alias{cusna_fran}
\title{Native simulation backend for RSiena's siena07()}
\usage{
cusna_fran(waves, effect_names, covariate = NULL, conditional = TRUE)
}
\arguments{
\item{waves}{a list of adjacency matrices (0/1), one per wave, matching the
  \code{RSiena} data object.}
\item{effect_names}{character vector of the included effect names, in the
  same order as the \code{RSiena} effects object rows (after the basic
  rates).}
\item{covariate}{optional actor-covariate values (numeric vector) used by the
  covariate effects (\code{egoX}/\code{altX}/\code{simX}/\code{sameX}/...).}
\item{conditional}{logical; must match the \code{cond} setting of the
  \code{RSiena} algorithm.}
}
\value{
A function suitable for assignment to \code{alg$FRAN}. Each call simulates
one chain per period at the current \code{theta} and returns, following
\pkg{RSiena}'s \code{simstats0c} contract, a list with \code{fra} and
\code{sc} (periods-by-parameters matrices of statistics and scores),
\code{ntim0} (one duration per period), \code{OK} and \code{feasible}. The
number of calls made so far is kept in \code{attr(fran, "counter")$i}.
}
\description{
Returns a \code{FRAN} closure that plugs the native SAOM simulator into the
unmodified \code{RSiena::siena07()} estimator. RSiena keeps its full
Robbins--Monro machinery (and hence its convergence behavior); only the
network simulation runs on the compiled engine. The full bridge (mapping an
arbitrary \code{siena07} effects object onto the native descriptors) is
under development; this closure covers the structural and single-covariate
effects named in \code{effect_names}.
}
\details{
The simulation seeds are drawn from R's random number stream, which
\code{siena07()} initialises from the algorithm's \code{seed}, so a run with
a fixed seed is reproducible. The closure checks that the parameter vector
passed by \code{siena07()} matches \code{effect_names} and
\code{conditional} and stops with an explanatory error otherwise.
Simulated networks are not returned to \pkg{RSiena}
(\code{returnDeps = TRUE} has no effect); use \code{\link{mom_estimate}}
with \code{return_nets = TRUE} for goodness-of-fit work on the native engine.
}
\seealso{
\code{\link{mom_estimate}} for the standalone native estimator.
}
\examples{
\donttest{
if (requireNamespace("RSiena", quietly = TRUE)) {
  # The s50 panel as an RSiena data object ...
  net <- RSiena::sienaDependent(array(unlist(s50$waves), dim = c(50, 50, 3)))
  dat <- RSiena::sienaDataCreate(net)
  eff <- RSiena::includeEffects(RSiena::getEffects(dat), "transTrip",
                                character = TRUE)
  # ... estimated by RSiena's own siena07(), simulating on the native engine.
  alg <- RSiena::sienaAlgorithmCreate(projname = NULL, nsub = 2, n3 = 200,
                                      seed = 123)
  alg$FRAN <- cusna_fran(s50$waves, c("density", "recip", "transTrip"))
  ans <- RSiena::siena07(alg, data = dat, effects = eff, batch = TRUE,
                         verbose = FALSE, silent = TRUE, useCluster = FALSE)
  ans
  attr(alg$FRAN, "counter")$i   # number of simulation calls
}
}
}
