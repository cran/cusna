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
A function suitable for assignment to \code{alg$FRAN}.
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
\seealso{
\code{\link{mom_estimate}} for the standalone native estimator.
}
\examples{
\dontrun{
alg <- RSiena::sienaAlgorithmCreate(projname = NULL, cond = FALSE)
alg$FRAN <- cusna_fran(waves = list(w1, w2, w3),
                       effect_names = c("density", "recip", "transTrip"),
                       conditional = FALSE)
ans <- RSiena::siena07(alg, data = dat, effects = eff, useCluster = FALSE)
}
}
