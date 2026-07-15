\name{ergm_term}
\alias{ergm_term}
\title{Construct an ERGM term}
\usage{
ergm_term(name, parameter = 0)
}
\arguments{
\item{name}{A single ERGM term name: one of \code{edges}, \code{mutual},
  \code{nodecov}, \code{nodematch}, \code{absdiff}, \code{gwesp}, \code{gwdsp}, or
  the temporal (tergm-only) terms \code{memory} (lagged tie) and \code{delrecip}
  (lagged reciprocity).}
\item{parameter}{Numeric term parameter (e.g. the gwesp/gwdsp decay); ignored by
  terms that take none.}
}
\value{
An \code{ergm_term} object.
}
\description{
Construct an ERGM term for use with \code{\link{ergm_stats}} and
\code{\link{ergm_simulate}}.
}
\examples{
ergm_term("edges")
ergm_term("gwesp", 0.25)
}
