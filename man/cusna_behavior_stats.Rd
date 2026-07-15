\name{cusna_behavior_stats}
\alias{cusna_behavior_stats}
\title{Host behavior evaluation statistics}
\usage{
cusna_behavior_stats(z, bty, xnet, bkdist, zbar, brange, bsimmean)
}
\arguments{
\item{z}{Integer M*n matrix of end behavior values (rows = chains).}
\item{bty}{Integer vector (kb) of behavior effect type ids
  (0 linear, 1 quad, 2 avAlt, 3 avSim, 4 totSim).}
\item{xnet}{Integer/logical n*n predictor network (start-of-period, less
  missing ties).}
\item{bkdist}{Integer vector (n): 1 if the actor counts (kmask_dist).}
\item{zbar, brange, bsimmean}{Numeric: overall mean, range and simMean.}
}
\value{
Numeric M*kb matrix (rows = chains).
}
\description{
Wrapper over the native \code{cusna_behavior_stats}, a port of
\code{behavior_stats} (RSiena convention 9: totSim subtracts nb*simMean, avSim
divides by nb, avAlt uses the full out-degree).
}
