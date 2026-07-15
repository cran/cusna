\name{cusna_network_stats}
\alias{cusna_network_stats}
\title{Host network evaluation statistics (RSiena conventions)}
\usage{
cusna_network_stats(a, type_ids, p1, p2, vmat, vmiss, dyn_flags, vb,
                    od0, id0, wtab, symflag = 0L)
}
\arguments{
\item{a}{Integer/logical n*n adjacency matrix (the masked statistic network
  \code{(A & keep) | force}), 0/1, zero diagonal.}
\item{type_ids}{Integer vector (k) of effect type ids.}
\item{p1, p2}{Numeric vectors (k) of effect parameters 1 and 2.}
\item{vmat}{Numeric k*n matrix of per-effect covariate values (rows = effects).}
\item{vmiss}{Integer/logical k*n matrix of covariate missing flags.}
\item{dyn_flags}{Integer vector (k): 1 if the effect covariate is behavior.}
\item{vb}{Numeric vector (n): centered start behavior (for dyn effects).}
\item{od0, id0}{Integer vectors (n): start out/in degrees.}
\item{wtab}{Numeric k*(n+1) matrix of gwesp/gwdsp weight tables (rows = effects).}
\item{symflag}{Integer: 1 for symmetric-network conventions (density halved).}
}
\value{
Numeric vector of length \code{k}.
}
\description{
Thin R wrapper over the native \code{cusna_network_stats}, a bit-exact port of
the reference \code{network_stats_ex}. On public datasets these statistics
reproduce the RSiena targets to machine precision.
}
