\name{cusna_gof_distribution}
\alias{cusna_gof_distribution}
\title{sienaGOF auxiliary distributions}
\usage{
cusna_gof_distribution(a, kind, levls)
}
\arguments{
\item{a}{Integer M*n*n array of end networks (first dim = chains).}
\item{kind}{Integer: 0 indegree, 1 outdegree, 2 geodesic.}
\item{levls}{Integer vector of cut levels.}
}
\value{
Numeric matrix with \code{M} rows and \code{length(levls)} columns (kind 0/1) or
\code{length(levls) + 1} columns (kind 2, the extra column is the unreachable
pair count).
}
\description{
Wrapper over the native \code{cusna_gof_distribution} (port of \code{gof.py}):
cumulative indegree/outdegree actor counts, or geodesic-distance ordered-pair
counts for a batch of simulated end networks.
}
\examples{
# Cumulative indegree counts and geodesic-distance pair counts of the three
# observed s50 waves, treated as a batch of three networks.
a <- aperm(simplify2array(s50$waves), c(3, 1, 2))   # 3 x 50 x 50
cusna_gof_distribution(a, kind = 0L, levls = 0:6)
cusna_gof_distribution(a, kind = 2L, levls = 1:5)
}
