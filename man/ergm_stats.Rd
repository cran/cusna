\name{ergm_stats}
\alias{ergm_stats}
\title{Observed ERGM statistics of a network}
\usage{
ergm_stats(x, terms, attr = NULL, directed = FALSE)
}
\arguments{
\item{x}{Integer/logical n*n adjacency matrix (0/1, zero diagonal).}
\item{terms}{A non-empty list of \code{\link{ergm_term}} objects.}
\item{attr}{Optional numeric node covariate of length n (for
  nodecov/nodematch/absdiff); \code{NULL} for none.}
\item{directed}{Logical: is the network directed (enables \code{mutual})?}
}
\value{
A named numeric vector of length \code{length(terms)}.
}
\description{
Native host statistics for the ERGM term list, a bit-exact port of the reference
\code{network_stats_ergm} (validated to machine zero in \code{native/test}).
}
\examples{
a <- matrix(0L, 4, 4)
a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L
ergm_stats(a, list(ergm_term("edges")))
}
