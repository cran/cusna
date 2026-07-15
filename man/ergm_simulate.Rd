\name{ergm_simulate}
\alias{ergm_simulate}
\title{Simulate networks from an ERGM (TNT sampler)}
\usage{
ergm_simulate(x, coef, terms, nsim = 100, attr = NULL, directed = FALSE,
              burnin = 16384L, interval = 1024L, seed = 1234L,
              return_nets = FALSE)
}
\arguments{
\item{x}{Integer/logical n*n start adjacency matrix (0/1, zero diagonal).}
\item{coef}{Numeric vector of natural parameters, one per term.}
\item{terms}{A non-empty list of \code{\link{ergm_term}} objects.}
\item{nsim}{Number of independent chains / simulated networks.}
\item{attr}{Optional numeric node covariate of length n; \code{NULL} for none.}
\item{directed}{Logical: is the network directed?}
\item{burnin, interval}{Metropolis burn-in and interval; raise the burn-in to
  ~10x the dyad count for large networks so chains forget the observed start.}
\item{seed}{Integer RNG seed.}
\item{return_nets}{Logical: if \code{TRUE}, also return the simulated networks.}
}
\value{
If \code{return_nets = FALSE} (default) a \code{nsim} x \code{length(terms)} matrix
of sampled statistics (named columns). If \code{return_nets = TRUE}, a list with
\code{stats} (that matrix) and \code{nets} (a list of \code{nsim} n*n integer
adjacency matrices, one per simulated network).
}
\description{
Draws \code{nsim} independent networks with the native TNT Metropolis sampler (the
CPU mirror, validated distributionally against \code{ergm::simulate}) and returns
their sufficient statistics. Full MCMC-MLE for ERGMs remains the domain of the
\code{ergm} package; cusna provides the simulator and MPLE.
}
\examples{
a <- matrix(0L, 6, 6)
a[1, 2] <- a[2, 1] <- a[2, 3] <- a[3, 2] <- 1L
sf <- ergm_simulate(a, coef = -1, terms = list(ergm_term("edges")),
                    nsim = 8, burnin = 2000L, seed = 1L)
colMeans(sf)
}
