\name{tergm_simulate}
\alias{tergm_simulate}
\title{Simulate from a temporal ERGM (TNT sampler with a lagged network)}
\usage{
tergm_simulate(x_start, lag, coef, terms = NULL, attr = NULL,
               mutual = TRUE, memory = TRUE, delrecip = TRUE,
               nsim = 100, burnin = NULL, seed = 1234L,
               return_nets = FALSE)
}
\arguments{
\item{x_start}{Integer/logical n*n start adjacency (0/1, zero diagonal), the
  chain's starting state.}
\item{lag}{Integer/logical n*n lagged network (the previous wave).}
\item{coef}{Numeric parameters, one per active term (matching \code{terms}, or
  the legacy order edges, [mutual], [memory], [delrecip] when \code{terms} is
  \code{NULL}).}
\item{terms}{A list of \code{\link{ergm_term}} objects (may include the temporal
  terms \code{memory} / \code{delrecip}). If \code{NULL} (default) the term list
  is assembled from the deprecated logical flags for backward compatibility.}
\item{attr}{Optional numeric node covariate of length n (for the
  nodecov / nodematch / absdiff cross-sectional terms); \code{NULL} for none.}
\item{mutual, memory, delrecip}{Deprecated: include the corresponding term. Used
  only when \code{terms} is \code{NULL}; the assembled order is edges, mutual,
  memory, delrecip.}
\item{nsim}{Number of independent chains / simulated networks.}
\item{burnin}{MCMC burn-in per chain (default scales with the dyad count).}
\item{seed}{Integer RNG seed.}
\item{return_nets}{Logical: if \code{TRUE}, also return the simulated networks.}
}
\value{
If \code{return_nets = FALSE} (default) a \code{nsim} x k matrix of sampled
statistics (named columns). If \code{return_nets = TRUE}, a list with \code{stats}
(that matrix) and \code{nets} (a list of \code{nsim} n*n integer adjacency
matrices).
}
\description{
Draws \code{nsim} directed networks from the conditional ERGM of the current wave
given the previous one, with the native CPU TNT sampler. The term list may mix the
cross-sectional ERGM terms (edges, mutual, nodecov, nodematch, absdiff) with the
temporal terms memory (lagged tie) and delrecip (lagged reciprocity). Python-free
(CPU engine). The pooled estimator is \code{\link{tergm_mple}}.
}
\examples{
\donttest{
set.seed(1)
lag <- matrix(as.integer(runif(20 * 20) < 0.2), 20, 20); diag(lag) <- 0L
sf <- tergm_simulate(lag, lag, coef = c(-2, 1, 1.5, 0.5), nsim = 20,
                     burnin = 3000L, seed = 1L)
colMeans(sf)
}
}
