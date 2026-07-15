\name{alaam_simulate}
\alias{alaam_simulate}
\title{Simulate the actor attribute from an ALAAM (Gibbs sampler)}
\usage{
alaam_simulate(net, coef, cov = NULL, y0 = NULL, contagion = "out",
               nsim = 100L, sweeps = 200L, seed = 1234L)
}
\arguments{
\item{net}{Integer/logical n*n adjacency of the fixed network.}
\item{coef}{Parameters in the order \code{density} (intercept), \code{contagion},
  then one per covariate column of \code{cov}.}
\item{cov}{Optional actor covariate(s) (length-n vector or n*p matrix) matching
  the covariate coefficients in \code{coef}.}
\item{y0}{Optional start attribute (length n); default all zero.}
\item{contagion}{\code{"out"} (default) uses out-neighbours; \code{"total"} adds
  in-neighbours (the joint-consistent contagion for a directed network).}
\item{nsim}{Number of independent chains.}
\item{sweeps}{Gibbs sweeps per chain (burn-in over the attribute vector).}
\item{seed}{Integer RNG seed.}
}
\value{
A list with \code{y} (a \code{nsim} x n matrix of simulated attributes) and
\code{stats} (a \code{nsim} x 2 matrix of \code{density} and \code{contagion}).
}
\description{
Draws \code{nsim} binary attribute vectors from an ALAAM at a fixed network with
the native Gibbs sampler: each full conditional is a logistic regression of
\code{y_i} on the contagion (the number of neighbours with \code{y = 1}, the
masked popcount used by \code{\link{alaam_mple}}, so estimator and simulator share
one statistic) plus the actor covariates. Python-free (host engine, CPU-only).
}
\examples{
\donttest{
set.seed(1)
net <- matrix(as.integer(runif(30 * 30) < 0.1), 30, 30)
net <- ((net + t(net)) > 0) * 1L        # undirected
diag(net) <- 0L
sim <- alaam_simulate(net, coef = c(density = -1, contagion = 0.3), nsim = 20L)
colMeans(sim$stats)
}
}
