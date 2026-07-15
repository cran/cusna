\name{alaam_mcmle}
\alias{alaam_mcmle}
\alias{print.alaam_mcmle}
\title{Maximum likelihood estimate of an ALAAM (MCMC-MLE)}
\usage{
alaam_mcmle(y, net, cov = NULL, contagion = "out", coef0 = NULL,
            nsim = 500L, sweeps = NULL, maxit = 30L, tol = 0.12,
            damp = 0.5, seed = 1234L, verbose = FALSE)
}
\arguments{
\item{y}{Binary actor attribute (0/1), length n.}
\item{net}{Integer/logical n*n adjacency of the fixed network.}
\item{cov}{Optional actor covariate(s) (length-n vector or n*p matrix).}
\item{contagion}{\code{"out"} (undirected neighbour count, default) or
  \code{"total"} (adds in-neighbours for a directed network).}
\item{coef0}{Optional starting parameters (default: the MPLE).}
\item{nsim}{Chains simulated per iteration (must exceed the parameter count).}
\item{sweeps}{Gibbs sweeps per chain (default scales with n).}
\item{maxit, tol, damp}{Outer iterations, convergence tolerance on the
  standardized deviations, and Hummel-style step damping.}
\item{seed}{Integer RNG seed.}
\item{verbose}{Print the per-iteration trajectory.}
}
\value{
An \code{alaam_mcmle} object with \code{coef}, \code{se} (from the inverse Fisher
information), \code{tstat}, \code{iterations}, \code{converged}, \code{sim_mean}
and \code{obs}.
}
\description{
Fits an ALAAM by MCMC maximum likelihood -- an external moment-matching loop over
the native Gibbs sampler, the attribute-DV twin of \code{\link{ergm_mcmle}}. The
autologistic model is a canonical exponential family, so its Fisher information is
the covariance of the sufficient statistics (density, contagion, covariate
products); each iteration simulates the attribute at the current parameters, forms
their mean and covariance, and takes a damped Newton step until the observed
statistics are matched. Unlike \code{\link{alaam_mple}} (a pseudo-likelihood), the
fitted model reproduces the observed moments. Pure R over the native engine (no
Python).
}
\examples{
\donttest{
set.seed(1)
net <- matrix(as.integer(runif(30 * 30) < 0.1), 30, 30)
net <- ((net + t(net)) > 0) * 1L
diag(net) <- 0L
y <- as.integer(runif(30) < 0.4)
fit <- alaam_mcmle(y, net, nsim = 300L, maxit = 20L)
fit
}
}
