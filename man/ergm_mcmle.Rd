\name{ergm_mcmle}
\alias{ergm_mcmle}
\alias{print.ergm_mcmle}
\title{Maximum likelihood estimate of an ERGM (MCMC-MLE)}
\usage{
ergm_mcmle(x, terms, attr = NULL, directed = FALSE, coef0 = NULL,
           nsim = 500L, burnin = NULL, interval = 1024L, maxit = 25L,
           tol = 0.1, damp = 0.5, seed = 1234L, verbose = FALSE)
}
\arguments{
\item{x}{Integer/logical n*n observed adjacency (0/1, zero diagonal).}
\item{terms}{A non-empty list of \code{\link{ergm_term}} objects.}
\item{attr}{Optional numeric node covariate of length n; \code{NULL} for none.}
\item{directed}{Logical: is the network directed?}
\item{coef0}{Optional starting parameters (default: edges from density, rest 0).}
\item{nsim}{Networks simulated per iteration (must exceed the number of terms).}
\item{burnin, interval}{MCMC schedule per iteration; chains start from the observed
  network, so a modest burn-in decorrelates near the MLE. Default burn-in scales
  with the dyad count.}
\item{maxit}{Maximum outer iterations.}
\item{tol}{Convergence: stop when every standardized deviation \code{|t|} is below it.}
\item{damp}{Step-length damping (Hummel-style partial stepping when far off).}
\item{seed}{Integer RNG seed.}
\item{verbose}{Print the per-iteration trajectory.}
}
\value{
An \code{ergm_mcmle} object with \code{coef}, \code{se} (from the inverse Fisher
information), \code{tstat}, \code{iterations}, \code{converged}, \code{sim_mean}
and \code{obs}.
}
\description{
Fits an ERGM by MCMC maximum likelihood: an external moment-matching loop over the
native TNT sampler. For the canonical exponential family the Fisher information is
the covariance of the sufficient statistics, so each iteration simulates from the
current parameters, forms the mean and covariance of the sufficient statistics, and
takes a damped Newton step until the observed statistics are matched. The loop is
pure R over the native engine (no Python). This closes the sampler + MPLE ->
maximum likelihood gap; full MCMC-MLE for arbitrary specifications remains the
domain of the \code{ergm} package.
}
\examples{
\donttest{
set.seed(1)
a <- matrix(as.integer(runif(18 * 18) < 0.3), 18, 18)
diag(a) <- 0L
fit <- ergm_mcmle(a, list(ergm_term("edges"), ergm_term("mutual")),
                  directed = TRUE, nsim = 200L, seed = 1L)
fit
}
}
