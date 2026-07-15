\name{stergm_cmle}
\alias{stergm_cmle}
\alias{print.stergm_cmle}
\title{Separable temporal ERGM by conditional MLE (STERGM CMLE)}
\usage{
stergm_cmle(nets, formation, dissolution = formation, attr = NULL,
            directed = TRUE, nsim = 500L, burnin = NULL, interval = 1024L,
            maxit = 30L, tol = 0.1, damp = 0.5, seed = 1234L, verbose = FALSE)
}
\arguments{
\item{nets}{A list of >= 2 n*n adjacency matrices (0/1, zero diagonal), in time
  order.}
\item{formation}{A non-empty list of [ergm_term()] objects for the formation
  model.}
\item{dissolution}{A list of [ergm_term()] objects for the dissolution
  (persistence) model; defaults to the formation terms.}
\item{attr}{Optional numeric node covariate of length n (for
  nodecov/nodematch/absdiff); \code{NULL} for none.}
\item{directed}{Logical: are the networks directed?}
\item{nsim}{Networks simulated per iteration per transition.}
\item{burnin, interval}{MCMC schedule (default burn-in scales with the dyad
  count).}
\item{maxit, tol, damp}{Outer iterations, convergence tolerance and step
  damping.}
\item{seed}{Integer RNG seed.}
\item{verbose}{Print the per-iteration trajectory.}
}
\value{
A \code{stergm_cmle} object with \code{formation} and \code{dissolution} sub-fits
(each a list with \code{coef}, \code{se}, \code{tstat}, \code{iterations},
\code{converged}, \code{sim_mean}, \code{obs}) and \code{transitions}.
}
\description{
Fits a STERGM (Krivitsky & Handcock 2014) to a network sequence by conditional
maximum likelihood. Each transition separates into two conditionally independent
ERGMs on discord-constrained networks: a formation model on the union
(only dyads absent at t-1 are free to toggle) and a dissolution model on the
intersection (only dyads present at t-1 are free), the latter in the persistence
parameterisation (positive coefficients favour tie persistence, matching
\code{tergm}'s \code{Persist()}). The CMLE is two ERGM MCMC-MLE fits (a reuse of
\code{\link{ergm_mcmle}}) with the native TNT sampler restricted to the free
dyads; statistics are pooled across transitions. Python-free (CPU engine).
}
\examples{
\donttest{
set.seed(1)
base <- matrix(as.integer(runif(20 * 20) < 0.2), 20, 20)
diag(base) <- 0L
nets <- list(base)
for (t in 2:3) {
  nx <- base
  flip <- matrix(runif(20 * 20) < 0.1, 20, 20)
  nx[flip] <- 1L - nx[flip]
  diag(nx) <- 0L
  nets[[t]] <- nx
  base <- nx
}
fit <- stergm_cmle(nets, formation = list(ergm_term("edges"), ergm_term("mutual")),
                   nsim = 200L, maxit = 20L)
fit
}
}
