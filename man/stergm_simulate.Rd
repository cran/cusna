\name{stergm_simulate}
\alias{stergm_simulate}
\title{Simulate a STERGM transition (constrained draw of the next wave)}
\usage{
stergm_simulate(y_lag, formation, dissolution = formation,
                coef_formation, coef_dissolution, attr = NULL,
                directed = TRUE, nsim = 100, burnin = NULL,
                interval = 1024L, seed = 1234L)
}
\arguments{
\item{y_lag}{Integer/logical n*n lagged network \eqn{Y^{t-1}} (0/1, zero
  diagonal).}
\item{formation}{A non-empty list of [ergm_term()] objects for the formation
  model.}
\item{dissolution}{A list of [ergm_term()] objects for the dissolution
  (persistence) model; defaults to the formation terms.}
\item{coef_formation}{Numeric formation parameters, one per formation term.}
\item{coef_dissolution}{Numeric dissolution parameters, one per dissolution
  term.}
\item{attr}{Optional numeric node covariate of length n (for
  nodecov/nodematch/absdiff); \code{NULL} for none.}
\item{directed}{Logical: are the networks directed?}
\item{nsim}{Number of independent transition draws.}
\item{burnin, interval}{MCMC schedule per chain (default burn-in scales with the
  dyad count).}
\item{seed}{Integer RNG seed (dissolution uses \code{seed + 500000L}).}
}
\value{
A list with \code{stats_formation} (a \code{nsim} x k matrix of the formation
sufficient statistics on the union space), \code{stats_dissolution} (the same for
the dissolution/persistence space) and \code{nets} (a list of \code{nsim} n*n
integer adjacency matrices, the composed \eqn{Y^t}).
}
\description{
Draws \code{nsim} networks from a separable temporal ERGM (Krivitsky & Handcock
2014) given the previous wave, the simulation counterpart of
\code{\link{stergm_cmle}}. The transition separates into two conditionally
independent constrained ERGMs, sampled with the native masked TNT sampler on
exactly the discord-constrained space the fit uses: a formation draw
\eqn{F \sim} ERGM(\code{coef_formation}) on \eqn{\{y \supseteq y_{lag}\}} (only
dyads absent at t-1 toggle) and a dissolution draw \eqn{D \sim}
ERGM(\code{coef_dissolution}) on \eqn{\{y \subseteq y_{lag}\}} (persistence
parameterisation, only dyads present at t-1 toggle). The two draws are
independent (separate seed streams, \code{seed} and \code{seed + 500000L}, as in
the fit) and composed into \eqn{y^t = (F \setminus y_{lag}) \cup D}: newly formed
ties where none existed, plus surviving ties where one did. Python-free (CPU
engine).
}
\examples{
\donttest{
set.seed(1)
y_lag <- matrix(as.integer(runif(20 * 20) < 0.2), 20, 20)
diag(y_lag) <- 0L
sim <- stergm_simulate(y_lag,
                       formation = list(ergm_term("edges"), ergm_term("mutual")),
                       coef_formation = c(-2, 0.5),
                       coef_dissolution = c(1.5, 0.2),
                       nsim = 10, burnin = 2000L)
sim$nets[[1]]
}
}
