\name{mom_estimate}
\alias{mom_estimate}
\alias{mom_control}
\title{Method-of-Moments estimation of a SAOM}
\usage{
mom_estimate(data, effects, beh_effects = NULL, rate_effects = NULL,
             conditional = FALSE, backend = c("cpu", "gpu"), maxdegree = 0,
             control = mom_control(), verbose = FALSE, return_nets = FALSE)

mom_control(firstg = 0.4, nsub = 4, n2 = c(40, 40, 60, 80),
            batch2 = 256, n1 = 1000, n3 = 4000, nD = 1000,
            diagonalize = 0.2, seed = 1234)
}
\arguments{
\item{data}{a \code{cusna_data} object from \code{\link{saom_data}}.}
\item{effects}{a non-empty list of \code{\link{cusna_effect}} objects
  (network evaluation/creation/endowment effects) and, optionally,
  \code{\link{cusna_interaction}} objects placed after their components.}
\item{beh_effects}{optional list of \code{\link{cusna_beh_effect}} objects for
  network--behavior co-evolution (requires \code{behavior} in the data;
  forces unconditional estimation).}
\item{rate_effects}{optional list of \code{\link{cusna_rate_effect}} objects
  (\code{outRate}/\code{inRate}/\code{RateX}); unconditional only.}
\item{conditional}{logical; \code{FALSE} (default) estimates the basic rate
  parameters jointly with the evaluation effects. \code{TRUE} conditions each
  period on the observed distance (RSiena's default for a single dependent
  network) and reports the rates separately.}
\item{backend}{\code{"cpu"} (default; the from-source CRAN build) or
  \code{"gpu"} (needs a CUDA build).}
\item{maxdegree}{optional maximum out-degree constraint (0 = none).}
\item{control}{a list of tuning parameters from \code{mom_control}.}
\item{verbose}{logical; if \code{TRUE}, print per-phase progress.}
\item{return_nets}{logical; if \code{TRUE}, keep the end networks of the
  phase-3 simulations (RSiena's \code{returnDeps}) in the returned fit as
  \code{sims}: a list with one integer array of dimension
  \code{c(n3, n, n)} per period (chains first, the layout
  \code{\link{cusna_gof_distribution}} takes). Co-evolution fits also return
  the simulated behavior as \code{sims_behavior} (one \code{n3 x n} matrix
  per period). The arrays hold \code{n3 * n^2} integers per period, so
  choose \code{n3} accordingly.}
\item{firstg}{initial Robbins--Monro gain.}
\item{nsub}{number of phase-2 subphases.}
\item{n2}{integer vector (length \code{>= nsub}) of phase-2 iteration counts
  per subphase.}
\item{batch2}{simulation batch size within phase-2 iterations.}
\item{n1}{number of phase-1 simulations (derivative estimation).}
\item{n3}{number of phase-3 simulations (standard errors and convergence).}
\item{nD}{simulations used to re-estimate the derivative each subphase.}
\item{diagonalize}{shrinkage of the derivative matrix towards its diagonal,
  in [0, 1]; higher is more robust on ill-conditioned models.}
\item{seed}{random seed for the simulation stream.}
}
\value{
\code{mom_estimate} returns a \code{\link{cusna_fit}} object.
\code{mom_control} returns a list of tuning parameters
(class \code{cusna_control}).
}
\description{
Fits a stochastic actor-oriented model by the RSiena-style method of moments
(Robbins--Monro), driving the compiled native simulator -- no Python. The
moment targets and per-period masks reproduce the reference implementation
bit-for-bit; the estimates agree with the reference estimator within
simulation standard errors.
}
\details{
For nearly collinear specifications the Robbins--Monro estimator may not
converge from a cold start; increase \code{diagonalize} or the phase counts,
or run \pkg{RSiena}'s own estimator on the native simulator with
\code{\link{cusna_fran}}.

The default schedule (65,320 simulated panels) is deliberately larger than
\pkg{RSiena}'s default, which gives phase-3 standard errors and convergence
ratios with less Monte Carlo error; on the CPU build a smaller schedule (for
example \code{n3 = 1000}) finishes correspondingly faster. The estimator
drives the CPU backend of the engine; \code{backend = "gpu"} is reserved for
a build that links the CUDA backend.

The simulated end networks returned with \code{return_nets = TRUE} are the
states at the end of each period (in conditional estimation, when the
observed distance is reached), before the statistic masks for missing,
structural and composition-change ties are applied.
}
\seealso{
\code{\link{cusna_fit}}, \code{\link{cusna_effect}}, \code{\link{saom_data}},
\code{\link{mom_estimate_multinet}}, \code{\link{cusna_fran}}.
}
\examples{
# The s50 panel: three waves of friendship among 50 girls.
dat <- saom_data(s50$waves)
effects <- list(cusna_effect("density"), cusna_effect("recip"),
                cusna_effect("transTrip"))
\donttest{
# A short schedule keeps the example quick; use the defaults in practice.
ctl <- mom_control(n1 = 200, nsub = 2, n2 = c(10, 10), batch2 = 50,
                   nD = 200, n3 = 300, seed = 1)
fit <- mom_estimate(dat, effects, conditional = TRUE, control = ctl,
                    return_nets = TRUE)
summary(fit)

# The phase-3 networks: indegree distribution of the simulated third wave.
dim(fit$sims[[2]])
indeg <- cusna_gof_distribution(fit$sims[[2]], kind = 0L, levls = 0:6)
colMeans(indeg)
}
}
