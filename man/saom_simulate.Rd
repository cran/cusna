\name{saom_simulate}
\alias{saom_simulate}
\alias{print.cusna_sim}
\title{Simulate a SAOM at given parameters}
\usage{
saom_simulate(data, effects, theta, nsim = 100, beh_effects = NULL,
              conditional = FALSE, maxdegree = 0, seed = NULL,
              return_nets = TRUE)

\method{print}{cusna_sim}(x, ...)
}
\arguments{
\item{data}{a \code{cusna_data} object from \code{\link{saom_data}}.}
\item{effects}{a non-empty list of \code{\link{cusna_effect}} (and
  \code{\link{cusna_interaction}}) objects, as for \code{\link{mom_estimate}}.}
\item{theta}{the parameter vector, in the layout of \code{coef()} of the
  matching fit (unconditional: the basic rates per period, then the
  effects; conditional: the effects only), or a \code{\link{cusna_fit}}.}
\item{nsim}{number of independent chains simulated for each period.}
\item{beh_effects}{optional list of \code{\link{cusna_beh_effect}} objects for
  network--behavior co-evolution.}
\item{conditional}{logical; if \code{TRUE}, each period runs until the
  observed distance is reached, as in conditional estimation.}
\item{maxdegree}{optional maximum out-degree constraint (0 = none).}
\item{seed}{optional seed; when given, \code{set.seed(seed)} is called before
  the per-period simulation seeds are drawn from R's random number stream.}
\item{return_nets}{logical; return the simulated end networks.}
\item{x}{a \code{cusna_sim} object.}
\item{...}{ignored.}
}
\value{
A \code{cusna_sim} list with elements \code{stats} (an \code{nsim}-by-\eqn{p}
matrix of the moment statistics, summed over periods, with the observed
\code{targets} alongside), \code{theta}, \code{durations} (conditional
simulation: an \code{nsim}-by-periods matrix of the simulated durations,
whose means estimate the rates), \code{nets} (when \code{return_nets}: one
\code{nsim x n x n} integer array of end networks per period, chains first),
\code{behavior} (co-evolution: one \code{nsim x n} matrix of end behavior per
period), \code{nsim}, \code{n_periods} and \code{conditional}. The
\code{print} method compares the simulated means of the statistics with the
observed targets.
}
\description{
Forward-simulates every period of the panel, starting each period from the
observed wave, at the parameter vector \code{theta} on the native engine.
The statistics are the moment statistics that \code{\link{mom_estimate}}
matches, so a simulation at the estimated parameters reproduces the targets
on average; the end networks serve goodness-of-fit checks and
counterfactual scenarios. This is the SAOM counterpart of
\code{\link{ergm_simulate}}.
}
\details{
Chains are independent and use counter-based random numbers, so a
simulation with a given seed is identical for any number of threads
(\code{\link{cusna_set_threads}}).
}
\seealso{
\code{\link{mom_estimate}} (whose \code{return_nets} keeps the phase-3
networks of an estimation run), \code{\link{cusna_gof_distribution}}.
}
\examples{
dat <- saom_data(s50$waves)
effects <- list(cusna_effect("density"), cusna_effect("recip"),
                cusna_effect("transTrip"))

# Conditional simulation at parameter values close to the s50 estimates.
sim <- saom_simulate(dat, effects, theta = c(-2.7, 2.4, 0.63), nsim = 200,
                     conditional = TRUE, seed = 1)
sim
colMeans(sim$durations)          # simulated durations ~ the rate estimates

# Outdegree distribution of the simulated third wave
outdeg <- cusna_gof_distribution(sim$nets[[2]], kind = 1L, levls = 0:6)
round(colMeans(outdeg), 1)
}
