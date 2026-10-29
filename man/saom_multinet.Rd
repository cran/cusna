\name{saom_multinet_data}
\alias{saom_multinet_data}
\alias{mom_estimate_multinet}
\title{Multi-network co-evolution}
\usage{
saom_multinet_data(waves_list, composition = NULL, cc_option = 1)

mom_estimate_multinet(data, effects_by_net, backend = c("cpu", "gpu"),
                      control = mom_control(), verbose = FALSE)
}
\arguments{
\item{waves_list}{a list with one entry per dependent network; each entry is
  itself a list of wave matrices (as in \code{\link{saom_data}}). All networks
  share the actor set and the number of waves.}
\item{composition}{optional composition change applied to the shared actor
  set (see \code{\link{saom_data}}).}
\item{cc_option}{the composition-change option (1--3).}
\item{data}{a \code{cusna_multinet_data} object.}
\item{effects_by_net}{a list with one effects list per network. Within-network
  effects are ordinary \code{\link{cusna_effect}} objects; cross-network
  effects (\code{"crprod"}, \code{"crprodRecip"}, \code{"crprodMutual"},
  \code{"from"}, \code{"to"}) additionally take \code{net_ref}, the 1-based
  index of the source network.}
\item{backend}{\code{"cpu"} (default) or \code{"gpu"}.}
\item{control}{tuning parameters from \code{\link{mom_control}}.}
\item{verbose}{logical; print per-phase progress.}
}
\value{
\code{saom_multinet_data} returns a \code{cusna_multinet_data} object.
\code{mom_estimate_multinet} returns a \code{\link{cusna_fit}}; its parameters
are laid out per network (rates then effects), following RSiena's convention.
}
\description{
Several dependent networks that co-evolve over one actor set, with
cross-network effects. \code{saom_multinet_data} builds the panel and
\code{mom_estimate_multinet} fits it by the method of moments on the native
multi-network simulator (unconditional only, as in RSiena for more than one
dependent variable).
}
\details{
Cross-network statistics read the other network at its period-start state;
the dynamics read the current co-evolving state. Multi-network models do not
support behavior co-evolution, interactions, or degree bounds.
}
\seealso{
\code{\link{mom_estimate}}, \code{\link{cusna_effect}}.
}
\examples{
# Friendship (the s50 waves) co-evolving with a second, synthetic network on
# the same 50 actors.
set.seed(11)
second <- lapply(s50$waves, function(w) {
  m <- w * (matrix(runif(2500), 50, 50) < 0.5); diag(m) <- 0L; m
})
dat <- saom_multinet_data(list(s50$waves, second))
dat
\donttest{
# A short schedule keeps the example quick; use the defaults in practice.
fit <- mom_estimate_multinet(dat, effects_by_net = list(
  list(cusna_effect("density"), cusna_effect("recip")),
  list(cusna_effect("density"), cusna_effect("crprod", net_ref = 1))
), control = mom_control(n1 = 150, nsub = 1, n2 = 10, batch2 = 40,
                         nD = 150, n3 = 200, seed = 5))
summary(fit)
}
}
