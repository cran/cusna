\name{cusna_effect}
\alias{cusna_effect}
\alias{cusna_beh_effect}
\alias{cusna_rate_effect}
\alias{cusna_interaction}
\title{Specify SAOM model effects}
\usage{
cusna_effect(name, covariate = NULL, parameter = NULL, dyn = FALSE,
             type = c("eval", "creation", "endow"), net_ref = NULL)

cusna_beh_effect(name)

cusna_rate_effect(name, covariate = NULL)

cusna_interaction(components)
}
\arguments{
\item{name}{effect name. For \code{cusna_effect}, a structural effect
  (e.g. \code{"density"}, \code{"recip"}, \code{"transTrip"}, the
  gwesp/gwdsp families), a covariate effect (\code{"egoX"}, \code{"altX"},
  \code{"simX"}, \code{"sameX"}, ...), or a dyadic-covariate effect
  (\code{"X"}). For \code{cusna_beh_effect}, one of \code{"linear"},
  \code{"quad"}, \code{"avAlt"}, \code{"avSim"}, \code{"totSim"}. For
  \code{cusna_rate_effect}, one of \code{"RateX"}, \code{"outRate"},
  \code{"inRate"}.}
\item{covariate}{covariate values. Actor covariate: a numeric vector of
  length \code{n} (constant) or an \code{n} by (waves) matrix (changing).
  Dyadic covariate: an \code{n} by \code{n} matrix, or \code{n} by \code{n}
  by (periods) for a changing one. Centered by the engine.}
\item{parameter}{optional internal effect parameter (RSiena's
  \code{parameter}); \code{NULL} uses the effect's default. The engine
  truncates this to an integer for gwesp/gwdsp/threshold effects, matching
  RSiena.}
\item{dyn}{logical; if \code{TRUE} the covariate is the co-evolving dependent
  behavior (choice dynamics use the current values, statistics use the
  start-of-period values).}
\item{type}{effect type: \code{"eval"} (default), \code{"creation"}, or
  \code{"endow"} (endowment).}
\item{net_ref}{for cross-network effects in multi-network models (see
  \code{\link{mom_estimate_multinet}}), the 1-based index of the source
  network; \code{NULL} otherwise.}
\item{components}{for \code{cusna_interaction}, a length-2 or length-3 vector
  of 1-based positions of the component effects in the effects list. The
  components must appear earlier in that list.}
}
\value{
A lightweight specification object of class \code{cusna_effect},
\code{cusna_beh_effect}, \code{cusna_rate_effect}, or
\code{cusna_interaction}, to be passed to \code{\link{mom_estimate}}.
}
\description{
Constructors for the terms of a stochastic actor-oriented model. They return
lightweight R objects; \code{\link{mom_estimate}} translates them into the
native engine's effect descriptors at fit time.
}
\seealso{
\code{\link{mom_estimate}}, \code{\link{saom_data}}.
}
\examples{
effects <- list(
  cusna_effect("density"),
  cusna_effect("recip"),
  cusna_effect("transTrip")
)

# A covariate effect and a user interaction of effects 1 and 2.
alcohol <- runif(50)
list(
  cusna_effect("egoX", covariate = alcohol),
  cusna_effect("altX", covariate = alcohol),
  cusna_interaction(c(1, 2))
)

# Behavior and rate effects.
cusna_beh_effect("linear")
cusna_rate_effect("outRate")
}
