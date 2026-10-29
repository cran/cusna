\name{saom_data}
\alias{saom_data}
\title{Create a SAOM data panel}
\usage{
saom_data(waves, behavior = NULL, composition = NULL, cc_option = 1)
}
\arguments{
\item{waves}{a list of at least two square adjacency matrices, one per
  observation wave, all on the same actor set. Missing ties are \code{NA}
  (or the code 9); structural zeros and ones are the codes 10 and 11.}
\item{behavior}{optional co-evolving dependent behavior: a numeric matrix
  with one row per actor and one column per wave (\code{NA} for missing).}
\item{composition}{optional composition change (joiners/leavers): a list
  with one entry per actor, each a numeric vector \code{c(start, end)} (or
  several stacked pairs) of 1-based wave numbers the actor is present for.
  Estimation is then unconditional, as in RSiena.}
\item{cc_option}{the \code{sienaCompositionChange} option (1--3) controlling
  how absent actors' entries are treated.}
}
\value{
An object of class \code{cusna_data} with fields \code{n} (actors),
\code{n_waves}, and \code{n_periods}, carrying the prepared engine data in
\code{$internal}.
}
\description{
Builds a longitudinal network panel for estimation, applying RSiena's
imputation and per-period masks (first-wave missings to 0, later waves
carried forward; structural values fixed; composition-change gating) in
native R code -- no Python. The construction is validated bit-for-bit
against the reference implementation on public datasets.
}
\seealso{
\code{\link{mom_estimate}}, \code{\link{cusna_effect}},
\code{\link{saom_multinet_data}}.
}
\examples{
dat <- saom_data(s50$waves)
dat

# with the co-evolving alcohol use as a dependent behavior
saom_data(s50$waves, behavior = s50$alcohol)
}
