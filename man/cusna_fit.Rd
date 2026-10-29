\name{cusna_fit}
\alias{cusna_fit}
\alias{print.cusna_fit}
\alias{summary.cusna_fit}
\alias{coef.cusna_fit}
\alias{vcov.cusna_fit}
\alias{as.data.frame.cusna_fit}
\title{Fitted SAOM model}
\usage{
\method{print}{cusna_fit}(x, ...)

\method{summary}{cusna_fit}(object, ...)

\method{coef}{cusna_fit}(object, ...)

\method{vcov}{cusna_fit}(object, ...)

\method{as.data.frame}{cusna_fit}(x, ...)
}
\arguments{
\item{x, object}{a \code{cusna_fit} object.}
\item{...}{ignored.}
}
\value{
\code{coef} returns the named parameter estimates; \code{vcov} the
delta-method covariance matrix (whose diagonal squares to the reported
standard errors); \code{as.data.frame} a data frame with columns
\code{effect}, \code{estimate}, \code{se}, and \code{t_conv}. \code{print}
and \code{summary} report the fit and return the object invisibly.

A \code{cusna_fit} is a list with elements \code{coefficients}, \code{se},
\code{tconv}, \code{estimates} (data frame), \code{rates} (conditional rate
estimates, or \code{NULL}), \code{tconv_max}, \code{n_sims}, \code{wall_time},
\code{backend}, \code{conditional}, \code{targets}, \code{sim_means},
\code{D} (derivative matrix), \code{msf} (moment covariance), and
\code{sims} / \code{sims_behavior}: the phase-3 end networks (one
\code{n3 x n x n} integer array per period) and, for co-evolution, end
behavior (one \code{n3 x n} matrix per period) when the fit was made with
\code{return_nets = TRUE}, otherwise \code{NULL}.
}
\description{
The object returned by \code{\link{mom_estimate}} and
\code{\link{mom_estimate_multinet}}, with the usual R accessors.
}
\seealso{
\code{\link{mom_estimate}}.
}
\examples{
\donttest{
# Conditional method-of-moments fit of the s50 panel (short schedule).
fit <- mom_estimate(
  saom_data(s50$waves),
  effects = list(cusna_effect("density"), cusna_effect("recip")),
  conditional = TRUE,
  control = mom_control(n1 = 200, nsub = 1, n2 = 15, batch2 = 50,
                        nD = 200, n3 = 300, seed = 2))
print(fit)
summary(fit)
coef(fit)
vcov(fit)
as.data.frame(fit)
}
}
