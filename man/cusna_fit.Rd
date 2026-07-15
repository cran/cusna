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
\code{D} (derivative matrix), and \code{msf} (moment covariance).
}
\description{
The object returned by \code{\link{mom_estimate}} and
\code{\link{mom_estimate_multinet}}, with the usual R accessors.
}
\seealso{
\code{\link{mom_estimate}}.
}
\examples{
\dontrun{
fit <- mom_estimate(dat, effects)
print(fit)
summary(fit)
coef(fit)
vcov(fit)
as.data.frame(fit)
}
}
