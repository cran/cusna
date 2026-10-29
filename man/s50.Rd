\name{s50}
\alias{s50}
\docType{data}
\title{Friendship and alcohol use of 50 girls (Teenage Friends and Lifestyle Study)}
\description{
A three-wave panel of friendship nominations among 50 randomly chosen girls
from the Teenage Friends and Lifestyle Study, with their alcohol use. It is
the standard small panel for stochastic actor-oriented models (no missing
values, quick to fit) and is also distributed with \pkg{RSiena} as
\code{s501}, \code{s502}, \code{s503} and \code{s50a}; the examples and
vignettes of this package use it.
}
\usage{s50}
\format{
A list with two elements:
\describe{
  \item{\code{waves}}{a list of three \eqn{50 \times 50}{50 x 50} integer
    adjacency matrices (0/1, zero diagonal), the friendship networks at time
    points 1--3; entry \code{[i, j]} is 1 when girl \code{i} names girl
    \code{j} as a friend.}
  \item{\code{alcohol}}{a \eqn{50 \times 3}{50 x 3} integer matrix of alcohol
    use at the three time points, coded 1--5 (higher values indicate higher
    consumption).}
}
}
\source{
\url{https://www.stats.ox.ac.uk/~snijders/siena/s50_data.htm}; identical to
the \pkg{RSiena} datasets \code{s501}, \code{s502}, \code{s503} and
\code{s50a}.
}
\references{
West, P. and Sweeting, H. (1995). \emph{Background Rationale and Design of the
West of Scotland 11-16 Study}. Working Paper No. 52, MRC Medical Sociology
Unit Glasgow.
}
\examples{
str(s50)
dat <- saom_data(s50$waves)
dat
# ties per wave
vapply(s50$waves, sum, numeric(1))
}
\keyword{datasets}
