\name{ergm_mple}
\alias{ergm_mple}
\title{Maximum pseudo-likelihood estimate (directed edges/mutual/ttriple demo)}
\usage{
ergm_mple(x)
}
\arguments{
\item{x}{Integer/logical n*n directed adjacency matrix (0/1, zero diagonal).}
}
\value{
A named numeric vector of pseudo-MLE coefficients (\code{edges}, \code{mutual},
\code{transitive}).
}
\description{
A logistic regression of tie presence on the dyadic change statistics the native
engine computes (Strauss & Ikeda 1990). This is the MPLE demo for the three
directed terms edges / mutual / transitive-triple; full MCMC-MLE stays with the
\code{ergm} package (framing M-3).
}
\examples{
set.seed(1)
a <- matrix(as.integer(runif(20 * 20) < 0.2), 20, 20)
diag(a) <- 0L
ergm_mple(a)
}
