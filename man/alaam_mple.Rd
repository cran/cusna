\name{alaam_mple}
\alias{alaam_mple}
\alias{print.alaam_mple}
\title{Autologistic actor attribute model by maximum pseudo-likelihood (ALAAM)}
\usage{
alaam_mple(y, net, cov = NULL, contagion = "out", directed = FALSE)
}
\arguments{
\item{y}{Binary actor attribute (0/1), length n -- the dependent variable.}
\item{net}{Integer/logical n*n adjacency of the fixed (exogenous) network.}
\item{cov}{Optional actor covariate(s): a length-n vector or an n*p matrix,
  entered as change statistics (the covariate value at each actor).}
\item{contagion}{Which contagion change statistic to include: \code{"out"}
  (out-neighbours with \code{y = 1}, the default and the neighbour count when
  undirected), \code{"in"}, \code{"both"}, or \code{FALSE} for none.}
\item{directed}{Logical: is \code{net} directed? (documentation only; the
  contagion columns are selected by \code{contagion}).}
}
\value{
An \code{alaam_mple} object with \code{coef}, \code{se}, \code{z}, \code{pval}
(from the logistic regression), \code{glm} (the fitted model) and \code{n}.
}
\description{
Fits an ALAAM (Robins, Pattison & Elliott 2001) -- a model for a binary actor
attribute at a fixed, exogenous network -- by maximum pseudo-likelihood, i.e. a
logistic regression of the attribute on the native contagion change statistic
(the number of network neighbours that also carry the attribute, a masked
popcount) plus any actor covariates. The dependent variable is the attribute,
not a tie; the network enters through the same bit-kernel primitive applied to
(network-row AND attribute-vector). Best suited to an undirected network, where
the contagion is the neighbour count with \code{y = 1}. Host C++ (no Python;
works in the CPU-only build).
}
\examples{
\donttest{
set.seed(1)
net <- matrix(as.integer(runif(30 * 30) < 0.1), 30, 30)
net <- ((net + t(net)) > 0) * 1L        # undirected
diag(net) <- 0L
y <- as.integer(runif(30) < 0.4)
fit <- alaam_mple(y, net)
fit
}
}
