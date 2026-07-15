\name{tergm_mple}
\alias{tergm_mple}
\alias{print.tergm_mple}
\title{Temporal ERGM by bootstrap pseudo-likelihood (btergm style)}
\usage{
tergm_mple(nets, mutual = TRUE, ttriple = TRUE, memory = TRUE,
           delrecip = TRUE, R = 200L, level = 0.95, seed = 1234L)
}
\arguments{
\item{nets}{A list of >= 2 directed n*n adjacency matrices (0/1, zero diagonal),
  in time order.}
\item{mutual, ttriple}{Include the cross-sectional mutual / transitive-triple terms
  (from the native change-statistic matrix).}
\item{memory, delrecip}{Include the temporal terms (lagged tie / lagged reciprocity).}
\item{R}{Bootstrap replications for the confidence intervals.}
\item{level}{Confidence level for the bootstrap intervals.}
\item{seed}{Integer RNG seed for the bootstrap.}
}
\value{
A \code{tergm_mple} object with \code{coef}, \code{ci} (bootstrap interval),
\code{boot} (replicate coefficients), \code{nobs} (pooled dyads) and
\code{transitions}.
}
\description{
Fits a temporal ERGM to a sequence of directed networks by pooled maximum
pseudo-likelihood -- the btergm estimator (Desmarais & Cranmer 2012). For every
transition the native engine builds the 0->1 change statistics of every dyad of
the current network (edges / mutual / transitive triple) and the temporal terms
are read off the lagged network (memory = lagged tie, delrecip = lagged
reciprocity); the rows are pooled across transitions and fit by weighted logistic
regression, with bootstrap confidence intervals. The design matrix is host C++ (no
Python; works in the CPU-only build).
}
\examples{
\donttest{
set.seed(1)
nets <- lapply(1:3, function(t) {
  m <- matrix(as.integer(runif(20 * 20) < 0.2), 20, 20)
  diag(m) <- 0L
  m
})
fit <- tergm_mple(nets, R = 50L)
fit
}
}
