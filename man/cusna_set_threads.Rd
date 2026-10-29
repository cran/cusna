\name{cusna_set_threads}
\alias{cusna_set_threads}
\title{Set the OpenMP thread count for the native CPU backend}
\usage{
cusna_set_threads(n)
}
\arguments{
\item{n}{Requested thread count (coerced to an integer, at least 1).}
}
\value{
Invisibly, the previous maximum thread count (\code{0L} if OpenMP is not
compiled in).
}
\description{
Explicitly caps the number of threads the CPU backend uses for its
OpenMP-parallel batch loops (SAOM simulation, ERGM sampling and MPLE).
Unlike \code{Sys.setenv(OMP_NUM_THREADS = ...)}, this takes effect
immediately via the OpenMP runtime API and does not depend on when the
environment variable would otherwise be read. Has no effect on a build
without OpenMP (e.g. the default macOS toolchain without \code{libomp}).
}
\examples{
old <- cusna_set_threads(1L)
cusna_set_threads(old)
}
