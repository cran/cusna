\name{cusna_abi_version}
\alias{cusna_abi_version}
\alias{cusna_has_cuda}
\title{Native cusna engine: ABI version and CUDA availability}
\usage{
cusna_abi_version()

cusna_has_cuda()
}
\value{
An integer (\code{cusna_abi_version}) or logical (\code{cusna_has_cuda}).
}
\description{
\code{cusna_abi_version()} returns the C ABI version compiled into the package.
\code{cusna_has_cuda()} returns \code{TRUE} only for a build that links the
CUDA backend of the engine and finds a usable device; the R package builds the
CPU backend, so it returns \code{FALSE}.
}
\examples{
cusna_abi_version()
cusna_has_cuda()
}
