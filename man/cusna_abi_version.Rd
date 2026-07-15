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
\code{cusna_has_cuda()} returns \code{TRUE} when the package was built with the
GPU path (a CUDA toolkit was detected at configure time) and a usable device is
present; the CRAN CPU-only build always returns \code{FALSE}.
}
