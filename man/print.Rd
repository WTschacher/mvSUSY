\name{print.mvsusy}
\alias{print.mvsusy}
\title{ mvsusy print method }
\description{
  Prints information about an \code{mvsusy} object.
}
\usage{
  \method{print}{mvsusy}(x, \dots)
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{\dots}{ Extra arguments passed to \code{print.data.frame} method. }
}
\value{
  Returns \code{x} invisibly. Display output to console as a side effect.
}
\seealso{\code{\link{mvsusy}}}
\examples{
set.seed(1)
data = as.data.frame(replicate(5, sample(10, 5000, TRUE)))
res = mvsusy(data, segment=10L, Hz=10L)
res
}
\keyword{ data }
