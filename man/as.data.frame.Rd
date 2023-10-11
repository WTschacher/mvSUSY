\name{as.data.frame.mvsusy}
\alias{as.data.frame.mvsusy}
\title{ mvsusy to data.frame conversion method }
\description{
  Turns \code{mvsusy} class object into a \code{data.frame}.
}
\usage{
\method{as.data.frame}{mvsusy}(x, row.names=NULL, optional=FALSE, \dots)
}
\arguments{
  \item{x}{ A \code{mvsusy} object. }
  \item{row.names}{ Ignored, only for consistency to generic \code{as.data.frame} method. }
  \item{optional}{ Ignored, only for consistency to generic \code{as.data.frame} method. }
  \item{\dots}{ Ignored. }
}
\value{
  Returns \code{data.frame}.
}
\seealso{\code{\link{mvsusy}}}
\examples{
data = as.data.frame(replicate(5, sample(10, 500, TRUE)))
res = mvsusy(data, segment=10L, Hz=10L, max_pseudo=100)
df = as.data.frame(res)
df
}
\keyword{ data }
