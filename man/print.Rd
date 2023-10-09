\name{print.mvsusy}
\alias{print.mvsusy}
\title{ mvsusy print method }
\description{
  Prints information about an \code{mvsusy} object.
}
\usage{
  \method{print}{mvsusy}(x, \dots, name_data="", piece="", participant="")
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{\dots}{ Extra arguments passed to \code{print.data.frame} method. }
  \item{name_data}{ Character. }
  \item{piece}{ Character. }
  \item{participant}{ Character. }
}
\value{
  Returns \code{x} invisibly. Display output to console as a side effect.
}
\seealso{\code{\link{mvsusy}}}
\examples{
data = as.data.frame(replicate(5, sample(10, 500, TRUE)))
res = mvsusy(data, segment=10L, Hz=10L, max_pseudo=100)
res
}
\keyword{ data }
