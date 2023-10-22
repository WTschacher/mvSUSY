\name{mvsusy}
\alias{mvsusy}
\title{ Multivariate Surrogate Synchrony }
\description{
  Multivariate Surrogate Synchrony estimates the synchrony within datasets that contain more than two time series. Multivariate SUSY was adapted and extended from SUCO and SUSY with respect to implementing the surrogate controls.
}
\usage{
mvsusy(x, segment, Hz, method=c("lambda_max","omega"), max_pseudo=1000, seed=1)
}
\arguments{
  \item{x}{ A data.frame of numeric columns. }
  \item{segment}{ Integer, size in seconds. Must not be smaller than \code{2 * maxlag}, must not be larger than half the the time series (\code{nrow(x)/2}). }
  \item{Hz}{ Integer, frames per second (sampling rate). }
  \item{method}{ Character, either \code{"lambda_max"} or \code{"omega"}. }
  \item{max_pseudo}{ Numeric, maximum number of permutations. }
  \item{seed}{ Numeric, random seed used to generate permutations. }
}
\details{
  Segments are non-overlapping, and the number of segments that fit into the time series may have a remainder (usually a few seconds at the end of the time series), which is not considered. Rows having missing values are removed.
}
\value{
  Object of class \code{mvsusy} is returned.
}
\seealso{\code{\link{plot.mvsusy}}, \code{\link{as.data.frame.mvsusy}}, \code{\link{print.mvsusy}}}
\examples{
set.seed(1)
data = as.data.frame(replicate(5, sample(10, 500, TRUE)))

## compute mvSUSY (lambda_max method)
res = mvsusy(data, segment=10, Hz=10, max_pseudo=100)
res
plot(res, type="eigenvalue")

## omega method
res = mvsusy(data, segment=10, Hz=10, max_pseudo=100, method="omega")
res
plot(res, type="density")

## export to flat file via data.frame and write.csv
df = as.data.frame(res)
df
}
\keyword{ data }
