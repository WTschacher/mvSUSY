\name{plot.mvsusy}
\alias{plot.mvsusy}
\title{ mvsusy plot method }
\description{
  Generate plot for a \code{mvsusy} object.
}
\usage{
  \method{plot}{mvsusy}(x,
    type=c("eigenvalue","density","free scale","segment-wise","time series"),
    \dots,
    plotly)
}
\arguments{
  \item{x}{ A \code{mvsusy} object. }
  \item{type}{ Character to specify type of the plot, one of \emph{eigenvalue, density, free scale, segment-wise, time series}. }
  \item{\dots}{ Ignored. }
  \item{plotly}{ Logical. Passing \code{FALSE} will force \emph{time series} plot not to use \code{plotly} even if installed. }
}
\details{
  Eigenvalue plot type works only for mvSUSY computed using \emph{lambda_max} method. Time series plot will by default use \code{plotly} package, if installed, to provide interactive plot.
}
\value{
  Returns \code{NULL} invisibly. Generate plot as a side effect.
}
\seealso{\code{\link{mvsusy}}}
\examples{
set.seed(1)
data = as.data.frame(replicate(5, sample(10, 500, TRUE)))
res = mvsusy(data, segment=10, Hz=10, max_pseudo=100)
plot(res, type="eigenvalue")
plot(res, type="density")
plot(res, type="free scale")
plot(res, type="segment-wise")
plot(res, type="time series")
plot(res, type="time series", plotly=FALSE)
}
\keyword{ data }
