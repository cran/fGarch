\name{fGarchData}
\alias{fGarchData}

\alias{dem2gbp}
\alias{sp500dge}

\concept{exchange rates data}
\concept{SP500 data}


\title{Time series datasets} 

\description{
  
  Datasets used in the examples, including DEM/GBP foreign exchange
  rates and data on SP500 index.

}
\format{

  \code{dem2gbp} is a data frame with one column \code{"DEM2GBP"} and
  1974 rows (observations).
  
  \code{sp500dge} is a data frame with one column \code{"SP500DGE"} and
  17055 rows (observations).

}
\details{
  The data represent retuns. No further details have been recorded.

  Further datasets are available in the packages that \pkg{fGarch}
  imports, see \code{\link[fBasics]{fBasicsData}} and
  \code{\link[timeSeries]{TimeSeriesData}}.
}

\seealso{
  \code{data(package = "fBasics")} and
  \code{data(package = "timeSeries")} for related datasets

}
\examples{
data(dem2gbp)
head(dem2gbp)
tail(dem2gbp)
str(dem2gbp)
plot(dem2gbp[[1]])

data(sp500dge)
head(sp500dge)
tail(sp500dge)
str(sp500dge)
plot(sp500dge[[1]])
}
\keyword{datasets}
\keyword{ts}
