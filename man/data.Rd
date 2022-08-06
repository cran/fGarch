\name{TimeSeriesData}

\alias{TimeSeriesData}

\alias{dem2gbp}
\alias{sp500dge}

\title{Time Series Data Sets} 

\description{
  
  Data sets used in the examples of the timeSeries packages, including
  DEM/GBP foreign exchange rates dataset and data on SP500 index.

}
\format{

  \code{dem2gbp} is a data frame with one column \code{"DEM2GBP"} and
  1974 rows (observations).
  
  \code{sp500dge} is a data frame with one column \code{"SP500DGE"} and
  17055 rows (observations).

}
\details{
  The data represent retuns. No further details have been recorded.
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
