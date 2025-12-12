\name{residuals-methods}
\alias{residuals-methods}

\docType{methods}

\alias{residuals}
\alias{residuals,fGARCH-method}


\title{Extract GARCH model residuals}

\description{
  Extracts residuals from a fitted GARCH object.
}

\usage{
\S4method{residuals}{fGARCH}(object, standardize = FALSE)
}

\arguments{
  \item{object}{

    an object of class \code{"fGARCH"} as returned by
    \code{\link{garchFit}}.

  }
  \item{standardize}{

    a logical, indicating if the residuals should be standardized.

  }
}

\details{

  The \code{"fGARCH"} method extracts the \code{@residuals} slot from an
  object of class \code{"fGARCH"} as returned by the function
  \code{garchFit} and optionally standardizes them, using conditional
  standard deviations.

}

\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\seealso{
  \code{\link{fitted}},
  \code{\link{predict}},
  \code{\link{garchFit}},
  class \code{\linkS4class{fGARCH}},
}

\examples{
stopifnot(require("timeSeries"))
## Swiss Pension fund Index
data(LPP2005REC, package = "timeSeries")
x <- as.timeSeries(LPP2005REC)

## Fit LPP40 Bechmark:
fit <- garchFit(LPP40 ~ garch(1, 1), data = 100*x, trace = FALSE)
fit

fitted <- fitted(fit)
head(fitted)
class(fitted)

res <- residuals(fit)
head(res)
class(res)
}

\keyword{models}
