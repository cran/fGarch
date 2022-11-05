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

    an object of class \code{"fGARCH"} as returned from
    \code{\link{garchFit}}.

  }
  \item{standardize}{

    a logical, indicating if the residuals should be standardized.

  }
}

\details{

  \code{\link{residuals}} is a generic function which extracts residual
  values from objects returned by modeling functions.

  The \code{"fGARCH"} method extracts the \code{@residuals} slot from an
  object of class \code{"fGARCH"} as returned by the function
  \code{garchFit} and optionally standardizes the residuals, using
  conditional standard deviations.

    The class of the returned value depends on the input to the
    function \code{garchFit} who created the object. The returned
    value is always of the same class as the input object to the
    argument \code{data} in the function \code{garchFit}, i.e. if
    you fit a \code{"timeSeries"} object, you will get back from
    the function \code{fitted} also a \code{"timeSeries"} object,
    if you fit an object of class \code{"zoo"}, you will get back
    again a \code{"zoo"} object. The same holds for a \code{"numeric"}
    vector, for a \code{"data.frame"}, and for objects of class
    \code{"ts", "mts"}.

    In contrast, the slot itself returns independent of the class
    of the data input always a numeric vector, i.e. the function
    call r\code{slot(object, "fitted")} will return a numeric vector.

}

\note{
  
  (GNB) Contrary to the description of the returned value of the
  \code{"fGARCH"} method, it is always \code{"numeric"}.

  TODO: either implement the documented behaviour or fix the
  documentation.

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
## Swiss Pension fund Index -
   stopifnot(require("timeSeries"))
   x = as.timeSeries(data(LPP2005REC))

## garchFit
   fit = garchFit(LPP40 ~ garch(1, 1), data = 100*x, trace = FALSE)
   fit

## residuals -
   res = residuals(fit)
   head(res)
   class(res)

## slot -
   identical(res, slot(fit, "residuals"))
}

\keyword{models}
