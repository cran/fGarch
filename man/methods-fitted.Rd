\name{fitted-methods}
\alias{fitted-methods}

\docType{methods}

\alias{fitted}
\alias{fitted,fGARCH-method}


\title{Extract GARCH model fitted values}

\description{

  Extracts fitted values from a fitted GARCH object.

}


\section{Methods}{
  Methods for \code{fitted} defined in package \pkg{fGarch}:
  \describe{

    \item{object = "fGARCH"}{
      Extractor function for fitted values.
    }

  }
}


\details{

    \code{fitted()} is a generic function which extracts fitted values
    from objects returned by modeling functions.  \code{methods(fitted)}
    will show the methods available in the current \R session.

    The method for \code{"fGARCH"} objects extracts the \code{@fitted}
    value slot from an object of class \code{"fGARCH"} as returned by
    the function \code{garchFit}.

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
  \code{\link{predict}},
  \code{\link{residuals}},
  \code{\link{garchFit}},
  class \code{\linkS4class{fGARCH}},

  \code{\link{plot}}
}

\examples{
stopifnot(require("timeSeries"))
## Swiss Pension fund Index -
   x = as.timeSeries(data(LPP2005REC))

## garchFit -
   # Fit LPP40 Bechmark:
   fit = garchFit(LPP40 ~ garch(1, 1), data = 100*x, trace = FALSE)
   fit

## fitted -
   # Fitted values are now a "timeSeries" oject:
   fitted = fitted(fit)
   head(fitted)
   class(fitted)

## slot -
   # The slot contains a numeric Vector:
   fitted = slot(fit, "fitted")
   head(fitted)
   class(fitted)
}

\keyword{models}

