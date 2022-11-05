\name{volatility-methods}

\docType{methods}

\alias{volatility}% S3 method now
\alias{volatility.fGARCH}% S3 method now


\title{Extract GARCH model volatility}

\description{

  Extracts volatility from a fitted GARCH object. 

}

\usage{
volatility(object, \dots)
\S3method{volatility}{fGARCH}(object, type = c("sigma", "h"), \dots)
}

\arguments{
  \item{object}{

    an object of class \code{"fGARCH"} as returned by
    \code{\link{garchFit}()}.

  }
  \item{type}{

    a character string denoting if the conditional standard deviations
    \code{"sigma"} or the variances \code{"h"} should be returned.

  }
  \item{\dots}{additional arguments to be passed.}
}

\details{

  \code{volatility} is an S3 generic function for computation of
  volatility, see \code{link[fBasics]{volatility}} for the default
  method.

  The method for \code{"fGARCH"} objects, described here, extracts the
  volatility from slot \code{@sigma.t} or \code{@h.t} of an
  \code{"fGARCH"} object usually obtained from the function
  \code{\link{garchFit}()}.

  The class of the returned value depends on the input to the function
  \code{garchFit} who created the object. The returned value is always
  of the same class as the input object to the argument \code{data} in
  the function \code{garchFit}, i.e. if you fit a \code{"timeSeries"}
  object, you will get back from the function \code{fitted} also a
  \code{"timeSeries"} object, if you fit an object of class
  \code{"zoo"}, you will get back again a \code{"zoo"} object. The same
  holds for a \code{"numeric"} vector, for a \code{"data.frame"}, and
  for objects of class \code{"ts", "mts"}.

  In contrast, the slot itself always contains a numeric vector,
  independently of the class of the input data input, i.e. the function
  call \code{slot(object, "fitted")} will return a numeric vector.

}

\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\note{
  
  (GNB) Contrary to the description of the returned value of the
  \code{"fGARCH"} method, it is always \code{"numeric"}.

  TODO: either implement the documented behaviour or fix the
  documentation.

}

\section{Methods}{
  Methods for \code{volatility} defined in package \pkg{fGarch}:
  \describe{
    \item{object = "fGARCH"}{
      Extractor function for volatility or standard deviation from
      an object of class \code{"fGARCH"}.
    }

  }
}

\seealso{
  \code{\link{garchFit}},
  class \code{\linkS4class{fGARCH}}
}

\examples{
## Swiss Pension fund Index -
   stopifnot(require("timeSeries")) # need package 'timeSeries'
   x = as.timeSeries(data(LPP2005REC, package = "timeSeries"))

## garchFit
   fit = garchFit(LPP40 ~ garch(1, 1), data = 100*x, trace = FALSE)
   fit

## volatility -
   # Standard Deviation:
   vola = volatility(fit, type = "sigma")
   head(vola)
   class(vola)
   # Variance:
   vola = volatility(fit, type = "h")
   head(vola)
   class(vola)

## slot -
   vola = slot(fit, "sigma.t")
   head(vola)
   class(vola)
   vola = slot(fit, "h.t")
   head(vola)
   class(vola)
}

\keyword{models}
