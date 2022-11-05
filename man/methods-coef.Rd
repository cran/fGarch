\name{coef-methods}
\alias{coef-methods}

\docType{methods}

\alias{coef}
\alias{coef,fGARCH-method}
\alias{coef,fGARCHSPEC-method}


\title{GARCH coefficients methods}

\description{
    Coefficients methods \code{coef()} for GARCH Models.
}

\section{Methods}{
  Methods for \code{coef} defined in package \pkg{fGarch}:
  \describe{
    \item{object = "fGARCH"}{
      Extractor function for coefficients from a fitted
      GARCH model.
    }
    \item{object = "fGARCHSPEC"}{
      Extractor function for coefficients from a GARCH
      specification structure.
    }
  }
}


\note{
  
  \code{\link{coef}} is a generic function which extracts coefficients
  from objects returned by modeling functions.

}


\author{
  Diethelm Wuertz for the Rmetrics \R-port
}


\examples{
## garchSpec -
   # Use default parameters beside alpha:
   spec = garchSpec(model = list(alpha = c(0.05, 0.05)))
   spec
   coef(spec)

## garchSim -
   # Simulate an univariate "timeSeries" series from specification 'spec':
   x = garchSim(spec, n = 2000)
   x = x[,1]

## garchFit --
   fit = garchFit( ~ garch(1, 1), data = x, trace = FALSE)

## coef -
   coef(fit)
}

\keyword{models}
