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

  The method for \code{"fGARCH"} objects extracts the \code{@fitted}
  value slot from an object of class \code{"fGARCH"} as returned by the
  function \code{garchFit}.

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
## see examples for 'residuals()'
}

\keyword{models}
