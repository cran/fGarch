\name{formula-methods}
\alias{formula-methods}

\docType{methods}

\alias{formula}
\alias{formula,fGARCH-method}


\title{Extract GARCH model formula}

\description{

  Extracts formula from a formula GARCH object. 
    
}


\section{Methods}{
  Methods for \code{formula} defined in package \pkg{fGarch}:
  \describe{

    \item{object = "fGARCH"}{
      Extractor function for formula expression.
    }
    
  }
}


\details{

  \code{formula} is a generic function which extracts the formula
  expression from objects returned by modeling functions.
  
  The \code{"fGARCH"} method extracts the \code{@formula} expression
  slot from an object of class \code{"fGARCH"} as returned by the
  function \code{garchFit}.
    
  The returned formula has always a left hand side. If the argument
  \code{data} was an univariate time series and no name was specified to
  the series, then the left hand side is assigned the name of the
  data.set. In the multivariate case the rectangular \code{data} object
  must always have column names, otherwise the fitting will be stopped
  with an error message
    
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
  \code{\link{garchFit}},
  class \code{\linkS4class{fGARCH}}
}

\examples{
## garchFit -
   fit = garchFit(~garch(1, 1), data = garchSim(), trace = FALSE)
   
## formula - 
   formula(fit)

## A Bivariate series and mis-specified formula:
   x = garchSim(n = 500)
   y = garchSim(n = 500)
   z = cbind(x, y)
   colnames(z)
   class(z)
   \dontrun{
   garchFit(z ~garch(1, 1), data = z, trace = FALSE)
   }
   # Returns:
   # Error in .garchArgsParser(formula = formula, data = data, trace = FALSE) :  
   #   Formula and data units do not match.
   
## Doubled column names in data set - formula can't fit:
   colnames(z) <- c("x", "x")
   z[1:6,]
   \dontrun{
   garchFit(x ~garch(1, 1), data = z, trace = FALSE)
   }
   # Again the error will be noticed:
   # Error in garchFit(x ~ garch(1, 1), data = z) : 
   #   Column names of data are not unique.

## Missing column names in data set - formula can't fit:
   z.mat <- as.matrix(z)
   colnames(z.mat) <- NULL
   z.mat[1:6,]
   \dontrun{
   garchFit(x ~ garch(1, 1), data = z.mat, trace = FALSE)
   }
   # Again the error will be noticed:
   # Error in .garchArgsParser(formula = formula, data = data, trace = FALSE) : 
   #   Formula and data units do not match
}

\keyword{models}
