\name{plot-methods}
\alias{plot-methods}

\docType{methods}

\alias{plot}
\alias{plot,fGARCH,missing-method}

\concept{GARCH diagnostics}
\concept{GARCH goodness=of-fit}


\title{GARCH plot methods}

\description{

  Plot methods for GARCH modelling. 

}


\usage{
\S4method{plot}{fGARCH,missing}(x, which = "ask", \dots)
}

\arguments{

  \item{x}{
    an object of class \code{"fGARCH"}.
  } 
  \item{which}{
    a character string or a vector of positive denoting which plot
    should be displayed, see section \sQuote{Details}.
  }
  \item{\dots}{
    optional arguments to be passed.
  }
 
}

\details{

  The \code{plot} method for \code{"fGARCH"} objects offers a selection
  of diagnostic, exploratory, and presentation plots from a
  menu. Argument \code{which} can be used to request specific plots,
  which is particularly useful in scripts.

  The generic function \code{plot} allows to display 13 graphs:

  \tabular{l}{
    Time SeriesPlot\cr
    Conditional Standard Deviation Plot\cr
    Series Plot with 2 Conditional SD Superimposed\cr
    Autocorrelation function Plot of Observations\cr
    Autocorrelation function Plot of Squared Observations\cr
    Cross Correlation Plot\cr
    Residuals Plot\cr
    Conditional Standard Deviations Plot\cr
    Standardized Residuals Plot\cr
    ACF Plot of Standardized Residuals\cr
    ACF Plot of Squared Standardized Residuals\cr
    Cross Correlation Plot between $r^2$ and r\cr
    Quantile-Quantile Plot of Standardized Residuals
  }        
}

\author{
  Diethelm Wuertz for the Rmetrics \R-port
}


\seealso{
  \code{\link[=tsdiag.fGARCH]{fGARCH}} method for \code{tsdiag},

  \code{\link{garchFit}},
  class \code{\linkS4class{fGARCH}},

  \code{\link{predict}},
  \code{\link{fitted}},
  \code{\link{residuals}}

  \code{\link{plot}}
}

\examples{
## garchSim -
   # Default Garch(1,1) Model:
   x = garchSim(n = 200)
   head(x) 

## garchFit - 
   fit = garchFit(formula = ~ garch(1, 1), data = x, trace = FALSE)
   
## Batch Plot:
   plot(fit, which = 3)
   
\dontrun{
## Plot:
   # Interactive Plot:
   plot(fit) 
}
}

\keyword{models}
