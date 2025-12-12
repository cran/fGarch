\name{plot-methods}
\alias{plot-methods}

\docType{methods}

\alias{plot}
\alias{plot,fGARCH,missing-method}

\concept{GARCH diagnostics}
\concept{GARCH goodness=of-fit}

\concept{VaR}
\concept{value-at-risk}
\concept{ES}
\concept{expected shortfall}


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

    a character string or a vector of positive integers specifying which
    plot(s) should be displayed, see section \sQuote{Details}.

  }
  \item{\dots}{
    optional arguments to be passed.
  }
 
}

\details{

  The \code{plot} method for \code{"fGARCH"} objects offers a selection
  of diagnostic, exploratory, and presentation plots from a menu.
  Argument \code{which} can be used to request specific plots.  This is
  particularly useful in scripts.

  If \code{which} is of length larger than one, all requested plots are
  produced. For this to be useful, the graphics window should be split
  beforehand in subwindows, e.g., using \code{par(mfrow = ...)},
  \code{par(mfcol = ...)}, or \code{layout()} (see section
  \sQuote{Examples}). If this is not done, then only the last plot will
  be visible.

  The following graphs are available:

  \tabular{rl}{
     1 \tab Time SeriesPlot\cr
     2 \tab Conditional Standard Deviation Plot\cr
     3 \tab Series Plot with 2 Conditional SD Superimposed\cr
     4 \tab Autocorrelation function Plot of Observations\cr
     5 \tab Autocorrelation function Plot of Squared Observations\cr
     6 \tab Cross Correlation Plot\cr
     7 \tab Residuals Plot\cr
     8 \tab Conditional Standard Deviations Plot\cr
     9 \tab Standardized Residuals Plot\cr
    10 \tab ACF Plot of Standardized Residuals\cr
    11 \tab ACF Plot of Squared Standardized Residuals\cr
    12 \tab Cross Correlation Plot between $r^2$ and r\cr
    13 \tab Quantile-Quantile Plot of Standardized Residuals\cr
    14 \tab Series with -VaR Superimposed\cr
    15 \tab Series with -ES Superimposed\cr
    16 \tab Series with -VaR & -ES Superimposed
  }        
}

\author{
  Diethelm Wuertz for the Rmetrics \R-port;

  VaR and ES graphs were added by Georgi N. Boshnakov in v4033.92
}

\seealso{
  \code{\link[=tsdiag.fGARCH]{fGARCH}} method for \code{tsdiag},

  \code{\link{garchFit}},
  class \code{\linkS4class{fGARCH}},

  \code{\link{predict}},
  \code{\link{fitted}},
  \code{\link{residuals}}
  \code{\link{VaR}}
  \code{\link{ES}}

  \code{\link{plot}}
}

\examples{
## simulate a Garch(1,1) time series
x <- garchSim(n = 200)
head(x) 

## fit GARCH(1,1) model
fit <- garchFit(formula = ~ garch(1, 1), data = x, trace = FALSE)

\dontrun{
## choose plots interactively
plot(fit)
}
   
## Batch Plot:
plot(fit, which = 3)

## a 2 by 2 matrix of plots
op <- par(mfrow = c(2,2))          # prepare 2x2 window
plot(fit, which = c(10, 11, 3, 16)) # plot
par(op)                            # restore the previous layout
}

\keyword{models}
