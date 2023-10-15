\name{garchSim}
\alias{garchSim}

\concept{GARCH model}
\concept{APARCH model}

\concept{simulate GARCH}
\concept{simulate APARCH}


\title{Simulate univariate GARCH/APARCH time series}

\description{

  Simulates univariate GARCH/APARCH time series.

}


\usage{
garchSim(spec = garchSpec(), n = 100, n.start = 100, extended = FALSE) 
}

\arguments{

  \item{spec}{
      
    a specification object of class \code{"\linkS4class{fGARCHSPEC}"} as
    returned by \code{\link{garchSpec}}. See also below for further
    details.

  }   
  \item{n}{

    length of the output series, an integer value, by default
    \code{n=100}.

  }
  \item{n.start}{

    length of \sQuote{burn-in} period, by default 100.

  }
  \item{extended}{
    
    logical parameter specifying what to return. If \code{FALSE}, return
    the univariate GARCH/APARCH time series. If \code{TRUE}, return a
    multivariate time series containing also the volatility and
    conditional innovations time series.
      
  }
}

\details{

  \code{garchSim} simulates an univariate GARCH or APARCH time series
  process as specified by argument \code{spec}.  The default model
  specifies Bollerslev's GARCH(1,1) model with normally distributed
  innovations.
  
  \code{spec} is an object of class \code{"fGARCHSPEC"} as returned by
  the function \code{\link{garchSpec}}. It comes with a slot
  \code{@model} which is a list of just the numeric parameter
  entries. These are recognized and extracted for use by the function
  \code{garchSim}.

  % (GB) commenting out since the return value is 'timeSeries'.
  %
  % By default the series will be returned as an object of class
  % \code{"ts"} or as a \code{"numeric"} vector. Having time/date
  % positions, e.g. from an empirical process the numeric vector can be
  % easily transformed into other time series objects like
  % \code{"timeSeries"} or \code{"zoo"}. So

  One can estimate the parameters of a GARCH process from empirical data
  using the function \code{garchFit} and then simulate statistically
  equivalent GARCH processes with the same set of model parameters using
  the function \code{garchSim}.
  
  % (GB) commenting out since there is no argument 'returnClass'
  %  (apparently it is superseded by 'extend')
  %
  % The third entry in the argument \code{returnClass="mts"} allows to
  % return a trivariate time series, where the first column contains the
  % simulated \code{"garch"} process, the second column the conditional
  % standard deviations \code{"h"}, and the last column the innovations
  % named \code{"eps"}.
}

\value{

  the simulated time series as an objects of
  class \code{"timeSeries"} with attribute \code{"spec"} containing the
  specification of the model.

  If \code{extended} is \code{TRUE}, then the time series is
  multivariate and contains also the volatility, \code{sigma}, and the
  conditional innovations, \code{eps}.
  
}

\note{
  
  An undocumented feature (so, it should not be relied on) is that the
  returned time series is timed so that the last observation is the day
  before the date when the function is executed. This probably should be
  controlled by an additional argument in \code{garchSim}.
  
}

\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\seealso{
  \code{\link{garchSpec}},
  \code{\link{garchFit}}
}

\examples{
## garchSpec -
   spec = garchSpec()
   spec

## garchSim -
   # Simulate a "timeSeries" object:
   x = garchSim(spec, n = 50)
   class(x)
   print(x) 
   
## More simulations ...

   # Default GARCH(1,1) - uses default parameter settings
   spec = garchSpec(model = list())
   garchSim(spec, n = 10)
   
   # ARCH(2) - use default omega and specify alpha, set beta=0!
   spec = garchSpec(model = list(alpha = c(0.2, 0.4), beta = 0))
   garchSim(spec, n = 10)
   
   # AR(1)-ARCH(2) - use default mu, omega
   spec = garchSpec(model = list(ar = 0.5, alpha = c(0.3, 0.4), beta = 0))
   garchSim(spec, n = 10)
   
   # AR([1,5])-GARCH(1,1) - use default garch values and subset ar[.]
   spec = garchSpec(model = list(mu = 0.001, ar = c(0.5,0,0,0,0.1)))
   garchSim(spec, n = 10)
   
   # ARMA(1,2)-GARCH(1,1) - use default garch values
   spec = garchSpec(model = list(ar = 0.5, ma = c(0.3, -0.3)))  
   garchSim(spec, n = 10)
   
   # GARCH(1,1) - use default omega and specify alpha/beta
   spec = garchSpec(model = list(alpha = 0.2, beta = 0.7))
   garchSim(spec, n = 10)
   
   # GARCH(1,1) - specify omega/alpha/beta
   spec = garchSpec(model = list(omega = 1e-6, alpha = 0.1, beta = 0.8))
   garchSim(spec, n = 10)
   
   # GARCH(1,2) - use default omega and specify alpha[1]/beta[2]
   spec = garchSpec(model = list(alpha = 0.1, beta = c(0.4, 0.4)))
   garchSim(spec, n = 10)
   
   # GARCH(2,1) - use default omega and specify alpha[2]/beta[1]
   spec = garchSpec(model = list(alpha = c(0.12, 0.04), beta = 0.08))
   garchSim(spec, n = 10)
   
   # snorm-ARCH(1) - use defaults with skew Normal
   spec = garchSpec(model = list(beta = 0, skew = 0.8), cond.dist = "snorm")
   garchSim(spec, n = 10)
   
   # sged-GARCH(1,1) - using defaults with skew GED
   model = garchSpec(model = list(skew = 0.93, shape = 3), cond.dist = "sged")
   garchSim(model, n = 10)
   
   # Taylor Schwert GARCH(1,1) - this belongs to the family of APARCH Models
   spec = garchSpec(model = list(delta = 1))
   garchSim(spec, n = 10)
   
   # AR(1)-t-APARCH(2, 1) - a little bit more complex specification ...
   spec = garchSpec(model = list(mu = 1.0e-4, ar = 0.5, omega = 1.0e-6, 
       alpha = c(0.10, 0.05), gamma = c(0, 0), beta = 0.8, delta = 1.8, 
       shape = 4, skew = 0.85), cond.dist = "sstd")
   garchSim(spec, n = 10)

   garchSim(spec, n = 10, extended = TRUE)
   
}

\keyword{models}
\keyword{ts}
