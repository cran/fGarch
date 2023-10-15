\name{fGARCH-class}
\alias{fGARCH-class}

\docType{class}

\alias{show,fGARCH-method}
\alias{update,fGARCH-method}

\concept{GARCH model}
\concept{APARCH model}
\concept{ARMA-GARCH model}
\concept{ARMA-APARCH model}
\concept{AR-GARCH model}
\concept{AR-APARCH model}
\concept{MA-GARCH model}
\concept{MA-APARCH model}



\title{Class "fGARCH"}

\description{ 

  The class 'fGARCH' represents a model of an heteroskedastic time
  series process.

}

\section{Objects from the Class}{

  Objects can be created by calls of the function \code{garchFit}.
  This object is a parameter estimate of an empirical GARCH process.
    
}

\section{Slots}{
  \describe{

    \item{\code{call}:}{Object of class \code{"call"}: 
      the call of the \code{garch} function.
    } 
    \item{\code{formula}:}{Object of class \code{"formula"}: 
      a formula object specifying the mean and variance equations.
    }
    \item{\code{method}:}{Object of class \code{"character"}: 
      a string denoting the optimization method, by default
      \code{"Max Log-Likelihood Estimation"}. 
    }
    \item{\code{data}:}{Object of class \code{"list"}: 
      a list with one entry named \code{x}, containing the data of
      the time series to be estimated, the same as given by the
      input argument \code{series}.
    }
    \item{\code{fit}:}{Object of class \code{"list"}: 
      a list with the results from the parameter estimation. The entries
      of the list depend on the selected algorithm, see below.
    }
    \item{\code{residuals}:}{Object of class \code{"numeric"}: 
      a numeric vector with the (raw, unstandardized) residual values.
    }
    \item{\code{fitted}:}{Object of class \code{"numeric"}: 
      a numeric vector with the fitted values.
    }
    \item{\code{h.t}:}{Object of class \code{"numeric"}: 

      a numeric vector with the conditional variances (\eqn{h_t =
      \sigma_t^\delta}{h.t = sigma.t^delta}).

    }
    \item{\code{sigma.t}:}{Object of class \code{"numeric"}: 
      a numeric vector with the conditional standard deviations.
    }
    \item{\code{title}:}{Object of class \code{"character"}: 
      a title string.
    }
    \item{\code{description}:}{Object of class \code{"character"}: 
      a string with a brief description.
    }
    
  }
}       


\section{Methods}{
  \describe{
    
    \item{plot}{\code{signature(x = "fGARCH", y = "missing")}: 
      plots an object of class \code{"fGARCH"}.
    }
    \item{show}{\code{signature(object = "fGARCH")}: 
      prints an object of class \code{"fGARCH"}.
    }      
    \item{summary}{\code{signature(object = "fGARCH")}: 
      summarizes an object of class \code{"fGARCH"}.
    }
    \item{predict}{\code{signature(object = "fGARCH")}: 
      forecasts mean and volatility from an object of class \code{"fGARCH"}.
    }
    \item{fitted}{\code{signature(object = "fGARCH")}: 
      extracts fitted values from an object of class \code{"fGARCH"}.
    }     
    \item{residuals}{\code{signature(object = "fGARCH")}: 
      extracts fresiduals from an object of class \code{"fGARCH"}.
    }
    \item{volatility}{\code{signature(object = "fGARCH")}: 
      extracts conditional volatility from an object of class \code{"fGARCH"}.
    }
    \item{coef}{\code{signature(object = "fGARCH")}: 
      extracts fitted coefficients from an object of class \code{"fGARCH"}.
    }
    \item{formula}{\code{signature(x = "fGARCH")}: 
      extracts formula expression from an object of class \code{"fGARCH"}.
    } 
    
  }
}

\author{
  Diethelm Wuertz and Rmetrics Core Team
}

\seealso{
  \code{\link{garchFit}},
  \code{\link{garchSpec}},
  \code{\link{garchFitControl}}
}

\examples{
## simulate a time series, fit a GARCH(1,1) model, and show it:
x <- garchSim( garchSpec(), n = 500)
fit <- garchFit(~ garch(1, 1), data = x, trace = FALSE)
fit # == print(fit) and also == show(fit)
}

\keyword{programming}
\keyword{ts}
