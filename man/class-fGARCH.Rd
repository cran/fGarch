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



\title{Class "fGARCH" - fitted ARMA-GARCH/APARCH models}

\description{ 

  Class 'fGARCH' represents models fitted to heteroskedastic time
  series, including ARCH, GARCH, APARCH, ARMA-GARCH and ARMA-APARCH
  models.

}

\section{Objects from the Class}{

  Objects from class \code{"fGARCH"} can be created by calls of the
  function \code{\link{garchFit}}.
    
}

\section{Slots}{
  \describe{

    \item{\code{call}:}{Object of class \code{"call"},
      the call used to fit the model and create the object.
    } 
    \item{\code{formula}:}{Object of class \code{"formula"},
      a formula object representing the mean and variance equations.
    }
    \item{\code{method}:}{Object of class \code{"character"},
      a string denoting the optimization method, by default
      \code{"Max Log-Likelihood Estimation"}.
    }
    \item{\code{data}:}{Object of class \code{"list"},
      a list with one entry, \code{x}, containing the data of
      the time series to which the model is fitted.
    }
    \item{\code{fit}:}{Object of class \code{"list"},
      a list with the results from the parameter estimation. The entries
      of the list depend on the selected algorithm, see below.
    }
    \item{\code{residuals}:}{Object of class \code{"numeric"},
      the raw, unstandardized residuals.
    }
    \item{\code{fitted}:}{Object of class \code{"numeric"},
      the fitted values.
    }
    \item{\code{h.t}:}{Object of class \code{"numeric"},
      the conditional variances (\eqn{h_t = \sigma_t^\delta}{h.t =
      sigma.t^delta}).
    }
    \item{\code{sigma.t}:}{Object of class \code{"numeric"},
      the conditional standard deviations.
    }
    \item{\code{title}:}{Object of class \code{"character"},
      a title string.
    }
    \item{\code{description}:}{Object of class \code{"character"},
      a string with a brief description.
    }
  }
}       


% \item{volatility}{\code{signature(object = "fGARCH")}: 
%   extracts conditional volatility from an object of class \code{"fGARCH"}.
% }
\section{Methods}{
  
  Besides the S4 methods described below, the are \code{"fGARCH"}
  methods (S3) for \code{tsdiag} (\code{\link{tsdiag}}), VaR
  (\code{\link{VaR}}), expected shortfall (\code{\link{ES}}), volatility
  (\code{\link{volatility}}), and maybe others.
  
  \describe{
    
    \item{plot}{\code{signature(x = "fGARCH", y = "missing")}: 
      plots an object of class \code{"fGARCH"}, see the
      \code{\link[=plot,fGARCH,missing-method]{help page}} of the method
      for details and options.
    }
    \item{show}{\code{signature(object = "fGARCH")}: 
      prints the object.
    }      
    \item{summary}{\code{signature(object = "fGARCH")}:
      summarizes the object. The \code{\link[=summary-methods]{help
      page}} of the \code{"fGARCH"} method gives details on the output,
      as well as interpretation of the results.
    }
    \item{predict}{\code{signature(object = "fGARCH")}:
      Computes forecasts of the mean and some measures of risk (such as
      volatility, value-at-risk and expected shortfall), see the
      method's \code{\link[=predict-methods]{help page}} for full
      details.
    }
    \item{fitted}{\code{signature(object = "fGARCH")}:
      extracts fitted values from the object
      (\code{\link[=fitted-methods]{help page}}).
    }     
    \item{residuals}{\code{signature(object = "fGARCH")}:
      returns residuals from the fitted model
      (\code{\link[=residuals-methods]{help page}}).
    }
    \item{coef}{\code{signature(object = "fGARCH")}:
      extracts the estimated coefficients
      (\code{\link[=coef-methods]{help page}}).
    }
    \item{formula}{\code{signature(x = "fGARCH")}:
      extracts the formula expression, see the method's
      \code{\link[=formula-methods]{help page}}.
    } 
    \item{update}{\code{signature(object = "fGARCH")}:
      ...
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
coef(fit)
summary(fit)
fit # == print(fit) and also == show(fit)
}

\keyword{programming}
\keyword{ts}
