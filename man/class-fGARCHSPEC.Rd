\name{fGARCHSPEC-class}
\alias{fGARCHSPEC-class}

\docType{class}

\alias{show,fGARCHSPEC-method}
\alias{update,fGARCHSPEC-method}

\concept{GARCH model}
\concept{APARCH model}
\concept{ARMA-GARCH model}
\concept{ARMA-APARCH model}


\title{Class "fGARCHSPEC"}

\description{

  Specification structure for an univariate GARCH time series model.

}


\section{Objects from the Class}{

  Objects can be created by calls of the function \code{garchSpec}.
  This object is specifies the parameters of an empirical GARCH process.
    
}


\section{Slots}{
  \describe{

    \item{\code{call}:}{Object of class \code{"call"}: 
      the call of the \code{garch} function.
    } 
    \item{\code{formula}:}{Object of class \code{"formula"}: 
      a list with two formula entries for the mean and variance 
      equation.
    }
    \item{\code{model}:}{Object of class \code{"list"}: 
      a list with the model parameters.
    }
    \item{\code{presample}:}{Object of class \code{"matrix"}: 
      a numeric matrix with presample values.
    }
    \item{\code{distribution}:}{Object of class \code{"character"}: 
      a character string with the name of the conditional distribution.
    }
    \item{\code{rseed}:}{Object of class \code{"numeric"}: 
      an integer with the random number generator seed.
    }
   
  }
}       


\section{Methods}{
  \describe{
    
    \item{show}{\code{signature(object = "fGARCHSPEC")}: 
      prints an object of class 'fGARCHSPEC'.
    }      
        
  }
}


\note{

  With Rmetrics Version 2.6.1 the class has been renamed from
  \code{"garchSpec"} to \code{"fGARCHSPEC"}.
    
}


\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\examples{
## garchSpec -
   spec = garchSpec()
   spec # print() or show() it
}

\keyword{programming}
