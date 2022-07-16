\name{show-methods}

\title{GARCH Modelling Show Methods}

\docType{methods}

\alias{show-methods}
\alias{show,ANY-method}
\alias{show,fGARCH-method}
\alias{show,fGARCHSPEC-method}

\description{
    Show methods for GARCH Modelling, i.e., methods for the
    \code{\link{show}()} generic function (which typically and also here
    are also used as \emph{the} \code{\link{print}} methods for S4 objects.
}


\section{Methods}{
  \describe{
    \item{object = "ANY"}{
        Generic function.
        }
    \item{object = "fGARCH"}{
        Print function for objects of class \code{"fGARCH"}.
        }
    \item{object = "fGARCHSPEC"}{
        Print function for objects of class \code{"fGARCH"}.
        }

    }
}


\author{
    Diethelm Wuertz for the Rmetrics \R-port.
}


\examples{
## garchSpec -
   spec = garchSpec()
   spec # print() or show() it

## garchSim -
   x = garchSim(spec, n = 500)

## garchFit -
   fit = garchFit(~ garch(1, 1), data = x)
   fit # == print(fit) and also == show(fit)
}
