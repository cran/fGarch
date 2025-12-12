\name{fUGARCHSPEC-class}
\alias{fUGARCHSPEC-class}
\alias{.ugarchFit}
\alias{.ugarchSpec}

\docType{class}

\title{Class 'fUGARCHSPEC'}
\description{Class 'fUGARCHSPEC'.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("fUGARCHSPEC", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"list"} ~~ }
    \item{\code{distribution}:}{Object of class \code{"list"} ~~ }
    \item{\code{optimization}:}{Object of class \code{"list"} ~~ }
    \item{\code{documentation}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
No methods defined with class \code{"fUGARCHSPEC"} in the signature.
}

\note{

  Was this class meant to replace \code{"fGARCHSPEC"}?
  
  (GNB) This class seems to be meant for internal use by the package.

  (GNB) Amendment: no, functions \code{.ugarchSpec} and
        \code{.ugarchFit} are exported. \code{.ugarchFit} fits the model
        from a spec, unlike \code{garchFit}.

  TODO: There is something unfinished here. Check and sort out. See also
        \code{\link{fUGARCHSPEC-class}}
  
}

\seealso{
  class \code{"\linkS4class{fGARCH}"}
}

\examples{
showClass("fUGARCHSPEC")
}

\keyword{classes}
