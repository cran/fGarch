\name{gedSlider}
\alias{gedSlider}


\title{Generalized error distribution slider}

\description{
    
  Displays interactively the dependence of the GED distribution on its
  parameters.
          
}


\usage{
gedSlider(type = c("dist", "rand"))
}

\arguments{

  \item{type}{

    a character string denoting which interactive plot should be
    displayed. Either a distribution plot \code{type = "dist"}, the
    default value, or a random variates plot, \code{type = "rand"}.

  }
}


\value{
  a Tcl object
}


\references{

Nelson D.B. (1991);
    \emph{Conditional Heteroscedasticity in Asset Returns: A New Approach},
    Econometrica, 59, 347--370.
    
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages. 
    
}


\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\seealso{
  \code{\link{ged}},
  \code{\link{gedFit}}
}
  
\examples{
\dontrun{  
## gedSlider -
   require(tcltk)
   gedSlider("dist")
   gedSlider("rand")
}
}

\keyword{distribution}

