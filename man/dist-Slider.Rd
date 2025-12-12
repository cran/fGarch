\name{dist-Slider}

\alias{snormSlider}

\alias{stdSlider}
\alias{sstdSlider}

\alias{gedSlider}
\alias{sgedSlider}


\concept{standardized GED distribution}
\concept{GED distribution}

\concept{standardized skew GED distribution}
\concept{skew GED distribution}
\concept{skewed GED distribution}

\concept{skew distribution}
\concept{skewed distribution}

\concept{skewed normal distribution}

\concept{skew distribution}
\concept{skewed distribution}

\concept{t-distribution}
\concept{Student-t distribution}
\concept{standardized Student t distribution}
\concept{standardized Student-t distribution}

\concept{standardized skew Student-t distribution}
\concept{skew Student-t distribution}
\concept{skewed Student-t distribution}



\title{
  Visualise skew normal, (skew) Student-t and (skew) GED distributions
}

\description{
    
  Displays interactively the dependence of distributions on their
  parameters. Package 'fGarch' provides sliders for the Student-t, GED,
  skew normal, skew Student-t and skew GED distributions,

}

\usage{
stdSlider(type = c("dist", "rand"))
gedSlider(type = c("dist", "rand"))

snormSlider(type = c("dist", "rand"))
sstdSlider(type = c("dist", "rand"))
sgedSlider(type = c("dist", "rand"))
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
  \code{\link{snorm}},
  \code{\link{std}},
  \code{\link{sstd}},
  \code{\link{ged}},
  \code{\link{sged}}
}
  
\examples{
\dontrun{  
require(tcltk)

snormSlider("dist")
snormSlider("rand")

stdSlider("dist")
stdSlider("rand")

sstdSlider("dist")
sstdSlider("rand")

gedSlider("dist")
gedSlider("rand")

sgedSlider("dist")
sgedSlider("rand")
}
}

\keyword{distribution}
