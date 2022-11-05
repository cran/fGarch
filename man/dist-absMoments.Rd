\name{absMoments}
\alias{absMoments}


\title{Absolute moments of GARCH distributions}

\description{
    
    Computes absolute moments of the standard normal, 
    standardized GED, and standardized skew Student-t 
    distributions.
    
}

\usage{
absMoments(n, density = c("dnorm", "dged", "dstd"), \dots)
}

\arguments{

  \item{n}{
    the order of the absolute moment, can be a vector to compute
    several absolute moments at once.
  } 
  \item{density}{
    a character string naming a symmetric density function.
  }
  \item{\dots}{
    parameters passed to the density function.
  }
}


\details{

  \code{absMoments} returns a numeric vector of length \code{n} with the
  values of the absolute moments, as specified by \code{n}, of the
  selected probability density function (pdf).

  If \code{density} names one of the densities in the signature of
  \code{absMoments}, the moments are calculated from known
  formulas.

  Otherwise, numerical integration is used and an attribute is attached
  to the results to report an estimate of the error. Note that the
  density is assumed symmetric wihtout a check.
  
}

\value{
  a numeric vector
}

\references{
    
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages. 
    
}


\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\seealso{
  \code{\link{ged}},
  \code{\link{std}}
}
\examples{
## absMoment -

absMoments(1, "dstd", nu = 6)
absMoments(1, "dstd", nu = 600)
absMoments(1, "dstd", nu = 60000)
absMoments(1, "dstd", nu = 600000)

absMoments(1, "dnorm")

## excess kurtosis of t_nu is  6/(nu - 4)
nu <- 6
absMoments(2*2, "dstd", nu = nu) / absMoments(2*1, "dstd", nu = nu)^2 - 3
6/(nu-4)

## 4th moment for t_4 is infinite
absMoments(4, "dstd", nu = 4)

absMoments(1, "dged", nu = 4)
}


\keyword{distribution}

