\name{sstd}
\alias{sstd}

\alias{dsstd}
\alias{psstd}
\alias{qsstd}
\alias{rsstd}


\title{Skew Student-t distribution}

\description{
    
  Functions to compute density, distribution function, quantile function
  and to generate random variates for the skew Student-t distribution.
       
}


\usage{
dsstd(x, mean = 0, sd = 1, nu = 5, xi = 1.5, log = FALSE)
psstd(q, mean = 0, sd = 1, nu = 5, xi = 1.5)
qsstd(p, mean = 0, sd = 1, nu = 5, xi = 1.5)
rsstd(n, mean = 0, sd = 1, nu = 5, xi = 1.5)
}


\arguments{
  \item{x, q}{
    a numeric vector of quantiles.
  }
  \item{p}{
    a numeric vector of probabilities.
  }
  \item{n}{
    number of observations to simulate.
  } 
  \item{mean}{
    location parameter.
  }
  \item{sd}{
    scale parameter.
  }
  \item{nu}{
    shape parameter (degrees of freedom).
  }
  \item{xi}{
    skewness parameter.
  }
  \item{log}{
    logical; if \code{TRUE}, densities are given as log densities.
  }
    
}


\details{

  \code{dsstd} computes the density,
  \code{psstd} the distribution function,
  \code{qsstd} the quantile function, and
  \code{rsstd} generates random deviates.

}

\value{
  numeric vector
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
  \code{\link{sstdFit}} (fit),
  \code{\link{sstdSlider}} (visualize)
}

\examples{
## sstd -
   par(mfrow = c(2, 2))
   set.seed(1953)
   r = rsstd(n = 1000)
   plot(r, type = "l", main = "sstd", col = "steelblue")
   
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
   box()
   x = seq(min(r), max(r), length = 201)
   lines(x, dsstd(x), lwd = 2)
   
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
   lines(x, psstd(x), lwd = 2)
   
   # Compute quantiles:
   round(qsstd(psstd(q = seq(-1, 5, by = 1))), digits = 6)
}

\keyword{distribution}
