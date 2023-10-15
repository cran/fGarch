\name{ged}
\alias{ged}

\alias{dged}
\alias{pged}
\alias{qged}
\alias{rged}


\concept{standardized GED distribution}
\concept{GED distribution}


\title{Standardized generalized error distribution}

\description{
    
  Functions to compute density, distribution function, quantile function
  and to generate random variates for the standardized generalized error
  distribution.
          
}


\usage{
dged(x, mean = 0, sd = 1, nu = 2, log = FALSE)
pged(q, mean = 0, sd = 1, nu = 2)
qged(p, mean = 0, sd = 1, nu = 2)
rged(n, mean = 0, sd = 1, nu = 2)
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
    shape parameter.
  }
  \item{log}{
    logical; if \code{TRUE}, densities are given as log densities.
  }
}

\details{

  The standardized GED is defined so that for a given \code{sd} it has
  the same variance, \code{sd^2}, for all values of the shape parameter,
  see the reference by Wuertz et al below.
  
  \code{dged} computes the density,
  \code{pged} the distribution function,
  \code{qged} the quantile function,
  and
  \code{rged} generates random deviates from the standardized-t
  distribution with the specified parameters.
  
}

\value{
  numeric vector
}

\references{

Nelson D.B. (1991);
    \emph{Conditional Heteroscedasticity in Asset Returns: A New Approach},
    Econometrica, 59, 347--370.
    
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages.

Wuertz D., Chalabi Y. and Luksan L. (????);
    \emph{Parameter estimation of ARMA  models with GARCH/APARCH errors: An R
      and SPlus software implementation},
    Preprint, 41 pages,
    \url{https://github.com/GeoBosh/fGarchDoc/blob/master/WurtzEtAlGarch.pdf}
}


\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\seealso{
  \code{\link{gedFit}},
  \code{\link{absMoments}},
  \code{\link{sged}} (skew GED),

  \code{\link{gedSlider}} for visualization
}

\examples{
## sged -
   par(mfrow = c(2, 2))
   set.seed(1953)
   r = rsged(n = 1000)
   plot(r, type = "l", main = "sged", col = "steelblue")
   
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
   box()
   x = seq(min(r), max(r), length = 201)
   lines(x, dsged(x), lwd = 2)
   
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
   lines(x, psged(x), lwd = 2)
   
   # Compute quantiles:
   round(qsged(psged(q = seq(-1, 5, by = 1))), digits = 6)
       
}

\keyword{distribution}
