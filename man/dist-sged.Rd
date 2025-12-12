\name{sged}
\alias{sged}

\alias{dsged}
\alias{psged}
\alias{qsged}
\alias{rsged}

\concept{standardized skew GED distribution}
\concept{skew GED distribution}
\concept{skewed GED distribution}

\concept{skew distribution}
\concept{skewed distribution}


\title{Skew generalized error distribution}

\description{
    
  Functions to compute density, distribution function, quantile function
  and to generate random variates for the skew generalized error
  distribution.
  Note that there are different ways to define a skew GED
  distribution, check section Details to see if this is the distribution
  you want.
          
}

\usage{
dsged(x, mean = 0, sd = 1, nu = 2, xi = 1.5, log = FALSE)
psged(q, mean = 0, sd = 1, nu = 2, xi = 1.5)
qsged(p, mean = 0, sd = 1, nu = 2, xi = 1.5)
rsged(n, mean = 0, sd = 1, nu = 2, xi = 1.5)
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
  \item{xi}{
    skewness parameter, a positive number. \code{xi = 1} gives a
    symmetric distribution (here normal).
  }
  \item{log}{
    logical; if \code{TRUE}, densities are given as log densities.
  }
}

\details{

  Skewness is based on the scheme by Fernandez and Steel (2000) but with
  different parameterisation as discussed by Wuertz et al (????), in
  that the skewness (\code{xi}), degrees of freedom (\code{nu}), 
  mean (\code{mean}) and standard deviation (\code{sd}) are independent
  parameters.
    
  \code{dsged} computes the density,
  \code{psged} the distribution function,
  \code{qsged} the quantile function, and
  \code{rsged} generates random deviates.

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
  \code{\link{sgedFit}} (fit),
  \code{\link{sgedSlider}} (visualize),

  \code{\link{snorm}} (skew normal),
  \code{\link{sstd}} (skew t)

  \code{\link{ged}} (symmetric GED)
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
