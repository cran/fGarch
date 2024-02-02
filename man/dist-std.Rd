\name{std}
\alias{std}

\alias{dstd}
\alias{pstd}
\alias{qstd}
\alias{rstd}

\concept{t-distribution}
\concept{Student-t distribution}
\concept{standardized Student t distribution}
\concept{standardized Student t-distribution}


\title{Standardized Student-t distribution}

\description{
    
  Functions to compute density, distribution function, quantile function
  and to generate random variates for the standardized Student-t
  distribution.
       
}

\usage{
dstd(x, mean = 0, sd = 1, nu = 5, log = FALSE)
pstd(q, mean = 0, sd = 1, nu = 5)
qstd(p, mean = 0, sd = 1, nu = 5)
rstd(n, mean = 0, sd = 1, nu = 5)
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
  \item{log}{
    logical; if \code{TRUE}, densities are given as log densities.
  }
}

\details{

  The standardized Student-t distribution is defined so that for a given
  \code{sd} it has the same variance, \code{sd^2}, for all degrees of
  freedom. For comparison, the variance of the usual Student-t
  distribution is \code{nu/(nu-2)}, where \code{nu} is the degrees of
  freedom.  The usual Student-t distribution is obtained by setting
  \code{sd = sqrt(nu/(nu - 2))}.

  Argument \code{nu} must be greater than 2. Although there is a default
  value for \code{nu}, it is rather arbitrary and relying on it is
  strongly discouraged.
  
  \code{dstd} computes the density,
  \code{pstd} the distribution function,
  \code{qstd} the quantile function,
  and
  \code{rstd} generates random deviates from the standardized-t
  distribution with the specified parameters.
  
}

\value{
  numeric vector
}


\references{
   
Fernandez C., Steel M.F.J. (2000); 
    \emph{On Bayesian Modelling of Fat Tails and Skewness},
    Preprint, 31 pages. 
    
Wuertz D., Chalabi Y. and Luksan L. (2006);
    \emph{Parameter estimation of ARMA  models with GARCH/APARCH errors: An R
      and SPlus software implementation},
    Preprint, 41 pages,
    \url{https://github.com/GeoBosh/fGarchDoc/blob/master/WurtzEtAlGarch.pdf}
}

\author{
  Diethelm Wuertz for the Rmetrics \R-port
}

\seealso{
  \code{\link{stdFit}} (fit).
  \code{\link{stdSlider}} (visualize),

  \code{\link{absMoments}}
}
\examples{
## std -

pstd(1, sd = sqrt(5/(5-2)), nu = 5) == pt(1, df = 5) # TRUE

   par(mfrow = c(2, 2))
   set.seed(1953)
   r = rstd(n = 1000)
   plot(r, type = "l", main = "sstd", col = "steelblue")
   
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
   box()
   x = seq(min(r), max(r), length = 201)
   lines(x, dstd(x), lwd = 2)
   
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
   lines(x, pstd(x), lwd = 2)
   
   # Compute quantiles:
   round(qstd(pstd(q = seq(-1, 5, by = 1))), digits = 6)
}

\keyword{distribution}
