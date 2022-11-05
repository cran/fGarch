\name{snorm}
\alias{snorm}

\alias{dsnorm}
\alias{psnorm}
\alias{qsnorm}
\alias{rsnorm}

\concept{skewed normal distribution}

\concept{skew distribution}
\concept{skewed distribution}


\title{Skew normal distribution}

\description{

  Functions to compute density, distribution function, quantile function
  and to generate random variates for the skew normal distribution.

}

\usage{
dsnorm(x, mean = 0, sd = 1, xi = 1.5, log = FALSE)
psnorm(q, mean = 0, sd = 1, xi = 1.5)
qsnorm(p, mean = 0, sd = 1, xi = 1.5)
rsnorm(n, mean = 0, sd = 1, xi = 1.5)
}

\arguments{
  \item{x, q}{a numeric vector of quantiles.}
  \item{p}{a numeric vector of probabilities.}
  \item{n}{the number of observations.} 
  \item{mean}{location parameter.}
  \item{sd}{scale parameter.}
  \item{xi}{skewness parameter.}
  \item{log}{a logical; if TRUE, densities are given as log densities.}
}


\details{

  \code{dsnorm} computed the density,
  \code{psnorm} the distribution function,
  \code{qsnorm} the quantile function,
  and
  \code{rsnorm} generates random deviates.

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
  \code{\link{snormFit}} (fit),
  \code{\link{snormSlider}} (visualize),

  \code{\link{sstd}} (skew Student-t),
  \code{\link{sged}} (skew GED)

}

\examples{
## snorm -
   # Ranbdom Numbers:
   par(mfrow = c(2, 2))
   set.seed(1953)
   r = rsnorm(n = 1000)
   plot(r, type = "l", main = "snorm", col = "steelblue")
   
   # Plot empirical density and compare with true density:
   hist(r, n = 25, probability = TRUE, border = "white", col = "steelblue")
   box()
   x = seq(min(r), max(r), length = 201)
   lines(x, dsnorm(x), lwd = 2)
   
   # Plot df and compare with true df:
   plot(sort(r), (1:1000/1000), main = "Probability", col = "steelblue",
     ylab = "Probability")
   lines(x, psnorm(x), lwd = 2)
   
   # Compute quantiles:
   round(qsnorm(psnorm(q = seq(-1, 5, by = 1))), digits = 6)
}

\keyword{distribution}
