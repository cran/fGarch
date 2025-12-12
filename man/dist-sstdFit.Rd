\name{sstdFit}
\alias{sstdFit}


\title{Skew Student-t distribution parameter estimation}

\description{
    
  Fits the parameters of the skew Student-t distribution.
       
}

\usage{
sstdFit(x, \dots)
}

\arguments{
  \item{x}{
    a numeric vector of quantiles.
  }
  \item{\dots}{
    parameters passed to the optimization function \code{nlm}. 
  }
}

\value{
  a list with the following components: 
  \item{par}{
    the best set of parameters found; 
  }
  \item{objective}{
    the value of objective corresponding to \code{par};
  }
  \item{convergence}{
    an integer code. 0 indicates successful convergence;
  }  
  \item{message}{
    a character string giving any additional information returned by the
    optimizer, or NULL. For details, see PORT documentation;
  }  
  \item{iterations}{
    number of iterations performed;
  } 
  \item{evaluations}{
    number of objective function and gradient function evaluations.
  }
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
  \code{\link{sstd}},
  \code{\link{stdFit}}
}

\examples{
set.seed(1953)
r <- rsstd(n = 1000)

sstdFit(r)
}

\keyword{distribution}
