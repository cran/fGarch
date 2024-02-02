\name{VaR}
\alias{VaR}

\alias{VaR.fGARCH}

\alias{ES}
\alias{ES.fGARCH}

\concept{VaR}
\concept{value-at-risk}
\concept{ES}
\concept{expected shortfall}


\title{Compute Value-at-Risk (VaR) and expected shortfall (ES)}

\description{

  Compute Value-at-Risk (VaR) and Expected Shortfall (ES) for a fitted
  GARCH-APARCH model.

}

\usage{
\method{VaR}{fGARCH}(dist, p_loss = 0.05, ..., tol)

\method{ES}{fGARCH}(dist, p_loss = 0.05, ...)
}

\arguments{
  \item{dist}{
    an object from class \code{"fGARCH"}, obtained from \code{garchFit()}.
  }

  \item{p_loss}{level, default is 0.05.}

  \item{...}{not used.}

  \item{tol}{tollerance}
}

\details{

  We provide methods for the generic functions \code{cvar::VaR} and
  \code{cvar::ES}.

}

\note{

  We use the traditional definition of VaR as the negated lower
  quantile. For example, if \eqn{X} are returns on an asset,
  VAR\eqn{{}_\alpha}{_a} = \eqn{-q_\alpha}{-q_a}, where
  \eqn{q_\alpha}{-q_a} is the lower \eqn{\alpha}{a} quantile of \eqn{X}.
  Equivalently, VAR\eqn{{}_\alpha}{_a} is equal to the lower
  \eqn{1-\alpha}{1-a} quantile of \eqn{-X}.

}

\seealso{

  \code{\link[cvar]{VaR}} and \code{\link[cvar]{ES}} in package \pkg{cvar}

}
\examples{
x <- garchSim( garchSpec(), n = 500)
fit <- garchFit(~ garch(1, 1), data = x, trace = FALSE)
head(VaR(fit))
head(ES(fit))
}
