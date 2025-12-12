\name{stats-tsdiag}
\alias{stats-tsdiag}

\alias{tsdiag}
\alias{tsdiag.fGARCH}


\title{Diagnostic plots and statistics for fitted GARCH models}

\description{

  Produce diagnostics for fitted GARCH/APARCH models. The method offers
  several tests, plots of autocorrelations and partial autocorrelations
  of the standardised conditional residuals, ability to control which
  graphs are produced (including interactively), as well as their
  layout.

}

\usage{
\method{tsdiag}{fGARCH}(object, gof.lag = NULL, ask = FALSE, \dots, plot = c(4L, 5L, 7L), 
       layout = NULL)
}

\arguments{
  \item{object}{
    an object from class \code{"fGARCH"}, as returned by
    \code{\link{garchFit}}. 
  }
  \item{gof.lag}{maximal lag for portmanteau tests.}
  \item{ask}{if \code{TRUE} present a menu of available plots, see Details.}
  \item{\dots}{not used.}
  \item{plot}{
    if \code{TRUE} all available plots; a vector of positive integers
    specifies a subset of the available plots.
  }
  \item{layout}{
    a list with arguments for \code{graphics::layout} for the plots. The
    default plots the autocorrelations of the standardised conditional
    residuals and their squares, as well as a QQ-plot for the fitted
    conditional distribution.
  }
}
\details{

  Compute and graph diagnostics for fitted ARMA-GARCH/APARCH models.

  % The method offers several portmanteau tests (including Ljung-Box,
  % Li-McLeod and Box-Pierce), plots of autocorrelations and partial
  % autocorrelations of the residuals, ability to control which graphs are
  % produced (including interactively), as well as their layout.
  % 
  % The method always makes a correction of the degrees of freedom of the
  % portmanteau tests (roughly, subtracting the number of estimated ARMA
  % parameters). Note that \code{stats::tsdiag} doesn't do that.

  \code{plot} can be \code{TRUE} to ask for all plots or a vector of
  positive integers specifying which plots to consider. Currently the
  following options are available:
  \Sexpr[stage=build,results=rd]{paste0("\\\\tabular{rl}{",
  paste0(1:length(fGarch:::.tsdiag_choices), " \\\\tab ",
  fGarch:::.tsdiag_choices, collapse = " \\\\cr\n "), "}\n")} The
  default produces plots of autocorrelations and partial
  autocorrelations of the standardised conditional residuals, as well as
  a QQ-plot for the fitted conditional distribution.  If \code{plot} is
  \code{TRUE}, you probably need also \code{ask = TRUE}.

  If argument \code{plot} is of length two the graphics window is split
  into 2 equal subwindows. Argument \code{layout} can still be used to
  change this.  If argument \code{plot} is of length one the graphics
  window is not split at all.
  
  In interactive sessions, if the number of requested graphs (as
  specified by argument \code{plot}) is larger than the number of graphs
  specified by the layout (by default 3), the function makes the first
  graph and then presents a menu of the requested plots. 

  Argument \code{layout} can be used to change the layout of the plot,
  for example to put two graphs per plot, see the examples. Currently it
  should be a list of arguments for \code{\link[graphics]{layout}}, see \code{?layout}.
  Don't call \code{layout} youself, as that will change the graphics
  device prematurely.
  
  The computed results are returned (invisibly). This is another
  difference from \code{stats::tsdiag} which doesn't return them.
}
\value{
  (experimental, may change)
  a list with components:
  \item{residuals}{standardised conditional residuals,}
  \item{gof}{

    goodness-of-fit tests, pretending parameters are known,

  }
  \item{gof_composite}{

    goodness-of-fit tests taking into account that the parameters are
    estimated.

  }

  Only components that are actually computed are included, the rest are
  NULL or absent.
}

\author{Georgi N. boshnakov}

\seealso{
  \code{\link[=plot,fGARCH,missing-method]{fGARCH}} method for \code{plot},

  \code{\link[stats]{tsdiag}}
}

\examples{
set.seed(20230612)
x <- garchSim(n = 200)
fit <- garchFit(formula = ~ garch(1, 1), data = x, trace = FALSE)
fit_test <- tsdiag(fit)
fit_test

## 2x2 matrix with acf of r, r^2 on diag, cor(r,r^2) below it, and qq-plot
tsdiag(fit, plot = c(4, 6, 7, 5), layout = list(matrix(1:4, nrow = 2)))
}

\keyword{ts}
\keyword{htest}
\concept{diagnostics}
\concept{diagnostic plots}
