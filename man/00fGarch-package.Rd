\name{fGarch-package}
\alias{fGarch-package}

\docType{package}

\alias{fGarch}

\concept{GARCH model}
\concept{APARCH model}

\concept{skew distribution}
\concept{fit skew distribution}

\concept{volatility}



\title{Modelling heterskedasticity in financial time series}

\description{

  The Rmetrics \pkg{fGarch} package is a collection of functions to
  analyze and model heteroskedastic behavior in financial time series.

}

\author{
Diethelm Wuertz [aut] (original code),
Yohan Chalabi [aut],
Tobias Setz [aut],
Martin Maechler [ctb] (<https://orcid.org/0000-0002-8685-9910>),
Chris Boudt [ctb]
Pierre Chausse [ctb],
Michal Miklovac [ctb],
Georgi N. Boshnakov [cre, ctb]

Maintainer: Georgi N. Boshnakov <georgi.boshnakov@manchester.ac.uk>
}

% \details{
%   
% }

\section{1 Introduction}{

  GARCH, Generalized Autoregressive Conditional Heteroskedastic, models 
  have become important in the analysis of time series data, particularly 
  in financial applications when the goal is to analyze and forecast 
  volatility. 
  
  For this purpose, the family of GARCH functions offers functions for
  simulating, estimating and forecasting various univariate GARCH-type
  time series models in the conditional variance and an ARMA
  specification in the conditional mean. The function
  \code{\link{garchFit}} is a numerical implementation of the maximum
  log-likelihood approach under different assumptions, Normal,
  Student-t, GED errors or their skewed versions.  The parameter
  estimates are checked by several diagnostic analysis tools including
  graphical features and hypothesis tests. Functions to compute n-step
  ahead forecasts of both the conditional mean and variance are also
  available.

  The number of GARCH models is immense, but the most influential models 
  were the first. Beside the standard ARCH model introduced by Engle [1982] 
  and the GARCH model introduced by Bollerslev [1986], the function
  \code{garchFit} also includes the more general class of asymmetric power
  ARCH models, named APARCH, introduced by Ding, Granger and Engle [1993]. 
  The APARCH models include as special cases the TS-GARCH model of 
  Taylor [1986] and Schwert [1989], the GJR-GARCH model of Glosten, 
  Jaganathan, and Runkle [1993], the T-ARCH model of Zakoian [1993], the 
  N-ARCH model of Higgins and Bera [1992], and the Log-ARCH model of
  Geweke [1986] and Pentula [1986].
  
  There exist a collection of review articles by Bollerslev, Chou and 
  Kroner [1992], Bera and Higgins [1993], Bollerslev, Engle and
  Nelson [1994], Engle [2001], Engle and Patton [2001], and Li, Ling 
  and McAleer [2002] which give a good overview of the scope of the 
  research.

}


\section{2 Time series simulation}{

  Functions to simulate artificial GARCH and APARCH time series
  processes.
    
  \tabular{ll}{ 
    \code{\link{garchSpec}} \tab specifies an univariate GARCH time series model\cr
    \code{\link{garchSim}}  \tab simulates a GARCH/APARCH process
  }
}


\section{3 Parameter estimation}{

  Functions to fit the parameters of GARCH and APARCH time series
  processes.
    
  \tabular{ll}{
    \code{\link{garchFit}} \tab fits the parameters of a GARCH process
  }
    
    
  \subsection{Extractor Functions:}{
    
    \tabular{ll}{
      \code{\link{residuals}}  \tab extracts residuals from a fitted \code{"fGARCH"} object\cr
      \code{\link{fitted}}     \tab extracts fitted values from a fitted \code{"fGARCH"} object\cr
      \code{\link{volatility}} \tab extracts conditional volatility from a fitted \code{"fGARCH"} object\cr
      \code{\link{coef}}       \tab extracts coefficients from a fitted \code{"fGARCH"} object\cr
      \code{\link{formula}}    \tab extracts formula expression from a fitted \code{"fGARCH"} object
    }
  }
}


\section{4 Forecasting}{

  Functions to forcecast mean and variance of GARCH and APARCH
  processes.
  
  \tabular{ll}{
    \code{\link{predict}} \tab forecasts from an object of class \code{"fGARCH"}
  }
  
}


\section{5 Standardized distributions}{

  This section contains functions to model standardized distributions.
    
  \subsection{Skew normal distribution:}{
    
    \tabular{ll}{
      \code{\link[=dnorm]{[dpqr]norm}}   \tab Normal distribution (base R)\cr
      \code{\link[=dsnorm]{[dpqr]snorm}} \tab Skew normal distribution\cr
      \code{\link[=snormFit]{snormFit}}  \tab fits parameters of Skew normal distribution
    }
  }
    
  \subsection{Skew generalized error distribution:}{
    
    \tabular{ll}{
      \code{\link[=dged]{[dpqr]ged}}   \tab Generalized error distribution\cr
      \code{\link[=psged]{[dpqr]sged}}  \tab Skew Generalized error  distribution\cr
      \code{\link[=gedFit]{gedFit}}     \tab fits parameters of Generalized error distribution\cr
      \code{\link[=sgedFit]{sgedFit}} \tab fits parameters of Skew generalized error distribution  }
  } 
    
  \subsection{Skew standardized Student-t distribution:}{
    
    \tabular{ll}{
      \code{\link[=dstd]{[dpqr]std}}  \tab Standardized Student-t distribution\cr
      \code{\link[=dsstd]{[dpqr]sstd}} \tab Skew standardized Student-t distribution\cr
      \code{\link[=stdFit]{stdFit}}  \tab fits parameters of Standardized Student-t distribution\cr
      \code{\link[=sstdFit]{sstdFit}}  \tab fits parameters of Skew standardized Student-t distribution
    }
  }
    
  \subsection{Absolute moments:}{
    \tabular{ll}{
      \code{\link{absMoments}}  \tab computes absolute moments of these distribution
    }
  }
  
}

\section{About Rmetrics}{

    The \code{fGarch} Rmetrics package is written for educational 
    support in teaching "Computational Finance and Financial Engineering" 
    and licensed under the GPL.
    
}

%  Diethelm Wuertz and Yohan Chalabi, with contributions from 
%  Michal Miklovac, Pierre Chausse, and Chris Boudt.

\keyword{package}
