% Convert this to NEWS.md:
%  1. inst/NEWS.Rd to NEWS.html:
%     tools:::Rd2HTML("inst/NEWS.Rd", "NEWS.html")
%
%  2. NEWS.html to NEWS.md (https://www.r-bloggers.com/2018/12/rstudio-pandoc-html-to-markdown/)
%     rmarkdown::pandoc_convert("NEWS.html", to = "markdown_strict", output = "NEWS.md")
%
%  3. Manually remove everything up to the first line containing:
%              ## CHANGES in fGarch VERSION ...
%
%  4. Bring the markdown sections one level up, i.e. remove one '#' from
%     each section level. For example using regex ^#\(#+\) with replacement value \1

\name{NEWS}
\title{News for \R Package \pkg{fGarch}}
\encoding{UTF-8}

\section{CHANGES in fGarch VERSION 4032.91 (2024-02-02, svn r6436--r????)}{
  \itemize{

    \item added computation of Value-at-Risk (VaR) and expected
    shortfall (ES) for fitted GARCH and APARCH models (in-sample and in
    the predict method). Just use something like `Var(fitted_object)`,
    `ES(fitted_object)` or `predict(fitted_object, ..., p_loss = 0.05)`.

  }
}


\section{CHANGES in fGarch VERSION 4022.90 (2023-10-15, svn r6333--r6435)}{
  \itemize{

    \item added \code{"fGARCH"} method for `stats::tsdiag`. The method
    produces diagnostic plot for fitted GARCH/APARCH models and computes
    some diagnostic tests. The plots can be chosen interactively and/or
    via arguments. The test results are in the returned value. The
    method is in development in that more plots may be made available
    and additional tests included in the returned value.

    \item refactored the \code{"fGARCH"} method for 'summary' to return
    an object from S3 class 'summary_fGARCH' equipped with a 'print'
    method. The printout is the same as before, except that now the
    numbers in the statistics column for the residual diagnostics are
    aligned on the decimal point (previously they were left-aligned due
    to a buglet).

    \item the \code{"fGARCH"} method for \code{fitted} was returning the
    data, not the fitted values. Fixes issue 6789, reported by Kouhei
    Hashinokuchi (hakoshie).

    \item the help pages for the \code{"fGARCH"} methods for
    \code{fitted()} and \code{residuals()} were stating that the
    returned results have the same class as the input time
    series. Actually, they return numeric vectors. (todo?: to make the
    returned values as previously documented, \code{garchFit()} would
    need to put the original data or the necessary information in the
    fitted object, e.g., \code{object@fit$data}.

    \item some tests were using deprecated \code{fBasics::.distCheck()}
    (notice the leading dot). Replaced such calls with the equivalent
    \code{fBasics::distCheck()}.

  }
}


\section{CHANGES in fGarch VERSION 4022.89 (2022-11-05, from svn r6316--r6326)}{
  \itemize{
    
    \item in \code{absMoments}, the absolute moments for the
    standardized Student-t distribution were wrong.

    \item in README, linked to the paper by Wuertz et al.

    \item substantially revised the documentation and filled gaps in it.

    \item removed the functions with suffix `_orig` which were kept
    temporarilly after the bug fix in v4021.87 since there were no
    reported problems with the fix.

  }
}


\section{CHANGES in fGarch VERSION 4021.88 (2022-09-28, svn r6276)}{
  \itemize{
    
    \item require Matrix (>= 1.5-0) to avoid problems for users who have
    earlier versions of Matrix on their devices (thanks to Mikael Jagan
    for checking for not strict enough dependency on Matrix and alerting
    the maintainer).

  }
}


\section{CHANGES in fGarch VERSION 4021.87 (2022-08-06, svn r6215--r6265)}{
  \subsection{NEW MAINTAINER}{
    \itemize{
      \item Georgi N. Boshnakov
    }
  }
  \subsection{VERSION NUMBERS}{
    \itemize{

      \item We continue to use the traditional Rmetrics scheme for
      numbering the versions of the package as Mmmr.n, where 'M' is the
      current major R version at the time of submission of the package
      to CRAN, 'mm' is the minor one and 'r' is the revision. 'n' is the
      sequential number of the CRAN submission of the package. For
      example, this release has version 4021.87 since it was released when R
      4.2.1 was current and 'n' in the previous version was 86.

    }
  }
  \subsection{BUG FIXES}{
    Fixed issue 6061 raised by William Scott, who also supplied
    examples.
    \itemize{

      \item The quantile function, \code{qsnorm}, was wrong around
      0.5. The error was in \code{.qsnorm}. For now its version before
      the fix is kept as \code{.qsnorm_orig}. Basically, branching was
      done w.r.t. \eqn{p = 0.5}, which is correct only for the symmetric
      case, \eqn{\xi = 1}, and should be \eqn{1/(1+\xi^2)} instead.  More
      details in the source code. The error was affecting the central
      part of the distrbution with the interval becoming larger for
      \eqn{\xi} further away from 1.

      \item The cdf, \code{psnorm}, had an error at a single point,
      coinciding with the wrong value for \eqn{p = 0.5} returned by
      \code{qsnorm(0.5)} before the fix. The result was that
      \code{psnorm(qsnorm(0.5))} was returning 0.5, falsely giving
      reassurance that \code{qsnorm(0.5)} was correct.

      \item Not mentioned in issue 6061 but the same problems held for
      the other skewed distributions: \code{qsstd}, \code{psstd},
      \code{qsged}, \code{psged}. The original versions of the relevant
      internal functions are kept for now with a suffix \code{_orig}, as
      above: \code{qsstd_orig}, \code{psstd_orig}, \code{qsged_orig},
      \code{psged_orig}.
      
    }
  }

  \subsection{Documentation}{
    \itemize{

      \item Edited the documentation of \code{"garchSpec"} and
      \code{garchSim}. It was somewhat incomplete and contained
      leftovers, apparently from old versions of the functions.

      \item Documented the datasets. Previously the help page for them
      was a placeholder, without the names of the available datasets.
      There is no information about the time span of the data or how the
      returns were calculated.

    }
  }

}

\section{CHANGES in fGarch VERSION 4021.86 (2022-06-23, svn r6188)}{
  \subsection{NEW MAINTAINER}{
    \itemize{
      \item Tobias Setz
    }
  }
  \subsection{Notes}{
    \itemize{
      \item This is a CRAN release of version 4001.1, with trivial
      changes in \file{DESCRIPTION}.
    }
  }
}



\section{CHANGES in fGarch VERSION 4001.1 (2022-06-23, svn r6184--r6185)}{
  \subsection{NEW MAINTAINER}{
    \itemize{
      \item ad interim: Martin Maechler
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item Packages \CRANpkg{timeSeries}, \CRANpkg{timeDate} and \CRANpkg{fBasics}
      are no longer in \code{Depends}, but only in \code{Imports} and hence no
      longer automatically attached to the \code{search()} path whenever
      \pkg{fGarch} is.

      This may require updates in your code, e.g., adding \preformatted{
   stopifnot(require("timeSeries"))}
      as it has been done in our own \pkg{fGarch}'s examples and tests.

      \item \code{.gogarchFit()} is at least \emph{mentioned} in the documentation.
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item Added registration of compiled functionality for speed up and
      as good practice.

      \item Removed all \code{Depends:} entries and checked exactly which
      parts of packages, notably \pkg{fBasics}, \pkg{timeDate}, and
      \pkg{timeSeries}, are needed and imported only these.

      \item Eliminated warning about 'length > 1' character formula in
      \code{garchFit()}, i.e., \code{.garchFit()}.

      \item Replaced the error-prone checking for `class()` equality
      by `inherits(*, <class>)`.
    }
  }

  \subsection{Misc}{
    \itemize{
      \item Exporting practically everything seems \dQuote{wrong} (according to MM):
      Several \code{.<some>} functions have \emph{no} documentation and
      hence should either be (renamed and) documented or no longer be exported.

      \item a \code{data} argument should never have a default: hence removed
      from \code{garchFit()}.
    }
  }
}


\section{CHANGES in fGarch, VERSION 3042.83.2 (2020-03-07, CRAN team)}{
  \subsection{Misc}{
    \itemize{

      % MM: add by comparing the directory trees (after "untar")
      % GNB: added by looking at the commit at https://github.com/cran/fGarch/
      %
      \item in \file{dist-norm.Rd}, removed the description of argument
      \code{...}, which is not in the argument list of any function
      described there.

    }
  }
}


%% CRAN Versions of  fGarch .../CRAN/src/contrib/Archive/fGarch -- 22.Jun 2022

%%     Size         Date  tarball-filename
%%   ------ ------------  ----------------------
%%   143455  7. Mar 2020  fGarch_3042.83.2.tar.gz  Packaged: 2020-03-07 10:25:26 UTC; hornik
%%   - - - - - - - - - -  - - - - - - - - - - - -
%%   145810 31. Jan 2019  fGarch_3042.83.1.tar.gz  Packaged: 2019-01-31 16:42:21 UTC; hornik
%%   160314 16. Nov 2017  fGarch_3042.83.tar.gz    Packaged: 2017-11-16 17:28:34 UTC; Tobias Setz
%%   160484 15. Aug 2016  fGarch_3010.82.1.tar.gz  Packaged: 2016-08-15 09:24:21 UTC; ripley
%%   152081  1. May 2013  fGarch_3010.82.tar.gz    Packaged: 2013-04-30 20:07:14 UTC; yankee (Yohan Ch)

\section{CHANGES in fGarch, VERSION 3042.83.1 (2019-01-31, CRAN team)}{
  \subsection{Misc}{
    % MM: add by comparing the directory trees (after "untar")
    % GNB: added by looking at the commit at https://github.com/cran/fGarch/
    %
    \itemize{

      \item in \file{NAMESPACE} and \file{R/methods-plot.R} renamed
      functions \code{.plot.garch.1}, \ldots, \code{.plot.garch.13} to
      \code{.plot.garch_1}, \ldots, \code{.plot.garch_13}.

      \item compressed datasets \file{data/dem2gbp.csv} and
      \file{data/sp500dge.csv} to \file{data/dem2gbp.csv.bz2}
      \file{data/sp500dge.csv.bz2}, respectively.
      
    }
  }
}

\section{CHANGES in fGarch, VERSION 3042.83 (2017-11-16, svn r...)}{
  \subsection{Misc}{
    \itemize{
      \item Startup message removed
      \item Incorporate fixes by CRAN team (Brian Ripley?)
      \item Checks and adaptions for R 3.4.2, e.g., \file{DESCRIPTION}, ...
    }
  }
}

\section{CHANGES in fGarch, VERSION 3010.82.1 (2016-08-14, CRAN team.)}{
  \subsection{Misc}{
    \itemize{

      \item in \file{NAMESPACE}, import (selectively) from \pkg{utils}.

      \item changed a couple of calls to \code{Matrix()} from package
      \pkg{Matrix} and \code{fastICA()} from \pkg{fastICA} to the fully
      qualified forms \code{Matrix::Matrix()} and
      \code{fastICA::fastICA}.

      \item removed some flag settings in \file{Makevars}.

      \item in \file{math.f}, move a \code{DATA} command out of the body
      of an \code{"if"} block putting it towards the beginning of the
      file.
      
    }
  }
}

\section{CHANGES in fGarch, VERSION 3010.82 (2013-04-30, svn r5509) -- and earlier}{
  \subsection{ChangeLog}{
    \itemize{
      \item Changes up to April 2013, by Yohan Chalabi, Diethelm Wuertz, Pierre
      Chausse and Martin Maechler
      are all in file \file{ChangeLog}. % = ../ChangeLog
    }
  }
}



%% CRAN-Archive for fGarch  -- older ones

%%     Size         Date  tarball-filename
%%   ------ ------------  ----------------------
%%   152081  1. Mai 2013  fGarch_3010.82.tar.gz
%%   151484 18. Sep 2012  fGarch_2150.81.tar.gz
%%   199196  4. Jun 2012  fGarch_2110.80.1.tar.gz
%%   198030 10. Nov 2009  fGarch_2110.80.tar.gz
%%   190133 28. Sep 2009  fGarch_2100.79.tar.gz
%%   183635 20. Apr 2009  fGarch_2100.78.tar.gz
%%   183616 16. Apr 2009  fGarch_290.77.tar.gz
%%   178881 28. Jan 2009  fGarch_290.76.tar.gz
%%   177761 27. Okt 2008  fGarch_280.75.tar.gz
%%   177703 25. Okt 2008  fGarch_280.74.tar.gz
%%   177695 23. Okt 2008  fGarch_280.73.tar.gz
%%   127399  8. Okt 2007  fGarch_260.72.tar.gz
%%   127396  2. Okt 2007  fGarch_260.71.tar.gz
