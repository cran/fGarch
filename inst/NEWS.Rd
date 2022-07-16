\name{NEWS}
\title{News for \R Package \pkg{fGarch}}
\encoding{UTF-8}

\section{CHANGES in fGarch VERSION 4001.1 (2022-06-23, svn r6184)}{
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
  \subsection{Author}{
    \itemize{
      \item ad interim: Martin Maechler
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
      \item . % add by comparing the directory trees (after "untar")
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
    \itemize{
      \item . % add by comparing the directory trees (after "untar")
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
      \item . % add by comparing the directory trees (after "untar")
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
