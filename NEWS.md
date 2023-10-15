## CHANGES in fGarch VERSION 4022.90 (20232-??-??, svn r6333–r????)

-   added `"fGARCH"` method for 'stats::tsdiag'. The method produces
    diagnostic plot for fitted GARCH/APARCH models and computes some
    diagnostic tests. The plots can be chosen iteractively and/or via
    arguments. The test results are in the returned value. The method is
    in development in that more plots may be made available and
    additional tests included in the returned value.

-   refactored the `"fGARCH"` method for 'summary' to return an object
    from S3 class 'summary\_fGARCH' equipped with a 'print' method. The
    printout is the same as before, except that now the numbers in the
    statistics column for the residual diagnostics are aligned on the
    decimal point (previously they were left-aligned due to a buglet).

-   the `"fGARCH"` method for `fitted` was returning the data, not the
    fitted values. Fixes issue 6789, reported by Kouhei Hashinokuchi
    (hakoshie).

-   the help pages for the `"fGARCH"` methods for `fitted()` and
    `residuals()` were stating that the returned results have the same
    class as the input time series. Actually, they return numeric
    vectors. (todo?: to make the returned values as previously
    documented, `garchFit()` would need to put the original data or the
    necessary information in the fitted object, e.g., `object@fit$data`.

-   some tests were using deprecated `fBasics::.distCheck()` (notice the
    leading dot). Replaced such calls with the equivalent
    `fBasics::distCheck()`.

## CHANGES in fGarch VERSION 4022.89 (2022-11-05, from svn r6316–r6326)

-   in `absMoments`, the absolute moments for the standardized Student-t
    distribution were wrong.

-   in README, linked to the paper by Wuertz et al.

-   substantially revised the documentation and filled gaps in it.

-   removed the functions with suffix '\_orig' which were kept
    temporarilly after the bug fix in v4021.87 since there were no
    reported problems with the fix.

## CHANGES in fGarch VERSION 4021.88 (2022-09-28, svn r6276)

-   require Matrix (&gt;= 1.5-0) to avoid problems for users who have
    earlier versions of Matrix on their devices (thanks to Mikael Jagan
    for checking for not strict enough dependency on Matrix and alerting
    the maintainer).

## CHANGES in fGarch VERSION 4021.87 (2022-08-06, svn r6215–r6265)

### NEW MAINTAINER

-   Georgi N. Boshnakov

### BUG FIXES

Fixed issue 6061 raised by William Scott, who also supplied examples.

-   The quantile function, `qsnorm`, was wrong around 0.5. The error was
    in `.qsnorm`. For now its version before the fix is kept as
    `.qsnorm_orig`. Basically, branching was done w.r.t. `p = 0.5`,
    which is correct only for the symmetric case, `\xi = 1`, and should
    be `1/(1+\xi^2)` instead. More details in the source code. The error
    was affecting the central part of the distrbution with the interval
    becoming larger for `\xi` further away from 1.

-   The cdf, `psnorm`, had an error at a single point, coinciding with
    the wrong value for `p = 0.5` returned by `qsnorm(0.5)` before the
    fix. The result was that `psnorm(qsnorm(0.5))` was returning 0.5,
    falsely giving reassurance that `qsnorm(0.5)` was correct.

-   Not mentioned in issue 6061 but the same problems held for the other
    skewed distributions: `qsstd`, `psstd`, `qsged`, `psged`. The
    original versions of the relevant internal functions are kept for
    now with a suffix `_orig`, as above: `qsstd_orig`, `psstd_orig`,
    `qsged_orig`, `psged_orig`.

### Documentation

-   Edited the documentation of `"garchSpec"` and `garchSim`. It was
    somewhat incomplete and contained leftovers, apparently from old
    versions of the functions.

-   Documented the datasets. Previously the help page for them was a
    placeholder, without the names of the available datasets. There is
    no information about the time span of the data or how the returns
    were calculated.

## CHANGES in fGarch VERSION 4021.86 (2022-06-23, svn r6188)

### NEW MAINTAINER

-   Tobias Setz

### Notes

-   This is a CRAN release of version 4001.1, with trivial changes in
    ‘<span class="file">DESCRIPTION</span>’.

## CHANGES in fGarch VERSION 4001.1 (2022-06-23, svn r6184–r6185)

### NEW MAINTAINER

-   ad interim: Martin Maechler

### NEW FEATURES

-   Packages [<span
    class="pkg">timeSeries</span>](https://CRAN.R-project.org/package=timeSeries),
    [<span
    class="pkg">timeDate</span>](https://CRAN.R-project.org/package=timeDate)
    and [<span
    class="pkg">fBasics</span>](https://CRAN.R-project.org/package=fBasics)
    are no longer in `Depends`, but only in `Imports` and hence no
    longer automatically attached to the `search()` path whenever <span
    class="pkg">fGarch</span> is.

    This may require updates in your code, e.g., adding

           stopifnot(require("timeSeries"))

    as it has been done in our own <span class="pkg">fGarch</span>'s
    examples and tests.

-   `.gogarchFit()` is at least *mentioned* in the documentation.

### BUG FIXES

-   Added registration of compiled functionality for speed up and as
    good practice.

-   Removed all `Depends:` entries and checked exactly which parts of
    packages, notably <span class="pkg">fBasics</span>, <span
    class="pkg">timeDate</span>, and <span
    class="pkg">timeSeries</span>, are needed and imported only these.

-   Eliminated warning about 'length &gt; 1' character formula in
    `garchFit()`, i.e., `.garchFit()`.

-   Replaced the error-prone checking for 'class()' equality by
    'inherits(\*, &lt;class&gt;)'.

### Misc

-   Exporting practically everything seems “wrong” (according to MM):
    Several `.<some>` functions have *no* documentation and hence should
    either be (renamed and) documented or no longer be exported.

-   a `data` argument should never have a default: hence removed from
    `garchFit()`.

## CHANGES in fGarch, VERSION 3042.83.2 (2020-03-07, CRAN team)

### Misc

-   in ‘<span class="file">dist-norm.Rd</span>’, removed the description
    of argument `...`, which is not in the argument list of any function
    described there.

## CHANGES in fGarch, VERSION 3042.83.1 (2019-01-31, CRAN team)

### Misc

-   in ‘<span class="file">NAMESPACE</span>’ and ‘<span
    class="file">R/methods-plot.R</span>’ renamed functions
    `.plot.garch.1`, ..., `.plot.garch.13` to `.plot.garch_1`, ...,
    `.plot.garch_13`.

-   compressed datasets ‘<span class="file">data/dem2gbp.csv</span>’ and
    ‘<span class="file">data/sp500dge.csv</span>’ to ‘<span
    class="file">data/dem2gbp.csv.bz2</span>’ ‘<span
    class="file">data/sp500dge.csv.bz2</span>’, respectively.

## CHANGES in fGarch, VERSION 3042.83 (2017-11-16, svn r...)

### Misc

-   Startup message removed

-   Incorporate fixes by CRAN team (Brian Ripley?)

-   Checks and adaptions for R 3.4.2, e.g., ‘<span
    class="file">DESCRIPTION</span>’, ...

## CHANGES in fGarch, VERSION 3010.82.1 (2016-08-14, CRAN team.)

### Misc

-   in ‘<span class="file">NAMESPACE</span>’, import (selectively) from
    <span class="pkg">utils</span>.

-   changed a couple of calls to `Matrix()` from package <span
    class="pkg">Matrix</span> and `fastICA()` from <span
    class="pkg">fastICA</span> to the fully qualified forms
    `Matrix::Matrix()` and `fastICA::fastICA`.

-   removed some flag settings in ‘<span class="file">Makevars</span>’.

-   in ‘<span class="file">math.f</span>’, move a `DATA` command out of
    the body of an `"if"` block putting it towards the beginning of the
    file.

## CHANGES in fGarch, VERSION 3010.82 (2013-04-30, svn r5509) – and earlier

### ChangeLog

-   Changes up to April 2013, by Yohan Chalabi, Diethelm Wuertz, Pierre
    Chausse and Martin Maechler are all in file ‘<span
    class="file">ChangeLog</span>’.
