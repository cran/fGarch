## CHANGES in fGarch VERSION 4021.89 (2022-??-??, from svn r6316--rXXXX)

-   in `absMoments`, the absolute moments for the standardized Student-t
    distribution were wrong.

-   substantially revised the documentation and filled gaps in it.

## CHANGES in fGarch VERSION 4021.88 (2022-09-28, svn r6276)

-   require Matrix (\>= 1.5-0) to avoid problems for users who have
    earlier versions of Matrix on their devices (thanks to Mikael Jagan
    for checking for not strict enough dependency on Matrix and alerting
    the maintainer).

## CHANGES in fGarch VERSION 4021.87 (2022-08-06, svn r6215--r6265)

### NEW MAINTAINER

-   Georgi N. Boshnakov

### BUG FIXES

Fixed issue 6061 raised by William Scott, who also supplied examples.

-   The quantile function, `qsnorm`, was wrong around 0.5. The error was
    in `.qsnorm`. For now its version before the fix is kept as
    `.qsnorm_orig`. Basically, branching was done w.r.t. *p = 0.5*,
    which is correct only for the symmetric case, *ξ = 1*, and should be
    *1/(1+ξ\^2)* instead. More details in the source code. The error was
    affecting the central part of the distrbution with the interval
    becoming larger for *ξ* further away from 1.

-   The cdf, `psnorm`, had an error at a single point, coinciding with
    the wrong value for *p = 0.5* returned by `qsnorm(0.5)` before the
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
    '[DESCRIPTION]{.file}'.

## CHANGES in fGarch VERSION 4001.1 (2022-06-23, svn r6184--r6185)

### NEW MAINTAINER

-   ad interim: Martin Maechler

### NEW FEATURES

-   Packages
    [[timeSeries]{.pkg}](https://CRAN.R-project.org/package=timeSeries),
    [[timeDate]{.pkg}](https://CRAN.R-project.org/package=timeDate) and
    [[fBasics]{.pkg}](https://CRAN.R-project.org/package=fBasics) are no
    longer in `Depends`, but only in `Imports` and hence no longer
    automatically attached to the `search()` path whenever
    [fGarch]{.pkg} is.

    This may require updates in your code, e.g., adding

           stopifnot(require("timeSeries"))

    as it has been done in our own [fGarch]{.pkg}\'s examples and tests.

-   `.gogarchFit()` is at least *mentioned* in the documentation.

### BUG FIXES

-   Added registration of compiled functionality for speed up and as
    good practice.

-   Removed all `Depends:` entries and checked exactly which parts of
    packages, notably [fBasics]{.pkg}, [timeDate]{.pkg}, and
    [timeSeries]{.pkg}, are needed and imported only these.

-   Eliminated warning about \'length \> 1\' character formula in
    `garchFit()`, i.e., `.garchFit()`.

-   Replaced the error-prone checking for \'class()\' equality by
    \'inherits(\*, \<class\>)\'.

### Misc

-   Exporting practically everything seems "wrong" (according to MM):
    Several `.<some>` functions have *no* documentation and hence should
    either be (renamed and) documented or no longer be exported.

-   a `data` argument should never have a default: hence removed from
    `garchFit()`.

## CHANGES in fGarch, VERSION 3042.83.2 (2020-03-07, CRAN team)

### Misc

-   in '[dist-norm.Rd]{.file}', removed the description of argument
    `...`, which is not in the argument list of any function described
    there.

## CHANGES in fGarch, VERSION 3042.83.1 (2019-01-31, CRAN team)

### Misc

-   in '[NAMESPACE]{.file}' and '[R/methods-plot.R]{.file}' renamed
    functions `.plot.garch.1`, \..., `.plot.garch.13` to
    `.plot.garch_1`, \..., `.plot.garch_13`.

-   compressed datasets '[data/dem2gbp.csv]{.file}' and
    '[data/sp500dge.csv]{.file}' to '[data/dem2gbp.csv.bz2]{.file}'
    '[data/sp500dge.csv.bz2]{.file}', respectively.

## CHANGES in fGarch, VERSION 3042.83 (2017-11-16, svn r\...)

### Misc

-   Startup message removed

-   Incorporate fixes by CRAN team (Brian Ripley?)

-   Checks and adaptions for R 3.4.2, e.g., '[DESCRIPTION]{.file}', \...

## CHANGES in fGarch, VERSION 3010.82.1 (2016-08-14, CRAN team.)

### Misc

-   in '[NAMESPACE]{.file}', import (selectively) from [utils]{.pkg}.

-   changed a couple of calls to `Matrix()` from package [Matrix]{.pkg}
    and `fastICA()` from [fastICA]{.pkg} to the fully qualified forms
    `Matrix::Matrix()` and `fastICA::fastICA`.

-   removed some flag settings in '[Makevars]{.file}'.

-   in '[math.f]{.file}', move a `DATA` command out of the body of an
    `"if"` block putting it towards the beginning of the file.

## CHANGES in fGarch, VERSION 3010.82 (2013-04-30, svn r5509) -- and earlier

### ChangeLog

-   Changes up to April 2013, by Yohan Chalabi, Diethelm Wuertz, Pierre
    Chausse and Martin Maechler are all in file '[ChangeLog]{.file}'.
