## Author: Georgi Boshnakov

## based on tsdiag.Sarima from package 'sarima'
.tsdiag_choices <- c(
    ## ## "classic (std. residuals, acf, portmanteau p-values)",
    ## "residuals",
    ## "acf of residuals",
    ## "p values for Ljung-Box statistic",
    ## "p values for Li-McLeod statistic",
    ## "p values for Box-Pierce statistic",
    ## "pacf of residuals"

    ## use 7-13; for consistency with the plot method we would like to keep the same numbers
    ##           but this seems not possible with menu().
    ##
    ##  .plot.garch.1               Plot Time Series
    ##  .plot.garch.2               Plot Conditional SD
    ##  .plot.garch.3               Plot Series with 2 Conditional SD Superimposed
    ##  .plot.garch.4               Plot ACF of Observations
    ##  .plot.garch.5               Plot ACF of Squared Observations
    ##  .plot.garch.6               Plot Cross Correlation
    
    "Residuals",                             #  .plot.garch.7
    "Conditional SDs",                       #  .plot.garch.8
    "Standardized Residuals",                #  .plot.garch.9
    "ACF of Standardized Residuals",         #  .plot.garch.10
    "ACF of Squared Standardized Residuals", #  .plot.garch.11
    "Cross Correlation between r^2 and r",   #  .plot.garch.12
    "QQ-Plot of Standardized Residuals"      #  .plot.garch.13
    ## TODO: pacf of r and r^2
)

## tsdiag.Sarima:
##
## { # 1:  "residuals",
##     plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
##     abline(h = 0)
##     #acf(err, main = "ACF of residuals from model", lag.max = lag.max)
##     #pacf(err, main = "PACF of residuals from model", lag.max = lag.max)
## },
## { # 2:  "ACF of residuals"
##     ## acf
##     acf(err, plot = TRUE, main = "ACF of Residuals", lag.max = lag.max, na.action = na.pass)
##         # acf(cdf, main = "", lag.max = lag.max)
##         # title("ACF of" ~U[t])
##         # pacf(cdf, main = "", lag.max = lag.max)
##         # title("PACF of" ~U[t])
## },
## { # 3: "Ljung-Box p-values"
##  acftest <- acfIidTest(sacf, npar = fitdf, nlags = 1:nlag, method = "LjungBox",
##                           interval = NULL)
##     res[["LjungBox"]] <- acftest
## },
## { # 4: "Li-McLeod p-values"
##  acftest <- acfIidTest(sacf, npar = fitdf, nlags = 1:nlag, method = "LiMcLeod",
##                           interval = NULL)
##     res[["LiMcLeod"]] <- acftest
## },
## { # 5: "Box-Pierce p-values"
##  acftest <- acfIidTest(sacf, npar = fitdf, nlags = 1:nlag, method = "BoxPierce",
##                           interval = NULL)
##     res[["BoxPierce"]] <- acftest
## },
## { # 6:  "PACF of residuals"
##     ## acf
##     pacf(err, plot = TRUE, main = "PACF of Residuals", lag.max = lag.max, na.action = na.pass)
## },
##     # { # 4: "ACF/Histogram of tau_residuals"
##     #     acf(err2, main = "ACF of tau_residuals", lag.max = lag.max)
##     #     hist(err2, freq = FALSE, main = "Histogram of tau_residuals", xlab  =  "",
##     #          ylim = c(0, 0.5))
##     #     lines(seq(-5, 5, .01), dnorm(seq(-5, 5, .01)), col = "red")
##     # }

tsdiag.fGARCH <- function(object, gof.lag = NULL, ask = FALSE, ..., plot = c(4L, 5L, 7L),
                          layout = NULL)
{
    ## Georgi Boshnakov
    n_per_page <- if(is.null(layout))
                      3
                  else
                      ## length(layout[[1]])
                      do.call("layout", layout)
    
    if(is.null(gof.lag))
        gof.lag <- 20  # :TODO: NOTE: arbitrary value
    else if(!is.numeric(gof.lag))
        stop("'gof.lag' must be numeric and contain positive integers")

    lag.max <- max(gof.lag)


    sres  <- residuals(object, standardize = TRUE)
    sres2 <- sres^2

    
    choices <- .tsdiag_choices
    chnum <- 1:length(choices)
    
    if(!isTRUE(plot)){                  # plot is typically numeric index here;
        choices <- choices[plot]        # FALSE or NULL give zero length result, so no plots
        chnum <- chnum[plot]

        if(anyNA(choices)){
            warning("'plot' should be TRUE/FALSE or vector of positive integers <= ",
                    length(.tsdiag_choices), ",\n", "ignoring non-existent values")
            chnum <- chnum[!is.na(choices)]
            choices <- choices[!is.na(choices)]
        }
    }

    if(length(choices) > 0){
	old.par <- par(no.readonly = TRUE)
	on.exit(par(old.par))     # restore graphics parameters before exiting.

        ask_user <- interactive() && (ask || length(choices) > n_per_page)

        ## adjust n_per_page if 'layout' is missing
        if(is.null(layout)) {
            n_per_page <- if(ask_user)
                              ## was: layout(matrix(1:3, ncol = 1))
                              layout(matrix(1:min(3, length(choices)), ncol = 1))
                          else
                              layout(matrix(1:min(3, length(choices)), ncol = 1))
        }
            
        choice_title <- "Select a plot number or 0 to exit"
        ch_index <- if(length(choices) == 1)
                        1
                    else if(ask)
                        menu(choices, title = choice_title)
                    else if(!identical(plot, FALSE))
                        1
                    else
                        integer(0)
        choice <- chnum[ch_index]

        ## ## precompute common stuff for portmanteau tests
        ## nlag <- gof.lag
        ## pval <- numeric(nlag)
        ## fitdf <- if(inherits(object, "Sarima"))
        ##              length(object$internal$nonfixed)
        ##          else if(inherits(object, "Arima"))
        ##              sum(object$arma[1:4]) # object$arma is: p, q, p_s. q_s, s, d, d_s
        ##          else
        ##              0
        ##             # for(i in 1L:nlag)
        ##             #     pval[i] <- Box.test(err, i, type="Ljung-Box", 
        ##             #                         fitdf = ifelse(i > fitdf, fitdf, i - 1))$p.value
        ## sacf <- autocorrelations(err, maxlag = nlag) # deal with NA's?

        
        res <- list(residuals = sres)
        while(length(choice) != 0){
            switch(choice,
                   .plot.garch.7 (object),  # "Residuals",                             
                   .plot.garch.8 (object),  # "Conditional SDs",                       
                   .plot.garch.9 (object),  # "Standardized Residuals",                
                   .plot.garch.10(object),  # "ACF of Standardized Residuals",         
                   .plot.garch.11(object),  # "ACF of Squared Standardized Residuals", 
                   .plot.garch.12(object),  # "Cross Correlation between r^2 and r",   
                   .plot.garch.13(object)   # "QQ-Plot of Standardized Residuals"      
                   )
                   
            if(length(chnum) == 1)  # length(choices) == 1
                break
            if(ask_user) { # was: interactive() && (ask || length(choices) > n_per_page)
                ch_index <- menu(choices, title = choice_title)
                choice <- chnum[ch_index]
            } else{
                ## just plot the next one
                ##  Note: this doesn't update ch_index
                chnum <- chnum[-1]
                choice <- chnum[1]
            }
        }
    }

    .f <- function(x)
        c(statistic = as.vector(x$statistic), p.value = x$p.value)

    if(requireNamespace("goftest", quietly = TRUE)) {
        gofargs <- .resid_with_dist(object)
        res$gof <- rbind(
            "Anderson-Darling" = .f(do.call(goftest::ad.test, gofargs)),
            "Cramer-vonMises" = .f(do.call(goftest::cvm.test, gofargs)) )

        gofargs$estimated <- FALSE
        res$gof_composite <- rbind(
            "Anderson-Darling" = .f(do.call(goftest::ad.test, gofargs)),
            "Cramer-vonMises" = .f(do.call(goftest::cvm.test, gofargs)) )
    } else{
        message("Please install package 'goftest' for additional tests")
    }

    class(res) <- "tsdiag_fGARCH"
    
    invisible(res)
}

print.tsdiag_fGARCH <- function(x, ...){
    ## for now just drop the values of the residuals
    x <- x[names(x) != "residuals"]
    for(s in names(x)) {
        cat(paste0( "\n", s, ":", "\n"))
        print(x[[s]])
    }
}



## GNB, based on .plot.garch.13
##      for tests in 'goftest'
.resid_with_dist <- function(x) {
    sres <- residuals(x, standardize = TRUE)

    cond_dist <- x@fit$params$cond.dist
    cond_cdf <- paste("p", cond_dist, sep = "")

    parNames <- names(x@fit$par)
    skew <-
        if ("skew" %in% parNames)
            x@fit$par["skew"]
        else
            x@fit$params$skew
    shape <-
        if ("shape" %in% parNames)
            x@fit$par["shape"]
        else
            x@fit$params$shape

    res <- list(x = sres)

    res$null <- if (cond_dist == "QMLE")
                    "pnorm"
                else
                    cond_cdf

    if (cond_dist == "std" || cond_dist == "ged")
        res$nu <- shape
    else if (cond_dist == "snorm")
        res$xi <- skew
    else if (cond_dist == "sstd" || cond_dist == "sged") {
        res$xi <- skew
        res$nu <- shape
    } else if (cond_dist == "snig") {
        res$rho <- skew
        res$zeta <- shape
    }

    res$estimated <- TRUE
    res$nullname <- cond_dist

    res
}
