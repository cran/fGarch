
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DESCRIPTION:
#  show.fGARCH               S4 Show method for an object of class 'fGARCH'
#  show.fGARCHSPEC           S4 Show method for an object of class 'fGARCHSPEC'
################################################################################

.prepare_GARCH_show <- function(object) {
    ## based on the original body of the show method
    common <- list(
        title = object@title,
        call = object@call,
        formula = object@formula,
        description = object@description
    )

    wrk <- if(as.character(object@call[1]) == ".gogarchFit") {
               ## multivariate
               list(
                   type = "multivariate",
                   cond_dist = object@fit[[1]]@fit$params$cond.dist,
                   number_of_margins = length(object@fit)
               )
           } else {
               ## univariate
               list(
                   type = "univariate",
                   cond_dist = object@fit$params$cond.dist,
                   par = object@fit$par,
                   se_method = if (object@fit$params$cond.dist == "QMLE")
                                   "robust"
                               else 
                                   "based on Hessian",
                   matcoef = object@fit$matcoef,
                   loglik = -object@fit$value,
                   normalized_loglik = -object@fit$value / NROW(object@data)
               )
           }
    c(common, wrk)
}

.print_title <- function(x) {
    cat("\nTitle:\n ")
    cat(x, "\n")
    invisible(NULL)
}

.print_call <- function(x) {
    cat("\nCall:\n ")
    cat(paste(deparse(x), sep = "\n", collapse = "\n"), "\n")
}

.print_cond_dist <- function(x) {
    cat("\nConditional Distribution:\n ")
    cat(x, "\n")
}

.print_coef <- function(x) {
    cat("\nCoefficient(s):\n")
    digits = max(5, getOption("digits") - 4)
    print.default(format(x, digits = digits), print.gap = 2, quote = FALSE)
}

.print_se_method <- function(x) {
    cat("\nStd. Errors:\n ")
    if (x == "QMLE")
        cat("robust", "\n")
    else
        cat("based on Hessian", "\n")
}

.print_error_analysis <- function(x) {
    digits = max(4, getOption("digits") - 5)
    signif.stars = getOption("show.signif.stars")
    cat("\nError Analysis:\n")
    printCoefmat(x, digits = digits, signif.stars = signif.stars)
}

.print_loglik <- function(x, nllh) {
    cat("\nLog Likelihood:\n ")
    cat(x, "   normalized: ", nllh, "\n")
}

.print_description <- function(x) {
    cat("\nDescription:\n ")
    cat(x, "\n")
    cat("\n")
}

.print_mean_var_eq <- function(formula) {
    cat("\nMean and Variance Equation:\n ")
    Name <- unclass(attr(formula, "data"))
    
    Formula <- formula
    attr(Formula, "data") <- NULL
    print(Formula)                # GNB: TODO: use arg. showEnv?

    cat(" [", Name, "]\n", sep = "")

    invisible(NULL)
}

.show_orig_body <- function(object, prepare = TRUE) {
    ## A function implemented by Diethelm Wuertz
    ## refactored and modified by GNB

    # Description:
    #   Print method for an object of class "fGARCH"

    # Arguments:
    #   object - an object of class 'fGARCH'

    # FUNCTION:

    res <- if(prepare)
               .prepare_GARCH_show(object)
           else
               object
    
    .print_title(res$title)
    .print_call(res$call)
    .print_mean_var_eq(res$formula)
    
    # Univariate or Multivariate Modeling ?
    if(res$type == "univariate") { # univariate Garch Models
        .print_cond_dist(res$cond_dist)
        .print_coef(res$par)
        .print_se_method(res$se_method)
        .print_error_analysis(res$matcoef)
        .print_loglik(res$loglik, res$normalized_loglik)

    } else {
        # For multivariate Garch Models ...
        #   extract information from first fitted instrument.
        object@fit[[1]]@fit$params$cond.dist
        
        # Conditional Distribution:
        cat("\nConditional Distribution:\n ")
        cat(object@fit[[1]]@fit$params$cond.dist, "\n")
        
        # Number of Margins:
        cat("\nNumber of Margins:\n ")
        cat(length(object@fit), "\n")
        
    }
    
    .print_description(res$description)

    invisible(res)
}

setMethod("show", "fGARCH",
          function(object) .show_orig_body(object)
)


# ------------------------------------------------------------------------------


setMethod(f = "show", signature(object = "fGARCHSPEC"), definition =
    function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   S4 Print Method for objects of class 'fGARCHSPEC'

    # Arguments:
    #   object - Object of class 'fGARCHSPEC'

    # FUNCTION:

    # Formula:
    x = object
    cat("\nFormula: \n ")
    cat(as.character(x@formula))

    # Model:
    cat("\nModel:")
    if (sum(abs(x@model$ar)) != 0)
        cat("\n ar:   ", x@model$ar)
    if (sum(abs(x@model$ma)) != 0)
        cat("\n ma:   ", x@model$ma)
    if (x@model$mu != 0)
        cat("\n mu:   ", x@model$mu)
    if (x@model$omega != 0)
        cat("\n omega:", x@model$omega)
    if (sum(abs(x@model$alpha)) != 0)
        cat("\n alpha:", x@model$alpha)
    if (sum(abs(x@model$gamma)) != 0)
        cat("\n gamma:", x@model$gamma)
    if (sum(abs(x@model$beta)) != 0)
        cat("\n beta: ", x@model$beta)
    if (x@model$delta != 2)
        cat("\n delta:", x@model$delta)

    # Distribution:
    cat("\nDistribution: \n ")
    cat(x@distribution)
    if (x@distribution != "norm") {
        if (x@distribution == "snorm") {
            cat("\nDistributional Parameters: \n")
            cat(" xi =", x@model$skew)
        }
        if (x@distribution == "ged" | x@distribution == "std") {
            cat("\nDistributional Parameter: \n")
            cat(" nu =", x@model$shape)
        }
        if (x@distribution == "sged" | x@distribution == "sstd") {
            cat("\nDistributional Parameters: \n")
            cat(" nu =", x@model$shape, " xi =", x@model$skew)
        }
    }

    # Seed:
    if (x@rseed != 0) {
        cat("\nRandom Seed: \n ")
        cat(x@rseed)
    }

    # Presample:
    cat("\nPresample: \n")
    n = -(length(x@presample[, 1])-1)
    time = n:0
    print(data.frame(cbind(time, x@presample)))

    # Return Value:
    invisible()
})


################################################################################
